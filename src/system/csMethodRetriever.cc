/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <string>
#include <fstream>
#include <cstring>
#include "csMethodRetriever.h"
#include "csException.h"
#include "csVector.h"
#include "geolib_string_utils.h"
#include "csModuleLists.h"
#include "geolib_platform_dependent.h"

#include "cseis_modules.h"
#include <dlfcn.h>

using namespace cseis_system;
using namespace cseis_geolib;

csMethodRetriever* csMethodRetriever::theinst_ = 0;

csMethodRetriever& csMethodRetriever::MethodRetriever(){
    if(!csMethodRetriever::theinst_){
        csMethodRetriever::theinst_ = new csMethodRetriever();
    }

    return *csMethodRetriever::theinst_;
}

csMethodRetriever::csMethodRetriever()
    : isOk_(false)
{
    init();
}
csMethodRetriever::~csMethodRetriever(){

}
void csMethodRetriever::init(){
    /*m_isOk = false;

    m_mlfnm = GetASModelListFileName();

    if(!File::exists(m_mlfnm.c_str())){
        m_errInfo.assign("Can't find file [");
        m_errInfo.append(m_mlfnm.c_str());
        m_errInfo.append("] !");
        return;
    }

    if(!m_modules.read(m_mlfnm.c_str(), NULL)){
        m_errInfo.assign("Load file [");
        m_errInfo.append(m_mlfnm.c_str());
        m_errInfo.append("] fails!");
        return;
    }

    m_nmodels = 0;
    for(int imodel = 0; imodel < m_modules.size(); imodel++){
        if(    (m_modules.getValue(imodel) == "EXE_SINGLE_TRACE")
            || (m_modules.getValue(imodel) == "EXE_MULTIE_TRACE")
            || (m_modules.getValue(imodel) == "EXE_FILE")
            || (m_modules.getValue(imodel) == "EXEC_TYPE_INPUT"))
            m_mnames[m_nmodels++] = m_modules.getKey(imodel).str();
    }

    m_isOk = true;*/
}

void csMethodRetriever::getParamInitMethod( std::string const& name, int verMajor, int verMinor, MParamPtr& param, MInitPtr& init ) {
  std::string nameLower = cseis_geolib::toLowerCase( name );
  
  char* soName     = new char[200];
  char const* ptr  = nameLower.c_str();
  sprintf(soName,"libas_%s.so.%d.%d",ptr,verMajor,verMinor);

  void* handle = dlopen( soName, RTLD_LAZY );
  const char *dlopen_error = dlerror();
  delete [] soName;

  if( dlopen_error ) {
    throw( cseis_geolib::csException("Error occurred while opening shared library. ...does module '%s' exist? Does version '%d.%d' exist?\nSystem message: %s\n",
                       name.c_str(), verMajor, verMinor, dlopen_error ) );
  }

  param = getParamMethod( nameLower, handle );
  init  = getInitMethod( nameLower, handle );

  //  dlclose(handle);
  //  fprintf(stderr,"Param method pointer found: %x\n", param);
}
//----------------------------------------------------------------------------
//
MParamPtr csMethodRetriever::getParamMethod( std::string const& nameLower, char const* soName ) {
  void* handle = dlopen( soName, RTLD_LAZY );
  const char *dlopen_error = dlerror();
  delete [] soName;
  if( dlopen_error ) {
    fprintf(stdout, "Error occurred while opening shared library. ...does module '%s' exist?\nSystem message: %s\n\n",
            nameLower.c_str(), dlopen_error );
    fflush(stdout);
    throw( cseis_geolib::csException("Error occurred while opening shared library. ...does module '%s' exist?\nSystem message: %s\n",
      				     nameLower.c_str(), dlopen_error ) );
  }
  return getParamMethod( nameLower, handle );
}
MParamPtr csMethodRetriever::getParamMethod( std::string const& name ) {
  std::string nameLower = cseis_geolib::toLowerCase( name );
  char* soName     = new char[200];
  char const* ptr  = nameLower.c_str();
  sprintf(soName,"libas_%s.so",ptr);
  return getParamMethod( nameLower, soName );
}
MParamPtr csMethodRetriever::getParamMethod( std::string const& name, int verMajor, int verMinor ) {
  std::string nameLower = cseis_geolib::toLowerCase( name );
  char* soName     = new char[200];
  char const* ptr  = nameLower.c_str();
  sprintf(soName,"libas_%s.so.%d.%d",ptr,verMajor,verMinor);
  return getParamMethod( nameLower, soName );
}
MParamPtr csMethodRetriever::getParamMethod( std::string const& name, std::string versionString ) {
  std::string nameLower = cseis_geolib::toLowerCase( name );
  char* soName     = new char[200];
  char const* ptr  = nameLower.c_str();
  sprintf(soName,"libas_%s.so.%s",ptr,versionString.c_str());
  return getParamMethod( nameLower, soName );
}
//----------------------------------------------------------------------------
//
MInitPtr csMethodRetriever::getInitMethod( std::string const& nameLower, void* handle  ) {
  char* methodName = new char[200];
  sprintf( methodName, "_init_mod_%s_", nameLower.c_str() );
  
  //  fprintf(stderr,"Init method name:  '%s'\n", methodName );

  //  MInitPtr method = reinterpret_cast<MInitPtr>( dlsym(handle,methodName) );
  MInitPtr method;
  void *ptr = dlsym(handle,methodName);
  memcpy(&method, &ptr, sizeof(void *));

  const char *dlsym_error = dlerror();
  delete [] methodName;
  if( dlsym_error ) {
    throw( cseis_geolib::csException("Cannot find init definition method. System message:\n%s\n", dlsym_error ) );
  }
  return method;
}
//----------------------------------------------------------------------------
//
MParamPtr csMethodRetriever::getParamMethod( std::string const& nameLower, void* handle ) {
  char* methodName = new char[200];
  sprintf( methodName, "_params_mod_%s_", nameLower.c_str() );

  //  MParamPtr method = reinterpret_cast<MParamPtr>( dlsym(handle,methodName) );
  MParamPtr method;
  void *ptr = dlsym(handle,methodName);
  memcpy(&method, &ptr, sizeof(void *));

  const char *dlsym_error = dlerror();

  delete [] methodName;
  if( dlsym_error ) {
    throw( cseis_geolib::csException("Cannot find parameter definition method. System message:\n%s\n", dlsym_error) );
  }

  return method;
}
//--------------------------------------------------------------------
//
void csMethodRetriever::getExecMethodSingleTrace( std::string const& name, int verMajor, int verMinor, MExecSingleTracePtr& exec ) {
  std::string nameLower = cseis_geolib::toLowerCase( name );

  char* soName     = new char[200];
  char const* ptr  = nameLower.c_str();
  sprintf(soName,"libas_%s.so.%d.%d",ptr,verMajor,verMinor);

  void* handle = dlopen( soName, RTLD_LAZY );
  const char *dlopen_error = dlerror();
  delete [] soName;

  if( dlopen_error ) {
    throw( cseis_geolib::csException("Error occurred while opening shared library. ...does module '%s' exist? Does version '%d.%d' exist?\nSystem message: %s\n",
                       name.c_str(), verMajor, verMinor, dlopen_error ) );
  }

  exec  = getExecMethodSingleTrace( nameLower, handle );

  //  myMethodNameList.insertEnd(nameLower);
  //  myExecList.insertEnd(exec);
}
//----------------------------------------------------------
//
void csMethodRetriever::getExecMethodMultiTrace( std::string const& name, int verMajor, int verMinor, MExecMultiTracePtr& exec ) {

  std::string nameLower = cseis_geolib::toLowerCase( name );

  char* soName     = new char[200];
  char const* ptr  = nameLower.c_str();
  sprintf(soName,"libas_%s.so.%d.%d",ptr,verMajor,verMinor);

  void* handle = dlopen( soName, RTLD_LAZY );
  const char *dlopen_error = dlerror();
  delete [] soName;

  if( dlopen_error ) {
    //    dlclose(handle);
    throw( cseis_geolib::csException("Error occurred while opening shared library. ...does module '%s' exist? Does version '%d.%d' exist?\nSystem message: %s\n",
                       name.c_str(), verMajor, verMinor, dlopen_error ) );
  }

  exec  = getExecMethodMultiTrace( nameLower, handle );

  //  myMethodNameList.insertEnd(nameLower);
  // myExecList.insertEnd(exec);
}
//
//---------------------------------------------------------
MExecSingleTracePtr csMethodRetriever::getExecMethodSingleTrace( std::string const& nameLower, void* handle  ) {
  char* methodName = new char[200];
  sprintf( methodName, "_exec_mod_%s_", nameLower.c_str() );
  
  //  fprintf(stderr,"Exec method name:  '%s'\n", methodName );

  //  MExecSingleTracePtr method = reinterpret_cast<MExecSingleTracePtr>( dlsym(handle,methodName) );
  MExecSingleTracePtr method;
  void *ptr = dlsym(handle,methodName);
  memcpy(&method, &ptr, sizeof(void *));

  const char *dlsym_error = dlerror();
  delete [] methodName;
  if( dlsym_error ) {
    throw( cseis_geolib::csException("Cannot find exec definition method. System message:\n%s\n", dlsym_error) );
  }
  return method;
}
//---------------------------------------------------------
MExecMultiTracePtr csMethodRetriever::getExecMethodMultiTrace( std::string const& nameLower, void* handle  ) {
  char* methodName = new char[200];
  sprintf( methodName, "_exec_mod_%s_", nameLower.c_str() );
  
  //  fprintf(stderr,"Exec method name:  '%s'\n", methodName );

  //  MExecMultiTracePtr method = reinterpret_cast<MExecMultiTracePtr>( dlsym(handle,methodName) );
  MExecMultiTracePtr method;
  void *ptr = dlsym(handle,methodName);
  memcpy(&method, &ptr, sizeof(void *));

  const char *dlsym_error = dlerror();
  delete [] methodName;
  if( dlsym_error ) {
    throw( cseis_geolib::csException("Cannot find exec definition method. System message:\n%s\n", dlsym_error) );
  }
  return method;
}


//--------------------------------------------------------------------
//
int csMethodRetriever::getNumStandardModules() {
  return N_METHODS;
}
//--------------------------------------------------------------------
//
std::string const* csMethodRetriever::getStandardModuleNames() {
  return NAMES;
}



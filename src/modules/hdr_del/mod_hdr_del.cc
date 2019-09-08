/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: HDR_DEL
 *
 * @author Bjorn Olofsson
 * @date   2007
 */

//*******************************************************************
//
//  Init phase
//
//*******************************************************************

void init_mod_hdr_del_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  edef->setVariables( NULL );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  csVector<std::string> valueList;
  int numLines = param->getNumLines("header");
  if( numLines < 1 ) {
    log->error("Required user parameter 'header' not found");
  }
  csVector<std::string> tmpList;
  for( int iline = 0; iline < numLines; iline++ ) {
    param->getAll( "header", &tmpList, iline );
    for( int i = 0; i < tmpList.size(); i++ ) {
      valueList.insertEnd( tmpList.at(i) );
    }
    tmpList.clear();
  }

  string text;
  bool isModeInclude = true;
  if( param->exists( "mode" ) ) {
    param->getString( "mode", &text );
    if( !text.compare("include") ) {
     isModeInclude = true;
    }
    else if( !text.compare("exclude") ) {
     isModeInclude = false;
    }
    else {
     log->error( "Mode option not recognised: '%s'", text.c_str() );
    }
  }

  tmpList.clear();
  for( int i = valueList.size()-1; i >= 0; i-- ) {
    string name = valueList.at(i);
    int length = name.length();
    if( !hdef->headerExists( name ) ) {
      bool found = false;
      // Check for wild card * at end:
      if( name.at(0) == '*' ) {
        for( int ihdr = 0; ihdr < hdef->numHeaders(); ihdr++ ) {
          int pos = hdef->headerName(ihdr).length() - (length-1);
          if( pos < 0 ) continue;
          if( !hdef->headerName(ihdr).substr(pos).compare(name.substr(1)) ) {
            tmpList.insertEnd(hdef->headerName(ihdr));
            found = true;
          }
        }
        valueList.remove(i);
      }
      else if( name.at(length-1) == '*' ) {
        for( int ihdr = 0; ihdr < hdef->numHeaders(); ihdr++ ) {
          if( !hdef->headerName(ihdr).substr(0,length-1).compare(name.substr(0,length-1)) ) {
            tmpList.insertEnd(hdef->headerName(ihdr));
            found = true;
          }
        }
        valueList.remove(i);
      }
      if( !found ) {
        log->line("Error: Unknown trace header '%s'", name.c_str() );
        env->addError();
      }
    }
  }
  for( int i = 0; i < tmpList.size(); i++ ) {
    valueList.insertEnd( tmpList.at(i) );
  }

  if( isModeInclude ) {
    for( int i = 0; i < valueList.size(); i++ ) {
      string hdrName = valueList.at(i);
      if( hdef->isSystemTraceHeader( hdrName ) ) {
        log->line("Trace header '%s' is a system trace header which cannot be deleted", hdrName.c_str());
        env->addError();
      }
      hdef->deleteHeader( hdrName );
    }
    log->line("List of trace headers to delete:");
  }
  else {
    for( int ihdr = 0; ihdr < hdef->numHeaders(); ihdr++ ) {
      bool deleteHdr = true;
      string hdrName = hdef->headerName( ihdr );
      if( hdef->isSystemTraceHeader(hdrName) ) continue;
      for( int i = 0; i < valueList.size(); i++ ) {
        string name = valueList.at(i);
        if( !name.compare( hdrName ) ) {
          deleteHdr = false;
          valueList.remove(i);
          break;
        }
      }
      if( deleteHdr ) {
        hdef->deleteHeader( hdrName );
      }
    }
    log->line("List of trace headers to keep:");
  }
  for( int i = 0; i < valueList.size(); i++ ) {
    log->line(" %s", valueList.at(i).c_str());
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_hdr_del_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  csExecPhaseDef* edef = env->execPhaseDef;

  if( edef->isCleanup()){
    return true;
  }

  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_hdr_del_( csParamDef* pdef ) {
  pdef->setModule( "HDR_DEL", "Delete trace headers", "Delete one or more trace headers" );

  pdef->addParam( "header", "List of header names to delete", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "List of trace header names", "Place a star  *  at the beginning or end of the header name as a wild card (* in arbitrary place not supported yet)" );

  pdef->addParam( "mode", "Mode of selection", NUM_VALUES_FIXED );
  pdef->addValue( "include", VALTYPE_OPTION );
  pdef->addOption( "include", "Kill specified trace headers" );
  pdef->addOption( "exclude", "Kill all trace headers EXCEPT the ones specified." );
}

extern "C" void _params_mod_hdr_del_( csParamDef* pdef ) {
  params_mod_hdr_del_( pdef );
}
extern "C" void _init_mod_hdr_del_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_hdr_del_( param, env, log );
}
extern "C" bool _exec_mod_hdr_del_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_hdr_del_( trace, port, env, log );
}


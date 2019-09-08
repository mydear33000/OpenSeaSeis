/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>
#include "csToken.h"
#include "csVector.h"
#include "csEquationSolver.h"

using std::string;
using namespace cseis_geolib;
using namespace cseis_system;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: HDR_MATH
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_hdr_math {
  struct VariableStruct {
    int nEquations;
    csEquationSolver* solver;
    int*          indexHdrs;
    type_t*       typeHdrs;
    std::string*  nameHdrs;
    int**         constHdrIndex;
    type_t**      constHdrType;
    std::string** constNames;
    int*          nHdrs;
    csVector<string>** constList;
  };
  static const type_t SPECIAL_TYPE_TRACE = 255;
}

using namespace mod_hdr_math;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_hdr_math_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  csVector<std::string> valueList;

  //---------------------------------------------------------
  // Initialize
  //
  vars->nEquations = 0;
  vars->solver     = NULL;
  vars->constList  = NULL;
  vars->indexHdrs  = NULL;
  vars->nameHdrs   = NULL;
  vars->typeHdrs   = NULL;
  vars->constHdrIndex  = NULL;
  vars->constNames     = NULL;
  vars->constHdrType   = NULL;
  vars->nHdrs          = NULL;
  
  //---------------------------------------------------------
  // Create new headers


  int nLines = param->getNumLines( "new" );
  for( int i = 0; i < nLines; i++ ) {
    valueList.clear();
    param->getAll( "new", &valueList, i );
    int numValues = valueList.size();
    try {
      if( numValues < 1 || numValues > 4 ) {
        log->line("Wrong number of arguments for user parameter 'new'. Expected: 1 to 4, found: %d.", numValues);
        for( int i = 0; i < numValues; i++ ) {
          log->line("Argument %d: %s\n", i+1, valueList.at(i).c_str() );
        }
        env->addError();
      }
      else {
        std::string name = valueList.at(0);
        if( hdef->headerExists( name ) ) {
          log->line("Cannot create new trace header '%s'. Header already exists.", name.c_str());
          env->addError();
        }
        else if( numValues == 1 ) {  // Standard header, name only
          hdef->addStandardHeader( name );
        }
        else {
          string typeText = valueList.at(1);
          type_t type = TYPE_UNKNOWN;
          if( !typeText.compare("int") ) {
            type = TYPE_INT;
          }
          else if( !typeText.compare("float") ) {
            type = TYPE_FLOAT;
          }
          else if( !typeText.compare("double") ) {
            type = TYPE_DOUBLE;
          }
          else if( !typeText.compare("int64") ) {
            type = TYPE_INT64;
          }
          else if( !typeText.compare("string") ) {
            type = TYPE_STRING;
          }
          else {
            log->line("Unknown trace header type: '%s', given trace header: '%s'", typeText.c_str(), name.c_str() );
            env->addError();
          }
          if( type == TYPE_STRING ) {
            int maxStrLen = -1;
            if( numValues > 2 ) {
              maxStrLen = atoi(valueList.at(2).c_str());
            }
            if( maxStrLen <= 0 || maxStrLen > 99 ) {
              log->line("For string headers, the maximum length of the string must be supplied\n(e.g.  new str_hdr_name string max_len \"String header description\" ).\n");
              env->addError();
            }
            else if( numValues > 3 ) {
              string desc = valueList.at(3);
              hdef->addHeader( type, name, desc, maxStrLen );
            }
            else {
              hdef->addHeader( type, name, maxStrLen );
            }
          }
          else if( numValues == 2 ) {  // No description given
            hdef->addHeader( type, name );
          }
          else {
            string desc = valueList.at(3);
            hdef->addHeader( type, name, desc );
          }
        }
      }
    }
    catch( csException& exc ) {
      log->line("Error: %s", exc.getMessage());
      log->line(" ...to add a standard header, use the following syntax:   new <standard_hdr_name>");
      env->addError();
    }
  }  // End for loop

  //---------------------------------------------------------
  // Prepare header equations
  //
  int nEquations = param->getNumLines( "equation" );

  vars->nEquations = nEquations;
  vars->solver     = new csEquationSolver[nEquations];
  vars->constList  = new csVector<string>*[nEquations];
  // Resultant headers:
  vars->indexHdrs  = new int[nEquations];
  vars->nameHdrs   = new string[nEquations];
  vars->typeHdrs   = new type_t[nEquations];
  // Equation 'constants':
  vars->constHdrIndex  = new int*[nEquations];
  vars->constNames     = new string*[nEquations];
  vars->constHdrType   = new type_t*[nEquations];
  vars->nHdrs          = new int[nEquations];

  for( int ieq = 0; ieq < nEquations; ieq++ ) {
    vars->constList[ieq] = NULL;
    vars->constHdrIndex[ieq] = NULL;
    vars->constNames[ieq]    = NULL;
    vars->constHdrType[ieq]  = NULL;
  }
  
  // List all header names. These may actually not really exist if they were deleted
  // --> Check later that these headers exist
  int numHeaders = hdef->numHeaders();
  string* allHeaderNames = new string[numHeaders];
  for( int ihdr = 0; ihdr < numHeaders; ihdr++ ) {
    allHeaderNames[ihdr] = hdef->headerName(ihdr);
  }
  int nHdrs;

  for( int ieq = 0; ieq < nEquations; ieq++ ) {
    param->getAll( "equation", &valueList, ieq );
    if( valueList.size() != 2 ) {
      log->line("Error: Wrong number of arguments for user parameter '%s'. Expected: 2, found: %d.", "equation", valueList.size());
      for( int i = 0; i < valueList.size(); i++ ) {
        log->line("Argument %d: %s", i+1, valueList.at(i).c_str() );
      }
      env->addError();
    }
    std::string name = valueList.at(0);
    if( hdef->headerExists( name ) ) {
      std::string equationText = valueList.at(1);
      if( equationText.size() == 0 ) {
        log->line("Error: Empty equation found for trace header '%s'", name.c_str());
        env->addError();
      }
      vars->nameHdrs[ieq]  = name;
      vars->typeHdrs[ieq]  = hdef->headerType(name.c_str());
      vars->indexHdrs[ieq] = hdef->headerIndex(name.c_str());
      if( vars->typeHdrs[ieq] == TYPE_STRING ) {
        vars->nHdrs[ieq] = 1;
        vars->constNames[ieq] = new string[1];
        vars->constNames[ieq][0] = equationText;
        continue;
      }
      vars->constList[ieq] = new csVector<string>;
      if( !vars->solver[ieq].prepare( equationText, allHeaderNames, numHeaders ) ) {
        log->error("Error occurred: %s", vars->solver[ieq].getErrorMessage().c_str() );
      }
      vars->solver[ieq].prepareUserConstants( vars->constList[ieq] );
      nHdrs = vars->constList[ieq]->size();
      vars->nHdrs[ieq] = nHdrs;
      vars->constNames[ieq]    = new string[nHdrs];
      vars->constHdrIndex[ieq] = new int[nHdrs];
      vars->constHdrType[ieq]  = new type_t[nHdrs];
      
      for( int ihdr = 0; ihdr < nHdrs; ihdr++ ) {
        vars->constNames[ieq][ihdr] = vars->constList[ieq]->at(ihdr);
      }
      // Final check if all identifiers really exist as header names
      for( int ihdr = 0; ihdr < vars->constList[ieq]->size(); ihdr++ ) {
        if( !hdef->headerExists( vars->constList[ieq]->at(ihdr) ) ) {
          log->line("Unknown identifier in equation: '%s'", vars->constList[ieq]->at(ihdr).c_str() );
          env->addError();
        }
        else {
          vars->constHdrIndex[ieq][ihdr] = hdef->headerIndex( vars->constNames[ieq][ihdr].c_str() );
          vars->constHdrType[ieq][ihdr]  = hdef->headerType( vars->constNames[ieq][ihdr].c_str() );
        }
      }
    }
    else {
      log->line("Error occurred: Trace header '%s' unknown. Please create first before setting this header.", name.c_str());
      env->addError();
    }
  }
  delete [] allHeaderNames;

  nLines = param->getNumLines( "delete" );
  csVector<std::string> nameList(2);
  for( int i = 0; i < nLines; i++ ) {
    valueList.clear();
    param->getAll( "delete", &valueList, i );
    for( int ihdr = 0; ihdr < valueList.size(); ihdr++ ) {
      std::string name = valueList.at(ihdr);
      std::string desc;
      bool contains = false;
      for( int ilist = 0; ilist < nameList.size(); ilist++ ) {
        if( !name.compare(nameList.at(ilist)) ) {
          contains = true;
          break;
        }
      }
      if( contains ) {
        log->warning("Trace header '%s' has already been deleted.", name.c_str());
      }
      else if( hdef->headerExists( name ) ) {
        hdef->deleteHeader( name );
      }
      else {
        log->line("Error occurred: Trace header '%s' unknown.", name.c_str());
        env->addError();
      }
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_hdr_math_(
                        csTrace* trace,
                        int* port,
                        csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;

  if( edef->isCleanup() ) {
    if( vars->solver != NULL ) {
      delete [] vars->solver; vars->solver = NULL;
    }
    // Resultant headers:
    if( vars->indexHdrs != NULL ) {
      delete [] vars->indexHdrs; vars->indexHdrs = NULL;
    }
    if( vars->nameHdrs != NULL ) {
      delete [] vars->nameHdrs;  vars->nameHdrs = NULL;
    }
    if( vars->typeHdrs != NULL ) {
      delete [] vars->typeHdrs; vars->typeHdrs = NULL;
    }
    // Equation 'constants':
    if( vars->nHdrs != NULL ) {
      delete [] vars->nHdrs; vars->nHdrs = NULL;
    }
    for( int ieq = 0; ieq < vars->nEquations; ieq++ ) {
      if( vars->constList[ieq] != NULL ) delete vars->constList[ieq];
      if( vars->constNames[ieq] != NULL ) delete [] vars->constNames[ieq];
      if( vars->constHdrIndex != NULL ) delete [] vars->constHdrIndex[ieq];
      if( vars->constHdrType != NULL ) delete [] vars->constHdrType[ieq];
    }
    if( vars->constList != NULL ) {
      delete [] vars->constList; vars->constList = NULL;
    }
    if( vars->constNames != NULL ) {
      delete [] vars->constNames;     vars->constNames  =  NULL;
    }
    if( vars->constHdrType != NULL ) {
      delete [] vars->constHdrType; vars->constHdrType = NULL;
    }
    if( vars->constHdrIndex != NULL ) {
      delete [] vars->constHdrIndex; vars->constHdrIndex = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }
  //----------------------------

  csTraceHeader* trcHdr = trace->getTraceHeader();

  try { 
    for( int ieq = 0; ieq < vars->nEquations; ieq++ ) {
      if( vars->typeHdrs[ieq] == TYPE_STRING ) {
        trcHdr->setStringValue( vars->indexHdrs[ieq], vars->constNames[ieq][0] );
        continue;
      }
      int nHdrs = vars->nHdrs[ieq];
      string hdrName;
      double* hdrValues = NULL;
      if( nHdrs > 0 ) {
        hdrValues = new double[nHdrs];
        for( int ihdr = 0; ihdr < nHdrs; ihdr++ ) {
          int indexHdr = vars->constHdrIndex[ieq][ihdr];
          type_t type  = vars->constHdrType[ieq][ihdr];
          if( type == TYPE_FLOAT ) {
            hdrValues[ihdr] = (double)trcHdr->floatValue(indexHdr);
          }
          else if( type == TYPE_DOUBLE ) {
            hdrValues[ihdr] = trcHdr->doubleValue( indexHdr );
          }
          else if( type == TYPE_INT ) {
            hdrValues[ihdr] = (double)trcHdr->intValue(indexHdr);
          }
          else if( type == TYPE_INT64 ) {
            hdrValues[ihdr] = (double)trcHdr->int64Value(indexHdr);
          }
        }
        vars->solver[ieq].setUserConstants( hdrValues, nHdrs );
        delete [] hdrValues;
      } // END if nHdrs
      double result = vars->solver[ieq].solve();
      type_t type   = vars->typeHdrs[ieq];
      if( type == TYPE_FLOAT ) {
        trcHdr->setFloatValue( vars->indexHdrs[ieq], (float)result );
      }
      else if( type == TYPE_DOUBLE ) {
        trcHdr->setDoubleValue( vars->indexHdrs[ieq], result );
      }
      else if( type == TYPE_INT ) {
        trcHdr->setIntValue( vars->indexHdrs[ieq], (int)result );
      }
      else if( type == TYPE_INT64 ) {
        trcHdr->setInt64Value( vars->indexHdrs[ieq], (csInt64_t)result );
      }
    }
  }
  catch( EquationException& e ) {
    log->line("Error occurred: %s", e.getMessage());
    env->addError();
  }

  return true;
}

//********************************************************************************
//
//
void params_mod_hdr_math_( csParamDef* pdef ) {
  pdef->setModule( "HDR_MATH", "Trace header computation",
                   "Set trace headers, perform mathematical computations based on trace headers, create trace headers" );

  pdef->addParam( "new", "Create new trace header", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "New trace header name" );
  pdef->addValue( "", VALTYPE_OPTION, "Trace header type" );
  pdef->addOption( "int", "Integer, 4 byte integer type" );
  pdef->addOption( "int64", "Integer, 8 byte integer type" );
  pdef->addOption( "float", "Float, 4 byte floating point type" );
  pdef->addOption( "double", "Double, 8 byte floating point type" );
  pdef->addOption( "string", "String type. The number of characters allocated for the string header is given in the next user parameter",
                   "Note: Once a string trace header has been created, its length cannot be increased.");
  pdef->addValue( "", VALTYPE_STRING, "Trace header description", "For string headers, specify the number of string characters here." );
  pdef->addValue( "", VALTYPE_STRING, "Trace header description (for string headers only)" );

  pdef->addParam( "equation", "Mathematical equation", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addValue( "", VALTYPE_STRING, "Mathematical equation (or constant text string for string header).",
                  "Constants: pi,e. Functions: abs,acos,asin,atan,atan2,ceil,cos,cosh,exp,floor,log,log10,max,min,mod,pow,int,round,sin,sinh,sqrt,tan,tanh,todegrees,toradians,sign");

  /*
  pdef->addParam( "tmp", "Create temporary variable", NUM_VALUES_VARIABLE, "Same as 'new' but only creating a temporary variable that can be used within this instance of the module HDR_MATH" );
  pdef->addValue( "", VALTYPE_STRING, "Name of temporary variable" );
  pdef->addValue( "", VALTYPE_OPTION, "Number type of temporary variable" );
  pdef->addOption( "int", "Integer, 4 byte integer type" );
  pdef->addOption( "int64", "Integer, 8 byte integer type" );
  pdef->addOption( "float", "Float, 4 byte floating point type" );
  pdef->addOption( "double", "Double, 8 byte floating point type" );
  */
}

extern "C" void _params_mod_hdr_math_( csParamDef* pdef ) {
  params_mod_hdr_math_( pdef );
}
extern "C" void _init_mod_hdr_math_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_hdr_math_( param, env, log );
}
extern "C" bool _exec_mod_hdr_math_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_hdr_math_( trace, port, env, log );
}



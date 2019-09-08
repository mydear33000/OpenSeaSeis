/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csSelectionManager.h"
#include "csSelection.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SELECT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_select {
  struct VariableStruct {
    csSelectionManager* selectionManager;
  };
}
using mod_select::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_select_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->selectionManager = NULL;

  //---------------------------------------------------------
  // Create new headers
  int nLines = param->getNumLines( "header" );
  if( nLines > 1 ) {
    log->line( "More than one line encountered for user parameter '%s'. Only one line is supported.", "header" );
    env->addError();
  }

  csVector<std::string> valueList;
  param->getAll( "header", &valueList );

  if( valueList.size() == 0 ) {
    log->line("Wrong number of arguments for user parameter '%s'. Expected: >0, found: %d.", "header", valueList.size());
    env->addError();
  }
  else {
    std::string text;
    param->getString( "select", &text );

    vars->selectionManager = NULL;
    try {
      vars->selectionManager = new csSelectionManager();
      vars->selectionManager->set( &valueList, &text, hdef );
    }
    catch( csException& e ) {
      vars->selectionManager = NULL;
      log->error( "%s", e.getMessage() );
    }
    if( edef->isDebug() ) vars->selectionManager->dump();
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_select_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;

  if( edef->isCleanup() ) {
    if( vars->selectionManager ) {
      delete vars->selectionManager; vars->selectionManager = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }
  
  if( vars->selectionManager->contains( trace->getTraceHeader() ) ) {
    return true;
  }
  else {
    return false;
  }

}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_select_( csParamDef* pdef ) {
  pdef->setModule( "SELECT", "Select traces", "Select traces that match specified header selection" );

  pdef->addParam( "header", "Names of trace headers used for trace selection", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "select", "Selection of header values", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING,
    "List of selection strings, one for each specified header. See documentation for more detailed description of selection syntax",
    "Examples: '10,13' '10-13'-->'10,11,12,13' '10-30(10)'-->'10,20,30'. Operators: <,<=,>=,>,!" );
}

extern "C" void _params_mod_select_( csParamDef* pdef ) {
  params_mod_select_( pdef );
}
extern "C" void _init_mod_select_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_select_( param, env, log );
}
extern "C" bool _exec_mod_select_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_select_( trace, port, env, log );
}


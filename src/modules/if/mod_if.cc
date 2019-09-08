/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csSelectionManager.h"
#include "csSelection.h"

using namespace cseis_geolib;
using namespace cseis_system;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: IF
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_if {
  struct VariableStruct {
    csSelectionManager* selectionManager;
  };
}
using mod_if::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//*************************************************************************************************
void init_mod_if_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  vars->selectionManager = NULL;

  csVector<std::string> valueList;

  param->getAll( "header", &valueList );
  if( valueList.size() == 0 ) {
    log->warning("%s: Wrong number of parameters for option 'HEADER'. Expected: > 0, found: %d.", edef->moduleName().c_str(), valueList.size());
    env->addError();
  }

  std::string text;
  param->getString( "select", &text );

  try {
    vars->selectionManager = new csSelectionManager();
    vars->selectionManager->set( &valueList, &text, hdef );
  }
  catch( csException& e ) {
    vars->selectionManager = NULL;
    log->error( "%s: %s", edef->moduleName().c_str(), e.getMessage() );
  }
  if( edef->isDebug() ) vars->selectionManager->dump();

  env->execPhaseDef->setExecType( EXEC_TYPE_SINGLETRACE );
}

//*************************************************************************************************
// Exec phase
//
//
//*************************************************************************************************
bool exec_mod_if_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;

  if( edef->isCleanup() ) {
    if( vars->selectionManager ) delete vars->selectionManager; vars->selectionManager = NULL;
    delete vars; vars = NULL;
    return true;
  }

  if( vars->selectionManager->contains( trace->getTraceHeader() ) ) {
    *port = 0;
  }
  else {
    *port = 1;
  }

  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_if_( csParamDef* pdef ) {
  pdef->setModule( "IF", "If statement", "Branch traces that match specified header selection" );

  pdef->addParam( "header", "Names of trace headers used for trace selection", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "select", "Selection of header values", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING,
      "List of selection strings, one for each specified header. See documentation for more detailed description of selection syntax" );
}

extern "C" void _params_mod_if_( csParamDef* pdef ) {
  params_mod_if_( pdef );
}
extern "C" void _init_mod_if_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_if_( param, env, log );
}
extern "C" bool _exec_mod_if_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_if_( trace, port, env, log );
}


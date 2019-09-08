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
 * Module: KILL
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_kill {
  struct VariableStruct {
    std::string* hdrNames;
    csSelectionManager* selectionManager;
    int nHeaders;
    bool isModeInclude;
    bool killZeroTraces;
    bool killOneZeroTrace;
  };
}
using mod_kill::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_kill_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->killZeroTraces  = false;
  vars->killOneZeroTrace = false;
  vars->selectionManager = NULL;
  vars->hdrNames = NULL;
  vars->nHeaders = 0;
  vars->isModeInclude = true;

  //---------------------------------------------------------
  std::string text;

  if( param->exists( "zero_traces" ) ) {
    param->getString( "zero_traces", &text );
    if( !text.compare("yes") ) {
      vars->killZeroTraces = true;
    }
    else if( !text.compare("no") ) {
      vars->killZeroTraces = false;
    }
    else if( !text.compare("one") ) {
      vars->killZeroTraces = true;
      vars->killOneZeroTrace = true;
    }
    else {
      log->error( "Option not recognised: '%s'", text.c_str() );
    }
  }

  if( !vars->killZeroTraces ) {
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
      vars->nHeaders = valueList.size();
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

  if( param->exists( "mode" ) ) {
    param->getString( "mode", &text );
    if( !text.compare("include") ) {
      vars->isModeInclude = true;
    }
    else if( !text.compare("exclude") ) {
      vars->isModeInclude = false;
    }
    else {
      log->error( "Mode option not recognised: '%s'", text.c_str() );
    }
  }

  if( edef->isDebug() ) log->line("Mode include: %d", vars->isModeInclude );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_kill_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup() ) {
    if( vars->selectionManager ) {
      delete vars->selectionManager; vars->selectionManager = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  if( !vars->killZeroTraces ) {
    if( vars->selectionManager->contains( trace->getTraceHeader() ) ) {
      return !vars->isModeInclude;
    }
    else {
      return vars->isModeInclude;
    }
  }
  else {
    float* samples = trace->getTraceSamples();
    if( !vars->killOneZeroTrace ) {
      for( int i = 0; i < shdr->numSamples; i++ ) {
        if( samples[i] != 0.0 ) return vars->isModeInclude;
      }
      return !vars->isModeInclude;
    }      
    else {
      for( int i = 0; i < shdr->numSamples; i++ ) {
        if( samples[i] == 0.0 ) {
          return !vars->isModeInclude;// false;
        }
      }
      return vars->isModeInclude;
    }
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//
//*************************************************************************************************
void params_mod_kill_( csParamDef* pdef ) {
  pdef->setModule( "KILL", "Kill traces", "Kill traces that match specified header selection" );

  pdef->addParam( "header", "Names of trace headers used for trace selection", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "select", "Selection of header values", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING,
      "List of selection strings, one for each specified header. See documentation for more detailed description of selection syntax" );

  pdef->addParam( "mode", "Mode of selection", NUM_VALUES_FIXED );
  pdef->addValue( "include", VALTYPE_OPTION );
  pdef->addOption( "include", "Kill traces matching the specified selection criteria" );
  pdef->addOption( "exclude", "Kill traces NOT matching the specified selection criteria" );

  pdef->addParam( "zero_traces", "Kill zero traces", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not kill zero traces" );
  pdef->addOption( "yes", "Kill zero traces" );
  pdef->addOption( "one", "Kill traces that contains at least one zero" );
}

extern "C" void _params_mod_kill_( csParamDef* pdef ) {
  params_mod_kill_( pdef );
}
extern "C" void _init_mod_kill_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_kill_( param, env, log );
}
extern "C" bool _exec_mod_kill_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_kill_( trace, port, env, log );
}


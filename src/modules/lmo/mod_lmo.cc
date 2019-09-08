/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_lmo {
  struct VariableStruct {
    int dummy;
  };
}
using mod_lmo::VariableStruct;

//*************************************************************************************************
// Init phase
//
void init_mod_lmo_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  //  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );
  log->error("Sorry this module has not been implemented yet");
}

//*************************************************************************************************
// Exec phase
//
bool exec_mod_lmo_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  //osSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    delete vars; vars = NULL;
    return true;
  }
  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_lmo_( csParamDef* pdef ) {
  pdef->setModule( "LMO", "Linear moveout correction (NOT IMPLEMENTED YET)", "...in the widest sense (also called HMO as in 'hyperbolic moveout')." );
}

extern "C" void _params_mod_lmo_( csParamDef* pdef ) {
  params_mod_lmo_( pdef );
}
extern "C" void _init_mod_lmo_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_lmo_( param, env, log );
}
extern "C" bool _exec_mod_lmo_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_lmo_( trace, port, env, log );
}


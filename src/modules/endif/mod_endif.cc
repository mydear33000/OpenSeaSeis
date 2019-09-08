/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: ELSEIF
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
using namespace cseis_system;

void init_mod_endif_( csParamManager* userParams, csInitPhaseEnv* env, csLogWriter* log ) {
  env->execPhaseDef->setExecType( EXEC_TYPE_SINGLETRACE );
}
bool exec_mod_endif_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return true;
}
void params_mod_endif_( csParamDef* pdef ) {
  pdef->setModule( "ENDIF", "Endif statement", "Marks the end of an If-elseif-else-endif block" );
}

extern "C" void _params_mod_endif_( csParamDef* pdef ) {
  params_mod_endif_( pdef );
}
extern "C" void _init_mod_endif_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_endif_( param, env, log );
}
extern "C" bool _exec_mod_endif_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_endif_( trace, port, env, log );
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: ELSE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */

using namespace cseis_system;

void init_mod_else_( csParamManager* userParams, csInitPhaseEnv* env, csLogWriter* log ) {
  env->execPhaseDef->setExecType( EXEC_TYPE_SINGLETRACE );
}
bool exec_mod_else_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return true;
}
void params_mod_else_( csParamDef* pdef ) {
  pdef->setModule( "ELSE", "Else statement", "Branch remaining traces from if-elseif..-endif block" );
}

extern "C" void _params_mod_else_( csParamDef* pdef ) {
  params_mod_else_( pdef );
}
extern "C" void _init_mod_else_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_else_( param, env, log );
}
extern "C" bool _exec_mod_else_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_else_( trace, port, env, log );
}


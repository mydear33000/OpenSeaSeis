/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: ENDSPLIT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */

using namespace cseis_system;

void init_mod_endsplit_( csParamManager* userParams, csInitPhaseEnv* env, csLogWriter* log ) {
  env->execPhaseDef->setExecType( EXEC_TYPE_SINGLETRACE );
}
bool exec_mod_endsplit_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  if( env->execPhaseDef->isCleanup() ) {
    return true;
  }
  return false;
}
void params_mod_endsplit_( csParamDef* pdef ) {
  pdef->setModule( "ENDSPLIT", "Endsplit statement", "Marks the end of a Split-endsplit block" );
}

extern "C" void _params_mod_endsplit_( csParamDef* pdef ) {
  params_mod_endsplit_( pdef );
}
extern "C" void _init_mod_endsplit_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_endsplit_( param, env, log );
}
extern "C" bool _exec_mod_endsplit_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_endsplit_( trace, port, env, log );
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cstddef>
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: FFT_2D
 *
 * @date   2011
 */

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_fft_2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;

  edef->setExecType( EXEC_TYPE_MULTITRACE );
  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  log->error("This is a dummy module which does not function");
}

//*******************************************************************************************
// Exec phase
//*******************************************************************************************

void exec_mod_fft_2d_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* log )
{

}

//*************************************************************************************************
// Parameter definition
//*************************************************************************************************
void params_mod_fft_2d_( csParamDef* pdef ) {
  pdef->setModule( "FFT_2D", "2D FFT - DUMMY MODULE");
  pdef->setVersion(0,5);
}

extern "C" void _params_mod_fft_2d_( csParamDef* pdef ) {
  params_mod_fft_2d_( pdef );
}
extern "C" void _init_mod_fft_2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_fft_2d_( param, env, log );
}
extern "C" void _exec_mod_fft_2d_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_fft_2d_( traceGather, port, numTrcToKeep, env, log );
}

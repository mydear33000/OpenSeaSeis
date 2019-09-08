/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: OVERLAP
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_overlap {
  struct VariableStruct {
    int hdrID_time_samp1_s;
    int hdrID_time_samp1_us;
    bool isAbsTime;
    bool isFirstCall;
    int overlap_ms;
    int overlap_numSamples;
    int numSamples_in;
  };
}
using namespace mod_overlap;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_overlap_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_MULTITRACE );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 3 );

  vars->hdrID_time_samp1_s  = -1;
  vars->hdrID_time_samp1_us = -1;
  vars->isAbsTime          = true;
  vars->isFirstCall         = true;
  vars->overlap_ms          = 0;
  vars->overlap_numSamples  = 0;
  vars->numSamples_in       = 0;

  //---------------------------------------------------------
  if( param->exists("absolute_time") ) {
    std::string text;
    param->getString( "absolute_time", &text );
    if( !text.compare( "yes" ) ) {
      vars->isAbsTime = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->isAbsTime = false;
    }
    else {
      log->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }
  //---------------------------------------------------------

  param->getInt( "overlap", &vars->overlap_ms );
  int traceLength_ms = (int)( shdr->numSamples * shdr->sampleInt + 0.5 );
  if( vars->overlap_ms <= 0 || vars->overlap_ms > traceLength_ms ) {
    log->line("Error in user parameter 'overlap': Inconsistent overlap specified: %d. Overlap has to be given in [ms]. Current trace length: %dms",
              vars->overlap_ms, traceLength_ms);
    env->addError();
  }

  vars->hdrID_time_samp1_s    = hdef->headerIndex( HDR_TIME_SAMP1.name );
  vars->hdrID_time_samp1_us = hdef->headerIndex( HDR_TIME_SAMP1_US.name );

  vars->overlap_numSamples = (int)( (float)vars->overlap_ms / shdr->sampleInt +0.5 );

  vars->numSamples_in = shdr->numSamples;
  shdr->numSamples += 2*vars->overlap_numSamples;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_overlap_(
                       csTraceGather* traceGather,
                       int* port,
                       int* numTrcToKeep,
                       csExecPhaseEnv* env,
                       csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    delete vars; vars = NULL;
    return;
  }

  int numTracesIn = traceGather->numTraces();
  if( traceGather->numTraces() != 3 ) {
    return;
  }
  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    float* samples1 = traceGather->trace(0)->getTraceSamples();
    float* samples2 = traceGather->trace(1)->getTraceSamples();
    float* samples3 = traceGather->trace(2)->getTraceSamples();
    for( int isamp = 0; isamp < vars->numSamples_in; isamp++ ) {
      samples1[isamp+vars->overlap_numSamples] = samples1[isamp];
      samples2[isamp+vars->overlap_numSamples] = samples2[isamp];
      samples3[isamp+vars->overlap_numSamples] = samples3[isamp];
    }
    for( int isamp = 0; isamp < vars->overlap_numSamples; isamp++ ) {
      samples1[isamp] = 0.0;
      samples2[isamp] = samples1[isamp+vars->overlap_numSamples];
      samples3[isamp] = samples2[isamp+vars->overlap_numSamples];
    }
    for( int isamp = 0; isamp < vars->overlap_numSamples; isamp++ ) {
      samples1[shdr->numSamples-isamp] = samples2[isamp];
      samples2[shdr->numSamples-isamp] = samples3[isamp];
    }
    traceGather->moveTraceTo( 0, traceGather, 3 );  // Move first trace to last position. --> Only last trace in gather is output
    *numTrcToKeep = 2;
  }
  else if( numTracesIn != 0 ) {
    //    float* samples1 = traceGather->trace(0)->getTraceSamples();
    float* samples2 = traceGather->trace(1)->getTraceSamples();
    float* samples3 = traceGather->trace(2)->getTraceSamples();
    for( int isamp = 0; isamp < vars->numSamples_in; isamp++ ) {
      samples3[isamp+vars->overlap_numSamples] = samples3[isamp];
    }
    for( int isamp = 0; isamp < vars->overlap_numSamples; isamp++ ) {
      samples3[isamp] = samples2[isamp+vars->overlap_numSamples];
    }
    for( int isamp = 0; isamp < vars->overlap_numSamples; isamp++ ) {
      samples2[shdr->numSamples-isamp] = samples3[isamp];
    }
    traceGather->moveTraceTo( 0, traceGather, 3 );  // Move first trace to last position. --> Only last trace in gather is output
    *numTrcToKeep = 2;
  }
  else {
    *numTrcToKeep = 0;
  }

  if( edef->isDebug() ) {
    for( int itrc = 0; itrc < traceGather->numTraces(); itrc++ ) {
      csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
      int time_samp1_s  = trcHdr->intValue( vars->hdrID_time_samp1_s );
      int time_samp1_us = trcHdr->intValue( vars->hdrID_time_samp1_us );
      log->line("Trace #-5d,  Time:  %12ds   %10dus", itrc, time_samp1_s, time_samp1_us);
    }
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_overlap_( csParamDef* pdef ) {
  pdef->setModule( "OVERLAP", "Create data overlap between adjacent traces", "Duplicates data from adjacent traces and pads it at start and end of trace" );

  pdef->addParam( "overlap", "Size of overlap [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Size of overlap [ms]", "Data of length 'overlap'ms from adjacent traces is added to start and end of trace" );

  pdef->addParam( "absolute_time", "Acknowledge absolute time?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Absolute time is acknowledged. No overlap will be created from adjacent traces that do not have adjacent absolute time stamps" );
  pdef->addOption( "no", "Absolute time is not acknowledged. Overlap will be created for all adjacent traces" );
}

extern "C" void _params_mod_overlap_( csParamDef* pdef ) {
  params_mod_overlap_( pdef );
}
extern "C" void _init_mod_overlap_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_overlap_( param, env, log );
}
extern "C" void _exec_mod_overlap_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_overlap_( traceGather, port, numTrcToKeep, env, log );
}


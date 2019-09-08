/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csGeolibUtils.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: RMS
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_rms {
  struct VariableStruct {
    int startSamp;
    int endSamp;
    int hdrId_rms;
    int hdrType_rms;
  };
}

using mod_rms::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_rms_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->startSamp   = 0;
  vars->endSamp     = 0;
  vars->hdrId_rms   = -1;
  vars->hdrType_rms = 0;

  //---------------------------------------------------------
  bool isTimeDomain = true;

  if( param->exists("domain") ) {
    std::string text;
    param->getString( "domain", &text );
    if( !text.compare( "sample" ) ) {
      isTimeDomain = false;
    }
    else if( !text.compare( "time" ) ) {
      isTimeDomain = true;
    }
    else {
      log->line("Domain option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }
  float startTime = 0.0;
  float endTime = 0.0;

  double start_in = 0.0;
  double end_in   = 0.0;

  if( param->exists("start") ) {
    param->getDouble( "start", &start_in );
  }
  if( param->exists("end") ) {
    param->getDouble( "end", &end_in );
  }

  if( isTimeDomain ) {
    startTime = start_in;
    endTime   = end_in;
    if( startTime < 0 ) {
      log->warning("Start time (%f) needs to be greater or equal to 0.0.", startTime);
      startTime = 0;
    }
    if( endTime > (shdr->numSamples-1)*shdr->sampleInt ) {
      log->warning("End time (%f) exceeds length of trace (%f).", endTime, (shdr->numSamples-1)*shdr->sampleInt );
      endTime = (shdr->numSamples-1)*shdr->sampleInt;
    }
    else if( endTime == 0 ) {
      endTime = (shdr->numSamples-1)*shdr->sampleInt;
    }
    if( startTime > endTime ) log->error("Start time (%f) needs to be smaller than end time (%f).", startTime, endTime);

    vars->startSamp = (int)(startTime / shdr->sampleInt);  // All in milliseconds
    vars->endSamp   = (int)(endTime / shdr->sampleInt);
  }
  else {
    //
    // NOTE: User input is '1' for first sample. Internally, '0' is used!!
    //
    vars->startSamp = (int)start_in;
    vars->endSamp   = (int)end_in;
    if( vars->startSamp < 1 ) log->error("Start sample (%d) needs to be greater or equal to 1.", vars->startSamp);
    if( vars->startSamp > vars->endSamp ) log->error("Start sample (%d) needs to be smaller than end sample (%d).", vars->startSamp, vars->endSamp);
    if( vars->endSamp > shdr->numSamples ) log->error("End sample (%d) exceeds number of samples (%d).", vars->endSamp, shdr->numSamples );
    else if( vars->endSamp == 0 ) {
      vars->endSamp = shdr->numSamples;
    }
    vars->startSamp -= 1;   // see note above..
    vars->endSamp   -= 1;
    startTime = (float)vars->startSamp * shdr->sampleInt;
    endTime   = (float)vars->endSamp * shdr->sampleInt;
  }

  //---------------------------------------------
  //
  std::string headerName_rms("rms");
  if( param->exists("hdr_rms") ) {
    param->getString("hdr_rms", &headerName_rms, 0);
  }
  if( !hdef->headerExists( headerName_rms ) ) {
    hdef->addHeader( TYPE_FLOAT, headerName_rms, "RMS value" );
  }

  vars->hdrId_rms = hdef->headerIndex(headerName_rms);
  vars->hdrType_rms = hdef->headerType(headerName_rms);
  if( vars->hdrType_rms != TYPE_FLOAT && vars->hdrType_rms != TYPE_DOUBLE ) {
    log->error("Trace header '%s' exists but has wrong type (%s). Type should be 'float' or 'double'.",
               cseis_geolib::csGeolibUtils::typeText(vars->hdrType_rms), headerName_rms.c_str());
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_rms_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    delete vars; vars = NULL;
    return true;
  }

  int nSamples = vars->endSamp - vars->startSamp + 1;
  float rms = compute_rms( &(trace->getTraceSamples()[vars->startSamp]), nSamples );

  // Set double value to account for different header types
  trace->getTraceHeader()->setDoubleValue( vars->hdrId_rms, rms );

  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_rms_( csParamDef* pdef ) {
  pdef->setModule( "RMS", "Compute RMS value in given time window", "Computed RMS value is stored in trace header 'rms'" );

  pdef->addParam( "domain", "Time or sample domain", NUM_VALUES_FIXED );
  pdef->addValue( "time", VALTYPE_OPTION );
  pdef->addOption( "time", "Window is specified in time [ms] (or frequency [Hz])" );
  pdef->addOption( "sample", "Window is specified in samples (1 for first sample)" );

  pdef->addParam( "start", "Start time/sample", NUM_VALUES_FIXED, "Start time or sample, this depends on the 'domain' setting" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Start time/sample of cross-correlation window" );

  pdef->addParam( "end", "End time/sample", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "End time/sample of cross-correlation window", "= 0 : Until end of trace" );

  pdef->addParam( "hdr_rms", "Name of trace header where RMS value shall be stored", NUM_VALUES_FIXED, "...may be new or existing trace header" );
  pdef->addValue( "rms", VALTYPE_STRING );
}

extern "C" void _params_mod_rms_( csParamDef* pdef ) {
  params_mod_rms_( pdef );
}
extern "C" void _init_mod_rms_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_rms_( param, env, log );
}
extern "C" bool _exec_mod_rms_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_rms_( trace, port, env, log );
}


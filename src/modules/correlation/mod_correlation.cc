/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "geolib_math.h"
#include "geolib_string_utils.h"
#include "geolib_methods.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: CORRELATION
 * 
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_correlation {
  struct VariableStruct {
    int hdrId_start;
    int hdrId_end;
    cseis_geolib::type_t hdrType_start;
    cseis_geolib::type_t hdrType_end;
    int hdrId_cross_lag;
    int hdrId_cross_amp;
    int startSamp; // Start sample index
    int endSamp;   // End sample index
    int maxLag_samples;    // Maximum lag time in samples
    float* buffer;
    int mode;
    int numSamplesOrig;
    int numSamplesBuffer;
    bool isTimeDomain;
    bool normalise;
    bool dampen;
    bool isCrossStacked;
  };
  static int const MODE_CROSS = 1;
  static int const MODE_AUTO  = 2;
  static int const MODE_AUTO_TWOSIDED  = 3;
  static int const MODE_CROSS_PILOT    = 4;
  static int const MODE_CROSS_PILOT_INPUT = 5;
}

//*************************************************************************************************
// Init phase
//
//
//*************************************************************************************************
void init_mod_correlation_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  mod_correlation::VariableStruct* vars = new mod_correlation::VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_MULTITRACE );

  vars->hdrId_cross_lag = -1;
  vars->hdrId_cross_amp = -1;
  vars->buffer          = NULL;
  vars->mode            = mod_correlation::MODE_CROSS;
  vars->startSamp       = 0;
  vars->endSamp         = shdr->numSamples-1;
  vars->maxLag_samples  = 0;
  vars->numSamplesBuffer = 0;

  vars->hdrId_start     = -1;
  vars->hdrId_end       = -1;
  vars->hdrType_start   = TYPE_UNKNOWN;
  vars->hdrType_end     = TYPE_UNKNOWN;
  vars->isTimeDomain    = true;
  vars->normalise = false;
  vars->dampen = false;
  vars->isCrossStacked = false;

  //---------------------------------------------------------
  // Set mode
  //
  std::string text;

  if( param->exists("mode") ) {
    param->getString( "mode", &text );
    if( !text.compare( "cross" ) ) {
      vars->mode = mod_correlation::MODE_CROSS;
    }
    else if( !text.compare( "auto" ) ) {
      vars->mode = mod_correlation::MODE_AUTO;
    }
    else if( !text.compare( "auto_twosided" ) ) {
      vars->mode = mod_correlation::MODE_AUTO_TWOSIDED;
    }
    else if( !text.compare( "cross_stacked" ) ) {
      vars->mode = mod_correlation::MODE_CROSS;
      vars->isCrossStacked = true;
    }
    else if( !text.compare( "cross_pilot" ) ) {
      vars->mode = mod_correlation::MODE_CROSS_PILOT;
    }
    else if( !text.compare( "cross_pilot_input" ) ) {
      vars->mode = mod_correlation::MODE_CROSS_PILOT_INPUT;
    }
    else {
      log->error("Unknown option: '%s'", text.c_str());
    }
  }
  if( vars->mode == mod_correlation::MODE_CROSS ) {
    edef->setTraceSelectionMode( TRCMODE_FIXED, 2 );
  }
  else if( vars->mode == mod_correlation::MODE_CROSS_PILOT || vars->mode == mod_correlation::MODE_CROSS_PILOT_INPUT ) {
    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }
  else {
    edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }
  
  if( param->exists("norm") ) {
    param->getString( "norm", &text );
    if( !text.compare( "yes" ) ) {
      vars->normalise = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->normalise = false;
    }
    else {
      log->error("Unknown option: '%s'", text.c_str());
    }
  }
  if( param->exists("damping") ) {
    param->getString( "damping", &text );
    if( !text.compare( "yes" ) ) {
      vars->dampen = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->dampen = false;
    }
    else {
      log->error("Unknown option: '%s'", text.c_str());
    }
  }

  //---------------------------------------------------------
  // Set time window
  //
  float startTime = 0.0;
  float endTime   = shdr->sampleInt * (float)vars->endSamp;

  vars->isTimeDomain = true;
  if( param->exists( "domain" ) ) {
    param->getString( "domain", &text );
    if( !text.compare( "sample" ) ) {
      vars->isTimeDomain = false;
    }
    else if( !text.compare( "time" ) ) {
      vars->isTimeDomain = true;
    }
    else {
      log->line("Domain not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }

  if( param->exists("start") ) {
    string textStart;
    param->getString( "start", &textStart );
    if( isDigit( textStart.at(0) ) ) {
      float start = (float)atof( textStart.c_str() );
      if( vars->isTimeDomain ) {
        startTime = start;
        vars->startSamp = (int)round( startTime / shdr->sampleInt );
      }
      else {
        vars->startSamp = (int)round( start ) - 1;  // User provides sample index starting at 1
        startTime = (float)vars->startSamp * shdr->sampleInt;
      }
    }
    else {
      vars->hdrId_start   = hdef->headerIndex(textStart);
      vars->hdrType_start = hdef->headerType(textStart);
    }
  }
  if( param->exists("end") ) {
    string textEnd;
    param->getString( "end", &textEnd );
    if( isDigit( textEnd.at(0) ) ) {
      float end = (float)atof( textEnd.c_str() );
      if( vars->isTimeDomain ) {
        endTime = end;
        vars->endSamp = (int)round( endTime / shdr->sampleInt );
      }
      else {
        vars->endSamp = (int)end - 1;  // User provides sample index starting at 1
        endTime = (float)vars->endSamp * shdr->sampleInt;
      }
    }
    else {
      vars->hdrId_end   = hdef->headerIndex(textEnd);
      vars->hdrType_end = hdef->headerType(textEnd);
    }
  }

  if( vars->isTimeDomain ) {
    if( startTime < 0.0 ) log->error("Start time (%f) needs to be greater or equal to 0.0.", startTime);
    if( startTime > endTime ) log->error("Start time (%f) needs to be smaller than end time (%f).", startTime, endTime);
    if( endTime >= (float)shdr->numSamples*shdr->sampleInt ) log->error("Specified end time (%fms) is greater than trace length (%fms).", endTime, (float)shdr->numSamples*shdr->sampleInt);
  }
  else {
    // NOTE: User input is '1' for first sample. Internally, '0' is used
    if( vars->startSamp < 0 ) log->error("Start sample (%d) needs to be greater or equal to 1.", vars->startSamp+1);
    if( vars->startSamp > vars->endSamp ) log->error("Start sample (%d) needs to be smaller than end sample (%d).", vars->startSamp+1, vars->endSamp+1);
    if( vars->endSamp > shdr->numSamples-1 ) log->error("End sample/time greater than number samples in trace (%d).", vars->endSamp+1, shdr->numSamples);
  }

  //---------------------------------------------------------
  // Set maximum lag
  //
  int nSampIn = vars->endSamp - vars->startSamp + 1;
  if( param->exists("max_lag") ) {
    float maxLag_ms;
    param->getFloat( "max_lag", &maxLag_ms );
    vars->maxLag_samples = (int)round( maxLag_ms / shdr->sampleInt );
    if( vars->maxLag_samples > nSampIn-1 ) {
      log->warning("Specified maximum lag time (%fms) exceeds number of samples (%d) in selected time window (%fms to %fms). Set maximum lag to %ms",
                   maxLag_ms, nSampIn, startTime, endTime, (float)((nSampIn-1)*shdr->sampleInt) );
      vars->maxLag_samples = nSampIn-1;
    }
  }
  else {
    vars->maxLag_samples = nSampIn - 1;
  }
  //---------------------------------------------------------
  // Other settings
  //
  vars->numSamplesOrig   = shdr->numSamples;
  vars->numSamplesBuffer = shdr->numSamples;
  if( vars->mode == mod_correlation::MODE_CROSS || vars->mode == mod_correlation::MODE_AUTO_TWOSIDED || vars->mode == mod_correlation::MODE_CROSS_PILOT ) {
    vars->numSamplesBuffer = 2*vars->maxLag_samples + 1;;
    shdr->numSamples       = vars->numSamplesBuffer;
  }
  else if( vars->mode == mod_correlation::MODE_AUTO ) {
    vars->numSamplesBuffer = vars->maxLag_samples + 1;
    shdr->numSamples       = vars->numSamplesBuffer;
  }
  else if( vars->mode == mod_correlation::MODE_CROSS_PILOT_INPUT ) {
    vars->numSamplesBuffer = 2*vars->maxLag_samples + 1;;
  }

  vars->buffer = new float[vars->numSamplesBuffer];

  if( !hdef->headerExists("cross_lag") ) {
    hdef->addHeader( TYPE_FLOAT, "cross_lag", "Cross-correlation lag time" );
  }
  vars->hdrId_cross_lag = hdef->headerIndex("cross_lag");
  type_t type = hdef->headerType("cross_lag");
  if( type != TYPE_FLOAT && type != TYPE_DOUBLE ) {
    log->error("Header 'cross_lag' exists but has wrong format. Must be FLOAT or DOUBLE.");
  }

  if( !hdef->headerExists("cross_amp") ) {
    hdef->addHeader( TYPE_FLOAT, "cross_amp", "Cross-correlation maximum amplitude" );
  }
  vars->hdrId_cross_amp = hdef->headerIndex("cross_amp");
  type = hdef->headerType("cross_amp");
  if( type != TYPE_FLOAT && type != TYPE_DOUBLE ) {
    log->error("Header 'cross_amp' exists but has wrong format. Must be FLOAT or DOUBLE.");
  }

  if( edef->isDebug() ) {
    log->line("time1: %f, time2: %f, sample1: %d, sample2: %d, sampInt: %f, max lag: %d samples\n", startTime, endTime, vars->startSamp, vars->endSamp, shdr->sampleInt, vars->maxLag_samples );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_correlation_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* log )
{
  mod_correlation::VariableStruct* vars = reinterpret_cast<mod_correlation::VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    if( vars->buffer != NULL ) {
      delete [] vars->buffer;
      vars->buffer = NULL;
    }
    delete vars; vars = NULL;
    return;
  }

  if( vars->mode == mod_correlation::MODE_CROSS && traceGather->numTraces()!= 2 ) {
    log->warning("Wrong number of input traces for cross-correlation. Expected: 2, found: %d", traceGather->numTraces() );
    traceGather->freeAllTraces();
    return;
  }

  float time_maxAmp = 0;
  float maxAmp = 0.0;

  int nSampIn = vars->endSamp - vars->startSamp + 1;

  //---------------------------------------------------
  // Read in window start/end times from trace header, if specified
  //
  if( vars->hdrId_start >= 0 ) {
    double start = 0.0;
    if( vars->hdrType_start == TYPE_INT ) {
      start = (double)traceGather->trace(0)->getTraceHeader()->intValue( vars->hdrId_start );
    }
    else {
      start = traceGather->trace(0)->getTraceHeader()->doubleValue( vars->hdrId_start );
    }
    if( vars->isTimeDomain ) {
      vars->startSamp = (int)round( start / shdr->sampleInt );
    }
    else {
      vars->startSamp = (int)round( start ) - 1;   // -1 because user specifies sample number starting at 1, not 0
    }
    if( vars->startSamp < 0 ) vars->startSamp = 0;
  }
  if( vars->hdrId_end >= 0 ) {
    double end = 0.0;
    if( vars->hdrType_end == TYPE_INT ) {
      end = (double)traceGather->trace(0)->getTraceHeader()->intValue( vars->hdrId_end );
    }
    else {
      end = traceGather->trace(0)->getTraceHeader()->doubleValue( vars->hdrId_end );
    }
    if( vars->isTimeDomain ) {
      vars->endSamp = (int)( end / shdr->sampleInt + 0.5 );
    }
    else {
      vars->endSamp = (int)( end + 0.5 ) - 1;   // -1 because user specifies sample number starting at 1, not 0
    }
    if( vars->endSamp > vars->numSamplesOrig-1 ) vars->endSamp = vars->numSamplesOrig - 1;
    if( vars->endSamp <= vars->startSamp )   vars->endSamp = vars->startSamp + vars->maxLag_samples;
  }
  // Check start/end sample versus max lag
  if( vars->hdrId_start >= 0 || vars->hdrId_end >= 0 ) {
    nSampIn = vars->endSamp - vars->startSamp + 1;
    if( vars->maxLag_samples > nSampIn-1 ) {
      log->warning("Specified maximum lag time (%fms) exceeds number of samples (%d) specified in input trace headers (time window: %fms to %fms)",
                   (float)vars->maxLag_samples * shdr->sampleInt, nSampIn, (float)vars->startSamp*shdr->sampleInt, (float)vars->endSamp*shdr->sampleInt );
      vars->endSamp = vars->startSamp + vars->maxLag_samples;
      if( vars->endSamp > vars->numSamplesOrig-1 ) {
        vars->endSamp = vars->numSamplesOrig-1;
        if( vars->endSamp-vars->startSamp < vars->maxLag_samples ) {
          vars->startSamp = vars->endSamp-vars->maxLag_samples;
          if( vars->startSamp < 0 ) {
            log->error("Inconsistent start/end samples found in input trace headers. Bailing out...");
          }
        }
      }
    }
  }

  //---------------------------------------------------
  // Cross-correlation
  //
  if( ( vars->mode == mod_correlation::MODE_CROSS ) ) {
    if( edef->isDebug() ) log->line("Cross-correlation input nSamples: %d, output nSamples: %d, start: %d, end: %d", nSampIn, vars->numSamplesBuffer, vars->startSamp, vars->endSamp );
    float* trace1Ptr = traceGather->trace(0)->getTraceSamples();
    float* trace2Ptr = traceGather->trace(1)->getTraceSamples();
    compute_twosided_correlation(
                          &trace1Ptr[vars->startSamp],
                          &trace2Ptr[vars->startSamp],
                          nSampIn,
                          vars->buffer,
                          vars->maxLag_samples,
			  vars->dampen );

    int sampleIndex_maxAmp = 0;
    if( vars->isCrossStacked ) {
      int numHalf = shdr->numSamples/2;
      for( int isamp = 0; isamp <= numHalf; isamp++ ) {
        trace1Ptr[numHalf-isamp] = 0.5f * ( vars->buffer[numHalf+isamp] + vars->buffer[numHalf-isamp] );
        trace1Ptr[numHalf+isamp] = trace1Ptr[numHalf-isamp];
        if( trace1Ptr[numHalf-isamp] > maxAmp ) {
          maxAmp = trace1Ptr[isamp];
          sampleIndex_maxAmp = isamp;
        }
      }
      traceGather->freeTrace(1);
    }
    else {
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        trace1Ptr[isamp] = vars->buffer[isamp];
        if( trace1Ptr[isamp] > maxAmp ) {
          maxAmp = trace1Ptr[isamp];
          sampleIndex_maxAmp = isamp;
        }
      }
    }

    float sampleIndex_maxAmp_float = getQuadMaxSample( trace1Ptr, sampleIndex_maxAmp, shdr->numSamples, &maxAmp );
    time_maxAmp = (sampleIndex_maxAmp_float-(float)vars->maxLag_samples)*shdr->sampleInt;
  }
  //---------------------------------------------------
  // Cross-correlation
  //
  else if( vars->mode == mod_correlation::MODE_CROSS_PILOT || vars->mode == mod_correlation::MODE_CROSS_PILOT_INPUT ) {
    float* trace1Ptr = traceGather->trace(0)->getTraceSamples();
    for( int itrc = 1; itrc < traceGather->numTraces(); itrc++ ) {
      float* trace2Ptr = traceGather->trace(itrc)->getTraceSamples();
      compute_twosided_correlation(
				   &trace1Ptr[vars->startSamp],
				   &trace2Ptr[vars->startSamp],
				   nSampIn,
				   vars->buffer,
				   vars->maxLag_samples,
				   vars->dampen );
      
      int sampleIndex_maxAmp = 0;
      maxAmp = 0.0;
      if( vars->mode == mod_correlation::MODE_CROSS_PILOT ) {
        memcpy( trace2Ptr, vars->buffer, vars->numSamplesBuffer*sizeof(float) );
      }
      for( int isamp = 0; isamp < vars->numSamplesBuffer; isamp++ ) {
	if( vars->buffer[isamp] > maxAmp ) {
	  maxAmp = vars->buffer[isamp];
	  sampleIndex_maxAmp = isamp;
	}
      }
      
      float sampleIndex_maxAmp_float = getQuadMaxSample( vars->buffer, sampleIndex_maxAmp, vars->numSamplesBuffer, &maxAmp );
      time_maxAmp = (sampleIndex_maxAmp_float-(float)vars->maxLag_samples)*shdr->sampleInt;
      csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
      trcHdr->setDoubleValue( vars->hdrId_cross_lag, time_maxAmp );
      trcHdr->setDoubleValue( vars->hdrId_cross_amp, maxAmp );
    }
    // Free pilot trace
    if( vars->mode == mod_correlation::MODE_CROSS_PILOT ) {
      traceGather->freeTrace( 0 );
    }
    else {
      traceGather->trace(0)->getTraceHeader()->setDoubleValue( vars->hdrId_cross_lag, 0 );
      traceGather->trace(0)->getTraceHeader()->setDoubleValue( vars->hdrId_cross_amp, 1 );
    }
  }
  //---------------------------------------------------
  // Autocorrelation
  //
  else {  // MODE_AUTO
    float* tracePtr = traceGather->trace(0)->getTraceSamples();
    float* bufferPtr;
    int sampleIndex_zeroLag = 0;
    if( vars->mode == mod_correlation::MODE_AUTO ) {
      bufferPtr = &vars->buffer[0];
    }
    else { // two-sided auto-correlation
      sampleIndex_zeroLag = (shdr->numSamples+1)/2 - 1;
      bufferPtr = &vars->buffer[sampleIndex_zeroLag];
    }
    compute_onesided_auto_correlation(
                                      &tracePtr[vars->startSamp],
                                      nSampIn,
                                      bufferPtr,
                                      vars->maxLag_samples,
				      vars->dampen );
    if( vars->normalise ) {
      for( int ilag = vars->maxLag_samples; ilag >= 0; ilag-- ) {
        bufferPtr[ilag] /= bufferPtr[0];
      }
    }
    if( vars->mode == mod_correlation::MODE_AUTO_TWOSIDED ) {
      for( int isamp = sampleIndex_zeroLag-1; isamp >= 0; isamp-- ) {
        vars->buffer[isamp] = vars->buffer[sampleIndex_zeroLag-isamp+sampleIndex_zeroLag];
      }
    }
    memcpy( tracePtr, vars->buffer, shdr->numSamples*sizeof(float) );
    maxAmp      = vars->buffer[sampleIndex_zeroLag];
    time_maxAmp = (float)sampleIndex_zeroLag * shdr->sampleInt;
  }

  if( vars->mode != mod_correlation::MODE_CROSS_PILOT && vars->mode != mod_correlation::MODE_CROSS_PILOT_INPUT ) {
    csTraceHeader* trcHdr = traceGather->trace(0)->getTraceHeader();
    //  trcHdr->setDoubleValue( vars->hdrId_cross_lag, (sampleIndex_maxAmp_float-nSamples+1)*shdr->sampleInt );
    trcHdr->setDoubleValue( vars->hdrId_cross_lag, time_maxAmp );
    trcHdr->setDoubleValue( vars->hdrId_cross_amp, maxAmp );
    
    if( vars->mode == mod_correlation::MODE_CROSS ) {
      // Free second trace, keep first trace
      traceGather->freeTrace( 1 );
    }
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_correlation_( csParamDef* pdef ) {
  pdef->setModule( "CORRELATION", "Cross-correlation between adjacent traces, or auto-correlation of same trace" );
  pdef->addDoc("Note that this version does not provide not a high-res correlation. Correlation lags are only accurate to ~1/2 the sample interval, not less");
  pdef->addDoc("The maximum cross-correlation time and amplitude are stored in new trace headers 'cross_lag' and 'cross_amp'");

  pdef->addParam( "mode", "Mode of correlation: Cross- or auto-correlation", NUM_VALUES_FIXED );
  pdef->addValue( "cross", VALTYPE_OPTION );
  pdef->addOption( "cross", "Cross-correlation of each consecutive trace pair" );
  pdef->addOption( "auto", "Auto-correlation of each single trace" );
  pdef->addOption( "auto_twosided", "Auto-correlation of each single trace", "Output two-sided auto-correlation" );
  pdef->addOption( "cross_stacked", "Cross-correlation of each consecutive trace pair", "Stack positive & negative side" );
  pdef->addOption( "cross_pilot", "Cross-correlation of each trace with first trace in ensemble", "Input one ensemble. First trace is treated as pilot trace. Output all traces except the pilot trace" );
  pdef->addOption( "cross_pilot_input", "Same as 'cross_pilot', but instead of outputting cross-correlation function, output input data with updated trace headers: cross_lag & cross_amp" );

  pdef->addParam( "domain", "Is correlation window given in time or in samples?", NUM_VALUES_FIXED );
  pdef->addValue( "time", VALTYPE_OPTION );
  pdef->addOption( "time", "Window is specified in time [ms]" );
  pdef->addOption( "sample", "Window is specified in samples (1 for first sample)" );

  pdef->addParam( "start", "Start time/sample", NUM_VALUES_FIXED, "Start time or sample, this depends on the 'domain' setting" );
  pdef->addValue( "0", VALTYPE_STRING, "Start time/sample of cross-correlation window", "Alternatively, name of trace header containing start time/sample" );

  pdef->addParam( "end", "End time/sample", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "End time/sample of cross-correlation window", "Alternatively, name of trace header containing end time/sample" );

  pdef->addParam( "max_lag", "Maximum cross-correlation lag [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum cross-correlation lag [ms]" );

  pdef->addParam( "norm", "Normalise output correlation function?", NUM_VALUES_FIXED, "Only applicable for auto-correlation" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Normalise zero lag amplitude to 1" );
  pdef->addOption( "no", "Do not normalise output" );

  pdef->addParam( "damping", "Dampen correlation function edges?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Apply damping to correlation edges" );
  pdef->addOption( "no", "Do not apply damping" );
}

extern "C" void _params_mod_correlation_( csParamDef* pdef ) {
  params_mod_correlation_( pdef );
}
extern "C" void _init_mod_correlation_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_correlation_( param, env, log );
}
extern "C" void _exec_mod_correlation_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_correlation_( traceGather, port, numTrcToKeep, env, log );
}


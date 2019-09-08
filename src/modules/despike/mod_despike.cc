/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csDespike.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: despike
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_despike {
  struct VariableStruct {
    int method;
    cseis_geolib::csDespike* despike;
    float* buffer;
    int output;
    int hdrID_num_spikes;
    int hdrID_num_spikes_samples;
  };
  static int const OUTPUT_APPLY = 1;
  static int const OUTPUT_DIFF  = 2;
}
using namespace mod_despike;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_despike_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->method  = csDespike::COSINE_TAPER;
  vars->buffer  = NULL;
  vars->despike = NULL;
  vars->output  = OUTPUT_APPLY;
  vars->hdrID_num_spikes = -1;
  vars->hdrID_num_spikes_samples = -1;

  std::string text;

  if( param->exists( "method" ) ) {
    param->getString( "method", &text );
    text = toLowerCase( text );
    if( !text.compare( "cos_taper" ) ) {
      vars->method = csDespike::COSINE_TAPER;
    }
    else if( !text.compare( "interpolation" ) ) {
      vars->method = csDespike::LINEAR_INTERPOLATION;
    }
    else if( !text.compare( "zero" ) ) {
      vars->method = csDespike::SET_TO_ZERO;
    }
    else {
      log->line("Unknown argument for user parameter 'method': '%s'.", text.c_str());
      env->addError();
    }
  }

  //---------------------------------------------------------
  if( param->exists("output") ) {
    param->getString( "output", &text );
    if( !text.compare( "apply" ) ) {
      vars->output = OUTPUT_APPLY;
    }
    else if( !text.compare( "diff" ) ) {
      vars->output = OUTPUT_DIFF;
    }
    else {
      log->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }

  //---------------------------------------------------------
  DespikeConfig config;
  if( shdr->domain == DOMAIN_XT || shdr->domain == DOMAIN_XD ) {
    csDespike::getDefaultTimeNoiseBurstConfig( config );
  }
  else {
    csDespike::getDefaultFrequencySpikeConfig( config );
  }

  int numValues = param->getNumValues("win_ref");
  if( numValues > 0 ) {
    param->getFloat( "win_ref", &config.widthRefWin, 0 );
    if( numValues > 1 ) {
      param->getFloat( "win_ref", &config.incWin, 1 );
    }
  }

  if( param->exists("win_spike") ) {
    param->getFloat( "win_spike", &config.minWidthSpikeWin, 0 );
  }

  if( param->exists("start") ) {
    param->getDouble( "start", &config.start );
  }
  if( param->exists("end") ) {
    param->getDouble( "end", &config.stop );
  }

  param->getFloat( "max_ratio", &config.maxRatio );
  //---------------------------------------------------------


  float traceLength = shdr->numSamples*shdr->sampleInt;
  if( config.widthRefWin < config.minWidthSpikeWin ) {
    log->warning("The despike window is chosen larger than the reference window: %f < %f. This is not good practise.", config.widthRefWin, config.minWidthSpikeWin );
  }
  if( config.widthRefWin >= traceLength ) {
    log->error("The reference window is chosen larger or equal to the current trace length: %f >= %f.", config.widthRefWin, traceLength );
  }
  if( (int)(config.minWidthSpikeWin/shdr->sampleInt) < 1 ) {
    log->error("The despike window is chosen too small to be effective: %f.", config.minWidthSpikeWin );
  }
  if( config.incWin > config.minWidthSpikeWin ) {
    log->error("The increment of the despike window is too large: %f. Maximum increment length: %f", config.incWin, config.minWidthSpikeWin );
  }
  if( (config.start >= config.stop || config.start >= traceLength) && (config.stop > 0 || config.start > config.stop)) {
    log->error("Inconsistent start/end points: %f %f", config.start, config.stop );
  }
  if( config.maxRatio <= 1.0 ) {
    log->error("Incorrect maximum ratio: %f. The threshold ratio must be greater than 1 to be effective.", config.maxRatio );
  }

  if( !hdef->headerExists("num_spikes") ) {
    hdef->addHeader( TYPE_INT, "num_spikes" );
  }
  else if( hdef->headerType("num_spikes") != TYPE_INT ) {
    log->error("Trace header 'num_spikes' exists but has wrong type: Should be INT");
  }
  if( !hdef->headerExists("nsamp_spikes") ) {
    hdef->addHeader( TYPE_INT, "nsamp_spikes" );
  }
  else if( hdef->headerType("nsamp_spikes") != TYPE_INT ) {
    log->error("Trace header 'nsamp_spikes' exists but has wrong type: Should be INT");
  }
  vars->hdrID_num_spikes = hdef->headerIndex("num_spikes");
  vars->hdrID_num_spikes_samples = hdef->headerIndex("nsamp_spikes");

  vars->despike = new csDespike( shdr->numSamples, shdr->sampleInt, config );
  vars->buffer = new float[shdr->numSamples];
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_despike_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    if( vars->despike != NULL ) {
      delete vars->despike;
      vars->despike = NULL;
    }
    if( vars->buffer != NULL ) {
      delete [] vars->buffer;
      vars->buffer = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  float* samples = trace->getTraceSamples();
  for( int i = 0; i < shdr->numSamples; i++ ) {
    vars->buffer[i] = samples[i];
  }

  int numSpikes = 0;
  int numSpikesSamples = 0;
  vars->despike->apply( vars->buffer, shdr->numSamples, numSpikes, numSpikesSamples );

  trace->getTraceHeader()->setIntValue( vars->hdrID_num_spikes, numSpikes );
  trace->getTraceHeader()->setIntValue( vars->hdrID_num_spikes_samples, numSpikesSamples );

  if( vars->output == OUTPUT_APPLY ) {
    if( vars->despike->method() == csDespike::COSINE_TAPER ) {
      for( int i = 0; i < shdr->numSamples; i++ ) {
        samples[i] *= vars->buffer[i];
      }
    }
    else {
      for( int i = 0; i < shdr->numSamples; i++ ) {
        samples[i] = vars->buffer[i];
      }
    }
  }
  else {  // OUTPUT_DIFF
    if( vars->despike->method() == csDespike::COSINE_TAPER ) {
      for( int i = 0; i < shdr->numSamples; i++ ) {
        samples[i] *= ( 1.0 - vars->buffer[i] );
      }
    }
    else {
      for( int i = 0; i < shdr->numSamples; i++ ) {
        samples[i] -= vars->buffer[i];
      }
    }
  }


  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_despike_( csParamDef* pdef ) {
  pdef->setModule( "DESPIKE", "Spike/noise burst removal" );

  pdef->addParam( "method", "Despike method. Specifies how to deal with the data window of the identified noise burst/spike", NUM_VALUES_FIXED );
  pdef->addValue( "cos_taper", VALTYPE_OPTION );
  pdef->addOption( "cos_taper", "Apply cosine taper over spike window");
  pdef->addOption( "interpolation", "Interpolate linearly between first and last sample of identified window" );
  pdef->addOption( "zero", "Set samples in spike window to zero" );

  pdef->addParam( "win_ref", "Reference window", NUM_VALUES_VARIABLE, "Width in units of trace, e.g. [ms] or [Hz]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Width of reference window in units of trace.", "Reference window over which median background value is computed" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Window increment in units of trace.", "Set to zero for sliding window, e.g. one sample interval" );

  pdef->addParam( "win_spike", "Width of despike window", NUM_VALUES_VARIABLE, "Width in units of trace, e.g. [ms] or [Hz]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Minimum width of despike window, in units of trace" );

  pdef->addParam( "start", "Start time/frequency for despike application", NUM_VALUES_FIXED,
                  "Despike operation will only be performed within the specified start/end window, given in the units of the trace , e.g. [ms] or [Hz]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Start of application window", "=0 : Use full trace from first sample" );

  pdef->addParam( "end", "End time/frequency for despike application", NUM_VALUES_FIXED,
                  "Despike operation will only be performed within the specified start/end window, given in the units of the trace , e.g. [ms] or [Hz]" );
  pdef->addValue( "0", VALTYPE_NUMBER, "End of application window", "=0 : Use full trace to last sample" );
  
  pdef->addParam( "max_ratio", "Maximum ratio for spike detection", NUM_VALUES_FIXED,
                  "Ratio is computed as follows: ratio[i] = sampleValue[i] / median( sampleValues[i +/- refWin/2] )" );
  pdef->addValue( "3", VALTYPE_NUMBER, "Max ratio for spike detection",
                  "Example: Ratio of 2 means that current sample value is twice the median value of the reference window" );

  pdef->addParam( "output", "Output options", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "apply", "Apply despike operation" );
  pdef->addOption( "diff", "Output detected spikes (difference)" );

}

extern "C" void _params_mod_despike_( csParamDef* pdef ) {
  params_mod_despike_( pdef );
}
extern "C" void _init_mod_despike_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_despike_( param, env, log );
}
extern "C" bool _exec_mod_despike_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_despike_( trace, port, env, log );
}


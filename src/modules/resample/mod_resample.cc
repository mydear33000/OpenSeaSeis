/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "geolib_methods.h"
#include "csFFTTools.h"
#include "csInterpolation.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: RESAMPLE
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_resample {
  struct VariableStruct {
    int   numSamplesOld;
    float sampleIntOld;  // [ms]
    float* buffer;
    float cutOffHz;
    bool debias;
    int filter;
    int hdrID_scalar;
    cseis_geolib::csInterpolation* interpol;
    int normOption;
    float normScalar;

    cseis_geolib::csFFTTools* fftTool;
  };
  static int const FILTER_NONE   = 0;
  static int const FILTER_FIR    = 11;
  static int const FILTER_BUTTER = 22;
  static int const NORM_YES   = 101;
  static int const NORM_NO    = 102;
  static int const NORM_RMS   = 103;
}
using namespace mod_resample;

void create_filter_coef( float cutOffFreqHz, int order, float* coef, int numCoef );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_resample_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->buffer = NULL;
  vars->debias = false;
  vars->filter = mod_resample::FILTER_NONE;
  vars->hdrID_scalar = -1;
  vars->interpol = NULL;
  vars->normOption = mod_resample::NORM_YES;
  vars->normScalar = 1.0f;

//---------------------------------------------
  vars->numSamplesOld = shdr->numSamples;
  vars->sampleIntOld  = shdr->sampleInt;

  float sampleIntNew;
  param->getFloat( "sample_int", &sampleIntNew );  // [ms]

  if( sampleIntNew <= 0 ) {
    log->error("Non-physical sample interval specified: %fms", sampleIntNew);
  }

  if( param->exists( "debias" ) ) {
    std::string text;
    param->getString( "debias", &text );
    text = toLowerCase( text );
    if( !text.compare( "yes" ) ) {
      vars->debias = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->debias = false;
    }
    else {
      log->line("Unknown argument for user parameter 'debias': '%s'.", text.c_str());
      env->addError();
    }
  }

  if( param->exists( "anti_alias" ) ) {
    std::string text;
    param->getString( "anti_alias", &text );
    text = toLowerCase( text );
    if( !text.compare( "yes" ) ) {
      vars->filter = mod_resample::FILTER_BUTTER;
    }
    else if( !text.compare( "fir" ) ) {
      vars->filter = mod_resample::FILTER_FIR;
    }
    else if( !text.compare( "no" ) ) {
      vars->filter = mod_resample::FILTER_NONE;
    }
    else {
      log->line("Unknown option: '%s'.", text.c_str());
      env->addError();
    }
  }

  float order = 20;
  if( param->exists("order") ) {
    param->getFloat("order", &order );
    if( order < 1 || order > 100 ) {
      log->error("Specified order = %.2f out of valid range = 1-100", order);
    }
    if( vars->filter == mod_resample::FILTER_NONE ) log->warning("Filter order specified but anti-alias filter is turned off");
  }
  if( param->exists("slope") ) {
    float slope;
    param->getFloat("slope", &slope );
    order = fabs(slope) / 6.0;
    if( vars->filter == mod_resample::FILTER_NONE ) log->warning("Filter slope specified but anti-alias filter is turned off");
  }
  float cutOffRatio = 0.8;
  if( param->exists( "cutoff" ) ) {
    param->getFloat( "cutoff", &cutOffRatio );
    if( cutOffRatio <= 0.0 || cutOffRatio >= 1.0 ) {
      log->error("Given cut-off ratio (=%f) outside of valid range (=0.0-1.0)", cutOffRatio);
    }
    if( vars->filter == mod_resample::FILTER_NONE ) log->warning("Filter cutoff specified but anti-alias filter is turned off");
  }

  bool isSinc = true;
  if( param->exists( "upsampling" ) ) {
    std::string text;
    param->getString( "upsampling", &text );
    if( !text.compare( "sinc" ) ) {
      isSinc = true;
    }
    else if( !text.compare( "quad" ) ) {
      isSinc = false;
    }
    else {
      log->line("Unknown option: '%s'.", text.c_str());
      env->addError();
    }
  }
  if( param->exists( "norm" ) ) {
    std::string text;
    param->getString( "norm", &text );
    if( !text.compare( "rms" ) ) {
      vars->normOption = mod_resample::NORM_RMS;
    }
    else if( !text.compare( "yes" ) ) {
      vars->normOption = mod_resample::NORM_YES;
    }
    else if( !text.compare( "no" ) ) {
      vars->normOption = mod_resample::NORM_NO;
    }
    else {
      log->line("Unknown option: '%s'.", text.c_str());
      env->addError();
    }
  }

  int numSamplesNew = shdr->numSamples;
  if( sampleIntNew < shdr->sampleInt ) {
    numSamplesNew = (int)( ((float)vars->numSamplesOld)*(shdr->sampleInt/sampleIntNew) );
    vars->buffer = new float[numSamplesNew];
    if( isSinc) vars->interpol = new csInterpolation( vars->numSamplesOld, vars->sampleIntOld, 8 );
  }
  else if( sampleIntNew > shdr->sampleInt ) {
    float freqNy      = 500.0/shdr->sampleInt;
    double ratio      = (int)((double)sampleIntNew / (double)shdr->sampleInt );
    int ratio_int = (int)round( ratio );
    int m;
    int twopm;
    csFFTTools::Powerof2( ratio_int, &m, &twopm );
    if( twopm != ratio_int ) {
      ratio = twopm * 2;
    }
    else {
      ratio = ratio_int;
    }
    if( fabs(sampleIntNew-shdr->sampleInt*ratio) > 0.001 ) {
      log->warning("Output sample interval was set to %f. Current implementation requires new_sample_int to be N*old_sample_int, where N=2^M", shdr->sampleInt*ratio);
    }
    sampleIntNew = shdr->sampleInt*ratio;
    float cutOffFreq  = ( freqNy / ratio ) * cutOffRatio;
    numSamplesNew     = (int)ceil((double)vars->numSamplesOld / ratio);  // Workaround BUGFIX 100504: Avoids clash for num samples = 2^N+1
    vars->fftTool     = new csFFTTools( vars->numSamplesOld, numSamplesNew, vars->sampleIntOld, sampleIntNew );
    vars->fftTool->setFilter( order, cutOffFreq, false );

    log->line("Sample int old/new: %f/%f ms\nNyquist: %f Hz\nCut-off frequency: %f Hz\nRatio: %f\nOrder: %f, #samples old/new/fft/fftout: %d/%d/%d/%d\n",
              vars->sampleIntOld, sampleIntNew, freqNy, cutOffFreq, ratio, order, vars->numSamplesOld, numSamplesNew, vars->fftTool->numFFTSamples(), vars->fftTool->numFFTSamplesOut() );
    if( edef->isDebug() ) {
      fprintf(stderr,"SampleInt: %f, Nyquist: %f, cutOff: %f, ratio: %f\n", shdr->sampleInt, freqNy, cutOffFreq, ratio );
    }
  }
  else {
    log->warning("Specified sample interval of %fms is the same as in the input data.\nData will stay unchanged.", sampleIntNew );
  }

  shdr->numSamples = numSamplesNew;
  shdr->sampleInt  = sampleIntNew;

  if( edef->isDebug() ) {
    float ratio = shdr->sampleInt/vars->sampleIntOld;
    int ratio_int = (int)round( ratio );
    log->line("nSamples (old,new): %d %d, sampleInt (old,new): %8.4f %8.4f, ratio: %f (%d)",
              vars->numSamplesOld, shdr->numSamples, vars->sampleIntOld, shdr->sampleInt, ratio, ratio_int);
    fprintf(stderr,"nSamples (old,new): %d %d, sampleInt (old,new): %8.4f %8.4f, ratio: %f (%d)\n",
            vars->numSamplesOld, shdr->numSamples, vars->sampleIntOld, shdr->sampleInt, ratio, ratio_int);
  }

  if( vars->normOption == mod_resample::NORM_RMS ) {
    if( !hdef->headerExists("resample_scalar") ) {
      hdef->addHeader( cseis_geolib::TYPE_FLOAT, "resample_scalar", "Inverse scalar applied during RESAMPLE" );
    }
    vars->hdrID_scalar = hdef->headerIndex("resample_scalar");
  }
  else if( vars->normOption == mod_resample::NORM_YES ) {
    float ratio = shdr->sampleInt / vars->sampleIntOld;
    if( ratio > 2.1f ) ratio /= 2.0f;   // ..no idea why this is necessary. Tested with 1ms and 2ms input data, resampled to 4ms, 8ms...
    vars->normScalar = 1.0/ratio;
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_resample_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    if( vars->buffer != NULL ) {
      delete [] vars->buffer;
      vars->buffer = NULL;
    }
    if( vars->fftTool != NULL ) {
      delete vars->fftTool;
      vars->fftTool = NULL;
    }
    if( vars->interpol != NULL ) {
      delete vars->interpol;
      vars->interpol = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  float* samples = trace->getTraceSamples();
  int numSamplesNew = shdr->numSamples;

  float resampleScalar = 1.0;
  if( shdr->sampleInt < vars->sampleIntOld ) {
    if( vars->interpol != NULL ) {
      for( int isamp = 0; isamp < numSamplesNew; isamp++ ) {
        float time_ms = (float)isamp * shdr->sampleInt;
        vars->buffer[isamp] = vars->interpol->valueAt( time_ms, samples );
      }
    }
    else {
      for( int isamp = 0; isamp < numSamplesNew; isamp++ ) {
        double sampleIndexOld = (double)(isamp*shdr->sampleInt) / vars->sampleIntOld;
        vars->buffer[isamp] = getQuadAmplitudeAtSample( samples, sampleIndexOld, vars->numSamplesOld );
      }
    }
    memcpy( samples, vars->buffer, numSamplesNew*sizeof(float) );
  }
  else if( shdr->sampleInt > vars->sampleIntOld ) {
    if( !vars->debias ) {
      resampleScalar = vars->fftTool->resample( samples, vars->filter != mod_resample::FILTER_NONE, vars->normOption == mod_resample::NORM_RMS );
    } // Remove DC bias first
    else {
      double sum = 0.0;
      int sampleCounter = 0;
      for( int isamp = 0; isamp < vars->numSamplesOld; isamp++ ) {
        if( samples[isamp] != 0 ) {
          sampleCounter += 1;
          sum += samples[isamp];
        }
      }
      float mean = sum / (double)sampleCounter;
      for( int isamp = 0; isamp < vars->numSamplesOld; isamp++ ) {
        if( samples[isamp] != 0 ) {
          samples[isamp] -= mean;
        }
      }

      resampleScalar = vars->fftTool->resample( samples, vars->filter != mod_resample::FILTER_NONE, vars->normOption == mod_resample::NORM_RMS );

      // Reapply DC bias afterwards
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        if( samples[isamp] != 0 ) {
          samples[isamp] += mean;
        }
      }
    }
    if( vars->normOption == mod_resample::NORM_YES ) {
      for( int isamp = 0; isamp < vars->numSamplesOld; isamp++ ) {
        samples[isamp] *= vars->normScalar;
      }      
    }
  }
  if( vars->normOption == mod_resample::NORM_RMS ) {
    if( resampleScalar > 0 ) resampleScalar = 1/resampleScalar;
    trace->getTraceHeader()->setFloatValue( vars->hdrID_scalar, resampleScalar );
  }
  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_resample_( csParamDef* pdef ) {
  pdef->setModule( "RESAMPLE", "Resample trace to different sample interval (LIMITED FUNCTIONALITY!)" );

  pdef->addParam( "sample_int", "New sample interval [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Sample interval [ms]" );

  pdef->addParam( "debias", "Remove DC bias before resampling?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Remove DC bias before resampling, reapply afterwards.", "This should be done to avoid FFT related artefacts due to DC bias" );
  pdef->addOption( "no", "Do not remove DC bias before resampling." );

  pdef->addParam( "anti_alias", "Apply anti-alias filter?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Apply anti-alias (Butterworth) filter. Specify input parameters 'cutoff' and 'slope' (or 'order')" );
  pdef->addOption( "no", "Do not apply anti-alias filter." );
  //  pdef->addOption( "fir", "Apply anti-alias filter. Specify input parameters 'cutoff' and 'slope' (or 'order')" );

  pdef->addParam( "slope", "Filter slope in dB/oct", NUM_VALUES_FIXED );
  pdef->addValue( "60", VALTYPE_NUMBER );

  pdef->addParam( "order", "Filter 'order'", NUM_VALUES_FIXED, "The filter 'order' defines the steepness of the filter taper. The higher the order, the shorter the taper." );
  pdef->addValue( "10", VALTYPE_NUMBER, "Filter order (1-100)" );

  pdef->addParam( "cutoff", "Cut-off (-3db) frequency", NUM_VALUES_FIXED );
  pdef->addValue( "0.8", VALTYPE_NUMBER, "Cut-off frequency, given as ratio of Nyquist (0.0-1.0)", "For example: 0.9 --> cut-off frequency is 90% of Nyquist" );

  pdef->addParam( "upsampling", "Interpolation method used for upsampling (output sample interval < input sample interval)", NUM_VALUES_FIXED );
  pdef->addValue( "sinc", VALTYPE_OPTION );
  pdef->addOption( "sinc", "Use sinc interpolation" );
  pdef->addOption( "quad", "Obsolete option. Use quadratic interpolation", "...for backward compatibility: This method was hardcoded in previous versions" );

  pdef->addParam( "norm", "Apply normalization?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Apply normalization" );
  pdef->addOption( "no", "Do not apply normalization" );
  pdef->addOption( "rms", "Obsolete option: Apply trace-by-trace RMS normalization. Inverse scalar is stored in output trace header 'resample_scalar'. Applying this scalar reverses the normalization" );

  //  pdef->addParam( "filter", "Anti-alias filter method to use", NUM_VALUES_FIXED );
  //  pdef->addValue( "butterworth", VALTYPE_OPTION );
  //  pdef->addOption( "butterworth", "Standard butterworth filter" );
  //  pdef->addOption( "fir", "Standard fir x-t filter" );
}


extern "C" void _params_mod_resample_( csParamDef* pdef ) {
  params_mod_resample_( pdef );
}
extern "C" void _init_mod_resample_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_resample_( param, env, log );
}
extern "C" bool _exec_mod_resample_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_resample_( trace, port, env, log );
}

/*
 * Create filter coefficients for standard FIR low-pass filter
 *
 * h(n) = w/pi  * sin(w*m)/(w*m)
 * m = (n - N/2)
 * N: Order
 * w = 2*pi*f
 * f = cut-off frequency
 */
void create_filter_coef( float cutOffFreqHz, int order, float* coef, int numCoef ) {
  double angularFreq = cutOffFreqHz * M_PI * 2.0;
  double orderHalf = (double)order / 2.0;

  for( int i = 0; i < order; i++ ) {
    double id = (double)i;
    coef[i] = (float)( sin( angularFreq * ( id - orderHalf ) ) / (M_PI * ( id - orderHalf )) );
  }
  for( int i = order+1; i < numCoef; i++ ) {
    double id = (double)i;
    coef[i] = (float)( sin( angularFreq * ( id - orderHalf ) ) / (M_PI * ( id - orderHalf )) );
  }
  coef[order] = (float)( 2.0*cutOffFreqHz );

  // Apply Hamming window:
  for( int i = 0; i < numCoef; i++ ) {
    coef[i] *= (float)( ( 0.54 - 0.46 * cos(2.0*M_PI*(double)i/(double)numCoef) ) );
  }
}



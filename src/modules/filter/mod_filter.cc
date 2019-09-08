/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csFFTTools.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: FILTER
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_filter {
  struct VariableStruct {
    float cutoffLow;
    float cutoffHigh;
    bool isLowPass;
    bool isHighPass;
    float order;
    csFFTTools* fftTool;
    csFFTTools** fftToolWin;
    int filterType;
    bool isPad;
    int paddedSamples;
    int numSamplesInclPad;
    float* bufferPaddedTrace;
    int output;
    float* bufferInput;
    int hdrID_both;

    // Windowed filter application
    int numWin; // Number of windows
    //    float* winTimesSamp; // Window end time in samples
    int* winStartSample; // Windows start samples (inclusive)
    int* winEndSample;   // Window end samples (inclusive)
    int* winNumSamples;  // Number of samples in windows
    int  winOverlapSamples; // Window overlap in number of samples 
    float* freqLowPass;
    float* freqHighPass;
    float* orderLowPass;
    float* orderHighPass;
    float* winBufferIn;
    float* winBufferOut;

    bool isNotchFilter;
    bool isNotchCosineTaper;
    float notchFreqHz;
    float notchWidthHz;
  };
  static int const UNIT_HZ   = 1;
  static int const UNIT_PERCENT = 2;
  static int const TYPE_BUTTER = 10;

  static int const OUTPUT_FILT = 41;
  static int const OUTPUT_DIFF = 42;
  static int const OUTPUT_BOTH = 43;
  static int const OUTPUT_IMPULSE = 44;
  static int const OUTPUT_NOTCH = 45;
}
using namespace mod_filter;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_filter_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  csTraceHeaderDef* hdef = env->headerDef;
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_MULTITRACE );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->cutoffLow     = 0;
  vars->cutoffHigh    = 0;
  vars->isLowPass     = false;
  vars->isHighPass    = false;
  vars->fftTool       = NULL;
  vars->order         = 4;
  vars->filterType = mod_filter::TYPE_BUTTER;

  vars->isPad             = false;
  vars->paddedSamples     = 0;
  vars->numSamplesInclPad = 0;
  vars->bufferPaddedTrace = NULL;
  vars->output = OUTPUT_FILT;
  vars->bufferInput = NULL;
  vars->hdrID_both = -1;

  vars->numWin = 0;
  vars->fftToolWin    = NULL;
  vars->winOverlapSamples = 0;
  vars->winStartSample = NULL;
  vars->winEndSample = NULL;
  vars->winNumSamples = NULL;
  vars->freqLowPass  = NULL;
  vars->freqHighPass = NULL;
  vars->orderLowPass  = NULL;
  vars->orderHighPass = NULL;
  vars->winBufferIn = NULL;
  vars->winBufferOut = NULL;

  vars->notchFreqHz = 0;
  vars->notchWidthHz = 0;
  vars->isNotchFilter = false;
  vars->isNotchCosineTaper = false;
//---------------------------------------------
//
  bool doRestore = false;
  std::string text;
  if( param->exists("mode") ) {
    param->getString("mode", &text );
    if( !text.compare("restore") ) {
      doRestore = true;
    }
    else if( !text.compare("apply") ) {
      doRestore = false;
    }
    else {
      log->error("Unknown option: %s", text.c_str() );
    }
  }

  if( param->exists("type") ) {
    param->getString("type", &text );
    if( !text.compare("butterworth") ) {
      vars->filterType = mod_filter::TYPE_BUTTER;
    }
    else {
      log->error("Unknown option: %s", text.c_str() );
    }
  }

//---------------------------------------------
//
  if( param->exists("pad") ) {
    float padLength;
    param->getFloat("pad", &padLength );
    vars->paddedSamples = (int)round( padLength / shdr->sampleInt );
    vars->numSamplesInclPad = 2 * vars->paddedSamples + shdr->numSamples;
    vars->bufferPaddedTrace = new float[vars->numSamplesInclPad];
    vars->isPad = true;
    if( shdr->numSamples/2 < vars->paddedSamples ) {
      log->error("Too much padding: Pad length cannot exceed 1/2 trace length. Maximum pad length = %fms", (float)shdr->numSamples*shdr->sampleInt/2.0f);
    }
  }
//---------------------------------------------
//

  if( param->exists("impulse") ) {
    std::string text;
    param->getString("impulse", &text);
    if( !text.compare("yes") ) {
      vars->output = mod_filter::OUTPUT_IMPULSE;
    }
    else if( !text.compare("no") ) {
      //        vars->outputImpulse = false;
    }
    else {
      log->error("Unknown option: '%s'", text.c_str() );
    }
  }
  
  int freqUnit = mod_filter::UNIT_HZ;
  if( param->exists("unit") ) {
    std::string text;
    param->getString("unit", &text);
    if( !text.compare("hz") ) {
      freqUnit = mod_filter::UNIT_HZ;
    }
    else if( !text.compare("percent") ) {
      freqUnit = mod_filter::UNIT_PERCENT;
    }
    else {
      log->error("Unknown option: '%s'", text.c_str() );
    }
  }

  if( param->exists("win_times") ) {
    if( vars->isPad ) log->error("Windowed filter application does not work in conjunction with padding. Remove user parameter 'pad' and retry.");
    if( doRestore ) log->error("Windowed filter application does not suppport  the 'mode restore' option. Remove user parameter 'mode' and retry.");
    if( vars->output == mod_filter::OUTPUT_IMPULSE ) log->error("Windowed filter application does not work in conjunction with impulse output option. Remove user parameter 'impulse' and retry.");
    if( freqUnit == mod_filter::UNIT_PERCENT ) log->error("Windowed filter application does not support percent option. Remove user parameter 'unit' and retry.");
    vars->numWin = param->getNumValues("win_times");
    int numLow  = param->getNumLines("lowpass");
    int numHigh = param->getNumLines("highpass");
    if( vars->numWin == 0 || (numLow!=vars->numWin && numHigh!=vars->numWin) ) {
      log->error("Inconsistent window definition. Need to specify same number of low (=%d) and/or high (=%d) pass points as window start times (=%d)", numLow, numHigh, vars->numWin);
    }
    if( numLow != 0 ) vars->isLowPass = true;
    if( numHigh != 0 ) vars->isHighPass = true;

    vars->winOverlapSamples = (int)(500 / shdr->sampleInt);
    if( param->exists("win_overlap") ) {
      float overlap_ms = 0;
      param->getFloat("win_overlap",&overlap_ms);
      vars->winOverlapSamples = (int)round( 0.5*overlap_ms / shdr->sampleInt ) * 2;
    }

    vars->winStartSample = new int[vars->numWin];
    vars->winEndSample   = new int[vars->numWin];
    vars->freqLowPass    = new float[vars->numWin];
    vars->freqHighPass   = new float[vars->numWin];
    vars->orderLowPass   = new float[vars->numWin];
    vars->orderHighPass  = new float[vars->numWin];
    float slope;
    for( int i = 0; i < vars->numWin; i++ ) {
      float time;
      param->getFloat("win_times", &time, i );
      vars->winStartSample[i] = (int)(round(time / shdr->sampleInt));
      if( vars->isLowPass ) {
        param->getFloatAtLine("lowpass", &vars->freqLowPass[i], i, 0 );
        param->getFloatAtLine("lowpass", &slope, i, 1 );
        vars->orderLowPass[i] = fabs(slope) / 6.0;
      }
      if( vars->isHighPass ) {
        param->getFloatAtLine("highpass", &vars->freqHighPass[i], i, 0 );
        param->getFloatAtLine("highpass", &slope, i, 1 );
        vars->orderHighPass[i] = fabs(slope) / 6.0;
      }
    }

    vars->winStartSample[0] = 0;
    vars->winEndSample[vars->numWin-1] = shdr->numSamples-1;
    for( int i = 1; i < vars->numWin; i++ ) {
      vars->winEndSample[i-1] = vars->winStartSample[i] + vars->winOverlapSamples/2;
      vars->winStartSample[i] -= vars->winOverlapSamples/2;
      if( vars->winStartSample[i] < 0 ) log->error("Window #%d, including overlap (=%.2fms), is too large.", i+1, vars->winOverlapSamples*shdr->sampleInt);
      if( vars->winEndSample[i-1] >= shdr->numSamples ) log->error("Window #%d, including overlap (=%.2fms), is too large.", i+1, vars->winOverlapSamples*shdr->sampleInt);
    }
    vars->fftToolWin = new csFFTTools*[vars->numWin];
    for( int i = 0; i < vars->numWin; i++ ) {
      int numSamplesWin = vars->winEndSample[i] - vars->winStartSample[i] + 1;
      vars->fftToolWin[i] = new csFFTTools( numSamplesWin, shdr->sampleInt );
    }
    vars->winBufferIn  = new float[shdr->numSamples];
    vars->winBufferOut = new float[shdr->numSamples];

    for( int iwin = 0; iwin < vars->numWin; iwin++ ) {
      log->line("Window #%d: Start time/sample %8.2fms/%-5d, end time/sample: %8.2fms/%-5d, number of samples: %d", iwin+1,
              shdr->sampleInt*vars->winStartSample[iwin], vars->winStartSample[iwin], shdr->sampleInt*vars->winEndSample[iwin], vars->winEndSample[iwin], vars->winEndSample[iwin]-vars->winStartSample[iwin]+1 );
    }
  } // END win_times

  //--------------------------------------------------------------------------------
  //
  if( vars->numWin == 0 && vars->filterType == mod_filter::TYPE_BUTTER ) {
    if( param->exists("highpass") ) {
      param->getFloat("highpass", &vars->cutoffHigh );
      vars->isHighPass = true;
    }

    if( param->exists("lowpass") ) {
      param->getFloat("lowpass", &vars->cutoffLow );
      vars->isLowPass = true;
    }

    if( param->exists("slope") ) {
      float slope;
      param->getFloat("slope", &slope );
      vars->order = fabs(slope) / 6.0;
    }
    if( param->exists("order") ) {
      if( param->exists("slope") ) {
	log->error("Specify only one of the following two user parameters: 'slope' and 'order'");
      }
      param->getFloat("order", &vars->order );
    }
    if( vars->order < 0.1 ) { // || vars->order > 100 ) {
      log->error("Specified order = %f out of valid range = 0.1-...", vars->order);
    }
    if( doRestore ) vars->order = -vars->order;

    if( param->exists("notch_filter") ) {
      param->getFloat("notch_filter", &vars->notchFreqHz, 0);
      param->getFloat("notch_filter", &vars->notchWidthHz, 1);
      std::string text;
      param->getString("notch_filter", &text, 2);
      if( !text.compare("cosine") ) {
        vars->isNotchCosineTaper = true;
      }
      else if( !text.compare("butter") ) {
        vars->isNotchCosineTaper = false;
      }
      else {
        log->error("Unknown option: '%s'", text.c_str() );
      }
      vars->isNotchFilter = true;
    }

    float freqNy = 500.0f/shdr->sampleInt;
    if( freqUnit == mod_filter::UNIT_PERCENT ) {
      if( vars->cutoffHigh < 0 || vars->cutoffHigh > 100 ) {
        log->error("Filter frequency specified as percent (%) exceeds valid range (0-100): %f", vars->cutoffHigh );
      }
      if( vars->cutoffLow < 0 || vars->cutoffLow > 100 ) {
        log->error("Filter frequency specified as percent (%) exceeds valid range (0-100): %f", vars->cutoffLow );
      }
      
      vars->cutoffHigh *= freqNy / 100.0f;
      vars->cutoffLow  *= freqNy / 100.0f;
    }
    if( vars->cutoffHigh < 0 || vars->cutoffHigh > freqNy ) {
      log->error("Filter frequency exceeds valid range (0-%f): %f", freqNy, vars->cutoffHigh );
    }
    if( vars->cutoffLow < 0 || vars->cutoffLow > freqNy ) {
      log->error("Filter frequency exceeds valid range (0-%f): %f", freqNy, vars->cutoffLow );
    }

    if( !vars->isLowPass && !vars->isHighPass && !vars->isNotchFilter ) {
      log->error("No filter option specified. Specify for user parameter 'lowpass' and/or 'highpass' and/or 'notch_filter'");
    }

    if( !vars->isPad ) {
      vars->fftTool = new csFFTTools( shdr->numSamples, shdr->sampleInt );
    }
    else {
      vars->fftTool = new csFFTTools( vars->numSamplesInclPad, shdr->sampleInt );
    }
  } // END: Setup bandpass filter

  if( param->exists("output") ) {
    if( param->exists("impulse") ) {
      log->error("Use parameter 'output' to specify output option. Remove obsolete parameter 'impulse' from flow");
    }
    param->getString("output", &text );
    if( !text.compare("filt") ) {
      vars->output = OUTPUT_FILT;
    }
    else if( !text.compare("diff") ) {
      vars->output = OUTPUT_DIFF;
    }
    else if( !text.compare("both") ) {
      vars->output = OUTPUT_BOTH;
    }
    else if( !text.compare("impulse") ) {
      vars->output = OUTPUT_IMPULSE;
    }
    else if( !text.compare("notch_filter") ) {
      vars->output = OUTPUT_NOTCH;
      shdr->numSamples = vars->fftTool->numFFTSamples();
      shdr->sampleInt  = vars->fftTool->sampleIntFreqHz();
      shdr->domain = DOMAIN_FX;
    }
    else {
      log->error("Unknown option: %s", text.c_str() );
    }
    if( vars->output == mod_filter::OUTPUT_DIFF || vars->output == mod_filter::OUTPUT_BOTH ) {
      vars->bufferInput = new float[shdr->numSamples];
    }
    if( vars->output == mod_filter::OUTPUT_BOTH ) {
      text = "trcno";
      if( param->getNumValues("output") > 1 ) {
        param->getString("output",&text,1);
      }
      vars->hdrID_both = hdef->headerIndex( text.c_str() );
    }
  }

  if( vars->isNotchFilter ) {
    vars->fftTool->setupNotchFilter( vars->notchFreqHz, vars->notchWidthHz, vars->order, vars->isNotchCosineTaper ); 
    //    for( int i = 0; i < vars->fftTool->numFFTSamples(); i++ ) {
    //   fprintf(stdout,"%f %f\n", vars->fftTool->sampleIntFreqHz()*i, ptr[i]);
    //  }
    // log->error("Whatever");
  }
}
//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_filter_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeaderDef const* hdef = env->headerDef;

  if( edef->isCleanup()){
    if( vars->fftTool != NULL ) {
      delete vars->fftTool;
      vars->fftTool = NULL;
    }
    if( vars->bufferPaddedTrace != NULL ) {
      delete [] vars->bufferPaddedTrace;
      vars->bufferPaddedTrace = NULL;
    }
    if( vars->bufferInput != NULL ) {
      delete [] vars->bufferInput;
      vars->bufferInput = NULL;
    }
    if( vars->winBufferIn != NULL ) {
      delete [] vars->winBufferIn;
      vars->winBufferIn = NULL;
    }
    if( vars->winBufferOut != NULL ) {
      delete [] vars->winBufferOut;
      vars->winBufferOut = NULL;
    }
    if( vars->fftToolWin != NULL ) {
      for( int iwin = 0; iwin < vars->numWin; iwin++ ) {
        delete vars->fftToolWin[iwin];
      }
      delete [] vars->fftToolWin;
      vars->fftToolWin = NULL;
    }
    if( vars->winStartSample != NULL ) {
      delete [] vars->winStartSample;
      vars->winStartSample= NULL;
    }
    if( vars->winEndSample != NULL ) {
      delete [] vars->winEndSample;
      vars->winEndSample = NULL;
    }
    if( vars->winNumSamples != NULL ) {
      delete [] vars->winNumSamples;
      vars->winNumSamples = NULL;
    }
    if( vars->freqLowPass  != NULL ) {
      delete [] vars->freqLowPass;
      vars->freqLowPass  = NULL;
    }
    if( vars->freqHighPass != NULL ) {
      delete [] vars->freqHighPass;
      vars->freqHighPass = NULL;
    }
    if( vars->orderLowPass != NULL ) {
      delete [] vars->orderLowPass;
      vars->orderLowPass  = NULL;
    }
    if( vars->orderHighPass != NULL ) {
      delete [] vars->orderHighPass;
      vars->orderHighPass = NULL;
    }
    delete vars; vars = NULL;
    return;
  }

  csTrace* trace = traceGather->trace(0);
  float* samplesInOut = trace->getTraceSamples();

  if( vars->output == mod_filter::OUTPUT_DIFF || vars->output == mod_filter::OUTPUT_BOTH ) {
    memcpy( vars->bufferInput, samplesInOut, shdr->numSamples*sizeof(float) );
  }

  // Padding: Pad samples and extrapolate data samples. Apply cosine taper to padded data.
  if( vars->isPad ) {
    for( int isamp = 0; isamp < vars->numSamplesInclPad; isamp++ ) {
      vars->bufferPaddedTrace[isamp] = 0;
    }
    // Copy input data
    memcpy( &vars->bufferPaddedTrace[vars->paddedSamples], samplesInOut, shdr->numSamples * sizeof(float) );
    // Extrapolate start of trace
    float value2x = 2 * samplesInOut[0];
    for( int isamp = vars->paddedSamples-1; isamp >= 0; isamp-- ) {
      float ratio = (float)isamp/(float)(vars->paddedSamples-1);
      float taper = 0.5 * ( 1 + cos( M_PI*(ratio-1.0) ) );
      vars->bufferPaddedTrace[isamp] = taper * ( value2x - samplesInOut[vars->paddedSamples-isamp] );
    }
    // Extrapolate end of trace
    value2x = 2 * samplesInOut[shdr->numSamples-1];
    int num = 2*shdr->numSamples+vars->paddedSamples-2;
    for( int isamp = shdr->numSamples+vars->paddedSamples; isamp < vars->numSamplesInclPad; isamp++ ) {
      float ratio = (float)(vars->numSamplesInclPad-isamp-1)/(float)(vars->paddedSamples-1);
      float taper = 0.5 * ( 1 + cos( M_PI*(ratio-1.0) ) );
      vars->bufferPaddedTrace[isamp] = taper * ( value2x - samplesInOut[num-isamp] );
    }
    // Redirect 'samples pointer' to padded trace
    samplesInOut = vars->bufferPaddedTrace;
    if( edef->isDebug() ) {
      for( int isamp = 0; isamp < vars->numSamplesInclPad; isamp++ ) {
        fprintf(stdout,"%f %f\n", isamp*shdr->sampleInt, samplesInOut[isamp]);
      }
    }
  }

  //----------------------------------------------------------------------
  // Windowed filter application
  if( vars->numWin != 0 ) {
    // Zero output buffer
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      vars->winBufferOut[isamp] = 0.0;
    }
    for( int iwin = 0; iwin < vars->numWin; iwin++ ) {
      int numSamplesWin = vars->winEndSample[iwin] - vars->winStartSample[iwin] + 1;
      // Copy windowed data to input buffer
      memcpy( vars->winBufferIn, &samplesInOut[vars->winStartSample[iwin]], numSamplesWin*sizeof(float) );
      // Apply low and/or high pass
      if( vars->isLowPass ) {
        vars->fftToolWin[iwin]->lowPass( vars->winBufferIn, vars->orderLowPass[iwin], vars->freqLowPass[iwin], (vars->output == OUTPUT_IMPULSE) );
      }
      if( vars->isHighPass ) {
        vars->fftToolWin[iwin]->highPass( vars->winBufferIn, vars->orderHighPass[iwin], vars->freqHighPass[iwin], (vars->output == OUTPUT_IMPULSE) );
      }
      // a) Copy non-overlap data to output buffer using memcpy
      int samp1Copy = vars->winStartSample[iwin] + vars->winOverlapSamples;
      int samp2Copy = vars->winEndSample[iwin]   - vars->winOverlapSamples;
      int sampFrom  = vars->winOverlapSamples;
      if( iwin == 0 ) {
        samp1Copy = 0;
        sampFrom = 0;
      }
      else if( iwin == vars->numWin-1 ) {
        samp2Copy = vars->winEndSample[iwin];
      }
      memcpy( &vars->winBufferOut[samp1Copy], &vars->winBufferIn[sampFrom], (samp2Copy-samp1Copy+1)*sizeof(float) );
      // b) Add overlap data to output buffer using linear taper
      
      if( iwin > 0 ) {
        for( int isamp = 0; isamp < vars->winOverlapSamples; isamp++ ) {
          float scalar = (float)(isamp) / (float)(vars->winOverlapSamples-1);
          vars->winBufferOut[isamp+vars->winStartSample[iwin]] += scalar*vars->winBufferIn[isamp];
        }
      }
      if( iwin < vars->numWin-1 ) {
        for( int isamp = 1; isamp <= vars->winOverlapSamples; isamp++ ) {
          float scalar = (float)(isamp-1) / (float)(vars->winOverlapSamples-1);
          vars->winBufferOut[vars->winEndSample[iwin]-isamp+1] += scalar*vars->winBufferIn[numSamplesWin-isamp];
        }
      }
      
    } // END: Windowed filter application
    memcpy( samplesInOut, vars->winBufferOut, shdr->numSamples*sizeof(float) );
  }
  //----------------------------------------------------------------------
  // Non-windowed filter application
  else if( vars->filterType == mod_filter::TYPE_BUTTER ) {
    if( vars->isLowPass ) {
      vars->fftTool->lowPass( samplesInOut, vars->order, vars->cutoffLow, (vars->output == OUTPUT_IMPULSE) );
    }
    if( vars->isHighPass ) {
      vars->fftTool->highPass( samplesInOut, vars->order, vars->cutoffHigh, (vars->output == OUTPUT_IMPULSE) );
    }
  }
  if( vars->isNotchFilter ) {
    if( vars->output == mod_filter::OUTPUT_NOTCH ) {
      double const* ptr = vars->fftTool->setupNotchFilter( vars->notchFreqHz, vars->notchWidthHz, vars->order, vars->isNotchCosineTaper );
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samplesInOut[isamp] = (float)ptr[isamp];
      }
      return;
    }
    vars->fftTool->notchFilter( samplesInOut, false );
  }

  // Padding: Remove padded samples
  if( vars->isPad ) {
    memcpy( trace->getTraceSamples(), &samplesInOut[vars->paddedSamples], shdr->numSamples * sizeof(float) );    
    samplesInOut = trace->getTraceSamples(); // Redirect pointer to actual trace
  }

  //--------------------------------------------------------------------------------
  //
  if( vars->output == mod_filter::OUTPUT_DIFF || vars->output == mod_filter::OUTPUT_BOTH ) {
    float* samplesOut = samplesInOut;
    if( vars->output == mod_filter::OUTPUT_BOTH ) {
      traceGather->createTrace( hdef, shdr->numSamples);
      traceGather->trace(1)->getTraceHeader()->copyFrom( trace->getTraceHeader() );
      samplesOut = traceGather->trace(1)->getTraceSamples();
      traceGather->trace(0)->getTraceHeader()->setIntValue( vars->hdrID_both,1 );
      traceGather->trace(1)->getTraceHeader()->setIntValue( vars->hdrID_both,2 );
    }
    // Compute difference:
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samplesOut[isamp] = vars->bufferInput[isamp] - samplesInOut[isamp];
    }
  }

}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_filter_( csParamDef* pdef ) {
  pdef->setModule( "FILTER", "Frequency filter" );

  pdef->addParam( "type", "Filter type", NUM_VALUES_FIXED );
  pdef->addValue( "butterworth", VALTYPE_OPTION );
  pdef->addOption( "butterworth", "Butterworth filter" );

  pdef->addParam( "lowpass", "Lowpass filter", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Cutoff frequency for low-pass filter [Hz]",
     "The cutoff frequency will be damped by -3db" );
  pdef->addValue( "12", VALTYPE_NUMBER, "Filter slope [dB/oct]" );

  pdef->addParam( "highpass", "Highpass filter", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Cutoff frequency for highpass filter [Hz]",
     "The cutoff frequency will be damped by -3db" );
  pdef->addValue( "12", VALTYPE_NUMBER, "Filter slope [dB/oct]" );

  pdef->addParam( "impulse", "Output filter impulse response", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Output filter impulse response, zero phase, placed at trace centre." );
  pdef->addOption( "no", "Do not output filter impulse response. (i.e. output filtered input data)" );

  pdef->addParam( "unit", "Unit of frequency values supplied in parameters", NUM_VALUES_FIXED );
  pdef->addValue( "hz", VALTYPE_OPTION );
  pdef->addOption( "hz", "Frequencies are specified in [Hz]" );
  pdef->addOption( "percent", "Frequencies are specified as percent of Nyquist" );

  pdef->addParam( "order", "Filter 'order'", NUM_VALUES_FIXED, "...can also be specified as slope/oct, see user parameter 'slope'. For backward compatibility, this parameter overrides the slope provided as the second value in the user parameters 'highpass' and 'lowpass'" );
  pdef->addValue( "4", VALTYPE_NUMBER, "Filter order (1-100)" );

  pdef->addParam( "slope", "Filter slope in dB/octave", NUM_VALUES_FIXED, "...can also be specified as filter 'order', see user parameter 'order'. For backward compatibility, this parameter overrides the slope provided as the second value in the user parameters 'highpass' and 'lowpass'" );
  pdef->addValue( "12", VALTYPE_NUMBER, "Filter slope [dB/oct]", "Will be converted into filter order (= 1/6 x slope)" );

  pdef->addParam( "pad", "Pad trace to avoid filter edge effects", NUM_VALUES_FIXED, "Pad and extrapolate trace at top and bottom before filter application. Remove padded samples after filter." );
  pdef->addValue( "0", VALTYPE_NUMBER, "Pad length at top and bottom, in units of trace (for example [ms])" );

  pdef->addParam( "mode", "Mode of operation: Apply specified filter or inverse filter", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "apply", "Apply filter" );
  pdef->addOption( "restore", "'Restore' filter: Apply inverse filter" );

  pdef->addParam( "output", "Output data", NUM_VALUES_VARIABLE );
  pdef->addValue( "filt", VALTYPE_OPTION );
  pdef->addOption( "filt", "Output filtered data" );
  pdef->addOption( "impulse", "Output filter impulse response, zero phase, placed at trace centre." );
  pdef->addOption( "diff", "Output difference between input and filtered data" );
  pdef->addOption( "both", "Output both filtered data and difference between input and filtered data" );
  pdef->addOption( "notch_filter", "Output notch filter" );
  pdef->addValue( "trcno", VALTYPE_STRING, "Optional, used for option 'both': Trace header to distinguish between filtered data (1) and difference (2)" );

  pdef->addParam( "win_times", "List of N window start times for windowed filter application", NUM_VALUES_VARIABLE, "Specify N lines with lowpass and/or highpass filter points, user parameters 'lowpass' and 'highpass'" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of times [ms]" );

  pdef->addParam( "win_overlap", "Overlap between windows", NUM_VALUES_FIXED, "Overlap is centred around window start times" );
  pdef->addValue( "500", VALTYPE_NUMBER, "Window time overlap [ms]" );

  pdef->addParam( "notch_filter", "Apply notch filter", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Notch frequency [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Width of notch filter [Hz]" );
  pdef->addValue( "butter", VALTYPE_OPTION, "Filter type" );
  pdef->addOption( "butter", "Use Butterworth filter. Specify user parameter 'slope' or 'order' for filter slope" );
  pdef->addOption( "cosine", "Use cosine taper" );
}


extern "C" void _params_mod_filter_( csParamDef* pdef ) {
  params_mod_filter_( pdef );
}
extern "C" void _init_mod_filter_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_filter_( param, env, log );
}
extern "C" void _exec_mod_filter_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_filter_( traceGather, port, numTrcToKeep, env, log );
}



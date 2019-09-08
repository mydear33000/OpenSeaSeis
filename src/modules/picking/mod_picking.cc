/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


#include "cseis_includes.h"
#include "geolib_methods.h"
#include "csSort.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: PICKING
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_picking {
  struct VariableStruct {
    int startSamp;
    int endSamp;
    //    float startSampFloat;
    //    float endSampFloat;
    float startTime;  // Start/end time may be different from start/end sample, i.e. they may be in between full samples
    float endTime;
    float threshold;
    float thresh_sign;
    int method;
    int modeSearch;
    int hdrId_timePick;
    int hdrId_amplitude;
    bool forcePick;

    bool isCross;
    int xcross_method;
    float xcross_thresholdTime;
    int xcross_startSamp;
    int xcross_endSamp;
    bool outputCross;
    float* xcrossBuffer;

    int hdrId_startTimeSample;
    int hdrId_endTimeSample;
    bool isSampleDomain;

    float valueTolerance;
  };
  static int const METHOD_PEAK_TROUGH   = 11;
  static int const METHOD_ZERO_CROSSING = 22;
  static int const METHOD_VALUE         = 23;
  
  static int const MODE_FIRST  = 101;
  static int const MODE_MAX    = 102;

  static int const XCROSS_ABSOLUTE = 201;
  static int const XCROSS_RELATIVE = 202;
}
using namespace mod_picking;

float cross_correlation( float const* s1, float const* s2, int firstSample, int numSamplesWindow, float* xcross, int numSamplesAll );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_picking_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_MULTITRACE );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->startSamp   = 0;
  vars->endSamp     = 0;
  vars->startTime   = 0.0;
  vars->endTime     = 0.0;
  vars->threshold   = 0;
  vars->thresh_sign = 0;
  vars->method      = 0;
  vars->modeSearch  = mod_picking::MODE_MAX;
  vars->hdrId_timePick = -1;
  vars->hdrId_amplitude = -1;
  vars->forcePick   = true;
  vars->isCross = false;
  vars->xcross_method = -1;
  vars->xcross_thresholdTime = 0;
  vars->xcross_startSamp = 0;
  vars->xcross_endSamp   = 0;
  vars->outputCross = false;
  vars->xcrossBuffer = NULL;
  vars->valueTolerance = 0.0f;
  vars->hdrId_startTimeSample = -1;
  vars->hdrId_endTimeSample = -1;
  vars->isSampleDomain = false;

  //---------------------------------------------------------
  std::string text;

  if( param->exists( "domain" ) ) {
    param->getString( "domain", &text );
    if( !text.compare( "sample" ) ) {
      vars->isSampleDomain = true;
    }
    else if( !text.compare( "time" ) || !text.compare( "freq" ) ) {
      vars->isSampleDomain = false;
    }
    else {
      log->line("Domain not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }

  csFlexNumber numberStart;
  csFlexNumber numberEnd;

  param->getString( "start", &text );
  if( !numberStart.convertToNumber( text ) ) {
    if( !hdef->headerExists(text) ) {
      log->error("Specified trace header name containing start time/sample does not exist: '%s'", text.c_str());
    }
    vars->hdrId_startTimeSample = hdef->headerIndex(text);
  }

  param->getString( "end", &text );
  if( !numberEnd.convertToNumber( text ) ) {
    if( !hdef->headerExists(text) ) {
      log->error("Specified trace header name containing end time/sample does not exist: '%s'", text.c_str());
    }
    vars->hdrId_endTimeSample = hdef->headerIndex(text);
  }

  if( vars->hdrId_startTimeSample != vars->hdrId_endTimeSample ) {
    if( vars->hdrId_startTimeSample < 0 || vars->hdrId_endTimeSample < 0 ) {
      log->error("Start & end time/sample: Both need to be specified as constant values or as input trace header");
    }
  }

  if( vars->isSampleDomain ) {
    //
    // NOTE: User input is '1' for first sample. Internally, '0' is used!!
    //
    if( vars->hdrId_startTimeSample < 0 ) {
      vars->startSamp = numberStart.intValue();
      if( vars->startSamp < 1 ) log->error("Start sample (%d) needs to be greater or equal to 1.", vars->startSamp);
      vars->startSamp -= 1;   // see note above..
      vars->startTime = (float)vars->startSamp * shdr->sampleInt;
    }
    if( vars->hdrId_endTimeSample < 0 ) {
      vars->endSamp = numberEnd.intValue();
      if( vars->endSamp > shdr->numSamples ) {
	log->warning("Specified end sample (%d) exceeds number of samples in trace (%d). Changed accordingly...", vars->endSamp, shdr->numSamples);
	vars->endSamp = shdr->numSamples;
      }
      if( vars->endSamp < vars->startSamp ) {
	log->error("End sample (%d) needs to be larger than start sample (%d).", vars->endSamp, vars->startSamp );
      }
      vars->endSamp   -= 1;
      vars->endTime   = (float)vars->endSamp * shdr->sampleInt;
    }
  }
  else {
    if( vars->hdrId_startTimeSample < 0 ) {
      vars->startTime = numberStart.floatValue();
      if( vars->startTime < 0.0 ) log->error("Start time (%f) needs to be greater or equal to 0.0.", vars->startTime);
      vars->startSamp = (int)(vars->startTime / shdr->sampleInt + 0.01);
    }
    if( vars->hdrId_endTimeSample < 0 ) {
      vars->endTime = numberEnd.floatValue();
      if( vars->endTime > shdr->numSamples*shdr->sampleInt ) {
	log->warning("Specified end time/freq (%f) exceeds trace length (%f). Changed accordingly...", vars->endTime, shdr->numSamples*shdr->sampleInt);
	vars->endTime = shdr->numSamples*shdr->sampleInt;
      }
      if( vars->startTime > vars->endTime ) log->error("Start time (%f) needs to be smaller than end time (%f).", vars->startTime, vars->endTime);
      vars->startSamp = (int)(vars->startTime / shdr->sampleInt + 0.01);
      vars->endSamp   = (int)(vars->endTime / shdr->sampleInt + 0.01);
    }
  }

  //-------------------------------------------------------------
  //
  param->getString( "method", &text );
  if( !text.compare("peak_trough") ) {
    vars->method = METHOD_PEAK_TROUGH;
  }
  else if( !text.compare("zero_crossing") ) {
    vars->method = METHOD_ZERO_CROSSING;
  }
  else if( !text.compare("value") ) {
    vars->method = METHOD_VALUE;
  }
  else {
    log->error("Option not recognised: '%s'", text.c_str());
  }

  //-------------------------------------------------------------
  // Cross-correlation
  //
  if( param->exists("xcross") ) {
    vars->isCross = true;
    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
    float threshold;
    param->getFloat( "xcross", &threshold, 0 );
    param->getString( "xcross", &text, 1 );
    if( !text.compare("absolute") ) {
      vars->xcross_method = XCROSS_ABSOLUTE;
    }
    else if( !text.compare("relative") ) {
      vars->xcross_method = XCROSS_RELATIVE;
    }
    else {
      log->error("Option not recognised: '%s'", text.c_str());
    }

    if( vars->isSampleDomain ) {
      vars->xcross_thresholdTime = threshold*shdr->sampleInt;
    }
    else {
      vars->xcross_thresholdTime = threshold;
    }
    if( vars->xcross_thresholdTime < 0 || vars->xcross_thresholdTime/shdr->sampleInt >= shdr->numSamples ) {
      log->error("Cross-correlation threshold out of range: %f  (unit of threshold depends on 'domain')", threshold );
    }
    vars->xcross_startSamp = vars->startSamp;
    vars->xcross_endSamp   = vars->endSamp;
    if( edef->isDebug() ) log->line("Cross-correlation window: start samp/end samp/number of samples: %d %d\n",
                                    vars->xcross_startSamp, vars->xcross_endSamp );
  }

  if( param->exists("xcross_win") ) {
    if( !vars->isCross ) {
      log->error("Cross-correlation window given (parameter 'xcross_win'), but parameter 'xcross' not specified.");
    }
    float start;
    float end;
    param->getFloat("xcross_win", &start, 0);
    param->getFloat("xcross_win", &end, 1);
    if( vars->isSampleDomain ) {
      if( start != 0 ) {
        vars->xcross_startSamp = (int)(start + 0.5);
      }
      if( end != 0 ) {
        vars->xcross_endSamp   = (int)(end + 0.5);
      }
    }
    else {
      if( start != 0 ) {
        vars->xcross_startSamp = (int)(start/shdr->sampleInt + 0.5);
      }
      if( end != 0 ) {
        vars->xcross_endSamp   = (int)(end/shdr->sampleInt + 0.5);
      }
    }
    if( vars->xcross_startSamp >= vars->xcross_endSamp ) {
      log->error("Inconsistent cross-correlation window:  start=%f  >=  end=%", vars->xcross_startSamp, vars->xcross_endSamp );
    }
    if( vars->xcross_endSamp-vars->xcross_startSamp >= shdr->numSamples-10 ) {
      log->error("Specified cross-correlation window is too large. Maximum +/- cross-correlation lag is constrained by number of samples before/after specified xcross window. Max lag ensures there are always %d samples in cross-correlation window. Number of samples in trace: %d", vars->xcross_endSamp-vars->xcross_startSamp+1, shdr->numSamples );
    }
    if( vars->xcross_startSamp < 0 ) {
      vars->xcross_startSamp = 0;
    }
    if( vars->xcross_endSamp >= shdr->numSamples ) {
      vars->xcross_endSamp = shdr->numSamples-1;
    }
  }

  if( param->exists("xcross_output") ) {
    param->getString( "xcross_output", &text, 1 );
    if( !text.compare("yes") ) {
      vars->outputCross = true;
    }
    else if( !text.compare("no") ) {
      vars->outputCross = false;
    }
    else {
      log->error("Option not recognised: '%s'", text.c_str());
    }
  }
  if( vars->isCross ) {
    int numSamplesWindow = vars->xcross_endSamp - vars->xcross_startSamp + 1;
    vars->xcrossBuffer = new float[numSamplesWindow];
    for( int iswin = 0; iswin < numSamplesWindow; iswin++ ) {
      vars->xcrossBuffer[iswin] = 0.0;
    }
  }

  //-------------------------------------------------------------
  //
  if( param->exists( "search_mode" ) ) {
    param->getString( "search_mode", &text );
    if( !text.compare("first") ) {
      vars->modeSearch = MODE_FIRST;
    }
    else if( !text.compare("max") ) {
      vars->modeSearch = MODE_MAX;
    }
    else {
      log->error("Option not recognised: '%s'", text.c_str());
    }
  }

  //-------------------------------------------------------------
  //
  if( param->exists("force_pick") ) {
    param->getString( "force_pick", &text );
    if( !text.compare("yes") ) {
      vars->forcePick = true;
    }
    else if( !text.compare("no") ) {
      vars->forcePick = false;
    }
    else {
      log->error("Option not recognised: '%s'", text.c_str());
    }
  }
  //-------------------------------------------------------------
  //
  if( param->exists("thresh") ) {
    param->getFloat( "thresh", &vars->threshold );
    if( vars->threshold < 0.0 ) {
      vars->thresh_sign = -1.0;
      vars->threshold *= -1.0;
    }
    else {
      vars->thresh_sign = 1.0;
    }
    if( vars->method == METHOD_VALUE ) {
      vars->threshold *= vars->thresh_sign; // Keep original value whether pos or neg
      if( param->getNumValues("thresh") > 1 ) {
	param->getFloat( "thresh", &vars->valueTolerance, 1 );
      }
    }
  }

  //-------------------------------------------------------------
  //
  string headerTimePick("time_pick");
  string headerAmplitude("amp_pick");
  if( param->exists( "header" ) ) {
    int numValues = param->getNumValues("header");     
    param->getString( "header", &headerTimePick, 0 );
    if( vars->method == METHOD_PEAK_TROUGH ) {
      if( numValues > 1 ) {
        param->getString( "header", &headerAmplitude, 1 );
      }
    }
  }
  if( hdef->headerExists( headerTimePick ) ) {
    type_t type = hdef->headerType( headerTimePick );
    if( type != TYPE_FLOAT && type != TYPE_DOUBLE ) {
      log->error("Trace header '%s' exists but has wrong type. Should be FLOAT or DOUBLE", headerTimePick.c_str() );
    }
  }
  else {
    hdef->addHeader( TYPE_FLOAT, headerTimePick, "Time pick [ms]" );
  }
  if( vars->method == METHOD_PEAK_TROUGH ) {
    if( hdef->headerExists( headerAmplitude ) ) {
      type_t type = hdef->headerType( headerAmplitude );
      if( type != TYPE_FLOAT && type != TYPE_DOUBLE ) {
	log->error("Trace header '%s' exists but has wrong type. Should be FLOAT or DOUBLE", headerAmplitude.c_str() );
      }
    }
    else {
      hdef->addHeader( TYPE_FLOAT, headerAmplitude, "Amplitude at time pick" );
    }
  }
  //--------------------------------------------------

  vars->hdrId_timePick  = hdef->headerIndex( headerTimePick );
  if( vars->method == METHOD_PEAK_TROUGH ) vars->hdrId_amplitude = hdef->headerIndex( headerAmplitude );

  if( edef->isDebug() ) {
    log->line("Start/end sample INDEX for picking: %d %d, start/end time/freq: %f %f", vars->startSamp, vars->endSamp, vars->startTime, vars->endTime );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_picking_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup() ) {
    if( vars->xcrossBuffer != NULL ) {
      delete [] vars->xcrossBuffer;
      vars->xcrossBuffer = NULL;
    }
    delete vars; vars = NULL;
    return;
  }

    int ntraces = traceGather->numTraces();
  //  if( !vars->isCross ) {
    for( int itrc = 0; itrc < ntraces; itrc++ ) {
      csTrace* trace = traceGather->trace(itrc);
      float* samples = trace->getTraceSamples();
      int samplePickedInt = -1;  // Sample index of picked sample, either min/max etc
      if( vars->hdrId_startTimeSample != -1 ) {
	if( vars->isSampleDomain ) {
	  vars->startSamp = trace->getTraceHeader()->intValue( vars->hdrId_startTimeSample ) - 1;
	  vars->startTime = (float)vars->startSamp * shdr->sampleInt;
	}
	else {
	  vars->startTime = trace->getTraceHeader()->floatValue( vars->hdrId_startTimeSample );
	  vars->startSamp = (int)(vars->startTime / shdr->sampleInt + 0.01);
	}
	if( vars->hdrId_endTimeSample < 0  && vars->startSamp > vars->endSamp ) {
	  vars->startSamp = vars->endSamp;
	  vars->startTime = (float)vars->startSamp * shdr->sampleInt;
	}
      }
      if( vars->hdrId_endTimeSample != -1 ) {
	if( vars->isSampleDomain ) {
	  vars->endSamp = trace->getTraceHeader()->intValue( vars->hdrId_endTimeSample ) - 1;
	  vars->endTime = (float)vars->endSamp * shdr->sampleInt;
	}
	else {
	  vars->endTime = trace->getTraceHeader()->floatValue( vars->hdrId_endTimeSample );
	  vars->endSamp = (int)(vars->endTime / shdr->sampleInt + 0.01);
	}
	if( vars->startSamp > vars->endSamp ) {
	  vars->endSamp = vars->startSamp;
	  vars->endTime = (float)vars->endSamp * shdr->sampleInt;
	}
      }

      if( vars->method == METHOD_VALUE ) {
        // Search for first sample value that matches threshold
	float threshMin = vars->threshold - vars->valueTolerance;
	float threshMax = vars->threshold + vars->valueTolerance;
        for( int i = vars->startSamp; i <= vars->endSamp; i++ ) {
          if( samples[i] >= threshMin && samples[i] <= threshMax ) {
            samplePickedInt = i;
            break;
          }
        }
        trace->getTraceHeader()->setFloatValue( vars->hdrId_timePick, (float)samplePickedInt*shdr->sampleInt );
      }
      else if( vars->modeSearch == MODE_FIRST ) {
        // Search for first sample value that exceeds the threshold
        for( int i = vars->startSamp; i <= vars->endSamp; i++ ) {
          if( vars->thresh_sign*samples[i] > vars->threshold ) {
            samplePickedInt = i;
            break;
          }
        }
      }
      else if( vars->modeSearch == MODE_MAX ) {
        float maxVal  = -1e9;
        for( int isamp = vars->startSamp; isamp <= vars->endSamp; isamp++ ) {
          float val = vars->thresh_sign*samples[isamp];
          if( val > maxVal ) {
            maxVal     = val;
            samplePickedInt = isamp;
          }
          //      log->line("%d %f", isamp, val);
        }
        if( edef->isDebug() ) log->line(" Maximum index/value: %d, %f (sign=%.0f)", samplePickedInt, maxVal, vars->thresh_sign );
      }
    
      if( samplePickedInt == -1 ) {
        trace->getTraceHeader()->setFloatValue( vars->hdrId_timePick, 0 );
        if( vars->method == METHOD_PEAK_TROUGH ) trace->getTraceHeader()->setFloatValue( vars->hdrId_amplitude, 0 );
        return;
      }
      
      if( vars->method == METHOD_PEAK_TROUGH ) {
        // a) Make sure that actual min/max is picked, even if it lays outside of user specified time window
        float maxValue = samples[samplePickedInt] * vars->thresh_sign;    
        int sampIndex  = samplePickedInt + 1;
        while( sampIndex < shdr->numSamples && samples[sampIndex]*vars->thresh_sign > maxValue ) {
          maxValue        = samples[sampIndex] * vars->thresh_sign;
          samplePickedInt = sampIndex;
          sampIndex++;
        }
        // b) Determine actual peak/trough using quadratic interpolation
        float ampPicked = 0.0;
        float samplePicked = getQuadMaxSample( samples, samplePickedInt, shdr->numSamples, &ampPicked );
        float timePicked = samplePicked * shdr->sampleInt;
        if( edef->isDebug() ) log->line( "Picked sample/time/amplitude:  %f  %f  %f", samplePicked, timePicked, ampPicked );
        // c) If actual pick is outside user specified time window, force pick to be either start or end time:
        if( vars->forcePick ) {
          if( timePicked > vars->endTime ) {
            timePicked = vars->endTime;
            ampPicked = getQuadAmplitudeAtSample( samples, timePicked/shdr->sampleInt, shdr->numSamples );
            if( edef->isDebug() ) log->line( "       ...forced to:  %f  %f  %f", timePicked/shdr->sampleInt, timePicked, ampPicked );
          }
          else if( timePicked < vars->startTime ) {
            timePicked = vars->startTime;
            ampPicked = getQuadAmplitudeAtSample( samples, timePicked/shdr->sampleInt, shdr->numSamples );
            if( edef->isDebug() ) log->line( "       ...forced to:  %f  %f  %f", timePicked/shdr->sampleInt, timePicked, ampPicked );
          }
        }
        trace->getTraceHeader()->setFloatValue( vars->hdrId_timePick, timePicked );
        trace->getTraceHeader()->setFloatValue( vars->hdrId_amplitude, ampPicked );
      }
      else if( vars->method == METHOD_ZERO_CROSSING ) {
        float samplePicked = (float)samplePickedInt;
        float maxValue = samples[samplePickedInt];
        int sampIndex  = samplePickedInt;
        while( sampIndex > 0 && samples[sampIndex]*maxValue >= 0.0 ) {
          samplePickedInt = sampIndex;
          sampIndex--;
        }
        if( sampIndex == vars->endSamp ) {
          samplePicked = vars->endSamp;
        }
        else {
          //      if( fabs(samples[samplePickedInt]) > fabs(samples[samplePickedInt+1]) ) {
          //    samplePickedInt += 1;
          //  }
          samplePicked = getQuadZeroSample( samples, samplePickedInt, shdr->numSamples );
        }
        trace->getTraceHeader()->setFloatValue( vars->hdrId_timePick, samplePicked*shdr->sampleInt );
      }
    }
    //  }
  //**************************************************
  // Cross-correlation 'picking'
  //
  if( vars->isCross ) {
    int firstSample = vars->xcross_startSamp;
    int numSamplesWindow = vars->xcross_endSamp - vars->xcross_startSamp + 1;
    float* timePick = new float[ntraces];
    float* result   = new float[ntraces];
    float* s1;
    float* s2 = traceGather->trace(0)->getTraceSamples();
    //    traceGather->trace(0)->getTraceHeader()->setFloatValue( vars->hdrId_amplitude, 0 );
    //    traceGather->trace(0)->getTraceHeader()->setFloatValue( vars->hdrId_timePick, 0 );
    float accumulated = 0.0;
    float timePickCurrent = traceGather->trace(0)->getTraceHeader()->floatValue( vars->hdrId_timePick );
    timePick[0] = timePickCurrent;
    result[0]   = 0;
    for( int itrc = 1; itrc < ntraces; itrc++ ) {
      csTrace* trace = traceGather->trace(itrc);
      timePickCurrent = traceGather->trace(itrc)->getTraceHeader()->floatValue( vars->hdrId_timePick );
      timePick[itrc]  = timePickCurrent;

      s1 = s2;
      s2 = trace->getTraceSamples();
      result[itrc] = cross_correlation( s1, s2, firstSample, numSamplesWindow, vars->xcrossBuffer, shdr->numSamples );
      if( vars->outputCross ) {
        for( int isamp = 0; isamp < numSamplesWindow; isamp++ ) {
          s1[isamp] = vars->xcrossBuffer[isamp];
        }
        for( int isamp = numSamplesWindow; isamp < shdr->numSamples; isamp++ ) {
          s1[isamp] = 0.0;
        }
      }
      accumulated += result[itrc];
      //trace->getTraceHeader()->setFloatValue( vars->hdrId_amplitude,  );
      //trace->getTraceHeader()->setFloatValue( vars->hdrId_timePick, accumulated*shdr->sampleInt );
    }

    // Relative cross-correlation times are used to fill in mispicks.
    // First, we need to find one time pick that matches the cross-correlation time, and use this as a seed point to check all other picks.
    if( vars->xcross_method == XCROSS_ABSOLUTE ) {
      int seedTrace = 1;
      float* sorted = new float[ntraces];
      int* index = new int[ntraces];
      memcpy(sorted,timePick,ntraces*sizeof(float));
      for( int itrc = 0; itrc < ntraces; itrc++ ) {
        index[itrc] = itrc;
      }
      csSort<float> sort;
      sort.simpleSortIndex( sorted, ntraces, index );
      seedTrace = index[ntraces/2];
/*      for( int itrc = 1; itrc < ntraces; itrc++ ) {
        if( fabs(timePick[itrc] - timePick[itrc-1] - result[itrc]) <= vars->xcross_thresholdTime ) {
          seedTrace = itrc;
          break;
        }
      }
      seedTrace = ntraces/2; */
      if( edef->isDebug() ) log->line("PICKING: Seed trace: %d", seedTrace);
      for( int itrc = seedTrace; itrc < ntraces; itrc++ ) {
        if( fabs(timePick[itrc]-timePick[itrc-1]-result[itrc]) > vars->xcross_thresholdTime ) {
          timePick[itrc] = result[itrc]*shdr->sampleInt + timePick[itrc-1];
          csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
          trcHdr->setFloatValue( vars->hdrId_timePick, timePick[itrc] );
          if( vars->method == METHOD_PEAK_TROUGH ) {
            trcHdr->setFloatValue( vars->hdrId_amplitude, getQuadAmplitudeAtSample( s2, timePick[itrc]/shdr->sampleInt, shdr->numSamples ) );
          }
        }
      }
      for( int itrc = seedTrace-1; itrc > 0; itrc-- ) {
        if( fabs(timePick[itrc+1]-timePick[itrc]-result[itrc+1]) > vars->xcross_thresholdTime ) {
          timePick[itrc] = timePick[itrc+1] - result[itrc+1]*shdr->sampleInt;
          csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
          trcHdr->setFloatValue( vars->hdrId_timePick, timePick[itrc] );
          if( vars->method == METHOD_PEAK_TROUGH ) {
            trcHdr->setFloatValue( vars->hdrId_amplitude, getQuadAmplitudeAtSample( s2, timePick[itrc]/shdr->sampleInt, shdr->numSamples ) );
          }
        }
      }
      
    }
    else {  // RELATIVE: Just scan through, no need for seed point
      for( int itrc = 1; itrc < ntraces; itrc++ ) {
        if( fabs(timePick[itrc]-timePick[itrc-1]) > vars->xcross_thresholdTime ) {
          timePick[itrc] = result[itrc]*shdr->sampleInt + timePick[itrc-1];
          csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
          trcHdr->setFloatValue( vars->hdrId_timePick, timePick[itrc] );
          if( vars->method == METHOD_PEAK_TROUGH ) {
            trcHdr->setFloatValue( vars->hdrId_amplitude, getQuadAmplitudeAtSample( s2, timePick[itrc]/shdr->sampleInt, shdr->numSamples ) );
          }
        }
      }
    }
    delete [] timePick;
    delete [] result;
  }

}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_picking_( csParamDef* pdef ) {
  pdef->setModule( "PICKING", "Pick first breaks or other event" );

  pdef->setVersion( 1, 0 );

  pdef->addParam( "domain", "Time or sample domain", NUM_VALUES_FIXED );
  pdef->addValue( "time", VALTYPE_OPTION );
  pdef->addOption( "time", "Window is specified in time [ms]" );
  pdef->addOption( "sample", "Window is specified in samples (1 for first sample)" );

  pdef->addParam( "start", "Start time/sample", NUM_VALUES_FIXED, "Start time or sample, this depends on the 'domain' setting" );
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "Start time/sample of picking window" );

  pdef->addParam( "end", "End time/sample", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_HEADER_NUMBER, "End time/sample of picking window" );

  pdef->addParam( "method", "Picking method", NUM_VALUES_FIXED );
  pdef->addValue( "peak_trough", VALTYPE_OPTION );
  pdef->addOption( "peak_trough", "Pick a peak or trough. Search for first sample value that exceeds the specified +/- threshold." );
  pdef->addOption( "zero_crossing", "Pick last zero-crossing before specified peak/trough (STILL EXPERIMENTAL)" );
  pdef->addOption( "value", "Pick sample time at first occurence of specific value", "Specify value to be picked with user parameter 'thresh'" );

  pdef->addParam( "search_mode", "Search mode", NUM_VALUES_FIXED );
  pdef->addValue( "max", VALTYPE_OPTION );
  pdef->addOption( "max", "Pick maximum value in specified time window (honoring the specified polarity)." );
  pdef->addOption( "first", "Pick first value in specified time window that exceeds the given threshold (honoring the specified polarity)." );

  pdef->addParam( "thresh", "Threshold value. Peak/trough will only be searched when exceeding threshold", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Threshold value. Can be positive (pick peak) or negative (pick trough)" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Tolerance (+/-). Only used for method 'value'" );

  pdef->addParam( "header", "Header where picked time is stored", NUM_VALUES_VARIABLE );
  pdef->addValue( "time_pick", VALTYPE_STRING, "Name of trace header where picked time is stored" );
  pdef->addValue( "amp_pick", VALTYPE_STRING, "Name of trace header where picked amplitude is stored (only for method 'peak_trough')" );

  pdef->addParam( "force_pick", "Enforce user specified start/end window", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Picked value is forced to be inside the user specified window" );
  pdef->addOption( "no", "Picked value may be sligthly outside user specified window",
                   "This may happen if the nearest interpolated maximum (or zero crossing..) of the time pick falls outside of the specified window" );

  pdef->addParam( "xcross", "Guide first break picking by cross-correlation of adjacent traces", NUM_VALUES_FIXED,
                  "For this method, input data is read in in ensemble order. In each ensemble, consecutive traces are cross-correlated to guide the time picking" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Threshold time/sample", "Depends on 'domain' setting");
  pdef->addValue( "absolute", VALTYPE_OPTION );
  pdef->addOption( "absolute", "If the time pick difference between two consecutive traces differs from the time shift determined by cross-correlation by more than the specified threshold, use the cross-correlation value instead" );
  pdef->addOption( "relative", "If the time pick difference between two consecutive traces exceeds the specified threshold, use the cross-correlation value instead" );

  pdef->addParam( "xcross_win", "Window for cross-corrrelation", NUM_VALUES_FIXED,
                  "Consecutive traces are cross-correlated in the specified window, to find relative time shift" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Start time/sample of picking window", "Depends on the 'domain' setting. Set to 0 to use same window as for picker" );
  pdef->addValue( "0", VALTYPE_NUMBER, "End time/sample of picking window", "Depends on the 'domain' setting. Set to 0 to use same window as for picker" );

  pdef->addParam( "xcross_output", "Output cross-correlation function?", NUM_VALUES_FIXED);
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not output cross-correlation function");
  pdef->addOption( "yes", "Overwrite input data with cross-correlation function. Set unused samples to zero");

  //  pdef->addParam( "value", "Exact value which shall be picked (method 'value')", NUM_VALUES_VARIABLE );
  //  pdef->addValue( "", VALTYPE_NUMBER, "Value" );
}

extern "C" void _params_mod_picking_( csParamDef* pdef ) {
  params_mod_picking_( pdef );
}
extern "C" void _init_mod_picking_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_picking_( param, env, log );
}
extern "C" void _exec_mod_picking_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_picking_( traceGather, port, numTrcToKeep, env, log );
}

//--------------------------------------------------------------------------------
//
float cross_correlation( float const* s1, float const* s2, int firstSample, int numSamplesWindow, float* xcross, int numSamplesAll ) {

  //  float* xcross = new float[numSamplesWindow];
  int half = numSamplesWindow / 2;

  float xmax   = 0.0;
  int indexMax = 0;
  int lastSample = firstSample + numSamplesWindow - 1;

  int firstWin = 0;
  int lastWin  = numSamplesWindow - 1;
  int firstSample2 = firstSample - half;
  int lastSample2 = lastSample + numSamplesWindow - half;
  int numWindows = lastWin - firstWin + 1;

  if( firstSample2 >= 0 && lastSample2 < numSamplesAll ) {
    for( int iswin = 0; iswin < numSamplesWindow; iswin++ ) {
      float sum = 0.0;
      for( int isamp = firstSample; isamp <= lastSample; isamp++ ) {
        sum += s1[isamp] * s2[isamp-half+iswin];
      }
      xcross[iswin] = sum/(float)numSamplesWindow;
      if( xcross[iswin] > xmax ) {
        xmax = xcross[iswin];
        indexMax = iswin;
      }
    }
  }
  else {
    if( firstSample2 < 0 ) {
      firstWin = -firstSample2;
    }
    if( lastSample2 >= numSamplesAll ) {
      lastWin = numSamplesWindow - 1 - ( lastSample2 - numSamplesAll + 1 );
    }
    int numWindows = lastWin - firstWin + 1;
    if( numWindows < 3 ) return 0;
    for( int iswin = firstWin; iswin <= lastWin; iswin++ ) {
      float sum = 0.0;
      for( int isamp = firstSample; isamp <= lastSample; isamp++ ) {
        sum += s1[isamp] * s2[isamp-half+iswin];
      }
      xcross[iswin] = sum/(float)numWindows;
      if( xcross[iswin] > xmax ) {
        xmax = xcross[iswin];
        indexMax = iswin;
      }
    }
  }

  float value = getQuadMaxSample( &xcross[firstWin], indexMax, numWindows, &xmax );

  return( value - (float)half );
}


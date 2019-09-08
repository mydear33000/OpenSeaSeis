/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csTableAll.h"
#include "csSort.h"
#include "csTimer.h"
#include <cmath>
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SELECT_TIME
 *
 * @author Bjorn Olofsson
 * @date   2007
 *
 * Bug fixes
 *  2009-03-25 - Correctly zero out new samples
 */
namespace mod_select_time {
  struct VariableStruct {
    int numSamplesOrig;
    int startSamp;
    int endSamp;
    int mode;
    int hdrID_time_samp1_s;
    int hdrID_time_samp1_us;
    int hdrID_key;
    cseis_geolib::type_t hdrType_key;
    int numLocations;
    double* keyValues;
    int* numLines;
    double** startTimes;
    double** endTimes;
    bool isTable;
    bool isDelTrace;
    float percentDelTrace_ratio;
    bool doTrim;
    bool isShiftTrace;
  };
  static int const MODE_ABSOLUTE = 1;
  static int const MODE_RELATIVE = 2;
}
using namespace mod_select_time;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_select_time_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->isShiftTrace = false;
  vars->startSamp = 0;
  vars->endSamp   = 0;
  vars->mode      = MODE_RELATIVE;
  vars->hdrID_time_samp1_s = -1;
  vars->hdrID_time_samp1_us = -1;
  vars->hdrID_key       = -1;
  vars->hdrType_key     = TYPE_UNKNOWN;
  vars->startTimes      = NULL;
  vars->endTimes        = NULL;
  vars->numLines        = NULL;
  vars->numLocations    = 0;
  vars->isTable         = false;
  vars->isDelTrace      = false;
  vars->percentDelTrace_ratio = 1.0;
        vars->numSamplesOrig = shdr->numSamples;
  vars->doTrim = false;

  std::string text;
  //---------------------------------------------------------
  if( param->exists("mode") ) {
    param->getString( "mode", &text );
    if( !text.compare( "relative" ) ) {
      vars->mode = MODE_RELATIVE;
    }
    else if( !text.compare( "absolute" ) ) {
      vars->mode = MODE_ABSOLUTE;
    }
    else {
      log->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }
  //---------------------------------------------------------
  if( param->exists("free_mem") ) {
    param->getString( "free_mem", &text );
    if( !text.compare( "yes" ) ) {
      vars->doTrim = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->doTrim = false;
    }
    else {
      log->line("Option not recognized: '%s'.", text.c_str());
     env->addError();
    }
  }
  if( param->exists("shift") ) {
    param->getString( "shift", &text );
    if( !text.compare( "yes" ) ) {
      vars->isShiftTrace = true;
      if( vars->mode == MODE_ABSOLUTE ) {
        log->error("Shift trace cannot be applied for absolute time selection");
      }
    }
    else if( !text.compare( "no" ) ) {
      vars->isShiftTrace = false;
    }
    else {
      log->line("Option not recognized: '%s'.", text.c_str());
     env->addError();
    }
  }
  //---------------------------------------------------------
  bool isTimeDomain = true;

  if( param->exists("domain") ) {
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

  //---------------------------------------------------------
  if( param->exists("del_traces") ) {
    param->getString( "del_traces", &text );
    if( !text.compare( "yes" ) ) {
      vars->isDelTrace = true;
      if( param->getNumValues("del_traces") > 1 ) {
        float percentDelTrace;
        param->getFloat( "del_traces", &percentDelTrace, 1 );
        if( percentDelTrace <= 0 || percentDelTrace > 100 ) {
          log->line("Error in user parameter 'del_traces': Percent threshold of trace must be in the range of ]0,100]. Given value: %f", percentDelTrace);
          env->addError();
        }
        vars->percentDelTrace_ratio = percentDelTrace / 100.0;
      }
    }
    else if( !text.compare( "no" ) ) {
      vars->isDelTrace = false;
    }
    else {
      log->line("Option not recognized: '%s'.", text.c_str());
      env->addError();
    }
  }

  //---------------------------------------------------------
  if( param->exists("table") ) {
    if( vars->isShiftTrace ) {
      log->error("Trace shift is not possile to perform when specifying start/end times in a table");
    }

    csSort<double> sortObj;

    vars->isTable = true;
    param->getString( "table", &text );
    csTimer timer;
    timer.start();
    csTableAll table( csTableAll::TABLE_TYPE_DUPLICATE_KEYS );
    try {
      table.initialize( text );
      int numKeys   = table.numKeys();
      int numValues = table.numValues();
      if( numKeys != 1 ) log->error("Number of key columns in input table file: %d. Only 1 key column is currently supported", numKeys);
      if( numValues != 2 ) log->error("Number of value columns in input table file: %d. Expected 2 value columns, i.e. 'start' and 'end'", numValues);
      if( table.valueName(0).compare("start") ) {
        log->warning("Name of first table value is '%s'. Correct/expected value name should be 'start'.", table.valueName(0).c_str() );
      }
      if( table.valueName(1).compare("end") ) {
        log->warning("Name of second table value is '%s'. Correct/expected value name should be 'end'.", table.valueName(1).c_str() );
      }
      vars->hdrID_key   = hdef->headerIndex( table.keyName(0) );
      vars->hdrType_key = hdef->headerType( table.keyName(0) );

/*      int numCols = 3;
      type_t* colTypes = new type_t[numCols];
      colTypes[0] = vars->hdrType_key;
      colTypes[1] = TYPE_DOUBLE;
      colTypes[2] = TYPE_DOUBLE;
      table.readTableContents( colTypes, numKeys+numValues );
      delete [] colTypes;
*/
      table.readTableContents();

      if( edef->isDebug() ) table.dump();
      int id_start = 0;
      int id_end   = 1;
      vars->numLocations = table.numLocations();
      vars->keyValues  = new double[vars->numLocations];
      vars->startTimes = new double*[vars->numLocations];
      vars->endTimes   = new double*[vars->numLocations];
      vars->numLines   = new int[vars->numLocations];
      int* indexList   = NULL;
      int maxLines = 0;
      double sampleIntInSeconds = (double)shdr->sampleInt / 1000.0;
      for( int iloc = 0; iloc < vars->numLocations; iloc++ ) {
        vars->keyValues[iloc] = table.getKeyValue( iloc, 0 );
        TableValueList const* tvl = table.getValues( &vars->keyValues[iloc] );
        int numLines = tvl->numLines();
        if( numLines > maxLines ) {
          if( indexList ) delete [] indexList;
          indexList = new int[numLines];
          maxLines  = numLines;
        }
        vars->numLines[iloc] = numLines+1;
        double* startTimes = new double[numLines+1];
        double* endTimes   = new double[numLines+1];
        for( int iline = 0; iline < numLines; iline++ ) {
          startTimes[iline] = tvl->get(id_start,iline);
          indexList[iline]  = iline;
        }
        sortObj.simpleSortIndex( startTimes, numLines, indexList );

        for( int iline = 0; iline < numLines; iline++ ) {
          endTimes[iline] = tvl->get(id_end,indexList[iline]);
        }

        if( edef->isDebug() ) {
          log->line("Good ranges:");
          for( int iline = 0; iline < numLines; iline++ ) {
            log->line("  #%-3d:  %f  - %f", iline, startTimes[iline], endTimes[iline]);
          }
        }

        // NOW, convert GOOD selection times into BAD selection times, for easier selection later on
        for( int iline = numLines; iline > 0; iline-- ) {
          endTimes[iline]   = startTimes[iline] - sampleIntInSeconds;          
          startTimes[iline] = endTimes[iline-1] + sampleIntInSeconds;
        }
        endTimes[0]        = startTimes[0] - sampleIntInSeconds;
        startTimes[0]      = 0;
        endTimes[numLines] = 1e40;
        vars->startTimes[iloc] = startTimes;
        vars->endTimes[iloc]   = endTimes;

        if( edef->isDebug() ) {
          log->line("Bad ranges:");
          for( int iline = 0; iline < numLines+1; iline++ ) {
            log->line("  #%-3d:  %f  - %f", iline, startTimes[iline], endTimes[iline]);
          }
        }
      }
      if( indexList ) delete [] indexList;

      log->line("Input table '%s'", table.tableName().c_str() ); 
      log->line("  Number of key locations:       %d", vars->numLocations);
      log->line("  CPU time for reading in table: %12.6f seconds\n\n", timer.getElapsedTime() );
    }
    catch( csException& exc ) {
      log->error("Error when initializing input table '%s':\n", text.c_str(), exc.getMessage() );
    }
  }
  else { // No table given
    // No table given
    //---------------------------------------------------------
    float startTime = 0.0;
    float endTime = 0.0;
    
    if( isTimeDomain ) {
      if( vars->mode == MODE_ABSOLUTE ) {
        vars->numLocations = 1;
        vars->numLines = new int[vars->numLocations];
        vars->numLines[0] = param->getNumValues( "start" );
        if( vars->numLines[0] == 0 ) log->error("No start values given");
        else if( param->getNumValues( "end" ) != vars->numLines[0] ) log->error("Unequal start times (%d) and end times (%d) given", vars->numLines[0],param->getNumValues( "end" ));
        vars->startTimes = new double*[vars->numLocations];
        vars->endTimes   = new double*[vars->numLocations];
        vars->startTimes[0] = new double[vars->numLines[0]];
        vars->endTimes[0]   = new double[vars->numLines[0]];
        for( int iline = 0; iline < vars->numLines[0]; iline++ ) {
          double startTime;
          double endTime;
          param->getDouble( "start", &startTime, iline );
          param->getDouble( "end", &endTime, iline );
          vars->startTimes[0][iline] = startTime * 1000.0; // Convert to milliseconds
          vars->endTimes[0][iline]   = endTime * 1000.0; // Convert to milliseconds
          if( edef->isDebug() ) {
            log->line("Start/end times[%d]: %15.3fs  %15.3fs", iline, startTime, endTime);
          }
        }
      }
      else {
        param->getFloat( "start", &startTime );
        endTime = (float)(shdr->numSamples-1) * shdr->sampleInt;
        if( param->exists("end") ) {
          param->getFloat( "end", &endTime );
        }
        if( startTime < 0.0 ) log->error("Start time (%f) needs to be greater or equal to 0.0.", startTime);
        if( startTime > endTime ) log->error("Start time (%f) needs to be smaller than end time (%f).", startTime, endTime);
        if( startTime >= shdr->numSamples*shdr->sampleInt ) log->error("Start time (%f) exceeds input trace end time (%f).", startTime, shdr->numSamples*shdr->sampleInt);
        vars->startSamp = (int)(startTime / shdr->sampleInt);  // All in milliseconds
        vars->endSamp   = (int)(endTime / shdr->sampleInt);
      }
    }
    else { // Sample domain
      if( vars->mode == MODE_ABSOLUTE ) {
        log->error("Absolute start/end times must be provided in seconds since 1-1-1970, not in samples");
      }
      //
      // NOTE: User input is '1' for first sample. Internally, '0' is used!!
      //
      param->getInt( "start", &vars->startSamp );
      vars->endSamp = shdr->numSamples;
      if( param->exists("end") ) {
        param->getInt( "end", &vars->endSamp );
      }
      if( vars->startSamp < 1 ) log->error("Start sample (%d) needs to be greater or equal to 1.", vars->startSamp);
      if( vars->startSamp > vars->endSamp ) log->error("Start sample (%d) needs to be smaller than end sample (%d).", vars->startSamp, vars->endSamp);
      vars->startSamp -= 1;   // see note above..
      vars->endSamp   -= 1;
      startTime = (float)vars->startSamp * shdr->sampleInt;
      endTime   = (float)vars->endSamp * shdr->sampleInt;
    }
    //    if( vars->endSamp > shdr->numSamples ) log->error("Specified end sample/time (%d/%.2f) is larger than number of input samples (%d)",
    //                                                vars->endSamp, vars->endSamp*shdr->sampleInt, shdr->numSamples);

    // Set new number of samples
    if( vars->mode != MODE_ABSOLUTE ) {
      if( !vars->isShiftTrace ) {
        shdr->numSamples = vars->endSamp + 1;
      }
      else {
        shdr->numSamples = vars->endSamp - vars->startSamp + 1;
      }
    }
    else {
      vars->hdrID_time_samp1_s  = hdef->headerIndex( HDR_TIME_SAMP1.name );
      vars->hdrID_time_samp1_us = hdef->headerIndex( HDR_TIME_SAMP1_US.name );
    }

    if( edef->isDebug() ) {
      log->line("time1: %f, time2: %f, sample1: %d, sample2: %d, sampInt: %f\n", startTime, endTime, vars->startSamp, vars->endSamp, shdr->sampleInt );
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_select_time_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  //  csTraceHeaderDef const* hdef = env->headerDef;

  if( edef->isCleanup()){
    if( vars->startTimes != NULL ) {
      for( int iloc = 0; iloc < vars->numLocations; iloc++ ) {
        delete [] vars->startTimes[iloc];
        delete [] vars->endTimes[iloc];
      }
      delete [] vars->startTimes;
      delete [] vars->endTimes;
      vars->startTimes = NULL;
      vars->endTimes = NULL;
    }
    if( vars->numLines != NULL ) {
      delete [] vars->numLines;
      vars->numLines = NULL;
    }
    if( vars->keyValues != NULL ) {
      delete [] vars->keyValues;
      vars->keyValues  = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }
  
  if( edef->isDebug() ) {
    log->line("SELECT_TIME: Number of samples in trace: %d, Super header: %d",
      trace->numSamples(), shdr->numSamples );
  }

  float* samples = trace->getTraceSamples();

  if( vars->isTable ) {
    double sampleIntInSeconds = (double)shdr->sampleInt / 1000.0;

    csTraceHeader* trcHdr = trace->getTraceHeader();
    int time_samp1_s  = trcHdr->intValue( vars->hdrID_time_samp1_s );
    int time_samp1_us = trcHdr->intValue( vars->hdrID_time_samp1_us );

    double keyValue = trace->getTraceHeader()->doubleValue( vars->hdrID_key );
    int location = 0;    
    while( location < vars->numLocations && vars->keyValues[location] != keyValue ) {
      location += 1;
    }
    if( location >= vars->numLocations ) {
      log->warning("Key value %d not found in input table. No time selection applied.", keyValue);
    }
    else {
      double time_firstSamp = (double)time_samp1_s + ((double)time_samp1_us)/1000000.0;   // Absolute time of first sample, in seconds since 01-01-1970
      double time_lastSamp  = time_firstSamp + (double)((shdr->numSamples-1)*shdr->sampleInt)/1000.0;    // Absolute time of last sample

      int iline = 0;
      int numLines = vars->numLines[location];
      double* startTimes = vars->startTimes[location];
      double* endTimes   = vars->endTimes[location];

      int numDeletedSamples = 0;

      if( edef->isDebug() ) log->line("--- KEY %f ---", keyValue );
      while( iline < numLines && endTimes[iline] < time_firstSamp ) iline += 1;
      while( iline < numLines && startTimes[iline] <= time_lastSamp ) {
        if( startTimes[iline] <= time_firstSamp && endTimes[iline] >= time_lastSamp ) {
          if( edef->isDebug() ) log->line("   ZERO'ing WHOLE TRACE!");
          memset( samples, 0, shdr->numSamples*4 );
          numDeletedSamples = shdr->numSamples;
          break;
        }
        int startSamp = (int)( ( MAX( startTimes[iline], time_firstSamp ) - time_firstSamp ) / sampleIntInSeconds + 0.5 );
        int endSamp   = MIN( (int)( ( MIN( endTimes[iline], time_lastSamp ) - time_firstSamp ) / sampleIntInSeconds + 0.5 ), shdr->numSamples-1 );
        if( edef->isDebug() ) log->line("   ZERO'ing samples: %10d to %10d,   %f  to %f  (%f  %f)", startSamp, endSamp, startTimes[iline], endTimes[iline], time_firstSamp, time_lastSamp );
        memset( &samples[startSamp], 0, (endSamp-startSamp+1)*4 );
        numDeletedSamples += (endSamp-startSamp+1);
        iline += 1;
      }
      if( vars->isDelTrace && (float)numDeletedSamples/(float)shdr->numSamples >= vars->percentDelTrace_ratio ) {
        return false;
      }
      return true;
    }
  }
  else if( vars->mode == MODE_ABSOLUTE ) {
    csTraceHeader* trcHdr = trace->getTraceHeader();
    double startTimeCurrent_ms = 1000.0 * (double)trcHdr->intValue( vars->hdrID_time_samp1_s ) + (double)trcHdr->intValue( vars->hdrID_time_samp1_us )/1000.0;
    double endTimeCurrent_ms   = (double)startTimeCurrent_ms + (double)( shdr->numSamples * shdr->sampleInt );

    if( edef->isDebug() ) {
      log->line("Start/end times: %15.3fs  %15.3fs", startTimeCurrent_ms/1000, endTimeCurrent_ms/1000);
    }
    for( int iline = 0; iline < vars->numLines[0]; iline++ ) {
      if( endTimeCurrent_ms <= vars->startTimes[0][iline] || startTimeCurrent_ms >= vars->endTimes[0][iline] ) {
        continue;
      }
      else {
        return true;  // At least some portion of trace is within specified time window
        break;
      }
    }
    return false;
  }
  else if( !vars->isShiftTrace ) {
    //
    // Zero out trace samples outside the selected time gate
    //
    if( vars->startSamp > 0 ) {
      memset( samples, 0, vars->startSamp*4 );
    }
    if( vars->numSamplesOrig < shdr->numSamples ) {
      memset( &samples[vars->numSamplesOrig], 0, (shdr->numSamples - vars->numSamplesOrig)*4 );
    }
  }
  else {
    for( int isamp = vars->startSamp; isamp <= vars->endSamp; isamp++ ) {
      int indexTo = isamp - vars->startSamp;
      samples[indexTo] = samples[isamp];
    }
  }
  // Free any left-over memory
  if( vars->doTrim ) {
    trace->trim();
  }

  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_select_time_( csParamDef* pdef ) {
  pdef->setModule( "SELECT_TIME", "Select time of interest", "Samples outside of the chosen time interval are either set to zero, or removed" );

  pdef->addDoc("For active data processing, only the relative time is usually of interest. For continuous data processing, absolute time selection is usually required.");

  pdef->addParam( "domain", "Time or sample domain", NUM_VALUES_FIXED );
  pdef->addValue( "time", VALTYPE_OPTION );
  pdef->addOption( "time", "Window is specified in time [ms]" );
  pdef->addOption( "sample", "Window is specified in samples (1 for first sample)" );

  pdef->addParam( "mode", "Time selection mode", NUM_VALUES_FIXED );
  pdef->addValue( "relative", VALTYPE_OPTION );
  pdef->addOption( "relative", "Time window(s) are specified in relative time or sample index" );
  pdef->addOption( "absolute", "Time window(s) are specified as absolute times" );

  pdef->addParam( "start", "List of start times/samples", NUM_VALUES_FIXED, "Depends on 'domain' parameter" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Start times/samples" );

  pdef->addParam( "end", "List of end times/samples", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "End times/samples" );

  pdef->addParam( "table", "Table with window start/end times or samples", NUM_VALUES_FIXED, "For absolute time mode, columns 'start_time' and 'end_time', given in [ms] since 01-Jan-1970, must be specified" );
  pdef->addValue( "", VALTYPE_STRING, "Full path name of table containing window start/end times");

  pdef->addParam( "del_traces", "Delete trace if more then specified amount of trace has been de-selected", NUM_VALUES_VARIABLE,
                  "This option really only makes sense when used in conjunction with an absolute time selection" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not delete any traces. Set de-selected samples to zero." );
  pdef->addOption( "yes", "Delete traces if more than X percent of trace has been de-selected. Otherwise set de-selected samples to zero." );
  pdef->addValue( "100", VALTYPE_NUMBER, "Threshold of removed data that triggers trace deletion, given in percent [%].");

  pdef->addParam( "free_mem", "Free unused data in case output trace is shorter than input trace", NUM_VALUES_VARIABLE,
                  "This can be used to boost memory performance" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "By default, the system decides whether to free any unused memory immediately or reuse it for later." );
  pdef->addOption( "yes", "Free extra memory." );

  pdef->addParam( "shift", "Shift first life sample to start of trace?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Do not shift trace. This means that output trace length will be (endSample+1) samples long" );
  pdef->addOption( "yes", "Shift first life sample to start of trace. This means that output trace length will be (endSample-startSample+1) samples long" );
}

extern "C" void _params_mod_select_time_( csParamDef* pdef ) {
  params_mod_select_time_( pdef );
}
extern "C" void _init_mod_select_time_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_select_time_( param, env, log );
}
extern "C" bool _exec_mod_select_time_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_select_time_( trace, port, env, log );
}

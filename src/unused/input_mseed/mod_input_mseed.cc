/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csSort.h"
#include "MSeedReader.h"
#include "FlexHeader.h"
#include "csFileUtils.h"
#include "csStandardHeaders.h"

extern "C" {
  #include "libmseed.h"
}

using namespace cseis_system;
using namespace cseis_geolib;
using namespace mseed;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: input_mseed
 *
 * @author Bjorn Olofsson
 * @date   2008
 */
namespace mod_input_mseed {
  struct VariableStruct {
    char** filenames;
    int numFiles;
    int currentFile;
    MSeedReader* reader;

    float* sampleValues;
    int numAllocatedSampleValues;
    
    int endTimePreviousRecord_s;  // End time in [s] of previous record
    int endTimePreviousRecord_us;  // End time fraction [us] of previous record

    int nTracesToRead;
    int traceCounter;

    bool atEOF;
    int* hdrIndex;
    type_t* hdrType;
    int hdrID_trcno;
    int hdrID_fileno;

    int hdrID_year;
    int hdrID_day;
    int hdrID_hour;
    int hdrID_min;
    int hdrID_sec;
    int hdrID_msec;
    int hdrID_time_samp1;
    int hdrID_time_samp1_us;
    int hdrID_nsamp;
    int hdrID_network;
    int hdrID_station;
    int hdrID_location;
    int hdrID_quality;
    int hdrID_channel;
    int hdrID_type;
  };
}
using namespace mod_input_mseed;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_mseed_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->filenames    = NULL;
  vars->numFiles     = 0;
  vars->currentFile  = 0;
  vars->reader       = NULL;

  vars->nTracesToRead  = 0;
  vars->traceCounter = 0;

  vars->atEOF        = false;
  vars->hdrIndex     = NULL;
  vars->hdrType      = NULL;
  vars->hdrID_trcno  = -1;
  vars->hdrID_fileno = -1;

  vars->hdrID_year     = -1;
  vars->hdrID_day      = -1;
  vars->hdrID_hour     = -1;
  vars->hdrID_min      = -1;
  vars->hdrID_sec      = -1;
  vars->hdrID_msec     = -1;
  vars->hdrID_time_samp1= -1;
  vars->hdrID_time_samp1_us= -1;
  vars->hdrID_nsamp    = -1;
  vars->hdrID_network  = -1;
  vars->hdrID_station  = -1;
  vars->hdrID_location = -1;
  vars->hdrID_quality  = -1;
  vars->hdrID_channel  = -1;
  vars->hdrID_type     = -1;

  //----------------------------------------------------
  // Read in parameters
  //
  if( param->exists( "ntraces" ) ) {
    param->getInt( "ntraces", &vars->nTracesToRead );
    if( vars->nTracesToRead < 0 ) {
      vars->nTracesToRead = 0;
    }
  }
  //------------------------------------
  int numSingleFiles = param->getNumLines( "filename" );
  csVector<string> fileList;
  if( param->exists("directory") ) {
    int numValues = param->getNumValues("directory");
    string directory;
    param->getString("directory", &directory, 0);
//    string extension = "_Z.seed";  // !!CHANGE!!
    string extension = "seed";
    if( numValues > 1 ) {
      param->getString("directory", &extension, 1);
    }
    csFileUtils::retrieveFiles( directory, extension, &fileList, true, log->getFile() );

    // Sort file names
    string* names = new string[fileList.size()];
    for( int i = 0; i < fileList.size(); i++ ) {
      names[i] = fileList.at(i);
    }
    csSort<std::string>().simpleSort( names, fileList.size() );
    for( int i = 0; i < fileList.size(); i++ ) {
      fileList.set( names[i], i );
    }
    delete [] names;
  }

  if( numSingleFiles <= 0 && fileList.size() == 0 ) {
    log->error("No filename specified.");
  }
  vars->numFiles = numSingleFiles + fileList.size();
  vars->filenames = new char*[vars->numFiles];

  for( int i = numSingleFiles; i < vars->numFiles; i++ ) {
    string s1 = fileList.at(i-numSingleFiles);
    int length = s1.size();
    vars->filenames[i] = new char[length+1];
    memcpy(vars->filenames[i],s1.c_str(),length);
    vars->filenames[i][length] = '\0';
  }

  for( int i = 0; i < numSingleFiles; i++ ) {
    string filename;
    param->getStringAtLine( "filename", &filename, i );
    int length = filename.length();
    vars->filenames[i] = new char[length+1];
    memcpy(vars->filenames[i],filename.c_str(),length);
    vars->filenames[i][length] = '\0';
  }

  for( int i = 0; i < vars->numFiles; i++ ) {
    log->line("Input file #%4d: %s", i, vars->filenames[i]);
  }
  log->line("");

  //----------------------------------------------------
  int nsamples_preset = 360000;
  if( param->exists("nsamples") ) {
    param->getInt("nsamples", &nsamples_preset);
    if( nsamples_preset < 0 ) {
      log->error("Incorrect number of samples specified: %d.", nsamples_preset );
    }
  }

  //-------------------------------------------------------------------------
  // Read first Mini Seed file
  //
  vars->reader = new MSeedReader();

  vars->currentFile = -1;
  bool failed;
  do {
    failed = false;
    vars->currentFile += 1;
    try {
      vars->reader->read( vars->filenames[vars->currentFile], NULL );
      failed = (vars->reader->sampleRate() <= 0 || vars->reader->numSamples() <= 0);
    }
    catch( Exception& exc ) {
      log->warning("Error occurred when reading input file '%s'. System message:\n%s", vars->filenames[vars->currentFile], exc.getMessage() );
      failed = true;
    }
    // Read until non-empty file has been found:
  } while( failed && vars->currentFile < vars->numFiles );

  if( failed ) {
    log->error("Failed to read any given input file.");
  }

  if( nsamples_preset == 0 ) {
    shdr->numSamples = vars->reader->numSamples();
  }
  else {
    if( nsamples_preset < vars->reader->numSamples() ) {
      log->warning("Number of samples found in first input file (%d) exceeds fixed number of samples (%d) (user parameter 'nsamples')",
        vars->reader->numSamples(), nsamples_preset );
    }
    shdr->numSamples = nsamples_preset;
  }
  if( vars->reader->sampleRate() <= 0 || shdr->numSamples <= 0 ) {
    log->error("Incorrect sample rate (=%f) or number of samples (=%d).", vars->reader->sampleRate(), shdr->numSamples);
  }

  shdr->sampleInt  = 1000.0 / vars->reader->sampleRate();

  //---------------------------------------------------------------------
  // Set trace headers
  //
  vars->hdrID_year  = hdef->addStandardHeader( HDR_TIME_YEAR.name );
  vars->hdrID_day   = hdef->addStandardHeader( HDR_TIME_DAY.name );
  vars->hdrID_hour  = hdef->addStandardHeader( HDR_TIME_HOUR.name );
  vars->hdrID_min   = hdef->addStandardHeader( HDR_TIME_MIN.name );
  vars->hdrID_sec   = hdef->addStandardHeader( HDR_TIME_SEC.name );
  vars->hdrID_msec  = hdef->addStandardHeader( HDR_TIME_MSEC.name );
  vars->hdrID_time_samp1    = hdef->addStandardHeader( HDR_TIME_SAMP1.name );
  vars->hdrID_time_samp1_us = hdef->addStandardHeader( HDR_TIME_SAMP1_US.name );
  vars->hdrID_nsamp         = hdef->addStandardHeader( HDR_NSAMP.name );

  vars->hdrID_network   = hdef->addHeader( TYPE_STRING, "network", "MSeed Network", 11 );
  vars->hdrID_station   = hdef->addHeader( TYPE_STRING, "station", "MSeed Station", 11 );
  vars->hdrID_location  = hdef->addHeader( TYPE_STRING, "location", "MSeed Location", 11 );
  vars->hdrID_channel   = hdef->addHeader( TYPE_STRING, "channel", "MSeed Channel", 11 );
  vars->hdrID_quality   = hdef->addHeader( TYPE_STRING, "quality", "MSeed Quality", 2 );
  vars->hdrID_type      = hdef->addHeader( TYPE_STRING, "sampletype", "MSeed Sample Type", 2 );

  vars->hdrID_trcno     = hdef->addStandardHeader( HDR_TRCNO.name );
  vars->hdrID_fileno    = hdef->addStandardHeader( HDR_FILENO.name );

  log->line("First input file:  %s", vars->filenames[0] );
  log->line("Sample interval:   %fms", shdr->sampleInt );
  log->line("Number of samples: %d", shdr->numSamples );
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_input_mseed_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const*  shdr = env->superHeader;
//  csTraceHeaderDef const*  hdef = env->headerDef;

  if( edef->isCleanup() ) {
    if( vars->reader != NULL ) {
      delete vars->reader;
      vars->reader = NULL;
    }
    if( vars->filenames != NULL ) {
      for( int i = 0; i < vars->numFiles; i++ ) {
        if( vars->filenames[i] != NULL ) delete [] vars->filenames[i];
      }
      delete [] vars->filenames;
      vars->filenames = NULL;
    }
    if( vars->hdrIndex != NULL ) {
      delete [] vars->hdrIndex;
      vars->hdrIndex = NULL;
    }
    if( vars->hdrType != NULL ) {
      delete [] vars->hdrType;
      vars->hdrType = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  if( vars->atEOF ) return false;

  if( vars->nTracesToRead > 0 && vars->nTracesToRead == vars->traceCounter ) {
    vars->atEOF = true;
    return false;
  }

  //-----------------------------------------------------------------------------------
  //
  MSeedReader* reader   = vars->reader;
  csTraceHeader* trcHdr = trace->getTraceHeader();
  float* samples        = trace->getTraceSamples();

  int numSamples_mseed = reader->numSamples();
  if( edef->isDebug() ) {
    fprintf(stderr,"Reading MSeed file #%3d (ns %d): %s\n", vars->traceCounter+1, numSamples_mseed, vars->filenames[vars->currentFile] );
    log->line("");
  }
  if( numSamples_mseed > shdr->numSamples ) {
    log->error("Unsupported inconsistency in number of samples in MiniSeed file, Actual > Expected: %d > %d", numSamples_mseed, shdr->numSamples);
  }

  reader->setData( samples, numSamples_mseed );
  for( int isamp = numSamples_mseed; isamp < shdr->numSamples; isamp++ ) {
    samples[isamp] = 0.0;
  }

  // Set headers
  trcHdr->setIntValue( vars->hdrID_year,   reader->headerValue(mseed::HDR_ID_YEAR)->intValue() );
  trcHdr->setIntValue( vars->hdrID_day,    reader->headerValue(mseed::HDR_ID_DAY)->intValue() );
  trcHdr->setIntValue( vars->hdrID_hour,   reader->headerValue(mseed::HDR_ID_HOUR)->intValue() );
  trcHdr->setIntValue( vars->hdrID_min,    reader->headerValue(mseed::HDR_ID_MIN)->intValue() );
  trcHdr->setIntValue( vars->hdrID_sec,    reader->headerValue(mseed::HDR_ID_SEC)->intValue() );
  trcHdr->setIntValue( vars->hdrID_msec,   reader->headerValue(mseed::HDR_ID_FRACT)->intValue() );

  trcHdr->setStringValue( vars->hdrID_network, reader->headerValue(mseed::HDR_ID_NETWORK)->stringValue() );
  trcHdr->setStringValue( vars->hdrID_location, reader->headerValue(mseed::HDR_ID_LOCATION)->stringValue() );
  trcHdr->setStringValue( vars->hdrID_station, reader->headerValue(mseed::HDR_ID_STATION)->stringValue() );
  trcHdr->setStringValue( vars->hdrID_channel, reader->headerValue(mseed::HDR_ID_CHANNEL)->stringValue() );
  trcHdr->setStringValue( vars->hdrID_quality, reader->headerValue(mseed::HDR_ID_QUALITY)->stringValue() );
  trcHdr->setStringValue( vars->hdrID_type, reader->headerValue(mseed::HDR_ID_TYPE)->stringValue() );

  int time_samp1_s  = (int)(reader->headerValue(mseed::HDR_ID_STARTTIME)->int64Value()/1000000);
  int time_samp1_us = (int)(reader->headerValue(mseed::HDR_ID_STARTTIME)->int64Value()%1000000);
  int endtime_s  = (int)(reader->headerValue(mseed::HDR_ID_ENDTIME)->int64Value()/1000000);
  int endtime_us = (int)(reader->headerValue(mseed::HDR_ID_ENDTIME)->int64Value()%1000000);
  int nsamp = (int)( ( (endtime_s-time_samp1_s)*1000 + (endtime_us-time_samp1_us)/1000 ) / shdr->sampleInt + 0.5);

  trcHdr->setIntValue( vars->hdrID_time_samp1,    time_samp1_s );
  trcHdr->setIntValue( vars->hdrID_time_samp1_us, time_samp1_us);
  trcHdr->setIntValue( vars->hdrID_nsamp,         nsamp );

  if( nsamp != numSamples_mseed ) {
    log->warning("Inconsistent number of samples in MSeed file '%s'.\n Actual number of data samples:  %d\n ...computed from MSeed header:  %d\n",
                 vars->filenames[vars->currentFile], numSamples_mseed, nsamp );
    trcHdr->setIntValue( vars->hdrID_nsamp,         numSamples_mseed );
  }
  if( edef->isDebug() ) {
    log->line("time_samp1 = %10ds (%dus), nsamp = %10d, numSamples_mseed = %10d (shdr: %d)",
	      time_samp1_s, time_samp1_us, nsamp, numSamples_mseed, shdr->numSamples);
  }

  vars->endTimePreviousRecord_s  = endtime_s;
  vars->endTimePreviousRecord_us = endtime_us;

  trcHdr->setIntValue( vars->hdrID_trcno, vars->traceCounter+1 );
  trcHdr->setIntValue( vars->hdrID_fileno, vars->currentFile+1 );

  //------------------------------------------------------------------------------------
  // Read in next file
  //
  vars->currentFile  += 1;
  vars->traceCounter += 1;
  if( vars->currentFile >= vars->numFiles ) {
    vars->atEOF = true;
    return true;
  }

  try {
    do { // Loop until non-empty file has been read in, or last file has been reached
      reader->read( vars->filenames[vars->currentFile], NULL );
      numSamples_mseed = reader->numSamples();
      if( numSamples_mseed <= 0 ) {
        if( edef->isDebug() ) {
          log->line(" Input file is empty: '%s'", vars->filenames[vars->currentFile] );
          fprintf(stderr," Input file is empty: '%s'\n", vars->filenames[vars->currentFile] );
        }
        vars->currentFile += 1;
        if( vars->currentFile >= vars->numFiles ) {
          vars->atEOF = true;
          return true;
        }
      }
    } while( numSamples_mseed <= 0 );
    if( numSamples_mseed > shdr->numSamples ) {
      log->error("Number of samples in input file '%s' (%d) exceeds number of samples in trace (%d)",
                 vars->filenames[vars->currentFile], numSamples_mseed, shdr->numSamples );
      //vars->atEOF = true;
      return true;
    }

    if( edef->isDebug() && numSamples_mseed < shdr->numSamples ) {
      log->warning("Fewer samples than expected in file '%s' (trace #%d):\n  Number of samples = %d < %d. Trace will be padded with zeros.",
                   vars->filenames[vars->currentFile], vars->traceCounter+1, numSamples_mseed, shdr->numSamples );
      
    }

    if( edef->isDebug() ) {
      int startTimeNextRecord_s = (int)(reader->headerValue(mseed::HDR_ID_STARTTIME)->int64Value()/1000000);
      int timeDiff = 1000*(startTimeNextRecord_s - vars->endTimePreviousRecord_s);
      log->line("Time difference: %d", timeDiff );
    }
  }
  catch( Exception& exc ) {
    log->warning("Error occurred when reading input file '%s'. System message:\n%s", vars->filenames[vars->currentFile], exc.getMessage() );
    vars->atEOF = true;
    return true;
  }

  //------------------------------------------------------------------------------------
  return true;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_mseed_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_MSEED", "Input Mini SEED file", "Reads data recorded in Mini SEED format" );

  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "directory", "Name of directory to search", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Name of directory", "All files in this directory and its subdirectories will be read" );
  pdef->addValue( "seed", VALTYPE_STRING, "File extension", "Only files with the given file extension will be read in" );

  pdef->addParam( "nsamples", "Manually set number of output samples", NUM_VALUES_FIXED );
  pdef->addValue( "360000", VALTYPE_NUMBER, "Number of output samples", "Specify 0 to set automatically. Program will abort if actual number of samples in input file is larger than specified number of samples" );

  pdef->addParam( "ntraces", "Number of traces to read in", NUM_VALUES_FIXED,
                  "Input of traces will stop when all traces have been read in, or if the number of traces specified has been reached. Traces will not be filled up to the specified range");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to read in" );
}

extern "C" void _params_mod_input_mseed_( csParamDef* pdef ) {
  params_mod_input_mseed_( pdef );
}
extern "C" void _init_mod_input_mseed_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_input_mseed_( param, env, log );
}
extern "C" bool _exec_mod_input_mseed_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_input_mseed_( trace, port, env, log );
}



/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "MSeedReader.h"
#include "FlexHeader.h"
#include <stdarg.h>
#include <cstring>
#include <exception>
#include <cmath>

using namespace mseed;
using namespace std;

static int MINIMUM( int val1, int val2 ) {
  return( val1 < val2 ? val1 : val2 );
}
static int MAXIMUM( int val1, int val2 ) {
  return( val1 > val2 ? val1 : val2 );
}

MSeedReader::MSeedReader() {
  myDataChunks   = new vector<DataChunk*>();
  myData         = NULL;
  myHeaderValues = new FlexHeader[NUM_MSEED_HEADERS];
  myNumStoredDataSamples    = 0;
  myNumAllocatedDataSamples = 0;
  myNumDataChunkSamples     = 0;
  mySampleRate = 0;
}

MSeedReader::~MSeedReader() {
  if( myData != NULL ) {
    delete [] myData;
    myData = NULL;
  }
  if( myDataChunks != NULL ) {
    freeMemory();
    delete myDataChunks;
    myDataChunks = NULL;
  }
  if( myHeaderValues != NULL ) {
    delete [] myHeaderValues;
    myHeaderValues = NULL;
  }
}

//--------------------------------------------------
void MSeedReader::freeMemory() {
  for( int i = 0; i < (int)myDataChunks->size(); i++ ) {
    DataChunk* chunk = myDataChunks->at(i);
    delete [] chunk->data;
    delete chunk;
  }
  myDataChunks->clear();
  myNumDataChunkSamples = 0;
}
//--------------------------------------------------
//
//
int MSeedReader::numSamples() const {
  return( myNumStoredDataSamples+myNumDataChunkSamples );
}
//--------------------------------------------------
//
//
FlexHeader const* MSeedReader::headerValue( int headerIndex ) const {
  return &myHeaderValues[headerIndex];
}
//--------------------------------------------------
//
//
void MSeedReader::setData( float* samples, int numSamplesToSet, int firstSampleIndex ) const {
  if( firstSampleIndex+numSamplesToSet > numSamples() ) {
    throw( Exception("MSeedReader::setData: Incorrect number of data samples: %d+%d > %d", firstSampleIndex, numSamplesToSet, numSamples()) );
  }
  int counterSamplesSet = 0;
  if( myNumStoredDataSamples > firstSampleIndex ) { // myData contains data samples --> copy these
    counterSamplesSet = MINIMUM( numSamplesToSet, myNumStoredDataSamples-firstSampleIndex );
    memcpy( samples, &myData[firstSampleIndex], counterSamplesSet*4 );
  }
  int counterNumSamples = myNumStoredDataSamples;

  // Difficult part: Copy correct values from data chunk objects to output 'samples' float array
  int ichunk = 0;
  //  while( ichunk < myDataChunks->size() && counterNumSamples+myDataChunks->at(ichunk).numSamples < firstSampleIndex ) {
  //  counterNumSamples += myDataChunks->at(ichunk).numSamples;
  //  ichunk += 1;
  // }
  while( counterSamplesSet < numSamplesToSet && ichunk < (int)myDataChunks->size() ) {
    DataChunk* chunk = myDataChunks->at(ichunk++);
    if( counterNumSamples+chunk->numSamples < firstSampleIndex ) {
      counterNumSamples += chunk->numSamples;
      continue;  // First sample index to copy has not been reached yet
    }
    int firstChunkIndex = MAXIMUM( firstSampleIndex - counterNumSamples, 0 );   // Index pointer into first value in data chunk to copy
    int numSamplesToCopy = MINIMUM( numSamplesToSet-counterSamplesSet, chunk->numSamples-firstChunkIndex );
    memcpy( &samples[counterSamplesSet], &chunk->data[firstChunkIndex], numSamplesToCopy*4 );
    counterSamplesSet += numSamplesToCopy;
    counterNumSamples += chunk->numSamples;
  }
}
//--------------------------------------------------
//
//
void MSeedReader::read( char* filename, FILE* debugStream ) {
  // If previous call to this method had data stored as data chunks:
  //   Try to (re-)allocate myData array so that it will store all data samples of the given file
  if( myDataChunks->size() != 0 ) {
    //    int numSampleCounter = 0;
    // for( int i = 0; i < myDataChunks->size(); i++ ) {
    //   DataChunk* chunk = myDataChunks->at(i);
    //  numSampleCounter += chunk->numSamples;
    // }
    if( myNumStoredDataSamples+myNumDataChunkSamples > myNumAllocatedDataSamples ) {
      if( myData != NULL ) {
        delete [] myData;
        myData = NULL;
      }
      myNumAllocatedDataSamples = myNumStoredDataSamples+myNumDataChunkSamples;
      myData = new float[myNumAllocatedDataSamples];
    }
    freeMemory();
  }
  myNumStoredDataSamples = 0;
  myNumDataChunkSamples  = 0;

  MSRecord* ms_record = msr_init( NULL );
  int numTotalSampleCounter = 0;  // Total number of samples read in so far
  int recordCounter         = 0;  // Counter of MS records read in so far
  int returnCode            = MS_NOERROR;

  int reclen = -1;
  off_t fpos = 0;
  int last   = 0;
  flag skipnotdata = 1; // Skip non-valid data? 0: no, 1: yes
  flag dataflag    = 1; // 1: True, unpack data correctly (?)
  flag verbose     = 0; // 1: Verbose; >1 : Very verbose
  if( debugStream != NULL ) verbose = 2;

  hptime_t endtime = 0;

  do {
    returnCode = ms_readmsr( &ms_record, filename, reclen, &fpos, &last, skipnotdata, dataflag, verbose );
    if( returnCode == MS_ENDOFFILE ) break;
    if( returnCode != MS_NOERROR ) {
      string message = "";
      switch( returnCode ) {
      case MS_GENERROR:
        message = "Generic unspecified error";
        break;
      case MS_NOTSEED:
        message = "Data is not SEED";
        break;
      case MS_WRONGLENGTH:
        message = "Length of data read was not correct";
        break;
      case MS_OUTOFRANGE:
        message = "SEED record length out of range";
        break;
      case MS_UNKNOWNFORMAT:
        message = "Unknown data encoding format";
        break;
      case MS_STBADCOMPFLAG:
        message = "Steim, invalid compression flag(s)";
        break;
      case MS_ENDOFFILE:
        message = "End of file reached";
        break;
      default:
        message = "Unknown error";
        break;
      }
      throw( Exception("Error when opening MiniSeed file %s. Returned error code: %d.\nMessage: %s", filename, returnCode, message.c_str()) );
    }
    recordCounter += 1;

    if( recordCounter == 1 ) { // Set trace headers
      BTime btime;
      ms_hptime2btime( ms_record->starttime, &btime );
      myHeaderValues[HDR_ID_YEAR] .setIntValue( (int)btime.year );
      myHeaderValues[HDR_ID_DAY]  .setIntValue( (int)btime.day );
      myHeaderValues[HDR_ID_HOUR] .setIntValue( (int)btime.hour );
      myHeaderValues[HDR_ID_MIN]  .setIntValue( (int)btime.min );
      myHeaderValues[HDR_ID_SEC]  .setIntValue( (int)btime.sec );
      myHeaderValues[HDR_ID_FRACT].setIntValue( (int)btime.fract );
      myHeaderValues[HDR_ID_STARTTIME].setInt64Value( (long long)ms_record->starttime );
      std::string text( ms_record->network );
      myHeaderValues[HDR_ID_NETWORK]  .setStringValue( text );
      text = ms_record->station;
      myHeaderValues[HDR_ID_STATION]  .setStringValue( text );
      text = ms_record->location;
      myHeaderValues[HDR_ID_LOCATION] .setStringValue( text );
      text = ms_record->channel;
      myHeaderValues[HDR_ID_CHANNEL]  .setStringValue( text );
      myHeaderValues[HDR_ID_QUALITY]  .setCharValue( ms_record->dataquality );
      myHeaderValues[HDR_ID_TYPE]     .setCharValue( ms_record->sampletype );
      myHeaderValues[HDR_ID_SAMPLERATE].setDoubleValue( ms_record->samprate );
      mySampleRate = ms_record->samprate;
    }
    endtime = ms_record->starttime + (hptime_t)ms_record->samplecnt*(hptime_t)(1000000.0/ms_record->samprate);

    float* floatDataPtr = NULL;
    // myData array is allocated and still has enough space
    if( (myNumAllocatedDataSamples != 0) && (numTotalSampleCounter <= myNumAllocatedDataSamples) &&
        ((numTotalSampleCounter+ms_record->numsamples) <= myNumAllocatedDataSamples) ) {
      floatDataPtr = &myData[numTotalSampleCounter];
      myNumStoredDataSamples += ms_record->numsamples;
    }
    // myData array does not have enough space (or hasn;t been allocated at all) --> Save data in data chunks
    else {
      DataChunk* chunk  = new DataChunk();
      chunk->numSamples = ms_record->numsamples;
      chunk->data       = new float[chunk->numSamples];
      floatDataPtr      = chunk->data;
      myDataChunks->push_back(chunk);
      myNumDataChunkSamples += chunk->numSamples;
    }
    numTotalSampleCounter += ms_record->numsamples;

    // Copy data elements from MS Record to float* array, pointed to by floatDataPtr
    // Set data sample values
    if( ms_record->sampletype == 'f' ){
      memcpy( floatDataPtr, ms_record->datasamples, ms_record->samplecnt*4 );
    }
    else if( ms_record->sampletype == 'i' ) {
      int32_t* dataPtr = (int32_t*)ms_record->datasamples;
      for( int isamp = 0; isamp < ms_record->numsamples; isamp++ ) {
        floatDataPtr[isamp] = (float)dataPtr[isamp];
      }
    }
    else if( ms_record->sampletype == 'd' ) {
      double* dataPtr = (double*)ms_record->datasamples;
      for( int isamp = 0; isamp < ms_record->numsamples; isamp++ ) {
        floatDataPtr[isamp] = (float)dataPtr[isamp];
      }
    }
    else {
      throw( Exception("Error, unrecognized sample type: '%c'\n", ms_record->sampletype) );
    }


    //----------------------------------------------------
    // Printout MS record contents
    if( debugStream != NULL ) {
      MSRecord* mr = ms_record;
      fprintf(debugStream,"\n------ MS record #%4d ----------------------\n", recordCounter);
      fprintf(debugStream,"Record length:          %d\n", mr->reclen );
      fprintf(debugStream,"Sequence number:        %d\n", mr->sequence_number );
      fprintf(debugStream,"Network:                '%s'\n", mr->network );
      fprintf(debugStream,"Station:                '%s'\n", mr->station );
      fprintf(debugStream,"Location:               '%s'\n", mr->location );
      fprintf(debugStream,"Channel:                '%s'\n", mr->channel );
      fprintf(debugStream,"Data quality:           '%c'\n", mr->dataquality );
      fprintf(debugStream,"Start time [us]:        %ld\n", (long int)mr->starttime );
      fprintf(debugStream,"Sample rate [Hz]:       %f\n", mr->samprate );
      fprintf(debugStream,"No. samples per record: %d\n", mr->samplecnt );
      fprintf(debugStream,"Encoding:               %d\n", mr->encoding );
      fprintf(debugStream,"Byte order:             %d\n", mr->byteorder );
      fprintf(debugStream,"Number of samples:      %d\n", mr->numsamples );
      fprintf(debugStream,"Sample type:            '%c'\n", mr->sampletype );

      int flag_details = 2;
      msr_print( ms_record, flag_details );
      fprintf(debugStream,"------ Other ----------------------\n");
      double sampleRate  = msr_samprate( ms_record );
      hptime_t startTime = msr_starttime( ms_record );
      int sampleByteSize = (int)get_samplesize( ms_record->sampletype );
      fprintf(debugStream,"OTHER: Sample rate: %f, start time: %ld, sample byte size: %d\n", sampleRate, (long int)startTime, sampleByteSize );
    } // isDebug
  } while( returnCode == MS_NOERROR );

  myHeaderValues[HDR_ID_ENDTIME].setInt64Value( (long long)endtime );
}
//--------------------------------------------------
//
//
void MSeedReader::dumpDataValues( std::FILE* stream ) {
  int nSamples = numSamples();
  DataChunk* cPtr = NULL;
  int counterChunkSamples = 0;
  int currentChunk = 0;
  float value = 0.0;
  for( int isamp = 0; isamp < nSamples; isamp++ ) {
    if( isamp < myNumStoredDataSamples ) {
      value = myData[isamp];
    }
    else {
      if( cPtr == NULL || counterChunkSamples == cPtr->numSamples ) {
        cPtr = myDataChunks->at(currentChunk);
        currentChunk += 1;
        counterChunkSamples = 0;
      }
      value = cPtr->data[counterChunkSamples];
      counterChunkSamples += 1;
    }
    fprintf(stream,"%10d %f\n", isamp, value );
  }
}
//--------------------------------------------------
//
//
void MSeedReader::dumpHeaderValues( std::FILE* stream ) {
  fprintf( stream, "Year:            %d\n", myHeaderValues[HDR_ID_YEAR].intValue() );
  fprintf( stream, "Day:             %d\n", myHeaderValues[HDR_ID_DAY].intValue() );
  fprintf( stream, "Hour:            %d\n", myHeaderValues[HDR_ID_HOUR].intValue() );
  fprintf( stream, "Minute:          %d\n", myHeaderValues[HDR_ID_MIN].intValue() );
  fprintf( stream, "Second:          %d\n", myHeaderValues[HDR_ID_SEC].intValue() );
  fprintf( stream, "Fraction:        %d\n", myHeaderValues[HDR_ID_FRACT].intValue() );
  fprintf( stream, "Start time[us]:  %lld\n", myHeaderValues[HDR_ID_STARTTIME].int64Value() );
  fprintf( stream, "End time[us]:    %lld\n", myHeaderValues[HDR_ID_ENDTIME].int64Value() );
  fprintf( stream, "Network:         '%s'\n", myHeaderValues[HDR_ID_NETWORK].stringValue() );
  fprintf( stream, "Station:         '%s'\n", myHeaderValues[HDR_ID_STATION].stringValue() );
  fprintf( stream, "Location:        '%s'\n", myHeaderValues[HDR_ID_LOCATION].stringValue() );
  fprintf( stream, "Channel:         '%s'\n", myHeaderValues[HDR_ID_CHANNEL].stringValue() );
  fprintf( stream, "Quality:         '%c'\n", myHeaderValues[HDR_ID_QUALITY].charValue() );
  fprintf( stream, "Sample type:     '%c'\n", myHeaderValues[HDR_ID_TYPE].charValue() );
  fprintf( stream, "Sample rate[Hz]: %f\n", myHeaderValues[HDR_ID_SAMPLERATE].doubleValue() );
}

//--------------------------------------------------
// Exception class
//
Exception::Exception() { myMessage = ""; }
Exception::Exception( std::string const& message ) { myMessage = message; }
Exception::Exception( int dummy, char const* text ) {
  myMessage = text;
}
Exception::Exception( char const* text, ... ) {
  char tmpString[1000];
  va_list argList;
  va_start( argList, text );
  vsprintf( tmpString, text, argList );
  myMessage = tmpString;
}
Exception::Exception( Exception const& obj ) {
  myMessage = obj.myMessage;
}
Exception::~Exception() {}
const char* Exception::getMessage() { return myMessage.c_str(); }


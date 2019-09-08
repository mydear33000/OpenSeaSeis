/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef MSEED_READER_H
#define MSEED_READER_H

#include <string>
#include <vector>
#include <cstdio>

extern "C" {
#include "libmseed.h"
}

namespace mseed {

static int const HDR_ID_YEAR  = 0;
static int const HDR_ID_DAY   = 1;
static int const HDR_ID_HOUR  = 2;
static int const HDR_ID_MIN   = 3;
static int const HDR_ID_SEC   = 4;
static int const HDR_ID_FRACT = 5;
static int const HDR_ID_STARTTIME = 6;
static int const HDR_ID_ENDTIME   = 7;
static int const HDR_ID_NETWORK   = 8;
static int const HDR_ID_STATION   = 9;
static int const HDR_ID_LOCATION  = 10;
static int const HDR_ID_CHANNEL   = 11;
static int const HDR_ID_QUALITY   = 12;
static int const HDR_ID_TYPE      = 13;
static int const HDR_ID_SAMPLERATE= 14;
static int const NUM_MSEED_HEADERS = 15;

struct DataChunk {
  int    numSamples;
  float* data;
};
 class FlexHeader;

/**
 * Mini SEED data reader
 * Usage:
 *  1) Create MSeedReader object
 *  2) Call method 'read()' once, for one MSeed file
 *  3) Allocate float* array of size numSamples() (number of samples that were read in)
 *  4) Call method 'setData()', passing the newly allocated float* array, to set all data sample values
 *  5) To read more files, go back to point 2)
 *
 * @author Bjorn Olofsson
 * @date   2008
 */
class MSeedReader {
 public:
  MSeedReader();
  ~MSeedReader();
  /**
   * Read single MSeed data file.
   * Throws Exception if error occurs.
   * @param filename    (i)  MSeed file name
   * @param debugStream (io) File stream. Pass NULL if no debug information shall be output
   */
  void read( char* filename, std::FILE* debugStream );
  /**
   * @return total number of data samples currently stored in MSeedReader object
   */
  int numSamples() const;
  /**
   * @return Sample rate [Hz]
   */
  double sampleRate() const { return mySampleRate; }
  /**
   * Set data samples
   * @param samples          (o) Data sample array which values shall be set
   * @param numSamplesToSet  (i) Number of samples to set in input array.
   * @param firstSampleIndex (i) Index of first sample to set (Exception is thrown if firstSampleIndex+numSamplesToSet > numSamples())
   */
  void setData( float* samples, int numSamplesToSet, int firstSampleIndex = 0 ) const;
  int numHeaderValues() const { return NUM_MSEED_HEADERS; }
  FlexHeader const* headerValue( int headerIndex ) const;
  FlexHeader const* headerValues() const { return myHeaderValues; }
  /**
   * Dump data sample values to specified file stream
   */
  void dumpDataValues( std::FILE* stream );
  /**
   * Dump header values to specified file stream
   */
  void dumpHeaderValues( std::FILE* stream );

 private:
  void freeMemory();
  /**
   * Array containing data samples.
   * Number of allocated samples in this array is myNumAllocatedDataSamples.
   * Number  of actually stored data samples in this array is myNumStoredDataSamples
   */
  float* myData;
  /**
   * Vector of 'data chunks'
   * All data that didn't fit into the pre-allocated myData array are stored as 'data chunks' in this vector object.
   * Each data chunk has a field giving the number of data samples, and a float array containing the data samples.
   */
  std::vector<DataChunk*>* myDataChunks;
  FlexHeader* myHeaderValues;
  double mySampleRate; // Hz

  /// Number of samples stored in myData array
  int  myNumStoredDataSamples;
  /// Number of samples allocated in myData array
  int  myNumAllocatedDataSamples;
  /// Total number of samples stored in data chunks
  int  myNumDataChunkSamples;
};

class Exception {
public:
  Exception();
  Exception( std::string const& message );
  Exception( int dummy, char const* text );
  Exception( char const* text, ... );
  Exception( Exception const& obj );
  ~Exception();
  const char* getMessage();
protected:
  std::string myMessage;
};


} // namespace

#endif



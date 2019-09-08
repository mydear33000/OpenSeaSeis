/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csRSFHeader.h"
#include "csRSFReader.h"
#include <cstring>
#include <cmath>

#include "csStandardHeaders.h"
#include "csFileUtils.h"
#include "csGeolibUtils.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace cseis_io;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: INPUT_RSF
 *
 */

namespace mod_input_rsf {
  struct VariableStruct {
    long traceCounter;
    csRSFReader* rsfReader;
    csRSFHeader* rsfHdr;
    int nTracesToRead;
    string filename;

    int dim1Counter;
    int dim2Counter;

    int hdrId_dim2;
    int hdrId_dim3;
    int hdrId_trcno;

    int hdrId_time_samp1;
    int hdrId_time_samp1_us;
    int hdrId_delay_time;
    double delayTime;
    bool atEOF;

    int count_dim2;
    int count_dim3;
  };
}
using mod_input_rsf::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_rsf_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  //----------------------------------------------
  // Initialise
  
  vars->traceCounter = 0;
  vars->rsfReader   = NULL;
  vars->nTracesToRead = 0;

  vars->hdrId_time_samp1    = -1;
  vars->hdrId_time_samp1_us = -1;
  vars->hdrId_dim2 = -1;
  vars->hdrId_dim3 = -1;
  vars->hdrId_trcno = -1;
  vars->hdrId_delay_time = -1;
  vars->delayTime = 0;
  vars->atEOF = false;
  vars->rsfHdr = new csRSFHeader();
  vars->count_dim2 = 0;
  vars->count_dim3 = 0;

  //------------------------------------------------
  
  std::string headerName;
  std::string yesno;
  int numTracesBuffer = 0;

  //------------------------------------

  param->getString("filename",&vars->filename);

  //----------------------------------------------------

  if( param->exists("ntraces") ) {
    param->getInt( "ntraces", &vars->nTracesToRead );
    if( vars->nTracesToRead < 0 ) {
      vars->nTracesToRead = 0;
    }
  }
  if( vars->nTracesToRead <= 0 ) vars->nTracesToRead = -1;  // Do not bother how many traces, read in all

  if( param->exists( "ntraces_buffer" ) ) {
    param->getInt( "ntraces_buffer", &numTracesBuffer );
    if( numTracesBuffer < 0 || numTracesBuffer > 9999999 ) {
      log->warning("Number of buffered traces out of range (=%d). Changed to default.", numTracesBuffer);
      numTracesBuffer = 0;
    }
  }

  bool rev_byte_order = false;
  if( param->exists("reverse_byte_order") ) {
    string text;
    param->getString( "reverse_byte_order", &text );
    if( !text.compare("yes") ) {
      rev_byte_order = true;
    }
    else if( !text.compare("no") ) {
      rev_byte_order = false;
    }
    else {
      log->error("Unknown option: %s", text.c_str());
    }
  }

  int numSamplesOut = 0;
  if( param->exists("nsamples") ) {
    param->getInt( "nsamples", &numSamplesOut );
  }

  //-------------------------------------------------------
  //

  try {
    vars->rsfReader = new csRSFReader( vars->filename, numTracesBuffer, rev_byte_order );
    vars->rsfReader->initialize( vars->rsfHdr );
  }
  catch( csException& e ) {
    vars->rsfReader = NULL;
    log->error("Error when opening RSF file '%s'.\nSystem message: %s", vars->filename.c_str(), e.getMessage() );
  }

  if( edef->isDebug() ) {
    vars->rsfReader->dump( log->getFile() );
  }

  vars->hdrId_dim2 = hdef->addHeader( cseis_geolib::TYPE_DOUBLE, "dim2", "RSF data dimension 2" );
  vars->hdrId_dim3 = hdef->addHeader( cseis_geolib::TYPE_DOUBLE, "dim3", "RSF data dimension 3" );
  vars->hdrId_trcno = hdef->addStandardHeader( HDR_TRCNO.name );
  vars->hdrId_delay_time = hdef->addStandardHeader( HDR_DELAY_TIME.name );

  if( numSamplesOut == 0 ) {
    shdr->numSamples = vars->rsfReader->numSamples();
  }
  else {
    shdr->numSamples = numSamplesOut;
  }
  shdr->sampleInt = vars->rsfReader->sampleInt();
  vars->delayTime = vars->rsfHdr->o1;

  log->line("");
  log->line( "  File name:            %s", vars->filename.c_str());
  log->line( "  Sample interval [ms]: %f", shdr->sampleInt );
  log->line( "  Number of samples:    %d", shdr->numSamples );
  log->line("");
  log->line(" RSF ASCII file dump:");
  vars->rsfReader->dump( log->getFile() );

  vars->traceCounter = 0;
  vars->count_dim2   = 1;
  vars->count_dim3   = 1;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_input_rsf_(
                          csTrace* trace,
                          int* port,
                          csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  csExecPhaseDef*         edef = env->execPhaseDef;

  if( edef->isCleanup() ) {
    if( vars->rsfReader != NULL ) {
      delete vars->rsfReader;
      vars->rsfReader = NULL;
    }
    if( vars->rsfHdr != NULL ) {
      delete vars->rsfHdr;
      vars->rsfHdr = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  if( vars->atEOF ) return false;

  csSuperHeader const*    shdr = env->superHeader;
  csTraceHeader* trcHdr = trace->getTraceHeader();
  float* samples = trace->getTraceSamples();

  if( (vars->nTracesToRead > 0 && vars->nTracesToRead == vars->traceCounter) ||
      !vars->rsfReader->getNextTrace( (byte_t*)samples, shdr->numSamples ) ) {
    if( vars->traceCounter == 0 ) {
      log->warning("RSF file '%s' does not contain any data trace.", vars->rsfReader->filename() );
    }
    vars->atEOF = true;
    return false;
  }
  
  double dim2 = vars->rsfHdr->o2 + (vars->count_dim2-1) * vars->rsfHdr->d2;
  double dim3 = vars->rsfHdr->o3 + (vars->count_dim3-1) * vars->rsfHdr->d3;
  trcHdr->setDoubleValue( vars->hdrId_dim2, dim2 );
  trcHdr->setDoubleValue( vars->hdrId_dim3, dim3 );
  trcHdr->setDoubleValue( vars->hdrId_delay_time, vars->delayTime );

  vars->traceCounter++;
  trcHdr->setIntValue( vars->hdrId_trcno, vars->traceCounter );

  if( vars->count_dim2 == vars->rsfHdr->n2 ) {
    vars->count_dim2 = 1;
    vars->count_dim3 += 1;
  }
  else {
    vars->count_dim2 += 1;
  }

  return true;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_rsf_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_RSF", "Input data in RSF format" );

  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "nsamples", "Number of samples to read in", NUM_VALUES_FIXED,
                  "If number of samples in input data set is smaller, traces will be filled with zeros. Set 0 to set number of samples from input data set.");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of samples to read in" );

  pdef->addParam( "ntraces", "Number of traces to read in", NUM_VALUES_FIXED,
                  "Input of traces will stop when all traces have been read in, or if the number of traces specified has been reached. Traces will not be filled up to the specified range");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces to read in" );

  pdef->addParam( "ntraces_buffer", "Number of traces to read into buffer at once", NUM_VALUES_FIXED,
                  "Reading in a large number of traces at once may enhance performance, but requires more memory" );
  pdef->addValue( "20", VALTYPE_NUMBER, "Number of traces to buffer" );

  pdef->addParam( "reverse_byte_order", "Reverse byte order of input file (endian byte swapping)", NUM_VALUES_VARIABLE );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Reverse byte order of input file" );
  pdef->addOption( "no", "Do not reverse byte order" );
}

extern "C" void _params_mod_input_rsf_( csParamDef* pdef ) {
  params_mod_input_rsf_( pdef );
}
extern "C" void _init_mod_input_rsf_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_input_rsf_( param, env, log );
}
extern "C" bool _exec_mod_input_rsf_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_input_rsf_( trace, port, env, log );
}


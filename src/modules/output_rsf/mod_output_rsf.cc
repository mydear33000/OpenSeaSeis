/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csRSFWriter.h"
#include "csRSFHeader.h"
#include "csFileUtils.h"
#include <cstring>
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace cseis_io;
using namespace std;

namespace mod_output_rsf {
  struct VariableStruct {
    long traceCounter;
    csRSFWriter* rsfWriter;
    int   nTracesRead;
    int   numSamplesOut;
    bool special;
    std::string hdrName_dim2;
    std::string hdrName_dim3;
    int hdrId_dim2;
    int hdrId_dim3;
    bool isFirstCall;
  };
}
using namespace mod_output_rsf;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_output_rsf_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->traceCounter = 0;
  vars->rsfWriter   = NULL;
  vars->nTracesRead  = 0;
  vars->numSamplesOut = shdr->numSamples;
  vars->special = true;
  vars->isFirstCall = true;
//--------------------------------------------------
  std::string headerName;
  std::string yesno;
  int numTracesBuffer = 20;
  std::string filename = "";
  //---------------------------------------------------------
  bool outputGrid = false;
  if( param->exists("output_grid") ) {
    string text;
    param->getString("output_grid", &text);
    if( !text.compare("yes") ) {
      outputGrid = true;
    }
    else if( !text.compare("no") ) {
      outputGrid = false;
    }
    else {
      log->error("Unknown option: %s", text.c_str());
    }
  }
  bool swapDim3 = false;
  /*
  if( param->exists("swap_dim3") ) {
    string text;
    param->getString("swap_dim3", &text);
    if( !text.compare("yes") ) {
      swapDim3 = true;
    }
    else if( !text.compare("no") ) {
      swapDim3 = false;
    }
    else {
      log->error("Unknown option: %s", text.c_str());
    }
  }
  if( swapDim3 ) {
    log->error("Swapping dimension 1 and 3 is not supported yet" );
  }
  */
  //--------------------------------------------------
  param->getString( "filename", &filename );
  if( filename.length() > 4 && !filename.substr(filename.length()-4,4).compare(".rsf") ) {
    filename = filename.substr(0,filename.length()-4);
    log->line("RSF output filename stump = '%s'\n", filename.c_str());
  }

  if( param->exists("nsamples" ) ) {
    param->getInt( "nsamples", &vars->numSamplesOut );
    if( vars->numSamplesOut <= 0 || vars->numSamplesOut > shdr->numSamples ) {
      log->error("Inconsistent number of samples specified: %d. Actual number of samples in trace: %d", vars->numSamplesOut, shdr->numSamples);
    }
  }

  double tolerance = 1e-5;
  if( param->exists("tolerance" ) ) {
    param->getDouble( "tolerance", &tolerance );
    tolerance = fabs(tolerance);
  }

  bool rev_byte_order = false;
  if( param->exists("reverse_byte_order") ) {
    string text;
    param->getString( "reverse_byte_order", &text );
    if( !text.compare("yes") ) {
      rev_byte_order = true;
    }
  }

  //---------------------------------------------------------
  if( param->exists( "ntraces_buffer" ) ) {
    param->getInt( "ntraces_buffer", &numTracesBuffer );
    if( numTracesBuffer <= 0 || numTracesBuffer > 999999 ) {
      log->warning("Number of buffered traces out of range (=%d). Changed to 20.", numTracesBuffer);
      numTracesBuffer = 20;
    }
  }
  else {
    numTracesBuffer = 20;
  }

  //--------------------------------------------------------------------------------
  // Set RSF header
  param->getString( "hdr_dim2", &vars->hdrName_dim2 );
  param->getString( "hdr_dim3", &vars->hdrName_dim3 );

  vars->hdrId_dim2 = hdef->headerIndex( vars->hdrName_dim2 );
  vars->hdrId_dim3 = hdef->headerIndex( vars->hdrName_dim3 );

  if( hdef->headerType( vars->hdrId_dim2 ) == TYPE_STRING ) {
    log->error("DIM2 trace header '%s' must be of numerical type", vars->hdrName_dim2.c_str());
  }

  cseis_io::csRSFHeader rsfHdr;

  /*
  param->getInt( "dim2", &rsfHdr.n2, 0 );
  param->getDouble( "dim2", &rsfHdr.o2, 1 );
  param->getDouble( "dim2", &rsfHdr.d2, 2 );
  if( param->getNumValues( "dim2" ) > 3 ) {
    param->getDouble( "dim2", &rsfHdr.e2, 3 );
  }
  else {
    rsfHdr.e2 = rsfHdr.o2 + rsfHdr.d2 * (rsfHdr.n2-1);
  }
  param->getInt( "dim3", &rsfHdr.n3, 0 );
  param->getDouble( "dim3", &rsfHdr.o3, 1 );
  param->getDouble( "dim3", &rsfHdr.d3, 2 );
  if( param->getNumValues( "dim3" ) > 3 ) {
    param->getDouble( "dim3", &rsfHdr.e3, 3 );
  }
  else {
    rsfHdr.e3 = rsfHdr.o3 + rsfHdr.d3 * (rsfHdr.n3-1);
  }
  */

  rsfHdr.n1 = shdr->numSamples;
  rsfHdr.o1 = 0;
  rsfHdr.d1 = shdr->sampleInt;
  rsfHdr.e1 = (shdr->numSamples-1) * shdr->sampleInt;

  if( param->exists("world_p1") ) {
    param->getDouble( "world_p1", &rsfHdr.world_x1, 0 );
    param->getDouble( "world_p1", &rsfHdr.world_y1, 1 );
    param->getDouble( "world_p1", &rsfHdr.il1, 2 );
    param->getDouble( "world_p1", &rsfHdr.xl1, 3 );

    param->getDouble( "world_p2", &rsfHdr.world_x2, 0 );
    param->getDouble( "world_p2", &rsfHdr.world_y2, 1 );
    param->getDouble( "world_p2", &rsfHdr.il2, 2 );
    param->getDouble( "world_p2", &rsfHdr.xl2, 3 );

    param->getDouble( "world_p3", &rsfHdr.world_x3, 0 );
    param->getDouble( "world_p3", &rsfHdr.world_y3, 1 );
    param->getDouble( "world_p3", &rsfHdr.il3, 2 );
    param->getDouble( "world_p3", &rsfHdr.xl3, 3 );

    param->getDouble( "ild_xld", &rsfHdr.ild, 0 );
    param->getDouble( "ild_xld", &rsfHdr.xld, 1 );
    outputGrid = true;
  }
  else if( outputGrid ) {
    log->error("User specified output of grid definition, but grid point 1 is not defined. Grid cannot be automatically defined.");
  }

  if( param->exists("data_format") ) {
    string text;
    param->getString( "data_format", &text );
    if( !text.compare("native_float") ) {
      rsfHdr.data_format = csRSFHeader::DATA_FORMAT_FLOAT;
      rsfHdr.esize = 4;
    }
    else {
      log->error("Unknown data sample format: '%s'", text.c_str());
    }
  }
  else {
    rsfHdr.data_format = csRSFHeader::DATA_FORMAT_FLOAT;
    rsfHdr.esize = 4;
  }

  if( param->exists("filename_bin") ) {
    string text;
    param->getString( "filename_bin", &text );
  }
  else {
    int length = filename.length();
    string binExt(".bin");
    if( length < 5 || filename.substr(length-5,4).compare(".rsf" ) ) {
      rsfHdr.filename_bin = filename;
      rsfHdr.filename_bin.append(".bin");
      filename.append(".rsf");
      rsfHdr.filename_bin_full_path = rsfHdr.filename_bin;
    }
    else {
      rsfHdr.filename_bin = filename.c_str();
      rsfHdr.filename_bin.append(".bin");
      rsfHdr.filename_bin_full_path = rsfHdr.filename_bin;
    }
  }

  //----------------------------------------------------
  //
  try {
    vars->rsfWriter = new csRSFWriter( filename, numTracesBuffer, rev_byte_order, swapDim3, outputGrid, tolerance );
  }
  catch( csException& e ) {
    vars->rsfWriter = NULL;
    log->error("Error when creating RSF writer.\nSystem message: %s", filename.c_str(), e.getMessage() );
  }

  //----------------------------------------------------
  try {
    vars->rsfWriter->initialize( &rsfHdr );
  }
  catch( csException& e ) {
    delete vars->rsfWriter;
    vars->rsfWriter = NULL;
    log->error("Error when initializing RSF writer.\nSystem message: %s", e.getMessage() );
  }

  log->line("");
  log->line("  File name:            %s", filename.c_str());
  log->line("  Sample interval [ms]: %f", shdr->sampleInt);
  log->line("  Number of samples:    %d", shdr->numSamples );
  log->line("  Sample data format:   %d", rsfHdr.data_format );
  log->line("");

  vars->traceCounter = 0;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_output_rsf_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  //  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup() ) {
    if( vars->rsfWriter != NULL ) {
      vars->rsfWriter->finalize();
      delete vars->rsfWriter;
      vars->rsfWriter = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }
  
  if( vars->isFirstCall ) {
    vars->isFirstCall = false;
    try {
      vars->rsfWriter->openFile();
    }
    catch( csException& e ) {
      log->error("Error when opening RSF output file.\nSystem message: %s", e.getMessage() );
    }
  }

  float* samples = trace->getTraceSamples();
  csTraceHeader* trcHdr = trace->getTraceHeader();

  double val_dim2 = trcHdr->doubleValue( vars->hdrId_dim2 );
  double val_dim3 = trcHdr->doubleValue( vars->hdrId_dim3 );

  try {
    vars->rsfWriter->writeNextTrace( (byte_t*)samples, vars->numSamplesOut, val_dim2, val_dim3 );
  }
  catch( csException& exc ) {
    log->error( "System message: %s" ,exc.getMessage() );
  }
  
  vars->traceCounter++;

  return true;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_output_rsf_( csParamDef* pdef ) {

  pdef->setModule( "OUTPUT_RSF", "Output data in RSF format" );
  pdef->addParam( "filename", "Output file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Output file name" );

  pdef->addParam( "hdr_dim2", "Trace header defining dimension 2", NUM_VALUES_FIXED,
                  "dimension 2 is the 'faster' dimension compared to dimension 3. In other words, dimension 2 trace header values are expected to change from traqce to trace, while dimension 3 trace header values should only change from section to section." );
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header defining dimension 2" );

  pdef->addParam( "hdr_dim3", "Trace header defining dimension 3", NUM_VALUES_FIXED, "dimension 3 is the 'slowest' dimension." );
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header defining dimension 3" );

  /*
  pdef->addParam( "dim2", "Dimension 2 definition", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "n2: Number of traces in dimension 2" );
  pdef->addValue( "", VALTYPE_NUMBER, "o2: Dimension 2 origin" );
  pdef->addValue( "", VALTYPE_NUMBER, "d2: Dimension 2 step/increment" );
  pdef->addValue( "", VALTYPE_NUMBER, "e2: Dimension 2 last coordinate" );

  pdef->addParam( "dim3", "Dimension 3 definition", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "n3: Number of traces in dimension 3" );
  pdef->addValue( "", VALTYPE_NUMBER, "o3: Dimension 3 origin" );
  pdef->addValue( "", VALTYPE_NUMBER, "d3: Dimension 3 step/increment" );
  pdef->addValue( "", VALTYPE_NUMBER, "e3: Dimension 3 last coordinate" );
  */

  pdef->addParam( "world_p1", "Grid definition: World point 1", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "X coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Y coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Inline number" );
  pdef->addValue( "", VALTYPE_NUMBER, "Crossline number" );

  pdef->addParam( "world_p2", "Grid definition: World point 2", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "X coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Y coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Inline number" );
  pdef->addValue( "", VALTYPE_NUMBER, "Crossline number" );

  pdef->addParam( "world_p3", "Grid definition: World point 3", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "X coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Y coordinate" );
  pdef->addValue( "", VALTYPE_NUMBER, "Inline number" );
  pdef->addValue( "", VALTYPE_NUMBER, "Crossline number" );

  pdef->addParam( "ild_xld", "Grid definition: Bin/cell size", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "ild: Bin/cell size in inline/row direction [m]" );
  pdef->addValue( "", VALTYPE_NUMBER, "xld: Bin/cell size in crossline/col direction [m]" );

  pdef->addParam( "data_format", "Outyput format for data samples", NUM_VALUES_FIXED);
  pdef->addValue( "native_float", VALTYPE_OPTION );
  pdef->addOption( "native_float", "Use native 32-byte floating point format" );
 
  pdef->addParam( "ntraces_buffer", "Number of traces to buffer before write operation", NUM_VALUES_FIXED,
                  "Writing a large number of traces at once enhances performance, but requires more memory" );
  pdef->addValue( "20", VALTYPE_NUMBER, "Number of traces to buffer before writing" );

  pdef->addParam( "reverse_byte_order", "Reverse byte order in output file (endian byte swapping)", NUM_VALUES_FIXED, "Setting this to 'yes' means that the output file will be written in Little endian byte order" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Reverse byte order for output file" );
  pdef->addOption( "no", "Do not reverse byte order" );

  pdef->addParam( "tolerance", "Tolerance value", NUM_VALUES_FIXED, "The tolerance is used to check dimension 2 and 3 trace header values" );
  pdef->addValue( "1e-5", VALTYPE_NUMBER, "Any difference between nominal and actual header values above the given tolerance will terminate the program" );

  //  pdef->addParam( "swap_dim3", "Swap dimension 3 and 1", NUM_VALUES_FIXED, "...create outputfile for Lysa modeling" );
  //  pdef->addValue( "no", VALTYPE_OPTION );
  //  pdef->addOption( "yes", "Swap dimension 1 and 3" );
  //  pdef->addOption( "no", "Do not swap dimensions" );

  pdef->addParam( "output_grid", "Output grid definition?", NUM_VALUES_FIXED,
		  "By default, grid definition is only output if three world points have been specified" );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Write grid definition to output file" );
  pdef->addOption( "no", "Do not output grid definition to output file" );  
}

extern "C" void _params_mod_output_rsf_( csParamDef* pdef ) {
  params_mod_output_rsf_( pdef );
}
extern "C" void _init_mod_output_rsf_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_output_rsf_( param, env, log );
}
extern "C" bool _exec_mod_output_rsf_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_output_rsf_( trace, port, env, log );
}


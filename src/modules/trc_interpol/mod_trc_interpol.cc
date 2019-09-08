/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFFTTools.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: TRC_INTERPOL
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_trc_interpol {
  struct VariableStruct {
    int method;
    int mode;
    int hdrID_in;
    int hdrID_out;
    double value;
    cseis_geolib::csFFTTools* fftTool1;
    cseis_geolib::csFFTTools* fftTool2;
    float* buffer;
  };
  static int const METHOD_SIMPLE_AVERAGE = 1;
  static int const METHOD_FX_AVERAGE     = 2;
  
  static int const MODE_ENSEMBLE = 11;
  static int const MODE_FIXED    = 12;
  static int const MODE_HEADER   = 13;
}

using namespace mod_trc_interpol;

//*******************************************************************
//
//  Init phase
//
//*******************************************************************

void init_mod_trc_interpol_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  vars->method = METHOD_SIMPLE_AVERAGE;
  vars->mode = MODE_FIXED;
  vars->hdrID_out = -1;
  vars->hdrID_in = -1;
  vars->value = 0.0;
  vars->fftTool1 = NULL;
  vars->fftTool2 = NULL;
  vars->buffer   = NULL;
  
  std::string text;
  if( param->exists("mode") ) {
    param->getString("mode",&text);
    if( !text.compare("fixed") ) {
      vars->mode = MODE_FIXED;
    }
    else if( !text.compare("ensemble") ) {
      vars->mode = MODE_ENSEMBLE;
    }
    else if( !text.compare("header") ) {
      vars->mode = MODE_HEADER;
      param->getString("header",&text);
      param->getDouble("header",&vars->value,1);
      vars->hdrID_in = hdef->headerIndex( text.c_str() );
    }
    else {
      log->error("Unknown option: %s", text.c_str() );
    }
  }

  if( param->exists("method") ) {
    param->getString("method",&text);
    if( !text.compare("simple_average") ) {
      vars->method = METHOD_SIMPLE_AVERAGE;
    }
    else if( !text.compare("fx_average") ) {
      vars->method = METHOD_FX_AVERAGE;
    }
    else {
      log->error("Unknown option: %s", text.c_str() );
    }
  }

  edef->setExecType( EXEC_TYPE_MULTITRACE );
  if( vars->mode == MODE_ENSEMBLE ) {
    edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }
  else if( vars->mode == MODE_FIXED ) {
    edef->setTraceSelectionMode( TRCMODE_FIXED, 2 );
  }
  else {
    edef->setTraceSelectionMode( TRCMODE_FIXED, 3 );
  }

  if( !hdef->headerExists("interpolated") ) {
    hdef->addHeader( TYPE_INT, "interpolated", "1: Trace is interpolated" );
  }
  vars->hdrID_out = hdef->headerIndex("interpolated");

  if( vars->method == METHOD_FX_AVERAGE ) {
    vars->fftTool1 = new cseis_geolib::csFFTTools( shdr->numSamples );
    vars->fftTool2 = new cseis_geolib::csFFTTools( shdr->numSamples );
    vars->buffer = new float[shdr->numSamples];
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************

void exec_mod_trc_interpol_(
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
    if( vars->fftTool1 != NULL ) {
      delete vars->fftTool1;
      vars->fftTool1 = NULL;
    }
    if( vars->fftTool2 != NULL ) {
      delete vars->fftTool2;
      vars->fftTool2 = NULL;
    }
    delete vars; vars = NULL;
    return;
  }

  if( edef->isDebug() ) log->line("Number of input traces: %d, num trc to keep: %d", traceGather->numTraces(), *numTrcToKeep);

  int numTraces = traceGather->numTraces();
  if( numTraces < 2 ) {  // This must be the last trace of the data set. Simply return, do nothing
    return;
  }

  int nSamples = shdr->numSamples;
  if( vars->mode == MODE_FIXED ) {
    // 1) Create new trace in between the two traces
    traceGather->createTraces( 1, 1, hdef, nSamples );

    // 2) Interpolate trace samples
    float* samples1 = traceGather->trace(0)->getTraceSamples();
    float* samples2 = traceGather->trace(2)->getTraceSamples();
    float* samplesNew = traceGather->trace(1)->getTraceSamples();
  
    if( vars->method == METHOD_SIMPLE_AVERAGE ) {
      for( int isamp = 0; isamp < nSamples; isamp++ ) {
	samplesNew[isamp] = (samples1[isamp] + samples2[isamp])*0.5;
      }
    }

    // 2) Interpolate trace header values
    csTraceHeader* trcHdr1 = traceGather->trace(0)->getTraceHeader();
    csTraceHeader* trcHdr2 = traceGather->trace(2)->getTraceHeader();
    csTraceHeader* trcHdrNew = traceGather->trace(1)->getTraceHeader();

    for( int ihdr = 0; ihdr < hdef->numHeaders(); ihdr++ ) {
      char type = hdef->headerType(ihdr);
      switch( type ) {
      case TYPE_INT:
        trcHdrNew->setIntValue( ihdr, (trcHdr1->intValue(ihdr)+trcHdr2->intValue(ihdr))/2 );
        break;
      case TYPE_FLOAT:
        trcHdrNew->setFloatValue( ihdr, (trcHdr1->floatValue(ihdr)+trcHdr2->floatValue(ihdr))/2 );
        break;
      case TYPE_DOUBLE:
        trcHdrNew->setDoubleValue( ihdr, (trcHdr1->doubleValue(ihdr)+trcHdr2->doubleValue(ihdr))/2 );
        break;
      case TYPE_INT64:
        trcHdrNew->setInt64Value( ihdr, (trcHdr1->int64Value(ihdr)+trcHdr2->int64Value(ihdr))/2 );
        break;
      case TYPE_STRING:
        trcHdrNew->setStringValue( ihdr, trcHdr1->stringValue(ihdr) );
        break;
      default:
        log->error("mod_trc_interpol: Unknown trace header type, code: %d", type);
      }
    }
    trcHdrNew->setIntValue( vars->hdrID_out, 1 );

    *numTrcToKeep = 1;  // Keep last trace for next pair of traces to interpolate in between
  } // END FIXED
  else if( vars->mode == MODE_ENSEMBLE ) {
    for( int itrc = numTraces-1; itrc > 0; itrc-- ) {
      int trcIndex1 = itrc-1;
      int trcIndex2 = itrc+1;
      // 1) Create new trace in between the two traces
      traceGather->createTrace( itrc, hdef, nSamples );

      // 2) Interpolate trace samples
      float* samples1 = traceGather->trace(trcIndex1)->getTraceSamples();
      float* samples2 = traceGather->trace(trcIndex2)->getTraceSamples();
      float* samplesNew = traceGather->trace(itrc)->getTraceSamples();
  
      if( vars->method == METHOD_SIMPLE_AVERAGE ) {
	for( int isamp = 0; isamp < nSamples; isamp++ ) {
	  samplesNew[isamp] = (samples1[isamp] + samples2[isamp])*0.5;
	}
      }
      else {
        bool success = vars->fftTool1->fft_forward( samples1 );
        if( !success ) {
          log->error("Error occurred in fft forward transform");
        }
        success = vars->fftTool2->fft_forward( samples2 );
        double* real1 = vars->fftTool1->getRealDataPointer();
        double* imag1 = vars->fftTool1->getImagDataPointer();
        double const* real2 = vars->fftTool2->realData();
        double const* imag2 = vars->fftTool2->imagData();
	for( int isamp = 0; isamp < vars->fftTool1->numFFTSamples(); isamp++ ) {
          real1[isamp] = (real1[isamp] + real2[isamp])*0.5;
          imag1[isamp] = (imag1[isamp] + imag2[isamp])*0.5;
          /*          float amp1 = 2.0 * sqrt(real1[isamp]*real1[isamp] + imag1[isamp]*imag1[isamp]);
          float amp2 = 2.0 * sqrt(real2[isamp]*real2[isamp] + imag2[isamp]*imag2[isamp]);
          float amp = 0.5 * ( amp1 + amp2 );
          float phase1 = atan2(-imag1[isamp],real1[isamp]);
          float phase2 = atan2(-imag2[isamp],real2[isamp]);
          float phase = atan2( sin(phase1)+sin(phase2), cos(phase1)+cos(phase2) );
          real1[isamp] = 0.5*amp * cos(phase);
          imag1[isamp] = -0.5*amp * sin(phase); */
        }
        success = vars->fftTool1->fft_inverse( );
        real1 = vars->fftTool1->getRealDataPointer();
	for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
	  samplesNew[isamp] = real1[isamp];
	}

      }

      // 2) Interpolate trace header values
      csTraceHeader* trcHdr1 = traceGather->trace(trcIndex1)->getTraceHeader();
      csTraceHeader* trcHdr2 = traceGather->trace(trcIndex2)->getTraceHeader();
      csTraceHeader* trcHdrNew = traceGather->trace(itrc)->getTraceHeader();

      for( int ihdr = 0; ihdr < hdef->numHeaders(); ihdr++ ) {
	char type = hdef->headerType(ihdr);
	switch( type ) {
	case TYPE_INT:
	  trcHdrNew->setIntValue( ihdr, (trcHdr1->intValue(ihdr)+trcHdr2->intValue(ihdr))/2 );
	  break;
	case TYPE_FLOAT:
	  trcHdrNew->setFloatValue( ihdr, (trcHdr1->floatValue(ihdr)+trcHdr2->floatValue(ihdr))/2 );
	  break;
	case TYPE_DOUBLE:
	  trcHdrNew->setDoubleValue( ihdr, (trcHdr1->doubleValue(ihdr)+trcHdr2->doubleValue(ihdr))/2 );
	  break;
	case TYPE_INT64:
	  trcHdrNew->setInt64Value( ihdr, (trcHdr1->int64Value(ihdr)+trcHdr2->int64Value(ihdr))/2 );
	  break;
	case TYPE_STRING:
	  trcHdrNew->setStringValue( ihdr, trcHdr1->stringValue(ihdr) );
	  break;
	default:
	  log->error("mod_trc_interpol: Unknown trace header type, code: %d", type);
	}
      }
      trcHdrNew->setIntValue( vars->hdrID_out, 1 );

    } // END for itrc
  } // END if MODE_ENSEMBLE
  else {
    if( numTraces < 3 ) {  // This must be the last trace of the data set. Simply return, do nothing
      return;
    }
    double value1 = traceGather->trace(1)->getTraceHeader()->doubleValue( vars->hdrID_in );
    double value2 = traceGather->trace(2)->getTraceHeader()->doubleValue( vars->hdrID_in );
    if( value1 == vars->value ) {
      // 2) Interpolate trace samples
      float* samples1 = traceGather->trace(0)->getTraceSamples();
      float* samples2 = traceGather->trace(2)->getTraceSamples();
      float* samplesMid = traceGather->trace(1)->getTraceSamples();
  
      if( vars->method == METHOD_SIMPLE_AVERAGE ) {
	for( int isamp = 0; isamp < nSamples; isamp++ ) {
	  samplesMid[isamp] = (samples1[isamp] + samples2[isamp])*0.5;
	}
      }
      traceGather->trace(1)->getTraceHeader()->setIntValue( vars->hdrID_out, 1 );
    }
    if( value2 == vars->value ) {
      *numTrcToKeep = 2;
    }
    else {
      *numTrcToKeep = 1;
    }
  } // END if MODE_HEADER
}


//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_trc_interpol_( csParamDef* pdef ) {
  pdef->setModule( "TRC_INTERPOL", "Interpolate traces (EXPERIMENTAL MODULE)", "Interpolate 1 trace between every two adjacent traces, e.g. if 10 traces are input, 19 traces will be output" );

  pdef->addParam( "method", "", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "simple_average", "Interpolate new trace by simply averaging adjacent traces. Header values are averaged." );
  pdef->addOption( "fx_average", "..." );

  pdef->addParam( "mode", "", NUM_VALUES_FIXED );
  pdef->addValue( "fixed", VALTYPE_OPTION );
  pdef->addOption( "fixed", "Input 2 sequential traces at a time" );
  pdef->addOption( "ensemble", "Input one ensemble at a time" );
  pdef->addOption( "header", "Trace header indicates traces which shall be replaced by an interpolated trace" );

  pdef->addParam( "header", "Trace header indicating trace which shall be interpolated", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Value indicating this trace shall be replaced/interpolated" );
}

extern "C" void _params_mod_trc_interpol_( csParamDef* pdef ) {
  params_mod_trc_interpol_( pdef );
}
extern "C" void _init_mod_trc_interpol_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_trc_interpol_( param, env, log );
}
extern "C" void _exec_mod_trc_interpol_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_trc_interpol_( traceGather, port, numTrcToKeep, env, log );
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csSort.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/*
 * This module has a known bug when padding:
 * If header values outside specified range exists, a segmentation error may occur
 * Workaround: delete all traces outside specified pad range before calling this module
*/

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: TRC_ADD_ENS
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_trc_add_ens {
  struct VariableStruct {
    int numTraces;
    float value;
    int method;
    int hdrID_pad;
    int hdrType_pad;
    int hdrID_padFlag;
    csFlexNumber start;
    csFlexNumber stop;
    csFlexNumber inc;
    bool deleteInconsistentTraces;
  };
  static const int METHOD_PAD     = 71;
  static const int METHOD_NTRACES = 72;
}
using mod_trc_add_ens::VariableStruct;

//*******************************************************************
//
//  Init phase
//
//*******************************************************************

void init_mod_trc_add_ens_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_MULTITRACE );
  env->execPhaseDef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->numTraces      = 0;
  vars->value          = 0;
  vars->start          = 0;
  vars->stop           = 0;
  vars->inc            = 0;
  vars->hdrID_padFlag  = -1;
  vars->hdrID_pad      = -1;
  vars->method         = mod_trc_add_ens::METHOD_NTRACES;
  vars->deleteInconsistentTraces = false;

  //-----------------------------------------------
  string text;
  if( param->exists("method") ) {
    param->getString( "method", &text );
    if( !text.compare("ntraces") ) {
      vars->method = mod_trc_add_ens::METHOD_NTRACES;
    }
    else if( !text.compare("pad") ) {
      vars->method = mod_trc_add_ens::METHOD_PAD;
      if( param->exists("delete") ) {
        string text;
        param->getString("delete",&text);
        if( !text.compare("yes") ) {
          vars->deleteInconsistentTraces = true;
        }
        else if( !text.compare("no") ) {
          vars->deleteInconsistentTraces = false;
        }
        else {
          log->error("Unknown option: %s", text.c_str());
        }
      }
    }
    else {
      log->error("Unknown option: %s", text.c_str());
    }
  }

  //-----------------------------------------------
  if( vars->method == mod_trc_add_ens::METHOD_NTRACES ) {
    param->getInt( "ntraces", &vars->numTraces );
  }
  else { // method PAD
    param->getString( "header", &text );
    if( !hdef->headerExists(text) ) {
      log->error("Trace header does not exist: %s", text.c_str() );
    }
    if( !hdef->headerExists("pad_flag") ) {
      hdef->addHeader(cseis_geolib::TYPE_INT, "pad_flag", "Padded trace? 0:no, 1:yes");
    }
    vars->hdrID_padFlag = hdef->headerIndex("pad_flag");
    vars->hdrID_pad   = hdef->headerIndex(text);
    vars->hdrType_pad = hdef->headerType(text);
    if( vars->hdrType_pad == TYPE_STRING || vars->hdrType_pad == TYPE_CHAR ) {
      log->error("Trace header does not have a number type. This is currently not supported.");
    }
    else if( vars->hdrType_pad != TYPE_INT ) {
      log->error("Trace header has floating point type, or 64bit integer. Only 32bit integer type trace headers are currently supported for trace padding.");
    }
    if( vars->hdrType_pad == TYPE_INT ) {
      int tmp3;
      param->getInt( "pad_inc", &tmp3 );
      vars->inc.setIntValue(tmp3);
      if( param->exists("pad") ) {
        int tmp1, tmp2;
        param->getInt( "pad", &tmp1, 0 );
        vars->start.setIntValue(tmp1);
        param->getInt( "pad", &tmp2, 1 );
        vars->stop.setIntValue(tmp2);
        vars->numTraces = (tmp2 - tmp1) / tmp3 + 1;
      }
    }
    else {
      double tmp3;
      param->getDouble( "pad_inc", &tmp3 );
      vars->inc.setDoubleValue(tmp3);
      if( param->exists("pad") ) {
        double tmp1, tmp2;
        param->getDouble( "pad", &tmp1, 0 );
        vars->start.setDoubleValue(tmp1);
        param->getDouble( "pad", &tmp2, 1 );
        vars->start.setDoubleValue(tmp2);
        vars->numTraces = (int)round((tmp2 - tmp1) / tmp3) + 1;
      }
    }
    if( param->exists("pad") ) log->line("Number of output traces per ensemble: %d", vars->numTraces );
  }

  if( param->exists("value") ) {
    param->getFloat( "value", &vars->value );
  }

  if( edef->isDebug() ) {
    log->line("Start/stop/inc: %f %f %f", vars->start.doubleValue(), vars->stop.doubleValue(), vars->inc.doubleValue() );
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_trc_add_ens_(
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
    delete vars; vars = NULL;
    return;
  }

  int numTracesOrig = traceGather->numTraces();
  if( numTracesOrig == 0 ) return;

  //------------------------------------------------
  if( vars->method == mod_trc_add_ens::METHOD_NTRACES ) {
    traceGather->createTraces( numTracesOrig, vars->numTraces, hdef, shdr->numSamples );

    csTraceHeader* trcHdrOrig = traceGather->trace(numTracesOrig-1)->getTraceHeader();
    for( int itrc = 0; itrc < vars->numTraces; itrc++ ) {
      int trcIndexNew = numTracesOrig + itrc;

      // Copy header data
      csTraceHeader* trcHdrCopy = traceGather->trace(trcIndexNew)->getTraceHeader();
      trcHdrCopy->copyFrom( trcHdrOrig );

      float* samples = traceGather->trace(trcIndexNew)->getTraceSamples();
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] = vars->value;
      }
    }
  }
  //------------------------------------------------
  else {  // method PAD
    int* indexBuffer;
    indexBuffer = new int[numTracesOrig];
    if( vars->hdrType_pad == TYPE_INT ) {
      int* hdrValueBuffer = new int[numTracesOrig];
      for( int itrc = 0; itrc < numTracesOrig; itrc++ ) {
        csTraceHeader* trcHdr = traceGather->trace(itrc)->getTraceHeader();
        trcHdr->setIntValue( vars->hdrID_padFlag, 0 );
        int value = trcHdr->intValue(vars->hdrID_pad);
        if( (value-vars->start.intValue()) % vars->inc.intValue() != 0 ) {
          log->error("Inconsistent trace header value: %d. Expected values should fit into the specified padding range from %d to %d, inc %d.",
            value, vars->start.intValue(), vars->stop.intValue(), vars->inc.intValue() );
        }
        hdrValueBuffer[itrc] = value;
        indexBuffer[itrc] = itrc;
        if( edef->isDebug() ) {
          log->line("Trace header value, trace #%-3d: %d", itrc+1, value);
        }
      }
      csSort<int>().simpleSortIndex( hdrValueBuffer, numTracesOrig, indexBuffer );
//      for( int itrc = 1; itrc < numTracesOrig; itrc++ ) {
//        if( indexBuffer[itrc] < indexBuffer[itrc-1] ) {
//        }
//      }
//      int minValue = hdrValueBuffer[0];
//      int maxValue = hdrValueBuffer[1];
//      if( vars->numTraces == 0 ) {
//      }

      if( numTracesOrig != vars->numTraces ) {
	//	fprintf(stdout,"Num traces orig: %d  / %d\n", numTracesOrig, vars->numTraces);
        traceGather->createTraces( numTracesOrig, vars->numTraces-numTracesOrig, hdef, shdr->numSamples );
        csTraceHeader* trcHdrOrig = traceGather->trace(0)->getTraceHeader();
        for( int itrc = numTracesOrig; itrc < vars->numTraces; itrc++ ) {
          traceGather->trace(itrc)->getTraceHeader()->copyFrom( trcHdrOrig );
          traceGather->trace(itrc)->getTraceHeader()->setIntValue( vars->hdrID_padFlag, 1 );
          float* samples = traceGather->trace(itrc)->getTraceSamples();
          for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
            samples[isamp] = vars->value;
          }
        }
        csTrace** tracePtr = new csTrace*[vars->numTraces];
        int currentTrcHdrValue = 0;
        int indexInputTrace    = 0;
        int indexNewTrace      = numTracesOrig;
	//	fprintf(stdout,"Num traces orig: %d\n", numTracesOrig);
        for( int itrc = 0; itrc < vars->numTraces; itrc++ ) {
          int value = vars->start.intValue() + itrc*vars->inc.intValue();
          if( indexInputTrace < numTracesOrig ) {
            currentTrcHdrValue = traceGather->trace( indexBuffer[indexInputTrace] )->getTraceHeader()->intValue(vars->hdrID_pad);
            if( edef->isDebug() ) {
              log->line("Current value: %d, computed value: %d", currentTrcHdrValue, value );
            }
            if( value == currentTrcHdrValue ) {
              tracePtr[itrc] = traceGather->trace( indexBuffer[indexInputTrace] );
              indexInputTrace += 1;
            }
            else if( value > currentTrcHdrValue ) {
              if( vars->deleteInconsistentTraces ) {
                log->warning("Inconsistent trace header value (%d > %d): This trace is duplicate and will be deleted. Index: %d %d (%d)",
                    value, currentTrcHdrValue, indexInputTrace, indexBuffer[indexInputTrace], numTracesOrig );
                tracePtr[itrc] = traceGather->trace( indexBuffer[indexInputTrace] );
                indexInputTrace += 1;
              }
              else {
                log->error("Inconsistent trace header value (%d > %d): This is a duplicate trace or a program bug in the module trc_add_ens. Index: %d %d (%d)",
                           value, currentTrcHdrValue, indexInputTrace, indexBuffer[indexInputTrace], numTracesOrig );
              }
            }
            else {
	      //	      fprintf(stdout," CHECK %d %d %d\n", itrc, indexNewTrace, traceGather->numTraces() );
	      //	      fflush(stdout);
              tracePtr[itrc] = traceGather->trace( indexNewTrace );
              tracePtr[itrc]->getTraceHeader()->setIntValue( vars->hdrID_pad, value );
              indexNewTrace += 1;
            }
          }
          else {
            tracePtr[itrc] = traceGather->trace( indexNewTrace );
            tracePtr[itrc]->getTraceHeader()->setIntValue( vars->hdrID_pad, value );
            indexNewTrace += 1;
          }
        }
        for( int itrc = 0; itrc < vars->numTraces; itrc++ ) {
          (*traceGather)[itrc] = tracePtr[itrc];
        }
        delete [] tracePtr;
      }
      delete [] hdrValueBuffer;
    }
    else {
//      double* hdrValueBuffer = new double[numTracesOrig];
//      delete [] hdrValueBuffer;
    }
    delete [] indexBuffer;
  }
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_trc_add_ens_( csParamDef* pdef ) {
  pdef->setModule( "TRC_ADD_ENS", "Ensemble trace adding" );

  pdef->addParam( "method", "Method for trace adding", NUM_VALUES_FIXED );
  pdef->addValue( "ntraces", VALTYPE_OPTION );
  pdef->addOption( "ntraces", "Add specified number of traces at beginning and end of ensemble" );
  pdef->addOption( "pad", "Pad traces as specified by trace header name and padding values" );

  pdef->addParam( "ntraces", "Number of traces to add at end of each ensemble", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER );

  pdef->addParam( "value", "Initialize trace samples to the given value", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER );

  pdef->addParam( "header", "Trace header name for trace padding", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "pad", "Pad traces following the specified trace header values", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "First value" );
  pdef->addValue( "", VALTYPE_NUMBER, "Last value" );

  pdef->addParam( "pad_inc", "Pad traces with the specified increment between trace header values", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Increment" );

  pdef->addParam( "delete", "Delete inconsistent traces (when trace padding)?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Delete inconsistent traces" );
  pdef->addOption( "no", "Do not delete inconsistent traces" );
}

extern "C" void _params_mod_trc_add_ens_( csParamDef* pdef ) {
  params_mod_trc_add_ens_( pdef );
}
extern "C" void _init_mod_trc_add_ens_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_trc_add_ens_( param, env, log );
}
extern "C" void _exec_mod_trc_add_ens_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_trc_add_ens_( traceGather, port, numTrcToKeep, env, log );
}


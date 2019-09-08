/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: DEBIAS
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_debias {
  struct VariableStruct {
    int mode;
    int hdrId_bias;
    bool isReapply;
    bool includeZeros;
    int startSample;
    int endSample;
//    int tmpCounter;
  };
  static int const MODE_ENSEMBLE = 11;
  static int const MODE_TRACE    = 12;
  static int const MODE_ALL      = 13;
}
using namespace mod_debias;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_debias_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader* shdr    = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  vars->hdrId_bias   = -1;
  vars->mode         = MODE_TRACE;
  vars->isReapply    = false;
  vars->includeZeros = true;
  vars->startSample = 0;
  vars->endSample = shdr->numSamples-1;

  edef->setExecType( EXEC_TYPE_MULTITRACE );

  if( param->exists( "mode" ) ) {
    std::string text;
    param->getString( "mode", &text );
    text = toLowerCase( text );
    if( !text.compare( "ensemble" ) ) {
      vars->mode = MODE_ENSEMBLE;
    }
    else if( !text.compare( "trace" ) ) {
      vars->mode = MODE_TRACE;
    }
    else {
      log->line("Unknown argument for user parameter 'mode': '%s'.", text.c_str());
      env->addError();
    }
  }
  if( param->exists( "reapply" ) ) {
    std::string text;
    param->getString( "reapply", &text );
    text = toLowerCase( text );
    if( !text.compare( "yes" ) ) {
      vars->isReapply = true;
    }
    else if( !text.compare( "no" ) ) {
      vars->isReapply = false;
    }
    else {
      log->line("Unknown argument for user parameter 'reapply': '%s'.", text.c_str());
      env->addError();
    }
  }
  if( param->exists( "zeros" ) ) {
    std::string text;
    param->getString( "zeros", &text );
    text = toLowerCase( text );
    if( !text.compare( "include" ) ) {
      vars->includeZeros = true;
    }
    else if( !text.compare( "exclude" ) ) {
      vars->includeZeros = false;
    }
    else {
      log->line("Unknown argument for user parameter 'zeros': '%s'.", text.c_str());
      env->addError();
    }
  }
  if( param->exists("window") ) {
    double start;
    double end;
    param->getDouble( "window", &start, 0 );
    param->getDouble( "window", &end, 1 );
    vars->startSample = (int)( start / shdr->sampleInt + 0.5 );
    vars->endSample   = (int)( end / shdr->sampleInt + 0.5 );

    if( vars->startSample < 0 ) vars->startSample = 0;
    if( vars->endSample >= shdr->numSamples ) vars->endSample = shdr->numSamples-1;
    if( vars->endSample <= vars->startSample ) log->error("Inconsistent start/end times in window: %f >=? %f\n", start, end);
  }

  if( vars->mode == MODE_ENSEMBLE ) {
    env->execPhaseDef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  }
  else if( vars->mode == MODE_TRACE ) {
    env->execPhaseDef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  }

  if( !hdef->headerExists("dc") ) {
    if( vars->isReapply ) {
      log->error("Trace header 'dc' not found: Required for removal of DC bias.");
    }
    hdef->addStandardHeader( "dc" );
  }
  vars->hdrId_bias = hdef->headerIndex("dc");
  
//  vars->tmpCounter = 0;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_debias_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    delete vars; vars = NULL;
    return;
  }

  int nTraces = traceGather->numTraces();
  
  if( !vars->isReapply ) {
    double mean_double = 0.0;
    for( int itrc = 0; itrc < nTraces; itrc++ ) {
      float* samples = traceGather->trace(itrc)->getTraceSamples();
      double sum = 0.0;
      int sampleCounter = 0;
      if( !vars->includeZeros ) {
        for( int isamp = vars->startSample; isamp <= vars->endSample; isamp++ ) {
          if( samples[isamp] != 0 ) {
            sampleCounter += 1;
            sum += (double)samples[isamp];
          }
        }
      }
      else {
        sampleCounter = vars->endSample - vars->startSample + 1;
        for( int isamp = vars->startSample; isamp <= vars->endSample; isamp++ ) {
          sum += (double)samples[isamp];
        }
      }
      if( sampleCounter != 0 ) mean_double += sum / (double)sampleCounter;
    }
    float mean_float = (float)( mean_double / (double)nTraces );
    
    for( int itrc = 0; itrc < nTraces; itrc++ ) {
      traceGather->trace(itrc)->getTraceHeader()->setFloatValue( vars->hdrId_bias, mean_float );
      float* samples = traceGather->trace(itrc)->getTraceSamples();
      if( !vars->includeZeros ) {
        for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
          if( samples[isamp] != 0 ) {
            samples[isamp] -= mean_float;
          }
        }
      }
      else {
        for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
          samples[isamp] -= mean_float;
        }
      }
    }
  }
  //-------------------------------------
  else {  // Re-apply DC value stored in trace header 'dc'
    float mean_float = 0.0;
    for( int itrc = 0; itrc < nTraces; itrc++ ) {
      mean_float = traceGather->trace(itrc)->getTraceHeader()->floatValue( vars->hdrId_bias );
      float* samples = traceGather->trace(itrc)->getTraceSamples();
      if( !vars->includeZeros ) {
        for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
          if( samples[isamp] != 0 ) {
            samples[isamp] += mean_float;
          }
        }
      }
      else {
        for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
          samples[isamp] += mean_float;
        }
      }
    }
  }

}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_debias_( csParamDef* pdef ) {
  pdef->setModule( "DEBIAS", "De-bias input data", "Remove average DC bias. Save DC bias value in trace header 'dc'.");
  
  pdef->addParam( "mode", "Mode of operation", NUM_VALUES_FIXED );
  pdef->addValue( "trace", VALTYPE_OPTION );
  pdef->addOption( "trace", "Remove average DC bias separately from each individual input trace." );
  pdef->addOption( "ensemble", "Remove average DC bias from whole input ensemble." );

  pdef->addParam( "reapply", "Re-apply DC bias?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Re-apply DC bias, e.g. add back DC bias stored in trace header 'dc'." );
  pdef->addOption( "no", "Do not re-apply DC bias (normal mode of operation)." );

  pdef->addParam( "zeros", "How shall zeros in data be handled?", NUM_VALUES_FIXED );
  pdef->addValue( "include", VALTYPE_OPTION );
  pdef->addOption( "exclude", "Exclude zeros from DC bias computation", "Zero values will not contribute to DC bias computation, and will remain unchanged in the output" );
  pdef->addOption( "include", "Include zeros in DC bias computation" );

  pdef->addParam( "window", "Computation window", NUM_VALUES_FIXED, "Full trace is used if not specified" );
  pdef->addValue( "", VALTYPE_NUMBER, "Start time" );
  pdef->addValue( "", VALTYPE_NUMBER, "End time" );
}

extern "C" void _params_mod_debias_( csParamDef* pdef ) {
  params_mod_debias_( pdef );
}
extern "C" void _init_mod_debias_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_debias_( param, env, log );
}
extern "C" void _exec_mod_debias_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_debias_( traceGather, port, numTrcToKeep, env, log );
}


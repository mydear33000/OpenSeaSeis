/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csTimeStretch.h"
#include "csTable.h"
#include <cmath>
#include <cstring>
#include <cstdio>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: TIME_STRETCH
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_time_stretch {
  struct VariableStruct {
    float* times;    // [s]
    float* values;  // [ms]
    csTimeStretch* timeStretchObj;
    int numTimes;
    float* traceBuffer;
    int mode;
    int valueType;

    int hdrId_stretch;

    csTable const* table; // Do not allocate or free!
    int numKeys;
    int* hdrId_keys;
    double* keyValueBuffer;
  };
  static int const TYPE_STRETCH = 1;
  static int const TYPE_TIME    = 2;
}
using mod_time_stretch::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_time_stretch_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

//---------------------------------------------
//
  vars->numTimes       = 0;
  vars->times          = NULL;
  vars->values          = NULL;
  vars->timeStretchObj = NULL;
  vars->traceBuffer    = NULL;
  vars->hdrId_stretch  = -1;
  vars->valueType      = mod_time_stretch::TYPE_STRETCH;

  vars->table          = NULL;
  vars->hdrId_keys     = NULL;
  vars->keyValueBuffer = NULL;
  vars->numKeys        = 0;
  vars->mode           = APPLY;

  csVector<std::string> valueList;

  //-------------------------------------
  std::string text;

  if( param->exists("mode") ) {
    param->getString( "mode", &text );
    if( !text.compare("apply") ) {
      vars->mode = APPLY;
    }
    else if( !text.compare("remove") ) {
      vars->mode = REMOVE;
    }
    else {
      log->error("Option not recognized: %s.", text.c_str());
    }
  }

  if( param->exists("value_type") ) {
    param->getString( "value_type", &text );
    if( !text.compare("stretch") ) {
      vars->valueType = mod_time_stretch::TYPE_STRETCH;
    }
    else if( !text.compare("time") ) {
      vars->valueType = mod_time_stretch::TYPE_TIME;
    }
    else {
      log->error("Option not recognized: %s.", text.c_str());
    }
  }

//---------------------------------------------
// Retrieve velocity table
//
/*
  if( param->exists("table") ) {
    param->getString("table", &text );
    vars->table = env->getTable( text );
    if( vars->table == NULL ) {
      log->error("No stretch table found with name '%s'. The table must be specified within the flow file with directive &table <tablename> <filename> {methodInterpolation}", text.c_str() );
    }
    vars->numKeys = vars->table->numKeys();
    if( vars->table->dimension() != csTable::TABLE_DIM_1D ) {
      log->error("Stretch table '%s' has a time column which makes it a 2D table.", text.c_str() );
    }
    vars->hdrId_keys     = new int[vars->numKeys];
    vars->keyValueBuffer = new double[vars->numKeys];
    for( int ikey = 0; ikey < vars->numKeys; ikey++ ) {
      if( !hdef->headerExists( vars->table->keyName(ikey) ) ) {
        log->error("No matching trace header found for table key '%s'", vars->table->keyName(ikey).c_str() );
      }
      vars->hdrId_keys[ikey] = hdef->headerIndex( vars->table->keyName(ikey) );
    }
  }
*/
//---------------------------------------------
// Retrieve time-stretch pairs
//
  if( param->exists( "time" ) ) {
    param->getAll( "time", &valueList );
    if( valueList.size() < 1 ){
      log->error("No times specified in user parameter 'time'");
    }
  }
  vars->numTimes = valueList.size();

  if( vars->numTimes > 0 ) {
    vars->times = new float[vars->numTimes];
    for( int i = 0; i < vars->numTimes; i++ ) {
      vars->times[i] = (float)atof( valueList.at(i).c_str() );
    }
  }
  else {
    vars->numTimes = 1;
    vars->times = new float[vars->numTimes];
    vars->times[0] = 0;
  }

  vars->values = new float[vars->numTimes];

  if( param->exists( "hdr_stretch" ) ) {
    if( vars->valueType != mod_time_stretch::TYPE_STRETCH ) {
      log->error("User parameter 'hdr_stretch' can only be specified in conjunction with value type 'stretch'");
    }
    param->getString( "hdr_stretch", &text );
    if( !hdef->headerExists( text ) ) {
      log->error("Unknown trace header: '%s'", text.c_str());
    }
    vars->hdrId_stretch = hdef->headerIndex( text );
  }
  else {  // No stretch header provided
    valueList.clear();
    param->getAll( "value", &valueList );
    int numValues = valueList.size();
    if( valueList.size() < 1 ) {
      log->error("No values supplied for user parameter 'value'");
    }
    else if( vars->valueType == mod_time_stretch::TYPE_STRETCH && numValues != vars->numTimes-1 ) {
      log->error("Non-matching number of stretch factors (=%d) and horizon times (=%d). Each stretch factor applies to one layer in between two horizon times, so the numberr of stretch factors should be one less",
        numValues, vars->numTimes );
    }
    else if( vars->valueType == mod_time_stretch::TYPE_TIME && numValues != vars->numTimes ) {
      log->error("Non-matching number of output (=%d) and input horizon times (=%d).", numValues, vars->numTimes );
    }
    float sign_scalar = ( vars->mode == APPLY ) ? 1.0 : -1.0;
    csFlexNumber number;
    for( int i = 0; i < numValues; i++ ) {
      if( !number.convertToNumber( valueList.at(i) ) ) {
        log->error("Specified value is not a valid number: '%s'", valueList.at(i).c_str() );
      }
      vars->values[i] = sign_scalar * number.floatValue();
      if( edef->isDebug() ) log->line("Value #%d: '%s' --> %f", i, valueList.at(i).c_str(), vars->values[i] );
    }
  }

  //---------------------------------------------
  // Set TimeStretch object
  //
  vars->timeStretchObj = new csTimeStretch( shdr->sampleInt, shdr->numSamples );
  vars->traceBuffer = new float[shdr->numSamples];

  if( vars->valueType == mod_time_stretch::TYPE_STRETCH ) {
    for( int i = 0; i < vars->numTimes-1; i++ ) {
      fprintf(stdout,"Horizon: %d  Top: %f  Bottom: %f   Stretch: %f\n", i, vars->times[i], vars->times[i+1], vars->values[i] );
    }
  }
  else {
    for( int i = 0; i < vars->numTimes; i++ ) {
      fprintf(stdout,"Horizon: %d  Input: %f   Output: %f\n", i, vars->times[i], vars->values[i] );
    }
  }
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_time_stretch_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    if( vars->times != NULL ) {
      delete [] vars->times;
      vars->times = NULL;
    }
    if( vars->values != NULL ) {
      delete [] vars->values;
      vars->values = NULL;
    }
    if( vars->timeStretchObj == NULL ) {
      delete vars->timeStretchObj;
      vars->timeStretchObj = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  float* samples = trace->getTraceSamples();
  memcpy( vars->traceBuffer, samples, shdr->numSamples*sizeof(float) );

  if( vars->hdrId_stretch > 0 ) {
    float s = trace->getTraceHeader()->floatValue( vars->hdrId_stretch );
    float sign_scalar = ( vars->mode == APPLY ) ? 1.0 : -1.0;
    for( int i = 0; i < vars->numTimes-1; i++ ) {
      vars->values[i] = sign_scalar * s;
      if( edef->isDebug() ) log->line("TIME_STRETCH DEBUG: Stretch #%d: %f", i, vars->values[i]);
    }
  }
  if( vars->valueType == mod_time_stretch::TYPE_STRETCH ) {
    vars->timeStretchObj->applyStretchFunction( vars->traceBuffer, vars->times, vars->values, vars->numTimes, samples );
  }
  else {
    vars->timeStretchObj->applyTimeInterval( vars->traceBuffer, vars->times, vars->values, vars->numTimes, samples );
  }
  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_time_stretch_( csParamDef* pdef ) {
  pdef->setModule( "TIME_STRETCH", "Stretch/squeeze data trace" );

  pdef->addDoc("To set up this module, specify either pairs of horizon input/output times, or otherwise the input horizon times and a stretch factor for each layer.");
  pdef->addDoc("Values can be specified either as ");

  pdef->addParam( "value_type", "Type of data provided to define stretch", NUM_VALUES_FIXED, "This value is provided by the usr parameter 'value', or through the trace header specified in 'hdr_stretch'" );
  pdef->addValue( "stretch", VALTYPE_OPTION );
  pdef->addOption( "stretch", "Provide stretch factor for each layer" );
  pdef->addOption( "time", "Provide output horizon time for each layer" );

  pdef->addParam( "time", "List of horizon times [ms]", NUM_VALUES_VARIABLE,
                  "For a model with N layers, specify N+1 time values. First value gives top of first layer, each subsequent value defines the bottom of one layer and the top of the next layer" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of horizon times [ms]." );

  pdef->addParam( "value", "List of stretch factors (or horizon output times) [ms]", NUM_VALUES_VARIABLE, "Number of stretch factors must be one less than number of times, one for each model layer. Number of output times must be the same as the number of input horizon times.");
  pdef->addValue( "", VALTYPE_NUMBER, "Stretch factor or output horizon time [ms]. Polarity for stretch factor: Stretch(+) or squeeze(-)" );

  pdef->addParam( "hdr_stretch", "Trace header containing stretch value", NUM_VALUES_FIXED, "This is an alternative to providing a list of values. This only applies to a single layer" );
  pdef->addValue( "", VALTYPE_STRING, "Name of trace header containing stretch value in [ms] (or [Hz])" );

  pdef->addParam( "mode", "Mode of application", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "apply", "Apply" );
  pdef->addOption( "remove", "Remove" );
}

extern "C" void _params_mod_time_stretch_( csParamDef* pdef ) {
  params_mod_time_stretch_( pdef );
}
extern "C" void _init_mod_time_stretch_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_time_stretch_( param, env, log );
}
extern "C" bool _exec_mod_time_stretch_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_time_stretch_( trace, port, env, log );
}


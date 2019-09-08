/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SCALING
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_scaling {
  struct VariableStruct {
    int option;
    int numTimes;
    float* sampIndex;
    float* scalar;
    int hdrId;
    bool isFile;
    bool isScalar;
    float* bufferScalar;
  };
}
using mod_scaling::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_scaling_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef  = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  
  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->option    = 0;
  vars->numTimes  = 0;
  vars->sampIndex = NULL;
  vars->scalar    = NULL;
  vars->hdrId     = -1;
  vars->bufferScalar = NULL;
  vars->isFile = false;
  vars->isScalar = false;

//---------------------------------------------
// Retrieve sampIndex and scalar
//
  csVector<std::string> valueList;

  if( param->exists( "filename" ) ) {
    std::string text;
    param->getString("filename",&text);
    csVector<float> timeList;
    csVector<float> scalarList;
    FILE* fin = fopen(text.c_str(),"r");
    if( fin == NULL ) log->error("Error occurred when opening input file %s.", text.c_str());
    char buffer[132];
    while( fgets(buffer,132,fin) != NULL ) {
      float time;
      float scalar;
      sscanf(buffer,"%f %f", &time, &scalar);
      timeList.insertEnd(time);
      scalarList.insertEnd(scalar);
      int num = timeList.size();
      if( num > 1 && timeList.at(num-1) <= timeList.at(num-2) ) log->error("Time #%d not monotonically increasing: %f(#%d) --> %f(#%d)",
                                                                           num, timeList.at(num-2), num-1, timeList.at(num-1), num);
    }
    fclose(fin);
    if( scalarList.size() == 0 ) log->error("Input file %s does not contain any data", text.c_str());
    int numSamples = scalarList.size();
    vars->bufferScalar = new float[shdr->numSamples];
    int indexCurrent = 0;
    float timeCurrent   = timeList.at(0);
    float timePrev      = timeCurrent;
    float scalarCurrent = scalarList.at(0);
    float scalarPrev    = scalarCurrent;
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      float time = (float)isamp * shdr->sampleInt;
      if( (indexCurrent == 0 && time <= timeCurrent) || (indexCurrent == numSamples-1) ) {
        vars->bufferScalar[isamp] = scalarCurrent;
      }
      else {
        if( time > timeCurrent ) {
          do {
            indexCurrent += 1;
            timeCurrent  = timeList.at(indexCurrent);  
          } while( time > timeCurrent && indexCurrent < numSamples-1 );
          timePrev   = timeList.at(indexCurrent-1);
          scalarPrev = scalarList.at(indexCurrent-1);
          scalarCurrent = scalarList.at(indexCurrent);
        }
        if( timePrev != timeCurrent ) {
          vars->bufferScalar[isamp] = scalarPrev + (time-timePrev)/(timeCurrent-timePrev) * (scalarCurrent - scalarPrev);
        }
        else {
          vars->bufferScalar[isamp] = scalarCurrent;
        }
      }
      if( edef->isDebug() ) fprintf(stdout,"%d %f %e\n", isamp, time, vars->bufferScalar[isamp]);
    }
    vars->isFile = true;

  }


  if( param->exists( "time" ) ) {
    param->getAll( "time", &valueList );
    if( valueList.size() < 1 ){
      log->error("No times specified in user parameter 'TIME'!");
    }
  }
  vars->numTimes = valueList.size();

  if( vars->numTimes > 0 ) {
    vars->sampIndex = new float[vars->numTimes];
    for( int i = 0; i < vars->numTimes; i++ ) {
      vars->sampIndex[i] = atof( valueList.at(i).c_str() ) / shdr->sampleInt;  // Convert to sample index
    }
  }
  else {
    vars->numTimes = 1;
    vars->sampIndex = new float[vars->numTimes];
    vars->sampIndex[0] = 0;
  }

  vars->scalar = new float[vars->numTimes];

  valueList.clear();

  if( param->exists("header") ) {
    string name;
    param->getString("header",&name);
    vars->hdrId = hdef->headerIndex( name );
  }
  else if( param->exists("scalar") ) {
    vars->isScalar = true; 
    param->getAll( "scalar", &valueList );

    if( valueList.size() < 1 ){
      log->error("Missing user parameter 'scalar'!");
    }
    else if( valueList.size() != vars->numTimes ) {
      log->error("Unequal number of scalar(%d) and times(%d)", valueList.size(), vars->numTimes );
    }
    csFlexNumber number;
    for( int i = 0; i < vars->numTimes; i++ ) {
//    vars->scalar[i] = atof( valueList.at(i).c_str() );
      if( !number.convertToNumber( valueList.at(i) ) ) {
        log->error("Specified scalar is not a valid number: '%s'", valueList.at(i).c_str() );
      }
      vars->scalar[i] = number.floatValue();
      if( edef->isDebug() ) log->line("Scalar #%d: '%s' --> %f", i, valueList.at(i).c_str(), vars->scalar[i] );
    }
  }
  else if( !vars->isFile ) {
    log->error("Need to specify at least one scaling parameter: 'header', 'scalar' or 'filename'");
  }
//----------------------------------------------
  
  std::string tableName;
  if( param->exists("table") ) {
    param->getString( "table", &tableName );
  }
//----------------------------------------------
  if( edef->isDebug() ) {
    for( int i = 0; i < vars->numTimes; i++ ) {
      log->line("sampleIndex: %8.2f , time: %8.2f, scalar: %8.4f", vars->sampIndex[i], vars->sampIndex[i]*shdr->sampleInt,
        vars->scalar[i] );
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_scaling_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    if( vars->sampIndex != NULL ) {
      delete [] vars->sampIndex; vars->sampIndex = NULL;
    }
    if( vars->scalar != NULL ) {
      delete [] vars->scalar; vars->scalar = NULL;
    }
    if( vars->bufferScalar != NULL ) {
      delete [] vars->bufferScalar;
      vars->bufferScalar = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  float* scalar = vars->scalar;
  float* sampIndex = vars->sampIndex;
  float* samples = trace->getTraceSamples();

  if( vars->hdrId >= 0 ) {
    double scalarValue = trace->getTraceHeader()->doubleValue(vars->hdrId);
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] *= scalarValue;
    }
  }
  else if( vars->isScalar ) {
    float sc;
    if( vars->numTimes == 1 ) {
      sc = scalar[0];
      for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] *= sc;
      }
    }
    else {
      int sampIndexFirst = (int)( sampIndex[0] ) + 1;
      int sampIndexLast  = (int)( sampIndex[vars->numTimes-1] + 0.01 );
      
      sc = scalar[0];
      for( int isamp = 0; isamp < sampIndexFirst; isamp++ ) {
        samples[isamp] *= sc;
      }
      sc = scalar[vars->numTimes-1];
      for( int isamp = sampIndexLast; isamp < shdr->numSamples; isamp++ ) {
        samples[isamp] *= sc;
      }
      
      int index = 1;
      float weight;
      float norm = 1.0;
      bool hasChanged = true;
      for( int isamp = sampIndexFirst; isamp < sampIndexLast; isamp++ ) {
        while( isamp > sampIndex[index] ) {
          hasChanged = true;
          index++;
        }
        if( hasChanged ) {
          norm = sampIndex[index] - sampIndex[index-1];
          hasChanged = false;
        }
        weight = ((float)isamp-sampIndex[index-1])/norm;
        sc = ( (1.0-weight) * scalar[index-1] + weight * scalar[index] );
        samples[isamp] *= sc;
      }
    }
  }
  if( vars->isFile ) {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] *= vars->bufferScalar[isamp];
    }
  }

  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_scaling_( csParamDef* pdef ) {
  pdef->setModule( "SCALING", "Scale trace data with linear function", "Apply time variant linear scaling function to seismic trace data" );

  pdef->setVersion( 1, 0 );

  pdef->addParam( "filename", "Input ASCII file name, containing table of pairs of time/scalar", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "time", "List of time value [ms]", NUM_VALUES_VARIABLE, "Time knee points at which specified scalars apply" );
  pdef->addValue( "", VALTYPE_NUMBER, "List of time values [ms]..." );

  pdef->addParam( "scalar", "List of scalars", NUM_VALUES_VARIABLE,
      "Scalar values that apply at specified time knee points. In between time knee points, scalar is linearly interpolated." );
  pdef->addValue( "", VALTYPE_NUMBER, "List of scalars..." );

  pdef->addParam( "header", "Trace header name containing scalar value", NUM_VALUES_FIXED);
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
}

extern "C" void _params_mod_scaling_( csParamDef* pdef ) {
  params_mod_scaling_( pdef );
}
extern "C" void _init_mod_scaling_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_scaling_( param, env, log );
}
extern "C" bool _exec_mod_scaling_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_scaling_( trace, port, env, log );
}


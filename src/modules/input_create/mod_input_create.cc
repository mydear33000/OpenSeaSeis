/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csStandardHeaders.h"
#include <cmath>
#include <cstring>
#include <ctime>    // For time
#include <cstdlib>  // For srand() and rand()

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: input_create
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_input_create {
  struct Receiver {
    int p;
    double x;
    double y;
  };
  struct VariableStruct {
    float value;
    int numTraces;
    bool atEOF;

    int numSpikes;
    int* samplesSpikes;
    float* values;

    float maxNoise;

    int hdrId_trcno;
    int hdrId_rcv;
    int hdrId_rec_x;
    int hdrId_rec_y;

    int traceCounter;

    // Plane wave:
    float pw_azimuth;  // in [rad]
    float pw_slowness; // s/m
    float pw_reftime_ms;
    float pw_x0;
    float pw_y0;
    int numPlaneWaves;

    // Point source:
    float* ps_x0;
    float* ps_y0;
    float ps_slowness; // s/m
    float ps_reftime_ms;
    float ps_spreadingFactor;
    int numPointSources;

    Receiver* receiver;
  };
}
using mod_input_create::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_input_create_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new mod_input_create::VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->hdrId_trcno  = -1;
  vars->hdrId_rcv    = -1;
  vars->hdrId_rec_x  = -1;
  vars->hdrId_rec_y  = -1;
  vars->atEOF        = false;
  vars->traceCounter = 0;
  vars->numTraces    = 1;
  vars->value        = 1.0f;
  vars->numSpikes = 0;
  vars->samplesSpikes = NULL;
  vars->values      = NULL;
  vars->maxNoise    = 0.0;

  vars->pw_azimuth  = 0.0;
  vars->pw_slowness = 0.0;
  vars->pw_reftime_ms = 0.0;
  vars->pw_x0       = 0.0;
  vars->pw_y0       = 0.0;
  vars->numPlaneWaves = 0;

  vars->ps_slowness = 0.0;
  vars->ps_reftime_ms = 0.0;
  vars->ps_x0       = NULL;
  vars->ps_y0       = NULL;
  vars->ps_spreadingFactor = 0;
  vars->numPointSources = 0;

  vars->receiver = NULL;

  //----------------------------------------------------
  double length;
  float sampleInt;

  param->getFloat( "sample_int", &sampleInt, 0 );
  shdr->sampleInt = sampleInt;

  string unit("ms");
  param->getDouble( "length", &length, 0 );
  if( param->getNumValues("length") > 1 ) {
    param->getString( "length", &unit, 1 );
  }
  bool isValue = false;
  if( param->exists("value") ) {
    param->getFloat( "value", &vars->value );
    isValue = true;
  }

  if( !unit.compare("samples") ) {
    shdr->domain     = DOMAIN_XT;
    shdr->numSamples = (int)(length + 0.5);
  }
  else if( !unit.compare("s") || !unit.compare("seconds") ) {
    shdr->domain     = DOMAIN_XT;
    shdr->numSamples = (int)( 1000*(length/shdr->sampleInt) + 0.5 );
  }
  else if( !unit.compare("ms") || !unit.compare("milliseconds") ) {
    shdr->domain     = DOMAIN_XT;
    shdr->numSamples = (int)( (length/shdr->sampleInt) + 0.5 );
  }
  else if( !unit.compare("m") || !unit.compare("meters") ) {
    shdr->domain     = DOMAIN_XD;
    shdr->numSamples = (int)( (length/shdr->sampleInt) + 0.5 );
  }
  else if( !unit.compare("km") || !unit.compare("kilometers") ) {
    shdr->domain     = DOMAIN_XD;
    shdr->numSamples = (int)( 1000*(length/shdr->sampleInt) + 0.5 );
  }
  else if( !unit.compare("Hz") ) {
    shdr->domain     = DOMAIN_FX;
    shdr->numSamples = (int)( (length/shdr->sampleInt) + 0.5 );
  }
  else {
    log->error("Unknown option: %s", unit.c_str() );
  }

  //-------------------------------------------------------------------------
  if( param->exists("noise") ) {
    param->getFloat("noise", &vars->maxNoise);
  }
  //-------------------------------------------------------------------------
  if( param->exists( "plane_wave" ) ) {
    int iv = 0;
    vars->numPlaneWaves = 1;
    param->getFloat( "plane_wave", &vars->pw_azimuth, iv++ );
    param->getFloat( "plane_wave", &vars->pw_slowness, iv++ );
    vars->pw_slowness /= 1000.0f;   // Convert from [s/km] to [s/m]
    param->getFloat( "plane_wave", &vars->pw_reftime_ms, iv++ );
    param->getFloat( "plane_wave", &vars->pw_x0, iv++ );
    param->getFloat( "plane_wave", &vars->pw_y0, iv++ );
    vars->pw_azimuth *= (float)(M_PI / 180.0);
    if( vars->value == 0.0 ) {
      log->error("When computing plane waves, a value (parameter 'value') other than 0 must be specified. Plane wave arrivals will be scaled by this value.");
    }
  }
  //-------------------------------------------------------------------------
  if( param->exists( "point_source" ) ) {
    vars->numPointSources = param->getNumLines("point_source_xy");
    if( vars->numPointSources == 0 ) {
      log->error("Specify at least one point source location with parameter 'point_source_xy'");
    }
    vars->ps_x0 = new float[vars->numPointSources];
    vars->ps_y0 = new float[vars->numPointSources];
    for( int i = 0; i < vars->numPointSources; i++ ) {
      param->getFloatAtLine( "point_source_xy", &vars->ps_x0[i], i, 0 );
      param->getFloatAtLine( "point_source_xy", &vars->ps_y0[i], i, 1 );
    }
    param->getFloat( "point_source", &vars->ps_slowness, 0 );
    vars->ps_slowness /= 1000.0f;   // Convert from [s/km] to [s/m]
    param->getFloat( "point_source", &vars->ps_reftime_ms, 1 );
    param->getFloat( "point_source", &vars->ps_spreadingFactor, 2 );
    if( vars->value == 0.0 ) {
      log->error("When computing point sources, a value (parameter 'value') other than 0 must be specified. Arrivals will be scaled by this value.");
    }
  }
  //-------------------------------------------------------------------------
  if( param->exists("spikes") ) {
    vars->numSpikes = param->getNumValues("spikes");
    vars->samplesSpikes = new int[vars->numSpikes];
    for( int i = 0; i < vars->numSpikes; i++ ) {
      float time;
      //      param->getInt("spikes", &vars->timesSpikes[i], i);
      param->getFloat("spikes", &time, i);
      vars->samplesSpikes[i] = (int)( time / shdr->sampleInt + 0.5);
    }

    vars->values = new float[vars->numSpikes];
    if( param->exists("values") ) {
      if( isValue ) log->error("Parameter 'values' cannot be specified in conjunction with parameter 'value'. Specify only one of these two");
      int numValues = param->getNumValues("values");
      if( numValues != vars->numSpikes ) log->error("Parameter 'values': Inconsistent number of values (=%d). Should match number of spikes (=%d)", numValues, vars->numSpikes);
      for( int i = 0; i < vars->numSpikes; i++ ) {
        param->getFloat("values", &vars->values[i], i);
      }
    }
    else {
      for( int i = 0; i < vars->numSpikes; i++ ) {
        vars->values[i] = vars->value;
      }
    }

  }

  //-------------------------------------------------------------------------
  if( param->exists("ntraces") ) {
    param->getInt( "ntraces", &vars->numTraces );
  }

  //-------------------------------------------------------------------------
  //
  if( param->exists("rec_geom") || vars->numPlaneWaves > 0 ) {
    std::string text;
    param->getString( "rec_geom", &text );
    FILE* fin = fopen(text.c_str(),"r");
    if( fin == NULL ) {
      log->error("Problems opening receiver geometry file: %s", text.c_str() );
    }
    csVector<mod_input_create::Receiver> receiverList;
    char line[132];
    while( fgets( line, 132, fin ) != NULL ) {
      mod_input_create::Receiver rcv;
      sscanf(line,"%d %lf %lf", &rcv.p, &rcv.x, &rcv.y);
      receiverList.insert(rcv);
    }
    vars->receiver = new mod_input_create::Receiver[receiverList.size()];
    for( int i = 0; i < receiverList.size(); i++ ) {
      vars->receiver[i].p = receiverList.at(i).p;
      vars->receiver[i].x = receiverList.at(i).x;
      vars->receiver[i].y = receiverList.at(i).y;
    }
    if( param->exists("ntraces") && vars->numTraces != receiverList.size() ) {
      log->warning("Receiver geometry overrules specified number of traces (=%d): Set number of traces to %d.", vars->numTraces, receiverList.size());
    }
    vars->numTraces = receiverList.size();
    vars->hdrId_rcv   = hdef->addHeader( &cseis_geolib::HDR_RCV );
    vars->hdrId_rec_x = hdef->addHeader( &cseis_geolib::HDR_REC_X );
    vars->hdrId_rec_y = hdef->addHeader( &cseis_geolib::HDR_REC_Y );

    if( edef->isDebug() ) {
      for( int i = 0; i < receiverList.size(); i++ ) {
        log->line("#%2d   %d %lf %lf", i+1, vars->receiver[i].p, vars->receiver[i].x, vars->receiver[i].y);
      }
    }
  }
  //-------------------------------------------------------------------------
  //
  vars->hdrId_trcno = hdef->addHeader( csStandardHeaders::get("trcno") );
  vars->traceCounter = 0;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_input_create_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const*  shdr = env->superHeader;
//  csTraceHeaderDef const*  hdef = env->headerDef;

  if( edef->isCleanup() ) {
    if( vars->samplesSpikes != NULL ) {
      delete [] vars->samplesSpikes; vars->samplesSpikes = NULL;
    }
    if( vars->receiver != NULL ) {
      delete [] vars->receiver;
      vars->receiver = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  if( vars->atEOF ) return false;

  if( vars->traceCounter == vars->numTraces-1 ) {
    vars->atEOF = true;
  }

  float* samples = trace->getTraceSamples();

  //-------------------------------------------------------------------------------------------
  // Add spikes
  //
  if( vars->samplesSpikes != NULL ) {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] = 0.0;
    }
    for( int ispike = 0; ispike < vars->numSpikes; ispike++ ) {
      if( vars->samplesSpikes[ispike] < shdr->numSamples && vars->samplesSpikes[ispike] >= 0 ) {
        samples[ vars->samplesSpikes[ispike] ] = vars->values[ispike];
      }
    }
  }
  else if( vars->numPlaneWaves > 0 || vars->numPointSources > 0 ) {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] = 0.0;
    }
    //-------------------------------------------------------------------------------------------
    // Add plane wave
    //
    if( vars->numPlaneWaves > 0 ) {
      double dx = vars->receiver[vars->traceCounter].x - vars->pw_x0;
      double dy = vars->receiver[vars->traceCounter].y - vars->pw_y0;
      double projection = dx * sin(vars->pw_azimuth) + dy * cos(vars->pw_azimuth);
      float delayTime_ms = (float)( projection * vars->pw_slowness ) * 1000.0; // [ms]
      int nearestSampleIndex = (int)round( (delayTime_ms+vars->pw_reftime_ms) / shdr->sampleInt );
      if( nearestSampleIndex >= 0 && nearestSampleIndex < shdr->numSamples ) {
        samples[nearestSampleIndex] += vars->value;
      }
      if( edef->isDebug() ) {
        log->line("Receiver %3d  %12f %12f  delay time  %10f ms (+%fms reference time)  sample_index %d",
          vars->receiver[vars->traceCounter].p, vars->receiver[vars->traceCounter].x, vars->receiver[vars->traceCounter].y,
          delayTime_ms, vars->pw_reftime_ms, nearestSampleIndex);
      }
    }
    //-------------------------------------------------------------------------------------------
    // Add point source
    //
    for( int ip = 0; ip < vars->numPointSources; ip++ ) {
      double dx = vars->receiver[vars->traceCounter].x - vars->ps_x0[ip];
      double dy = vars->receiver[vars->traceCounter].y - vars->ps_y0[ip];
      double offset = sqrt( dx * dx + dy * dy );
      float delayTime_ms = (float)( offset * vars->ps_slowness * 1000.0 ); // [ms]
      int nearestSampleIndex = (int)round( (delayTime_ms+vars->ps_reftime_ms) / shdr->sampleInt );
      if( nearestSampleIndex >= 0 && nearestSampleIndex < shdr->numSamples ) {
        samples[nearestSampleIndex] += vars->value / (float)max(1.0,pow(offset,(double)vars->ps_spreadingFactor));
      }
      if( edef->isDebug() ) {
        log->line("Receiver %3d  %12f %12f  delay time  %10f ms (+%fms reference time)  sample_index %d, offset: %f",
          vars->receiver[vars->traceCounter].p, vars->receiver[vars->traceCounter].x, vars->receiver[vars->traceCounter].y,
          delayTime_ms, vars->ps_reftime_ms, nearestSampleIndex, offset);
      }
    }
  }
  //-------------------------------------------------------------------------------------------
  // Set to constant value
  //
  else {
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      samples[isamp] = vars->value;
    }
  }
  //-------------------------------------------------------------------------------------------
  // Add noise
  if( vars->maxNoise != 0.0 ) {
    srand( (vars->traceCounter+1)*shdr->numSamples );
    for( int isamp = 0; isamp < shdr->numSamples; isamp++ ) {
      double random  = 2.0 * ( (double)rand() - 0.5*(double)RAND_MAX ) / (double)RAND_MAX;
      samples[isamp] += (float)random * vars->maxNoise;
    }
  }
  csTraceHeader* trcHdr = trace->getTraceHeader();
  trcHdr->setIntValue( vars->hdrId_trcno, vars->traceCounter+1 );

  if( vars->hdrId_rcv >= 0 ) {
    trcHdr->setIntValue( vars->hdrId_rcv, vars->receiver[vars->traceCounter].p );
    trcHdr->setDoubleValue( vars->hdrId_rec_x, vars->receiver[vars->traceCounter].x );
    trcHdr->setDoubleValue( vars->hdrId_rec_y, vars->receiver[vars->traceCounter].y );
  }

  vars->traceCounter += 1;

  return true;
}
//********************************************************************************
// Parameter definition
//
//
//********************************************************************************
void params_mod_input_create_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_CREATE", "Create traces", "Simple synthetic trace generator" );

  pdef->addParam( "length", "Trace length", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Length of traces" );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );
  pdef->addOption( "samples", "Samples" );
  pdef->addOption( "m", "Meters" );
  pdef->addOption( "km", "Kilometers" );

  pdef->addParam( "ntraces", "Number of traces", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Number of traces" );

  pdef->addParam( "sample_int", "Sample interval", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Sample interval, unit specified under user parameter 'length'" );

  pdef->addParam( "value", "Constant value to set trace samples", NUM_VALUES_FIXED );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Set trace samples to this value" );

  pdef->addParam( "noise", "Add Gaussian random noise", NUM_VALUES_FIXED );
  pdef->addValue( "1.0", VALTYPE_NUMBER, "Maximum amplitude of random noise" );

  pdef->addParam( "spikes", "Add spikes", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "List of times/depths/frequencies (unit see user parameter 'length') where a spike shall be added (rest of trace is set to 0)" );

  pdef->addParam( "values", "List of values corresponding to spikes", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Set each spike to the corresponding value (use parameter 'value' to set each spike to the same value)" );

  pdef->addParam( "plane_wave", "Generate plane wave (spikes only)", NUM_VALUES_FIXED, "Requires the definiton of a receiver geometry" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Azimuth [deg]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Slowness [s/km]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Reference time [ms]: Time to add to computed delay time" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Origin X coordinate [m] for delay time computation" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Origin Y coordinate [m] for delay time computation" );

  pdef->addParam( "point_source", "Generate point source (spikes only)", NUM_VALUES_FIXED, "Requires the definiton of a receiver geometry" );
//  pdef->addValue( "0.0", VALTYPE_NUMBER, "Point source Z coordinate [m]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Slowness [s/km]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Reference time [ms]: Time to add to computed travel time" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Geometric spreading 1/r^N" );

  pdef->addParam( "point_source_xy", "Point source location", NUM_VALUES_FIXED );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Point source X coordinate [m]" );
  pdef->addValue( "0.0", VALTYPE_NUMBER, "Point source Y coordinate [m]" );

  pdef->addParam( "rec_geom", "Receiver geometry definition", NUM_VALUES_FIXED, "Specify ASCII file name containing listing of receiver XY positions" );
  pdef->addValue( "", VALTYPE_STRING, "ASCII file name. Format:  rcv  rec_x  rec_y" );
}

extern "C" void _params_mod_input_create_( csParamDef* pdef ) {
  params_mod_input_create_( pdef );
}
extern "C" void _init_mod_input_create_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_input_create_( param, env, log );
}
extern "C" bool _exec_mod_input_create_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_input_create_( trace, port, env, log );
}



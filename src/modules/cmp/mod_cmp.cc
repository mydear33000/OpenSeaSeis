/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "geolib_methods.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: CMP
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_cmp {
  struct VariableStruct {
    float depth;
    float bin_scalar;
    int method;
    int hdrId_offset;
    int hdrId_cmp;
    int hdrId_rcv;
    int hdrId_source;
    int hdrId_sou_z;
    int hdrId_rec_elev;
  };
  static int const METHOD_NORMAL   = 1;
  static int const METHOD_OBC_ASYM = 2;
}
using namespace mod_cmp;

//*************************************************************************************************
// Init phase
//
//
//*************************************************************************************************
void init_mod_cmp_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->depth      = 0;
  vars->bin_scalar = 2.0;
  vars->method     = METHOD_NORMAL;
  vars->hdrId_offset   = -1;
  vars->hdrId_cmp      = -1;
  vars->hdrId_rcv      = -1;
  vars->hdrId_source   = -1;
  vars->hdrId_sou_z    = -1;
  vars->hdrId_rec_elev = -1;

  //-------------------------------------------------------------
  //
  std::string text;
  if( param->exists("method") ) {
    param->getString("method",&text);
    if( !text.compare("normal") ) {
      vars->method = METHOD_NORMAL;
    }
    else if( !text.compare("obc_asym") ) {
      vars->method = METHOD_OBC_ASYM;
      param->getFloat( "depth", &vars->depth );
    }
    else {
      log->error("Unknown option: '%s'", text.c_str());
    }
  }

  if( param->exists("bin_scalar") ) {
    param->getFloat("bin_scalar",&vars->bin_scalar);
  }

  //-------------------------------------------------------------
  // Headers
  //
  if( !hdef->headerExists( "source" ) ) {
    log->error("Trace header 'source' does not exist.");
  }
  if( !hdef->headerExists( "rcv" ) ) {
    log->error("Trace header 'rcv' does not exist.");
  }

  if( vars->method == METHOD_OBC_ASYM ) {
    if( !hdef->headerExists( "offset" ) ) {
      log->error("Trace header 'offset' does not exist.");
    }
    else if( hdef->headerType( "offset" ) != TYPE_FLOAT && hdef->headerType( "offset" ) != TYPE_DOUBLE ) {
      log->error("Trace header 'offset' exists but has the wrong number type. Should be FLOAT.");
    }
    if( !hdef->headerExists( "rec_elev" ) ) {
      log->error("Trace header 'rec_elev' does not exist.");
    }
    if( !hdef->headerExists( "sou_z" ) ) {
      log->error("Trace header 'sou_z' does not exist.");
    }
    vars->hdrId_offset = hdef->headerIndex( "offset" );
    vars->hdrId_rec_elev = hdef->headerIndex( "rec_elev" );
    vars->hdrId_sou_z    = hdef->headerIndex( "sou_z" );
  }
  vars->hdrId_rcv    = hdef->headerIndex( "rcv" );
  vars->hdrId_source = hdef->headerIndex( "source" );

  if( !hdef->headerExists( "cmp" ) ) {
    hdef->addStandardHeader( "cmp" );
  }
  vars->hdrId_cmp = hdef->headerIndex( "cmp" );

  // cmp scalar = 1: cmp numbering follows that of source and receiver stations
  // cmp scalar = 2: cmp numbering is twice that of source and receiver stations
  // ...
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_cmp_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    delete vars; vars = NULL;
    return true;
  }

  csTraceHeader* trcHdr = trace->getTraceHeader();
  double rcv     = (double)trcHdr->intValue( vars->hdrId_rcv );
  double source  = (double)trcHdr->intValue( vars->hdrId_source );
  int cmp;

  if( vars->method == METHOD_NORMAL ) {
    cmp = (int)( vars->bin_scalar * 0.5 * ( source + rcv ) + 0.5 );
  }
  else {
    double rec_elev = trcHdr->doubleValue( vars->hdrId_rec_elev );
    double sou_z    = trcHdr->doubleValue( vars->hdrId_sou_z );
    double offset   = trcHdr->doubleValue( vars->hdrId_offset );

    double dz       = rec_elev - sou_z;
    if( dz > vars->depth ) {
      log->warning("Specified target depth (%f) is shallower than receiver depth (%f), rcv: %d, source: %d",
        vars->depth, rec_elev, rcv, source );
    }
    double offset_pside = offset / ( 2.0 - dz/vars->depth );
    cmp = (int)( vars->bin_scalar * ( source + (offset_pside/offset)*( rcv - source ) ) + 0.5 );
  }

  trcHdr->setIntValue( vars->hdrId_cmp, cmp );

  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_cmp_( csParamDef* pdef ) {
  pdef->setModule( "CMP", "Perform CMP binning" );

  pdef->addParam( "method", "CMP binning method", NUM_VALUES_FIXED );
  pdef->addValue( "normal", VALTYPE_OPTION );
  pdef->addOption( "normal", "'Normal' CMP binning. Set CMP at centre between source and receiver." );
  pdef->addOption( "obc_asym", "Asymmetric CMP binning for OBC geometry before redatuming. This method requires a target depth at which CMP binning will be correct" );

  pdef->addParam( "bin_scalar", "Numbering of CMP: Apply scalar to 'normal' CMP number", NUM_VALUES_FIXED );
  pdef->addValue( "2", VALTYPE_NUMBER, "CMP numbering is computed as follows: cmp = bin_scalar * (source + rcv)/2" );

  pdef->addParam( "depth", "Target depth at which cmp binning shall be performed", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Target depth [m]" );
}

extern "C" void _params_mod_cmp_( csParamDef* pdef ) {
  params_mod_cmp_( pdef );
}
extern "C" void _init_mod_cmp_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_cmp_( param, env, log );
}
extern "C" bool _exec_mod_cmp_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_cmp_( trace, port, env, log );
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"
#include "csSortManager.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SORT
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_sort {
  struct VariableStruct {
    int numHeaders;
    int* indexHdr;
    cseis_geolib::type_t* hdrTypes;
    std::string* hdrNames;
    int* sortDir;
    int traceCounter;
    cseis_geolib::csSortManager* sortManager;
  };
  static int const INCREASING = 1;
  static int const DECREASING = 2;
}
using namespace mod_sort;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_sort_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
//  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_MULTITRACE );
  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->numHeaders = 0;
  vars->indexHdr = NULL;
  vars->hdrTypes  = NULL;
  vars->hdrNames = NULL;
  vars->sortDir  = NULL;
  vars->traceCounter = 0;
  vars->sortManager  = NULL;

  //-----------------------------------------

  int nLines = param->getNumLines( "header" );
//  if( nLines > 1 ) {
//    log->error( "More than one line encountered for user parameter HEADER. Only one line is supported at the moment." );
//  }
  vars->numHeaders = nLines;

  csVector<string> valueList;
  param->getAll( "header", &valueList, 0 );
  if( valueList.size() < 1 ) {
    log->line("Wrong number of arguments for user parameter 'header'. Expected: 2, found: %d.", valueList.size());
    env->addError();
  }
  else {
    vars->sortDir   = new int[vars->numHeaders];
    vars->indexHdr  = new int[vars->numHeaders];
    vars->hdrTypes = new type_t[vars->numHeaders];
    vars->hdrNames  = new std::string[vars->numHeaders];
    valueList.clear();
    for( int ihdr = 0; ihdr < vars->numHeaders; ihdr++ ) {
      param->getAll( "header", &valueList, ihdr );
      string name = valueList.at(0);
      vars->hdrNames[ihdr] = name;
      if( valueList.size() > 1 ) {
        string sortDirText = valueList.at(1);
        if( !sortDirText.compare("increasing") ) {
          vars->sortDir[ihdr] = INCREASING;
        }
        else {
          vars->sortDir[ihdr] = DECREASING;
        }
      }
      else {
        vars->sortDir[ihdr] = INCREASING;
      }
      if( hdef->headerExists( name.c_str() ) ) {
        vars->indexHdr[ihdr]  = hdef->headerIndex( name.c_str() );
        vars->hdrTypes[ihdr] = hdef->headerType( name.c_str() );
        if( vars->hdrTypes[ihdr] == TYPE_STRING ) {
          log->line("String headers are not supported: %s", name.c_str());
          env->addError();
        }
      }
      else {
        log->line("Unknown trace header: %s", name.c_str());
        env->addError();
      }
    }
  }
  int sortMethod = csSortManager::SIMPLE_SORT;
  if( param->exists("method") ) {
    std::string text;
    param->getString("method", &text);
    if( !text.compare("simple") ) {
      sortMethod = csSortManager::SIMPLE_SORT;
    }
    else if( !text.compare("tree") ) {
      sortMethod = csSortManager::TREE_SORT;
    }
    else {
      log->error("Unknown option: '%s'", text.c_str() );
    }
  }

  vars->sortManager = new csSortManager( vars->numHeaders, sortMethod );

  vars->traceCounter = 0;
}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_sort_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
//  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    if( vars->indexHdr ) {
      delete [] vars->indexHdr;
      vars->indexHdr = NULL;
    }
    if( vars->hdrTypes ) {
      delete [] vars->hdrTypes;
      vars->hdrTypes = NULL;
    }
    if( vars->hdrNames ) {
      delete [] vars->hdrNames;
      vars->hdrNames = NULL;
    }
    if( vars->sortDir ) {
      delete [] vars->sortDir;
      vars->sortDir = NULL;
    }
    if( vars->sortManager != NULL ) {
      delete vars->sortManager;
      vars->sortManager = NULL;
    }
    delete vars; vars = NULL;
    return;
  }

//-------------------------------------------
  int nTraces = traceGather->numTraces();
  if( edef->isDebug() ) log->line("Number of input traces: %d", nTraces );
  vars->traceCounter += nTraces;

  vars->sortManager->resetValues( nTraces );
  for( int ihdr = 0; ihdr < vars->numHeaders; ihdr++ ) {
    int sign = ( vars->sortDir[ihdr] == INCREASING ) ? 1 : -1;
    if( vars->hdrTypes[ihdr] == TYPE_INT ) {
      for( int itrc = 0; itrc < nTraces; itrc++ ) {
        int value = sign*traceGather->trace(itrc)->getTraceHeader()->intValue(vars->indexHdr[ihdr]);
        vars->sortManager->setValue( itrc, vars->numHeaders-ihdr-1, csFlexNumber(value) );
      }
    }
    else if( vars->hdrTypes[ihdr] == TYPE_DOUBLE ) {
      for( int itrc = 0; itrc < nTraces; itrc++ ) {
        double value = (double)sign*traceGather->trace(itrc)->getTraceHeader()->doubleValue(vars->indexHdr[ihdr]);
        vars->sortManager->setValue( itrc, vars->numHeaders-ihdr-1, csFlexNumber(value) );
      }
    }
    else { // INT64
      for( int itrc = 0; itrc < nTraces; itrc++ ) {
        csInt64_t value = sign*traceGather->trace(itrc)->getTraceHeader()->int64Value(vars->indexHdr[ihdr]);
        vars->sortManager->setValue( itrc, vars->numHeaders-ihdr-1, csFlexNumber(value) );
      }
    }
  }
  //  fprintf(stdout,"START sorting %d ( =? %d ) traces...\n", nTraces, vars->sortManager->numValues());
  //  vars->sortManager->dump();
  //  fflush(stdout);
  vars->sortManager->sort();
  // fprintf(stdout,"END sort...\n");
  // fflush(stdout);

  csTrace** tracePtr = new csTrace*[nTraces];
  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    tracePtr[itrc] = traceGather->trace( vars->sortManager->sortedIndex(itrc) );
  }
  for( int itrc = 0; itrc < nTraces; itrc++ ) {
    (*traceGather)[itrc] = tracePtr[itrc];
  }
  delete [] tracePtr;
}


//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_sort_( csParamDef* pdef ) {
  pdef->setModule( "SORT", "Sort traces", "Very simple implementation. Requires to hold all traces to be sorted in memory at once." );

  pdef->addParam( "mode", "Sort mode", NUM_VALUES_FIXED );
  pdef->addValue( "ensemble", VALTYPE_OPTION );
  pdef->addOption( "ensemble", "Sort all traces in input ensemble" );

  pdef->addParam( "header", "Trace header to sort on", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addValue( "increasing", VALTYPE_OPTION );
  pdef->addOption( "increasing", "Sort in increasing order" );
  pdef->addOption( "decreasing", "Sort in decreasing order" );

  pdef->addParam( "method", "Sort method", NUM_VALUES_FIXED );
  pdef->addValue( "simple", VALTYPE_OPTION );
  pdef->addOption( "simple", "Simplest sort method. Fastest for small and partially pre-sorted data sets" );
  pdef->addOption( "tree", "Tree sorting method. Most efficient for large, totally un-sorted data sets" );
}


extern "C" void _params_mod_sort_( csParamDef* pdef ) {
  params_mod_sort_( pdef );
}
extern "C" void _init_mod_sort_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_sort_( param, env, log );
}
extern "C" void _exec_mod_sort_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_sort_( traceGather, port, numTrcToKeep, env, log );
}


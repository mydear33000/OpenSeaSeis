/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"

#ifdef PLATFORM_LINUX
#include "csTimer.h"
#include "libmatlabTest.h"
#endif

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: MATLAB_TEST
 *
 * @author Bjorn Olofsson
 * @date   2008
 */
namespace mod_matlab_test {
  struct VariableStruct {
#ifdef PLATFORM_LINUX
    mxArray* matlab_data_in;
    mxArray* matlab_data_out;
#endif
    std::string method_name;
  };
}
using mod_matlab_test::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_matlab_test_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*  shdr   = env->superHeader;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

#ifdef PLATFORM_LINUX
  vars->matlab_data_in  = NULL;
  vars->matlab_data_out = NULL;

  //--------------------------------------
  //
  //  param->getString("method_name", &vars->method_name);
  

  log->line(" Start up Matlab Runtime Environment...");
  fprintf(stdout," Start up Matlab Runtime Environment...\n");
  csTimer timer;
  timer.start();
  // Initialize Matlab Compiler Runtime (MCR) Environment
  mclmcrInitialize();

  log->line(" Time taken for Matlab startup (Runtime Environment):   %12.6fs\n", timer.getElapsedTime() );
  fprintf(stdout," Time taken for Matlab startup (Runtime Environment):   %12.6fs\n", timer.getElapsedTime() );

  // Initialize linked Matlab method...
  if( !libmatlabTestInitialize() ) {
    log->error("Could not initialize Matlab library.");
  }

  log->line(" Time taken for Matlab startup (Runtime + Module Init): %12.6fs\n", timer.getElapsedTime() );
  fprintf(stdout," Time taken for Matlab startup (Runtime + Module Init): %12.6fs\n", timer.getElapsedTime() );

  // Create Matlab objects
  vars->matlab_data_in = mxCreateNumericMatrix( shdr->numSamples, 1, mxSINGLE_CLASS, mxREAL );  // Create float array (mxSINGLE_CLASS)
  //  vars->matlab_data_in  = mxCreateDoubleMatrix( shdr->numSamples, 1, mxREAL );  // double
#endif

#ifdef PLATFORM_WINDOWS
  log->error("Matlab test module only available on Linux platform");
#endif

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_matlab_test_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

#ifdef PLATFORM_LINUX
  if( edef->isCleanup()){
    if( vars->matlab_data_in != NULL ) {
      mxDestroyArray( vars->matlab_data_in );
      vars->matlab_data_in = NULL;
    }
    if( vars->matlab_data_out != NULL ) {
      mxDestroyArray( vars->matlab_data_out );
      vars->matlab_data_out = NULL;
    }
    mclTerminateApplication();
    delete vars; vars = NULL;
    return true;
  }

  float* samples = trace->getTraceSamples();
  int numSamples = shdr->numSamples;

  memcpy( mxGetData(vars->matlab_data_in), samples, numSamples*sizeof(float) );
  mlfMatlabTest( 1, &vars->matlab_data_out, vars->matlab_data_in );
  memcpy( samples, mxGetData(vars->matlab_data_out), numSamples*sizeof(float) );
#endif

  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_matlab_test_( csParamDef* pdef ) {
  pdef->setModule( "MATLAB_TEST", "Matlab test module (EXPERIMENTAL MODULE)" );

  pdef->addParam( "method_name", "Matlab method name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Name of Matlab method to call.", "The shared object file must have the name lib<method_name>.so>" );

}

extern "C" void _params_mod_matlab_test_( csParamDef* pdef ) {
  params_mod_matlab_test_( pdef );
}
extern "C" void _init_mod_matlab_test_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_matlab_test_( param, env, log );
}
extern "C" bool _exec_mod_matlab_test_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_matlab_test_( trace, port, env, log );
}


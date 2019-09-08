/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include <cmath>
#include "csToken.h"
#include "csVector.h"
#include "csEquationSolver.h"
#include "csSort.h"
#include "cseis_curveFitting.h"
#include "csSelectionManager.h"
#include "geolib_methods.h"

using std::string;
using namespace cseis_geolib;
using namespace cseis_system;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: HDR_MATH_ENS
 *
 * @author Bjorn Olofsson
 * @date   2007
 */
namespace mod_hdr_math_ens {
  struct VariableStruct {
    int  nHeaders;
    int  hdrId1;
    int  hdrId2;
    int  hdrStddevId1;
    int  hdrStddevId2;
    char hdrType1;
    char hdrType2;
    int  method;
    int  nTraces;
    int  output_option;
    int  nEquations;
    bool computeStddev;
    cseis_geolib::csSort<double>* sortObj;
    csSelectionManager* selectionManager;
    int selectFailOption;
    int hdrId_selectFail;
    cseis_geolib::type_t hdrType_selectFail;
    double selectFailValue;
  };
  static const int MEDIAN  = 1;
  static const int MEAN    = 2;
  static const int XCOR_COS2 = 3;
  static const int LINEAR   = 4;
  static const int MAXIMUM  = 5;
  static const int MINIMUM  = 6;
  static const int MEAN_ANGLE = 7;

  static const int OUTPUT_ALL     = 0;
  static const int OUTPUT_FIRST   = 1;
  static const int OUTPUT_LAST    = 2;

  static const int SELECT_FAIL_NONE   = 0;
  static const int SELECT_FAIL_ALL    = 1;
  static const int SELECT_FAIL_HEADER = 2;
  static const int SELECT_FAIL_VALUE  = 3;
}

#define EPSILON 0.000001

void ls_linear_fit( double const* xValuesIn, double const* yValuesIn, int numValues, double* valuesOut );

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_hdr_math_ens_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
	//  csSuperHeader*    shdr = env->superHeader;
  mod_hdr_math_ens::VariableStruct* vars = new mod_hdr_math_ens::VariableStruct();
  edef->setVariables( vars );
  edef->setExecType( EXEC_TYPE_MULTITRACE );
  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE, 2 );

  vars->nHeaders = 1;
  vars->hdrId1   = -1;
  vars->hdrId2   = -1;
  vars->hdrStddevId1 = -1;
  vars->hdrStddevId2 = -1;
  vars->hdrType1 = 0;
  vars->hdrType2 = 0;
  vars->method   = mod_hdr_math_ens::MEAN;
  vars->nTraces  = 0;
  vars->output_option = mod_hdr_math_ens::OUTPUT_ALL;
  vars->nEquations = 0;
  vars->computeStddev = false;
  vars->sortObj     = NULL;
  vars->selectionManager = NULL;
  vars->selectFailOption = mod_hdr_math_ens::SELECT_FAIL_NONE;
  vars->hdrId_selectFail = -1;
  vars->hdrType_selectFail = TYPE_UNKNOWN;
  vars->selectFailValue = 0;

  //----------------------------------------------------------------
  std::string text;
  std::string headerName1("");
  std::string headerName2("");

  param->getString( "method", &text );
  text = toLowerCase( text );
  if( !text.compare("median") ) {
    vars->method = mod_hdr_math_ens::MEDIAN;
    param->getString( "header1", &headerName1 );
  }
  else if( !text.compare("mean") ) {
    vars->method = mod_hdr_math_ens::MEAN;
    param->getString( "header1", &headerName1 );
  }
  else if( !text.compare("minimum") ) {
    vars->method = mod_hdr_math_ens::MINIMUM;
    param->getString( "header1", &headerName1 );
  }
  else if( !text.compare("maximum") ) {
    vars->method = mod_hdr_math_ens::MAXIMUM;
    param->getString( "header1", &headerName1 );
  }
  else if( !text.compare("mean_angle") ) {
    vars->method = mod_hdr_math_ens::MEAN_ANGLE;
    param->getString( "header1", &headerName1 );
  }
  else if( !text.compare("linear") ) {
    vars->method = mod_hdr_math_ens::LINEAR;
    param->getString( "header1", &headerName1 );
    param->getString( "header2", &headerName2 );
    if( hdef->headerExists( headerName2.c_str() ) ) {
      vars->hdrId2   = hdef->headerIndex( headerName2.c_str() );
      vars->hdrType2 = hdef->headerType( vars->hdrId2 );
      vars->nHeaders = 2;
      if( vars->hdrType2 != TYPE_FLOAT && vars->hdrType2 != TYPE_DOUBLE ) {
        log->error("For method 'linear', trace header '%s' must be of floating point type", headerName2.c_str() );
      }
    }
    else {
      log->line(" Trace header '%s' does not exist.", headerName2.c_str() );
      env->addError();
    }
  }
  else if( !text.compare("xcor_cos2") ) {
    vars->method = mod_hdr_math_ens::XCOR_COS2;
    param->getString( "header1", &headerName1 );
    param->getString( "header2", &headerName2 );
    if( hdef->headerExists( headerName2.c_str() ) ) {
      vars->hdrId2   = hdef->headerIndex( headerName2.c_str() );
      vars->hdrType2 = hdef->headerType( vars->hdrId2 );
      vars->nHeaders = 2;
      if( vars->hdrType2 == TYPE_STRING ) {
        log->error("Header '%s' is of type string. This module requires trace headers to be of number type", headerName2.c_str() );
      }
    }
    else {
      log->line(" Trace header '%s' does not exist.", headerName2.c_str() );
      env->addError();
    }
  }
  else {
    log->line("Method not recognized: %s.", text.c_str());
    env->addError();
  }

  //----------------------------------------------------------------
  if( hdef->headerExists( headerName1.c_str() ) ) {
    vars->hdrId1   = hdef->headerIndex( headerName1.c_str() );
    vars->hdrType1 = hdef->headerType( vars->hdrId1 );
    if( vars->hdrType1 == TYPE_STRING ) {
      log->error("Header '%s' is of type string. This module requires trace headers to be of number type", headerName1.c_str() );
    }
  }
  else {
    log->line(" Trace header '%s' does not exist.", headerName1.c_str() );
    env->addError();
  }
  //----------------------------------------------------------------
  if( param->exists("output") ) {
    std::string outputText;
    param->getString( "output", &outputText );
    outputText = toLowerCase( outputText );
    if( !outputText.compare("all") ) {
      vars->output_option = mod_hdr_math_ens::OUTPUT_ALL;
    }
    else if( !outputText.compare("first") ) {
      vars->output_option = mod_hdr_math_ens::OUTPUT_FIRST;
    }
    else if( !outputText.compare("last") ) {
      vars->output_option = mod_hdr_math_ens::OUTPUT_LAST;
    }
    else {
      log->line("Output option not recognized: %s.", outputText.c_str());
      env->addError();
    }
  }
  //----------------------------------------------------------------
  std::string yesno;
  if( param->exists( "stddev" ) ) {
    param->getString( "stddev", &yesno );
    yesno = toLowerCase( yesno );
    if( !yesno.compare("yes") ) {
      vars->computeStddev = true;
      std::string headerName_stddev = std::string(headerName1 + "_stddev");
      if( !hdef->headerExists( headerName_stddev ) ) {
        hdef->addHeader( vars->hdrType1, headerName_stddev, "Standard deviation" );
      }
      else if( vars->hdrType1 != hdef->headerType( headerName_stddev ) ) {
        log->error("Trace header '%s' already exists but has wrong type. Should have same type as header '%s'.",
                   headerName_stddev.c_str(), headerName1.c_str() );
      }
      vars->hdrStddevId1 = hdef->headerIndex( headerName_stddev );
    }
    else if( !yesno.compare("no") ) {
      vars->computeStddev = false;
    }
    else {
      log->line("Option not recognized: %s. Should be 'yes' or 'no'", yesno.c_str());
      env->addError();
    }
  }

  if( vars->method == mod_hdr_math_ens::MEDIAN ) {
    vars->sortObj = new csSort<double>();
  }

  //---------------------------------------------------------
  // Create new headers
  if( param->exists("header_select") ) {
    if( vars->method == mod_hdr_math_ens::LINEAR ) {
      log->error("Sorry, header selection (header_select) is currently not supported for linear fitting method");
    }
    csVector<std::string> valueList;
    param->getAll( "header_select", &valueList );

    if( valueList.size() == 0 ) {
      log->line("Wrong number of arguments for user parameter 'header_select'. Expected: >0, found: %d.", valueList.size());
      env->addError();
    }
    else {
      std::string text;
      param->getString( "select", &text );

      vars->selectionManager = NULL;
      try {
        vars->selectionManager = new csSelectionManager();
        vars->selectionManager->set( &valueList, &text, hdef );
      }
      catch( csException& e ) {
        vars->selectionManager = NULL;
        log->error( "%s", e.getMessage() );
      }
      if( edef->isDebug() ) vars->selectionManager->dump();
    }
    if( param->exists("select_fail") ) {
      param->getString( "select_fail", &text );
      if( !text.compare("none") ) {
        vars->selectFailOption = mod_hdr_math_ens::SELECT_FAIL_NONE;
      }
      else if( !text.compare("all") ) {
        vars->selectFailOption = mod_hdr_math_ens::SELECT_FAIL_ALL;
      }
      else if( !text.compare("value") ) {
        vars->selectFailOption = mod_hdr_math_ens::SELECT_FAIL_VALUE;
        param->getString( "select_fail", &text, 1 );
        vars->selectFailValue = atof(text.c_str());
      }
      else if( !text.compare("header") ) {
        vars->selectFailOption = mod_hdr_math_ens::SELECT_FAIL_HEADER;
        param->getString( "select_fail", &text, 1 );
        if( !hdef->headerExists(text) ) {
          log->line("Trace header does not exist: %s.", text.c_str());
          env->addError();
        }
        else {
          vars->hdrId_selectFail   = hdef->headerIndex(text);
          vars->hdrType_selectFail = hdef->headerType(vars->hdrId_selectFail);
        }
      }
      else {
        log->line("Output option not recognized: %s.", text.c_str());
        env->addError();
      }
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
void exec_mod_hdr_math_ens_(
                            csTraceGather* traceGather,
                            int* port,
                            int* numTrcToKeep,
                            csExecPhaseEnv* env,
                            csLogWriter* log )
{
  mod_hdr_math_ens::VariableStruct* vars = reinterpret_cast<mod_hdr_math_ens::VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;

  if( edef->isCleanup() ) {
    if( vars->sortObj != NULL ) {
      delete vars->sortObj;
      vars->sortObj = NULL;
    }
    if( vars->selectionManager ) {
      delete vars->selectionManager; vars->selectionManager = NULL;
    }
    delete vars; vars = NULL;
    return;
  }
  //----------------------------

  int nTracesIn = traceGather->numTraces();

  //-----------------------------------------
  // Read in header values into value arrays
  //
  bool* isSelected = new bool[nTracesIn];
  int nTracesSelected = 0;
  if( vars->selectionManager == NULL ) {
    for( int itrc = 0; itrc < nTracesIn; itrc++) {
      isSelected[itrc] = true;
    }
    nTracesSelected = nTracesIn;
  }
  else {
    for( int itrc = 0; itrc < nTracesIn; itrc++) {
      if( vars->selectionManager->contains( traceGather->trace(itrc)->getTraceHeader() ) ) {
        isSelected[itrc] = true;
        nTracesSelected += 1;
      }
      else {
        isSelected[itrc] = false;
      }
    }
  }
  double* values  = NULL;
  double* values2 = NULL;

  if( nTracesSelected == 0 ) {
    log->warning("Module HDR_MATH_ENS: No traces match selection criteria");
    nTracesSelected = nTracesIn;
    for( int itrc = 0; itrc < nTracesIn; itrc++) {
      isSelected[itrc] = true;
    }
    if( vars->selectFailOption == mod_hdr_math_ens::SELECT_FAIL_ALL ) {
      // Nothing
    }
//    else if( vars->selectFailOption = SELECT_FAIL_NONE ) {
//    }
    else {  // if( vars->selectFailOption = SELECT_FAIL_HEADER ) {
      double value = vars->selectFailValue;
      for( int itrc = 0; itrc < nTracesIn; itrc++) {
        if( vars->selectFailOption == mod_hdr_math_ens::SELECT_FAIL_HEADER ) {
          value = traceGather->trace(itrc)->getTraceHeader()->doubleValue(vars->hdrId_selectFail);
        }
        if( vars->hdrType1 == TYPE_FLOAT ) {
          traceGather->trace(itrc)->getTraceHeader()->setFloatValue( vars->hdrId1, (float)value );
        }
        else if( vars->hdrType1 == TYPE_DOUBLE ) {
          traceGather->trace(itrc)->getTraceHeader()->setDoubleValue( vars->hdrId1, value );
        }
        else if( vars->hdrType1 == TYPE_INT ) {
          traceGather->trace(itrc)->getTraceHeader()->setIntValue( vars->hdrId1, (int)value );
        }
        else if( vars->hdrType1 == TYPE_INT64 ) {
          traceGather->trace(itrc)->getTraceHeader()->setInt64Value( vars->hdrId1, (csInt64_t)value );
        }
        if( vars->nHeaders == 2 ) {
          if( vars->hdrType2 == TYPE_FLOAT ) {
            traceGather->trace(itrc)->getTraceHeader()->setFloatValue( vars->hdrId2, (float)value );
          }
          else if( vars->hdrType2 == TYPE_DOUBLE ) {
            traceGather->trace(itrc)->getTraceHeader()->setDoubleValue( vars->hdrId2, value );
          }
          else if( vars->hdrType2 == TYPE_INT ) {
            traceGather->trace(itrc)->getTraceHeader()->setIntValue( vars->hdrId2, (int)value );
          }
          else if( vars->hdrType2 == TYPE_INT64 ) {
            traceGather->trace(itrc)->getTraceHeader()->setInt64Value( vars->hdrId2, (csInt64_t)value );
          }
        }
      }
    }
  }

//-------------------------------------------------------------------------------
    values = new double[nTracesSelected];
    int counter = 0;
    if( vars->hdrType1 == TYPE_FLOAT ) {
      for( int itrc = 0; itrc < nTracesIn; itrc++) {
        if( isSelected[itrc] ) values[counter++] = (double)traceGather->trace(itrc)->getTraceHeader()->floatValue(vars->hdrId1);
      }
    }
    else if( vars->hdrType1 == TYPE_DOUBLE ) {
      for( int itrc = 0; itrc < nTracesIn; itrc++) {
        if( isSelected[itrc] ) values[counter++] = traceGather->trace(itrc)->getTraceHeader()->doubleValue(vars->hdrId1);
      }
    }
    else if( vars->hdrType1 == TYPE_INT ) {
      for( int itrc = 0; itrc < nTracesIn; itrc++) {
        if( isSelected[itrc] ) values[counter++] = (double)traceGather->trace(itrc)->getTraceHeader()->intValue(vars->hdrId1);
      }
    }
    else if( vars->hdrType1 == TYPE_INT64 ) {
      for( int itrc = 0; itrc < nTracesIn; itrc++) {
        if( isSelected[itrc] ) values[counter++] = (double)traceGather->trace(itrc)->getTraceHeader()->int64Value(vars->hdrId1);
      }
    }
    else {
      log->error("mod_hdr_math_ens: Unsupported trace header type, code: %d", vars->hdrType1 );
    }
    if( vars->nHeaders == 2 ) {
      int counter = 0;
      values2 = new double[nTracesSelected];
      if( vars->hdrType2 == TYPE_FLOAT ) {
        for( int itrc = 0; itrc < nTracesIn; itrc++) {
          if( isSelected[itrc] ) values2[counter++] = (double)traceGather->trace(itrc)->getTraceHeader()->floatValue(vars->hdrId2);
        }
      }
      else if( vars->hdrType2 == TYPE_DOUBLE ) {
        for( int itrc = 0; itrc < nTracesIn; itrc++) {
          if( isSelected[itrc] ) values2[counter++] = traceGather->trace(itrc)->getTraceHeader()->doubleValue(vars->hdrId2);
        }
      }
      else if( vars->hdrType2 == TYPE_INT ) {
        for( int itrc = 0; itrc < nTracesIn; itrc++) {
          if( isSelected[itrc] ) values2[counter++] = (double)traceGather->trace(itrc)->getTraceHeader()->intValue(vars->hdrId2);
        }
      }
      else if( vars->hdrType2 == TYPE_INT64 ) {
        for( int itrc = 0; itrc < nTracesIn; itrc++) {
          if( isSelected[itrc] ) values2[counter++] = (double)traceGather->trace(itrc)->getTraceHeader()->int64Value(vars->hdrId2);
        }
      }
    }

  //----------------------------------------------------
  // Compute function
  //
  double result = 0.0;
  double result2 = 0.0;
  double stddev = 0.0;
  if( vars->method == mod_hdr_math_ens::MEDIAN ) {
    vars->sortObj->treeSort( values, nTracesSelected );
    result = values[nTracesSelected/2];
  }
  else if( vars->method == mod_hdr_math_ens::MEAN ) {
    for( int itrc = 0; itrc < nTracesSelected; itrc++) {
      result += values[itrc];
    }
    result /= nTracesSelected;
  }
  else if( vars->method == mod_hdr_math_ens::MEAN_ANGLE ) {
    // Assumes input values are angles in degrees
    // --> Take average of angular vector, then compute back average angle
    for( int itrc = 0; itrc < nTracesSelected; itrc++) {
      double an_rad = values[itrc] * M_PI / 180.0;
      result  += sin(an_rad);
      result2 += cos(an_rad);
    }
    result  /= nTracesSelected;
    result2 /= nTracesSelected;
    // Convert vector back into angle:
    result  = fmod( atan2(result,result2) * 180.0 / M_PI + 360.0, 360.0 );
  }
  else if( vars->method == mod_hdr_math_ens::MAXIMUM ) {
    result = values[0];
    for( int itrc = 1; itrc < nTracesSelected; itrc++) {
      if( values[itrc] > result ) result = values[itrc];
    }
  }
  else if( vars->method == mod_hdr_math_ens::MINIMUM ) {
    result = values[0];
    for( int itrc = 1; itrc < nTracesSelected; itrc++) {
      if( values[itrc] < result ) result = values[itrc];
    }
  }
  else if( vars->method == mod_hdr_math_ens::LINEAR ) {
    double coefficients[2];
    polynom_fit( values, values2, nTracesSelected, 1, coefficients );

    if( edef->isDebug() ) log->line("HDR_MATH_ENS: Linear polynomial coefficients: %e x^0  %e x^1\n", coefficients[0], coefficients[1]);
    //    ls_linear_fit( values, values2, nTraces, values2 );
    for( int itrc = 0; itrc < nTracesSelected; itrc++) {
      values2[itrc] = values[itrc] * coefficients[1] + coefficients[0];
    }

    if( vars->hdrType2 == TYPE_FLOAT ) {
      for( int itrc = 0; itrc < nTracesSelected; itrc++) {
        traceGather->trace(itrc)->getTraceHeader()->setFloatValue( vars->hdrId2, (float)values2[itrc] );
      }
    }
    else if( vars->hdrType2 == TYPE_DOUBLE ) {
      for( int itrc = 0; itrc < nTracesSelected; itrc++) {
        traceGather->trace(itrc)->getTraceHeader()->setDoubleValue( vars->hdrId2, values2[itrc] );
      }
    }
    if( values != NULL ) delete [] values;
    if( values2 != NULL ) delete [] values2;
    if( isSelected != NULL )  delete [] isSelected;

    return;
  }
  //--------------------------------------------------------------------------------
  else if( vars->method == mod_hdr_math_ens::XCOR_COS2 ) {
    // For shear-wave splitting analysis, array 'values' contains azimuth angle in radians, 'values' contains cross-correlation time lag

    int periodicity = 2;
    for( int itrc = 0; itrc < nTracesSelected; itrc++) {
      values[itrc] = DEG2RAD( fmod(values[itrc],360.0) );
    }
    computeXcorCos( values, values2, nTracesSelected, periodicity, vars->computeStddev, result, stddev, result2 );

/*
    result  = values[0];
    result2 = values2[0];
    int maxtrc = 0;
    for( int itrc = 1; itrc < nTracesSelected; itrc++) {
//      fprintf(stdout,"%d   %f %f %f\n", itrc, values[itrc], RAD2DEG(values[itrc]), values2[itrc]);
      if( values2[itrc] > result2 ) {
        result2 = values2[itrc];
        result = values[itrc];
        maxtrc = itrc;
      //  fprintf(stdout,"  %d %f %f %f\n", itrc, result, RAD2DEG(result), result2);
      }
    }
    //fprintf(stdout,"\n");
    
    
    float val2[3];
    if( maxtrc == 0 ) {
      values[2] = values[1];
      values[1] = values[0];
      values[0] = values[nTracesSelected-1];
      val2[2] = values2[1];
      val2[1] = values2[0];
      val2[0] = values2[nTracesSelected-1];
      maxtrc = 1;
    }
    else if( maxtrc == nTracesSelected-1 ) {
      values[2] = values[0];
      values[1] = values[nTracesSelected-1];
      values[0] = values[nTracesSelected-2];
      val2[2] = values2[0];
      val2[1] = values2[nTracesSelected-1];
      val2[0] = values2[nTracesSelected-2];
      maxtrc = 1;
    }
    else {
      values[2] = values[maxtrc+1];
      values[1] = values[maxtrc];
      values[0] = values[maxtrc-1];
      val2[2] = values2[maxtrc+1];
      val2[1] = values2[maxtrc];
      val2[0] = values2[maxtrc-1];
      maxtrc = 1;
    }
    float result3;
    float sampleIndex_maxAmp_float = getQuadMaxSample( val2, maxtrc, nTracesSelected, &result3 );
    result2 = (double)result3;
    int sampleIndex_int = (int)sampleIndex_maxAmp_float;
    if( sampleIndex_int == 2 ) sampleIndex_int = 1;
    result = (sampleIndex_maxAmp_float-sampleIndex_int) * ( values[sampleIndex_int+1]-values[sampleIndex_int] ) + values[sampleIndex_int];

     stddev = 0;
*/
  }
  // Standard deviation
  if( vars->computeStddev && vars->method != mod_hdr_math_ens::XCOR_COS2 ) {
    for( int itrc = 0; itrc < nTracesSelected; itrc++) {
      stddev += CS_SQR( values[itrc] - result );
    }
    stddev = sqrt( stddev/(double)nTracesSelected );
  }

  if( values != NULL ) delete [] values;
  if( values2 != NULL ) delete [] values2;
  if( isSelected != NULL ) delete [] isSelected;
  //--------------------------------------------------------
  // Free traces that shall not be output from this module
  //
  if( nTracesIn > 1 && vars->output_option != mod_hdr_math_ens::OUTPUT_ALL ) {
    if( vars->output_option == mod_hdr_math_ens::OUTPUT_LAST ) {
      traceGather->freeTraces( 0, nTracesIn-1 );
    }
    else if( vars->output_option == mod_hdr_math_ens::OUTPUT_FIRST ) {
      traceGather->freeTraces( 1, nTracesIn-1 );
    }
  }

  int nTracesOut = traceGather->numTraces();
    
  //--------------------------------------------------------
  // Put results back into trace headers
  //
  if( vars->hdrType1 == TYPE_FLOAT ) {
    for( int itrc = 0; itrc < nTracesOut; itrc++) {
      traceGather->trace(itrc)->getTraceHeader()->setFloatValue( vars->hdrId1, (float)result );
    }
    if( vars->computeStddev ) {
      for( int itrc = 0; itrc < nTracesOut; itrc++) {
        traceGather->trace(itrc)->getTraceHeader()->setFloatValue( vars->hdrStddevId1, (float)stddev );
      }
    }
  }
  else if( vars->hdrType1 == TYPE_DOUBLE ) {
    for( int itrc = 0; itrc < nTracesOut; itrc++) {
      traceGather->trace(itrc)->getTraceHeader()->setDoubleValue( vars->hdrId1, result );
    }
    if( vars->computeStddev ) {
      for( int itrc = 0; itrc < nTracesOut; itrc++) {
        traceGather->trace(itrc)->getTraceHeader()->setDoubleValue( vars->hdrStddevId1, stddev );
      }
    }
  }
  else if( vars->hdrType1 == TYPE_INT ) {
    for( int itrc = 0; itrc < nTracesOut; itrc++) {
      traceGather->trace(itrc)->getTraceHeader()->setIntValue( vars->hdrId1, (int)result );
    }
    if( vars->computeStddev ) {
      for( int itrc = 0; itrc < nTracesOut; itrc++) {
        traceGather->trace(itrc)->getTraceHeader()->setIntValue( vars->hdrStddevId1, (int)stddev );
      }
    }
  }
  else if( vars->hdrType1 == TYPE_INT64 ) {
    for( int itrc = 0; itrc < nTracesOut; itrc++) {
      traceGather->trace(itrc)->getTraceHeader()->setInt64Value( vars->hdrId1, (csInt64_t)result );
    }
    if( vars->computeStddev ) {
      for( int itrc = 0; itrc < nTracesOut; itrc++) {
        traceGather->trace(itrc)->getTraceHeader()->setInt64Value( vars->hdrStddevId1, (csInt64_t)stddev );
      }
    }
  }

  if( vars->method == mod_hdr_math_ens::XCOR_COS2 ) {
    for( int itrc = 0; itrc < nTracesOut; itrc++) {
      traceGather->trace(itrc)->getTraceHeader()->setDoubleValue( vars->hdrId2, result2 );
    }
  }
  
}

//********************************************************************************
//
//
void params_mod_hdr_math_ens_( csParamDef* pdef ) {
  pdef->setModule( "HDR_MATH_ENS", "Multi-trace header computation", "Apply multi-trace mathematical equation" );

  pdef->addParam( "method", "", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_OPTION );
  pdef->addOption( "median", "Compute median value of header1" );
  pdef->addOption( "mean", "Compute mean value of header1" );
  pdef->addOption( "linear", "Fit linear polynom to header1/header2 value pairs" );
  pdef->addOption( "xcor_cos2", "Cross-correlate with cosine function (2*theta variation)", "header2: S1 azimuth [deg], header2: Time lag from cross-correlation" );
  pdef->addOption( "minimum", "Take minimum value of header1" );
  pdef->addOption( "maximum", "Take maximum value of header1" );
  pdef->addOption( "mean_angle", "Compute mean value of header1 (contains angle in degrees)" );

  pdef->addParam( "header1", "Trace header 1 to be used in computation", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "header2", "Trace header 2 to be used in computation", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "output", "", NUM_VALUES_FIXED );
  pdef->addValue( "all", VALTYPE_OPTION );
  pdef->addOption( "all", "Output all traces. Each trace contains the result in the specified header" );
  pdef->addOption( "last", "Output last trace only" );
  pdef->addOption( "first", "Output first trace only" );

  pdef->addParam( "stddev", "", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Compute standard deviation.", "A new trace header with extension _stddev will be created" );
  pdef->addOption( "no", "Do not compute standard deviation" );

  pdef->addParam( "header_select", "Names of trace headers used for trace selection", NUM_VALUES_VARIABLE,
                  "Only selected traces will form input to specified analysis method. All output traces will be updated with the analysis result" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "select", "Selection of header values", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING,
                  "List of selection strings, one for each specified header. See documentation for more detailed description of selection syntax");

  pdef->addParam( "select_fail", "What to do when selection fails?", NUM_VALUES_VARIABLE, "This parameter gives several options what will be done when the given selection criteria in parameters 'header_select' and 'select' does not match any input trace" );
  pdef->addValue( "all", VALTYPE_OPTION );
  pdef->addOption( "all", "Perform specified header operation on all traces in input gather" );
  pdef->addOption( "none", "Do not perform operation on any input traces. Leave specified trace header unchanged" );
  pdef->addOption( "value", "Set specified trace header to constant value" );
  pdef->addOption( "header", "Set specified trace header to value from another trace header" );
  pdef->addValue( "", VALTYPE_STRING, "For option 'value': Constant value; for option 'header': Name of trace header; otherwise blank");
}

extern "C" void _params_mod_hdr_math_ens_( csParamDef* pdef ) {
  params_mod_hdr_math_ens_( pdef );
}
extern "C" void _init_mod_hdr_math_ens_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_hdr_math_ens_( param, env, log );
}
extern "C" void _exec_mod_hdr_math_ens_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* log ) {
  exec_mod_hdr_math_ens_( traceGather, port, numTrcToKeep, env, log );
}

void ls_linear_fit( double const* xValuesIn, double const* yValuesIn, int numValues, double* valuesOut ) {

  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_xx = 0.0;
  double sum_xy = 0.0;
  for( int ival = 0; ival < numValues; ival++ ) {
    double xval = xValuesIn[ival];
    double yval = yValuesIn[ival];

    sum_x += xval;
    sum_y += yval;

    sum_xx += xval*xval;
    sum_xy += xval*yval;
  }

  double tmp1 = sum_x*sum_x - numValues*sum_xx;

  if( fabs(tmp1) != 0.0 ) {
    double aa = ((sum_x * sum_xy) - (sum_xx * sum_y)) / tmp1;
    double bb = ((sum_x * sum_y) - (numValues * sum_xy)) / tmp1;
    for( int ival = 0; ival < numValues; ival++ ) {
      valuesOut[ival] = aa * xValuesIn[ival] + bb;
    }
  }

}



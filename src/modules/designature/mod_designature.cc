/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csFlexNumber.h"
#include "csFFTTools.h"
#include "csFFTDesignature.h"
#include "csASCIIFileReader.h"
#include "csFileUtils.h"
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: DESIGNATURE
 *
 * @author Bjorn Olofsson
 * @date   2011
 */
namespace mod_designature {
  struct VariableStruct {
    csFFTDesignature* fftDesig;
    int asciiFormat;
    std::string filenameOut;
    bool isFileOutWavelet;
  };
  static int const UNIT_MS = 3;
  static int const UNIT_S  = 4;
}
using mod_designature::VariableStruct;

//*************************************************************************************************
// Init phase
//
//
//
//*************************************************************************************************
void init_mod_designature_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->asciiFormat     = -1;
  vars->fftDesig        = NULL;
  vars->filenameOut     = "";
  vars->isFileOutWavelet = false;

//---------------------------------------------
//

  std::string filename;
  std::string text;
  int unitTime_ascii = mod_designature::UNIT_MS;

  bool overrideSampleInt = false;
  if( param->exists("override_sample_int") ) {
    param->getString("override_sample_int", &text);
    if( !text.compare("no") ) {
      overrideSampleInt = false;
    }
    else if( !text.compare("yes") ) {
      overrideSampleInt = true;
    }
    else {
      log->error("Unknown option: '%s'", text.c_str());
    }
  }
  param->getString("input_wavelet", &filename);
  if( param->getNumValues("input_wavelet") > 1 ) {
    param->getString("input_wavelet", &text, 1);
    if( !text.compare("ms") ) {
      unitTime_ascii = mod_designature::UNIT_MS;
    }
    else if( !text.compare("s") ) {
      unitTime_ascii = mod_designature::UNIT_S;
    }
    else {
      log->error("Unknown option: '%s'", text.c_str());
    }
  }
  param->getString("format",&text,0);
  if( !text.compare("columns") ) {
    vars->asciiFormat = cseis_io::csASCIIFileReader::FORMAT_COLUMNS;
  }
  else if( !text.compare("signature") ) {
    vars->asciiFormat = cseis_io::csASCIIFileReader::FORMAT_NUCLEUS_SIGNATURE;
  }
  else {
    log->error("Unknown option: %s", text.c_str() );
  }

  int colIndexTime  = 0;
  int colIndexValue = 1;
  int colIndexTrace = -1;
  if( vars->asciiFormat == cseis_io::csASCIIFileReader::FORMAT_COLUMNS ) {
    if( param->exists("format") && param->getNumValues("format") > 1 ) {
      param->getInt( "format", &colIndexTime, 1 );
      colIndexTime  -= 1;
      if( param->getNumValues("format") > 2 ) {
        param->getInt( "format", &colIndexValue, 2 );
        colIndexValue -= 1;
        if( param->getNumValues("format") > 3 ) {
          param->getInt( "format", &colIndexTrace, 3 );
          colIndexTrace -= 1;
        }
      }
    }
  }

  float percWhiteNoise = 0.01f;
  if( param->exists("white_noise") ) {
    param->getFloat("white_noise",&percWhiteNoise);
  }

  //----------------------------------------------------------------------
  // Read in signature from external ASCII file
  //
  cseis_io::ASCIIParam asciiParam;
  try {
    cseis_io::csASCIIFileReader asciiFileReader( filename, vars->asciiFormat );
    bool success = asciiFileReader.initialize( &asciiParam, colIndexTime, colIndexValue, colIndexTrace );
    if( !success ) log->error("Unknown error occurred during initialization of signature input file. Incorrect or unsupported format?");
    success = asciiFileReader.readNextTrace( &asciiParam );
    if( !success ) log->error("Unknown error occurred when reading in samples from signature input file. Incorrect or unsupported format?");
  }
  catch( csException& e ) {
    log->error("Error occurred when initializing input ASCII file: %s", e.getMessage() );
  }

  //
  //----------------------------------------------------------------------

  float sampleInt_ms = asciiParam.sampleInt;
  float timeZero_ms  = -(float)asciiParam.timeFirstSamp;
  // ^^ Minus sign in case time of first sample is negative...or so. See Nucleus format
  // ..this is required because the csDesignature class simply assumes time of first sample = 0

  if( unitTime_ascii == mod_designature::UNIT_S ) {
    sampleInt_ms *= 1000.0;
    timeZero_ms  *= 1000.0;
  }
  if( sampleInt_ms != shdr->sampleInt ) {
    if( !overrideSampleInt ) {
      log->error("ASCII input file has different sample interval (=%f ms) than input data (=%f ms). Unsupported case.", sampleInt_ms, shdr->sampleInt);
    }
    else {
      log->warning("ASCII input file has different sample interval (=%f ms) than input data (=%f ms). Ignored.", sampleInt_ms, shdr->sampleInt);
    }
  }
  if( asciiParam.numSamples() > shdr->numSamples ) {
    log->error("ASCII input file has more samples (=%d) than input data (=%d). Unsupported case.", asciiParam.numSamples(), shdr->numSamples);
  }
  if( asciiParam.numSamples() <= 0 ) {
    log->error("Could not read in any sample values from input file. Unsupported or incorrect file format?");
  }

  float* buffer = new float[shdr->numSamples];
  for( int isamp = 0; isamp < asciiParam.numSamples(); isamp++ ) {
    buffer[isamp] = asciiParam.sample(isamp);
  }
  // Pad signature file with zeros...
  for( int isamp = asciiParam.numSamples(); isamp < shdr->numSamples; isamp++ ) {
    buffer[isamp] = 0.0;
  }
  

  //----------------------------------------------------------------------
  // Read in output wavelet from external ASCII file
  //
  std::string filenameOut;
  int unitTimeOut_ascii = unitTime_ascii;
  float* bufferOut = NULL;
  if( param->exists("output_wavelet") ) {
    param->getString("output_wavelet", &filenameOut);
    if( param->getNumValues("output_wavelet") > 1 ) {
      param->getString("output_wavelet", &text, 1);
      if( !text.compare("ms") ) {
        unitTimeOut_ascii = mod_designature::UNIT_MS;
      }
      else if( !text.compare("s") ) {
        unitTimeOut_ascii = mod_designature::UNIT_S;
      }
      else {
        log->error("Unknown option: '%s'", text.c_str());
      }
    }

    cseis_io::ASCIIParam asciiParamOut;
    try {
      cseis_io::csASCIIFileReader asciiFileReader( filenameOut, vars->asciiFormat );
      bool success = asciiFileReader.initialize( &asciiParamOut, colIndexTime, colIndexValue, colIndexTrace );
      if( !success ) log->error("Unknown error occurred during initialization of output wavelet file. Incorrect or unsupported format?");
      success = asciiFileReader.readNextTrace( &asciiParamOut );
      if( !success ) log->error("Unknown error occurred when reading in samples from output wavelet file. Incorrect or unsupported format?");
    }
    catch( csException& e ) {
      log->error("Error occurred when initializing input ASCII file: %s", e.getMessage() );
    }
    if( asciiParamOut.numSamples() != asciiParam.numSamples() ) {
      log->error("Output wavelet contains %d number of samples, not matching the number of samples in input wavelet = %d.\nThis is not supported",
                 asciiParamOut.numSamples() != asciiParam.numSamples() );
    }
    bufferOut = new float[shdr->numSamples];
    for( int isamp = 0; isamp < asciiParamOut.numSamples(); isamp++ ) {
      bufferOut[isamp] = asciiParamOut.sample(isamp);
    }
    // Pad signature file with zeros...
    for( int isamp = asciiParamOut.numSamples(); isamp < shdr->numSamples; isamp++ ) {
      bufferOut[isamp] = 0.0;
    }
  }
  if( unitTimeOut_ascii != unitTime_ascii ) {
    log->error("Different time units specified for input & output wavelet. This is not supported.");
  }
  //--------------------------------------------------

  if( param->exists("zero_time") ) {
    param->getFloat("zero_time",&timeZero_ms);
    if( param->getNumValues("zero_time") > 1 ) {
      string unit;
      param->getString( "zero_time", &unit, 1 );
      if( !unit.compare("samples") ) {
        timeZero_ms *= asciiParam.sampleInt;
      }
      else if( !unit.compare("s") ) {
        timeZero_ms *= 1000.0;
      }
      else if( !unit.compare("ms") ) {
        // Nothing
      }
    }
  }
  if( edef->isDebug() ) {
    log->line("Zero time = %f ms", timeZero_ms);
  }

  vars->fftDesig = new csFFTDesignature( shdr->numSamples, shdr->sampleInt, buffer, timeZero_ms/1000.0, percWhiteNoise, bufferOut );
  delete [] buffer;
  if( bufferOut != NULL ) delete [] bufferOut;

  if( param->exists("option") ) {
    int filterType = csFFTDesignature::AMP_PHASE;
    param->getString("option", &text );
    if( !text.compare("phase_only") ) {
      filterType = csFFTDesignature::PHASE_ONLY;
    }
    else if( !text.compare("amp_only") ) {
      filterType = csFFTDesignature::AMP_ONLY;
    }
    else if( !text.compare("amp_phase") ) {
      filterType = csFFTDesignature::AMP_PHASE;
    }
    else {
      log->error("Unknown option: %s", text.c_str() );
    }
    vars->fftDesig->setDesigFilterType( filterType );
  }

  float freqNy = 500.0f/asciiParam.sampleInt;
  float orderDefault = 5;
  if( param->exists("highpass") ) {
    float cutOff;
    float order = orderDefault;
    param->getFloat("highpass", &cutOff );
    if( cutOff < 0 || cutOff > freqNy ) {
      log->error("Filter frequency exceeds valid range (0-%fHz): %fHz", freqNy, cutOff );
    }
    if( param->getNumValues("highpass") > 1 ) {
      float slope;
      param->getFloat("highpass", &slope, 1 );
      order = fabs(slope) / 6.0;
    }
    vars->fftDesig->setDesigHighPass( cutOff, order );
  }
  if( param->exists("lowpass") ) {
    float cutOff;
    float order = orderDefault;
    param->getFloat("lowpass", &cutOff );
    if( param->getNumValues("lowpass") > 1 ) {
      float slope;
      param->getFloat("lowpass", &slope, 1 );
      order = fabs(slope) / 6.0;
    }
    if( cutOff < 0 || cutOff > freqNy ) {
      log->error("Filter frequency exceeds valid range (0-%fHz): %fHz", freqNy, cutOff );
    }
    vars->fftDesig->setDesigLowPass( cutOff, order );
  }

  if( param->exists("high_set") ) {
    float freqHighSet = 0;
    param->getFloat("high_set", &freqHighSet);
    vars->fftDesig->setDesigHighEnd( freqHighSet );
  }


  //--------------------------------------------------
  if( param->exists("notch_suppression") ) {
    float notchFreq;
    float notchWidth;
    int numLines = param->getNumLines("notch_suppression");
    for( int iline = 0; iline < numLines; iline++ ) {
      param->getFloatAtLine("notch_suppression", &notchFreq, iline, 0);
      param->getFloatAtLine("notch_suppression", &notchWidth, iline, 1);
      vars->fftDesig->setNotchSuppression( notchFreq, notchWidth );
    }
  }


  // Write filter to output file
  if( param->exists("filename_output") ) {
    param->getString("filename_output", &vars->filenameOut);
    if( !csFileUtils::createDoNotOverwrite( vars->filenameOut ) ) {
      log->error("Unable to open filter output file %s. Wrong path name..?", vars->filenameOut.c_str() );
    }
    if( param->getNumValues("filename_output") > 1 ) {
      param->getString("filename_output", &text, 1);
      if( !text.compare("wavelet") ) {
	vars->isFileOutWavelet = true;
      }
      else if( !text.compare("spectrum") ) {
	vars->isFileOutWavelet = false;
      }
      else {
	log->error("Unknown option: '%s'", text.c_str());
      }
    }
  }


  if( edef->isDebug() ) {
    //    log->line("---------------ASCII file dump-------------------");
    for( int isamp = 0; isamp < asciiParam.numSamples(); isamp++ ) {
      //      log->line("%d %f", isamp, asciiParam.sample(isamp) );
      fprintf( stdout, "%d %f\n", isamp, asciiParam.sample(isamp) );
    }
  }

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_designature_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;

  if( edef->isCleanup()){
    if( vars->fftDesig != NULL ) {
      delete vars->fftDesig;
      vars->fftDesig = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  if( vars->filenameOut.size() > 0 ) {
    FILE* fout = fopen(vars->filenameOut.c_str(), "w");
    if( vars->isFileOutWavelet ) {
      vars->fftDesig->dump_wavelet(fout,true);
    }
    else {
      vars->fftDesig->dump_spectrum(fout);
    }
    fclose(fout);
    fout = NULL;
    vars->filenameOut = "";
  }
  float* samples = trace->getTraceSamples();
  vars->fftDesig->applyFilter( samples, shdr->numSamples );

  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_designature_( csParamDef* pdef ) {
  pdef->setModule( "DESIGNATURE", "Designature filter operation" );
  pdef->addDoc("This module generates a filter (in the frequency domain) from a specified signature wavelet, and applies this filter to all input traces.");
  pdef->addDoc("Currently supported options: Create zero-phasing filter that removes the specified signature; in other words, apply the transfer function that converts the specified signature wavelet into a zero-phase spike.");
  pdef->addDoc("Add specified white noise to signature spectrum before computing transfer function. Finally, the designature filter spectrum can be further bandpass filtered.");

  pdef->addDoc("Note that this is a naive implementation, using spectral division. The user has the choice to limit the frequency range, add white noise, and to fill notches to alleviate artefacts.");


  pdef->addParam( "input_wavelet", "Name of file containing input signature/response wavelet", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit of time values found in input file, if any" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );

  pdef->addParam( "format", "Input wavelet ASCII file format", NUM_VALUES_VARIABLE);
  pdef->addValue( "columns", VALTYPE_OPTION );
  pdef->addOption( "signature", "Read in source signature from Nucleus ASCII file" );
  pdef->addOption( "columns", "Simple file format with 2 or 3 columns:  Time[ms]  Amplitude  (Trace number)",
    "Trace number column is optional." );
  pdef->addValue( "1", VALTYPE_NUMBER, "For 'column' format: Column number specifying time (or other sample unit)" );
  pdef->addValue( "2", VALTYPE_NUMBER, "For 'column' format: Column number specifying sample value" );
  pdef->addValue( "-1", VALTYPE_NUMBER, "Optional: Column number specifying trace number. -1: Trace number not applicable" );

  pdef->addParam( "output_wavelet", "Name of file containing output wavelet", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit of time values found in input file, if any" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );

  pdef->addParam( "zero_time", "Zero time in input wavelet", NUM_VALUES_VARIABLE,
                  "Specifying this parameter overrides the zero time that is found in the ASCII signature file" );
  pdef->addValue( "0", VALTYPE_NUMBER );
  pdef->addValue( "ms", VALTYPE_OPTION, "Unit" );
  pdef->addOption( "ms", "Milliseconds" );
  pdef->addOption( "s", "Seconds" );
  pdef->addOption( "samples", "Samples" );

  pdef->addParam( "white_noise", "Amount of white noise (in percent of maximum amplitude) to add to the signature/response spectrum",
                  NUM_VALUES_FIXED );
  pdef->addValue( "0.01", VALTYPE_NUMBER );

  pdef->addParam( "notch_suppression", "Suppress notches in signature wavelet", NUM_VALUES_FIXED,
                  "Apply cosine taper around notch frequency to filter" );
  pdef->addValue( "", VALTYPE_NUMBER, "Notch frequency [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "Width of suppression filter [Hz]" );

  pdef->addParam( "option", "Filter type option", NUM_VALUES_FIXED);
  pdef->addValue( "amp_phase", VALTYPE_OPTION );
  pdef->addOption( "amp_phase", "Create amplitude & phase designature filter" );
  pdef->addOption( "amp_only", "Create amplitude only designature filter" );
  pdef->addOption( "phase_only", "Create phase only designature filter" );

  pdef->addParam( "lowpass", "Lowpass filter to apply to designature filter before application", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Cutoff frequency for low-pass filter [Hz]",
     "The cutoff frequency will be damped by -3db" );
  pdef->addValue( "30", VALTYPE_NUMBER, "Filter slope (db/oct)" );

  pdef->addParam( "highpass", "Highpass filter to apply to designature filter before application", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Cutoff frequency for highpass filter [Hz]",
     "The cutoff frequency will be damped by -3db" );
  pdef->addValue( "30", VALTYPE_NUMBER, "Filter slope (db/oct)" );

  pdef->addParam( "high_set", "Set high-end of filter constant", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Frequency [Hz]", "Above this frequency, keep amplitude and phase constant");

  pdef->addParam( "filename_output", "Name of file where designature filter shall be written to", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addValue( "spectrum", VALTYPE_OPTION );
  pdef->addOption( "spectrum", "Output filter as amplitude/phase spectrum" );
  pdef->addOption( "wavelet", "Output filter as wavelet" );

  pdef->addParam( "override_sample_int", "Override sample interval?", NUM_VALUES_FIXED);
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "" );
  pdef->addOption( "yes", "Ignore sample interval of input wavelet. Assume it is the same as the input data" );
}

extern "C" void _params_mod_designature_( csParamDef* pdef ) {
  params_mod_designature_( pdef );
}
extern "C" void _init_mod_designature_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_designature_( param, env, log );
}
extern "C" bool _exec_mod_designature_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_designature_( trace, port, env, log );
}

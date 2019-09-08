/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csFFTDesignature.h"
#include "csException.h"
#include "geolib_math.h"
#include "geolib_defines.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace cseis_geolib;

/*
 * Designature filter explained:
 * numSamples (input wavelet)
 */

csFFTDesignature::csFFTDesignature( int numSamples, float sampleInt, float const* input_wavelet, float timeZero_s, float percWhiteNoise, float const* output_wavelet ) :
  csFFTTools( numSamples, sampleInt )
{
  //  myNotchIndexFirstRed = 0;
  //  myNotchIndexLastRed  = 0;
  //  myIsNotchSuppression = false;
  myDesigAmpFilter  = NULL;
  myDesigPhaseShift = NULL;
  myAmpSpecIn   = NULL;
  myPhaseSpecIn = NULL;
  myFilterType = AMP_PHASE;
  initDesignature( input_wavelet, timeZero_s, percWhiteNoise, output_wavelet );
}
csFFTDesignature::~csFFTDesignature() {
  if( myDesigAmpFilter != NULL ) {
    delete [] myDesigAmpFilter;
    myDesigAmpFilter  = NULL;
  }
  if( myDesigPhaseShift != NULL ) {
    delete [] myDesigPhaseShift;
    myDesigPhaseShift  = NULL;
  }
  if( myAmpSpecIn != NULL ) {
    delete [] myAmpSpecIn;
    myAmpSpecIn = NULL;
  }
  if( myPhaseSpecIn != NULL ) {
    delete [] myPhaseSpecIn;
    myPhaseSpecIn = NULL;
  }
}

//--------------------------------------------------------------------------------
//
//
//
void csFFTDesignature::initDesignature( float const* input_wavelet, float timeZero_s, float percWhiteNoise,
                                        float const* output_wavelet ) {
  myAmpSpecIn       = new float[myNumFFTSamplesIn/2+1];
  myPhaseSpecIn     = new float[myNumFFTSamplesIn/2+1];
  myDesigAmpFilter  = new float[myNumFFTSamplesIn/2+1];
  myDesigPhaseShift = new float[myNumFFTSamplesIn/2+1];

  float* ampSpecOut   = NULL;
  float* phaseSpecOut = NULL;

  // Forward transform wavelet
  fft_forward( input_wavelet, myAmpSpecIn, myPhaseSpecIn );
  float df = 1000.0 / ( (float)myNumFFTSamplesIn*mySampleIntIn );

  // Compute maximum amplitude, for application of percentage white noise
  float maxAmp = myAmpSpecIn[0];
  for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
    if( myAmpSpecIn[i] > maxAmp ) maxAmp = myAmpSpecIn[i];
    //  float df = 1.0 / sampleIntOut; // [Hz]
    //  float freq = df * (float)(i);
  }

  float whiteNoise = (percWhiteNoise/100.0) * maxAmp;
  if( whiteNoise <= 0.0 ) whiteNoise = 1e-50;

  // Spiking filter:
  if( output_wavelet == NULL ) {
    // Compute inverse (designature) filter = Inverse amplitude spectrum after adding white noise
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      myDesigAmpFilter[i] = maxAmp / (myAmpSpecIn[i] + whiteNoise);
    }

    // Compute phase shift for zero-phasing
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      float freq          = df * (float)(i);
      float phaseShift    = 2.0 * M_PI * freq * timeZero_s;
      myDesigPhaseShift[i] = fmod( phaseShift - myPhaseSpecIn[i] , 2.0*M_PI );
      if( myDesigPhaseShift[i] > M_PI ) {
        myDesigPhaseShift[i] -= 2.0*M_PI;
      }
    }
  }
  // 'Transfer' filter from input to output wavelet
  else {
    ampSpecOut       = new float[myNumFFTSamplesIn/2+1];
    phaseSpecOut     = new float[myNumFFTSamplesIn/2+1];
    fft_forward( output_wavelet, ampSpecOut, phaseSpecOut );

    float maxAmpOut = ampSpecOut[0];
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      if( ampSpecOut[i] > maxAmpOut ) maxAmpOut = ampSpecOut[i];
    }

    // Compute inverse (designature) filter = Inverse amplitude spectrum after adding white noise
    //    float ratio = (maxAmp/maxAmpOut);
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      myDesigAmpFilter[i] =  ampSpecOut[i] / (myAmpSpecIn[i] + whiteNoise);
    }

    // Compute phase shift
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      //   float freq          = df * (float)(i);
      //   float phaseShift    = 2.0 * M_PI * freq * timeZero_s;
      myDesigPhaseShift[i] = fmod( phaseSpecOut[i] - myPhaseSpecIn[i] , 2.0*M_PI );
      if( myDesigPhaseShift[i] > M_PI ) {
        myDesigPhaseShift[i] -= 2.0*M_PI;
      }
    }
  }
}
//--------------------------------------------------------------------------------
//
// TODO: Apply optional taper to data before FFT transform
//
bool csFFTDesignature::applyFilter( float* samples, int numSamples ) {
  if( numSamples != myNumSamplesIn ) {
    throw(csException("csFFTDesignature::applyFilter: Inconsistent number of samples: %d != %d", numSamples, myNumSamplesOut));
  }

  // Forward transform input data and compute amplitude & phase spectrum
  bool success = fft_forward( samples, myAmpSpecIn, myPhaseSpecIn );
  if( !success ) return success;

  if( myFilterType == csFFTDesignature::AMP_PHASE ) {
    // Apply inverse filter & zero-phasing to amplitude & phase spectrum
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      myAmpSpecIn[i]   *= myDesigAmpFilter[i];
      myPhaseSpecIn[i] += myDesigPhaseShift[i];
    }
  }
  else if( myFilterType == csFFTDesignature::AMP_ONLY ) {
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      myAmpSpecIn[i]   *= myDesigAmpFilter[i];
    }
  }
  else if( myFilterType == csFFTDesignature::PHASE_ONLY ) {
    for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
      myPhaseSpecIn[i] += myDesigPhaseShift[i];
    }
  }

  /*  if( myIsNotchSuppression ) {
    for( int isamp = myNotchIndexFirstRed; isamp <= myNotchIndexLastRed; isamp++ ) {
      double phase  = 2.0*( ( (double)(isamp-myNotchIndexFirst) / (double)myNotchWidth ) - 1.0 ) * M_PI;
      double scalar = 0.5 * (cos(phase) + 1.0);
      myAmpSpecIn[isamp] *= scalar;
    }
    } */

  // Transform back to X-T. Do not perform normalisation in fft (not sure how when normalisation is required..?)
  success = fft_inverse( myAmpSpecIn, myPhaseSpecIn, false );
  for( int i = 0; i < numSamples; i++ ) {
    samples[i] = myBufferReal[i];
  }

  return success;
}

void csFFTDesignature::setDesigFilterType( int filterType ) {
  myFilterType = filterType;
}

void csFFTDesignature::setDesigLowPass( float cutOffHz, float order ) {
  double G0 = 1.0;
  double power = order * 2.0;
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );

  for( int is = 0; is <= myNumFFTSamplesIn/2; is++ ) {
    double freq = (double)is * df;
    double dampG = sqrt( G0 / (1.0 + pow(freq/cutOffHz,power) ) );
    myDesigAmpFilter[is] *= dampG;
  }
}

void csFFTDesignature::setDesigHighPass( float cutOffHz, float order ) {
  double G0 = 1.0;
  double power = order * 2.0;
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );

  myDesigAmpFilter[0] = 0.0;

  for( int is = 1; is <= myNumFFTSamplesIn/2; is++ ) {
    double freq = (double)is * df;
    double dampG;
    dampG = sqrt( G0 / (1.0 + pow(cutOffHz/freq,power) ) );
    myDesigAmpFilter[is] *= dampG;
  }
}

void csFFTDesignature::setDesigHighEnd( float freq ) {
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );
  int index = (int)round( freq / df );
  if( index <= 0 || index > myNumFFTSamplesIn/2 ) return;
  float amplitude = myDesigAmpFilter[index];
  float phase     = myDesigPhaseShift[index];

  for( int is = index+1; is <= myNumFFTSamplesIn/2; is++ ) {
    myDesigAmpFilter[is]  = amplitude;
    myDesigPhaseShift[is] = phase;
  }
}

// Apply cosine taper around notch frequency
void csFFTDesignature::setNotchSuppression( float notchFreq, float notchWidth ) {
  double df = 1.0 / ( (double)myNumFFTSamplesIn*mySampleIntIn/1000.0 );
  int indexFirst = (int)round( ( notchFreq - 0.5*notchWidth ) / df );
  int indexLast  = (int)round( ( notchFreq + 0.5*notchWidth ) / df );
  int width = indexLast - indexFirst;

  int indexFirstRed = std::max(0,indexFirst);
  int indexLastRed  = std::min(myNumFFTSamplesIn,indexLast);
  /*
  myIsNotchSuppression = true;
  myNotchIndexFirst    = indexFirst;
  myNotchWidth = width;
  myNotchIndexFirstRed = indexFirstRed;
  myNotchIndexLastRed  = indexLastRed;
  */
  for( int isamp = indexFirstRed; isamp <= indexLastRed; isamp++ ) {
    double phase = 2.0*( ( (double)(isamp-indexFirst) / (double)width ) - 1.0 ) * M_PI;
    double scalar   = 0.5 * (cos(phase) + 1.0);
    myDesigAmpFilter[isamp] = myDesigAmpFilter[isamp] * scalar;
  }
}

//--------------------------------------------------------------------------------

void csFFTDesignature::dump_spectrum( FILE* stream ) const {
  float df = 1000.0 / ( (float)myNumFFTSamplesIn*mySampleIntIn );
  for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
    float freq = df * (float)(i);
    fprintf(stream,"%.6f  %.10e %.10e\n", freq, myDesigAmpFilter[i], myDesigPhaseShift[i] );
  }
}

void csFFTDesignature::dump_wavelet( FILE* stream, bool doNormalize ) {
  float* desigAmpFilter   = new float[myNumFFTSamplesIn/2+1];
  float* desigPhaseShift = new float[myNumFFTSamplesIn/2+1];

  float phaseScalar = 0.5 * myNumFFTSamplesIn * mySampleIntIn / 1000.0; // Linear phase shift to center main wavelet energy
  float df = 1000.0 / ( (float)myNumFFTSamplesIn*mySampleIntIn );
  for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
    float freq = df * (float)(i);
    float phaseShift    = freq * M_PI * phaseScalar;
    desigAmpFilter[i]  = myDesigAmpFilter[i];
    desigPhaseShift[i] = myDesigPhaseShift[i] + phaseShift;
  }
  convertFromAmpPhase( desigAmpFilter, desigPhaseShift );
  // bool success = 
  fft( csFFTTools::INVERSE, myTwoPowerIn, myBufferReal, myBufferImag, false );
  double const* real = realData();
  float normScalar = 1.0f;
  if( doNormalize ) normScalar = 1.0f / (float)(myNumFFTSamplesIn/2);
  for( int i = 0; i <= myNumFFTSamplesIn/2; i++ ) {
    float time = mySampleIntIn * (float)i;
    fprintf(stream,"%.6f  %.10e\n", time, normScalar*real[i]  );
  }
  delete [] desigAmpFilter;
  delete [] desigPhaseShift;
}

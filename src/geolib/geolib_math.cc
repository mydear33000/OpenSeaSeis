/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cmath>
#include <cstdio>
#include <limits>
#include "geolib_math.h"
#include "geolib_methods.h"

namespace cseis_geolib {

//---------------------------------------------
  int compute_correlation_length( int maxlag ) {
    return( 2*maxlag+1 );
  }

  /*
   *
   * auto[f(i)] = a[tau] = a[dt * i] = 1 / (m - i + 1) * SUM_j=0-m[ f(j) * f(i+j) ]
   */
  void compute_onesided_auto_correlation( float const* samples, int nSampIn, float* autocorr ) {
    compute_onesided_auto_correlation( samples, nSampIn, autocorr, nSampIn-1, false );
  }
  void compute_onesided_auto_correlation( float const* samples, int nSampIn, float* autocorr, int maxlag_in_num_samples ) {
    compute_onesided_auto_correlation( samples, nSampIn, autocorr, maxlag_in_num_samples, false );
  }
  void compute_onesided_auto_correlation( float const* samples, int nSampIn, float* autocorr, int maxlag_in_num_samples,
					  bool dampen ) {
    for( int ilag = 0; ilag <= maxlag_in_num_samples; ilag++ ) {
      float sum = 0;
      int sampEnd = nSampIn-ilag;
      for( int isamp = 0; isamp < sampEnd; isamp++ ) {
        sum += samples[isamp]*samples[isamp+ilag];
      }
      if( !dampen ) {
	autocorr[ilag] = sum;
      }
      else {
	autocorr[ilag] = (float)sampEnd * sum / (float)nSampIn;
      }
    }
  }

  
  void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight,
                                     int nSampIn, float* corr ) {
    compute_twosided_correlation( samplesLeft, samplesRight, nSampIn, corr, nSampIn-1, false );
  }
  void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight,
				     int nSampIn, float* corr, int maxlag_in_num_samples ) {
    compute_twosided_correlation( samplesLeft, samplesRight, nSampIn, corr, maxlag_in_num_samples, false );
  }
  void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight,
				     int nSampIn, float* corr, int maxlag_in_num_samples, bool dampen ) {

    if( !dampen ) {
      //---------------------------------------
      // Compute negative lags
      //
      int sampEnd    = nSampIn;
      for( int ilag = -maxlag_in_num_samples; ilag < 0; ilag++ ) {
        int sampStart  = -ilag;
        float sum = 0;
        for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
          sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
        }
        corr[ilag+maxlag_in_num_samples] = sum;
      }
      int sampStart = 0;
      //---------------------------------------
      // Compute positive lags
      //
      for( int ilag = 0; ilag <= maxlag_in_num_samples; ilag++ ) {
        int sampEnd    = nSampIn-ilag;
        float sum = 0;
        for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
          sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
        }
        corr[ilag+maxlag_in_num_samples] = sum; 
      }
    }
    else { // dampen
      //---------------------------------------
      // Compute negative lags
      //
      int sampEnd    = nSampIn;
      for( int ilag = -maxlag_in_num_samples; ilag < 0; ilag++ ) {
        int sampStart  = -ilag;
        float sum = 0;
        for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
          sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
        }
        int nSamp = sampEnd-sampStart;
        corr[ilag+maxlag_in_num_samples] = (float)nSamp * sum / (float)nSampIn;
      }
      int sampStart  = 0;
      //---------------------------------------
      // Compute positive lags
      //
      for( int ilag = 0; ilag <= maxlag_in_num_samples; ilag++ ) {
        int sampEnd    = nSampIn-ilag;
        float sum = 0;
        for( int isamp = sampStart; isamp < sampEnd; isamp++ ) {
          sum += samplesLeft[isamp]*samplesRight[isamp+ilag];
        }
        int nSamp = sampEnd-sampStart;
        corr[ilag+maxlag_in_num_samples] = (float)nSamp * sum / (float)nSampIn;
      }
    }
  }

  void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight,
                                     int nSampIn, float* corr, int maxlag_in_num_samples, float* sampleIndex_maxAmp, float* maxAmp ) {

    compute_twosided_correlation( samplesLeft, samplesRight, nSampIn, corr, maxlag_in_num_samples );

    // Determine maximum cross-correlation lag time & amplitude
    int nSampCorr = maxlag_in_num_samples*2 + 1;
    int sampleIndex_int = 0;
    *maxAmp = corr[sampleIndex_int];
    for( int isamp = 0; isamp < nSampCorr; isamp++ ) {
      if( corr[isamp] > *maxAmp ) {
        *maxAmp = corr[isamp];
        sampleIndex_int = isamp;
      }
    }
    *sampleIndex_maxAmp = getQuadMaxSample( corr, sampleIndex_int, nSampCorr, maxAmp );
  }

  float compute_rms( float const* samples, int nSamples ) {
    if( nSamples > 0 ) {
      double sum_sqr = 0.0;
      for( int isamp = 0; isamp < nSamples; isamp++ ) {
        sum_sqr += (double)samples[isamp]*(double)samples[isamp];
      }
      return( sqrt( (float)sum_sqr/(float)nSamples ) );
    }
    else {
      return 0.0;
    }
  }

  int factor_2357( int *my_in ) {
    int  test,  n2, n3, n5, n7, count; 
    
    int in = *my_in;
    int out; 
    
    out = -1;
    
    if ( in < 0 || in > (std::numeric_limits<int>::max()-1) ) return out; // Fail if input bad  
    if ( in <=1 ){ out = 1; return out; } // Skip the loop for the easy one. 
    
    // Factor-out small primes from the input number. It should reduce to 1 once 
    // it is a perfect factor of 2,3,5 and 7. If not, increment the input and repeat.
    // Fail if the number of iterations or the test value gets to large. 
    test  = 0;
    count = 0;
    while ( test != 1 && count < std::numeric_limits<int>::max() ){
      test = in+count;
      if ( test == std::numeric_limits<int>::max() ) break;
      count++;     
      
      n2=n3=n5=n7=0;
      while ( test%7 == 0 ){ test = test/7; n7++; }
      while ( test%5 == 0 ){ test = test/5; n5++; }
      while ( test%3 == 0 ){ test = test/3; n3++; }
      while ( test%2 == 0 ){ test = test/2; n2++; }
      
    }
    if ( test != 1 ) return out;
    out = (int)round( pow(2.0,n2) * pow(3.0,n3) * pow(5.0,n5) * pow(7.0,n7) );
    
    return out;
  }

}


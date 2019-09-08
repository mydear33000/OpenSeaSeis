/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

/**
 * Math 'library' methods
 *
 * @author Bjorn Olofsson
 * @date 2007
 */

#ifndef GEOLIB_MATH_H
#define GEOLIB_MATH_H

namespace cseis_geolib {

#define ERROR    0
#define SUCCESS  1

#define REMOVE 11
#define APPLY  12

#define CS_SQR(a) ((a)*(a))

#define SPECIAL_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

// The following macros are here for historical reasons
#define IMIN(a,b) ((a) < (b) ? (a) : (b))

template<class A> A DEG2RAD( A val ) {
  return( val / 57.29577951308232286 );
}
template<class A> A RAD2DEG( A val ) {
  return( val * 57.29577951308232286 );
}
template<class A> A SIGN( A val ) {
  return( val >= 0 ? static_cast<A>(1) : static_cast<A>(-1) );
}
template<class A> int NINT( A val ) {
  return( (val >= 0) ? (int)( (double)val + 0.55555555 ) : ( (int)( (double)val - 0.55555555 ) ) );
}
template<class A> static A MAX( A val1, A val2 ) {
  return( val1 > val2 ? val1 : val2 );
}
template<class A> static A MIN( A val1, A val2 ) {
  return( val1 < val2 ? val1 : val2 );
}

 /**
 * Compute correlation length for specified maximum time lag
 * This function shall guarantee consistency of length of correlation function when using the following methods
 * @return correlation length
 */
 int compute_correlation_length( int maxlag );
 /**
  * Compute correlation between two series (left & right)
  * Number of output samples is determined from maximum lag:
  *   nSamplesOut = 2*maxlag_in_num_samples+1
  *
  * @param samplesLeft  (i) 'Left' series
  * @param samplesRight (i) 'Right' series
  * @param nSampIn      (i) Number of values in input series
  * @param corr         (o) Ouput correlation function
  * @param maxlag_in_num_samples  (i) Maximum time lag to compute, in number of samples
  */
 void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight,
                                    int nSampIn, float* corr );
 void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight,
                                    int nSampIn, float* corr, int maxlag_in_num_samples );
 void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight,
                                    int nSampIn, float* corr, int maxlag_in_num_samples, bool dampen );
 /**
  * Compute correlation between two series (left & right)
  * Number of output samples is determined from maximum lag:
  *   nSamplesOut = 2*maxlag_in_num_samples+1
  *
  * @param samplesLeft  (i) 'Left' series
  * @param samplesRight (i) 'Right' series
  * @param nSampIn      (i) Number of values in input series
  * @param corr         (o) Ouput correlation function
  * @param maxlag_in_num_samples  (i) Maximum time lag to compute, in number of samples
  * @param sampleIndex_maxAmp (o) Interpolated maximum correlation lag (in samples)
  * @param maxAmp             (o) Cross-correlation function amplitude at maximum correlation lag
  */
 void compute_twosided_correlation( float const* samplesLeft, float const* samplesRight,
                                    int nSampIn, float* corr, int maxlag_in_num_samples, float* sampleIndex_maxAmp, float* maxAmp );

 void compute_onesided_auto_correlation( float const* samples, int nSampIn, float* autocorr );
 void compute_onesided_auto_correlation( float const* samples, int nSampIn, float* autocorr, int maxlag_in_num_samples );
 void compute_onesided_auto_correlation( float const* samples, int nSampIn, float* autocorr, int maxlag_in_num_samples, bool dampen );

 /**
 * Compute RMS value for given data series
 * @param samples  Input samples
 * @param nSamples Number of input samples
 * @return RMS value
 */
 float compute_rms( float const* samples, int nSamples );

////////////////////////////////////////////////////////////////////////////////
// Return the first multiple of small primes (2,3,5 and 7) greater than or 
// equal to the input value. 
//
// The return value is -1 if the value can't be reduced to a factor of small
// primes. That is, if it is less than 0 or too large (*).
//
// (*)NOTE: The largest 32-bit integer that is still a factor of small primes 
// appears to be 
//
//          INT_MAX - 3330622  = 2147483647 - 3330622
//                             = 2144153025
//
// Numbers larger than this are still valid integers up to INT_MAX but can't 
// be factored.
// 
// It is up to the calling program to decide what to do in that case. The return
// value will still be -1.
//
////////////////////////////////////////////////////////////////////////////////
 int factor_2357( int *my_in );

} // namespace

#endif


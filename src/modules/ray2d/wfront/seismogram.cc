/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

void addRickerWavelet( float freq_hz, float dampingFactor, float phaseShift_rad, float sampleInt_s, float timeSamp1_s,
                       float eventTime_s, float eventPhase_rad, float eventAmplitude, float* samples, int nSamples )
{
  if( dampingFactor < 0.0001f ) {
    dampingFactor = 1.0f;
  }

  float omegaHalf    = (float)( M_PI * freq_hz );
  float exponentDamped_sq = omegaHalf * omegaHalf / (dampingFactor * dampingFactor);  // = og
  //  float omegaHalf_sq = omegaHalf * omegaHalf;

  float waveletWidth_s = (float)( dampingFactor * 3.0 * sqrt(6) / M_PI ) / freq_hz;
  float waveletWidthHalf_s = 0.5f * waveletWidth_s;

  float time1 = eventTime_s - waveletWidthHalf_s - timeSamp1_s;
  float time2 = eventTime_s + waveletWidthHalf_s - timeSamp1_s;

  /*
  float maxTime = (float)nSamples * sampleInt_s;
  if( time1 >= maxTime || time2 <= 0.0 ) {
    fprintf(stderr,"Outside of valid window. Max time: %.2fms, time1: %.2fms, time2: %.2f\n", maxTime * 1000.0, time1*1000, time2*1000);
    return;
  }
  */

  int index1 = std::max( (int)round( time1 / sampleInt_s ), 0 );
  int index2 = std::min( (int)round( time2 / sampleInt_s ), nSamples-1 );
    fprintf(stderr,"Start/end window: %d %d\n", index1, index2);

  for( int kk = index1; kk <= index2; kk++ ) {
    float time_s = timeSamp1_s + sampleInt_s * (float)(kk);
    float td = timeSamp1_s + sampleInt_s * (float)(kk) - eventTime_s;
    float phase     = phaseShift_rad - eventPhase_rad + 2.0f * omegaHalf * td;
    float amplitude = eventAmplitude * cos(phase) * exp( -exponentDamped_sq * td * td );
    samples[kk] += amplitude;

    //    float tmp_sq  = omegaHalf_sq * td * td;
    //    float tmp2_sq   = exponentDamped_sq * td * td;
    //    float amplitude2 = eventAmplitude * (1.0 - 2.0 * tmp_sq ) * exp( -tmp2_sq );
    //    fprintf(stderr," %f %e %f %f\n", time_s, amplitude, cos(phase), (1.0 - 2.0 * tmp_sq ));
  }
 
}



#ifndef CS_MODULES_H
#define CS_MODULES_H


const int N_METHODS_SINGLE = 47;
const int N_METHODS_MULTI  = 36;
const int N_METHODS = N_METHODS_MULTI + N_METHODS_SINGLE;

std::string NAMES[N_METHODS] = {
  std::string("ATTRIBUTE"),
  std::string("BEAM_FORMING"),
  std::string("BIN"),
  std::string("CCP"),
  std::string("CMP"),
  std::string("CONCATENATE"),
  std::string("CONVOLUTION"),
  std::string("CORRELATION"),
  std::string("DEBIAS"),
  std::string("DESIGNATURE"),
  std::string("DESPIKE"),
  std::string("ELSE"),
  std::string("ELSEIF"),
  std::string("ENDIF"),
  std::string("ENDSPLIT"),
  std::string("ENS_DEFINE"),
  std::string("FFT"),
  std::string("FFT_2D"),
  std::string("FILTER"),
  std::string("FXDECON"),
  std::string("GAIN"),
  std::string("GEOTOOLS"),
  std::string("HDR_DEL"),
  std::string("HDR_MATH"),
  std::string("HDR_MATH_ENS"),
  std::string("HDR_PRINT"),
  std::string("HDR_SET"),
  std::string("HISTOGRAM"),
  std::string("HODOGRAM"),
  std::string("IF"),
  std::string("IMAGE"),
  std::string("INPUT"),
  std::string("INPUT_ASCII"),
  std::string("INPUT_CREATE"),
  std::string("INPUT_RSF"),
  std::string("INPUT_SEGD"),
  std::string("INPUT_SEGY"),
  std::string("INPUT_SINEWAVE"),
  std::string("KILL"),
  std::string("KILL_ENS"),
  std::string("LMO"),
  std::string("MIRROR"),
  std::string("MUTE"),
  std::string("NMO"),
  std::string("OFF2ANGLE"),
  std::string("ORIENT"),
  std::string("ORIENT_CONVERT"),
  std::string("OUTPUT"),
  std::string("OUTPUT_RSF"),
  std::string("OUTPUT_SEGY"),
  std::string("OVERLAP"),
  std::string("P190"),
  std::string("PICKING"),
  std::string("POSCALC"),
  std::string("PZ_SUM"),
  std::string("RAY2D"),
  std::string("READ_ASCII"),
  std::string("REPEAT"),
  std::string("RESAMPLE"),
  std::string("RESEQUENCE"),
  std::string("RMS"),
  std::string("ROTATE"),
  std::string("SCALING"),
  std::string("SELECT"),
  std::string("SELECT_TIME"),
  std::string("SEMBLANCE"),
  std::string("SORT"),
  std::string("SPLIT"),
  std::string("SPLITTING"),
  std::string("STACK"),
  std::string("STATICS"),
  std::string("SUMODULE"),
  std::string("TEST"),
  std::string("TEST_MULTI_ENSEMBLE"),
  std::string("TEST_MULTI_FIXED"),
  std::string("TIME_SLICE"),
  std::string("TIME_STRETCH"),
  std::string("TRC_ADD_ENS"),
  std::string("TRC_INTERPOL"),
  std::string("TRC_MATH"),
  std::string("TRC_MATH_ENS"),
  std::string("TRC_PRINT"),
  std::string("TRC_SPLIT")
};

#endif

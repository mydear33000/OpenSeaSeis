/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "cseis_includes.h"
#include "csVector.h"
#include "csSortManager.h"
#include "csFlexNumber.h"
#include "csTime.h"
#include <cstring>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: READ_ASCII
 *
 * @author Bjorn Olofsson
 * @date   2007
 *
 * Bug fixes and updates
 *  2007-Oct-08 - Convert key position to C++ (start at 0). Correct key & header position during input
 *  2009-Apr-05 - Change second position given in parameter 'header' and 'key' to position (used to be length, wrongly documented)
 */
namespace mod_read_ascii {
  struct VariableStruct {
    int npos;             // Number of lines in input ASCII file
    int currentPos;       // Index of line from which last value was retrieved
    double** keyValues;   // Key values as they appear in input ASCII file
    double** headerValues;
    double* keyDoubleBuffer;  // Key values of current trace
    int* keyIntBuffer;        // Key values of current trace
    std::string* keyStringBuffer;
    cseis_geolib::csVector<int>* headerIndexList;
    cseis_geolib::csVector<int>* headerTypeList;
    cseis_geolib::csVector<int>* keyIndexList;
    cseis_geolib::csVector<int>* keyTypeList;
    char  ignoreChar;
    char  selectChar;
    int numKeys;
    int numHeaders;
    bool dropUnmatchedTraces;
    bool showWarnings;

    bool      isTimeKey;
    bool      isTimeKey_us;
//    csDate_t  date;
    int       timeKeyYear;
    int       hdrId_time_key;
  };

// Do not change index number! ...gives number of parameters for each method
  static int const METHOD_COLUMNS   = 2;
  static int const METHOD_POSITIONS = 3;

struct Pos {
  int start;
  int length;
};

}
using namespace mod_read_ascii;

//*************************************************************************************************
// Init phase
//
//
//*************************************************************************************************
void init_mod_read_ascii_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_SINGLETRACE );

  vars->npos       = 0;
  vars->currentPos = 0;
  vars->keyValues       = NULL;
  vars->headerValues    = NULL;
  vars->keyDoubleBuffer = NULL;
  vars->keyIntBuffer    = NULL;
  vars->keyStringBuffer = NULL;
  vars->headerIndexList = NULL;
  vars->headerTypeList  = NULL;
  vars->keyIndexList    = NULL;
  vars->keyTypeList     = NULL;
  vars->ignoreChar   = ' ';
  vars->selectChar   = ' ';
  vars->numKeys      = 0;
  vars->numHeaders   = 0;
  vars->dropUnmatchedTraces = false;
  vars->showWarnings = false;

  vars->isTimeKey    = false;
  vars->isTimeKey_us = false;
  vars->hdrId_time_key = -1;

  //----------------------------------------------------

  std::string filename;
  std::string methodText;
  csVector<std::string> valueList;
  csVector<std::string> headerNameList;
  csVector<int> headerColumnList;
  csVector<Pos> headerPosList;

  csVector<std::string> keyNameList;
  csVector<int> keyColumnList;
  csVector<Pos> keyPosList;

  std::string headerName;
  std::string keyName;
  int column;

  FILE* f_in;

  vars->headerIndexList = new csVector<int>();
  vars->headerTypeList  = new csVector<int>();
  vars->keyIndexList    = new csVector<int>();
  vars->keyTypeList     = new csVector<int>();
  //---------------------------------------------------------
  // Read in general parameters
  //
  int method = -1;

  param->getString( "filename", &filename );
  param->getString( "method", &methodText );
  methodText = toLowerCase( methodText );

  if( !methodText.compare("columns") ) {
    method = METHOD_COLUMNS;
  }
  else if( !methodText.compare("positions") ) {
    method = METHOD_POSITIONS;
  }
  else {
    log->line("Option for parameter METHOD not recognised: %s", methodText.c_str());
    env->addError();
  }

  std::string text;
  vars->showWarnings = false;
  if( param->exists("warnings") ) {
    param->getString( "warnings", &text );
    if( !text.compare("yes") ) {
      vars->showWarnings = true;
    }
    else if( !text.compare("no") ) {
      vars->showWarnings = false;
    }
    else {
      log->line("Option for parameter 'warnings' not recognised: %s", text.c_str());
    }
  }

  bool checkConsistency = false;
  if( param->exists("check") ) {
    param->getString( "check", &text );
    if( !text.compare("yes") ) {
      checkConsistency = true;
    }
    else if( !text.compare("no") ) {
      checkConsistency = false;
    }
    else {
      log->line("Option not recognised: %s", text.c_str());
    }
  }

//-------------------------------------------
//
  Pos timeKeyPos;
  if( param->exists("key_sps_time") ) {
    param->getInt( "time_year", &vars->timeKeyYear );
    string timeKeyName;
    vars->isTimeKey = true;
    int numValues = param->getNumValues("key_sps_time");
    if( numValues < 2 ) {
      log->error("Parameter key_sps_time expects two or three values: Header name, column or start position, and length");
    }
    param->getString( "key_sps_time", &timeKeyName, 0 );
    param->getInt( "key_sps_time", &timeKeyPos.start, 1 );
    timeKeyPos.start -= 1;  // Convert to C-type position index
    if( method == METHOD_POSITIONS ) {
      if( numValues < 3 ) {
        log->error("Parameter key_sps_time expects three values: Header name, column or start position, and length");
      }
      param->getInt( "key_sps_time", &timeKeyPos.length, 2 );
      if( timeKeyPos.start < 0 || timeKeyPos.length <= 0 ) {
        log->line("Error: Inconsistent position/length given for time key %s: start=%d, length=%d",
                  timeKeyName.c_str(), timeKeyPos.start+1, timeKeyPos.length);
        env->addError();
      }
    }
    vars->hdrId_time_key = hdef->headerIndex( timeKeyName.c_str() );
    vars->keyIndexList->insertEnd( hdef->headerIndex(timeKeyName) );
    vars->keyTypeList->insertEnd( hdef->headerType(timeKeyName) );

    if( param->exists("key_sps_time_us") ) {
      param->getString( "key_sps_time_us", &timeKeyName, 0 );
      vars->isTimeKey_us = true;
      vars->keyIndexList->insertEnd( hdef->headerIndex(timeKeyName) );
      vars->keyTypeList->insertEnd( hdef->headerType(timeKeyName) );
    }

//    if( !text.compare("jjjhhmmss") ) {
//      vars->showWarnings = true;
//    }
//    else {
//      log->line("Specified time key format not supported: %s. Use jjjhhmmss instead", text.c_str());
//    }
//    vars->dateFormat = new char[20];
//    sprintf(vars->dateFormat,"%%3d%%2d%%2d%%2d");
  }

//-------------------------------------------
// Ignore/select certain lines... not implemented yet
  bool isIgnore = false;
  bool isSelect = false;
  if( param->exists("ignore_line") ) {
    valueList.clear();
    param->getAll( "ignore_line", &valueList );
    if( valueList.size() > 1 ) log->error("Only one 'ignore' character supported...");
    isIgnore = true;
    vars->ignoreChar = valueList.at(0).at(0);
    if( edef->isDebug() ) log->line("Ignore character: '%c'", vars->ignoreChar );
  }

//-------------------------------------------

  if( param->exists("select_line") ) {
    valueList.clear();
    param->getAll( "select_line", &valueList );
    if( valueList.size() > 1 ) log->error("Only one 'select' character supported...");
    isSelect = true;
    vars->selectChar = valueList.at(0).at(0);
  }

//---------------------------------------------------------------
// Read in user parameters for headers
//
  int maxNumColumn = 0;
  int maxPosition  = 0;

  vars->numHeaders = param->getNumLines( "header" );
        //   = nLines;   // !CHANGE! Make sure numHeaders is used in all instances of the headers, not list.size() or similar...
  if( vars->numHeaders == 0 ) {
    log->line("User parameter 'header' missing. No headers specified which shall be read in.");
    env->addError();
  }
  for( int iLine = 0; iLine < vars->numHeaders; iLine++ ) {
    valueList.clear();
    param->getAll( "header", &valueList, iLine );
    if( valueList.size() != method ) {
      log->line("Wrong number of arguments for user parameter 'header'. Expected: %d, found: %d.",
                method, valueList.size());
      env->addError();
    }
    else {
      headerName = valueList.at(0);
      headerNameList.insertEnd( headerName );
      if( method == METHOD_COLUMNS ) {
        column = atoi(valueList.at(1).c_str()) - 1;  // C style column index (starting with 0)
        headerColumnList.insertEnd( column );
        if( column > maxNumColumn ) maxNumColumn = column;
      }
      else {
        Pos pos;
        pos.start  = atoi(valueList.at(1).c_str()) - 1;  // C style index (starting with 0)
        pos.length = atoi(valueList.at(2).c_str());
        if( pos.start < 0 || pos.length <= 0 ) {
          log->line("Error: Inconsistent positions given for header %s: start=%d, length=%d",
                    headerName.c_str(), pos.start+1, pos.length);
          env->addError();
        }
        headerPosList.insertEnd( pos );
        if( pos.start + pos.length > maxPosition ) maxPosition = pos.start + pos.length;
      }
    }
  }
  vars->numHeaders = headerNameList.size();
  if( vars->numHeaders == 0 ) log->error("No (valid) trace headers specified.");


//---------------------------------------------------------------
// Read in user parameters for keys
//
  vars->numKeys = param->getNumLines( "key" );
  if( vars->numKeys <= 0 && !vars->isTimeKey ) {
    log->error("User parameter KEY missing. At least one key needs to be specified.");
  }
  else if( vars->numKeys > 0 && vars->isTimeKey ) {
    log->error("When time key is used, no other key is supported.");
  }
  else if( vars->numKeys > 0 ) {
    vars->keyDoubleBuffer = new double[vars->numKeys];
    vars->keyIntBuffer    = new int[vars->numKeys];
//    vars->keyStringBuffer = new std::string[vars->numKeys];
    for( int ikey = 0; ikey < vars->numKeys; ikey++ ) {
      vars->keyIntBuffer[ikey] = 0;
      vars->keyDoubleBuffer[ikey] = 0.0;
//      vars->keyStringBuffer[ikey] = "";
    }

    for( int ikey = 0; ikey < vars->numKeys; ikey++ ) {
      if( edef->isDebug() ) log->line("Reading parameters for key %d", ikey+1);
      valueList.clear();
      param->getAll( "key", &valueList, ikey );
      if( valueList.size() != method ) {
        log->line("Wrong number of arguments for user parameter 'key'. Expected: %d, found: %d.",
                  method, valueList.size());
        env->addError();
      }
      else {
        keyName = valueList.at(0);
        keyNameList.insertEnd( keyName );
        if( method == METHOD_COLUMNS ) {
          column = atoi(valueList.at(1).c_str()) - 1;  // C style column index (starting with 0)
          keyColumnList.insertEnd( column );
          if( column > maxNumColumn ) maxNumColumn = column;
        }
        else {
          Pos pos;
          pos.start  = atoi(valueList.at(1).c_str()) - 1;  //  C style column index (starting with 0)
          pos.length = atoi(valueList.at(2).c_str());
          if( pos.start < 0 || pos.length <= 0 ) {
            log->line("Error: Inconsistent positions given for key %s: start=%d, length=%d",
                      keyName.c_str(), pos.start+1, pos.length);
            env->addError();
          }
          keyPosList.insertEnd( pos );
        }
      }
    }
  }


//---------------------------------------------------------------
  vars->dropUnmatchedTraces = false;
  if( param->exists("drop_traces") ) {
    std::string text;
    param->getString( "drop_traces", &text );
    if( !text.compare("yes") ) {
      vars->dropUnmatchedTraces = true;
    }
    else if( !text.compare("no") ) {
      vars->dropUnmatchedTraces = false;
    }
    else {
      log->error("Unknown option '%s'", text.c_str());
    }
  }
//---------------------------------------------------------------
// Check if specified headers exist. param->get header index & type
//
  for( int i = 0; i < vars->numHeaders; i++ ) {
    std::string headerName = headerNameList.at(i);
                bool headerExists = hdef->headerExists( headerName );
    if( headerExists || csStandardHeaders::isStandardHeader(headerName) ) {
                        if( !headerExists ) {
                                hdef->addStandardHeader( headerName );
                        }
      vars->headerIndexList->insertEnd( hdef->headerIndex(headerName.c_str()) );
      type_t type = hdef->headerType(headerName.c_str());
      if( type == TYPE_STRING ) {
        log->error("String headers are currently not supported for this module: Trace header %s", headerName.c_str());
      }
      vars->headerTypeList->insertEnd( type );
    }
    else {
      log->line("Trace header %s does not exist, and is not a standard header.", headerName.c_str());
      env->addError();
    }
  }
  for( int i = 0; i < vars->numKeys; i++ ) {
    std::string keyName = keyNameList.at(i);
    if( edef->isDebug() ) log->line("Checking parameters for key %d (%s)", i+1, keyName.c_str());
    if( hdef->headerExists( keyName ) ) {
      for( int k = 0; k < headerNameList.size(); k++ ) {
        if( !keyName.compare( headerNameList.at(k) ) ) {
          log->line("Trace header %s specified both as header and as key.", headerName.c_str());
          env->addError();
        }
      }
      vars->keyIndexList->insertEnd( hdef->headerIndex(keyName.c_str()) );
      vars->keyTypeList->insertEnd( hdef->headerType(keyName.c_str()) );
    }
    else {
      log->line("Trace header %s does not exist.", keyName.c_str());
      env->addError();
    }
  }
  if( env->errorCount() != 0 ) {
    return;
  }

//---------------------------------------------------------------
// Read in ASCII file
//
  if( (f_in = fopen( filename.c_str(), "r" )) == (FILE*) NULL ) {
    log->error("Could not open file: '%s'", filename.c_str());
  }

  char buffer[1024];
  int counterLines = 0;
  while( fgets( buffer, 1024, f_in ) != NULL ) {
    if( strlen(buffer) <= 1 ) continue;  // Ignore empty lines (Allow 1 character for trailing newline)
    if( isIgnore && buffer[0] == vars->ignoreChar ) continue;
    if( isSelect && buffer[0] != vars->selectChar ) continue;
    counterLines++;
  }
  vars->npos = counterLines;
  rewind( f_in );

  if( vars->numKeys > 0 ) {
    vars->keyValues = new double*[vars->numKeys];
    for( int i = 0; i < vars->numKeys; i++ ) {
      vars->keyValues[i] = new double[vars->npos];
    }
  }
  else {  // Time key:
    if( vars->isTimeKey_us ) {
      vars->keyValues    = new double*[2];
      vars->keyValues[1] = new double[vars->npos];  // Allocate second key value 'just in case' microseconds are required
    }
    else {
      vars->keyValues    = new double*[1];
    }
    vars->keyValues[0] = new double[vars->npos];
  }
  vars->headerValues = new double*[vars->numHeaders];
  for( int i = 0; i < vars->numHeaders; i++ ) {
    vars->headerValues[i] = new double[vars->npos];
  }

  csDate_t date;
  date.year = vars->timeKeyYear;


  counterLines = 0;

  if( method == METHOD_COLUMNS ) {
    if( vars->isTimeKey ) {
      if( timeKeyPos.start > maxNumColumn ) {
        maxNumColumn = timeKeyPos.start;
      }
    }
    while( fgets( buffer, 1024, f_in ) != NULL ) {
      int bufferLength = strlen(buffer);
      if( bufferLength <= 1 ) continue;  // Ignore empty lines (Allow 1 character for trailing newline)
      if( isIgnore && buffer[0] == vars->ignoreChar ) continue;
      if( isSelect && buffer[0] != vars->selectChar ) continue;
      if( edef->isDebug() ) log->line("ASCII file, line #%3d: %s", counterLines, buffer);
      valueList.clear();
      tokenize( buffer, valueList );
      if( valueList.size() < maxNumColumn+1 ) {
        log->line("Input file contains too few columns. Number of columns found: %d. Maximum column number for key/header: %d. First line:\n%s", valueList.size(), maxNumColumn+1, buffer);
        env->addError();
        break;
      }
      for( int i = 0; i < vars->numKeys; i++ ) {
        column = keyColumnList.at(i);
        vars->keyValues[i][counterLines] = atof(valueList.at(column).c_str());
      }
      for( int ihdr = 0; ihdr < vars->numHeaders; ihdr++ ) {
        column = headerColumnList.at(ihdr);
        vars->headerValues[ihdr][counterLines] = atof(valueList.at(column).c_str());
      }
      if( vars->isTimeKey ) {
        string dateString( valueList.at(timeKeyPos.start) );
        date.julianDay = atoi( dateString.substr(0,3).c_str() );
        date.hour      = atoi( dateString.substr(3,2).c_str() );
        date.min       = atoi( dateString.substr(5,2).c_str() );
        date.sec       = atoi( dateString.substr(7,2).c_str() );
        vars->keyValues[0][counterLines] = (double)date.unixTime();
        if( vars->isTimeKey_us ) {
          date.usec = atoi( dateString.substr(10,6).c_str() );
          vars->keyValues[1][counterLines] = (double)date.usec;
        }
        if( edef->isDebug() ) log->line("Date: %s", date.getString() );
      }
      counterLines++;
    }
  }
  else if( method == METHOD_POSITIONS ) {
    while( fgets( buffer, 1024, f_in ) != NULL ) {
      int bufferLength = strlen(buffer);
      if( bufferLength <= 1 ) continue;  // Ignore empty lines (Allow 1 character for trailing newline)
      if( isIgnore && buffer[0] == vars->ignoreChar ) continue;
      if( isSelect && buffer[0] != vars->selectChar ) continue;

      valueList.clear();
      std::string bufferString = buffer;
      if( edef->isDebug() ) {
        log->line("ASCII file, line #%3d: %s", counterLines, bufferString.c_str());
      }
      if( bufferLength < maxPosition ) {
        log->line("Input line contains too few characters. Expected, according to specified key/header positions: %d, found: %d\nLine: %s",
                  maxPosition, bufferLength, bufferString.c_str() );
        env->addError();
        break;
      }
      for( int i = 0; i < vars->numKeys; i++ ) {
        Pos pos = keyPosList.at(i);
        if( pos.start+pos.length > bufferLength ) {
          log->line("Error: Key %s, ASCII file line #%d: Start position/length exceeds line length", keyNameList.at(i).c_str(), counterLines+1 );
          env->addError();
        }
        vars->keyValues[i][counterLines] = atof(bufferString.substr(pos.start,pos.length).c_str());
      }
      for( int ihdr = 0; ihdr < vars->numHeaders; ihdr++ ) {
        Pos pos = headerPosList.at(ihdr);
        if( pos.start+pos.length > bufferLength ) {
          log->line("Error: Header %s, ASCII file line #%d: Start position/length exceeds line length", headerNameList.at(ihdr).c_str(), counterLines+1 );
          env->addError();
        }
        vars->headerValues[ihdr][counterLines] = atof(bufferString.substr(pos.start,pos.length).c_str());
      }
      if( vars->isTimeKey ) {
        if( timeKeyPos.start+timeKeyPos.length > bufferLength ) {
          log->line("Error: Time key, ASCII file line #%d: Start position/length exceeds line length", counterLines+1 );
          env->addError();
        }
        string dateString = bufferString.substr(timeKeyPos.start,timeKeyPos.length);
        date.julianDay = atoi( dateString.substr(0,3).c_str() );
        date.hour      = atoi( dateString.substr(3,2).c_str() );
        date.min       = atoi( dateString.substr(5,2).c_str() );
        date.sec       = atoi( dateString.substr(7,2).c_str() );
        vars->keyValues[0][counterLines] = (double)date.unixTime();
        if( vars->isTimeKey_us ) {
          date.usec = atoi( dateString.substr(10,6).c_str() );
          vars->keyValues[1][counterLines] = (double)date.usec;
        }
        if( edef->isDebug() ) log->line("Date: %s   %d", date.getString(), date.unixTime() );
      }
      counterLines++;
    }
  }

  fclose( f_in );

//-----------------------------------------------------------------------------
// Set time key values
//
  if( vars->isTimeKey ) {
    vars->numKeys = 1;
    if( vars->isTimeKey_us ) vars->numKeys += 1;
    vars->keyDoubleBuffer = new double[vars->numKeys];
    vars->keyIntBuffer    = new int[vars->numKeys];
  }

//-----------------------------------------------------------------------------
// Sort key values
//
        bool doSort = false;
        if( param->exists("sort") ) {
                string text;
                param->getString("sort",&text);

                if( !text.compare("yes") ) {
                        doSort = true;
                }
                else if( !text.compare("no") ) {
                        doSort = false;
                }
                else {
                        log->error("Option not recognised: %s", text.c_str());
                }
        }

        if( doSort ) {

                cseis_geolib::csSortManager sortManager( vars->numKeys, csSortManager::TREE_SORT );
                sortManager.resetValues( vars->npos );
  
                for( int ikey = 0; ikey < vars->numKeys; ikey++ ) {
                        for( int ipos = 0; ipos < vars->npos; ipos++ ) {
                                double value = vars->keyValues[ikey][ipos];
                                sortManager.setValue( ipos, vars->numKeys-ikey-1, csFlexNumber(value) );
                        }
                }
  
                sortManager.sort();
  
                for( int ikey = 0; ikey < vars->numKeys; ikey++ ) {
                        double* sortedValues = new double[vars->npos];
                        for( int ipos = 0; ipos < vars->npos; ipos++ ) {
                                sortedValues[ipos] = vars->keyValues[ikey][ sortManager.sortedIndex(ipos) ];
                        }
                        delete [] vars->keyValues[ikey];
                        vars->keyValues[ikey] = sortedValues;
                }

                for( int ihdr = 0; ihdr < vars->numHeaders; ihdr++ ) {
                        double* sortedValues = new double[vars->npos];
                        for( int ipos = 0; ipos < vars->npos; ipos++ ) {
                                sortedValues[ipos] = vars->headerValues[ihdr][ sortManager.sortedIndex(ipos) ];
                        }
                        delete [] vars->headerValues[ihdr];
                        vars->headerValues[ihdr] = sortedValues;
                }
        }
//-----------------------------------------------------------------------------
//
  if( checkConsistency ) {
    bool abort = false;
    for( int ipos = 1; ipos < vars->npos; ipos++ ) {
      bool isSame = true;
      for( int ikey = 0; ikey < vars->numKeys; ikey++ ) {
        isSame = isSame && ( vars->keyValues[ikey][ipos] == vars->keyValues[ikey][ipos-1] );
      }
      if( isSame ) {
        log->warning("Duplicate key(s) in input file:");
        for( int ikey = 0; ikey < vars->numKeys; ikey++ ) {
          log->write("Key #%d: %f  ", ikey+1, vars->keyValues[ikey][ipos] );
        }
        log->write("\n");
        abort = true;
      }
    }
    if( abort ) {
      log->error("Inconsistencies found in input file.");
    }
  }

//-----------------------------------------------------------------------------
// Debug dump
//
  if( edef->isDebug() ) {
    for( int i = 0; i < vars->keyIndexList->size(); i++ ) {
      log->write( "%15s ", keyNameList.at(i).c_str() );
    }
    log->write(" ### ");
    for( int i = 0; i < vars->headerIndexList->size(); i++ ) {
      log->write( "%15s ", headerNameList.at(i).c_str() );
    }
    if( vars->isTimeKey ) {
      log->write( "%15s ", "Time key" );
    }
    log->line("");
    double value;
    for( int iPos = 0; iPos < vars->npos; iPos++ ) {
      for( int i = 0; i < vars->keyIndexList->size(); i++ ) {
        value = vars->keyValues[i][iPos];
        log->write( "%15.3f ", value );
      }
      log->write(" ### ");
      for( int i = 0; i < vars->headerIndexList->size(); i++ ) {
        value = vars->headerValues[i][iPos];
        log->write( "%15.3f ", value );
      }
      if( vars->isTimeKey ) {
        log->write( "%12d %s ", (int)vars->keyValues[0][iPos], csGeolibUtils::UNIXsec2dateString((int)vars->keyValues[0][iPos]).c_str() );
        if( vars->isTimeKey_us ) {
          log->write( "(%6d) ", (int)vars->keyValues[1][iPos] );
        }
      }
      log->line("");
    }
  }

  vars->currentPos = 0;

}

//*************************************************************************************************
// Exec phase
//
//
//
//*************************************************************************************************
bool exec_mod_read_ascii_(
  csTrace* trace,
  int* port,
  csExecPhaseEnv* env, csLogWriter* log )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;

  if( edef->isCleanup()){

    if( vars->keyIntBuffer ) {
      delete [] vars->keyIntBuffer;
      vars->keyIntBuffer = NULL;
    }
    if( vars->keyStringBuffer ) {
      delete [] vars->keyStringBuffer;
      vars->keyStringBuffer = NULL;
    }
    if( vars->keyDoubleBuffer ) {
      delete [] vars->keyDoubleBuffer;
      vars->keyDoubleBuffer = NULL;
    }
    if( vars->keyIndexList ) {
      if( vars->keyValues ) {
        for( int i = 0; i < vars->keyIndexList->size(); i++ ) {
          delete [] vars->keyValues[i];
        }
        delete [] vars->keyValues;
      }
      delete vars->keyIndexList;
    }
    if( vars->headerIndexList ) {
      if( vars->headerValues ) {
        for( int i = 0; i < vars->headerIndexList->size(); i++ ) {
          delete [] vars->headerValues[i];
        }
        delete [] vars->headerValues;
      }
      delete vars->headerIndexList;
    }
    if( vars->keyTypeList ) {
      delete vars->keyTypeList;
      vars->keyTypeList = NULL;
    }
    if( vars->headerTypeList ) {
      delete vars->headerTypeList;
      vars->headerTypeList = NULL;
    }
    delete vars; vars = NULL;
    return true;
  }

  csTraceHeader* trcHdr = trace->getTraceHeader();
//------------------------------------------------------------------
// Get key values from trace header
//
  for( int iKey = 0; iKey < vars->numKeys; iKey++ ) {
    int type = vars->keyTypeList->at(iKey);
    int keyIndex = vars->keyIndexList->at(iKey);
    switch( type ) {
      case TYPE_FLOAT:
      case TYPE_DOUBLE:
        vars->keyDoubleBuffer[iKey] = trcHdr->doubleValue( keyIndex );
        break;
      case TYPE_INT:
        vars->keyIntBuffer[iKey] = trcHdr->intValue( keyIndex );
        break;
      case TYPE_INT64:
        vars->keyDoubleBuffer[iKey] = (double)trcHdr->int64Value( keyIndex );
        break;
      case TYPE_STRING:
        vars->keyStringBuffer[iKey] = trcHdr->stringValue( keyIndex );
        break;
      default:
        log->line("ERROR!!!");
    }
  }

  if( edef->isDebug() ) {
    log->write("Keys: ");
    for( int iKey = 0; iKey < vars->numKeys; iKey++ ) {
      switch( vars->keyTypeList->at(iKey) ) {
        case TYPE_FLOAT:
        case TYPE_DOUBLE:
        case TYPE_INT64:
          log->write("#%d,float: %f, ", iKey, vars->keyDoubleBuffer[iKey] );
        break;
      case TYPE_INT:
          log->write("#%d,int: %d, ", iKey, vars->keyIntBuffer[iKey] );
        break;
      case TYPE_STRING:
          log->write("#%d,string: %s , ", iKey, vars->keyStringBuffer[iKey].c_str() );
      }
    }
    log->write("\n");
  }

  int counterKey = 0;

  vars->currentPos = 0;  // TEMP ...to make it work in a brute way

  while( counterKey < vars->numKeys ) {
    int type = vars->keyTypeList->at(counterKey);

    if( type == TYPE_FLOAT || type == TYPE_DOUBLE || type == TYPE_INT64 ) {
      while( vars->currentPos < vars->npos && vars->keyValues[counterKey][vars->currentPos] != vars->keyDoubleBuffer[counterKey] ) {
        vars->currentPos++;
      }
      if( vars->currentPos == vars->npos ) {
        vars->currentPos = 0;
      }
      while( vars->currentPos < vars->npos && vars->keyValues[counterKey][vars->currentPos] != vars->keyDoubleBuffer[counterKey] ) {
        vars->currentPos++;
      }
      if( vars->currentPos == vars->npos ) {
        if( vars->showWarnings ) {
          log->write("Unable to find matching line in ASCII file for ");
          for( int i = 0; i < vars->numKeys; i++ ) {
            log->write("key #%d: %f  ", i+1, vars->keyDoubleBuffer[i] );
          }
          log->write("\n");
        }
        return !vars->dropUnmatchedTraces;
      }
      else {
        for( int jKey = counterKey-1; jKey >= 0; jKey-- ) {
          if( vars->keyValues[jKey][vars->currentPos] != vars->keyDoubleBuffer[jKey] ) {
            if( vars->showWarnings ) {
              log->write("Unable to find matching line in ASCII file for ");
              for( int i = 0; i < vars->numKeys; i++ ) {
                log->write("key #%d: %f  ", i+1, vars->keyDoubleBuffer[i] );
              }
              log->write("\n");
            }
            return !vars->dropUnmatchedTraces;
          }
        }
      }

    }
    // TYPE_INT
    else if( type == TYPE_INT ) {
      while( vars->currentPos < vars->npos && vars->keyValues[counterKey][vars->currentPos] != (double)vars->keyIntBuffer[counterKey] ) {
        vars->currentPos++;
      }
      if( vars->currentPos == vars->npos ) {
        vars->currentPos = 0;
      }
      while( vars->currentPos < vars->npos && vars->keyValues[counterKey][vars->currentPos] != (double)vars->keyIntBuffer[counterKey] ) {
        vars->currentPos++;
      }
      if( vars->currentPos == vars->npos ) {
        if( vars->showWarnings ) {
          log->write("Unable to find matching line in ASCII file for ");
          for( int i = 0; i < vars->numKeys; i++ ) {
            log->write("key #%d: %d  ", i+1, vars->keyIntBuffer[i] );
          }
          log->write("\n");
        }
        return !vars->dropUnmatchedTraces;
      }
      else {
        for( int jKey = counterKey-1; jKey >= 0; jKey-- ) {
          if( vars->keyValues[jKey][vars->currentPos] != (double)vars->keyIntBuffer[jKey] ) {
            if( vars->showWarnings ) {
              log->write("Unable to find matching line in ASCII file for ");
              for( int i = 0; i < vars->numKeys; i++ ) {
                log->write("key #%d: %d  ", i+1, vars->keyIntBuffer[i] );
              }
              log->write("\n");
            }
            return !vars->dropUnmatchedTraces;
          }
        }
      }

    }
    // TYPE_STRING
    else {
      // Not implemented yet
    }
    counterKey += 1;
  }
//  if( vars->isTimeKey ) {
//    int currentTime_s = trcHdr->intValue( vars->hdrId_time_key );
//    if( currentTime_s
//  }

  for( int i = 0; i < vars->numHeaders; i++ ) {
    int type = vars->headerTypeList->at(i);
    int index = vars->headerIndexList->at(i);
    if( type == TYPE_INT ) {
      trcHdr->setIntValue( index, (int)vars->headerValues[i][vars->currentPos] );
    }
    else if( type == TYPE_INT64 ) {
      trcHdr->setInt64Value( index, (csInt64_t)vars->headerValues[i][vars->currentPos] );
    }
    else {
      trcHdr->setDoubleValue( index, vars->headerValues[i][vars->currentPos] );
    }
  }

  return true;
}

//*************************************************************************************************
// Parameter definition
//
//
//*************************************************************************************************
void params_mod_read_ascii_( csParamDef* pdef ) {
  pdef->setModule( "READ_ASCII", "Read trace header values from ASCII file" );

  pdef->setVersion( 1, 0 );

  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "method", "Method to specify header locations in ASCII file", NUM_VALUES_FIXED );
  pdef->addValue( "columns", VALTYPE_OPTION );
  pdef->addOption( "columns", "Location of key/header values is given by column number", "Columns are separated by white spaces" );
  pdef->addOption( "positions", "Location of key/header values is given by character positions from the start of each line", "First character in input line is at position 1" );

  pdef->addParam( "key", "Key trace header used to match specified trace header values", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name of key header" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number/Start position", "Depends on setting for user parameter 'method'" );
  pdef->addValue( "", VALTYPE_NUMBER, "Length", "Only used for method 'positions'" );

  pdef->addParam( "header", "Trace header to be read in", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number/Start position", "Depends on setting for user parameter 'method'" );
  pdef->addValue( "", VALTYPE_NUMBER, "Length", "Only used for method 'positions'" );

  pdef->addParam( "ignore_line", "Ignore lines starting with specified characters", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Ignore lines starting with specified characters" );

  pdef->addParam( "select_line", "Select lines starting with specified characters", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Select only lines starting with the specified characters" );

  pdef->addParam( "key_sps_time", "Use time (precision = 1s) as the key", NUM_VALUES_VARIABLE,
                  "Time in the ASCII file is expected in SPS format (dddhhmmss)" );
        //  pdef->addValue( "jjjhhmmss", VALTYPE_STRING, "Time format expected in input file. Valid identifiers: j(Julian day, starting at 1 for first day of year), h(hour), m(minute), s(second)", "Example: jjjhhmmss (SPS format)" );
        pdef->addValue( "", VALTYPE_STRING, "Trace header containing time in UNIX seconds [s] (for example time_samp1)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Column number/Start position", "Depends on setting of user parameter METHOD" );
  pdef->addValue( "", VALTYPE_NUMBER, "Length", "Only used for method 'positions'");

  pdef->addParam( "key_sps_time_us", "..to be used in conjunction with parameter 'key_sps_time'. Increase precision of match to microseconds",
                  NUM_VALUES_VARIABLE, "Time in the ASCII file is expected in EXTENDED SPS format (dddhhmmss.uuuuuu)" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header containing shot time microsecond [us] fraction (for example time_samp1_us)" );

  pdef->addParam( "time_year", "Year. Required in case a 'time' key is used.", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Year, e.g. 2009." );

  pdef->addParam( "warnings", "Show warnings?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Show warnings");
  pdef->addOption( "no", "Do not show warnings" );

  pdef->addParam( "check", "Check input file for consistency?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Check input file for consistency");
  pdef->addOption( "no", "Do not perform consistency check." );

  pdef->addParam( "sort", "Sort key values from input ASCII file?", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Sort key values in input ASCII file");
  pdef->addOption( "no", "Do not sort key values. WARNING: This assumes key values are already sorted on input, otherwise the method will not work" );
  
  pdef->addParam( "drop_traces", "Drop unmatched traces?", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "yes", "Drop traces for which no match could be found in input ASCII file");
  pdef->addOption( "no", "Do not drop unmatched traces" );
}

extern "C" void _params_mod_read_ascii_( csParamDef* pdef ) {
  params_mod_read_ascii_( pdef );
}
extern "C" void _init_mod_read_ascii_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* log ) {
  init_mod_read_ascii_( param, env, log );
}
extern "C" bool _exec_mod_read_ascii_( csTrace* trace, int* port, csExecPhaseEnv* env, csLogWriter* log ) {
  return exec_mod_read_ascii_( trace, port, env, log );
}


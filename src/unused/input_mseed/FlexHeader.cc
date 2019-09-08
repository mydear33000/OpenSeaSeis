/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "FlexHeader.h"

using namespace mseed;
using namespace std;

const int FlexHeader::MIN_BYTE_SIZE = 8;

FlexHeader::FlexHeader() {
  init();
  myByteSize    = MIN_BYTE_SIZE;
  myValue       = new char[myByteSize];
  myType        = TYPE_INT;
}
void FlexHeader::init() {  
  myValue    = NULL;
  myByteSize = 0;
}
void FlexHeader::init( int byteSize ) {  
  myByteSize    = byteSize;
  myValue       = new char[myByteSize];
}
FlexHeader::FlexHeader( FlexHeader const& obj ) {
  init( obj.myByteSize );
  myType     = obj.myType;
  memcpy( myValue, obj.myValue, myByteSize );
}
FlexHeader::FlexHeader( double d ) {
  init( MIN_BYTE_SIZE );
//  fprintf(stdout,"Constructor double %x %f\n", myValue, d);
  myType        = TYPE_DOUBLE;
  *((double*)myValue) = d;
}
FlexHeader::FlexHeader( float f ) {
  init( MIN_BYTE_SIZE );
//  fprintf(stdout,"Constructor float %x %f\n", myValue, f);
  myType        = TYPE_FLOAT;
  *((float*)myValue) = f;
}
FlexHeader::FlexHeader( int i ) {
  init( MIN_BYTE_SIZE );
//  fprintf(stdout,"Constructor int   %x %d\n", myValue, i );
  myType        = TYPE_INT;
  *((int*)myValue) = i;
}
FlexHeader::FlexHeader( char c ) {
  init( MIN_BYTE_SIZE );
//  fprintf(stdout,"Constructor int   %x %d\n", myValue, i );
  myType        = TYPE_CHAR;
  myValue[0] = c;
  myValue[1] = '\0';
}
FlexHeader::FlexHeader( int64_t i ) {
  init( MIN_BYTE_SIZE );
  myType        = TYPE_INT64;
  *((int64_t*)myValue) = i;
}
FlexHeader::FlexHeader( std::string& value ) {
//  fprintf(stdout,"Constructor int   %x %d\n", myValue, i );
  int newByteSize  = value.length()+1;
  init( max( newByteSize, MIN_BYTE_SIZE ) );
  myType = TYPE_STRING;

  memcpy( myValue, value.c_str(), newByteSize-1 );
  myValue[newByteSize-1] = '\0';
}
FlexHeader::~FlexHeader() {
//  fprintf(stdout,"Destructor %x\n", myValue);
  myType = TYPE_UNKNOWN;
  if( myValue ) {
    delete [] myValue;
    myValue = NULL;
  }
  myByteSize = 0;
}
//-------------------------------------------------------------
int FlexHeader::stringSize() {
  if( myType != TYPE_STRING ) {
    return 0;
  }
  else {
    return strlen(myValue);
  }
}

//-------------------------------------------------------------
FlexHeader& FlexHeader::operator=( const FlexHeader& obj ) {
//  fprintf(stdout,"Operator=            %x\n", myValue );
  if( obj.myType == TYPE_STRING && obj.myByteSize > myByteSize ) {
    freeMemory();
    myByteSize = obj.myByteSize;
    myValue    = new char[myByteSize];
  }
  myType = obj.myType;
  memcpy( myValue, obj.myValue, obj.myByteSize );
  return *this;
}
FlexHeader& FlexHeader::operator=( int i ) {
//  fprintf(stdout,"Operator=(int)       %x\n", myValue );
  setIntValue( i );
  return *this;
}
FlexHeader& FlexHeader::operator=( char c ) {
//  fprintf(stdout,"Operator=(int)       %x\n", myValue );
  setCharValue( c );
  return *this;
}
FlexHeader& FlexHeader::operator=( int64_t i ) {
//  fprintf(stdout,"Operator=(int)       %x\n", myValue );
  setInt64Value( i );
  return *this;
}
FlexHeader& FlexHeader::operator=( double d ) {
//  fprintf(stdout,"Operator=(double)    %x\n", myValue );
  setDoubleValue( d );
  return *this;
}
FlexHeader& FlexHeader::operator=( std::string& text ) {
  setStringValue( text );
  return *this;
}
//---------------------------------------------------------------------
//
bool FlexHeader::operator==( const FlexHeader& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    if( obj.myType != TYPE_STRING ) {
      return( doubleValue_internal() == obj.doubleValue() );
    }
    break;
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() == obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_INT || obj.myType == TYPE_FLOAT ) {
      return( floatValue_internal() == obj.floatValue() );
    }
    break;
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)intValue_internal() == obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (float)intValue_internal() == obj.floatValue_internal() );
    }
    else if( obj.myType == TYPE_INT ) {
      return( intValue_internal() == obj.intValue() );
    }
    break;
  case TYPE_STRING:
    if( obj.myType == TYPE_STRING ) {
      return( !strcmp( myValue, obj.myValue)  );
    }
    break;
  }
  return false;
}
//---------------------------------------------------------------------
bool FlexHeader::operator!=( const FlexHeader& obj ) const {
  switch( myType ) {
  case TYPE_DOUBLE:
    return( doubleValue_internal() != obj.doubleValue() );
  case TYPE_FLOAT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)floatValue_internal() != obj.doubleValue_internal() );
    }
    else {
      return( floatValue_internal() != obj.floatValue() );
    }
  case TYPE_INT:
    if( obj.myType == TYPE_DOUBLE ) {
      return( (double)intValue_internal() != obj.doubleValue_internal() );
    }
    else if( obj.myType == TYPE_FLOAT ) {
      return( (float)intValue_internal() != obj.floatValue_internal() );
    }
    else {
      return( intValue_internal() != obj.intValue() );
    }
  }
  return false;
}
//---------------------------------------------------------------------
//
void FlexHeader::dump() const {

  fprintf(stdout," Type: %2d: ", myType );
  if( myType == TYPE_INT ) {
    fprintf(stdout," valueInt:     %d\n", *((int*)myValue) );
  }
  else if( myType == TYPE_FLOAT ) {
    fprintf(stdout," valueFloat:   %f\n", *((float*)myValue) );
  }
  else if( myType == TYPE_DOUBLE ) {
    fprintf(stdout," valueDouble:  %f\n", *((double*)myValue) );
  }
  else if( myType == TYPE_INT64 ) {
    fprintf(stdout," valueInt64:   %lld\n", *((int64_t*)myValue) );
  }
  else if( myType == TYPE_CHAR ) {
    fprintf(stdout," valueChar:   '%c'\n", myValue[0] );
  }
  else { // if( myType == TYPE_STRING ) {
    fprintf(stdout," valueString: '%s'\n", myValue);
  }
}



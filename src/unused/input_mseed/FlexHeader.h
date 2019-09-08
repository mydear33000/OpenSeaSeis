/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef FLEX_HEADER_H
#define FLEX_HEADER_H

#include <algorithm>
#include <string>

namespace mseed {

  /**
   * Flexible header class
   * This class stores header value in a memory efficient way.
   * Supported formats are
   *  int    (32bit)
   *  float  (32bit)
   *  double (64bit)
   *  int64  (64bit)
   *  string (variable size)
   *
   * @author Bjorn Olofsson
   * @date   2007
   */

class FlexHeader {
 public:
  typedef unsigned char type_t;
  typedef signed long long int64_t;

  /// Trace and parameter types
  static const type_t TYPE_UNKNOWN = 255;
  static const type_t TYPE_EMPTY   = 0;
  static const type_t TYPE_INT     = 1;  // 32bit
  static const type_t TYPE_FLOAT   = 2;  // 32bit
  static const type_t TYPE_DOUBLE  = 3;  // 64bit
  static const type_t TYPE_CHAR    = 4;  // 8bit
  static const type_t TYPE_STRING       = 5;  // = ARRAY_CHAR
  static const type_t TYPE_INT64        = 6;  // 64bit
  static const type_t TYPE_ARRAY_INT    = 7;
  static const type_t TYPE_ARRAY_FLOAT  = 8;
  static const type_t TYPE_ARRAY_DOUBLE = 9;
  static const type_t TYPE_OPTION = 10;
  static const type_t TYPE_FLOAT_DOUBLE = 23;

public:
  FlexHeader();
  FlexHeader( double d );
  FlexHeader( float f );
  FlexHeader( int i );
  FlexHeader( char c );
  FlexHeader( int64_t i );
  FlexHeader( std::string& text );
  FlexHeader( FlexHeader const& obj );
  void init();
  void init( int byteSize );
  ~FlexHeader();
  FlexHeader& operator=( const FlexHeader& obj );
  FlexHeader& operator=( double d );
  FlexHeader& operator=( int i );
  FlexHeader& operator=( char c );
  FlexHeader& operator=( int64_t i );
  FlexHeader& operator=( std::string& text );
  /// @return type of number
  inline type_t type() const { return myType; }
  inline void setStringValue( std::string& value ) {
    int newByteSize = std::max( MIN_BYTE_SIZE, (int)value.length()+1 );
    if( myByteSize < newByteSize ) {
      freeMemory();
      myValue = new char[newByteSize];
    }
    myByteSize = newByteSize;
    memcpy( myValue, value.c_str(), value.length() );
    myValue[value.length()] = '\0';
    myType = TYPE_STRING;
  }
  inline void setFloatValue( float value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    *((float*)myValue) = value;
    myType = TYPE_FLOAT;
  }
  inline void setDoubleValue( double value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    *((double*)myValue) = value;
    myType = TYPE_DOUBLE;
  }
  inline void setIntValue( int value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    *((int*)myValue) = value;
    myType = TYPE_INT;
  }
  inline void setCharValue( char value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    myValue[0] = value;
    myValue[1] = '\0';
    myType = TYPE_CHAR;
  }
  inline void setInt64Value( int64_t value ) {
    if( myType == TYPE_STRING ) {
      freeMemory();
      myValue = new char[MIN_BYTE_SIZE];
      myByteSize = MIN_BYTE_SIZE;
    }
    *((int64_t*)myValue) = value;
    myType = TYPE_INT64;
  }

  inline char charValue() const {
    switch( myType ) {
    case TYPE_CHAR:
    case TYPE_STRING:
      return myValue[0];
    default:
      return ' ';
    }
  }
  inline int intValue() const {
    switch( myType ) {
    case TYPE_INT:
      return *((int*)myValue);
    case TYPE_FLOAT:
      return (int)( *((float*)myValue) );
    case TYPE_DOUBLE:
      return (int)( *((double*)myValue) );
    case TYPE_INT64:
      return (int)( *((int64_t*)myValue) );
    default:
      return 0;
    }
  }
  inline int64_t int64Value() const {
    switch( myType ) {
    case TYPE_INT64:
      return *((int64_t*)myValue);
    case TYPE_INT:
      return (int64_t)( *((int*)myValue) );
    case TYPE_FLOAT:
      return (int64_t)( *((float*)myValue) );
    case TYPE_DOUBLE:
      return (int64_t)( *((double*)myValue) );
    default:
      return 0;
    }
  }
  inline float floatValue() const {
    switch( myType ) {
    case TYPE_FLOAT:
      return *((float*)myValue);
    case TYPE_DOUBLE:
      return (float)( *((double*)myValue) );
    case TYPE_INT:
      return (float)( *((int*)myValue) );
    case TYPE_INT64:
      return (float)( *((int64_t*)myValue) );
    default:
      return 0.0;
    }
  }
  inline double doubleValue() const {
    switch( myType ) {
    case TYPE_DOUBLE:
      return *((double*)myValue);
    case TYPE_FLOAT:
      return (double)( *((float*)myValue) );
    case TYPE_INT:
      return (double)( *((int*)myValue) );
    case TYPE_INT64:
      return (double)( *((int64_t*)myValue) );
    default:
      return 0.0;
    }
  }
  inline char const* stringValue() const {
    switch( myType ) {
    case TYPE_STRING:
    case TYPE_CHAR:
      return myValue;
    default:
      return " ";
    }
  }
  int stringSize();
  
  inline bool isUndefined() const { return( myType == TYPE_UNKNOWN ); }

  bool operator==( const FlexHeader& obj ) const;
  bool operator!=( const FlexHeader& obj ) const;

  void dump() const;
private:
  static const int MIN_BYTE_SIZE;
  char*    myValue;
  type_t   myType;
  int      myByteSize;

// Internal methods, without any checks
  inline char charValue_internal() const {
    return myValue[0];
  }
  inline int intValue_internal() const {
    return *((int*)myValue);
  }
  inline int64_t int64Value_internal() const {
    return *((int64_t*)myValue);
  }
  inline float floatValue_internal() const {
    return *((float*)myValue);
  }
  inline double doubleValue_internal() const {
    return *((double*)myValue);
  }
  inline char const* stringValue_internal() const {
    return myValue;
  }
  inline void freeMemory() {
    if( myValue != NULL ) {
      delete [] myValue;
      myValue = NULL;
    }
  }

};

}

#endif


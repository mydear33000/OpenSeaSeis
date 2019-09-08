/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_POINT_H
#define CS_POINT_H

#include <cstdlib>
#include <cstdio>

namespace cseis_geolib {

/**
 * Point representation
 *
 *
 */
class csPoint {
 public:
 csPoint() : x(0), y(0), z(0) {}
 csPoint( csPoint const& obj ) : x(obj.x), y(obj.y), z(obj.z) {}
 csPoint( double xin, double yin, double zin ) : x(xin), y(yin), z(zin) {}
  ~csPoint() {}
  csPoint& operator=( csPoint const& obj ) {
    x = obj.x;
    y = obj.y;
    z = obj.z;
    return *this;
  }
  void dump( std::FILE* stream ) const {
    if( stream == NULL ) stream = stdin;
    fprintf(stream,"%f %f %f\n", x, y, z);
  }
 public:
  double x;
  double y;
  double z;
};

} // end namespace

#endif


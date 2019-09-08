/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_LINE_3D_H
#define CS_LINE_3D_H

#include <cstdlib>
#include "csPoint.h"

namespace cseis_geolib {

/**
 * 3D Line representation.
 * Line is defined by two end points
 */
class csLine3D {
 public:
  csLine3D();
  csLine3D( csPoint const& p1, csPoint const& p2 );
  virtual ~csLine3D();
  virtual csPoint pointAt( double lamda ) const;

  void set( csPoint const& p1, csPoint const& p2 );

  csPoint p1() const { return myP1; }
  csPoint p2() const { return myP2; }
  csPoint delta() const { return myDelta; }
  /// Offset between end points (= horizontal distance)
  double offset() const { return myOffset; }
  /// Distance between end points (=3D distance)
  double distance() const { return myDistance; }
  /// Azimuth of line
  double azimuth() const { return myAzimuth; }

 protected:
  void compute_internal();

  csPoint myP1;
  csPoint myP2;
  csPoint myDelta;
  double myOffset;
  double myDistance;
  double myAzimuth;
};


} // end namespace

#endif


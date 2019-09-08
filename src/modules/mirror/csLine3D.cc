/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csLine3D.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>

using namespace std;
using namespace cseis_geolib;

csLine3D::csLine3D() : myP1(), myP2(), myDelta(), myOffset(0), myDistance(0), myAzimuth(0) {
}

csLine3D::csLine3D( csPoint const& p1, csPoint const& p2 ) : myP1(p1), myP2(p2), myDelta(), myOffset(0), myDistance(0), myAzimuth(0) {
  compute_internal();
}
void csLine3D::set( csPoint const& p1, csPoint const& p2 ) {
  myP1 = p1;
  myP2 = p2;
  compute_internal();
}
void csLine3D::compute_internal() {
  double dx = myP2.x - myP1.x;
  double dy = myP2.y - myP1.y;
  double dz = myP2.z - myP1.z;
  myAzimuth = atan2( dx, dy );
  myOffset  = sqrt( dx*dx + dy*dy );
  myDistance= sqrt( dx*dx + dy*dy + dz*dz );
  myDelta.x = dx;
  myDelta.y = dy;
  myDelta.z = dz;
}
//--------------------------------------------------
csLine3D::~csLine3D() {
}
//--------------------------------------------------
csPoint csLine3D::pointAt( double lamda ) const {
  return csPoint( myP1.x + lamda*myDelta.x, myP1.y + lamda*myDelta.y, myP1.z + lamda*myDelta.z );
}


/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csMatrixFStyle.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace cseis_geolib;
using namespace std;

int wfront_check_params2(
                         FILE* f_out,
                         csMatrixFStyle<float> const& x_int,
                         int* npoints_int,
                         int nreclines,
                         int* nrec,
                         float x_source, float z_source,
                         csMatrixFStyle<float> const& xz_rec,
                         int flag_stack )
{

  //------
  // Check source
  if( flag_stack != 0 ) {
    if( x_source != xz_rec(0,0) ) {
      x_source = xz_rec(0,0);
      fprintf(f_out,"Changed source X coordinate to %f\n", x_source);
    }
    z_source = 0.0;
  }

  //-------------
  // Check receivers

  for( int i = 0; i < nreclines; i++ ) {
    if( xz_rec(0,i) < x_int(0,0) ) {
      fprintf(stderr,"Receiver line %d, receiver 1 is outside the model\n", i+1);
      return -1;
    }
    else if( xz_rec(nrec[i]-1,i) > x_int(npoints_int[0]-1,0) ) {
      fprintf(stderr,"Receiver line %d: Receivers are located outside the model\n", i+1);
      return -1;
    }
  }
  return 0;
}

int wfront_check_params(
          int n_int,
          csMatrixFStyle<float> const& x_int,
          int* npoints_int,
          int* index_grid,
          csMatrixFStyle<float> const& x_grid,
          csMatrixFStyle<float> const& z_grid,
          int* nx_grid,
          int* nz_grid,
          int lay_source, int int_source,
          float anglein, float angleout,
          int nray_source,
          int nreclines,
          int* lay_rec,
          int* int_rec,
          int f_wfronts,
          FILE* f_out,
          int step_wfronts,
          int ncodes,
          csMatrixFStyle<int> const& code,
          int* codesteps,
          int flag_2point,
          int flag_stack )
{
  // Check interfaces, model
  //
  for( int i = 0; i < n_int; i++ ) {
    if( x_int(0,i) != x_int(0,0) ) {
      fprintf(stderr,"X-coord. of interface no.%d= %f or int no.1= %f, not at model border\n", i, x_int(0,i), x_int(0,0) );
      return -1;
    }
    else if( x_int(npoints_int[i]-1,i) != x_int(npoints_int[0]-1,0) ) {
      fprintf(stderr,"X-coord. of interface no.%d= %f or int no.1= %f, not at model border\n",
              i, x_int(npoints_int[i]-1,i), x_int(npoints_int[0]-1,0));
      return -1;
    }
  }
  //--------------
  // Check velocity distribution and ptos
  
  for( int i = 0; i < n_int-1; i++ ) {
    if( index_grid[i] == 0 ) {
      for( int j = 1; j < nx_grid[i]; j++ ) {
        if( x_grid(j,i) <= x_grid(j-1,i) ) {
          fprintf(stderr," Velocity grid points of layer %d not in ascending order." \
                  " X-coord. %d: %f , X-coord. %d: %f\n", i, j-1, x_grid(j-1,i), j, x_grid(j,i) );
          return -1;
        }
      }
      for( int j = 1; j < nz_grid[i]; j++ ) {
        if( z_grid(j,i) <= z_grid(j-1,i) ) {
          fprintf(stderr," Velocity grid points of layer %d not in ascending order." \
                  " Z-coord. %d: %f , Z-coord. %d: %f\n", i,j-1,z_grid(j-1,i),j,z_grid(j,i) );
          return -1;
        }
      }
    }
  }
  //------
  // Check source
  if( flag_stack == 0 ) {
    if( lay_source < 1 || lay_source > (n_int-1) ) {
      fprintf(stderr,"Incorrect layer specified for source: %d\n", lay_source);
      return -1;
    }
    else if( int_source > 1 ) {
      fprintf(stderr,"Source must be on first interface or within a layer.\n");
      fprintf(stderr,"Input setting:\n");
      fprintf(stderr,"Source   layer / interface =  %d / %d\n", lay_source, int_source);
      fprintf(stderr,"Receiver layer / interface =  %d / %d\n", lay_rec[0], int_rec[0]);
      return -1;
    }
    else if( fabs(angleout - anglein) > 2*M_PI ) {
      fprintf(stderr,"The beam is too big (more than 360 degrees)\n");
      return -1;
    }
    else if( nray_source < 2 ) {
      fprintf(stderr, "At least two rays must start from the source...\n");
      return -1;
    }
  }
  else {
    if( lay_source != lay_rec[0] ) {
      lay_source = lay_rec[0];
      fprintf(f_out, "Changed source layer to %d\n", lay_source);
    }
    if( int_source != int_rec[0] ) {
      int_source = int_rec[0];
      fprintf(f_out, "Changed source interface to %d\n", int_source);
    }
  }

  //-------------
  // Check receivers

  for( int i = 0; i < nreclines; i++ ) {
    if( lay_rec[i] != int_rec[i] && lay_rec[i] != (int_rec[i] - 1) ) {
      fprintf(stderr,"Receiver line %d: Receivers are not located in a layer beneath or above the interface (%d). Set layer number (%d) correctly\n", i+1, int_rec[i], lay_rec[i]);
      return -1;
    }
  }
  //-------------
  // Check codes
  //

  if( int_source == 0 ) {
    for( int i = 0; i < ncodes; i++ ) {
      if( code(0,i) == 0 ) {
        fprintf(stderr,"Source inside layer: Specify if the ray goes up (-1) or down (1) in code no. %d\n",i+1);
        return -1;
      }
      else if( fabs(code(1,i)) != lay_source ) {
        fprintf(stderr," Code no. %d step 1 is wrong: %d. Should be +/-%d (source layer index)\n", i+1, code(1,i), lay_source);
        return -1;
      }
    }
  }
  else {
    int lay1 = lay_source;
    int lay2 = 0;
    if( int_source == lay_source ) {
      lay2 = lay_source-1;
    }
    else {
      lay2 = lay_source+1;
    }
    for( int i = 0; i < ncodes; i++ ) {
      if( abs(code(1,i)) != lay1 && abs(code(1,i)) != lay2 ) {
        fprintf(stderr," Code no. %d step 1 is wrong.\n", i+1);
        fprintf(stderr," Code no. %d step 1 is wrong: %d. Should be +/-%d or +/-%d (layer index above or below source interface)\n", i+1, code(1,i), lay1, lay2);
        return -1;
      }
    }
  }
  /*
  for( int i = 0; i < ncodes; i++ ) {
    for( int j = 1; j <= codesteps[i]; j++ ) {
      int codeCurrent = abs(code(j,i));
      fprintf(stdout,"Code #%d, step #%d = %d\n", i+1, j, codeCurrent);
    }
  }
  */
  int updownflag = 0;
  for( int i = 0; i < ncodes; i++ ) {
    int codePrev = abs(code(1,i));
    if( codePrev == 0 ) {
      fprintf(stderr,"Code #%d step #%d is zero\n", i+1, 1);
      return -1;
    }
    else if( codePrev >= n_int ) {
      fprintf(stderr," Code #%d step #%d (=%d) exceeds number of layers (%d)\n", i+1, 1, codePrev, n_int-1 );
      return -1;
    }
    if( int_source == 0 ) {
      updownflag = code(0,i);
    }
    else {
      if( codePrev == int_source ) {
        updownflag = 1;
      }
      else {
        updownflag = -1;
      }
    }
    for( int j = 2; j <= codesteps[i]; j++ ) {
      int codeCurrent = abs( code(j,i) );
      //      fprintf(stdout,"Code  #%d #%d = %d   (prev = %d) num = %d\n", i+1, j, code(j,i), codePrev, codesteps[i] );
      if( codeCurrent == 0 ) {
        fprintf(stderr,"Code #%d step #%d is zero\n", i+1, j+1);
        return -1;
      }
      else if( codeCurrent >= n_int ) {
        fprintf(stderr," Code #%d step #%d (=%d) exceeds number of layers (=%d)\n", i+1, j, codeCurrent, n_int-1 );
        return -1;
      }
      else if( codePrev == (n_int-1) ) {
        if( codeCurrent == (n_int-1) && updownflag > 0 ) {
          fprintf(stderr," Code #%d step #%d (=%d): Reflection at last interface not possible\n", i+1, j, codeCurrent );
          return -1;
        }
      }
      if( updownflag < 0 ) {
        if( codeCurrent != codePrev && codeCurrent != (codePrev-1) ) {
          fprintf(stderr," Code #%d step #%d (=%d): Expected (+/-)%d or (+/-)%d \n", i+1, j, codeCurrent, codePrev, (codePrev-1) );
          return -1;
        }
      }
      else if( codeCurrent != codePrev && codeCurrent != (codePrev+1) ) {
          fprintf(stderr," Code #%d step #%d (=%d): Expected (+/-)%d or (+/-)%d \n", i+1, j, codeCurrent, codePrev, (codePrev+1) );
          return -1;
      }
      if( codePrev == codeCurrent ) {
        updownflag = -updownflag;
      }
      codePrev = codeCurrent;
    } // END for j
  } // END for i

  //-------------
  // Check stacking specialities

  if( flag_stack != 0 ) {
    if( f_wfronts != 0 ) {
      fprintf(stderr,"WARNING: You specified to plot out wavefronts during a STACK calculation\n");
      return -1;
    }
  }

  //-------------
  // Check two-point ray-tracing specialities

  if( flag_2point != 0 ) {
    if( nreclines > 1 ) {
      fprintf(stderr,"For two-point ray-tracing only a single receiver line is supported\n");
      return -1;
    }
  }
  
  return 0;
}


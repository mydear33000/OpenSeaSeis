c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm GET_MIN
c 'Berechnet' Minimalwert (z-Komponente) eines kubischen Splines
c
      subroutine get_min(x_int,z_int,d_int,c_int,b_int,npoints_int,
     &     zmin,
     &     MAXPOINTS_INT)

      implicit none

      integer MAXPOINTS_INT
      integer npoints_int
      real x_int(MAXPOINTS_INT), z_int(MAXPOINTS_INT)
      real d_int(MAXPOINTS_INT), c_int(MAXPOINTS_INT), b_int(MAXPOINTS_INT)
      real zmin

      real xstep,x,x0,z
      integer i, actp, nsteps

      nsteps = MAXPOINTS_INT*10
      
      xstep = (x_int(npoints_int) - x_int(1)) / float(nsteps-1)
      zmin = 1000000.
      x0 = x_int(1)

      do i = 1, nsteps
         x = x0 + xstep*float(i-1)
         call get_z(x,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &        z,actp,
     &        MAXPOINTS_INT)
         if (z .lt. zmin) zmin = z
      end do

      do i = 1, npoints_int
         if (z_int(i) .lt. zmin) zmin = z_int(i)
      end do

      end
c********************************************************************
c Unterprogramm GET_MAX
c 'Berechnet' Minimalwert (z-Komponente) eines kubischen Splines
c
      subroutine get_max(x_int,z_int,d_int,c_int,b_int,npoints_int,
     &     zmax,
     &     MAXPOINTS_INT)

      implicit none

      integer MAXPOINTS_INT
      integer npoints_int
      real x_int(MAXPOINTS_INT), z_int(MAXPOINTS_INT)
      real d_int(MAXPOINTS_INT), c_int(MAXPOINTS_INT), b_int(MAXPOINTS_INT)
      real zmax

      real xstep,x,x0,z
      integer i, nsteps, actp

      nsteps = MAXPOINTS_INT*10
      
      xstep = (x_int(npoints_int) - x_int(1)) / float(nsteps-1)
      zmax = -1000000.
      x0 = x_int(1)

      do i = 1, nsteps
         x = x0 + xstep*float(i-1)
         call get_z(x,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &        z,actp,
     &        MAXPOINTS_INT)
         if (z .gt. zmax) zmax = z
      end do

      do i = 1, npoints_int
         if (z_int(i) .gt. zmax) zmax = z_int(i)
      end do

      end




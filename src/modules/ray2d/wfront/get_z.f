c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm GET_Z / Hauptprogramm WFRONT
c Berechnet die z-Komponente zur entsprechenden x-Komponente einer Grenzflaeche
c XIN  : x-Komponente des Punktes auf Grenzflaeche
c ZOUT : z-Komponente
c ACTP : Index des Stuetzpunktes auf der Grenzflaeche, der zur Berechnung von ZOUT
c        benutzt wird
c
      subroutine get_z(xin,x_int,z_int,d_int,c_int,b_int,npoints,
     &                zout,actp,
     &                MAXPOINTS_INT)

      implicit none

      integer MAXPOINTS_INT
      integer npoints
      real xin, x_int(MAXPOINTS_INT), z_int(MAXPOINTS_INT)
      real d_int(MAXPOINTS_INT), c_int(MAXPOINTS_INT), b_int(MAXPOINTS_INT)
      real zout
      integer actp
      
      real dx
      integer i

      do 1 i = 1, npoints - 1
         actp = i
         if (xin .lt. x_int(i+1)) goto 2
 1    continue
 2    dx = xin - x_int(actp)
      zout = ((d_int(actp)*dx + c_int(actp))*dx + b_int(actp))*dx + z_int(actp)

      end

c****************************************************************************************************
c Unterprogramm GET_Z1 / Hauptprogramm WFRONT
c Berechnet die z-Komponente zur entsprechenden x-Komponente einer Grenzflaeche sowie
c die Steigung der Grenzflaeche dz/dx.
c XIN  : x-Komponente des Punktes auf Grenzflaeche
c ZOUT : z-Komponente
c DZDX : Steigung
c ACTP : Index des Stuetzpunktes auf der Grenzflaeche, der zur Berechnung von Z und DZDX
c        benutzt wird
c
      subroutine get_z1(xin,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &     zout,dzdx,actp,
     &     MAXPOINTS_INT)

      implicit none

      integer MAXPOINTS_INT
      integer npoints_int
      real d_int(MAXPOINTS_INT),c_int(MAXPOINTS_INT)
      real b_int(MAXPOINTS_INT),z_int(MAXPOINTS_INT)
      real x_int(MAXPOINTS_INT)
      integer actp
      real xin
      real zout, dzdx

      real dx
      integer actp0,i
      
c      if (actp .gt. npoints_int) actp = 1

      if (xin .lt. x_int(actp)) then
         do i = actp-1,1,-1
            actp0 = i
            if (xin .gt. x_int(actp0)) goto 1
         end do
      else if (actp .ne. npoints_int .and. xin .gt. x_int(actp+1)) then
         do i = actp+1, npoints_int-1
            actp0 = i
            if (xin .lt. x_int(actp0+1)) goto 1
         end do
      else
         actp0 = actp
      end if


 1    dx = xin - x_int(actp0)
      zout = ((d_int(actp0)*dx + c_int(actp0))*dx + b_int(actp0))*dx + z_int(actp0)
      dzdx = (3.*d_int(actp0)*dx + 2.*c_int(actp0))*dx + b_int(actp0)

      actp = actp0

      end
      
c****************************************************************************************************
c Unterprogramm GET_Z2 / Hauptprogramm WFRONT
c Berechnet die z-Komponente zur entsprechenden x-Komponente einer Grenzflaeche,
c die Steigung der Grenzflaeche dz/dx sowie die zweite Ableitung d2z/dx2.
c X      : x-Komponente des Punktes auf Grenzflaeche
c Z      : z-Komponente
c DZDX   : Steigung
c D2ZDX2 : zweite Ableitung
c ACTP   : Index des Stuetzpunktes auf der Grenzflaeche, der zur Berechnung von Z und DZDX
c          benutzt wird
c
      subroutine get_z2(x,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &     z,dzdx,d2zdx2,actp,
     &     MAXPOINTS_INT)

      implicit none

      integer MAXPOINTS_INT
      integer npoints_int
      real d_int(MAXPOINTS_INT),c_int(MAXPOINTS_INT)
      real b_int(MAXPOINTS_INT),z_int(MAXPOINTS_INT)
      real x_int(MAXPOINTS_INT)
      integer actp
      real x
      real z, dzdx, d2zdx2

      real dx
      integer actp0,i

c Sicherheit
      if (actp .lt. 1) then
         actp = 1
      else if (actp .ge. npoints_int) then
         actp = npoints_int-1
      end if

      if (x .lt. x_int(actp)) then
         do i = actp-1,1,-1
            actp0 = i
            if (x .gt. x_int(actp0)) goto 1
         end do
      else if (actp .ne. npoints_int .and. x .gt. x_int(actp+1)) then
         do i = actp+1,npoints_int-1
            actp0 = i
            if (x .lt. x_int(actp0+1)) goto 1
         end do
      else
         actp0 = actp
      end if

 1    dx = x - x_int(actp0)
      z = ((d_int(actp0)*dx + c_int(actp0))*dx + b_int(actp0))*dx + z_int(actp0)
      dzdx = (3.*d_int(actp0)*dx + 2.*c_int(actp0))*dx + b_int(actp0)
      d2zdx2 = 6.*d_int(actp0)*dx + 2.*c_int(actp0)

      actp = actp0

      end
      
c****************************************************************************************************
c Unterprogramm GET_ACTP
c Bestimmt den Index ACTP des Grenzflaechenstueckes, auf dem der Punkt (XIN,) liegt
c
      subroutine get_actp(xin,interface,x_int,npoints,
     &                actp,
     &                MAX_INT,MAXPOINTS_INT)

      implicit none

      integer MAX_INT, MAXPOINTS_INT
      integer interface, npoints(MAXPOINTS_INT)
      real xin, x_int(MAXPOINTS_INT,MAX_INT)

      integer actp
      
      integer i

      do 1 i = 1, npoints(interface) - 1
         actp = i
         if (xin .lt. x_int(i+1,interface)) goto 2
 1    continue

 2    end

c****************************************************************************************************
c Unterprogramm GET_Z1_AT_ACTP / Hauptprogramm WFRONT
c Berechnet die z-Komponente zur entsprechenden x-Komponente einer Grenzflaeche
c Index ACTP ist zwingend vorgeschrieben
c XIN  : x-Komponente des Punktes auf Grenzflaeche
c ZOUT : z-Komponente
c ACTP : Index des Stuetzpunktes auf der Grenzflaeche, der zur Berechnung von ZOUT
c        benutzt wird
c
      subroutine get_z1_at_actp(xin,x_int,z_int,d_int,c_int,b_int,
     &                zout,dzdxout)

      implicit none

      real xin, x_int, z_int
      real d_int, c_int, b_int
      real zout,dzdxout
      
      real dx

      dx = xin - x_int
      zout = ((d_int*dx + c_int)*dx + b_int)*dx + z_int
      dzdxout = (3.*d_int*dx + 2.*c_int)*dx + b_int

      end

c****************************************************************************************************
c Unterprogramm GET_DZDX / Hauptprogramm WFRONT
c Berechnet die z-Komponente zur entsprechenden x-Komponente einer Grenzflaeche sowie
c die Steigung der Grenzflaeche dz/dx.
c X    : x-Komponente des Punktes auf Grenzflaeche
c DZDX : Steigung
c ACTP : Index des Stuetzpunktes auf der Grenzflaeche, der zur Berechnung von Z und DZDX
c        benutzt wird
c
      subroutine get_dzdx(x,x_int,d_int,c_int,b_int,npoints_int,
     &     dzdx,actp,
     &     MAXPOINTS_INT)

      implicit none

      integer MAXPOINTS_INT
      integer npoints_int
      real d_int(MAXPOINTS_INT),c_int(MAXPOINTS_INT)
      real b_int(MAXPOINTS_INT)
      real x_int(MAXPOINTS_INT)
      integer actp
      real x
      real dzdx

      real dx
      integer actp0,i
      
c      if (actp .gt. npoints_int) actp = 1

      if (x .lt. x_int(actp)) then
         do i = actp-1,1,-1
            actp0 = i
            if (x .gt. x_int(actp0)) goto 1
         end do
      else if (actp .ne. npoints_int .and. x .gt. x_int(actp+1)) then
         do i = actp+1, npoints_int-1
            actp0 = i
            if (x .lt. x_int(actp0+1)) goto 1
         end do
      else
         actp0 = actp
      end if


 1    dx = x - x_int(actp0)
      dzdx = (3.*d_int(actp0)*dx + 2.*c_int(actp0))*dx + b_int(actp0)

      actp = actp0

      end
      


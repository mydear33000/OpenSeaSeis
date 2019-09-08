c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm VELOCITY / Hauptprogramm WFRONT
c Dieses Programm stammt vom Programmpaket SEIS88, Originalname VELOC aus MODSP.F
c Determination of velocity and derivatives from bicubic polynomial approximation
c Ausgabe:
c vel(1) = v
c vel(2) = dv/dx
c vel(3) = dv/dz
c vel(4) = d2v/dx2
c vel(5) = d2v/dxdz
c vel(6) = d2v/dz2
c
C*************************************************************
C
      subroutine velocity(x0, z0, nx_grid, nz_grid, x_grid, z_grid, v_grid,
     &     veltype, vpvs_ratio, actlayer,
     &     vel,
     &     MAXP_XGRID, MAXP_ZGRID)

      implicit none

      integer MAXP_XGRID, MAXP_ZGRID
      real x0, z0, vpvs_ratio
      integer nx_grid, nz_grid, veltype, actlayer
      real x_grid(MAXP_XGRID), z_grid(MAXP_ZGRID), v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real vel(6)

      integer i, ix1,iz1, layer
      real xlast, zlast
      real xa, xb, hx, hz, za, zb, ax, az
      real b00,b02,b20,b22,c00,c02,c20,c22,d00,d02,d20,d22,e00,e02,e20,e22
      real aux1,aux2,aux3,aux4
      
      integer ix2,iz2
      real b01,b03,b10,b11,b12,b13,b21,b23,b30,b31,b32,b33
      real d10,d12,d30,d32

      common layer,b01,b03,b10,b11,b12,b13,b21,b23,b30,b31,b32,b33,b00,b02,b20,b22,ix2,iz2,xlast,zlast

C     DETERMINATION OF VELOCITY AND ITS DERIVATIVES
C     FOR BICUBIC POLYNOMIAL APPROXIMATION

      do i = 2, nx_grid
         ix1 = i
         if (x0 .lt. x_grid(ix1)) goto 2
      end do
    2 continue
c Punkt ist ausserhalb des velocity network (z-Wert zu klein)
c wird aber trotzdem berechnet (Extrapolation des velocity grids)
      if (z0 .lt. z_grid(1)) then
         iz1 = 2
         goto 4
      end if
      do i = 2, nz_grid
         iz1 = i
         if (z0 .lt. z_grid(iz1)) goto 4
      end do
    4 continue

c-----------------
c Bestimmung der Parameter (nur, wenn nicht letztes Mal dieselben gebraucht wurden)
c
      if (layer .ne. actlayer .or. ix1 .ne. ix2 .or. iz1 .ne. iz2) then
         ix2 = ix1
         iz2 = iz1
         ix1 = ix1 - 1
         iz1 = iz1 - 1

         if(ix2-1 .lt. 1) goto 999
         xlast = x_grid(ix2-1)
         zlast = z_grid(iz2-1)

         hx = x_grid(ix2) - xlast
         hz = z_grid(iz2) - zlast
         XA = 3.*HX
         ZA = 3.*HZ
         XB = HX/3.
         ZB = HZ/3.

         B20 = v_grid(2,iz1,ix1)
         B02 = v_grid(3,iz1,ix1)
         B22 = v_grid(4,iz1,ix1)
         B00 = v_grid(1,iz1,ix1)

         C20 = v_grid(2,iz1,ix1+1)
         C02 = v_grid(3,iz1,ix1+1)
         C22 = v_grid(4,iz1,ix1+1)
         C00 = v_grid(1,iz1,ix1+1)

         d20 = v_grid(2,iz1+1,ix1)
         d02 = v_grid(3,iz1+1,ix1)
         d22 = v_grid(4,iz1+1,ix1)
         d00 = v_grid(1,iz1+1,ix1)

         e20 = v_grid(2,iz1+1,ix1+1)
         e02 = v_grid(3,iz1+1,ix1+1)
         e22 = v_grid(4,iz1+1,ix1+1)
         e00 = v_grid(1,iz1+1,ix1+1)

         D32 = (e22-d22) / XA
         D30 = (e20-d20) / XA
         B30 = (c20-b20) / XA
         B32 = (c22-b22) / XA
         D12 = (e02-d02) / HX - XB * (e22 + 2.*d22)
         D10 = (e00-d00) / HX - XB * (e20 + 2.*d20)
         B10 = (c00-b00) / HX - XB * (C20 + 2.*B20)
         B12 = (c02-b02) / HX - XB * (C22 + 2.*B22)
         B03 = (d02-b02) / ZA
         B13 =(D12-B12) / ZA
         B23 = (d22-b22) / ZA
         B33 = (D32-B32) / ZA
         B01 = (D00-B00) / HZ - ZB * (D02 + 2.*B02)
         B11 = (D10-B10) / HZ - ZB * (D12 + 2.*B12)
         B21 = (D20-B20) / HZ - ZB * (D22 + 2.*B22)
         B31 = (D30-B30) / HZ - ZB * (D32 + 2.*B32)

         layer = actlayer
      end if
      ax = x0 - xlast
      az = z0 - zlast

      AUX1 = ((B33*AZ + B32)*AZ + B31)*AZ + B30
      AUX2 = ((B23*AZ + B22)*AZ + B21)*AZ + B20
      AUX3 = ((B13*AZ + B12)*AZ + B11)*AZ + B10
      AUX4 = ((B03*AZ + B02)*AZ + B01)*AZ + B00
      vel(1) = ((AUX1*AX + AUX2)*AX + AUX3)*AX + AUX4
      vel(2) = (3.*AUX1*AX + 2.*AUX2)*AX + AUX3
      vel(4) = 6.*AUX1*AX + 2.*AUX2

      AUX1 = (3.*B33*AZ + 2.*B32)*AZ + B31
      AUX2 = (3.*B23*AZ + 2.*B22)*AZ + B21
      AUX3 = (3.*B13*AZ + 2.*B12)*AZ + B11
      AUX4 = (3.*B03*AZ + 2.*B02)*AZ + B01
      vel(3) = ((AUX1*AX + AUX2)*AX + AUX3)*AX + AUX4
      vel(5) = (3.*AUX1*AX + 2.*AUX2)*AX + AUX3

      AUX1 = 3.*B33*AZ + B32
      AUX2 = 3.*B23*AZ + B22
      AUX3 = 3.*B13*AZ + B12
      AUX4 = 3.*B03*AZ + B02
      vel(6) = 2.*(((AUX1*AX + AUX2)*AX + AUX3)*AX + AUX4)

      if (veltype .lt. 0) then
         do i = 1, 6
            vel(i) = vel(i) / vpvs_ratio
         end do
      end if

 999  end

c********************************************************************
c Unterprogramm VELONLY
c Errechnung NUR der Geschwindigkeit an einem Punkt (X0,Z0)
c
      subroutine velonly(x0, z0, nx_grid, nz_grid, x_grid, z_grid, v_grid,
     &     veltype, vpvs_ratio, actlayer,
     &     vel,
     &     MAXP_XGRID, MAXP_ZGRID)

      implicit none

      integer MAXP_XGRID, MAXP_ZGRID
      real x0, z0, vpvs_ratio
      integer nx_grid, nz_grid, veltype, actlayer
      real x_grid(MAXP_XGRID), z_grid(MAXP_ZGRID), v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real vel

      integer i, ix1,iz1
      real xlast, zlast
      real xa, xb, hx, hz, za, zb, ax, az
      real b00,b02,b20,b22,c00,c02,c20,c22,d00,d02,d20,d22,e00,e02,e20,e22
      real aux1,aux2,aux3,aux4
      
      integer ix2,iz2,layer
      real b01,b03,b10,b11,b12,b13,b21,b23,b30,b31,b32,b33
      real d10,d12,d30,d32

      common layer,b01,b03,b10,b11,b12,b13,b21,b23,b30,b31,b32,b33,b00,b02,b20,b22,ix2,iz2,xlast,zlast

C     DETERMINATION OF VELOCITY AND ITS DERIVATIVES
C     FOR BICUBIC POLYNOMIAL APPROXIMATION

      do i = 2, nx_grid
         ix1 = i
         if (x0 .lt. x_grid(ix1)) goto 2
      end do
    2 continue
c Punkt ist ausserhalb des velocity network (z-Wert zu klein)
      if (z0 .lt. z_grid(1)) then
         iz1 = 2
         goto 4
c         write(*,*) "Oops! Checking point outside the network! ---"
c         goto 999
      end if
      do i = 2, nz_grid
         iz1 = i
         if (z0 .lt. z_grid(iz1)) goto 4
      end do
    4 continue

      if (layer .ne. actlayer .or. ix1 .ne. ix2 .or. iz1 .ne. iz2) then
         ix2 = ix1
         iz2 = iz1
         ix1 = ix1 - 1
         iz1 = iz1 - 1

         if (ix1 .lt. 1) then
            write(*,*) 'a'
         end if

         xlast = x_grid(ix2-1)
         zlast = z_grid(iz2-1)

         B20 = v_grid(2,iz1,ix1)
         B02 = v_grid(3,iz1,ix1)
         B22 = v_grid(4,iz1,ix1)
         B00 = v_grid(1,iz1,ix1)

         C20 = v_grid(2,iz1,ix1+1)
         C02 = v_grid(3,iz1,ix1+1)
         C22 = v_grid(4,iz1,ix1+1)
         C00 = v_grid(1,iz1,ix1+1)

         d20 = v_grid(2,iz1+1,ix1)
         d02 = v_grid(3,iz1+1,ix1)
         d22 = v_grid(4,iz1+1,ix1)
         d00 = v_grid(1,iz1+1,ix1)

         e20 = v_grid(2,iz1+1,ix1+1)
         e02 = v_grid(3,iz1+1,ix1+1)
         e22 = v_grid(4,iz1+1,ix1+1)
         e00 = v_grid(1,iz1+1,ix1+1)

         hx = x_grid(ix2) - xlast
         hz = z_grid(iz2) - zlast
         XA = 3.*HX
         ZA = 3.*HZ
         XB = HX/3.
         ZB = HZ/3.

         D32 = (e22-d22) / XA
         D30 = (e20-d20) / XA
         B30 = (c20-b20) / XA
         B32 = (c22-b22) / XA
         D12 = (e02-d02) / HX - XB * (e22 + 2.*d22)
         D10 = (e00-d00) / HX - XB * (e20 + 2.*d20)
         B10 = (c00-b00) / HX - XB * (C20 + 2.*B20)
         B12 = (c02-b02) / HX - XB * (C22 + 2.*B22)
         B03 = (d02-b02) / ZA
         B13 =(D12-B12) / ZA
         B23 = (d22-b22) / ZA
         B33 = (D32-B32) / ZA
         B01 = (D00-B00) / HZ - ZB * (D02 + 2.*B02)
         B11 = (D10-B10) / HZ - ZB * (D12 + 2.*B12)
         B21 = (D20-B20) / HZ - ZB * (D22 + 2.*B22)
         B31 = (D30-B30) / HZ - ZB * (D32 + 2.*B32)
         layer = actlayer
      end if

      ax = x0 - xlast
      az = z0 - zlast

      AUX1 = ((B33*AZ + B32)*AZ + B31)*AZ + B30
      AUX2 = ((B23*AZ + B22)*AZ + B21)*AZ + B20
      AUX3 = ((B13*AZ + B12)*AZ + B11)*AZ + B10
      AUX4 = ((B03*AZ + B02)*AZ + B01)*AZ + B00
      vel = ((AUX1*AX + AUX2)*AX + AUX3)*AX + AUX4

      if (veltype .lt. 0) then
         vel = vel / vpvs_ratio
      end if

      end

c********************************************************************
c Unterprogramm VEL_INITIALIZE
c Initialisierung fuer die Routinen VELOCITY und VELONLY
c Die Variablen des COMMON-Blocks werden hier gesetzt!
c
      subroutine vel_initialize(x0, z0, nx_grid, nz_grid, x_grid, z_grid, v_grid,
     &     actlayer,
     &     MAXP_XGRID, MAXP_ZGRID)

      implicit none

      integer MAXP_XGRID, MAXP_ZGRID
      real x0, z0
      integer nx_grid, nz_grid, actlayer
      real x_grid(MAXP_XGRID), z_grid(MAXP_ZGRID), v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      integer i, ix1,iz1
      real xlast, zlast
      real xa, xb, hx, hz, za, zb
      real b00,b02,b20,b22,c00,c02,c20,c22,d00,d02,d20,d22,e00,e02,e20,e22
      
      integer ix2,iz2,layer
      real b01,b03,b10,b11,b12,b13,b21,b23,b30,b31,b32,b33
      real d10,d12,d30,d32

      common layer,b01,b03,b10,b11,b12,b13,b21,b23,b30,b31,b32,b33,b00,b02,b20,b22,ix2,iz2,xlast,zlast

      do i = 2, nx_grid
         ix1 = i
         if (x0 .lt. x_grid(ix1)) goto 2
      end do
    2 continue
c Punkt ist ausserhalb des velocity network (z-Wert zu klein)
      if (z0 .lt. z_grid(1)) then
         iz1 = 2
         goto 4
      end if
      do i = 2, nz_grid
         iz1 = i
         if (z0 .lt. z_grid(iz1)) goto 4
      end do
    4 continue

      ix2 = ix1
      iz2 = iz1
      ix1 = ix1 - 1
      iz1 = iz1 - 1

      xlast = x_grid(ix2-1)
      zlast = z_grid(iz2-1)

      B20 = v_grid(2,iz1,ix1)
      B02 = v_grid(3,iz1,ix1)
      B22 = v_grid(4,iz1,ix1)
      B00 = v_grid(1,iz1,ix1)
      
      C20 = v_grid(2,iz1,ix1+1)
      C02 = v_grid(3,iz1,ix1+1)
      C22 = v_grid(4,iz1,ix1+1)
      C00 = v_grid(1,iz1,ix1+1)

      d20 = v_grid(2,iz1+1,ix1)
      d02 = v_grid(3,iz1+1,ix1)
      d22 = v_grid(4,iz1+1,ix1)
      d00 = v_grid(1,iz1+1,ix1)

      e20 = v_grid(2,iz1+1,ix1+1)
      e02 = v_grid(3,iz1+1,ix1+1)
      e22 = v_grid(4,iz1+1,ix1+1)
      e00 = v_grid(1,iz1+1,ix1+1)

      hx = x_grid(ix2) - xlast
      hz = z_grid(iz2) - zlast
      XA = 3.*HX
      ZA = 3.*HZ
      XB = HX/3.
      ZB = HZ/3.

      D32 = (e22-d22) / XA
      D30 = (e20-d20) / XA
      B30 = (c20-b20) / XA
      B32 = (c22-b22) / XA
      D12 = (e02-d02) / HX - XB * (e22 + 2.*d22)
      D10 = (e00-d00) / HX - XB * (e20 + 2.*d20)
      B10 = (c00-b00) / HX - XB * (C20 + 2.*B20)
      B12 = (c02-b02) / HX - XB * (C22 + 2.*B22)
      B03 = (d02-b02) / ZA
      B13 =(D12-B12) / ZA
      B23 = (d22-b22) / ZA
      B33 = (D32-B32) / ZA
      B01 = (D00-B00) / HZ - ZB * (D02 + 2.*B02)
      B11 = (D10-B10) / HZ - ZB * (D12 + 2.*B12)
      B21 = (D20-B20) / HZ - ZB * (D22 + 2.*B22)
      B31 = (D30-B30) / HZ - ZB * (D32 + 2.*B32)
      layer = actlayer

c      write(*,*) layer,b01,b03,b10,b11,b12,b13,b21,b23,b30,b31,b32,b33,b00,b02,b20,b22,ix2,iz2,xlast,zlast

      end

c********************************************************************
c Unterprogramm VELONLY_LINEAR
c Diese Routine entstammt dem Programmpaket SEIS88, Originalname VELOC aus MODIS.F
c Veraenderungen sind vorgenommen worden
c
c     DETERMINATION OF VELOCITY FROM VELOCITY ISOLINES
c
      subroutine velonly_linear(x, z, layer, npoints_int, x_int, z_int,
     &     b_int, c_int, d_int, v_bottom, v_top, veltype, vpvs_ratio,
     &     vel,
     &     MAX_INT, MAXPOINTS_INT)

      implicit none

      integer MAX_INT, MAXPOINTS_INT
      integer npoints_int(MAX_INT), layer, veltype
      real x_int(MAXPOINTS_INT,MAX_INT), z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT), c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT)
      real v_bottom(MAX_INT), v_top(MAX_INT), vpvs_ratio(MAX_INT)
      real x, z, vel
      integer i, j, nlay, n
      real zz(2)
      real aux, aux1, aux2
      
      n = 0
      do nlay = layer, layer + 1
         n = n + 1
         do i = 1, npoints_int(nlay) - 1
            j = i
            if (x .lt. x_int(i+1,nlay)) goto 3
         end do
 3       aux = x - x_int(j,nlay)
         zz(n) = ((d_int(j,nlay)*aux + c_int(j,nlay))*aux + b_int(j,nlay))*aux + z_int(j,nlay)
      end do
      aux = v_bottom(layer) - v_top(layer)
      aux1 = zz(2) - zz(1)
      aux2 = z - zz(1)

      if (abs(aux1) .gt. .000001) then
         vel = (aux / aux1) * aux2 + v_top(layer)
      else
         vel = v_top(layer)
      end if

      if (veltype .gt. 0) goto 999
      vel = vel / vpvs_ratio(layer)
    
 999  end


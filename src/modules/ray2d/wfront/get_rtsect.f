c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm GET_RTSECT / Hauptprogramm WFRONT
c Aufruf des Unterprogramms VELOCITY
c Berechnet Hilfsgroessen fuer Reflection/Transmission
c am Schnittpunkt (SECT(1),SECT(2)):
c  rthalt(1) = xint_normal   : Einheitsvektor der Normalen auf der Grenzflaeche
c  rthalt(2) = zint_normal   :  -"-
c  rthalt(3) = px_unit : Einheitsvektor der Slowness (in lok. kart. Koordinaten)
c  rthalt(4) = pz_unit :  -"-
c  rthalt(5) = radcurv : Radius der Kruemmung der Grenzflaeche
c  rthalt(6) = akapa1
c  rthalt(7) = vs1
c  rthalt(8) = vel(1)  : Geschwindigkeit
c
      subroutine get_rtsect(sect,dzdx,d2zdx2,
     &     nx_grid,nz_grid,x_grid,z_grid,
     &     v_grid,ptos,veltype,actlayer,
     &     rtsect,
     &     N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer N_RTPARAM,MAXP_XGRID,MAXP_ZGRID
      integer veltype,actlayer
      real sect(4)
      real dzdx,d2zdx2,ptos
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real rtsect(N_RTPARAM)

      real xint_normal,zint_normal,px_unit,pz_unit,radcurv,akapa1,vs1,vel(6)
c      real p

      call velocity(sect(1), sect(2), nx_grid, nz_grid,
     &     x_grid, z_grid, v_grid,
     &     veltype, ptos,actlayer,
     &     vel,
     &     MAXP_XGRID, MAXP_ZGRID)

c Statt vel(1) zu benuetzen, muesste man vielleicht besser |p| nehmen !!!

c      p = sqrt(sect(3)*sect(3) + sect(4)*sect(4))
      zint_normal = 1.0/sqrt(1.0 + dzdx*dzdx)
      xint_normal = -dzdx * zint_normal
      pz_unit = (xint_normal*sect(3) + zint_normal*sect(4))*vel(1)
c      pz_unit = (xint_normal*sect(3) + zint_normal*sect(4))/p

      if (pz_unit .gt. 0.0) then
         xint_normal = -xint_normal
         zint_normal = -zint_normal
         pz_unit = -pz_unit
      end if

      px_unit = (xint_normal*sect(4) - zint_normal*sect(3))*vel(1)
c      px_unit = (xint_normal*sect(4) - zint_normal*sect(3))/p
      radcurv = d2zdx2*zint_normal*zint_normal*zint_normal

      akapa1 = vel(2)*sect(4) - vel(3)*sect(3)
c      vs1 = (vel(2)*sect(3) + vel(3)*sect(4))/p
      vs1 = (vel(2)*sect(3) + vel(3)*sect(4))*vel(1)

      rtsect(1) = xint_normal
      rtsect(2) = zint_normal
      rtsect(3) = px_unit
      rtsect(4) = pz_unit
      rtsect(5) = radcurv
      rtsect(6) = akapa1
      rtsect(7) = vs1
      rtsect(8) = vel(1)
c      rtsect(8) = 1./p

      end


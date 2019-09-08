c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm GET_FLAGRAY / Hauptprogramm WFRONT
c Wird aufgerufen von DO_INTERSECTION
c Aufruf der Unterprogramme DERY, RUNKUTTA2
c Bestimmt, ob der Strahl RAY unter- oder oberhalb am Punkt (XINT,ZINT) vorbeigeht
c FLAGRAY=+1: Strahl verlaeuft unterhalb des Punktes (ZRAY > ZINT)
c FLAGRAY=-1: Strahl verlaeuft oberhalb ...(ZRAY < ZINT)
c FLAGRAY= 0: normal...
c
c DTNEW: Hilfe fuer INTERSECT_RAY: Zeitschritt
c
      subroutine get_flagray(xint,zint,ray,dtray,nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     flagray,dtnew,
     &     MAXP_XGRID,MAXP_ZGRID,N_PARAM)

      implicit none

      integer MAXP_XGRID,MAXP_ZGRID,N_PARAM,veltype,actlayer
      real ray(N_PARAM), raynew(N_PARAM), ptos, dtray, xint, zint
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      integer flagray
      real dtnew

      integer idiv,i, count, MAX_COUNTS
      parameter (MAX_COUNTS = 20)
      real divide, dtnewhalf,dt1,dt2,dx1,dz1,dx2,dz2,dxnew
      real kh1(N_PARAM)

      idiv = 20
      divide = 1./ float(idiv)

      dtnew = dtray * divide
      dtnewhalf = 0.5*dtnew

      call dery(ray,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &     kh1,
     &     4,MAXP_XGRID,MAXP_ZGRID)

      call runkutta2(ray,dtnew,dtnewhalf,kh1,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &     raynew,
     &     4,MAXP_XGRID,MAXP_ZGRID)

      dx1 = ray(1) - xint
      dz1 = ray(2) - zint
      dx2 = raynew(1) - xint
      dz2 = raynew(2) - zint

      do i = 2, idiv
         if (dx1*dx2 .le. 0.0) goto 10
         dtnew = dtray*divide*float(i)
         dtnewhalf = 0.5*dtnew
         call runkutta2(ray,dtnew,dtnewhalf,kh1,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &        raynew,
     &        4,MAXP_XGRID,MAXP_ZGRID)
         dx1 = dx2
         dz1 = dz2
         dx2 = raynew(1) - xint
         dz2 = raynew(2) - zint
      end do
      if (dx1*dx2 .gt. 0.0) then
c         write(*,*) "Error in GET_FLAG!!! Program will perhaps end on error!"
         dtnew = dtray
         flagray = 0
         goto 999
      end if

 10   continue

      dt2 = dtnew
      dt1 = dtnew - dtray*divide
      count = 0

      do while (dz1*dz2 .lt. 0.0 .and. count .lt. MAX_COUNTS)
         count = count + 1
         dtnew = 0.5 * (dt1 + dt2)
         dtnewhalf = 0.5*dtnew
         call runkutta2(ray,dtnew,dtnewhalf,kh1,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &        raynew,
     &        4,MAXP_XGRID,MAXP_ZGRID)
         dxnew = raynew(1) - xint
         if (dx1*dxnew .gt. 0.0) then
            dx1 = dxnew
            dz1 = raynew(2) - zint
            dt1 = dtnew
         else
            dx2 = dxnew
            dz2 = raynew(2) - zint
            dt2 = dtnew
         end if
      end do

      if (dz1 .gt. 0.0 .or. dz2 .gt. 0.0) then
         flagray = 1
      else
         flagray = -1
      end if

c Um sicher zu gehen, dass DTNEW wirklich in die neue Schicht laeuft!
c Abstand DX2 muss klein sein, damit man sicher sein kann, dass der Strahl bis DTNEW die Grenzflaeche
c nicht vielleicht zweimal schneidet! (bei fast senkrechter Grenzflaeche im Extremfall)
c
c ist nicht ganz ausgereift! (oder doch!)
      count = 0
      do while (abs(dx2) .gt. 0.0001 .and. count .lt. MAX_COUNTS)
         count = count + 1
         dtnew = 0.5 * (dt1 + dt2)
         dtnewhalf = 0.5*dtnew
         call runkutta2(ray,dtnew,dtnewhalf,kh1,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &        raynew,
     &        4,MAXP_XGRID,MAXP_ZGRID)
         dxnew = raynew(1) - xint
         if (dx1*dxnew .gt. 0.0) then
            dx1 = dxnew
            dt1 = dtnew
         else
            dx2 = dxnew
            dt2 = dtnew
         end if
      end do

c Diskussion, ob hier DT1 oder DT2 zu nehmen ist:
c DT2 pro: Man ist sicher, dass die Grenzflaeche wirklich ueberschritten wurde.
c DT1 pro: Man ist sicher, dass nicht DT2 direkt in INTERSECT_RAY genommen wird, was der Fall sein kann, wenn DZ2 < EPSILON ist.
c DT1 : Auch mit DT1 ist Grenzflaeche wahrscheinlich bereits ueberschritten;
c   naemlich immer, wenn count < MAX_COUNTS (obere Schleife endet mit DZ2*DZ1 > 0.0). ABER: Auch ist's moeglich,
c   dass der Strahl wie eine Schale den Punkt umschliesst, d.h. mit Anfangs- sowie Endpunkt unter oder
c   ueber dem Eckpunkt liegt. Dann ist es moeglich, dass trotz count < MAX_COUNTS DT1 noch nicht bis in die neue Schicht
c   laeuft.
c FAZIT: Es muss die Moeglichkeit geschaffen werden, dass DT1 genommen werden kann... ODER:
c !!!!! In INTERSECT_RAY zwar abs(chkmid) < epsilon als Treffer nehmen, aber RAYNEW muss in der alten Schicht liegen !!!!!
c Dann wird DT2 beibelassen.

      dtnew = dt2

 999  end


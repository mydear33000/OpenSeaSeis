c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm GET_BOUNDRAY / Hauptprogramm WFRONT
c Wird aufgerufen von WFRONT
c Aufruf der Unterprogramme RUNKUTTA, RUNKUTTA2, DERY, INTERSECT_WFRONT, INTERPOL_BOUND, GET_RTSECT, INTERSECT_RAY, GET_Z
c Berechnet den genauen boundary ray. Dieser liegt (normalerweise) zwischen RAY1 und RAY2, und
c wird zu einem Zeitpunkt zwischen ACTTIME und ACTTIME-DTRAY angesetzt.
c Vorgehensweise:
c RAY1 liegt in der neuen Schicht, RAY2 in der alten.
c Wie in WFRONT wird ein Checkwert zwischen zwei Strahlen RAY1 und RAY2 berechnet (CHKBOUND).
c (Checkwert: inwieweit ist Wellenfront senkrecht zu Grenzflaeche am Schnittpunkt zwischen
c der Grenzflaeche und einer Geraden durch RAY1 und RAY2?)
c Die Strahlen RAY1 und RAY2 werden um DTRAY zurueckpropagiert.
c Von dort aus wird nun iterativ durch jeweiliges Halbieren von DTRAY ein Checkwert
c bestimmt, bis dieser kleiner als EPSILON ist.
c FLAG:
c   =  0 : boundary ray liegt zwischen Strahlen RAY1 und RAY2
c   = -1 : boundary ray liegt ausserhalb von Strahl RAY2, d.h. RAY2 schneidet
c          die Grenzflaeche zweimal und muss geloescht werden/neu bearbeitet werden!
c   = +1 : Abbruch wegen caustic zone!
c INF_PART: infinitesimale Entfernung von Grenzflaeche; an dieser Stelle wird der
c    boundary ray angesetzt.
c    Genauer: Der bound.ray wird zwischen ray1 und ray2 interpoliert. Die Entfernung zwischen
c    ray1 und ray2 sei 1.0; Entfernung zwischen ray1 und dem eigentlichen boundary ray sei
c    PART. Dann wird der boundary ray im Abstand PART+INF_PART von ray1 angesetzt. 
c
c Ausgabe:
c BOUNDTMP: Der (vorbeilaufende) boundary ray an dem Ort, an dem er interpoliert wird
c BOUNDRAY: Der (vorbeilaufende) boundary ray auf der aktuellen Wellenfront
c .....
c
      subroutine get_boundray(ray1,ray2,amplitude1,amplitude2,phase1,phase2,angle1,angle2,
     &     chkbound0,acttime,kmah,dtray, dthalf, interface,cross1,actlayer,
     &     npoints_int,x_int,z_int,b_int,c_int,d_int,iii,iii2, veltype, ptos,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     boundtmp,ttmp,boundray,amplboundray,phaseboundray,angleboundray,kmahboundray,flag,
     &     crossect,sect,tsect,rtsect,fictsect,amplsect,phasesect,kmahsect,anglesect,
     &     N_PARAM, MAX_INT, MAXPOINTS_INT, MAXP_XGRID,MAXP_ZGRID, N_RTPARAM, INF_DIST)

      implicit none

      integer MAX_INT, MAXPOINTS_INT, N_PARAM, MAXP_XGRID, MAXP_ZGRID, N_RTPARAM, MAX_ALLINT
      integer npoints_int(MAX_INT),iii(MAXPOINTS_INT,MAX_INT),iii2(MAXPOINTS_INT,MAX_INT)
      real x_int(MAXPOINTS_INT,MAX_INT), z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT), c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT)
      real ptos
      integer veltype, interface, cross1, actlayer
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID),v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real ray1(N_PARAM),ray2(N_PARAM),amplitude1,amplitude2,phase1,phase2,angle1,angle2
      real chkbound0, acttime, dtray, dthalf
      integer kmah
 
      real boundtmp(N_PARAM),boundray(N_PARAM),amplboundray,phaseboundray,angleboundray,ttmp
      real sect(N_PARAM),rtsect(N_RTPARAM),tsect,amplsect,phasesect,anglesect
      integer kmahboundray,flag,crossect,fictsect,kmahsect

      integer MAX_COUNTS
      parameter (MAX_COUNTS = 30)
      real ray1old(N_PARAM),ray2old(N_PARAM),ray1new(N_PARAM),ray2new(N_PARAM)
      real chkright, chkleft, chkmid, dtleft, dtright, dtnew, dtnewhalf, raysect(4)
      real kh1one(N_PARAM), kh1two(N_PARAM)

      real s, sbound, part, INF_DIST, INF_PART
      real epsilon, dzdx, tnew2, dz, dtnew2, zint, z2int, d2zdx2, dummy
      integer i, count, actp, actp2, abbruch

c      if (NPARAM_LOCAL .ne. N_PARAM) then
c         write(*,*) "GET_BOUNDRAY: WARNING! Parameter NPARAM_LOCAL not equal N_PARAM!"
c      end if

      abbruch = 0
      do i = 1, 4
         raysect(i) = 0.0
      end do
      MAX_ALLINT = MAX_INT
c      INF_PART = 0.01
      epsilon = 0.0001
      flag = 0
      chkright = chkbound0

      if (abs(chkright) .lt. epsilon) then
         do i = 1, N_PARAM
            ray1new(i) = ray1(i)
            ray2new(i) = ray2(i)
         end do
         ttmp = acttime
      else
c-------------
c Zurueckpropagation um DTRAY, Bestimmung der alten Werte RAYOLD:
c
         dtnew = -dtray
         dtnewhalf = -dthalf
         call runkutta(ray1,dtnew,dtnewhalf,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        ray1old,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         call runkutta(ray2,dtnew,dtnewhalf,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        ray2old,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
c--------------
c 'linker' (alter) Checkwert (CHKLEFT) wird bestimmt:
c
         call intersect_wfront(ray1old,ray2old,interface,cross1,
     &        npoints_int,x_int,z_int,b_int,c_int,d_int,iii,iii2,veltype,ptos,actlayer,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        dzdx,raysect,abbruch,
     &        MAX_INT,MAXPOINTS_INT,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         if (abbruch .ne. 0) goto 990
         chkleft = dzdx*raysect(3) - raysect(4)
         dtleft = 0.0
         dtright = dtray
c
c---------------
c 'mittlerer' Checkwert wird bestimmt
c
         dtnew = 0.5 * dtray
         dtnewhalf = 0.5 * dtnew
         call dery(ray1old,nx_grid,nz_grid,x_grid,z_grid,
     &        v_grid,ptos,veltype,actlayer,
     &        kh1one,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         call dery(ray2old,nx_grid,nz_grid,x_grid,z_grid,
     &        v_grid,ptos,veltype,actlayer,
     &        kh1two,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         call runkutta2(ray1old,dtnew,dtnewhalf,kh1one,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        ray1new,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         call runkutta2(ray2old,dtnew,dtnewhalf,kh1two,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        ray2new,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)

         call intersect_wfront(ray1new,ray2new,interface,cross1,
     &        npoints_int,x_int,z_int,b_int,c_int,d_int,iii,iii2,veltype,ptos,actlayer,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        dzdx,raysect,abbruch,
     &        MAX_INT,MAXPOINTS_INT,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         if (abbruch .ne. 0) goto 990
         chkmid = dzdx*raysect(3) - raysect(4)
         
c-----------------
c Iterative Schleife
         count = 0
c count > MAX_COUNTS sollte nur vorkommen, falls Interface eine scharfe Ecke hat (iii(.,.) = -1!)
c dann ist dies auch OK!
         do while (abs(chkmid) .gt. epsilon .and. count .lt. MAX_COUNTS)
            count = count + 1
            if ((chkright*chkmid) .lt. 0.0) then
               dtleft = dtnew
               chkleft = chkmid
            else
               dtright = dtnew
               chkright = chkmid
            end if
            dtnew = 0.5 * (dtleft + dtright)
            dtnewhalf = 0.5 * dtnew
            call runkutta2(ray1old,dtnew,dtnewhalf,kh1one,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           ptos,veltype,actlayer,
     &           ray1new,
     &           N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            call runkutta2(ray2old,dtnew,dtnewhalf,kh1two,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           ptos,veltype,actlayer,
     &           ray2new,
     &           N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            call intersect_wfront(ray1new,ray2new,interface,cross1,
     &           npoints_int,x_int,z_int,b_int,c_int,d_int,iii,iii2,veltype,ptos,actlayer,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           dzdx,raysect,abbruch,
     &           MAX_INT,MAXPOINTS_INT,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            if (abbruch .ne. 0) goto 990
            chkmid = dzdx*raysect(3) - raysect(4)
         end do
c------------------

         ttmp = acttime - dtray + dtnew
      end if

c      s = (ray1new(1)-ray2new(1))**2 + (ray1new(2)-ray2new(2))**2
c      sbound = (ray1new(1)-raysect(1))**2 + (ray1new(2) - raysect(2))**2
c      part = sqrt(sbound/s) + INF_PART
c-----------

      s = sqrt((ray1new(1)-ray2new(1))**2 + (ray1new(2)-ray2new(2))**2)
      sbound = sqrt((ray1new(1)-raysect(1))**2 + (ray1new(2) - raysect(2))**2)
      part = sbound/s
c keine Interpolation, wenn ein Strahl RAY2 die Grenzflaeche zweimal schneidet
      if (part .gt. 1.0) then
         flag = -1
         goto 999
      end if

      INF_PART = INF_DIST/s
      part = part + INF_PART
c keine Interpolation, wenn Strahl RAY2 bereits nah genug an Grenzflaeche vorbeilaeuft:
c Eigentlich muesste jetzt der zweite boundary ray trotzdem interpoliert werden...
      if (part .gt. 1.0) then
         flag = 1
         goto 999
      end if

c------------------------
c Erster Boundray (der, der vorbeilaeuft)
c (BOUNDTMP wird bestimmt)
      call interpol_bound(ray1new,ray2new,amplitude1,amplitude2,phase1,phase2,angle1,angle2,
     &     acttime,kmah,part,ttmp,veltype, ptos, dtnew,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,actlayer,
     &     boundtmp,amplboundray,phaseboundray,kmahboundray,angleboundray,abbruch,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)
      if (abbruch .ne. 0) goto 990

c Propagieren zur aktuellen Wellenfront
      dtnewhalf = 0.5*dtnew
      call runkutta(boundtmp,dtnew,dtnewhalf,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     boundray,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

c========================================================================================
c Zweiter Boundray (der, der die Grenzflaeche gerade noch schneidet)
c
c nur, wenn RAY1 nicht bereits im Genauigkeitsgebiet liegt!
      if (part-2.*INF_PART .gt. 0.0) then
         part = part - 2.*INF_PART
         
         call interpol_bound(ray1new,ray2new,amplitude1,amplitude2,phase1,phase2,angle1,angle2,
     &        acttime,kmah,part,ttmp,veltype, ptos, dtnew,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,actlayer,
     &        ray2old,amplsect,phasesect,kmahsect,anglesect,abbruch,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         if (abbruch .ne. 0) goto 990

c Falls der neue Strahl nicht in der neuen Schicht liegt, wird er geloescht:
c Bedeutet: Zu spitzer Winkel in Modell (geringer als 90 Grad)
         call get_z(ray2old(1),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &        b_int(1,interface),npoints_int(interface),
     &        z2int,actp2,
     &        MAXPOINTS_INT)
         if (cross1 .lt. 0 .and. z2int .lt. ray2old(2)) then
            crossect = 0
            goto 999
         else if (cross1 .gt. 0 .and. z2int .gt. ray2old(2)) then
            crossect = 0
            goto 999
         end if
c+++++++++++++++++++++
c ray2old muss nun zurueckpropagiert werden und mit der Grenzflaeche geschnitten werden...:
         dtnew = dtnew - dtray
         dtnewhalf = 0.5*dtnew
         call runkutta(ray2old,dtnew,dtnewhalf,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        ray2new,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)

         dz = ray2old(2) - z2int
         call get_z(ray2new(1),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &        b_int(1,interface),npoints_int(interface),
     &        zint,actp,
     &        MAXPOINTS_INT)
         tnew2 = acttime - dtray
         dtnew2 = -dtnew
         dtnew = -dtray
         dtnewhalf = 0.5*dtnew
         count = 0
         do while ((ray2new(2)-zint)*dz .gt. 0.0)
            count = count + 1
            if (count .eq. 10) then
               crossect = 0
               goto 999
            end if
            tnew2 = tnew2 - dtray
            dtnew2 = dtray
            actp2 = actp
            z2int = zint
            do i = 1, N_PARAM
               ray2old(i) = ray2new(i)
            end do
            call runkutta(ray2old,dtnew,dtnewhalf,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           ptos,veltype,actlayer,
     &           ray2new,
     &           N_PARAM,MAXP_XGRID,MAXP_ZGRID)

            call get_z(ray2new(1),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &           b_int(1,interface),npoints_int(interface),
     &           zint,actp,
     &           MAXPOINTS_INT)
         end do
c+++++++++++++++++++++++++++++
c jetzt liegt RAY2NEW in der alten Schicht, RAY2OLD in der neuen
         call intersect_ray(ray2new,ray2old,actp2,cross1,z2int,tnew2,
     &        0,dtnew2,1,dummy,
     &        x_int,d_int,c_int,b_int,z_int,npoints_int,iii,iii2,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        actlayer,ptos,veltype,0,
     &        sect,tsect,fictsect,dzdx,d2zdx2,
     &        N_PARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)

         call get_rtsect(sect,dzdx,d2zdx2,
     &        nx_grid,nz_grid,x_grid,z_grid,
     &        v_grid,ptos,veltype,actlayer,
     &        rtsect,
     &        N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
         crossect = cross1
      end if

      goto 999
 990  flag = 1

 999  end

c************************************************************************
c Unterprogramm INTERPOL_BOUND / Hauptprogramm WFRONT
c Wird aufgerufen von GET_BOUNDRAY
c Interpolation eines boundary rays.
c Dieser liegt zwischen RAY1 und RAY2 (PART >= 0.0)
c oder ausserhalb von RAY2 (PART < 0.0)
c TTMP ist der totale Zeitpunkt, zu dem der boundary ray loslaeuft.
c qdynnew,pdynnew werden durch einfache Mittelung bestimmt.
c
      subroutine interpol_bound(ray1,ray2,amplitude1,amplitude2,phase1,phase2,angle1,angle2,
     &     acttime,kmah,part,ttmp,veltype, ptos, dtnew,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,actlayer,
     &     raynew,amplitudenew,phasenew,kmahnew,anglenew,abbruch,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer N_PARAM,MAXP_XGRID,MAXP_ZGRID,actlayer
      real ray1(N_PARAM),ray2(N_PARAM),amplitude1,amplitude2,phase1,phase2,angle1,angle2
      real acttime,ttmp
      real part
      real raynew(N_PARAM), amplitudenew, phasenew, anglenew
      integer kmah,kmahnew,abbruch
      real tnew, ttmp2
c--------
      real vel(6),dtnew
      integer i
      integer veltype
      real ptos
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)
c---------

      abbruch = 0

      do i = 1, 2
         raynew(i) = part * (ray2(i) - ray1(i)) + ray1(i)
      end do

      ttmp2 = ttmp

c Ueberlegung: Man sollte wie bei GET_RAYONINT die beiden unteren Ergebnisse 'quadratisch' mitteln,
c oder vielleicht direkt GET_RAYONINT aufrufen?

      if (part .lt. 0.5) then
         if (abs(ray1(5)) .lt. 0.00001) then
            write(*,*) 'INTERPOL_BOUND: CAUSTIC!'
            abbruch = 1
            goto 999
         end if
         call velocity(ray1(1), ray1(2), nx_grid, nz_grid,
     &        x_grid, z_grid, v_grid,veltype, ptos,actlayer,
     &        vel,
     &        MAXP_XGRID, MAXP_ZGRID)
         call paraxial(ray1,ttmp,vel,
     &        raynew,tnew,N_PARAM)
      else
         if (abs(ray2(5)) .lt. 0.00001) then
            write(*,*) 'INTERPOL_BOUND: CAUSTIC!'
            abbruch = 1
            goto 999
         end if
         call velocity(ray2(1), ray2(2), nx_grid, nz_grid,
     &        x_grid, z_grid, v_grid, veltype, ptos,actlayer,
     &        vel,
     &        MAXP_XGRID, MAXP_ZGRID)
         call paraxial(ray2,ttmp,vel,
     &        raynew,tnew,N_PARAM)
      end if

      do i = 5, 6
         raynew(i) = part * (ray2(i) - ray1(i)) + ray1(i)
      end do
      amplitudenew = part * (amplitude2 - amplitude1) + amplitude1
      phasenew = part * (phase2 - phase1) + phase1
      anglenew = part * (angle2 - angle1) + angle1

      kmahnew = kmah

c Sollten die Strahlen konvergieren (kann eigentlich nicht vorkommen),
c das heisst tnew > acttime, gilt das folgende trotzdem:
      dtnew = acttime - tnew

 999  end


c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm GET_CORNERRAYS / Hauptprogramm WFRONT
c Wird aufgerufen von SORTSECT
c Aufruf der Unterprogramme INTERPOLATE2, DO_INTERSECTION, RUNKUTTA
c Interpoliert zwei Strahlen in die Ecke einer Blockstruktur
c RAY1 : Strahl, der obere Grenzflaeche schnitt
c RAY2 : Strahl, der untere Grenzflaeche schnitt
c
c Ausgabe:
c RAYUP  : interpolierter Strahl, der obere Granzflaeche schneidet
c FLAGUP : =0: RAYUP ist gleich RAY1; sonst: neuer Strahl!
c RAYDOWN  : interpolierter Strahl, der untere Grenzflaeche schneidet
c FLAGDOWN : =0: RAYDOWN ist gleich RAY2; sonst: neuer Strahl!
c ...
c
      subroutine get_cornerrays(ray1,time1,amplitude1,phase1,angle1,kmah,
     &     ray2,time2,amplitude2,phase2,angle2,
     &     dtray,actlayer,left_border,right_border,f_out,
     &     x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,veltype,ptos,
     &     rayup,timeup,rtup,amplitudeup,phaseup,angleup,kmahup,fictup,flagup,
     &     raydown,timedown,rtdown,amplitudedown,phasedown,angledown,kmahdown,fictdown,flagdown,
     &     error,
     &     N_PARAM,N_RTPARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer N_PARAM,N_RTPARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID
      real ray1(N_PARAM),ray2(N_PARAM),amplitude1,amplitude2,phase1,phase2,angle1,angle2
      real time1,time2,dtray
      integer kmah,actlayer,veltype,ihaltup,ihaltdown,f_out
      real left_border,right_border
      real x_int(MAXPOINTS_INT,MAX_INT),z_int(MAXPOINTS_INT,MAX_INT),d_int(MAXPOINTS_INT,MAX_INT)
      real c_int(MAXPOINTS_INT,MAX_INT),b_int(MAXPOINTS_INT,MAX_INT), ptos
      integer npoints_int(MAX_INT),iii(MAXPOINTS_INT,MAX_INT),iii2(MAXPOINTS_INT,MAX_INT),nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID),v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real rayup(N_PARAM),rtup(N_RTPARAM),timeup,raydown(N_PARAM),rtdown(N_RTPARAM),timedown
      real amplitudeup,amplitudedown,phaseup,phasedown,angleup,angledown
      integer flagup,flagdown,kmahup,kmahdown,fictup,fictdown,error
      
      integer N_INTERPOLATIONS, N_RTPARAM_LOCAL
      parameter(N_INTERPOLATIONS = 10, N_RTPARAM_LOCAL = 8)
      real raymid(N_PARAM),rayup0(N_PARAM),raydown0(N_PARAM),amplitudemid,phasemid,anglemid
      real rayold(N_PARAM),raynew(N_PARAM),raysect(N_PARAM),tsect,rtsect(N_RTPARAM_LOCAL)
      integer crossect, kmahsect, fictsect, kmahmid
      real acttime,time0,dtnew,dtnewhalf,chkbounddummy,spread_dummy,cos_dummy,dtnew0
      real part,partup,partdown
      integer crossdummy,flaginterpol,flagintersect,i,n,iloop,nloops

      if (N_RTPARAM .gt. N_RTPARAM_LOCAL) goto 990

      spread_dummy = 0.0
      cos_dummy = 2.0
      flaginterpol = 0
      flagintersect = 0
      ihaltdown = 1
      ihaltup = 1

      flagup = 0
      flagdown = 0


c==========================================================
c Rueckpropagation des Strahls mit der groesseren Zeit
c
      if (time1 .lt. time2) then
         dtnew = time1-time2
         do i = 1, N_PARAM
            rayold(i) = ray2(i)
         end do
c----
c Der Zeitschritt muss kleinergleich DTRAY sein!
         if (abs(dtnew) .gt. dtray) then
            nloops = abs(int(dtnew/dtray))
            dtnew = mod(dtnew,dtray)
            
            dtnew0 = -dtray
            dtnewhalf = 0.5*dtnew0
            do iloop = 1, nloops
               call runkutta(rayold,dtnew0,dtnewhalf,
     &              nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &              ptos,veltype,actlayer,
     &              raydown0,
     &              N_PARAM,MAXP_XGRID,MAXP_ZGRID)
               do i = 1, N_PARAM
                  rayold(i) = raydown0(i)
               end do
            end do
         end if
         dtnewhalf = 0.5*dtnew
         call runkutta(rayold,dtnew,dtnewhalf,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        raydown0,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
 
         amplitudedown = amplitude2
         phasedown = phase2
         angledown = angle2
         do i = 1, N_PARAM
            rayup0(i) = ray1(i)
         end do
         amplitudeup = amplitude1
         phaseup = phase1
         angleup = angle1
         time0 = time1
c---------------------------------------
      else
         dtnew = time2-time1
         do i = 1, N_PARAM
            rayold(i) = ray1(i)
         end do
c----
c Der Zeitschritt muss kleinergleich DTRAY sein!
         if (abs(dtnew) .gt. dtray) then
            nloops = abs(int(dtnew/dtray))
            dtnew = mod(dtnew,dtray)
            
            dtnew0 = -dtray
            dtnewhalf = 0.5*dtnew0
            do iloop = 1, nloops
               call runkutta(rayold,dtnew0,dtnewhalf,
     &              nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &              ptos,veltype,actlayer,
     &              rayup0,
     &              N_PARAM,MAXP_XGRID,MAXP_ZGRID)
               do i = 1, N_PARAM
                  rayold(i) = rayup0(i)
               end do
            end do
         end if
         dtnewhalf = 0.5*dtnew
         call runkutta(rayold,dtnew,dtnewhalf,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        rayup0,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         amplitudeup = amplitude1
         phaseup = phase1
         angleup = angle1
         do i = 1, N_PARAM
            raydown0(i) = ray2(i)
         end do
         amplitudedown = amplitude2
         phasedown = phase2
         angledown = angle2
         time0 = time2
      end if
c
c Rueckpropagation Ende
c================================================================

      partup = 0.0
      partdown = 1.0
      dtnew = 0.25 * dtray
      dtnewhalf = 0.5*dtnew
c
c N Interpolationen zwischen RAYUP und RAYDOWN
c
      do n = 1, N_INTERPOLATIONS
         part = 0.5 * (partdown + partup)
         call interpolate2(rayup0,amplitudeup,phaseup,kmah,angleup,
     &        raydown0,amplitudedown,phasedown,angledown,part,
     &        raymid,amplitudemid,phasemid,kmahmid,anglemid,
     &        N_PARAM)
         do i = 1, N_PARAM
            raynew(i) = raymid(i)
         end do
         crossect = 0
         acttime = time0
         do while (crossect .eq. 0)
            acttime = acttime + dtnew
            do i = 1, N_PARAM
               rayold(i) = raynew(i)
            end do
            call runkutta(rayold,dtnew,dtnewhalf,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           ptos,veltype,actlayer,
     &           raynew,
     &           N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            call do_intersection(rayold,actlayer,(acttime-dtnew),dtnew,ihaltup,ihaltdown,veltype,flagintersect,
     &           x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,ptos,left_border,right_border,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           raynew,kmah,chkbounddummy,crossdummy,
     &           raysect,tsect,rtsect,kmahsect,crossect,fictsect,
     &           MAX_INT,MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
            if (crossect .eq. 0 .and. (raynew(1) .le. left_border .or. raynew(1) .ge. right_border)) then
               goto 500
            end if
         end do
         if (crossect .lt. 0) then
            partup = part
            do i = 1, N_PARAM
c Muss dat oder nich?
c               rayup0(i) = raymid(i)
               rayup(i) = raysect(i)
            end do
            amplitudeup = amplitudemid
            phaseup = phasemid
            angleup = anglemid
            flagup = 1
            fictup = fictsect
            timeup = tsect
            kmahup = kmahsect
            if (ihaltup .ne. 0) then
               do i = 1, N_RTPARAM
                  rtup(i) = rtsect(i)
               end do
            end if
         else
            partdown = part
            do i = 1, N_PARAM
c Muss dat oder nich?
c               raydown0(i) = raymid(i)
               raydown(i) = raysect(i)
            end do
            amplitudedown = amplitudemid
            phasedown = phasemid
            angledown = anglemid
            flagdown = 1
            timedown = tsect
            kmahdown = kmahsect
            fictdown = fictsect
            if (ihaltdown .ne. 0) then
               do i = 1, N_RTPARAM
                  rtdown(i) = rtsect(i)
               end do
            end if
         end if
      end do
      goto 999

c--------------------------------------------------------
c Not corner ray, but small layer:
c Keine Interpolation, aber moeglichst Warnung ausgeben
c Ausnahme: Corner Point liegt auf Modellrand. Dann keine Warnung ausgeben!
c
 500  continue

      if (ray1(3) .ge. 0.0) then
         if ((z_int(actlayer,npoints_int(actlayer)) .eq. z_int(actlayer+1,npoints_int(actlayer+1))) .and.
     &        (z_int(actlayer,npoints_int(actlayer)-1) .ne. z_int(actlayer+1,npoints_int(actlayer+1)-1))) then
c nothing
         else
            write(*,1000,err=980) actlayer, ray1(1),ray1(2),ray2(1),ray2(2)
            write(f_out,1000,err=980) actlayer, ray1(1),ray1(2),ray2(1),ray2(2)
         end if
      else if ((z_int(actlayer,1) .eq. z_int(actlayer+1,1)) .and.
     &        (z_int(actlayer,2) .ne. z_int(actlayer+1,2))) then
c nothing
      else
         write(*,1000,err=980) actlayer, ray1(1),ray1(2),ray2(1),ray2(2)
         write(f_out,1000,err=980) actlayer, ray1(1),ray1(2),ray2(1),ray2(2)
      end if

      flagup = 0
      flagdown = 0

      goto 999

 1000 format("WARNING: SORTSECT: Layer no. ",i2," is very small.",/,
     &     "Two neighbour rays intersect both interfaces at points (",f6.2,",",f6.2,") and (",f6.2,",",f6.2,")",
     &     /,"Dimine the input variable SPREAD_MAX, the minimum distance between two rays.",/)

 980  write(*,*) "GET_CORNERRAYS: Error while writing in output file no.",f_out
      error = 1
      goto 999
 990  write(*,*) "GET_CORNERRAYS: PARAMETER 'N_RTPARAM_LOCAL' too low!"
      write(*,*) "Change it to at least ",N_RTPARAM
      error = 1

 999  end


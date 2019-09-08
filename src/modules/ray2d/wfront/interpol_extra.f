c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm INTERPOL_EXTRA / Hauptprogramm WFRONT
c Wird aufgerufen von SORTSECT
c Aufruf der Unterprogramme INTERPOLATE, DO_INTERSECTION, RUNKUTTA
c Interpoliert zwei Strahlen, damit zwischen ihnen paraxial approximiert werden darf (Abstand!)
c und zwar fuer den Fall, dass zwischen RAY1 und RAY2 ein Stuetzpunkt liegt.
c
c FICT1 : Index des Stuetzpunktes von RAY1
c
c Ausgabe:
c RAYLEFT  : interpolierter Strahl, der auf der Seite von RAY1 ist
c FLAGLEFT : =0: RAYLEFT ist gleich RAY1; sonst: neuer Strahl!
c RAYRIGHT  : interpolierter Strahl, der auf der Seite von RAY2 ist
c FLAGRIGHT : =0: RAYRIGHT ist gleich RAY2; sonst: neuer Strahl!
c ...
c
      subroutine interpol_extra(ray1,time1,amplitude1,phase1,angle1,kmah,fict1,
     &     ray2,time2,amplitude2,phase2,angle2,
     &     dtray,actlayer,left_border,right_border,
     &     x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,veltype,ptos,
     &     rayleft,timeleft,rtleft,amplitudeleft,phaseleft,angleleft,kmahleft,fictleft,flagleft,
     &     rayright,timeright,rtright,amplituderight,phaseright,angleright,kmahright,fictright,flagright,
     &     error,
     &     N_PARAM,N_RTPARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer N_PARAM,N_RTPARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID
      real ray1(N_PARAM),ray2(N_PARAM),amplitude1,amplitude2,phase1,phase2,angle1,angle2
      real time1,time2,dtray
      integer kmah,actlayer,veltype,ihaltup,ihaltdown,fict1
      real left_border,right_border
      real x_int(MAXPOINTS_INT,MAX_INT),z_int(MAXPOINTS_INT,MAX_INT),d_int(MAXPOINTS_INT,MAX_INT)
      real c_int(MAXPOINTS_INT,MAX_INT),b_int(MAXPOINTS_INT,MAX_INT), ptos
      integer npoints_int(MAX_INT),iii(MAXPOINTS_INT,MAX_INT),iii2(MAXPOINTS_INT,MAX_INT),nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID),v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real rayleft(N_PARAM),rtleft(N_RTPARAM),timeleft,rayright(N_PARAM),rtright(N_RTPARAM),timeright
      real amplitudeleft,amplituderight,phaseleft,phaseright,angleleft,angleright
      integer flagleft,flagright,kmahleft,kmahright,fictleft,fictright,error
      
      integer N_INTERPOLATIONS, N_RTPARAM_LOCAL
      parameter(N_INTERPOLATIONS = 5, N_RTPARAM_LOCAL = 8)
      real raymid(N_PARAM),rayleft0(N_PARAM),rayright0(N_PARAM),amplitudemid,phasemid,anglemid
      real rayold(N_PARAM),raynew(N_PARAM),raysect(N_PARAM),tsect,rtsect(N_RTPARAM_LOCAL)
      integer crossect, kmahsect, fictsect, kmahmid
      real acttime,time0,dtnew,dtnewhalf,chkbounddummy,spread_dummy,cos_dummy,dtnew0
      integer crossdummy,flaginterpol,flagintersect,i,n,iloop,nloops

      if (N_RTPARAM .gt. N_RTPARAM_LOCAL) goto 990

      spread_dummy = 0.0
c wichtig: cos_min muss 2.0 sein, damit auf jeden Fall in INTERPOLATE interpoliert wird!
      cos_dummy = 2.0
      flaginterpol = 0
      flagintersect = 0
      ihaltdown = 1
      ihaltup = 1

      flagleft = 0
      flagright = 0


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
     &              rayright0,
     &              N_PARAM,MAXP_XGRID,MAXP_ZGRID)
               do i = 1, N_PARAM
                  rayold(i) = rayright0(i)
               end do
            end do
         end if
         dtnewhalf = 0.5*dtnew
         call runkutta(rayold,dtnew,dtnewhalf,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        rayright0,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
 
         amplituderight = amplitude2
         phaseright = phase2
         angleright = angle2
         do i = 1, N_PARAM
            rayleft0(i) = ray1(i)
         end do
         amplitudeleft = amplitude1
         phaseleft = phase1
         angleleft = angle1
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
     &              rayleft0,
     &              N_PARAM,MAXP_XGRID,MAXP_ZGRID)
               do i = 1, N_PARAM
                  rayold(i) = rayleft0(i)
               end do
            end do
         end if
         dtnewhalf = 0.5*dtnew
         call runkutta(rayold,dtnew,dtnewhalf,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        rayleft0,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         amplitudeleft = amplitude1
         phaseleft = phase1
         angleleft = angle1
         do i = 1, N_PARAM
            rayright0(i) = ray2(i)
         end do
         amplituderight = amplitude2
         phaseright = phase2
         angleright = angle2
         time0 = time2
      end if
c
c Rueckpropagation Ende
c================================================================

      dtnew = 0.5 * dtray
      dtnewhalf = 0.5*dtnew
c
c N Interpolationen zwischen RAYLEFT und RAYRIGHT
c
      do n = 1, N_INTERPOLATIONS
         call interpolate(rayleft0,amplitudeleft,phaseleft,kmah,angleleft,
     &        rayright0,amplituderight,phaseright,angleright,
     &        spread_dummy, cos_dummy, 1,
     &        raymid,amplitudemid,phasemid,kmahmid,anglemid,flaginterpol,
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
         end do
         if (fictsect .eq. fict1) then
            do i = 1, N_PARAM
               rayleft0(i) = raymid(i)
               rayleft(i) = raysect(i)
            end do
            amplitudeleft = amplitudemid
            phaseleft = phasemid
            angleleft = anglemid
            flagleft = 1
            do i = 1, N_RTPARAM
               rtleft(i) = rtsect(i)
            end do
            timeleft = tsect
            kmahleft = kmahsect
            fictleft = fictsect
         else
            do i = 1, N_PARAM
               rayright0(i) = raymid(i)
               rayright(i) = raysect(i)
            end do
            amplituderight = amplitudemid
            phaseright = phasemid
            angleright = anglemid
            flagright = 1
            do i = 1, N_RTPARAM
               rtright(i) = rtsect(i)
            end do
            timeright = tsect
            kmahright = kmahsect
            fictright = fictsect
         end if
      end do

      goto 999

 990  write(*,*) "INTERPOL_EXTRA: PARAMETER 'N_RTPARAM_LOCAL' too low!"
      write(*,*) "Change it to at least ",N_RTPARAM
      error = 1

 999  end


c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c SUBROUTINE PROPAGATE / Hauptprogramm WFRONT
c wird aufgerufen von WFRONT
c Aufruf der Unterprogramme DO_INTERSECTION, INTERSECT_RAY, GET_Z, GET_RTSECT, RUNKUTTA2, DERY, DERY_EXTRA
c
c Ein Strahl (RAYLAST) wird um DTNEW weiterpropagiert --> RAY.
c Ausgabe: RAY, THESECT, RTSECT ...
c FLAG: =-1: Strahl ist gerade von Grenzflaeche gestartet
c FLAG: =0 : sonst
c
      subroutine propagate(raylast,actlayer,acttime,dtray,dthalf,ihaltup,ihaltdown,veltype,flag,
     &     n_int,x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,ptos,left_border,right_border,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ray,kmah,chkbound,iraylast,iray,nextraylast,nextray,
     &     thesect,tsect,rtsect,kmahsect,crossect,fictsect,
     &     leftboresect,tleftboresect,kmahlboresect,
     &     rightboresect,trightboresect,kmahrboresect,halfsteps_count,
     &     n_allint,MAXPOINTS_INT,N_PARAM,N_RTPARAM,n_boreholes,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer n_int,n_allint,MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID
      integer iray,iraylast,actlayer,veltype,ihaltdown,ihaltup,flag
      integer halfsteps_count

      integer nextraylast,nextray,crossect,kmah,kmahsect,fictsect
      real raylast(N_PARAM),ray(N_PARAM),thesect(N_PARAM),rtsect(N_RTPARAM)
      real tsect,chkbound,dtray,dthalf,acttime

      integer npoints_int(n_allint),iii(MAXPOINTS_INT,n_allint),iii2(MAXPOINTS_INT,n_int)
      real x_int(MAXPOINTS_INT,n_allint), z_int(MAXPOINTS_INT,n_allint)
      real b_int(MAXPOINTS_INT,n_allint), c_int(MAXPOINTS_INT,n_allint)
      real d_int(MAXPOINTS_INT,n_allint)
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID),v_grid(4,MAXP_ZGRID,MAXP_XGRID)
      real ptos,left_border, right_border     
      integer n_boreholes,kmahlboresect(n_boreholes),kmahrboresect(n_boreholes)
      real leftboresect(N_PARAM,n_boreholes),tleftboresect(n_boreholes)
      real rightboresect(N_PARAM,n_boreholes),trightboresect(n_boreholes)

      integer actp,flagray,newcross,ibore,hit,loop,nloops,i,idummy
      integer crossbore, actbore
      real dzdx,d2zdx2,x1out,x2out,dummy
      real ray2(N_PARAM),raylast2(N_PARAM),dtnew,dtnewhalf,epsilon6,diff1,diff2
      real vel(6),kh1(N_PARAM),maxgrad

      epsilon6 = 0.000001
      newcross = 0


c============================================================
c Propagation von RAYLAST um DTRAY auf RAY

c Voreinstellung von KH1 (Hauptsaechlich geht es um das Erhalten von VEL() zum Checken)
c
      call dery_extra(raylast,nx_grid,nz_grid,x_grid,z_grid,
     &     v_grid,ptos,veltype,actlayer,
     &     kh1,vel,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

      maxgrad = max(abs(vel(2)),abs(vel(3)))
c-------------------
c Zeitschritt muss halbiert werden... (Genauigkeit!)
c 0.06 ist dabei ein empirisch bestimmter Checkwert
c
      if (maxgrad*dtray .gt. 0.06) then
         nloops = 1
         dtnew = dtray
         dtnewhalf = dthalf
         do while (maxgrad*dtnew .gt. 0.06)
            halfsteps_count = halfsteps_count + 1
            nloops = 2 * nloops
            dtnew = dtnew * 0.5
            dtnewhalf = dtnew * 0.5
         end do
         do i = 1, N_PARAM
            raylast2(i) = raylast(i)
         end do
         do 1 loop = 1, nloops
            call runkutta2(raylast2,dtnew,dtnewhalf,kh1,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           ptos,veltype,actlayer,
     &           ray,
     &           N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            if (loop .ne. nloops) then
               do i = 1, N_PARAM
                  raylast2(i) = ray(i)
               end do
               call dery(raylast2,nx_grid,nz_grid,x_grid,z_grid,
     &              v_grid,ptos,veltype,actlayer,
     &              kh1,
     &              N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            end if
 1       continue
c-------------------
c Voller Zeitschritt...
c
      else
         call runkutta2(raylast,dtray,dthalf,kh1,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        ray,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
      end if
c
c Ende Propagation um DTRAY
c==============================================================================================
c Pruefe auf Schnittpunkt
c Im Notfall Schnittpunkt bestimmen
c
      if (crossect .eq. 0) then
         call do_intersection(raylast,actlayer,(acttime-dtray),dtray,ihaltup,ihaltdown,veltype,flag,
     &        x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,ptos,left_border,right_border,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ray,kmah,chkbound,newcross,
     &        thesect,tsect,rtsect,kmahsect,crossect,fictsect,
     &        n_int,MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
      end if
c-------------------------------------
c Strahl ueberschritt Modellrand:
c
      if (ray(1) .le. left_border .or. ray(1) .ge. right_border) then
         nextraylast = -abs(nextray)
         nextray = 0
c-------------------------------------
c Normal
c
      else
         iraylast = iray
         if (raylast(5)*ray(5) .lt. 0.0) then
            kmah = kmah + 1
         end if
      end if

c==============================================================================================
c Boreholes
c Denke dran: x_int(,) sind die z-Koordinaten, z_int(,) die x-Koordinaten!
c
      ibore = n_int
      hit = 0
      flagray = 0
      do while (ibore .lt. n_allint)
         ibore = ibore + 1
         if (crossect .eq. 0) then
            if (raylast(2) .le. x_int(npoints_int(ibore),ibore) .and.
     &           ray(2) .le. x_int(npoints_int(ibore),ibore)) then
               call get_z(raylast(2),x_int(1,ibore),z_int(1,ibore),d_int(1,ibore),c_int(1,ibore),b_int(1,ibore),npoints_int(ibore),
     &              x1out,actp,
     &              MAXPOINTS_INT)
               call get_z(ray(2),x_int(1,ibore),z_int(1,ibore),d_int(1,ibore),c_int(1,ibore),b_int(1,ibore),npoints_int(ibore),
     &              x2out,actp,
     &              MAXPOINTS_INT)
               diff1 = raylast(1) - x1out
               diff2 = ray(1) - x2out
c Falls Bohrloch geschnitten wurde, der erste Punkt des Strahls aber nicht auf dem Bohrloch liegt...
c               if ((diff1*diff2 .lt. 0.0 .and. abs(diff1) .gt. epsilon6) .or.
c     &              abs(diff2) .lt. epsilon6) then
c Neu: geschnitten oder RAY liegt genau auf Bohrloch
               if ((diff1*diff2) .lt. 0.0 .or. diff2 .eq. 0.0) then
                  do i = 1, N_PARAM
                     ray2(i) = ray(i)
                  end do
                  hit = 1
                  dtnew = dtray
               end if               
            end if
c Strahl hat gerade eine Grenzflaeche ueberschritten:
         else if (newcross .ne. 0) then
            if (raylast(2) .le. x_int(npoints_int(ibore),ibore) .and.
     &           thesect(2) .le. x_int(npoints_int(ibore),ibore)) then 
               call get_z(raylast(2),x_int(1,ibore),z_int(1,ibore),d_int(1,ibore),c_int(1,ibore),b_int(1,ibore),npoints_int(ibore),
     &              x1out,actp,
     &              MAXPOINTS_INT)
               call get_z(thesect(2),x_int(1,ibore),z_int(1,ibore),d_int(1,ibore),c_int(1,ibore),b_int(1,ibore),npoints_int(ibore),
     &              x2out,actp,
     &              MAXPOINTS_INT)
               diff1 = raylast(1) - x1out
               diff2 = thesect(1) - x2out
c     Falls Bohrloch geschnitten wurde, der erste Punkt des Strahls aber nicht auf dem Bohrloch liegt...
c               if ((diff1*diff2 .lt. 0.0 .and. abs(diff1) .gt. epsilon6) .or.
c     &              abs(diff2) .lt. epsilon6) then
c Neu:
               if ((diff1*diff2) .lt. 0.0 .or. diff2 .eq. 0.0) then
                  do i = 1, N_PARAM
                     ray2(i) = thesect(i)
                  end do
                  hit = 1
                  dtnew = tsect + dtray - acttime
               end if               
            end if
         end if

c---------------------------------
c Bohrloch getroffen:
         if (hit .eq. 1) then
            actbore = ibore - n_int
c Fall 'diff1 = 0.0' kann nicht auftreten (vorher abgefragt)
c-----------
c Strahl kommt von links:
            if (diff1 .gt. 0.0) then
c Nur neuen Schnittpunkt suchen, wenn noch keiner vorhanden (--> DAS ist eingeplante die Unzulaenglichkeit bei Bohrloechern)
               if (tleftboresect(actbore) .eq. 0.0) then
                  crossbore = -ibore
                  call intersect_ray(raylast,ray2,actp,crossbore,x2out,(acttime-dtray),
     &                 flag,dtnew,2,dummy,
     &                 x_int,d_int,c_int,b_int,z_int,npoints_int,iii,iii2,
     &                 nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &                 actlayer,ptos,veltype,flagray,
     &                 leftboresect(1,actbore),tleftboresect(actbore),idummy,dzdx,d2zdx2,
     &                 N_PARAM,n_allint,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)
                  if (raylast(5)*leftboresect(5,actbore) .lt. 0.0) then
                     kmahlboresect(actbore) = kmah + 1
                  else
                     kmahlboresect(actbore) = kmah
                  end if
c                  write(86,*) leftboresect(1,actbore),leftboresect(2,actbore)
               end if
c-----------
c Strahl kommt von rechts:
            else
               if (trightboresect(actbore) .eq. 0.0) then
                  crossbore = ibore
                  call intersect_ray(raylast,ray2,actp,crossbore,x2out,(acttime-dtray),
     &                 flag,dtnew,2,dummy,
     &                 x_int,d_int,c_int,b_int,z_int,npoints_int,iii,iii2,
     &                 nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &                 actlayer,ptos,veltype,flagray,
     &                 rightboresect(1,actbore),trightboresect(actbore),idummy,dzdx,d2zdx2,
     &                 N_PARAM,n_allint,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)
                  if (raylast(5)*rightboresect(5,actbore) .lt. 0.0) then
                     kmahrboresect(actbore) = kmah + 1
                  else
                     kmahrboresect(actbore) = kmah
                  end if
c                  write(86,*) rightboresect(1,actbore),rightboresect(2,actbore)
               end if
            end if
         
         end if
      end do
c
c Ende Schleife ueber alle Bohrloecher
c----------------------------------------

      end


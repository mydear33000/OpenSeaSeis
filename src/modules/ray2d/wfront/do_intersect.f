c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm DO_INTERSECTION / Hauptprogramm WFRONT
c wird aufgerufen von: WFRONT, PROPAGATE, GET_CORNERRAY
c Aufruf der Unterprogramme INTERSECT_RAY, GET_Z, GET_RTSECT, GET_FLAGRAY
c
c Ueberpruefung auf Schnittpunkt eines Strahls mit Grenzflaeche
c RAYOLD: Anfangswerte des Strahls
c RAY: Endwerte des Strahls
c Ausgabe: SECT, RTSECT ...
c FLAG: =-1: Strahl ist gerade von Grenzflaeche gestartet; =0: sonst.
c
      subroutine do_intersection(rayold,actlayer,dtime,dtray,ihaltup,ihaltdown,veltype,flag,
     &     x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,ptos,left_border,right_border,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ray,kmah,chkbound,newcross,
     &     sect,tsect,rtsect,kmahsect,crossect,fictsect,
     &     MAX_INT,MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer MAX_INT,MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID
      integer actlayer,veltype,ihaltdown,ihaltup,flag

      integer crossect,kmah,kmahsect,fictsect
      real rayold(N_PARAM),ray(N_PARAM),sect(N_PARAM),rtsect(N_RTPARAM)
      real tsect,chkbound,dtray,dtime

      integer npoints_int(MAX_INT),iii(MAXPOINTS_INT,MAX_INT),iii2(MAXPOINTS_INT,MAX_INT)
      real x_int(MAXPOINTS_INT,MAX_INT), z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT), c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT)
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID),v_grid(4,MAXP_ZGRID,MAXP_XGRID)
      real ptos,left_border, right_border     

      integer lastactp,firstactp,actp,flagray,newcross,i,istep,iactp,actptmp
      real zint,dzdx,d2zdx2
      real dtnew

      flagray = 0

c------------------------
c Pruefe auf Schnittpunkt
c
      call get_z(ray(1),x_int(1,actlayer),z_int(1,actlayer),d_int(1,actlayer),c_int(1,actlayer),
     &     b_int(1,actlayer),npoints_int(actlayer),
     &     zint,lastactp,
     &     MAXPOINTS_INT)
c--------------
c Strahl hat die obere Grenzflaeche ueberschritten:
      if (ray(2) .le. zint) then
         call get_actp(rayold(1),actlayer,x_int,npoints_int,
     &        firstactp,
     &        MAX_INT,MAXPOINTS_INT)
c--------------
c mind. ein Stuetzpunkt liegt zwischen ray und rayold!
c Ueberpruefung, ob Strahl in Ecke laeuft und dabei evt. die andere Grenzflaeche zuerst schneidet:
         if (firstactp .ne. lastactp) then
            if (ray(1) .lt. rayold(1)) then
               firstactp = firstactp - 1
               istep = -1
            else
               firstactp = firstactp + 1
               istep = 1
            end if
            do iactp = firstactp, lastactp, istep
               do i = 1, npoints_int(actlayer+1)                  
                  if (iii(i,actlayer+1) .eq. iactp) then
                     if (ray(1) .lt. rayold(1)) then
                        actptmp = iactp + 1
                     else
                        actptmp = iactp
                     end if
                     flagray = -1
                     call get_flagray(x_int(actptmp,actlayer),z_int(actptmp,actlayer),rayold,dtray,
     &                    nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &                    flagray,dtnew,
     &                    MAXP_XGRID,MAXP_ZGRID,N_PARAM)
                     if (flagray .gt. 0) then
                        call get_z(ray(1),x_int(1,actlayer+1),z_int(1,actlayer+1),d_int(1,actlayer+1),
     &                       c_int(1,actlayer+1),b_int(1,actlayer+1),npoints_int(actlayer+1),
     &                       zint,lastactp,
     &                       MAXPOINTS_INT)
                        goto 20
                     end if
                     goto 10
                  end if
               end do
            end do
         end if
         goto 10
c---------------------
c Pruefe auf Schnittpunkt mit unterer Grenzflaeche
      else
         call get_z(ray(1),x_int(1,actlayer+1),z_int(1,actlayer+1),d_int(1,actlayer+1),
     &        c_int(1,actlayer+1),b_int(1,actlayer+1),npoints_int(actlayer+1),
     &        zint,lastactp,
     &        MAXPOINTS_INT)
c--------------
c Strahl hat die untere Grenzflaeche ueberschritten:
         if (ray(2) .ge. zint) then
            call get_actp(rayold(1),actlayer+1,x_int,npoints_int,
     &           firstactp,
     &           MAX_INT,MAXPOINTS_INT)
c--------------
c mind. ein Stuetzpunkt liegt zwischen ray und rayold!
c Ueberpruefung, ob Strahl in Ecke laeuft und dabei evt. die andere Grenzflaeche zuerst schneidet:
            if (firstactp .ne. lastactp) then
               if (ray(1) .lt. rayold(1)) then
                  firstactp = firstactp - 1
                  istep = -1
               else
                  firstactp = firstactp + 1
                  istep = 1
               end if
               do iactp = firstactp, lastactp, istep
                  if (iii(iactp,actlayer+1) .gt. 0) then
                     if (ray(1) .lt. rayold(1)) then
                        actptmp = iii(iactp,actlayer+1) + 1
                     else
                        actptmp = iii(iactp,actlayer+1)
                     end if
                     flagray = 1
                     call get_flagray(x_int(actptmp,actlayer),z_int(actptmp,actlayer),rayold,dtray,
     &                    nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &                    flagray,dtnew,
     &                    MAXP_XGRID,MAXP_ZGRID,N_PARAM)
                     if (flagray .lt. 0) then
                        call get_z(ray(1),x_int(1,actlayer),z_int(1,actlayer),d_int(1,actlayer),
     &                       c_int(1,actlayer),b_int(1,actlayer),npoints_int(actlayer),
     &                       zint,lastactp,
     &                       MAXPOINTS_INT)
                        goto 10
                     end if
                     goto 20
                  end if
               end do
            end if
            goto 20
         end if
      end if
      goto 30

c--------------------------------------------------
c Strahl hat die obere Grenzflaeche ueberschritten:
c
 10   continue
      crossect = -actlayer
      call intersect_ray(rayold,ray,lastactp,crossect,zint,dtime,
     &     flag,dtray,1,dtnew,
     &     x_int,d_int,c_int,b_int,z_int,npoints_int,iii,iii2,
     &     nx_grid,nz_grid,x_grid,z_grid,
     &     v_grid,
     &     actlayer,ptos,veltype,flagray,
     &     sect,tsect,fictsect,dzdx,d2zdx2,
     &     N_PARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)

c Ueberpruefung auf Konsistenz bei corner ray:
c im Notfall wird der Schnittpunkt des corner rays in die Naehe des Eckpunktes gesetzt
      if (flagray .ne. 0) then
         call get_actp(sect(1),actlayer,x_int,npoints_int,
     &        actp,
     &        MAX_INT,MAXPOINTS_INT)
c Strahl laeuft nach links:
         if (ray(1) .lt. rayold(1) .and. actp .lt. actptmp) then
            sect(1) = x_int(actptmp,actlayer) + 0.000001
            call get_z(sect(1),x_int(1,actlayer),z_int(1,actlayer),d_int(1,actlayer),c_int(1,actlayer),
     &           b_int(1,actlayer),npoints_int(actlayer),
     &           sect(2),actp,
     &           MAXPOINTS_INT)
            fictsect = actp
c            sect(2) = z_int(actptmp,actlayer)
c Strahl laeuft nach rechts:
         else if (ray(1) .ge. rayold(1) .and. actp .ge. actptmp) then
            if (actp .eq. actptmp .and. abs(x_int(actp,actlayer)-sect(1)) .lt. 0.000001) then
            else
               sect(1) = x_int(actptmp,actlayer) - 0.000001
               call get_z(sect(1),x_int(1,actlayer),z_int(1,actlayer),d_int(1,actlayer),c_int(1,actlayer),
     &              b_int(1,actlayer),npoints_int(actlayer),
     &              sect(2),actp,
     &              MAXPOINTS_INT)
               fictsect = actp
c               sect(2) = z_int(actptmp,actlayer)
            end if
         end if
      end if

      if (sect(1) .le. left_border .or. sect(1) .ge. right_border) then
         crossect = 0
      else
         chkbound = dzdx*sect(3) - sect(4)
c Aufgrund Rechenungenauigkeit muss hier gepfuscht werden:
         if (abs(chkbound) .lt. 0.005) then
            chkbound = 0.1
         end if
         if (rayold(5)*sect(5) .lt. 0.0) then
            kmahsect = kmah + 1
         else
            kmahsect = kmah
         end if
         newcross = -1
         if (ihaltup .ne. 0) then
            call get_rtsect(sect,dzdx,d2zdx2,
     &           nx_grid,nz_grid,x_grid,z_grid,
     &           v_grid,ptos,veltype,actlayer,
     &           rtsect,
     &           N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
         end if
      end if
      goto 30

c---------------------------------------------------
c Strahl hat die untere Grenzflaeche ueberschritten:
c
 20   continue
      crossect = actlayer+1
      call intersect_ray(rayold,ray,lastactp,crossect,zint,dtime,
     &     flag,dtray,1,dtnew,
     &     x_int,d_int,c_int,b_int,z_int,npoints_int,iii,iii2,
     &     nx_grid,nz_grid,x_grid,z_grid,
     &     v_grid,
     &     actlayer,ptos,veltype,flagray,
     &     sect,tsect,fictsect,dzdx,d2zdx2,
     &     N_PARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)

c Ueberpruefung auf Konsistenz bei corner ray:
c im Notfall wird der Schnittpunkt des corner rays in die Naehe des Eckpunktes gesetzt
      if (flagray .ne. 0) then
         call get_actp(sect(1),actlayer,x_int,npoints_int,
     &        actp,
     &        MAX_INT,MAXPOINTS_INT)
c Strahl laeuft nach links:
         if (ray(1) .lt. rayold(1) .and. actp .lt. actptmp) then
            sect(1) = x_int(actptmp,actlayer) + 0.000001
            call get_z(sect(1),x_int(1,crossect),z_int(1,crossect),d_int(1,crossect),c_int(1,crossect),
     &           b_int(1,crossect),npoints_int(crossect),
     &           sect(2),actp,
     &           MAXPOINTS_INT)
            fictsect = actp
c            sect(2) = z_int(actptmp,actlayer)
c Strahl laeuft nach rechts:
         else if (ray(1) .ge. rayold(1) .and. actp .ge. actptmp) then
            if (actp .eq. actptmp .and. abs(x_int(actp,actlayer)-sect(1)) .lt. 0.000001) then
            else
               sect(1) = x_int(actptmp,actlayer) - 0.000001
               call get_z(sect(1),x_int(1,crossect),z_int(1,crossect),d_int(1,crossect),c_int(1,crossect),
     &              b_int(1,crossect),npoints_int(crossect),
     &              sect(2),actp,
     &              MAXPOINTS_INT)
               fictsect = actp
c               sect(2) = z_int(actptmp,actlayer)
            end if
         end if
      end if

      if (sect(1) .le. left_border .or. sect(1) .ge. right_border) then
         crossect = 0
      else
         chkbound = dzdx*sect(3) - sect(4)
c Aufgrund Rechenungenauigkeit muss hier gepfuscht werden:
         if (abs(chkbound) .lt. 0.005) then
            chkbound = -0.1
         end if
         if (rayold(5)*sect(5) .lt. 0.0) then
            kmahsect = kmah + 1
         else
            kmahsect = kmah
         end if
         newcross = 1
         if (ihaltdown .ne. 0) then
            call get_rtsect(sect,dzdx,d2zdx2,
     &           nx_grid,nz_grid,x_grid,z_grid,
     &           v_grid,ptos,veltype,actlayer,
     &           rtsect,
     &           N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
         end if
      end if

 30   continue

      end


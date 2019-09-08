c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm SORTSECT / Hauptprogramm WFRONT
c Wird aufgerufen von WFRONT
c Aufruf der Unterprogramme GET_CORNERRAY,INTERPOL_EXTRA,GET_RAYONINT,GET_Z2,GET_RTSECT,SET_FICT,FIND_EDGE
c
c - loescht Strahlen, die den Modellrand ueberschritten haben
c - loescht Strahlen, die eine Grenzflaeche ueberschritten haben, die nicht mehr gebraucht wird
c - Interpoliert neue Strahlen an folgenden Stellen:
c   a) in Ecken von Blockstrukturen
c   b) zwischen unterschiedlichen Teilstuecken einer Grenzflaeche
c   c) wenn zwischen zwei Schnittpunkten mehrere Stuetzpunkte liegen!
c      Das darf nicht so bleiben, weil man so nicht genau weiss, auf welchen Teilen man sich befindet
c Punkt d) ist nicht mehr noetig, seit GET_RAYONINT verbessert wurde:
c   d) bei zu grossem Abstand zur Erhoehung der Genauigkeit, bevor paraxial approximiert wird
c
      subroutine sortsect(sect,tsect,rtsect,amplitude,phase,angle,kmahsect,crossect,ihaltup,ihaltdown,
     &     dtray,actlayer,left_border,right_border,nray,spread_max2,INF_DIST,
     &     x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,f_out,f_xrays,flag_smooth,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,veltype,ptos,
     &     nextsect,fictsect,error,
     &     N_PARAM,N_RTPARAM,MAX_RAYS,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer N_PARAM,N_RTPARAM,MAX_RAYS,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID
      real sect(N_PARAM,MAX_RAYS),tsect(MAX_RAYS),amplitude(MAX_RAYS),phase(MAX_RAYS),angle(MAX_RAYS)
      integer kmahsect(MAX_RAYS),nextsect(0:MAX_RAYS),crossect(MAX_RAYS),fictsect(MAX_RAYS)
      real rtsect(N_RTPARAM,MAX_RAYS),dtray,left_border,right_border,spread_max2
      integer actlayer,veltype,nray,f_out,flag_smooth,ihaltup,ihaltdown
      real x_int(MAXPOINTS_INT,MAX_INT),z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT),c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT),ptos
      integer npoints_int(MAX_INT),iii(MAXPOINTS_INT,MAX_INT),iii2(MAXPOINTS_INT,MAX_INT),nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID),v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      integer error, f_xrays

      integer flag,flagup,flagdown,iplus,iplus2,actptmp,iray0
      integer irayold,iray,next,actp0,actp1,actp2,actint,found,i,actp3,int3
      integer n_interpolations,n
      real part, delta_part
      real dx,x0,z0,dzdx,d2zdx2,dzdx1,dzdx2,d2zdx21,d2zdx22,epsilon,INF_DIST,distance

      epsilon = 0.0001

c==============================================================================
c Schleife ueber alle Strahlen
c
      iray = 0
      irayold = 0
      next = abs(nextsect(iray))
      do while (next .ne. 0)
         iray = next
         next = abs(nextsect(iray))
c===================================================
c Ebene (1)
c Strahl ueberschritt Modellrand und wird geloescht:
c
         if (crossect(iray) .eq. 0) then
            nextsect(irayold) = -next
c===================================================
c Ebene (1)
c Strahl liegt innerhalb eines coinciding Grenzflaechenpakets, darf nicht angefasst werden:
c
         else if (abs(fictsect(iray)) .gt. 999) then
            irayold = iray
c===================================================
c Ebene (1)
c Strahlen IRAY und NEXT sind benachbart und haben gleiche KMAH-Indizes
c Fall 'FICTSECT(next) .gt. 999' kann nur auftreten, wenn Modell geglaettet wurde (FLAG_SMOOTH .ne. 0)
c
         else if (nextsect(iray) .gt. 0 .and. kmahsect(iray) .eq. kmahsect(next) .and. fictsect(next) .lt. 1000) then
c===========================================
c Ebene (2)
c Strahl IRAY traf andere Grenzflaeche als Strahl NEXT
c
            if (crossect(iray) .ne. crossect(next)) then
c---------------------------------------------------
c Ebene (3)
c Strahl NEXT lief aus dem Modell... STrahl IRAY wird eingebaut:
c
               if (crossect(next) .eq. 0) then
                  actp1 = fictsect(iray)
                  actint = abs(crossect(iray))
                  call set_fict(iii(actp1,actint),iii2(actp1,actint),fictsect(iray))
                  irayold = iray
c---------------------------------------------------
c Ebene (3)
c Welle laeuft in Ecke einer Blockstruktur: 2 Strahlen interpolieren:
c
               else
                  actint = abs(crossect(iray))
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Interpolation zweier Strahlen in Ecken von Blockstrukturen:
c----------------------------
c Ebene (4)
c Erster Fall: IRAY schneidet untere Grenzflaeche:
c Aendern: fuer den Fall, dass mehrere Stuetzpunkte dazwischen liegen!
                  if (crossect(iray) .gt. 0) then
                     iplus = nray + 1
                     iplus2 = nray + 2
                     if (iplus2 .gt. MAX_RAYS) goto 990
                     call get_cornerrays(sect(1,next),tsect(next),amplitude(next),phase(next),angle(next),kmahsect(next),
     &                    sect(1,iray),tsect(iray),amplitude(iray),phase(iray),angle(iray),
     &                    dtray,actlayer,left_border,right_border,f_out,
     &                    x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,
     &                    nx_grid,nz_grid,x_grid,z_grid,v_grid,veltype,ptos,
     &                    sect(1,iplus2),tsect(iplus2),rtsect(1,iplus2),amplitude(iplus2),phase(iplus2),angle(iplus2),
     &                    kmahsect(iplus2),fictsect(iplus2),flagup,
     &                    sect(1,iplus),tsect(iplus),rtsect(1,iplus),amplitude(iplus),phase(iplus),angle(iplus),
     &                    kmahsect(iplus),fictsect(iplus),flagdown,error,
     &                    N_PARAM,N_RTPARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)
                     if (error .ne. 0) goto 999
c Neuer Strahl, der untere Grenzflaeche traf:
                     if (flagdown .ne. 0) then
                        if (f_xrays .ne. 0) write(f_xrays,*) sect(1,iplus),sect(2,iplus),' corner ray'
                        crossect(iplus) = crossect(iray)
                        nextsect(iray) = iplus
                        nextsect(iplus) = -next
                        nray = iplus
                     end if
c Neuer Strahl, der obere Grenzflaeche traf:
                     if (flagup .ne. 0) then
                        if (f_xrays .ne. 0) write(f_xrays,*) sect(1,iplus2),sect(2,iplus2),' corner ray'
                        if (flagdown .eq. 0) then
                           do i = 1, N_PARAM
                              sect(i,iplus) = sect(i,iplus2)
                           end do
                           do i = 1, N_RTPARAM
                              rtsect(i,iplus) = rtsect(i,iplus2)
                           end do
                           tsect(iplus) = tsect(iplus2)
                           amplitude(iplus) = amplitude(iplus2)
                           phase(iplus) = phase(iplus2)
                           kmahsect(iplus) = kmahsect(iplus2)
                           fictsect(iplus) = fictsect(iplus2)
                           angle(iplus) = angle(iplus2)
                           crossect(iplus) = crossect(next)
                           nextsect(iray) = -iplus
                           nextsect(iplus) = next
                           nray = iplus
                        else
                           crossect(iplus2) = crossect(next)
                           nextsect(iplus) = -iplus2
                           nextsect(iplus2) = next
                           nray = iplus2
                        end if
                     end if
c----------------------------
c Ebene (4)
c Zweiter Fall: IRAY schneidet obere Grenzflaeche:
                  else
                     iplus = nray + 1
                     iplus2 = nray + 2
                     if (iplus2 .gt. MAX_RAYS) goto 990
                     call get_cornerrays(sect(1,iray),tsect(iray),amplitude(iray),phase(iray),angle(iray),kmahsect(iray),
     &                    sect(1,next),tsect(next),amplitude(next),phase(next),angle(next),
     &                    dtray,actlayer,left_border,right_border,f_out,
     &                    x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,
     &                    nx_grid,nz_grid,x_grid,z_grid,v_grid,veltype,ptos,
     &                    sect(1,iplus),tsect(iplus),rtsect(1,iplus),amplitude(iplus),phase(iplus),angle(iplus),
     &                    kmahsect(iplus),fictsect(iplus),flagup,
     &                    sect(1,iplus2),tsect(iplus2),rtsect(1,iplus2),amplitude(iplus2),phase(iplus2),angle(iplus2),
     &                    kmahsect(iplus2),fictsect(iplus2),flagdown,error,
     &                    N_PARAM,N_RTPARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)
                     if (error .ne. 0) goto 999
c Neuer Strahl, der obere Grenzflaeche traf:
                     if (flagup .ne. 0) then
                        if (f_xrays .ne. 0) write(f_xrays,*) sect(1,iplus),sect(2,iplus),' corner ray'
                        crossect(iplus) = crossect(iray)
                        nextsect(iray) = iplus
                        nextsect(iplus) = -next
                        nray = iplus
                     end if
c Neuer Strahl, der untere Grenzflaeche traf:
                     if (flagdown .ne. 0) then
                        if (f_xrays .ne. 0) write(f_xrays,*) sect(1,iplus2),sect(2,iplus2),' corner ray'
                        if (flagup .eq. 0) then
                           do i = 1, N_PARAM
                              sect(i,iplus) = sect(i,iplus2)
                           end do
                           do i = 1, N_RTPARAM
                              rtsect(i,iplus) = rtsect(i,iplus2)
                           end do
                           tsect(iplus) = tsect(iplus2)
                           amplitude(iplus) = amplitude(iplus2)
                           phase(iplus) = phase(iplus2)
                           kmahsect(iplus) = kmahsect(iplus2)
                           fictsect(iplus) = fictsect(iplus2)
                           angle(iplus) = angle(iplus2)
                           crossect(iplus) = crossect(next)
                           nextsect(iray) = -iplus
                           nextsect(iplus) = next
                           nray = iplus
                        else
                           crossect(iplus2) = crossect(next)
                           nextsect(iplus) = -iplus2
                           nextsect(iplus2) = next
                           nray = iplus2
                        end if
                     end if
                  end if
c------------
c Ebene (4)
c keine einziger neuer Strahl:
                  if (flagup .eq. 0 .and. flagdown .eq. 0) then                     
                     irayold = iray
                  else
c Noch mal denselben Strahl durchlaufen lassen, damit fictsect gesetzt werden kann, bei Stuetzpunkten
c dazwischen interpoliert werden kann usw...:
                     next = iray
                  end if
c
c Ende Ebene (4)
c Ende Interpolation von Strahlen in Ecke einer Blockstruktur (corner rays)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               end if
c
c Ende Ebene (3)
c-------------------------------------------

c============================================
c Ebene (2)
c Schnittpunkt wird nicht weiter gebraucht, wird also geloescht:
c
            else if ((crossect(iray) .gt. 0 .and. ihaltdown .eq. 0) .or.
     &              (crossect(iray) .lt. 0 .and. ihaltup .eq. 0)) then
               nextsect(irayold) = -next
c=============================================================================================================
c Ebene (2)
c (mindestens ein) Stuetzpunkt liegt zwischen Strahl IRAY und NEXT:
c     
            else if (fictsect(iray) .ne. fictsect(next)) then

               actp1 = fictsect(iray)
               actp2 = fictsect(next)
               actint = abs(crossect(iray))

c ...wird nach einer Verbesserung von GET_RAYONINT nicht mehr gebraucht:
c               distance = (sect(1,iray)-sect(1,next))**2 + (sect(2,iray)-sect(2,next))**2
               distance = 0.0
c---------------------------------------------------
c NICHT MEHR: Extra Interpolation, falls Abstand zu gross ist:
c ...aber wenn mehr als ein Punkt dazwischen ist
c               if (distance .gt. spread_max2 .or. abs(actp1-actp2) .gt. 1) then
               if (abs(actp1-actp2) .gt. 1) then
                  iplus = nray + 1
                  iplus2 = nray + 2
                  if (iplus2 .gt. MAX_RAYS) goto 990
                  call interpol_extra(sect(1,iray),tsect(iray),amplitude(iray),phase(iray),angle(iray),kmahsect(iray),
     &                 fictsect(iray),sect(1,next),tsect(next),amplitude(next),phase(next),angle(next),
     &                 dtray,actlayer,left_border,right_border,
     &                 x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,
     &                 nx_grid,nz_grid,x_grid,z_grid,v_grid,veltype,ptos,
     &                 sect(1,iplus),tsect(iplus),rtsect(1,iplus),amplitude(iplus),phase(iplus),angle(iplus),
     &                 kmahsect(iplus),fictsect(iplus),flagup,
     &                 sect(1,iplus2),tsect(iplus2),rtsect(1,iplus2),amplitude(iplus2),phase(iplus2),angle(iplus2),
     &                 kmahsect(iplus2),fictsect(iplus2),flagdown,error,
     &                 N_PARAM,N_RTPARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)
                  if (error .ne. 0) goto 999

                  if (flagup .ne. 0) then
                     if (f_xrays .ne. 0) write(f_xrays,*) sect(1,iplus),sect(2,iplus),' extra interpol'
                     actp1 = fictsect(iray)
                     call set_fict(iii(actp1,actint),iii2(actp1,actint),fictsect(iray))
                     crossect(iplus) = crossect(iray)
                     nextsect(iray) = iplus
                     nextsect(iplus) = next
                     nray = iplus
                     irayold = iray
                     iray = iplus
                  end if
                  if (flagdown .ne. 0) then
                     if (f_xrays .ne. 0) write(f_xrays,*) sect(1,iplus2),sect(2,iplus2),' extra interpol'
                     if (flagup .eq. 0) then
                        do i = 1, N_PARAM
                           sect(i,iplus) = sect(i,iplus2)
                        end do
                        do i = 1, N_RTPARAM
                           rtsect(i,iplus) = rtsect(i,iplus2)
                        end do
                        tsect(iplus) = tsect(iplus2)
                        amplitude(iplus) = amplitude(iplus2)
                        phase(iplus) = phase(iplus2)
                        kmahsect(iplus) = kmahsect(iplus2)
                        fictsect(iplus) = fictsect(iplus2)
                        angle(iplus) = angle(iplus2)
                        crossect(iplus) = crossect(next)
                        nextsect(iray) = iplus
                        nextsect(iplus) = next
                        nray = iplus
                        next = iplus
                     else
                        crossect(iplus2) = crossect(next)
                        nextsect(iray) = iplus2
                        nextsect(iplus2) = next
                        nray = iplus2
                        next = iplus2
                     end if
                     actp2 = fictsect(next)
                  end if
                  distance = (sect(1,iray)-sect(1,next))**2 + (sect(2,iray)-sect(2,next))**2
               end if

c-----------------------------------------
c Ist der Abstand noch immer (bzw: ist er j e t z t) zu gross,
c handelt es sich ziemlich sicher um einen nicht erkannten boundary ray
c (nicht erkannt aufgrund Aufloesungsschwierigkeiten) -> kein neuer Strahl wird interpoliert
c Ausserdem: Wenn noch immer mehrere Stuetzpunkte dazwischen liegen: wie oben...
c
               if (distance .gt. spread_max2 .or. abs(actp1-actp2) .gt. 1) then
                  call set_fict(iii(actp1,actint),iii2(actp1,actint),fictsect(iray))
                  nextsect(iray) = -next
               else
                  if (actp2 .gt. actp1) then
                     actp0 = actp2
                  else if (actp1 .gt. actp2) then
                     actp0 = actp1
                  else
                     write(*,*) "SORTSECT: PROGRAM ERROR!"
                     error = 1
                     goto 999
                  end if
c---------------------------------------------------------------------
c Suche: FOUND=0: Keine Interpolation, !=0: Interpolation muss erfolgen.
c Eckpunkt dazwischen:
                  call find_edge(actp0,actint,crossect(iray),iii,iii2,found,
     &                 MAX_INT,MAXPOINTS_INT)

                  call set_fict(iii(actp1,actint),iii2(actp1,actint),fictsect(iray))
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Veschiedene Interface-Stuecke: 2 Strahlen interpolieren:
c
                  if (found .ne. 0) then
                     iray0 = iray
c--------------------------------
c     interpoliere ersten Strahl:
c
                     dx = sect(1,iray)-x_int(actp0,actint)
                     if (abs(dx) .gt. INF_DIST) then
                        iplus = nray + 1
                        if (iplus .gt. MAX_RAYS) goto 990
                        sect(1,iplus) = x_int(actp0,actint) + sign(1.0,dx)*INF_DIST
                        actptmp = actp0
                        
                        call get_z2(sect(1,iplus),x_int(1,actint),z_int(1,actint),
     &                       d_int(1,actint),c_int(1,actint),b_int(1,actint),npoints_int(actint),
     &                       sect(2,iplus),dzdx,d2zdx2,actptmp,
     &                       MAXPOINTS_INT)
                        call get_rayonint(sect(1,next),amplitude(next),phase(next),tsect(next),angle(next),
     &                       sect(1,iray),amplitude(iray),phase(iray),tsect(iray),angle(iray),
     &                       tsect(iplus),veltype,ptos,actlayer,dtray,
     &                       nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &                       sect(1,iplus),amplitude(iplus),phase(iplus),angle(iplus),flag,
     &                       N_PARAM,MAXP_XGRID,MAXP_ZGRID)
                        if (flag .eq. 0) then
                           call get_rtsect(sect(1,iplus),dzdx,d2zdx2,
     &                          nx_grid,nz_grid,x_grid,z_grid,
     &                          v_grid,ptos,veltype,actlayer,
     &                          rtsect(1,iplus),
     &                          N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
                           if (f_xrays .ne. 0) write(f_xrays,*) sect(1,iplus),sect(2,iplus)," different parts"
                           kmahsect(iplus) = kmahsect(iray)
                           fictsect(iplus) = fictsect(iray)
                           crossect(iplus) = crossect(iray)
                           nextsect(iray) = iplus
c     nextsect(iplus) = next
c     next = iplus
                           nray = iplus
                           iray = iplus
                        end if
                     end if
c--------------------------------
c     interpoliere zweiten Strahl:
c     
                     dx = sect(1,next)-x_int(actp0,actint)
                     if (abs(dx) .gt. INF_DIST) then
                        iplus = nray + 1
                        if (iplus .gt. MAX_RAYS) goto 990
                        sect(1,iplus) = x_int(actp0,actint) + sign(1.0,dx)*INF_DIST
                        actptmp = actp0
                        call get_z2(sect(1,iplus),x_int(1,actint),z_int(1,actint),
     &                       d_int(1,actint),c_int(1,actint),b_int(1,actint),npoints_int(actint),
     &                       sect(2,iplus),dzdx,d2zdx2,actptmp,
     &                       MAXPOINTS_INT)
                        call get_rayonint(sect(1,next),amplitude(next),phase(next),tsect(next),angle(next),
     &                       sect(1,iray0),amplitude(iray0),phase(iray0),tsect(iray0),angle(iray0),
     &                       tsect(iplus),veltype,ptos,actlayer,dtray,
     &                       nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &                       sect(1,iplus),amplitude(iplus),phase(iplus),angle(iplus),flag,
     &                       N_PARAM,MAXP_XGRID,MAXP_ZGRID)
                        if (flag .eq. 0) then
                           call get_rtsect(sect(1,iplus),dzdx,d2zdx2,
     &                          nx_grid,nz_grid,x_grid,z_grid,
     &                          v_grid,ptos,veltype,actlayer,
     &                          rtsect(1,iplus),
     &                          N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
                           if (f_xrays .ne. 0) write(f_xrays,*) sect(1,iplus),sect(2,iplus)," different parts"
                           kmahsect(iplus) = kmahsect(next)
                           fictsect(iplus) = actp2
                           crossect(iplus) = crossect(next)
                           nextsect(iplus) = next
                           next = iplus
                           nray = iplus
                        end if
                     end if

c--------------------------------------------------
c Modell wird 'geglaettet' bei FLAG_SMOOTH != 0, d.h. die Wellenfront wird nicht getrennt.
c Hier wird auch noch geprueft, ob diese Stelle nicht extra auf 'Ecke' geschaltet worden ist, indem
c III in der Eingabedatei explicit auf -1 gesetzt wurde
c
                     if (flag_smooth .ne. 0) then
                        actp3 = actp0
                        int3 = actint
                        if (crossect(iray) .lt. 0) then
                           do i = 1, iii2(actp0,actint)-1
                              actp3 = iii(actp3,int3)
                              int3 = int3 - 1
                           end do
                        end if
                        if (iii(actp3,int3) .eq. -1) then
c ...wird dennoch getrennt
                           nextsect(iray) = -next
                        else
c ...klebt zusammen
                           nextsect(iray) = next
                        end if
                     else
                        nextsect(iray) = -next
c--------------------------------------------------
c Im Falle, dass die Grenzflaeche am Punkt ACTP0 glatt ist (und nicht fictious ist),
c wird FICTSECT += 10 gesetzt und NEXTSECT positiv belassen.
c Bei einer Transmission wird dann NEXTSECT = -NEXTSECT gesetzt. (siehe REFLTRANS)
c
                        x0 = x_int(actp0,actint)
                        actptmp = actp0
                        call get_z2(x0,x_int(1,actint),z_int(1,actint),d_int(1,actint),c_int(1,actint),
     &                       b_int(1,actint),npoints_int(actint),
     &                       z0,dzdx1,d2zdx21,actptmp,
     &                       MAXPOINTS_INT)
                        actptmp = actp0 - 1
                        dx = x0 - x_int(actptmp,actint)
                        dzdx2 = (3.*d_int(actptmp,actint)*dx + 2.*c_int(actptmp,actint))*dx + b_int(actptmp,actint)
                        d2zdx22 = 6.*d_int(actptmp,actint)*dx + 2.*c_int(actptmp,actint)
                        if (abs(dzdx1-dzdx2) .lt. epsilon) then
                           if (abs(d2zdx21-d2zdx22) .lt. epsilon) then
                              fictsect(iray) = fictsect(iray) + 10
                              nextsect(iray) = next
                           end if
                        end if
                     end if
                  end if
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               end if
               irayold = iray
c==============================================================
c Ebene (2)
c Normal: Strahl IRAY ist abgecheckt und darf normal eingebunden werden
c
            else

               actp1 = fictsect(iray)
               actint = abs(crossect(iray))
               call set_fict(iii(actp1,actint),iii2(actp1,actint),fictsect(iray))
               distance = (sect(1,iray)-sect(1,next))**2 + (sect(2,iray)-sect(2,next))**2
c---------------------------------------------------
c Extra Interpolation(en), falls Abstand zu gross ist:
               if (distance .gt. spread_max2) then
                  n_interpolations = int(0.3*distance/spread_max2)
                  delta_part = 1.0 / float(n_interpolations+1)
                  actptmp = actp1
                  iray0 = iray

c========================================================================================
c 'n_interpolations' Interpolationen, alle zwischen Originalstrahlen 'IRAY' und 'NEXT':
c
                  do 10 n = 1, n_interpolations
                     part = delta_part*float(n)
                     iplus = nray + 1

                     if (iplus .gt. MAX_RAYS) goto 990
                     sect(1,iplus) = sect(1,iray0) + part * (sect(1,next) - sect(1,iray0))
                     call get_z2(sect(1,iplus),x_int(1,actint),z_int(1,actint),
     &                    d_int(1,actint),c_int(1,actint),b_int(1,actint),npoints_int(actint),
     &                    sect(2,iplus),dzdx,d2zdx2,actptmp,
     &                    MAXPOINTS_INT)
                     call get_rayonint(sect(1,next),amplitude(next),phase(next),tsect(next),angle(next),
     &                    sect(1,iray0),amplitude(iray0),phase(iray0),tsect(iray0),angle(iray0),
     &                    tsect(iplus),veltype,ptos,actlayer,dtray,
     &                    nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &                    sect(1,iplus),amplitude(iplus),phase(iplus),angle(iplus),flag,
     &                    N_PARAM,MAXP_XGRID,MAXP_ZGRID)
                     if (flag .eq. 0) then
                        call get_rtsect(sect(1,iplus),dzdx,d2zdx2,
     &                       nx_grid,nz_grid,x_grid,z_grid,
     &                       v_grid,ptos,veltype,actlayer,
     &                       rtsect(1,iplus),
     &                       N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
                        if (f_xrays .ne. 0) write(f_xrays,*) sect(1,iplus),sect(2,iplus)," distance interpolation"
                        kmahsect(iplus) = kmahsect(iray0)
                        fictsect(iplus) = fictsect(iray0)
                        crossect(iplus) = crossect(iray0)
                        nextsect(iray) = iplus

                        nextsect(iplus) = next

                        nray = iplus
                        iray = iplus
                     end if
 10               continue
               end if

c               actp1 = fictsect(iray)
c               actint = abs(crossect(iray))
c               call set_fict(iii(actp1,actint),iii2(actp1,actint),fictsect(iray))
               irayold = iray
            end if
c
c Ende Ebene (2)
c==============================================================
c Ebene (1)
c Keine Interpolationen etc. zwischen IRAY und NEXT moeglich,
c Strahl IRAY wird normal eingebunden:
c
         else
            actp1 = fictsect(iray)
            actint = abs(crossect(iray))
            call set_fict(iii(actp1,actint),iii2(actp1,actint),fictsect(iray))
            irayold = iray
         end if
c
c Ende Ebene (1)
c==============================================================
      end do
c
c Ende Schleife ueber alle Strahlen
c====================================================================================

      goto 999
 990  write(*,*) "SORTSECT: Too many rays!"
      write(*,*) "Enlarge PARAMETER 'MAX_RAYS' in WFRONT"
      error = 1

 999  end

c************************************************************************************************************
c Unterprogramm SET_FICT
c Setzt FICTNEW
c
      subroutine set_fict(index1,index2,fictnew)

      implicit none

      integer index1,index2,fictnew

c----------------------------
c coinciding  u n d  fictious
      if (index2 .lt. 0) then
         fictnew = 1000*index2
c----------------------------
c nur fictious
      else if (index1 .eq. -2) then
         fictnew = -1000
c----------------------------
c normal
      else if (index2 .eq. 0) then
         fictnew = 0
c----------------------------
c nur coinciding
      else
         fictnew = 1000*index2 + 1
      end if

      end

c************************************************************************************************************
c Unterprogramm FIND_EDGE
c Wird aufgerufen von SORTSECT, INTERSECT_RAY, INTERSECT_WFRONT
c Ueberprueft, ob der Stuetzpunkt ACTP0 ein Eckpunkt ist (die Grenzschicht dort nicht glatt ist)
c Suche: FOUND=0: Keine Interpolation, =1: Interpolation muss erfolgen.
c
      subroutine find_edge(actp0,actint,crossect,iii,iii2,found,
     &     MAX_INT,MAXPOINTS_INT)

      implicit none

      integer MAX_INT,MAXPOINTS_INT
      integer iii(MAXPOINTS_INT,MAX_INT),iii2(MAXPOINTS_INT,MAX_INT)
      integer actint,actp0,crossect

      integer found

      integer i, int3, actp3


      if (iii(actp0,actint) .eq. -1) then
         found = 1
c Die naechsten Zeilen sind fuer die Vollstaendigkeit der Einsaetze auskommentiert:
c Auch wenn sich Wellenfront nicht trennt, kan sie sich aber schneiden, und die
c neu interpolierten Strahlen werden die wichtigen Randstrahlen
c      else if (flag_smooth .ne. 0) then
c         found = 0
c Beide Stuecke (links und rechts des Mittelpunkts) sind mit gleich vielen anderen Interfaces coinciding:
      else if (iii2(actp0,actint) .eq. iii2(actp0-1,actint)) then
c fictious interface (beide Stuecke):
         if (iii2(actp0,actint) .lt. 0) then
            found = 0
         else
            actp3 = actp0
            int3 = actint
c Strahl verlaeuft nach oben und coinciding interfaces:
            if (crossect .lt. 0 .and. iii2(actp0,actint) .ne. 0) then
               do i = 1, iii2(actp0,actint)-1
                  actp3 = iii(actp3,int3)
                  int3 = int3 - 1
               end do
            end if
c Eckpunkt dazwischen (fuer 'Strahl laeuft nach unten' oder keine coinciding interfaces doppelte Abfrage):
            if (iii(actp3,int3) .eq. -1) then
               found = 1
c Stuetzpunkt ACTP3 ist glatt und keine Seite ist fictious:
            else if (iii(actp3,int3) .eq. 0 .and. iii(actp3-1,int3) .gt. -2) then
               found = 0
c beide Seiten des Steutzpunktes ACTP3 sind fictious:
            else if (iii(actp3,int3) .eq. iii(actp3-1,int3)) then
               found = 0
            else
               found = 1
            end if
         end if
c unterschiedliche Anzahl coinciding interfaces auf beiden Seiten:
      else
         found = 1
      end if

      end


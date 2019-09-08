c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm REFSTACK / Hauptprogramm WFSTACK
c Wird aufgerufen von WFRONT
c Aufruf des Unterprogramme GET_CRITRAY, RTCALC
c Reflektiert/transmittiert alle eingetroffenen Strahlen des Halts (HALT).
c Ausgabe:
c RAY0, TRAY0, KMAH0, NRAY0: reflektierte, transmittierte Strahlen.
c NEXTSECT : siehe WFRONT.
c
c RTFLAGNEW = 0  : aktueller Strahl wurde reflektiert/transmittiert; != 0 : kein R/T fand statt
c RTFLAGLAST : dasselbe fuer den vorherigen Strahl 
c

      subroutine refstack(halt,thalt,rthalt,amplhalt,phasehalt,anglehalt,inthalt,nrayhalt,kmahhalt,
     &     ficthalt,nexthalt,signrt,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &     ptos,actlayer,oldlayer,veltype,oldveltype,inf_angle,inf_dist,updownflag,flag_surface,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,flag_stack,
     &     ray0,tray0,amplitude0,phase,angle,kmah0,nextsect,nray0,error,f_xrays,
     &     crossect,fictsect,process_2point,ccp,MAX_RECEIVERS,
     &     MAX_INT,MAXPOINTS_INT,N_PARAM,MAX_RAYS,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID,reflpoint)

      implicit none

      integer MAX_INT, MAXPOINTS_INT,N_PARAM,MAX_RAYS,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID
      real halt(N_PARAM,MAX_RAYS), thalt(MAX_RAYS), rthalt(N_RTPARAM,MAX_RAYS), amplhalt(MAX_RAYS)
      real phasehalt(MAX_RAYS), anglehalt(MAX_RAYS)
      integer inthalt, nrayhalt, kmahhalt(MAX_RAYS), nexthalt(0:MAX_RAYS), ficthalt(MAX_RAYS),signrt
      integer npoints_int(MAX_INT), actlayer, oldlayer, veltype, oldveltype, updownflag
      integer flag_surface
      real x_int(MAXPOINTS_INT,MAX_INT), z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT), c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT)
      real ptos(MAX_INT), rho1(MAX_INT), rho2(MAX_INT), inf_angle, inf_dist, dtray
      integer nx_grid(MAX_INT),nz_grid(MAX_INT), flag_stack
      real x_grid(MAXP_XGRID,MAX_INT),z_grid(MAXP_ZGRID,MAX_INT)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID,MAX_INT)
      integer f_xrays

      real ray0(N_PARAM,MAX_RAYS),tray0(MAX_RAYS), amplitude0(MAX_RAYS), phase(MAX_RAYS),angle(MAX_RAYS)
      integer kmah0(MAX_RAYS),nextsect(0:MAX_RAYS),fictsect(MAX_RAYS),crossect(MAX_RAYS),nray0

      integer N_RTPARAM_LOCAL
      parameter (N_RTPARAM_LOCAL = 8)
      real critray1(N_PARAM),traycrit1, amplcrit1, phasecrit1, anglecrit1
      real critray2(N_PARAM),traycrit2, amplcrit2, phasecrit2, anglecrit2
      integer iray,irayold,iray0,iray0old, rtflagnew, rtflaglast, iplus,iplus2,error
      integer iray1,iray2, flag_crit1, flag_crit2, otherlayer, layer2, i,ficttmp, lastfict, flag

      real raynorm(N_PARAM), traynorm, amplnorm, phasenorm, anglenorm
      real rtraynorm(N_RTPARAM_LOCAL)
      real chknormold, chknorm, dzdx
      integer actp
      integer MAX_RECEIVERS
      logical process_2point
      real ccp(2,MAX_RECEIVERS)
      real reflpoint(2,MAX_RAYS)
      real xreflsave, zreflsave

c      rthalt(1) = xint_normal
c      rthalt(2) = zint_normal
c      rthalt(3) = px_unit
c      rthalt(4) = pz_unit
c      rthalt(5) = radcurv
c      rthalt(6) = akapa1
c      rthalt(7) = vs1
c      rthalt(8) = vel(1) bzw 1.0/|p|

      if (flag_stack .ne. 0 .and. signrt .lt. 0) then
         actp = 1
c Checkwert bestimmen, der die Orientierung des Strahls zur Oberflaeche angibt
         call get_dzdx(halt(1,1),x_int(1,inthalt),d_int(1,inthalt),c_int(1,inthalt),
     &        b_int(1,inthalt),npoints_int(inthalt),
     &        dzdx,actp,
     &        MAXPOINTS_INT)
c Skalarprodukt: parallel oder antiparallel?
         chknormold = halt(3,1) + dzdx*halt(4,1)
      end if

      irayold = 0
      iray0 = 1
      iray0old = 0
      rtflaglast = 0
      rtflagnew = 0
      lastfict = 0
c==============================================================
c Schleife ueber alle Strahlen des Halts...
c
      iray = 0
      do while (iray .ne. nrayhalt)
         iray = iray + 1
         rtflaglast = rtflagnew
c--------------------------------
c Neue Ueberpruefung von FICTHALT
c-----------
c Ganz normal
         if (ficthalt(iray) .eq. 0 .or. ficthalt(iray) .eq. 10) then
            fictsect(iray0) = 0
            layer2 = actlayer
            if (signrt .lt. 0) then
               otherlayer = oldlayer - updownflag
            else
               otherlayer = oldlayer
            end if
c-----------
c fictious Interface:
         else if (ficthalt(iray) .lt. 0) then
c-----------
c Transmission
            if (signrt .gt. 0 .and. veltype .eq. oldveltype) then
               fictsect(iray0) = ficthalt(iray) + 1000
               do i = 1, N_PARAM
                  ray0(i,iray0) = halt(i,iray)
               end do
               amplitude0(iray0) = amplhalt(iray)
               phase(iray0) = phasehalt(iray)
               tray0(iray0) = thalt(iray)
               angle(iray0) = anglehalt(iray)
c rtsect eigentlich setzen, wird aber eh' nich gebraucht
               goto 2
c-----------
c Nichts, gehe zum Ende:
            else
               goto 1
            end if
         else
c-----------
c Erste Reflexion/Transmission bei coinciding interfaces:
            if (mod(ficthalt(iray),1000) .ne. 0) then
               ficttmp = (ficthalt(iray) - 1) / 1000
               if (signrt .lt. 0) then
                  fictsect(iray0) = 0
                  otherlayer = oldlayer - updownflag*ficttmp
                  layer2 = actlayer
               else
                  ficttmp = ficttmp - 1
                  fictsect(iray0) = ficttmp*1000
                  layer2 = actlayer + updownflag*ficttmp
                  otherlayer = oldlayer
c-----------
c Strahl erreicht bald die Erdoberflaeche:
c Keine Transmission, sondern direkt rueberschreiben:
                  if (layer2 .eq. 0) then
c Fall: Konversion auf letztem Stueck zur Erdoberflaeche -> nicht moeglich!
                     if (oldveltype .ne. veltype) then
                        goto 1
                     end if
                     do i = 1, N_PARAM
                        ray0(i,iray0) = halt(i,iray)
                     end do
                     amplitude0(iray0) = amplhalt(iray)
                     phase(iray0) = phasehalt(iray)
                     tray0(iray0) = thalt(iray)
                     angle(iray0) = anglehalt(iray)
                     goto 2
                  end if
               end if
c-----------
c Transmission innerhalb coinciding interfaces
            else if (signrt .gt. 0 .and. veltype .eq. oldveltype) then
               fictsect(iray0) = ficthalt(iray) - 1000
               do i = 1, N_PARAM
                  ray0(i,iray0) = halt(i,iray)
               end do
               amplitude0(iray0) = amplhalt(iray)
               phase(iray0) = phasehalt(iray)
               tray0(iray0) = thalt(iray)
               angle(iray0) = anglehalt(iray)
c rtsect eigentlich setzen, wird aber eh' nich gebraucht
               goto 2
c-----------
c Nichts, gehe zum Ende:
            else
               goto 1
            end if
         end if

c Bei Reflexion und stack: rote Lampen leuchten auf!
         if (flag_stack .ne. 0 .and. signrt .lt. 0) then
c Checkwert bestimmen, der die Orientierung des Strahls zur Oberflaeche angibt
            call get_dzdx(halt(1,iray),x_int(1,inthalt),d_int(1,inthalt),c_int(1,inthalt),
     &           b_int(1,inthalt),npoints_int(inthalt),
     &           dzdx,actp,
     &           MAXPOINTS_INT)
            chknorm = halt(3,iray) + dzdx*halt(4,iray)
c hier muss geprueft werden, dass sich der aktuelle Schnittpunkt INMITTEN eines Grenzflaechenstueckes
c befindet, an dem reflektiert wird. Liegt er am RANDE, muss der Wert CHKNORMOLD neu gesetzt werden
            if (iray .gt. 1 .and. ficthalt(iray-1) .ne. 0) then
               if (mod(ficthalt(iray-1),1000) .eq. 0) then
                  chknormold = chknorm
               end if
            end if
c 90 Grad ('normaler') incident ray ist vorhanden, er liegt zwischen diesem und dem letzten Strahl 
            if (chknorm*chknormold .lt. 0) then
               iray1 = iray - 1
               iray2 = iray
c Korrigiert ptos --> ptos(oldlayer)...hoffe das ist korrekt
               call get_normray(halt(1,iray1),thalt(iray1),rthalt(1,iray1),halt(1,iray2),thalt(iray2),
     &              rthalt(1,iray2),amplhalt(iray1),amplhalt(iray2),phasehalt(iray1),phasehalt(iray2),
     &              anglehalt(iray1),anglehalt(iray2),
     &              oldlayer,oldveltype,inthalt,
     &              npoints_int,x_int,z_int,b_int,c_int,d_int,ptos(oldlayer),
     &              nx_grid(oldlayer),nz_grid(oldlayer),x_grid(1,oldlayer),z_grid(1,oldlayer),
     &              v_grid(1,1,1,oldlayer),
     &              chknormold,chknorm,
     &              raynorm,amplnorm,phasenorm,anglenorm,traynorm,rtraynorm,flag,error,
     &              MAX_INT,MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
c               flag = 1
               if (error .ne. 0) goto 999
c caustic: Setze einfach einen der beiden Strahlen
               if (flag .ne. 0) then
                  do i = 1, N_PARAM
                     raynorm(i) = halt(i,iray)
                  end do
                  do i = 1, N_RTPARAM
                     rtraynorm(i) = rthalt(i,iray)
                  end do
                  phasenorm = phasehalt(iray)
                  amplnorm = amplhalt(iray)
                  traynorm = thalt(iray)
                  anglenorm = anglehalt(iray)
               end if
               call rtcalc(raynorm,rtraynorm,amplnorm,phasenorm,oldlayer,layer2,
     &              otherlayer,oldveltype,veltype,inthalt,signrt,ptos,flag_surface,
     &              nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,
     &              ray0(1,iray0),amplitude0(iray0),phase(iray0),rtflagnew,
     &              MAX_INT,N_RTPARAM,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
               tray0(iray0) = traynorm
c               angle(iray0) = anglenorm
c OBS!! Nur zum Test gedacht, sollte spaeter geaendert werden!
               angle(iray0) = 2. * traynorm
               chknormold = chknorm
               goto 2
            end if
            chknormold = chknorm
            goto 1
c normal, kein stack oder stack ohne Reflexion
         else
            call rtcalc(halt(1,iray),rthalt(1,iray),amplhalt(iray),phasehalt(iray),oldlayer,layer2,
     &           otherlayer,oldveltype,veltype,inthalt,signrt,ptos,flag_surface,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,
     &           ray0(1,iray0),amplitude0(iray0),phase(iray0),rtflagnew,
     &           MAX_INT,N_RTPARAM,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            if( process_2point .and. signrt .lt. 0 .and. rtflagnew .eq. 0 ) then
c               write(*,*) "Reflection ",ray0(1,iray0),ray0(2,iray0)

c               if( iray0 .le. MAX_RECEIVERS ) then
c                  ccp(1,iray0) = ray0(1,iray0)
c                  ccp(2,iray0) = ray0(2,iray0)
c               else
c                  write(*,*) "Error in REFLSTACK: More rays than receivers..."
c               end if
            else if( signrt .lt. 0 ) then
c Speicherung des Reflexionspunktes
               ray0(7,iray0) = halt(1,iray)
               ray0(8,iray0) = halt(2,iray)
c               write(*,*) "RAY ", ray0(7,iray0), ray0(8,iray0), rtflagnew, process_2point
            end if
         end if

c========================================================================================================
c Kritischen Winkel suchen
c
c Aber: falls vom vorherigen Strahl FICTHALT = 10, 1011, 2011... -> keine Suche!
c Fall 'FICTHALT(iray) .ne. FICTHALT(iray-1)' kann nur auftreten, wenn Modell geglaettet wurde (FLAG_SMOOTH .ne. 0)
c
         if (rtflagnew .ne. rtflaglast .and. (lastfict .eq. 0 .and. nexthalt(iray-1) .eq. iray .and.
     &        kmahhalt(iray) .eq. kmahhalt(iray-1) .and. ficthalt(iray) .eq. ficthalt(iray-1))) then
            xreflsave = ray0(7,iray0)
            zreflsave = ray0(8,iray0)
c IRAY1 ist der Strahl, der noch transmittiert wurde
c IRAY2 ist der ueberkritische Strahl
c--------------------------------------
c Fall: Strahl IRAY ist ueberkritisch
c
            if (rtflagnew .ne. 0) then
               iray1 = iray-1
               iray2 = iray
               iplus = iray0 + 1
               if (iplus .gt. MAX_RAYS) goto 990
               call get_critray( halt(1,iray1), thalt(iray1), rthalt(1,iray1), halt(1,iray2), thalt(iray2),
     &              amplhalt(iray1),amplhalt(iray2),phasehalt(iray1),phasehalt(iray2),
     &              anglehalt(iray1),anglehalt(iray2),
     &              inf_angle,inf_dist,oldlayer,layer2,otherlayer,oldveltype,veltype,inthalt,signrt,
     &              npoints_int,x_int,z_int,b_int,c_int,d_int,ptos,
     &              nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,
     &              ray0(1,iplus),amplitude0(iplus),phase(iplus),angle(iplus),tray0(iplus),flag_crit1,
     &              ray0(1,iray0),amplitude0(iray0),phase(iray0),angle(iray0),tray0(iray0),flag_crit2,error,
     &              MAX_INT,MAXPOINTS_INT,N_RTPARAM,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
               if (error .ne. 0) goto 999
               if (flag_crit1 .ne. 0) then
                  goto 1
               else if (flag_crit2 .ne. 0) then
                  do i = 1, N_PARAM
                     ray0(i,iray0) = ray0(i,iplus)
                  end do
                  amplitude0(iray0) = amplitude0(iplus)
                  phase(iray0) = phase(iplus)
                  angle(iray0) = angle(iplus)
                  tray0(iray0) = tray0(iplus)
               else
                  if (f_xrays .ne. 0) write(f_xrays,*) ray0(1,iray0),ray0(2,iray0)," critical ray (2)"
                  kmah0(iray0) = kmahhalt(iray)
                  ray0(7,iray0) = xreflsave
                  ray0(8,iray0) = zreflsave
                  crossect(iray0) = 0
                  nextsect(iray0-1) = iray0
                  iray0 = iplus
                  fictsect(iray0) = fictsect(iray0-1)
               end if
               if (f_xrays .ne. 0) write(f_xrays,*) ray0(1,iray0),ray0(2,iray0)," critical ray"
               kmah0(iray0) = kmahhalt(iray)
               crossect(iray0) = 0
               nextsect(iray0-1) = iray0
               ray0(7,iray0) = xreflsave
               ray0(8,iray0) = zreflsave
               iray0 = iray0 + 1
               goto 1
c--------------------------------------
c Fall: Strahl IRAY-1 ist ueberkritisch
c
            else
               iray1 = iray
               iray2 = iray-1
               halt(7,iray2) = halt(7,iray1)
               halt(8,iray2) = halt(8,iray1)
               call get_critray(halt(1,iray1),thalt(iray1),rthalt(1,iray1),halt(1,iray2),thalt(iray2),
     &              amplhalt(iray1),amplhalt(iray2),phasehalt(iray1),phasehalt(iray2),
     &              anglehalt(iray1),anglehalt(iray2),
     &              inf_angle,inf_dist,oldlayer,layer2,otherlayer,oldveltype,veltype,inthalt,signrt,
     &              npoints_int,x_int,z_int,b_int,c_int,d_int,ptos,
     &              nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,
     &              critray1,amplcrit1,phasecrit1,anglecrit1,traycrit1,flag_crit1,
     &              critray2,amplcrit2,phasecrit2,anglecrit2,traycrit2,flag_crit2,error,
     &              MAX_INT,MAXPOINTS_INT,N_RTPARAM,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
               if (error .ne. 0) goto 999
               
               if (flag_crit1 .ne. 0) then
                  tray0(iray0) = thalt(iray)
                  angle(iray0) = anglehalt(iray)
                  kmah0(iray0) = kmahhalt(iray)
                  crossect(iray0) = 0
                  nextsect(iray0-1) = -iray0
                  ray0(7,iray0) = xreflsave
                  ray0(8,iray0) = zreflsave
                  iray0 = iray0 + 1
                  irayold = iray
                  goto 1
               else if (flag_crit2 .ne. 0) then
                  iplus = iray0 + 1
                  if (iplus .gt. MAX_RAYS) goto 990
                  
c Schreibe transmittierten Originalstrahl in IPLUS
                  do i = 1, N_PARAM
                     ray0(i,iplus) = ray0(i,iray0)
                  end do
                  amplitude0(iplus) = amplitude0(iray0)
                  phase(iplus) = phase(iray0)
                  fictsect(iplus) = fictsect(iray0)
                  crossect(iplus) = 0
                  angle(iplus) = anglehalt(iray)
                  tray0(iplus) = thalt(iray)
                  kmah0(iplus) = kmahhalt(iray)

c Schreibe kritischen Strahl in IRAY0
                  do i = 1, N_PARAM
                     ray0(i,iray0) = critray1(i)
                  end do
                  amplitude0(iray0) = amplcrit1
                  phase(iray0) = phasecrit1
                  angle(iray0) = anglecrit1
                  tray0(iray0) = traycrit1
                  crossect(iray0) = 0
                  kmah0(iray0) = kmahhalt(iray)
                  
                  if (f_xrays .ne. 0) write(f_xrays,*) ray0(1,iray0),ray0(2,iray0)," critical ray"

                  ray0(7,iplus) = xreflsave
                  ray0(8,iplus) = zreflsave
                  ray0(7,iray0) = xreflsave
                  ray0(8,iray0) = zreflsave
                  nextsect(iray0-1) = -iray0
                  iray0 = iplus
                  nextsect(iray0-1) = iray0
                  iray0 = iray0 + 1
                  irayold = iray
                  goto 1
               else
                  iplus = iray0 + 1
                  iplus2 = iray0 + 2
                  if (iplus2 .gt. MAX_RAYS) goto 990

c Schreibe transmittierten Originalstrahl in IPLUS2
                  do i = 1, N_PARAM
                     ray0(i,iplus2) = ray0(i,iray0)
                  end do
                  amplitude0(iplus2) = amplitude0(iray0)
                  phase(iplus2) = phase(iray0)
                  fictsect(iplus2) = fictsect(iray0)
                  crossect(iplus2) = 0
                  angle(iplus2) = anglehalt(iray)
                  tray0(iplus2) = thalt(iray)
                  kmah0(iplus2) = kmahhalt(iray)

c Schreibe kritischen Strahl in IRAY0
                  do i = 1, N_PARAM
                     ray0(i,iray0) = critray1(i)
                  end do
                  amplitude0(iray0) = amplcrit1
                  phase(iray0) = phasecrit1
                  angle(iray0) = anglecrit1
                  tray0(iray0) = traycrit1
                  crossect(iray0) = 0
                  kmah0(iray0) = kmahhalt(iray)
                  if (f_xrays .ne. 0) write(f_xrays,*) ray0(1,iray0),ray0(2,iray0)," critical ray"

c Schreibe Extra-Strahl in IPLUS
                  do i = 1, N_PARAM
                     ray0(i,iplus) = critray2(i)
                  end do
                  amplitude0(iplus) = amplcrit2
                  phase(iplus) = phasecrit2
                  angle(iplus) = anglecrit2
                  tray0(iplus) = traycrit2
                  fictsect(iplus) = fictsect(iray0)
                  crossect(iplus) = 0
                  kmah0(iplus) = kmahhalt(iray)
                  if (f_xrays .ne. 0) write(f_xrays,*) ray0(1,iplus),ray0(2,iplus)," critical ray (2)"

                  nextsect(iray0-1) = -iray0
                  nextsect(iray0) = iplus
                  nextsect(iplus) = iplus2
                  
                  ray0(7,iray0) = xreflsave
                  ray0(8,iray0) = zreflsave
                  ray0(7,iplus) = xreflsave
                  ray0(8,iplus) = zreflsave
                  ray0(7,iplus2) = xreflsave
                  ray0(8,iplus2) = zreflsave
                  iray0 = iplus2 + 1
                  irayold = iray
                  goto 1
               end if
            end if
c
c Ende kritischen Winkel suchen
c========================================================================================================
         else if (rtflagnew .ne. 0) then
            goto 1
         else
            tray0(iray0) = thalt(iray)
            angle(iray0) = anglehalt(iray)
         end if
            
 2       kmah0(iray0) = kmahhalt(iray)
         crossect(iray0) = 0

         if ((iray-irayold) .eq. 1) then
            nextsect(iray0-1) = isign(1,nexthalt(irayold))*iray0
         else
            nextsect(iray0-1) = -iray0
         end if
 
         if (lastfict .ne. 0) then
            if (signrt .gt. 0) then
               nextsect(iray0-1) = -iray0
            end if
         end if

         irayold = iray
         iray0 = iray0 + 1

 1       if (iray0 .gt. MAX_RAYS) goto 990
         lastfict = (ficthalt(iray) - (ficthalt(iray)/1000)*1000) / 10
      end do

      nray0 = iray0-1
      nextsect(0) = -1
      nextsect(nray0) = 0

c setze alle next's auf -1, das heisst, keine Verwandtschaft zwischen Strahlen
      if (flag_stack .ne. 0 .and. signrt .lt. 0) then
         do iray = 1, nray0 -1
            nextsect(iray) = -abs(nextsect(iray))
         end do
      end if

      goto 999
 990  write(*,*) "REFLTRANS: Too many rays!"
      write(*,*) "Enlarge PARAMETER 'MAX_RAYS' in WFRONT"
      error = 1

 999  end

c***************************************************************************************************
c Unterprogramm RTCALC / Hauptprogramm WFRONT
c Wird aufgerufen von REFLTRANS, GET_RECORD, GET_CRITRAY
c Aufruf der Unterprogramme COEF8, VELOCITY
c Berechnet Reflexion/Transmission eines Strahls
c FLAG = -1  : Keine Reflexion/Transmission erfolgt
c FLAG = 0   : R/T erfolgt
c SIGNRT = -1: Reflexion
c              OLDLAYER  : Schicht, aus der Strahl kommt
c              ACTLAYER  : Schicht, in die Strahl reflektiert wird (= OLDLAYER)
c              OTHERLAYER: Schicht auf der anderen Seite der Grenzflaeche
c SIGNRT = +1: Transmission
c              OLDLAYER  : Schicht, aus der Strahl kommt
c              ACTLAYER  : Schicht, in die Strahl transmittiert wird
c              OTHERLAYER: Schicht auf der anderen (alten) Seite der Grenzflaeche (=OLDLAYER)
c
      subroutine rtcalc(rayold,rtray,amplold,phaseold,oldlayer,actlayer,
     &     otherlayer,oldveltype,veltype,actint,signrt,ptos,flag_surface,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,
     &     raynew,amplnew,phasenew,flag,
     &     MAX_INT,N_RTPARAM,N_PARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer MAX_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID
      integer actlayer, otherlayer, veltype, oldveltype,signrt, oldlayer, actint, flag_surface
      real rayold(N_PARAM), rtray(N_RTPARAM), amplold, phaseold
      real ptos(MAX_INT),rho1(MAX_INT),rho2(MAX_INT)
      integer nx_grid(MAX_INT),nz_grid(MAX_INT)
      real x_grid(MAXP_XGRID,MAX_INT),z_grid(MAXP_ZGRID,MAX_INT)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID,MAX_INT)

      real raynew(N_PARAM), amplnew, phasenew
      integer flag

      real vel(6), vel1, root, extra, akapa2, pz2_unit, vs2, s1
      integer i
      real p,vp1,vshear1,vp2,vshear2,ro1,ro2,r,phs
      integer nc,nd


      flag = 0

      call velocity(rayold(1), rayold(2), nx_grid(actlayer), nz_grid(actlayer),
     &     x_grid(1,actlayer), z_grid(1,actlayer), v_grid(1,1,1,actlayer),
     &     veltype, ptos(actlayer),actlayer,
     &     vel,
     &     MAXP_XGRID, MAXP_ZGRID)

      root = 1.0/(vel(1)*vel(1)) +
     &     (rtray(4)*rtray(4) - 1.0) / (rtray(8)*rtray(8))

      if (root .le. 0.0) then
         flag = -1
         goto 999
      end if

c=====================
c neue slowness px,pz:
c
      root = signrt*sqrt(root)
      extra = rtray(4)/rtray(8) + root
      
      raynew(3) = rayold(3) - extra*rtray(1)
      raynew(4) = rayold(4) - extra*rtray(2)
      
c================
c neue Q, P:
c
      pz2_unit = -root*vel(1)
      akapa2 = vel(2)*raynew(4) - vel(3)*raynew(3)
      vs2 = (vel(2)*raynew(3) + vel(3)*raynew(4))*vel(1)

      s1 = (rtray(7)-vs2)*(rtray(3)*rtray(3)/(rtray(8)*rtray(8)))
      s1 = s1 + rtray(5)*(pz2_unit/vel(1) - rtray(4)/rtray(8))
      s1 = s1 + 2.*(akapa2*pz2_unit - rtray(4)*rtray(6))*rtray(3)/rtray(8)

      raynew(5) = rayold(5)*pz2_unit/rtray(4)
      raynew(6) = (rayold(6)*rtray(4) - rayold(5)*s1/rtray(4))/pz2_unit

c=======================================
c x-,z-Koordinaten, Zeit bleiben gleich:
      do i = 1, 2
         raynew(i) = rayold(i)
      end do
c Reflexionspukt bleibt gleich
      do i = 7, 8
         raynew(i) = rayold(i)
      end do

c==========================================================
c Amplitude : Anew = Aold * SQRT( V2*|sina1| / V1*|sina2| )

c-----------------------------------------
c Setzen der Werte der ersten Schicht (oldlayer):
c
      if (oldveltype .gt. 0) then
         vp1 = rtray(8)
         vshear1 = vp1 / ptos(oldlayer)
         if (signrt .lt. 0) then
            if (veltype .gt. 0) then
               nc = 1
            else
               nc = 2
            end if
         else
            if (veltype .gt. 0) then
               nc = 3
            else
               nc = 4
            end if
         end if            
      else
         vshear1 = rtray(8)
         vp1 = vshear1 * ptos(oldlayer)
         if (signrt .lt. 0) then
            if (veltype .gt. 0) then
               nc = 5
            else
               nc = 6
            end if
         else
            if (veltype .gt. 0) then
               nc = 7
            else
               nc = 8
            end if
         end if         
      end if
      ro1 = rho1(oldlayer) + rho2(oldlayer)*vp1

c--------------------------------------------
c Setzen der Werte der zweiten Schicht: (R: otherlayer, T: actlayer)
c
c Reflexion:
      if (signrt .lt. 0) then
         if (actint .eq. 1 .or. otherlayer .eq. 0) then
c            write(*,*) "Reflection at earth's surface!!"
c FLAG_SURFACE != 0: Keine Reflexion an Erdoberflaeche:
            if (flag_surface .ne. 0) then
               flag = -1
               goto 999
            end if
            ro2 = 0.0
         else
            call velonly(rayold(1), rayold(2), nx_grid(otherlayer), nz_grid(otherlayer),
     &           x_grid(1,otherlayer), z_grid(1,otherlayer), v_grid(1,1,1,otherlayer),
     &           1,ptos(otherlayer),otherlayer,
     &           vel1,
     &           MAXP_XGRID, MAXP_ZGRID)
            vp2 = vel1
            vshear2 = vp2 / ptos(otherlayer)
            ro2 = rho1(otherlayer) + rho2(otherlayer)*vp2
         end if
c Transmission:
      else
         if (veltype .gt. 0) then
            vp2 = vel(1)
            vshear2 = vp2 / ptos(actlayer)
         else
            vshear2 = vel(1)
            vp2 = vshear2 * ptos(actlayer)
         end if
         ro2 = rho1(actlayer) + rho2(actlayer)*vp2
      end if

      p = rtray(3) / rtray(8)
      if (p .lt. 0.0) then
         nd = 1
         p = abs(p)
      else
         nd = 0
      end if

c==============================================================
c Berechnung des Reflexions/Transmissionskoeffizienten 'R' (sowie Phasenverschiebung 'PHS'):
c
      call coef8(p,vp1,vshear1,ro1,vp2,vshear2,ro2,nc,nd,r,phs)
c
c==============================================================

c Bei Reflexion bleibt Dichte gleich
      if (signrt .lt. 0) then
         amplnew = amplold * sqrt(vel(1) * abs(pz2_unit/(rtray(8)*rtray(4)))) * r
      else
         if (oldveltype .lt. 0) then
            ro1 = rho1(oldlayer) + rho2(oldlayer)*rtray(8)*ptos(oldlayer)
         else
            ro1 = rho1(oldlayer) + rho2(oldlayer)*rtray(8)
         end if
         if (veltype .lt. 0) then
            ro2 = rho1(actlayer) + rho2(actlayer)*vel(1)*ptos(actlayer)
         else
            ro2 = rho1(actlayer) + rho2(actlayer)*vel(1)
         end if
         amplnew = amplold * sqrt(abs( (vel(1)*ro2*pz2_unit) / (rtray(8)*rtray(4)*ro1) )) * r
      end if

      phasenew = phaseold + phs

 999  end

c**********************************************************************************
c Unterprogramm GET_CRITRAY / Hauptprogramm WFRONT
c Wird aufgerufen von REFLTRANS
c Aufruf der Unterprogramme GET_RAYONINT, RTCALC, GET_RTSECT, PARAXIAL_P, GET_Z1, GET_Z, VELOCITY, VELONLY
c Bestimmung eines Strahles RAYNEW0, der gerade noch transmittiert wird (critical ray)
c Ausserdem: Bestimmung eines zweiten Strahles RAYNEW2, der zwischen dem zuletzt transmittierten
c Strahl und dem critical ray verlaeuft (dies zur besseren Aufloesung des Transmissionskoeffizienten)
c RAY1 ist der Strahl, der noch transmittiert wurde
c RAY2 ist bereits ueberkritisch
c FLAG = 0: critical ray gefunden
c FLAG !=0: critical ray nicht gefunden
c
      subroutine get_critray(ray1,tray1,rtray1,ray2,tray2,
     &     ampl1,ampl2,phase1,phase2,angle1,angle2,
     &     inf_angle,INF_DIST,oldlayer,actlayer,otherlayer,oldveltype,veltype,actint,signrt,
     &     npoints_int,x_int,z_int,b_int,c_int,d_int,ptos,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,
     &     raynew0,amplnew0,phasenew0,anglenew0,tnew0,flag1,
     &     raynew2,amplnew2,phasenew2,anglenew2,tnew2,flag2,error,
     &     MAX_INT,MAXPOINTS_INT,N_RTPARAM,N_PARAM,MAXP_XGRID,MAXP_ZGRID)

c Pruefen: ASIN und ACOS pruefen!

      implicit none

      integer MAX_INT, MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID
      real ray1(N_PARAM),ray2(N_PARAM),rtray1(N_RTPARAM),tray1,tray2
      real ampl1,phase1,angle1,ampl2,phase2,angle2,dtray
      integer oldlayer,actlayer,otherlayer,oldveltype,veltype,actint,signrt
      integer npoints_int(MAX_INT)
      real x_int(MAXPOINTS_INT,MAX_INT), z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT), c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT)
      real ptos(MAX_INT), rho1(MAX_INT), rho2(MAX_INT), inf_angle, INF_DIST
      integer nx_grid(MAX_INT),nz_grid(MAX_INT)
      real x_grid(MAXP_XGRID,MAX_INT),z_grid(MAXP_ZGRID,MAX_INT)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID,MAX_INT)

      real raynew0(N_PARAM),amplnew0,phasenew0,tnew0,anglenew0
      real raynew2(N_PARAM),amplnew2,phasenew2,tnew2,anglenew2
      integer flag1,flag2,error

      integer N_RTPARAM_LOCAL, MAX_COUNTS
      parameter (N_RTPARAM_LOCAL = 8, MAX_COUNTS = 30)
      real raynew(N_PARAM),amplnew,phasenew,rtnew(N_RTPARAM_LOCAL)
      real velact1,velactnew,vold1(6),vold2(6),ptmp(2),v1_v2,dzdx,dx,d2zdx2
      real totaldist,part,x1,x2,dummy,inf_angle2
      integer actp,found,count, flag_surface

      double precision pi,sin_inf,phi,phi_crit,phi_suff,sin_crit

      if (N_RTPARAM_LOCAL .ne. N_RTPARAM) goto 980
c      if (N_PARAM_LOCAL .ne. N_PARAM) goto 990

      if (abs(ray1(5)) .lt. 0.00001 .or. abs(ray2(5)) .lt. 0.00001) then
         write(*,*) 'GET_CRITRAY: CAUSTIC!'
         flag1 = -1
         goto 999
      end if

      pi = 3.141592653589
      flag_surface = 0
c---------------
c Aendern: INF_ANGLE abhaengig von d2zdx2 (Kruemmung der Oberflaeche) machen! Und vom Geschw.-Gradienten!

      inf_angle2 = inf_angle
      sin_inf = sin(pi*0.5 - inf_angle2)
      sin_crit = sin(0.5*(pi - inf_angle2))
      actp = 1

      call velonly(ray1(1), ray1(2), nx_grid(actlayer), nz_grid(actlayer),
     &     x_grid(1,actlayer), z_grid(1,actlayer), v_grid(1,1,1,actlayer),
     &     veltype, ptos(actlayer),actlayer,
     &     velact1,
     &     MAXP_XGRID, MAXP_ZGRID)

      v1_v2 = rtray1(8) / velact1
c Es kann keinen kritischen Strahl geben, weil Geschwindigkeit abnimmt!
      if (abs(v1_v2) .gt. 1.0) then
         flag1 = -1
         goto 999
      end if

      phi = pi - acos(rtray1(4))
      phi_crit = asin(v1_v2 * sin_crit)
      phi_suff = asin(v1_v2 * sin_inf)

c-----------------------
c Einfallswinkel von Strahl 1 liegt nicht im gesuchten Genauigkeitsbereich
c
 100  continue
      if (phi .lt. phi_suff) then
c Geschwindigkeiten in der alten Schicht:
         call velocity(ray1(1), ray1(2), nx_grid(oldlayer), nz_grid(oldlayer),
     &        x_grid(1,oldlayer), z_grid(1,oldlayer), v_grid(1,1,1,oldlayer),
     &        oldveltype, ptos(oldlayer),oldlayer,
     &        vold1,
     &        MAXP_XGRID, MAXP_ZGRID)
         call velocity(ray2(1), ray2(2), nx_grid(oldlayer), nz_grid(oldlayer),
     &        x_grid(1,oldlayer), z_grid(1,oldlayer), v_grid(1,1,1,oldlayer),
     &        oldveltype,ptos(oldlayer),oldlayer,
     &        vold2,
     &        MAXP_XGRID, MAXP_ZGRID)
c Neuer Punkt mitten zwischen RAY1 und RAY2
         raynew(1) = 0.5*(ray1(1) + ray2(1))
         call get_z1(raynew(1),x_int(1,actint),z_int(1,actint),d_int(1,actint),c_int(1,actint),b_int(1,actint),
     &        npoints_int(actint),
     &        raynew(2),dzdx,actp,
     &        MAXPOINTS_INT)
c---------
c par. Approximation der slowness sowie Mittelung
         call paraxial_p(ray1,vold1,
     &        raynew,N_PARAM)
         ptmp(1) = raynew(3)
         ptmp(2) = raynew(4)
         call paraxial_p(ray2,vold2,
     &        raynew,N_PARAM)
         raynew(3) = 0.5 * (ptmp(1) + raynew(3))
         raynew(4) = 0.5 * (ptmp(2) + raynew(4))
c---------
         call get_rtsect(raynew,dzdx,dummy,
     &        nx_grid(oldlayer),nz_grid(oldlayer),x_grid(1,oldlayer),z_grid(1,oldlayer),
     &        v_grid(1,1,1,oldlayer),ptos(oldlayer),oldveltype,oldlayer,
     &        rtnew,
     &        N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
         call velonly(raynew(1), raynew(2), nx_grid(actlayer), nz_grid(actlayer),
     &        x_grid(1,actlayer), z_grid(1,actlayer), v_grid(1,1,1,actlayer),
     &        veltype, ptos(actlayer),actlayer,
     &        velactnew,
     &        MAXP_XGRID, MAXP_ZGRID)
         v1_v2 = rtnew(8) / velactnew
c Es kann keinen kritischen Strahl geben, weil Geschwindigkeit abnimmt!
         if (abs(v1_v2) .gt. 1.0) then
            flag1 = -1
            goto 999
         end if

         phi = pi - acos(rtnew(4))
         phi_crit = asin(v1_v2 * sin_crit)
         phi_suff = asin(v1_v2 * sin_inf)

         x1 = ray1(1)
         x2 = ray2(1)
         totaldist = abs(x1-x2)
         if (totaldist .lt. 0.000001) totaldist = 0.000001
c-----------------------------------------
c Iteration:
c
         found = 0
         count = 0
         do while (found .eq. 0 .and. count .lt. MAX_COUNTS)
            count = count + 1
            if (phi .lt. phi_suff) then
               x1 = raynew(1)
            else if (phi .gt. phi_crit) then
               x2 = raynew(1)
            else
               found = 1
            end if
            if (found .eq. 0) then
               raynew(1) = 0.5 * (x1 + x2)
               call get_z1(raynew(1),x_int(1,actint),z_int(1,actint),d_int(1,actint),c_int(1,actint),b_int(1,actint),
     &              npoints_int(actint),
     &              raynew(2),dzdx,actp,
     &              MAXPOINTS_INT)
               call paraxial_p(ray1,vold1,
     &              raynew,N_PARAM)
               ptmp(1) = raynew(3)
               ptmp(2) = raynew(4)
               call paraxial_p(ray2,vold2,
     &              raynew,N_PARAM)
               part = abs(ray1(1)-raynew(1)) / totaldist
               raynew(3) = ptmp(1) + part * (raynew(3) - ptmp(1))
               raynew(4) = ptmp(2) + part * (raynew(4) - ptmp(2))
               call get_rtsect(raynew,dzdx,dummy,
     &              nx_grid(oldlayer),nz_grid(oldlayer),x_grid(1,oldlayer),z_grid(1,oldlayer),
     &              v_grid(1,1,1,oldlayer),ptos(oldlayer),oldveltype,oldlayer,
     &              rtnew,
     &              N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
               call velonly(raynew(1), raynew(2), nx_grid(actlayer), nz_grid(actlayer),
     &              x_grid(1,actlayer), z_grid(1,actlayer), v_grid(1,1,1,actlayer),
     &              veltype, ptos(actlayer),actlayer,
     &              velactnew,
     &              MAXP_XGRID, MAXP_ZGRID)
               v1_v2 = rtnew(8) / velactnew
c Es kann keinen kritischen Strahl geben, weil Geschwindigkeit abnimmt!
               if (abs(v1_v2) .gt. 1.0) then
                  flag1 = -1
                  goto 999
               end if
c               phi_suff = asin(v1_v2 * sin_inf)
c               phi_crit = 0.5*(asin(v1_v2) + phi_suff)

               phi = pi - acos(rtnew(4))
               phi_crit = asin(v1_v2 * sin_crit)
               phi_suff = asin(v1_v2 * sin_inf)
            end if
         end do

c Berechne ALLE Parameter des neuen Punktes
         call get_rayonint(ray1,ampl1,phase1,tray1,angle1,ray2,ampl2,phase2,tray2,angle2,
     &        tnew0,oldveltype,ptos(oldlayer),oldlayer,dtray,
     &        nx_grid(oldlayer),nz_grid(oldlayer),x_grid(1,oldlayer),z_grid(1,oldlayer),
     &        v_grid(1,1,1,oldlayer),
     &        raynew,amplnew,phasenew,anglenew0,flag1,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         if (flag1 .ne. 0) goto 999

c Berechne ALLE R/T-Paramter
         dx = raynew(1) - x_int(actp,actint)
         d2zdx2 = 6.*d_int(actp,actint)*dx + 2.*c_int(actp,actint)
         call get_rtsect(raynew,dzdx,d2zdx2,
     &        nx_grid(oldlayer),nz_grid(oldlayer),x_grid(1,oldlayer),z_grid(1,oldlayer),
     &        v_grid(1,1,1,oldlayer),ptos(oldlayer),oldveltype,oldlayer,
     &        rtnew,
     &        N_RTPARAM_LOCAL,MAXP_XGRID,MAXP_ZGRID)

         call rtcalc(raynew,rtnew,amplnew,phasenew,oldlayer,actlayer,
     &        otherlayer,oldveltype,veltype,actint,signrt,ptos,flag_surface,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,
     &        raynew0,amplnew0,phasenew0,flag1,
     &        MAX_INT,N_RTPARAM_LOCAL,N_PARAM,MAXP_XGRID,MAXP_ZGRID)

         if (flag1 .ne. 0) then
c Im Falle dass: Diese Zeilen muessten die Loesung sein:
            inf_angle2 = 2.*inf_angle2
            v1_v2 = rtray1(8) / velact1
            sin_inf = sin(pi*0.5 - inf_angle2)
            phi = pi - acos(rtray1(4))
            phi_crit = asin(v1_v2 * sin_crit)
            phi_suff = asin(v1_v2 * sin_inf)
            goto 100
         end if

c==============================================================================================
c Bestimmung eines weiteren Strahls nahe dem kritischen Strahl (Genauigkeit der dyn. Parameter!)
c
         if (ray1(1) .gt. ray2(1)) then
            raynew(1) = raynew0(1) + INF_DIST
            if (raynew(1) .ge. ray1(1)) then
               flag2 = -1
               goto 999
            end if
         else
            raynew(1) = raynew0(1) - INF_DIST
            if (raynew(1) .le. ray1(1)) then
               flag2 = -1
               goto 999
            end if
         end if
         call get_z(raynew(1),x_int(1,actint),z_int(1,actint),d_int(1,actint),c_int(1,actint),b_int(1,actint),
     &        npoints_int(actint),
     &        raynew(2),actp,
     &        MAXPOINTS_INT)

         call get_rayonint(ray1,ampl1,phase1,tray1,angle1,ray2,ampl2,phase2,tray2,angle2,
     &        tnew2,oldveltype,ptos(oldlayer),oldlayer,dtray,
     &        nx_grid(oldlayer),nz_grid(oldlayer),x_grid(1,oldlayer),z_grid(1,oldlayer),
     &        v_grid(1,1,1,oldlayer),
     &        raynew,amplnew,phasenew,anglenew2,flag2,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
         if (flag2 .ne. 0) goto 999

c Berechne ALLE R/T-Paramter
         dx = raynew(1) - x_int(actp,actint)
         d2zdx2 = 6.*d_int(actp,actint)*dx + 2.*c_int(actp,actint)
         call get_rtsect(raynew,dzdx,d2zdx2,
     &        nx_grid(oldlayer),nz_grid(oldlayer),x_grid(1,oldlayer),z_grid(1,oldlayer),
     &        v_grid(1,1,1,oldlayer),ptos(oldlayer),oldveltype,oldlayer,
     &        rtnew,
     &        N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)

         call rtcalc(raynew,rtnew,amplnew,phasenew,oldlayer,actlayer,
     &        otherlayer,oldveltype,veltype,actint,signrt,ptos,flag_surface,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,
     &        raynew2,amplnew2,phasenew2,flag2,
     &        MAX_INT,N_RTPARAM,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
c================================================================================================
      else
c         write(*,*) 'ray already critical'
         flag1 = -1
      end if

      goto 999
 980  write(*,*) "GET_CRITRAY : Parameter N_RTPARAM_LOCAL not equal N_RTPARAM from WFRONT"
      write(*,*) "Change N_RTPARAM_LOCAL to the amount of N_RTPARAM"
      error = 22
      goto 999
 990  write(*,*) "GET_CRITRAY : Parameter N_PARAM_LOCAL not equal N_PARAM from WFRONT"
      write(*,*) "Change N_PARAM_LOCAL to the amount of N_PARAM"
      error = 22

 999  end


c**********************************************************************************
c Unterprogramm GET_NORMRAY / Hauptprogramm WFSTACK
c Wird aufgerufen von REFSTACK
c Aufruf der Unterprogramme PARAXIAL, GET_Z1, GET_DZDX, VELOCITY, VELONLY
c Bestimmung eines Strahles RAYNEW, der senkrecht auf die Grenzflaeche einfaellt
c (NORMaler Strahl).
c Indikator fuer den senkrechten Einfall ist ein Kreuzprodukt aus Slowness
c und einem tangentialen Vektor auf der Grenzflaeche
c
      subroutine get_normray(ray1,tray1,rtray1,ray2,tray2,rtray2,
     &     ampl1,ampl2,phase1,phase2,angle1,angle2,
     &     actlayer,veltype,actint,
     &     npoints_int,x_int,z_int,b_int,c_int,d_int,ptos,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     chknorm1, chknorm2,
     &     raynew,amplnew,phasenew,anglenew,tnew,rtnew,flag,error,
     &     MAX_INT,MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer MAX_INT, MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID
      real ray1(N_PARAM),ray2(N_PARAM),tray1,tray2,rtray1(N_RTPARAM),rtray2(N_RTPARAM)
      real ampl1,phase1,angle1,ampl2,phase2,angle2
      integer actlayer,veltype,actint
      integer npoints_int(MAX_INT)
      real x_int(MAXPOINTS_INT,MAX_INT), z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT), c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT), ptos
      real chknorm1, chknorm2
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID),v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real raynew(N_PARAM), amplnew, phasenew, tnew, anglenew, rtnew(N_RTPARAM)
      integer error, flag

      integer MAX_COUNTS
      parameter (MAX_COUNTS = 30)
      real dzdx, dx, d2zdx2, part, part1, part2,ttmp, pxtmp, pztmp, vel1(6), vel2(6)
      real chknormnew, epsilon
      integer actp, count, i
c      double precision pi, threshold, cosnew, xint_normal, zint_normal

c-------------------------------------------------------------

c      if (N_PARAM_LOCAL .ne. N_PARAM) goto 990

      if (abs(ray1(5)) .lt. 0.00001 .or. abs(ray2(5)) .lt. 0.00001) then
         write(*,*) 'GET_NORMRAY: CAUSTIC!'
         flag = -1
         goto 999
      end if

      epsilon = 0.0001
      count = 0
      flag = 0
c      pi = 3.141592653589
c      threshold = cos(inf_angle*0.01)
c      threshold = 1.0

c      if (abs(rtray1(4)) .gt. threshold) then
      if (abs(chknorm1) .lt. epsilon) then
         do i = 1, N_PARAM
            raynew(i) = ray1(i)
         end do
         amplnew = ampl1
         tnew = tray1
         anglenew = angle1
         phasenew = phase1
         do i = 1, N_RTPARAM
            rtnew(i) = rtray1(i)
         end do
         goto 999
c      else if (abs(rtray2(4)) .gt. threshold) then
      else if (abs(chknorm2) .lt. epsilon) then
         do i = 1, N_PARAM
            raynew(i) = ray2(i)
         end do
         amplnew = ampl2
         tnew = tray2
         anglenew = angle2
         phasenew = phase2
         do i = 1, N_RTPARAM
            rtnew(i) = rtray2(i)
         end do
         goto 999
      end if

c      cosnew = 0.0
      part1 = 0.0
      part2 = 1.0
      actp = 1
      chknormnew = chknorm1
      part = 0.5
      call velocity(ray1(1), ray1(2),nx_grid,nz_grid,
     &     x_grid,z_grid,v_grid,veltype,ptos,actlayer,
     &     vel1,
     &     MAXP_XGRID, MAXP_ZGRID)
      call velocity(ray2(1), ray2(2), nx_grid, nz_grid,
     &     x_grid,z_grid,v_grid,veltype,ptos,actlayer,
     &     vel2,
     &     MAXP_XGRID, MAXP_ZGRID)
c-------------------------------------------------------------
c      do while (cosnew .lt. threshold)
      do while (abs(chknormnew) .gt. epsilon)
         count = count + 1
         if (count .gt. MAX_COUNTS) goto 950
         raynew(1) = ray1(1) + part * (ray2(1) - ray1(1))

         call get_z1(raynew(1),x_int(1,actint),z_int(1,actint),d_int(1,actint),
     &        c_int(1,actint),b_int(1,actint),npoints_int(actint),
     &        raynew(2),dzdx,actp,
     &        MAXPOINTS_INT)

         call paraxial(ray1,tray1,vel1,
     &        raynew,ttmp,N_PARAM)
         pxtmp = raynew(3)
         pztmp = raynew(4)
         call paraxial(ray2,tray2,vel2,
     &        raynew,tnew,N_PARAM)
         raynew(3) = pxtmp + part * (raynew(3) - pxtmp)
         raynew(4) = pztmp + part * (raynew(4) - pztmp)

c         call velonly(raynew(1), raynew(2), nx_grid, nz_grid,
c     &        x_grid, z_grid, v_grid,
c     &        veltype, ptos, actlayer,
c     &        vel,
c     &        MAXP_XGRID, MAXP_ZGRID)
c         zint_normal = 1.0/dsqrt(1.0 + dble(dzdx*dzdx))
c         xint_normal = dble(-dzdx) * zint_normal
c         cosnew = dabs((xint_normal*dble(raynew(3)) + zint_normal*dble(raynew(4))) * dble(vel))

         chknormnew = raynew(3) + dzdx*raynew(4)
         if (chknormnew*chknorm1 .lt. 0.0) then
            part2 = part
         else
            part1 = part
         end if
         part = 0.5 * (part1 + part2)
      end do
c-------------------------------------------------------------
      tnew = ttmp + part * (tnew - ttmp)

      do i = 5, 8
         raynew(i) = ray1(i) + part*(ray2(i) - ray1(i))
      end do
      amplnew = ampl1 + part*(ampl2 - ampl1)
      phasenew = phase1 + part*(phase2 - phase1)
      anglenew = angle1 + part*(angle2 - angle1)

      dx = raynew(1) - x_int(actp,actint)
      d2zdx2 = 6.*d_int(actp,actint)*dx + 2.*c_int(actp,actint)
      call get_rtsect(raynew,dzdx,d2zdx2,
     &     nx_grid,nz_grid,x_grid,z_grid,
     &     v_grid,ptos,veltype,actlayer,
     &     rtnew,
     &     N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)

c-------------------------------------------------------------

      goto 999
 950  tnew = ttmp + part * (tnew - ttmp)
      do i = 5, 6
         raynew(i) = ray1(i) + part*(ray2(i) - ray1(i))
      end do
      amplnew = ampl1 + part*(ampl2 - ampl1)
      phasenew = phase1 + part*(phase2 - phase1)
      anglenew = angle1 + part*(angle2 - angle1)
      dx = raynew(1) - x_int(actp,actint)
      d2zdx2 = 6.*d_int(actp,actint)*dx + 2.*c_int(actp,actint)
      call get_rtsect(raynew,dzdx,d2zdx2,
     &     nx_grid,nz_grid,x_grid,z_grid,
     &     v_grid,ptos,veltype,actlayer,
     &     rtnew,
     &     N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
      write(*,*) "WARNING :"
      write(*,*) "GET_NORMRAY : More than MAX_COUNTS iterations... is there a bug?"
      goto 999
 990  write(*,*) "GET_NORMRAY : Parameter N_PARAM_LOCAL not equal N_PARAM from WFRONT"
      write(*,*) "Change N_PARAM_LOCAL to the amount of N_PARAM"
      error = 22
 999  end



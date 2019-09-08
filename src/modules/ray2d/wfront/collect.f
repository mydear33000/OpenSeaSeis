c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm COLLECT_HALT / Hauptprogramm WFRONT
c Ordnet die Strahlen von SECT* in die entsprechenden Felder HALT*
c Genauer: Schnittpunkte mit der oberen und/oder unteren Grenzflaeche werden
c gesammelt.
c Ausserdem werden an Randstellen (fictious interface, coinciding interfaces) neue Strahlen interpoliert!
c CROSSECT > 0: untere Grenzflaeche
c CROSSECT < 0: obere Grenzflaeche
c CROSSECT = 0: kein Schnittpunkt, Strahl verliess aeusseren Modellrand
c IHALTUP, IHALTDOWN: Index des Halts. =0: Nicht aufsammeln
c

      subroutine collect_halt(sect,tsect,rtsect,amplitude,phase,angle,kmahsect,nextsect,crossect,fictsect,actlayer,
     &     code,actcode,actcodestep,nactcodestep,actcompounds,actncomp,hold,
     &     thehalt,thalt,rthalt,amplhalt,phasehalt,anglehalt,kmahhalt,nexthalt,ficthalt,nrayhalt,ncodehalt,codestepshalt,
     &     inthalt,compoundshalt,ncomphalt,ihaltup,ihaltdown,ihaltcomp,error,
     &     N_PARAM,MAX_RAYS,MAX_HALTS,MAX_CODES,MAX_CODESTEPS,N_RTPARAM)

      implicit none

      integer N_PARAM, MAX_RAYS, MAX_HALTS, MAX_CODES, MAX_CODESTEPS, N_RTPARAM
      real sect(N_PARAM,MAX_RAYS), tsect(MAX_RAYS), rtsect(N_RTPARAM,MAX_RAYS), amplitude(MAX_RAYS)
      real phase(MAX_RAYS), angle(MAX_RAYS)
      integer nextsect(0:MAX_RAYS),crossect(MAX_RAYS),kmahsect(MAX_RAYS),fictsect(MAX_RAYS)
      integer ihaltup,ihaltdown,actlayer,ihaltcomp
      integer code(0:MAX_CODESTEPS,MAX_CODES),actcode,actcodestep(MAX_CODESTEPS),nactcodestep
      integer hold(0:MAX_CODESTEPS,MAX_CODES),actncomp,actcompounds(MAX_CODESTEPS)

c TEMP
c      integer N_PARAM_HALT
c      parameter( N_PARAM_HALT=7 )
c END TEMP
      real thehalt(N_PARAM,MAX_RAYS,MAX_HALTS), thalt(MAX_RAYS,MAX_HALTS)
      real rthalt(N_RTPARAM,MAX_RAYS,MAX_HALTS), amplhalt(MAX_RAYS,MAX_HALTS)
      real phasehalt(MAX_RAYS,MAX_HALTS), anglehalt(MAX_RAYS,MAX_HALTS)
      integer nexthalt(0:MAX_RAYS,MAX_HALTS), kmahhalt(MAX_RAYS,MAX_HALTS), ficthalt(MAX_RAYS,MAX_HALTS)
      integer nrayhalt(MAX_HALTS),ncodehalt(MAX_HALTS)
      integer codestepshalt(MAX_HALTS),inthalt(MAX_HALTS)
      integer ncomphalt(MAX_HALTS),compoundshalt(MAX_CODESTEPS,MAX_HALTS),error

      integer iray,iraydown,irayup,i
      integer int1,int2,next


      iraydown = 0
      irayup = 0
      iray = 0
      next = abs(nextsect(iray))
      int1 = actlayer
      int2 = actlayer + 1
c Muessen die Schnittpunkte des oberen und/oder unteren Interfaces gesammelt werden?
      if (ihaltdown .ne. 0) then
         nexthalt(0,ihaltdown) = -1
         if (ihaltup .ne. 0) nexthalt(0,ihaltup) = -1
      else if (ihaltup .ne. 0) then
         nexthalt(0,ihaltup) = -1
      else
c Keine Aufsammlung, nichts passiert:
         goto 999
      end if      
c=====================================================================================================
c Durchlauf aller Schnittpunkte
c
      do while (next .ne. 0)
         iray = next
         next = abs(nextsect(iray))
c         if (crossect(iray) .eq. 0) then
c            do while (crossect(iray) .eq. 0 .and. next .ne. 0)
c               iray = next
c               next = abs(nextsect(iray))
c            end do
c            if (ihaltdown .ne. 0) nexthalt(iraydown,ihaltdown) = -(iraydown+1)
c            if (ihaltup .ne. 0) nexthalt(irayup,ihaltup) = -(irayup+1)
c         end if

c-------------------------------------
c Schnitt mit der unteren Grenzflaeche:
c
         if (crossect(iray) .gt. 0 .and. ihaltdown .ne. 0) then
c Richtiges Setzen von NEXTHALT:
            if (nexthalt(iraydown,ihaltdown) .eq. iray) then
               nexthalt(iraydown,ihaltdown) = iraydown + 1
            else
               nexthalt(iraydown,ihaltdown) = -(iraydown+1)
            end if

            iraydown = iraydown + 1
            do i = 1, N_PARAM
               thehalt(i,iraydown,ihaltdown) = sect(i,iray)
            end do
            do i = 1, N_RTPARAM
               rthalt(i,iraydown,ihaltdown) = rtsect(i,iray)
            end do
            thalt(iraydown,ihaltdown) = tsect(iray)
            amplhalt(iraydown,ihaltdown) = amplitude(iray)
            phasehalt(iraydown,ihaltdown) = phase(iray)
            anglehalt(iraydown,ihaltdown) = angle(iray)
            kmahhalt(iraydown,ihaltdown) = kmahsect(iray)
            ficthalt(iraydown,ihaltdown) = fictsect(iray)
c            nexthalt(iraydown,ihaltdown) = isign(1,nextsect(iray)) * (iraydown+1)
c temporaeres Setzen:
            nexthalt(iraydown,ihaltdown) = nextsect(iray)
c-------------------------------------------------------------------------------------
c Schnitt mit der oberen Grenzflaeche:
c
         else if (crossect(iray) .lt. 0 .and. ihaltup .ne. 0) then
c Richtiges Setzen von NEXTHALT:
            if (nexthalt(irayup,ihaltup) .eq. iray) then
               nexthalt(irayup,ihaltup) = irayup + 1
            else
               nexthalt(irayup,ihaltup) = -(irayup+1)
            end if

            irayup = irayup + 1
            do i = 1, N_PARAM
               thehalt(i,irayup,ihaltup) = sect(i,iray)
            end do
            do i = 1, N_RTPARAM
               rthalt(i,irayup,ihaltup) = rtsect(i,iray)
            end do
            thalt(irayup,ihaltup) = tsect(iray)
            amplhalt(irayup,ihaltup) = amplitude(iray)
            phasehalt(irayup,ihaltup) = phase(iray)
            anglehalt(irayup,ihaltup) = angle(iray)
            kmahhalt(irayup,ihaltup) = kmahsect(iray)
            ficthalt(irayup,ihaltup) = fictsect(iray)
c            nexthalt(irayup,ihaltup) = isign(1,nextsect(iray)) * (irayup+1)
c temporaeres Setzen
            nexthalt(irayup,ihaltup) = nextsect(iray)
         end if
      end do
c
c Ende Durchlauf aller Schnittpunkte
c=================================================================================
c Setzen des Codes etc.:
c

c Untere Grenzflaeche:
c
      if (iraydown .ne. 0) then
         if (iraydown .gt. MAX_RAYS) goto 980
         nexthalt(iraydown,ihaltdown) = 0
         nrayhalt(ihaltdown) = iraydown
         inthalt(ihaltdown) = int2
         ncomphalt(ihaltdown) = actncomp
         do i = 1, actncomp
            compoundshalt(i,ihaltdown) = actcompounds(i)
         end do
c--------------
c Fall: erster Codeschritt; source liegt innerhalb Schicht; code wird spaeter nochmal gebraucht:
         if (actcodestep(nactcodestep) .eq. 1 .and. hold(0,actcode) .ne. 0) then
            ncodehalt(ihaltdown) = hold(0,actcode)
            codestepshalt(ihaltdown) = 1
c-------------
c Fall: normal
         else
            if (ihaltdown .eq. ihaltcomp) then
               ncomphalt(ihaltdown) = ncomphalt(ihaltdown) + 1
               compoundshalt(ncomphalt(ihaltdown),ihaltdown) = actcodestep(nactcodestep)
               codestepshalt(ihaltdown) = actcodestep(nactcodestep) + 1
            else
               codestepshalt(ihaltdown) = actcodestep(nactcodestep)
            end if
c--------------
c compound element; code wird spaeter nochmal gebraucht:
            if (hold(actcodestep(nactcodestep),actcode) .gt. 1) then
               ncodehalt(ihaltdown) = hold(actcodestep(nactcodestep),actcode)
            else
               ncodehalt(ihaltdown) = actcode
            end if
         end if
      else
         ihaltdown = 0
      end if

c------------------------------------------------------------------------------------------------------
c Obere Grenzflaeche
c
      if (irayup .ne. 0) then
         if (irayup .gt. MAX_RAYS) goto 980
         nexthalt(irayup,ihaltup) = 0
         nrayhalt(ihaltup) = irayup
         inthalt(ihaltup) = int1
         ncomphalt(ihaltup) = actncomp
         do i = 1, actncomp
            compoundshalt(i,ihaltup) = actcompounds(i)
         end do
c--------------
c Fall: erster Codeschritt; source liegt innerhalb Schicht; code wird spaeter nochmal gebraucht:
         if (actcodestep(nactcodestep) .eq. 1 .and. hold(0,actcode) .ne. 0 .and. code(0,actcode) .gt. 0) then
            ncodehalt(ihaltup) = hold(0,actcode)
            codestepshalt(ihaltup) = 1
c-------------
c Fall: normal
         else
c--------------
c compound element; code wird spaeter nochmal gebraucht:
            if (hold(actcodestep(nactcodestep),actcode) .gt. 1) then
               ncodehalt(ihaltup) = hold(actcodestep(nactcodestep),actcode)
            else
               ncodehalt(ihaltup) = actcode
            end if
c Falls der letzte Strahlenweg ein compound element war...
            if (ihaltup .eq. ihaltcomp) then
               ncomphalt(ihaltup) = ncomphalt(ihaltup) + 1
               compoundshalt(ncomphalt(ihaltup),ihaltup) = actcodestep(nactcodestep)
               codestepshalt(ihaltup) = actcodestep(nactcodestep) + 1
            else
               codestepshalt(ihaltup) = actcodestep(nactcodestep)
            end if
         end if
      else
         ihaltup = 0
      end if

c 1000 format(/"WARNING: Interfaces no.",i2," and no.",i2," seem to be too close to each other",/,
c     &     "         near the point (",f7.3,",",f7.3,"). The actual resolution can't manage it.",/,
c     &     "         Reduce the input variable 'spread_max' (maximal distance between rays."/)
c 1010 format(/,"COLLECT_HALT: Write error in file no.",i2,/)
 1020 format(/,"COLLECT_HALT: Too many rays!",/," Enlarge PARAMETER 'MAX_RAYS' in WFRONT.F",/)

      goto 999
 980  write(*,1020)
      error = 2
      goto 999

 999  end

c******************************************************************************************************
c Unterprogramm COLLECT_BORE / Hauptprogramm WFRONT
c Ordnet die Strahlen von SECT* in ein entsprechendes Feld HALT*
c Genauer: Schnittpunkte mit einem Bohrloch werden gesammelt.
c TSECT = 0.0 : kein Schnittpunkt
c sonst       : Schnittpunkt
c
      subroutine collect_bore(leftsect,tleftsect,kmahleftsect,rightsect,trightsect,kmahrightsect,
     &     nextsect,ibore,amplitude,phase,angle,
     &     actcode,actcodestep,nactcodestep,actcompounds,actncomp,
     &     halt,thalt,amplhalt,phasehalt,anglehalt,kmahhalt,nexthalt,nrayhalt,ncodehalt,codestepshalt,
     &     compoundshalt,ncomphalt,flag,
     &     N_PARAM,MAX_RAYS,MAX_CODESTEPS,n_boreholes)

      implicit none

      integer N_PARAM,MAX_RAYS,MAX_CODESTEPS, n_boreholes, ibore
      real leftsect(N_PARAM,n_boreholes,MAX_RAYS), tleftsect(n_boreholes,MAX_RAYS), amplitude(MAX_RAYS)
      real rightsect(N_PARAM,n_boreholes,MAX_RAYS), trightsect(n_boreholes,MAX_RAYS)
      real angle(MAX_RAYS), phase(MAX_RAYS)
      integer nextsect(0:MAX_RAYS),kmahleftsect(n_boreholes,MAX_RAYS),kmahrightsect(n_boreholes,MAX_RAYS)
      integer actcode,actcodestep(MAX_CODESTEPS),nactcodestep
      integer actncomp,actcompounds(MAX_CODESTEPS)
      real halt(N_PARAM,MAX_RAYS), thalt(MAX_RAYS), amplhalt(MAX_RAYS), phasehalt(MAX_RAYS), anglehalt(MAX_RAYS)
      integer nexthalt(0:MAX_RAYS), kmahhalt(MAX_RAYS)
      integer nrayhalt, ncodehalt, codestepshalt
      integer ncomphalt, compoundshalt(MAX_CODESTEPS), flag

      integer iray, iraynew, i, next

      iraynew = 0
c-----------------------------------------------------------------      
c Durchlauf aller evt. Schnittpunkte
c
c Einsaetze von links einsammeln:
c
      iray = 0
      next = 1

      do while (nextsect(iray) .ne. 0)
         if (isign(1,nextsect(iray)) .eq. -1) then
            next = -1
         end if
         iray = abs(nextsect(iray))
         if (tleftsect(ibore,iray) .gt. 0.0) then
            nexthalt(iraynew) = next*(iraynew+1)
            iraynew = iraynew + 1            
            do i = 1, N_PARAM
               halt(i,iraynew) = leftsect(i,ibore,iray)
            end do
            thalt(iraynew) = tleftsect(ibore,iray)
            amplhalt(iraynew) = amplitude(iray)
            phasehalt(iraynew) = phase(iray)
            anglehalt(iraynew) = angle(iray)
            kmahhalt(iraynew) = kmahleftsect(ibore,iray)
            next = isign(1,nextsect(iray))
c         else if (trightsect(ibore,iray) .ne. 0.0) then
c            next = -1
         end if
      end do

c Einsaetze von rechts einsammeln:
c
      iray = 0
      next = -1

      do while (nextsect(iray) .ne. 0)
         if (isign(1,nextsect(iray)) .eq. -1) then
            next = -1
         end if
         iray = abs(nextsect(iray))
         if (trightsect(ibore,iray) .gt. 0.0) then
            nexthalt(iraynew) = next*(iraynew+1)
            iraynew = iraynew + 1            
            do i = 1, N_PARAM
               halt(i,iraynew) = rightsect(i,ibore,iray)
            end do
            thalt(iraynew) = trightsect(ibore,iray)
            amplhalt(iraynew) = amplitude(iray)
            phasehalt(iraynew) = phase(iray)
            anglehalt(iraynew) = angle(iray)
            kmahhalt(iraynew) = kmahrightsect(ibore,iray)
            next = isign(1,nextsect(iray))
c         else if (tleftsect(ibore,iray) .ne. 0.0) then
c            next = -1
         end if
      end do
c
c Ende Durchlauf aller Schnittpunkte
c---------------------------------------------------------------------------------

c Setzen des Codes etc.:
      if (iraynew .ne. 0) then
         nexthalt(iraynew) = 0
         nrayhalt = iraynew
         ncomphalt = actncomp
         do i = 1, actncomp
            compoundshalt(i) = actcompounds(i)
         end do
         ncodehalt = actcode
         codestepshalt = actcodestep(nactcodestep)
c Vorsicht: compound elements koennen innerhalb des Borhlochs variieren!
c Deshalb wird die letzte Anzahl genommen, ohne die evt. neue Aenderung zu beruecksichtigen
         flag = 0
      else
         flag = -1
      end if

c      do i = 1, nrayhalt
c         write(45,*) halt(1,i),halt(2,i)
c         write(46,*) halt(1,i),thalt(i)
c         write(47,*) halt(1,i),amplhalt(i)
c         if (nexthalt(i) .lt. 0) then
c            write(45,'()')
c            write(46,'()')
c            write(47,'()')
c         end if
c      end do
c      write(45,'()')
c      write(45,'()')
c      write(46,'()')
c      write(46,'()')
c      write(47,'()')
c      write(47,'()')

      end



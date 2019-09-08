c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm GET_ACTHALT / Hauptprogramm WFRONT
c Wird aufgerufen von WFRONT
c Aufruf der Unterprogramme GET_Z
c Bestimmt den naechsten Codeschritt...:
c ACTCODESTEP und ACTCODE werden bestimmt.
c ACTHALT: Index des Halts, von dem aus weiterpropagiert wird
c ACTHALT = 0: Propagation von der Quelle
c FLAG = -1: letzter Code, Programm beenden!
c FLAG =  0: sonst
c
      subroutine get_acthalt(code,codesteps,ncodehalt,ncodes,ihaltup,ihaltdown,
     &     ncomphalt,compoundshalt,hold,codestepshalt,flag_compounds,all_reclines,
     &     actcode,actcodestep,nactcodestep,acthalt,actualhalt,ihaltnorm,ihaltcomp,actncomp,
     &     actcompounds,icompout,flag,
     &     MAX_HALTS, MAX_CODES, MAX_CODESTEPS,MAX_RECLINES)

      implicit none

      integer MAX_HALTS, MAX_CODES, MAX_CODESTEPS, MAX_RECLINES
      integer code(0:MAX_CODESTEPS,MAX_CODES),codesteps(MAX_CODES), hold(0:MAX_CODESTEPS,MAX_CODES)
      integer ncodehalt(MAX_HALTS),codestepshalt(MAX_HALTS), nactcodestep, flag_compounds
      integer ncodes,ihaltup,ihaltdown
      integer ncomphalt(MAX_HALTS),compoundshalt(MAX_CODESTEPS,MAX_HALTS)
      integer actcode, actcodestep(MAX_CODESTEPS), actualhalt(MAX_CODESTEPS), acthalt, ihaltnorm, ihaltcomp
      integer flag
      integer actncomp,actcompounds(MAX_CODESTEPS),icompout(MAX_RECLINES),all_reclines

      integer nactcodesteptmp
      integer i,j,ihalt

      flag = 0
      nactcodesteptmp = nactcodestep
 1    if (actcode .eq. ncodes) then
         do while (actcodestep(nactcodestep) .eq. codesteps(ncodes))
            if (nactcodestep .eq. 1) goto 900
            nactcodestep = nactcodestep - 1
         end do
      end if

c-----------------------------------------------------------
c Fall: Dieser Code ist abgearbeitet, der naechste ist dran
c
      if (actcodestep(nactcodestep) .eq. codesteps(actcode)) then
c die letzte Verzweigung aufgrund eines compound elements ist dran...
         if (nactcodestep .ne. 1) then
            nactcodestep = nactcodestep - 1
            goto 1
         end if
         actcode = actcode + 1
         do i = 1, all_reclines
            icompout(i) = 1
         end do
         actcodestep(1) = -2
c Ist ein Teilstueck dieses Strahlencodes bereits berechnet worden?
         do ihalt = 1, MAX_HALTS
            if (codestepshalt(ihalt) .ne. 0) then
               j = 0
               if (codestepshalt(ihalt) .le. codesteps(actcode)) then
c Pruefe, wieviele Steps uebereinstimmen:
                  do while (j .le. codestepshalt(ihalt) .and. 
     &                 code(j,ncodehalt(ihalt)) .eq. code(j,actcode))
                     j = j + 1
                  end do
               end if
               j = j - 1
               if (j .gt. actcodestep(1)) then
                  if (j-actcodestep(1) .eq. 1 .and. ncodehalt(ihalt) .eq. actcode) then
                     nactcodestep = nactcodestep + 1
                     actcodestep(nactcodestep) = j
                     actualhalt(nactcodestep) = ihalt
                  else
                     nactcodestep = 1
                     actcodestep(1) = j
                     actualhalt(1) = ihalt
                  end if
               else if (j .eq. actcodestep(1)) then
                  nactcodestep = nactcodestep + 1
                  actcodestep(nactcodestep) = j
                  actualhalt(nactcodestep) = ihalt
c letzter Schritt von nactcodestep war compound, und hold(,) > 1:
               else if (actcodestep(1) - j .eq. 1 .and. ncodehalt(actualhalt(1)) .eq. actcode) then
                  nactcodestep = nactcodestep + 1
                  actcodestep(nactcodestep) = j
                  actualhalt(nactcodestep) = ihalt                  
               end if
            end if
         end do
c Keine Uebereinstimmung gefunden... vom Anfang:
         if (actcodestep(1) .lt. 0) then
            acthalt = 0
            actcodestep(1) = 1
            actncomp = 0
         else
            acthalt = actualhalt(nactcodestep)
            actcodestep(nactcodestep) = actcodestep(nactcodestep) + 1
            actncomp = ncomphalt(acthalt)
            do i = 1, actncomp
               actcompounds(i) = compoundshalt(i,acthalt)
            end do
         end if
c-----------------------------------------------------------
c Fall: naechster Code-Schritt ist dran
c
      else if (nactcodestep .eq. nactcodesteptmp) then
c 1. Fall: compound element bei ihaltdown:
         if (ihaltcomp .ne. 0 .and. ihaltdown .eq. ihaltcomp .and. codestepshalt(ihaltdown) .ne. 0) then
            if (codestepshalt(ihaltdown) .eq. codesteps(actcode)) then
               ihaltcomp = 0
               if (hold(codesteps(actcode),actcode) .eq. 0) then
                  codestepshalt(ihaltdown) = 0
               end if
               goto 1
            end if
            if (codestepshalt(ihaltnorm) .ne. 0) then
               nactcodestep = nactcodestep + 1
            end if
            acthalt = ihaltdown
            actualhalt(nactcodestep) = acthalt
            actcodestep(nactcodestep) = codestepshalt(ihaltdown) + 1
c 2. Fall: compound element bei ihaltup:
         else if (ihaltcomp .ne. 0 .and. ihaltup .eq. ihaltcomp .and. codestepshalt(ihaltup) .ne. 0) then
            if (codestepshalt(ihaltup) .eq. codesteps(actcode)) then
               ihaltcomp = 0
               if (hold(codesteps(actcode),actcode) .eq. 0) then
                  codestepshalt(ihaltup) = 0
               end if
               goto 1
            end if
            if (codestepshalt(ihaltnorm) .ne. 0) then
               nactcodestep = nactcodestep + 1
            end if
            acthalt = ihaltup
            actualhalt(nactcodestep) = acthalt
            actcodestep(nactcodestep) = codestepshalt(ihaltup) + 1
c 3. Fall: normal weiter bei ihaltdown
         else if (ihaltdown .eq. ihaltnorm .and. codestepshalt(ihaltdown) .ne. 0) then
            actcodestep(nactcodestep) = actcodestep(nactcodestep) + 1
            acthalt = ihaltdown
            actualhalt(nactcodestep) = acthalt
c 4. Fall: normal weiter bei ihaltup
         else if (ihaltup .eq. ihaltnorm .and. codestepshalt(ihaltup) .ne. 0) then
            actcodestep(nactcodestep) = actcodestep(nactcodestep) + 1
            acthalt = ihaltup
            actualhalt(nactcodestep) = acthalt
c 5. Fall: Kein Halt mit dem gesuchten code wurde gefunden... die Strahlen muessen versiegt sein
c          (oder die vorherige Verzweigung aufgrund eines compound elements ist dran)
         else
            if (nactcodestep .ne. 1) then
               nactcodestep = nactcodestep - 1
            else
               write(*,*) "Wellenfront versiegt... Code Nr. ",actcode
               actcodestep(nactcodestep) = codesteps(actcode)
            end if
            goto 1
         end if
         actncomp = ncomphalt(acthalt)
         do i = 1, actncomp
            actcompounds(i) = compoundshalt(i,acthalt)
         end do
c----------
c Der letzte (NACTCODESTEP) aufgezeichnete Schritt ist dran...
      else
         actcodestep(nactcodestep) = actcodestep(nactcodestep) + 1
         acthalt = actualhalt(nactcodestep)
         actncomp = ncomphalt(acthalt)
         do i = 1, actncomp
            actcompounds(i) = compoundshalt(i,acthalt)
         end do
      end if

c------------
c Freimachen des aktuellen Halts
c
      if (acthalt .ne. 0 .and. hold(codestepshalt(acthalt),actcode) .eq. 0) then
         codestepshalt(acthalt) = 0
      end if
c-----------------
c Setzen von IHALTNORM und IHALTCOMP
c
      ihaltnorm = 1
      do while (codestepshalt(ihaltnorm) .ne. 0)
         ihaltnorm = ihaltnorm + 1
         if (ihaltnorm .gt. MAX_HALTS) then
            write(*,1010)
            ihaltnorm = 1
            codesteps(ihaltnorm) = 0
         end if
      end do
      actualhalt(nactcodestep) = ihaltnorm

c Falls eine Reflexion folgt oder dieser code-Anfang nochmal fuer einen anderen
c up/down-flag gebraucht wird (nur, wenn Quelle innerhalb Schicht liegt)

      if ((actcodestep(nactcodestep) .ne. codesteps(actcode) .and.
     &     code(actcodestep(nactcodestep),actcode) .eq. code(actcodestep(nactcodestep)+1,actcode)) .or.
     &     (actcodestep(nactcodestep) .eq. 1 .and. hold(0,actcode) .ne. 0) .or.
     &     hold(actcodestep(nactcodestep),actcode) .gt. 1 ) then
         if (flag_compounds .ne. 0) then
            ihaltcomp = 0
         else
            ihaltcomp = ihaltnorm + 1
            if (ihaltcomp .gt. MAX_HALTS) then
               write(*,1010)
               ihaltcomp = 1
               codesteps(ihaltcomp) = 0
            end if
            do while (codestepshalt(ihaltcomp) .ne. 0)
               ihaltcomp = ihaltcomp + 1
               if (ihaltcomp .gt. MAX_HALTS) then
                  write(*,1010)
                  ihaltcomp = ihaltnorm + 1
                  codesteps(ihaltcomp) = 0
               end if
            end do
         end if
      else
         ihaltcomp = 0
      end if

      goto 999

 1010 format("WARNING!!"/"GET_ACTHALT : Too many halts needed!"/"Enlarge parameter MAX_HALTS in WFRONT")

 900  flag = -1
      goto 999

 999  end


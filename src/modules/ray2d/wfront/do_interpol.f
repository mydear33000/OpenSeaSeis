c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm DO_INTERPOL / Hauptprogramm WFRONT
c letzte Aenderung 6.02.97  /  getestet
c Wird nur aufgerufen von WFRONT
c Aufruf des Unterprogramms GET_Z
c Interpoliert neue Strahlen, falls Abstand oder Winkel zwischen zwei benachbarten
c Strahlen zu gross werden. Es werden eventuell mehrere Strahlen gleichzeitig
c neu angesetzt.
c
      subroutine do_interpol(actlayer,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &     spread_max, cos_min, angle_max, next,
     &     ray,amplitude,phase,angle,kmah,
     &     crossect,fictsect,kmahsect,nextray,nextsect,irayold,iray,nray,
     &     nextboresect,tleftboresect,trightboresect,error,
     &     N_PARAM,MAX_INT,MAXPOINTS_INT,n_boreholes,MAX_RAYS)

      implicit none

      integer N_PARAM,MAX_INT,MAXPOINTS_INT,MAX_RAYS
      real ray(N_PARAM,MAX_RAYS),amplitude(MAX_RAYS)
      real phase(MAX_RAYS),angle(MAX_RAYS),spread_max,cos_min,angle_max
      integer actlayer, npoints_int(MAX_INT)
      real x_int(MAXPOINTS_INT,MAX_INT), z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT), c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT)
      integer crossect(MAX_RAYS),nextray(0:MAX_RAYS),nextsect(0:MAX_RAYS),kmah(MAX_RAYS)
      integer fictsect(MAX_RAYS), kmahsect(MAX_RAYS)
      integer n_boreholes,nextboresect(0:MAX_RAYS)
      real tleftboresect(n_boreholes,MAX_RAYS), trightboresect(n_boreholes,MAX_RAYS)
      integer irayold,iray,next,nray,error

      integer n_interpolations,i,n,nraynew,iraysave
      real epsilon,dx,dz,distance,p1,p2,p1p2,p1_x_p2,lamda,my,xvirtual,zvirtual,r1,r2,pnew,ptmp
      real part,delta_part,actangle

      real zint,extrasign
      integer actp,ibore

c Sicherheitsepsilon, um Fall 1/0 abzufangen
      epsilon = 0.000001
      iraysave = iray

      dx = ray(1,next) - ray(1,iray)
      dz = ray(2,next) - ray(2,iray)
      distance = dx*dx + dz*dz
      p1 = sqrt(ray(3,iray)*ray(3,iray) + ray(4,iray)*ray(4,iray))
      p2 = sqrt(ray(3,next)*ray(3,next) + ray(4,next)*ray(4,next))
      p1p2 = ray(3,iray)*ray(3,next) + ray(4,iray)*ray(4,next)

c----------------------------
c Nur Interpolation, wenn Abstand zwischen Strahlen groesser spread_max oder
c Winkel zwischen p1 und p2 kleiner cos_min (dies ist cos(max. Winkel))
c
c

      if (distance .lt. spread_max .and. abs(p1p2/(p1*p2)) .gt. cos_min) goto 900

      p1_x_p2 = ray(4,iray)*ray(3,next) - ray(3,iray)*ray(4,next)
c----------------------------
c Strahlen sind nicht parallel (sonst einfache Mittelung)
c
      if (abs(p1_x_p2) .gt. epsilon) then
         lamda = (dz*ray(3,next) - dx*ray(4,next))/p1_x_p2
c----------------------------
c Keine Interpolation bei Konvergenz, und Strahlen nah beieinander:
c
         if (lamda .gt. 0.0) then
            if (distance .lt. spread_max) goto 900
            extrasign = -1.0
         else
            extrasign = 1.0
         end if
         my = (dz*ray(3,iray) - dx*ray(4,iray))/p1_x_p2
         xvirtual = ray(1,iray) + lamda*ray(3,iray)
         zvirtual = ray(2,iray) + lamda*ray(4,iray)
         r1 = abs(lamda*p1)
         r2 = abs(my*p2)
      end if

      n_interpolations = int(sqrt(distance/spread_max))
c pruefen, ob nicht IMAX
      actangle = p1p2/(p1*p2)
      if (actangle .lt. 0.0) then
         if (distance .lt. spread_max) goto 900
         actangle = 0.0
      else if (actangle .gt. 1.0) then
         actangle = 0.0
      else
         actangle = acos(actangle)
      end if
c hier geaendert:
      n_interpolations = max(n_interpolations,int(actangle/angle_max))
      if (n_interpolations .eq. 0) n_interpolations = 1
      delta_part = 1.0 / float(n_interpolations+1)

c========================================================================================
c 'n_interpolations' Interpolationen, alle zwischen Originalstrahlen 'IRAYSAVE' und 'NEXT':
c
      do 10 n = 1, n_interpolations
         part = delta_part*float(n)
         nraynew = nray + 1
         if (nraynew .gt. MAX_RAYS) goto 990
c Versuch einer neuen Interpolation der Slowness:

c lineare Interpolation des Betrages der Slowness
         pnew = p1 + part*(p2 - p1)

c lineare Interpolation der Slowness ...es geht aber noch weiter...
         ray(3,nraynew) = ray(3,iraysave) + part*(ray(3,next) - ray(3,iraysave))
         ray(4,nraynew) = ray(4,iraysave) + part*(ray(4,next) - ray(4,iraysave))

c ...naemlich damit, dass der Einheitsvektor der neuen Slowness genommen wird und
c mit dem neuen Betrag der Slowness multipliziert wird 
         ptmp = pnew / sqrt(ray(3,nraynew)*ray(3,nraynew) + ray(4,nraynew)*ray(4,nraynew))
         ray(3,nraynew) = ray(3,nraynew) * ptmp
         ray(4,nraynew) = ray(4,nraynew) * ptmp

c Bestimmung des Ortes
         if (abs(p1_x_p2) .gt. epsilon) then
            ray(1,nraynew) = xvirtual + extrasign*(r1 + part*(r2 - r1))*(ray(3,nraynew)/pnew)
            ray(2,nraynew) = zvirtual + extrasign*(r1 + part*(r2 - r1))*(ray(4,nraynew)/pnew)
         else
            do i = 1, 2
               ray(i,nraynew) = ray(i,iraysave) + part*(ray(i,next) - ray(i,iraysave))
            end do
         end if

c--------------
c Pruefe auf Schnittpunkt mit oberer Grenzflaeche:
         call get_z(ray(1,nraynew),x_int(1,actlayer),z_int(1,actlayer),d_int(1,actlayer),c_int(1,actlayer),
     &        b_int(1,actlayer),npoints_int(actlayer),
     &        zint,actp,
     &        MAXPOINTS_INT)
         if (ray(2,nraynew) .le. zint) goto 10

c--------------
c Pruefe auf Schnittpunkt mit unterer Grenzflaeche:
         call get_z(ray(1,nraynew),x_int(1,actlayer+1),z_int(1,actlayer+1),d_int(1,actlayer+1),
     &        c_int(1,actlayer+1),b_int(1,actlayer+1),npoints_int(actlayer+1),
     &        zint,actp,
     &        MAXPOINTS_INT)
         if (ray(2,nraynew) .ge. zint) goto 10

         nray = nraynew

         ray(5,nray) = ray(5,iraysave) + part*(ray(5,next) - ray(5,iraysave))
         ray(6,nray) = ray(6,iraysave) + part*(ray(6,next) - ray(6,iraysave))
         ray(7,nray) = ray(7,iraysave) + part*(ray(7,next) - ray(7,iraysave))
         ray(8,nray) = ray(8,iraysave) + part*(ray(8,next) - ray(8,iraysave))
      
         amplitude(nray) = amplitude(iraysave) + part*(amplitude(next) - amplitude(iraysave))
         phase(nray) = phase(iraysave) + part*(phase(next) - phase(iraysave))
         angle(nray) = angle(iraysave) + part*(angle(next) - angle(iraysave))
         kmah(nray) = kmah(iraysave)

         crossect(nray) = 0
         fictsect(nray) = 0
         kmahsect(nray) = 0
         nextray(nray) = next
         nextsect(nray) = next
         nextray(iray) = nray
         nextsect(iray) = nray
         nextboresect(nray) = next
         nextboresect(iray) = nray
         do ibore = 1, n_boreholes
            if (tleftboresect(ibore,iray) .ne. 0.0 .or. tleftboresect(ibore,next) .ne. 0.0) then
               tleftboresect(ibore,nray) = -1.0
            else
               tleftboresect(ibore,nray) = 0.0
            end if
            if (trightboresect(ibore,iray) .ne. 0.0 .or. trightboresect(ibore,next) .ne. 0.0) then
               trightboresect(ibore,nray) = -1.0
            else
               trightboresect(ibore,nray) = 0.0
            end if
         end do
         iray = nray

 10   continue
c
c========================================================================================
 
 900  irayold = iraysave

      goto 999
 990  write(*,*) "DO_INTERPOL: Too many rays!"
      write(*,*) "Enlarge PARAMETER 'MAX_RAYS' in WFRONT"
      error = 1
 999  end


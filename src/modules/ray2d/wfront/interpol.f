c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Interpolation zwischen benachbarten Strahlen
c wird aufgerufen von DO_INTERPOL, INTERPOL_EXTRA
c pnew,qdynnew,pdynnew werden durch einfache Mittelung bestimmt.
c xnew,znew werden durch Winkelmittelung sowie Mittelung
c des Abstandes zum virtuellen Schnittpunkt berechnet.
c FLAG = 0 : Interpolation erfolgt
c FLAG = -1: keine Interpolation

      subroutine interpolate(ray1,amplitude1,phase1,kmah1,angle1,
     &     ray2,amplitude2,phase2,angle2,
     &     spread_max, cos_min, cross,
     &     raynew,amplitudenew,phasenew,kmahnew,anglenew,flag,
     &     N_PARAM)

      implicit none

      integer N_PARAM, kmah1, kmahnew, cross
      real ray1(N_PARAM),ray2(N_PARAM),amplitude1,amplitude2,phase1,phase2,angle1,angle2
      real spread_max, cos_min

      real raynew(N_PARAM),amplitudenew,phasenew,anglenew
      integer flag

      real epsilon,dx,dz,d,p1,p2,p1p2,p1_x_p2,lamda,my,xvirtual,zvirtual,r1,r2,pnew

c Sicherheitsepsilon, um Fall 1/0 abzufangen
      epsilon = 0.000001

      flag = 0
      dx = ray2(1) - ray1(1)
      dz = ray2(2) - ray1(2)
      d = dx*dx + dz*dz
      p1 = sqrt(ray1(3)*ray1(3) + ray1(4)*ray1(4))
      p2 = sqrt(ray2(3)*ray2(3) + ray2(4)*ray2(4))
      p1p2 = ray1(3)*ray2(3) + ray1(4)*ray2(4)

c Nur Interpolation, wenn Abstand zwischen Strahlen groesser spread_max oder
c Winkel zwischen p1 und p2 kleiner cos_min (dies ist cos(max.Winkel))
      if (d .lt. spread_max .and. abs(p1p2/(p1*p2)) .gt. cos_min) then
         if (cross .eq. 0) then
            flag = -1
            goto 999
         end if
      end if

      p1_x_p2 = ray1(4)*ray2(3) - ray1(3)*ray2(4)

c neue Slowness :
      raynew(3) = 0.5*(ray1(3) + ray2(3))
      raynew(4) = 0.5*(ray1(4) + ray2(4))
      
c Berechnung der neuen x-, z-Koordinate :

c Sind die Strahlen annaehernd parallel, wird einfach gemittelt
      if (abs(p1_x_p2) .lt. epsilon) then
         raynew(1) = 0.5*(ray1(1) + ray2(1))
         raynew(2) = 0.5*(ray1(2) + ray2(2))
      else
         lamda = (dz*ray2(3) - dx*ray2(4))/p1_x_p2
c Bei Konvergenz der Strahlen keine Interpolation
c Allerdings wird in INTERPOL_EXTRA auch diese Interpolation gebraucht:
         if (lamda .gt. 0.0) then
            if (cos_min .gt. 1.5) then
               my = (dz*ray1(3) - dx*ray1(4))/p1_x_p2
               pnew = sqrt(raynew(3)*raynew(3) + raynew(4)*raynew(4))
               xvirtual = ray1(1) + lamda*ray1(3)
               zvirtual = ray1(2) + lamda*ray1(4)
               r1 = abs(lamda*p1)
               r2 = abs(my*p2)
               raynew(1) = xvirtual - 0.5*(r1 + r2)*(raynew(3)/pnew)
               raynew(2) = zvirtual - 0.5*(r1 + r2)*(raynew(4)/pnew)
            else
               flag = -1
               goto 999
            end if
         else
            my = (dz*ray1(3) - dx*ray1(4))/p1_x_p2
            pnew = sqrt(raynew(3)*raynew(3) + raynew(4)*raynew(4))
            xvirtual = ray1(1) + lamda*ray1(3)
            zvirtual = ray1(2) + lamda*ray1(4)
            r1 = abs(lamda*p1)
            r2 = abs(my*p2)
            raynew(1) = xvirtual + 0.5*(r1 + r2)*(raynew(3)/pnew)
            raynew(2) = zvirtual + 0.5*(r1 + r2)*(raynew(4)/pnew)
         end if
      end if

c--------------------------------------------------------------------------
c neue dyn. Parameter, Amplitude, Phase, KMAH-Index, initialer Winkel
c werden linear interpoliert
c
      raynew(5) = 0.5*(ray1(5) + ray2(5))
      raynew(6) = 0.5*(ray1(6) + ray2(6))
      amplitudenew = 0.5 * (amplitude1 + amplitude2)
      phasenew = 0.5 * (phase1 + phase2)
      kmahnew = kmah1
      anglenew = 0.5*(angle1 + angle2)
c Reflektionspunkt wird linear interpoliert
      raynew(7) = 0.5*(ray1(7) + ray2(7))
      raynew(8) = 0.5*(ray1(8) + ray2(8))

 999  end

c******************************************************************************************
c Interpolation zwischen benachbarten Strahlen
c wird aufgerufen von GET_CORNERRAYS
c pnew,qdynnew,pdynnew werden durch einfache Mittelung bestimmt.
c xnew,znew werden durch Winkelmittelung sowie Mittelung
c des Abstandes zum virtuellen Schnittpunkt berechnet.

      subroutine interpolate2(ray1,amplitude1,phase1,kmah1,angle1,
     &     ray2,amplitude2,phase2,angle2,part,
     &     raynew,amplitudenew,phasenew,kmahnew,anglenew,
     &     N_PARAM)

      implicit none

      integer N_PARAM, kmah1, kmahnew
      real ray1(N_PARAM),ray2(N_PARAM),amplitude1,amplitude2,phase1,phase2,angle1,angle2,part

      real raynew(N_PARAM),amplitudenew,phasenew,anglenew

      real epsilon,dx,dz,p1,p2,p1_x_p2,lamda,my,xvirtual,zvirtual,r1,r2,pnew

c Sicherheitsepsilon, um Fall 1/0 abzufangen
      epsilon = 0.000001

      dx = ray2(1) - ray1(1)
      dz = ray2(2) - ray1(2)
      p1 = sqrt(ray1(3)*ray1(3) + ray1(4)*ray1(4))
      p2 = sqrt(ray2(3)*ray2(3) + ray2(4)*ray2(4))
      p1_x_p2 = ray1(4)*ray2(3) - ray1(3)*ray2(4)

c neue Slowness :
      raynew(3) = ray1(3) + part*(ray2(3) - ray1(3))
      raynew(4) = ray1(4) + part*(ray2(4) - ray1(4))

c Berechnung der neuen x-, z-Koordinate :

c Sind die Strahlen annaehernd parallel, wird einfach gemittelt
      if (abs(p1_x_p2) .lt. epsilon) then
         raynew(1) = ray1(1) + part*(ray2(1) - ray1(1))
         raynew(2) = ray1(2) + part*(ray2(2) - ray1(2))
      else
         lamda = (dz*ray2(3) - dx*ray2(4))/p1_x_p2
c Bei Konvergenz der Strahlen keine Interpolation
c Allerdings wird in INTERPOL_EXTRA auch diese Interpolation gebraucht:
         if (lamda .gt. 0.0) then
            my = (dz*ray1(3) - dx*ray1(4))/p1_x_p2
            pnew = sqrt(raynew(3)*raynew(3) + raynew(4)*raynew(4))
            xvirtual = ray1(1) + lamda*ray1(3)
            zvirtual = ray1(2) + lamda*ray1(4)
            r1 = abs(lamda*p1)
            r2 = abs(my*p2)
            raynew(1) = xvirtual - (r1 + part*(r2 - r1))*(raynew(3)/pnew)
            raynew(2) = zvirtual - (r1 + part*(r2 - r1))*(raynew(4)/pnew)
         else
            my = (dz*ray1(3) - dx*ray1(4))/p1_x_p2
            pnew = sqrt(raynew(3)*raynew(3) + raynew(4)*raynew(4))
            xvirtual = ray1(1) + lamda*ray1(3)
            zvirtual = ray1(2) + lamda*ray1(4)
            r1 = abs(lamda*p1)
            r2 = abs(my*p2)
            raynew(1) = xvirtual + (r1 + part*(r2 - r1))*(raynew(3)/pnew)
            raynew(2) = zvirtual + (r1 + part*(r2 - r1))*(raynew(4)/pnew)
         end if
      end if

c--------------------------------------------------------------------------
c neue dyn. Parameter, Amplitude, Phase, KMAH-Index, initialer Winkel
c werden linear interpoliert
c
      raynew(5) = ray1(5) + part*(ray2(5) - ray1(5))
      raynew(6) = ray1(6) + part*(ray2(6) - ray1(6))
      amplitudenew = amplitude1 + part * (amplitude2 - amplitude1)
      phasenew = phase1 + part * (phase2 - phase1)
      kmahnew = kmah1
      anglenew = angle1 + part*(angle2 - angle1)

      raynew(7) = ray1(7) + part*(ray2(7) - ray1(7))
      raynew(8) = ray1(8) + part*(ray2(8) - ray1(8))
 
      end


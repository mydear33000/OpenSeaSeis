c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm PARAXIAL / Hauptprogramm WFRONT
c Wird aufgerufen von GET_RAYONINT, INTERPOL_BOUND
c Berechnung der Laufzeit und Slowness im Punkt raynew(1),raynew(2)
c durch paraxiales raytracing.
c
c pnew(i)    = pold(i) + W(i,j)*delta(j)
c tnew       = told + 0.5 * delta(i) * (pnew(i) + pold(i))
c
c delta(i)   = xnew(i) - xold(i)
c w(i,j)     = omega(k,i) * a(k,m) * omega(m,j)
c omega(k,i) = siehe unten...Tangente
c a(k,m)     = siehe unten...
c
c tangent(i) = pold(i) / pold
c 
      subroutine paraxial(rayold,told,vel,
     &     raynew,tnew,N_PARAM)

      implicit none

      integer N_PARAM
      real raynew(N_PARAM),rayold(N_PARAM),tnew,told
      real delta(2),tangent(2)
      real omega(2,2),a(2,2),aux(2,2),w(2,2)
      real vn,vs,v2,p
      real vel(6)

      integer i,j,k

c p koennte man durch 1/vel(1) vertauschen ...Schnelligkeit
      p = sqrt(rayold(3)*rayold(3) + rayold(4)*rayold(4))
      do i = 1, 2
         delta(i) = raynew(i) - rayold(i)
         tangent(i) = rayold(i+2)/p
      end do
      omega(1,1) = tangent(2)
      omega(1,2) = -tangent(1)
      omega(2,1) = tangent(1)
      omega(2,2) = tangent(2)
      vn = vel(2)*tangent(2) - vel(3)*tangent(1)
      vs = vel(2)*tangent(1) + vel(3)*tangent(2)
      v2 = vel(1)*vel(1)

c Idee: Pruefe, ob mit v2 = 1./(p*p) bessere Ergebnisse zustande kommen
c Weiterhin pruefe mit nur          t1 = t1 + rayold(i+2)*delta(i)

      a(1,1) = rayold(6)/rayold(5)
      a(1,2) = -vn/v2
      a(2,1) = -vn/v2
      a(2,2) = -vs/v2

      do 1 i = 1, 2
         do 1 j = 1, 2
            aux(i,j) = 0.0
            do 1 k = 1, 2
               aux(i,j) = aux(i,j) + a(i,k)*omega(k,j)
 1    continue      
      do 2 i = 1, 2
         do 2 j = 1, 2
            w(i,j) = 0.0
            do 2 k = 1, 2
               w(i,j) = w(i,j) + omega(k,i)*aux(k,j)
 2    continue
c---------
c Zeit alt:
c
c      t2 = 0.0
c      t1 = 0.0
c      do 3 i = 1, 2
c         t1 = t1 + rayold(i+2)*delta(i)
c         do 3 j = 1, 2
c            t2 = t2 + delta(i)*delta(j)*w(i,j)
c 3    continue
c      tnew = told + t1 + 0.5*t2

c-------------------
c slowness:
c
      do 4 i = 3, 4
         raynew(i) = rayold(i)
         do 5 j = 3, 4
            raynew(i) = raynew(i) + w(i-2,j-2)*delta(j-2)
 5       continue
 4    continue

c-------------------
c Zeit:
c      
      tnew = told
      do i = 3, 4
         tnew = tnew + 0.5*delta(i-2)*(raynew(i) + rayold(i))
      end do
      
      end

c*********************************************************************************************
c Unterprogramm PARAXIAL_P / Hauptprogramm WFRONT
c Wird aufgerufen von INTERSECT_WFRONT, GET_CRITRAY
c Berechnung der Slowness im Punkt raynew(1),raynew(2)
c durch paraxiales raytracing.
c
c pnew(i)    = pold(i) + W(i,j)*delta(j)
c
c delta(i)   = xnew(i) - xold(i)
c w(i,j)     = omega(k,i) * a(k,m) * omega(m,j)
c omega(k,i) = siehe unten...Tangente
c a(k,m)     = siehe unten...
c
c tangent(i) = pold(i) / pold
c 
      subroutine paraxial_p(rayold,vel,
     &     raynew,N_PARAM)

      implicit none

      integer N_PARAM
      real raynew(4),rayold(N_PARAM)
      real delta(2),tangent(2)
      real omega(2,2),a(2,2),aux(2,2),w(2,2)
      real vn,vs,v2,p

      real vel(6)

      integer i,j,k

c p koennte man durch 1/vel(1) vertauschen ...Schnelligkeit
      p = sqrt(rayold(3)*rayold(3) + rayold(4)*rayold(4))
      do i = 1, 2
         delta(i) = raynew(i) - rayold(i)
         tangent(i) = rayold(i+2)/p
      end do
      omega(1,1) = tangent(2)
      omega(1,2) = -tangent(1)
      omega(2,1) = tangent(1)
      omega(2,2) = tangent(2)
      vn = vel(2)*tangent(2) - vel(3)*tangent(1)
      vs = vel(2)*tangent(1) + vel(3)*tangent(2)
      v2 = vel(1)*vel(1)

      a(1,1) = rayold(6)/rayold(5)
      a(1,2) = -vn/v2
      a(2,1) = -vn/v2
      a(2,2) = -vs/v2

      do 1 i = 1, 2
         do 1 j = 1, 2
            aux(i,j) = 0.0
            do 1 k = 1, 2
               aux(i,j) = aux(i,j) + a(i,k)*omega(k,j)
 1    continue      
      do 2 i = 1, 2
         do 2 j = 1, 2
            w(i,j) = 0.0
            do 2 k = 1, 2
               w(i,j) = w(i,j) + omega(k,i)*aux(k,j)
 2    continue

c-------------------
c slowness:
c
      do 4 i = 3, 4
         raynew(i) = rayold(i)
         do 5 j = 1, 2
            raynew(i) = raynew(i) + w(i-2,j)*delta(j)
 5       continue
 4    continue


      end


c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm RUNKUTTA / Hauptprogramm WFRONT
c Aufruf des Unterprogramms DERY
c Raytracing. Propagiert einen Strahl um den Zeitschritt DTRAY anhand
c eines Runge-Kutta-Verfahrens 4. Ordnung
c Ausgabe:
c RAYNEW : Neue Strahlenparameter
c
      subroutine runkutta(ray,dtray,dthalf,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     raynew,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer N_PARAM, MAXP_XGRID, MAXP_ZGRID
      integer veltype,actlayer
      real ray(N_PARAM), raynew(N_PARAM), ptos, dtray, dthalf
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      integer N_PARAM_KH
      parameter( N_PARAM_KH=6 )
      real kh1(N_PARAM_KH), kh2(N_PARAM_KH), kh3(N_PARAM_KH), kh4(N_PARAM_KH)
      real raytmp(N_PARAM)
      integer i
      

      call dery(ray,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &     kh1,
     &     N_PARAM_KH,MAXP_XGRID,MAXP_ZGRID)

      do 1 i = 1, N_PARAM_KH
         raytmp(i) = ray(i) + dthalf*kh1(i)
 1    continue

      call dery(raytmp,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &     kh2,
     &     N_PARAM_KH,MAXP_XGRID,MAXP_ZGRID)

      do 2 i = 1, N_PARAM_KH
         raytmp(i) = ray(i) + dthalf*kh2(i)
 2    continue

      call dery(raytmp,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &     kh3,
     &     N_PARAM_KH,MAXP_XGRID,MAXP_ZGRID)

      do 3 i = 1, N_PARAM_KH
         raytmp(i) = ray(i) + dtray*kh3(i)
 3    continue

      call dery(raytmp,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &     kh4,
     &     N_PARAM_KH,MAXP_XGRID,MAXP_ZGRID)

      do 4 i = 1, N_PARAM_KH
         raynew(i) = ray(i) + dtray/6.0*(kh1(i) + kh4(i) + 2.0*(kh2(i) + kh3(i)))
 4    continue

      do i = 7, 8
         raynew(i) = ray(i)
      end do

      end

c*****************************************************************************
c Unterprogramm DERY_EXTRA
c Aufruf des Unterprogramms VELOCITY
c Bestimmung der Hilfswerte K_OVER_H fuer das Runge-Kutta-Verfahren
c Ausserdem: Rueckgabe der Geschwindigkeiten VEL(6)
c
      subroutine dery_extra(ray,nx_grid,nz_grid,x_grid,z_grid,
     &     v_grid,ptos,veltype,actlayer,
     &     k_over_h,vel,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)
      
      implicit none

      integer MAXP_XGRID, MAXP_ZGRID, N_PARAM
      integer veltype, actlayer
      real ray(N_PARAM), ptos
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)
      real k_over_h(N_PARAM)
      real vel(6),vv

      call velocity(ray(1), ray(2), nx_grid, nz_grid,
     &     x_grid, z_grid, v_grid,
     &     veltype, ptos, actlayer,
     &     vel,
     &     MAXP_XGRID, MAXP_ZGRID)

      vv = vel(1) * vel(1)
      k_over_h(1) = vv * ray(3)
      k_over_h(2) = vv * ray(4)
      k_over_h(3) = -vel(2)/vel(1)
      k_over_h(4) = -vel(3)/vel(1)

      if (N_PARAM .eq. 4) goto 999

      k_over_h(5) = vv * ray(6)
      k_over_h(6) = -(vel(4)*ray(4)*ray(4) - 2.0*vel(5)*ray(3)*ray(4) + vel(6)*ray(3)*ray(3)) * vel(1)*ray(5)

c vel(6) ist eh 0! (vel(6) = d2v/dz2)
c      k_over_h(6) = -(vel(4)*ray(4)*ray(4) - 2.0*vel(5)*ray(3)*ray(4)) * vel(1)*ray(5)

  999 end

      
c*****************************************************************************
c Unterprogramm DERY
c Aufruf des Unterprogramms VELOCITY
c Bestimmung der Hilfswerte K_OVER_H fuer das Runge-Kutta-Verfahren
c
      subroutine dery(ray,nx_grid,nz_grid,x_grid,z_grid,
     &     v_grid,ptos,veltype,actlayer,
     &     k_over_h,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer MAXP_XGRID, MAXP_ZGRID, N_PARAM
      integer veltype,actlayer
      real ray(N_PARAM), ptos
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)
      real k_over_h(N_PARAM)
      real vel(6),vv

      call velocity(ray(1), ray(2), nx_grid, nz_grid,
     &     x_grid, z_grid, v_grid,
     &     veltype, ptos, actlayer,
     &     vel,
     &     MAXP_XGRID, MAXP_ZGRID)

      vv = vel(1) * vel(1)
      k_over_h(1) = vv * ray(3)
      k_over_h(2) = vv * ray(4)
      k_over_h(3) = -vel(2)/vel(1)
      k_over_h(4) = -vel(3)/vel(1)

      if (N_PARAM .eq. 4) goto 999

      k_over_h(5) = vv * ray(6)
      k_over_h(6) = -(vel(4)*ray(4)*ray(4) - 2.0*vel(5)*ray(3)*ray(4) + vel(6)*ray(3)*ray(3)) * vel(1)*ray(5)

c vel(6) ist eh 0! (vel(6) = d2v/dz2)
c      k_over_h(6) = -(vel(4)*ray(4)*ray(4) - 2.0*vel(5)*ray(3)*ray(4)) * vel(1)*ray(5)

 999  end
 

c***********************************************************************************************************
c Unterprogramm RUNKUTTA2 / Hauptprogramm WFRONT
c Aufruf des Unterprogramms DERY
c Raytracing. Propagiert einen Strahl um den Zeitschritt DTRAY anhand
c eines Runge-Kutta-Verfahrens 4. Ordnung
c Unterschied zu RUNKUTTA : Die Hilfsgroessen KH1() sind bereits bekannt!
c Ausgabe:
c RAYNEW : Neue Strahlenparameter
c
      subroutine runkutta2(ray,dtray,dthalf,kh1,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     raynew,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer N_PARAM, MAXP_XGRID, MAXP_ZGRID
      integer veltype,actlayer
      real ray(N_PARAM), raynew(N_PARAM), ptos, dtray, dthalf
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      integer N_PARAM_KH
      parameter( N_PARAM_KH=6 )
      real kh1(N_PARAM_KH), kh2(N_PARAM_KH), kh3(N_PARAM_KH), kh4(N_PARAM_KH)
      real raytmp(N_PARAM)
      integer i, nparam_min

c Manchmal wird diese Funktion mit N_PARAM=4 aufgerufen
c Stelle sicher dass keine Feldgrenzen ueberschritten werden:
      nparam_min = N_PARAM
      if ( N_PARAM_KH .lt. N_PARAM ) then
         nparam_min = N_PARAM_KH
      end if

      do 1 i = 1, nparam_min
         raytmp(i) = ray(i) + dthalf*kh1(i)
 1    continue

      call dery(raytmp,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &     kh2,
     &     nparam_min,MAXP_XGRID,MAXP_ZGRID)

      do 2 i = 1, nparam_min
         raytmp(i) = ray(i) + dthalf*kh2(i)
 2    continue

      call dery(raytmp,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &     kh3,
     &     nparam_min,MAXP_XGRID,MAXP_ZGRID)

      do 3 i = 1, nparam_min
         raytmp(i) = ray(i) + dtray*kh3(i)
 3    continue

      call dery(raytmp,nx_grid,nz_grid,x_grid,z_grid,v_grid,ptos,veltype,actlayer,
     &     kh4,
     &     nparam_min,MAXP_XGRID,MAXP_ZGRID)

      do 4 i = 1, nparam_min
         raynew(i) = ray(i) + dtray/6.0*(kh1(i) + kh4(i) + 2.0*(kh2(i) + kh3(i)))
 4    continue

      if ( N_PARAM .gt. 6 ) then
         do i = 7, 8
            raynew(i) = ray(i)
         end do
      end if

      end



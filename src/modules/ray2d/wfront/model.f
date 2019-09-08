c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm MODEL / Hauptprogramm WFRONT
c Die zugrunde liegenden Routinen entstammen demm Programmpaket SEIS88
c Wird aufgerufen von WFRONT
c Aufruf des Unterprogramms SPLIN, MODEL_COMMON, GET_MAX, GET_MIN, GET_Z, BIAP,
c VELONLY_LINEAR, VEL_INITIALIZE
c Erstellt das Modell, d.h. bestimmt die Schichtgrenzen anhand
c von kubischen Splines mit Interpolationspunkten x_int(), z_int().
c
      subroutine model(n_int,n_boreholes,n_allint,npoints_int,x_int,z_int,iii,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     index_grid,v_top,v_bottom,vpvs_ratio,f_out,flag_smooth,
     &     b_int,c_int,d_int,
     &     left_border,right_border,
     &     iii2,error,
     &     MAXPOINTS_INT, MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer n_allint, MAXPOINTS_INT, MAXP_XGRID,MAXP_ZGRID
      integer n_int, npoints_int(n_allint), n_boreholes
      integer iii(MAXPOINTS_INT,n_allint),iii2(MAXPOINTS_INT,n_int)
      real x_int(MAXPOINTS_INT,n_allint), z_int(MAXPOINTS_INT,n_allint)
      integer nx_grid(n_int),nz_grid(n_int)
      real x_grid(MAXP_XGRID,n_int),z_grid(MAXP_ZGRID,n_int),v_grid(4,MAXP_ZGRID,MAXP_XGRID,n_int)
      integer f_out,flag_smooth
      integer index_grid(n_int)
      real v_top(n_int),v_bottom(n_int),vpvs_ratio(n_int),left_border,right_border

      real b_int(MAXPOINTS_INT,n_allint), c_int(MAXPOINTS_INT,n_allint)
      real d_int(MAXPOINTS_INT,n_allint)

      real zmin, zmax
      integer nc, i, j, ilay, iiisave
      
      integer MAX_INT_LOCAL,error
      parameter (MAX_INT_LOCAL = 50)
      real dz,dx

      integer actp, begin, actint, actp2, actint2, count
      real ztmp,x_bore,delta,x0,z0,x1,z1

c      integer found, i2, n, nx_velout, nz_velout
c      real zout
c      integer MAXX_VELOUT,MAXZ_VELOUT, f_velout
c      parameter (MAXX_VELOUT = 500, MAXZ_VELOUT = 500)
c      real x_output(MAXX_VELOUT), z_output(MAXZ_VELOUT), v_output(MAXZ_VELOUT)
c      real zint(MAXX_VELOUT,MAX_INT_LOCAL)
c      character*30 fname_velout

      if (MAX_INT_LOCAL .lt. n_int) goto 980

c==========================================================================================
c Interfaces
c
c left and right border:
      left_border = x_int(1,1)
      right_border = x_int(npoints_int(1),1)
c-------------------------------
c Initialisierung:
      do actint = 1, n_int
         do actp = 1, npoints_int(actint)
            iii2(actp,actint) = 0
         end do
      end do
c------
c Pruefe erste und letzte Grenzflaeche:
      do actp = 1, npoints_int(1)
         if (iii(actp,1) .ne. 0 .and. iii(actp,1) .ne. -1) then
            actint = 1
            goto 920
         else if (actp .gt. 1 .and. actp .lt. npoints_int(1)) then
            if (iii(actp,1) .eq. -1 .and. index_grid(1) .eq. 1) then
               actint = 1
               goto 925
            end if
         end if
      end do
      do actp = 1, npoints_int(n_int)
         if (iii(actp,n_int) .ne. 0 .and. iii(actp,n_int) .ne. -1) then
            actint = n_int
            goto 920
         else if (actp .gt. 1 .and. actp .lt. npoints_int(n_int)) then
            if (iii(actp,n_int) .eq. -1 .and. index_grid(n_int-1) .eq. 1) then
               actint = n_int
               goto 925
            end if
         end if
      end do

c------
c Setze III2 und weitere Ueberpruefung:
c
      do actint = n_int-1, 2, -1
         if (iii(1,actint) .eq. 0) then
            iii(1,actint) = -1
            write(f_out,1340,err=990) actint,1
         end if
         if (iii(npoints_int(actint),actint) .ne. -1) then
            iii(npoints_int(actint),actint) = -1
            write(f_out,1340,err=990) actint,npoints_int(actint)
         end if
         do actp = 1, npoints_int(actint)-1
            if (x_int(actp+1,actint) .le. x_int(actp,actint)) then
               write(*,*) x_int(actp+1,actint),x_int(actp,actint),actp, actint
               goto 940
            end if
            if (iii2(actp,actint) .eq. 0) then
               if (iii(actp,actint) .gt. 0) then
                  if (iii(actp+1,actint) .eq. 0 .and. flag_smooth .eq. 0) then
                     iii(actp+1,actint) = -1
                     write(f_out,1340,err = 990) actint,(actp+1)
                  end if
                  actp2 = iii(actp,actint)
                  actint2 = actint-1
                  x0 = x_int(actp,actint)
                  z0 = z_int(actp,actint)
                  x1 = x_int(actp+1,actint)
                  z1 = z_int(actp+1,actint)
                  count = 2
                  do while (iii(actp2,actint2) .gt. 0)
                     if (x0 .ne. x_int(actp2,actint2) .or. z0 .ne. z_int(actp2,actint2)) goto 930
                     if (x1 .ne. x_int(actp2+1,actint2) .or. z1 .ne. z_int(actp2+1,actint2)) goto 935
                     iii2(actp2,actint2) = 1
                     actp2 = iii(actp2,actint2)
                     actint2 = actint2-1
                     count = count + 1
                  end do
                  if (x0 .ne. x_int(actp2,actint2) .or. z0 .ne. z_int(actp2,actint2)) goto 930
                  if (x1 .ne. x_int(actp2+1,actint2) .or. z1 .ne. z_int(actp2+1,actint2)) goto 935
                  if (iii(actp2,actint2) .eq. -2) then
                     iii2(actp,actint) = -count
                     iii2(actp2,actint2) = -count
                  else
                     iii2(actp,actint) = count
                     iii2(actp2,actint2) = count
                  end if
c Setzen von III aller coinciding interfaces auf ein- und denselben Wert:
                  if (flag_smooth .ne. 0) then
                     iiisave = iii(actp2+1,actint2)
                     actp2 = actp
                     actint2 = actint
                     if (iiisave .eq. -2) iiisave = -1
                     if (iiisave .gt. 0) iiisave = -1
                     do i = 2, count - 1
                        if (iii(actp2+1,actint2) .lt. 1) then
                           iii(actp2+1,actint2) = iiisave
                        end if
                        actp2 = iii(actp2,actint2)
                        actint2 = actint2 - 1
                     end do
                  end if
               end if
            end if
         end do
      end do

      call model_common(n_int,npoints_int,x_int,z_int,iii,
     &     b_int,c_int,d_int,error,
     &     MAXPOINTS_INT)
      if (error .ne. 0) goto 999
c==========================================================================================
c Boreholes
c

c check borehole boundary points
      begin = n_int + 1
      do i = begin, n_boreholes + n_int
         x_bore = z_int(1,i)
         call get_z(x_bore,x_int(1,1),z_int(1,1),d_int(1,1),c_int(1,1),b_int(1,1),npoints_int(1),
     &        ztmp,actp,
     &        MAXPOINTS_INT)
         if (ztmp .ne. x_int(1,i)) then
            x_int(1,i) = ztmp
            write(f_out,1100,err = 990) i,x_int(1,i)
         end if
         nc = npoints_int(i)
         x_bore = z_int(nc,i)
         call get_z(x_bore,x_int(1,n_int),z_int(1,n_int),d_int(1,n_int),c_int(1,n_int),b_int(1,n_int),npoints_int(n_int),
     &        ztmp,actp,
     &        MAXPOINTS_INT)
         if (ztmp .lt. x_int(nc,i)) then
            x_int(nc,i) = ztmp
            write(f_out,1110,err = 990) nc,i,x_int(nc,i)
         end if
      end do

c x- und z-Koordinate werden vertauscht!
      call model_common(n_boreholes,npoints_int(begin),x_int(1,begin),z_int(1,begin),iii(1,begin),
     &     b_int(1,begin),c_int(1,begin),d_int(1,begin),error,
     &     MAXPOINTS_INT)
      if (error .ne. 0) goto 999

c==========================================================================================
c Geschwindigkeitsmodell:
c

c Genauigkeitsparameter (in Meter)... Wird zum Ergebnis der Bestimmung der Maximal- und Minimalausdehnung
c von Schichten hinzuaddiert (-subtrahiert)
      delta = 0.01
      do ilay = 1, n_int-1
c Geschwindigkeitsgitter muss erst berechnet werden:
         if (index_grid(ilay) .ne. 0) then
            nx_grid(ilay) = 2
            nz_grid(ilay) = 2
            call get_min(x_int(1,ilay),z_int(1,ilay),d_int(1,ilay),c_int(1,ilay),b_int(1,ilay),npoints_int(ilay),
     &           zmin,
     &           MAXPOINTS_INT)
            call get_max(x_int(1,ilay+1),z_int(1,ilay+1),d_int(1,ilay+1),c_int(1,ilay+1),b_int(1,ilay+1),npoints_int(ilay+1),
     &           zmax,
     &           MAXPOINTS_INT)
c Zur Sicherheit, damit Schicht wirklich innerhalb des grids liegt
            zmin = zmin - delta
            zmax = zmax + delta
            dz = (zmax-zmin) / float(nz_grid(ilay)-1)
            dx = (right_border-left_border) / float(nx_grid(ilay)-1)
            do i = 1, nx_grid(ilay)
               x_grid(i,ilay) = dx*float(i-1) + left_border
            end do
            do i = 1, nz_grid(ilay)
               z_grid(i,ilay) = dz*float(i-1) + zmin
               do j = 1, nx_grid(ilay)
                  call velonly_linear(x_grid(j,ilay),z_grid(i,ilay),ilay,
     &                 npoints_int,x_int,z_int,b_int,c_int,d_int,v_bottom,v_top,1,vpvs_ratio,
     &                 v_grid(1,i,j,ilay),
     &                 n_allint,MAXPOINTS_INT)
               end do
            end do
         end if
         call biap(nx_grid(ilay),nz_grid(ilay),x_grid(1,ilay),z_grid(1,ilay),
     &        v_grid(1,1,1,ilay),error,
     &        MAXP_XGRID, MAXP_ZGRID)
         if (error .ne. 0) goto 999
      end do
c-----------------
c Initialisierung der routine VELOCITY:
c
      ilay = 1
      call vel_initialize( x_int(1,ilay), z_int(1,ilay), nx_grid(ilay), nz_grid(ilay),
     &     x_grid(1,ilay), z_grid(1,ilay), v_grid(1,1,1,ilay), ilay,
     &     MAXP_XGRID, MAXP_ZGRID)

c---------------------------------
c Ausgabe des Modells in Datei f_intout     
c ...geloescht...

c Ausgabe der Bohrloecher
c ...geloescht...

c-----------------------------------------------------
c Ausgabe des Geschwindigkeitsmodells in F_VELOUT
c

c$$$      f_velout = 28
c$$$      fname_velout = "vel.out"
c$$$      nx_velout = 100
c$$$      nz_velout = 50
c$$$      if (f_velout .ne. 0) then
c$$$         write(*,*) "Plotting model..."
c$$$         if (nx_velout .gt. MAXX_VELOUT) then
c$$$            write(f_out,*,err=990) "WARNING: Too many grid points for velocity model plot!"
c$$$            write(f_out,*,err=990) "         NX_VELOUT was set to ",MAXX_VELOUT
c$$$            write(f_out,*,err=990) "         Enlarge PARAMETER 'MAXX_VELOUT' in WMODEL.F"
c$$$            nx_velout = MAXX_VELOUT
c$$$         end if
c$$$         if (nz_velout .gt. MAXX_VELOUT) then
c$$$            write(f_out,*,err=990) "WARNING: Too many grid points for velocity model plot!"
c$$$            write(f_out,*,err=990) "         NZ_VELOUT was set to ",MAXZ_VELOUT
c$$$            write(f_out,*,err=990) "         Enlarge PARAMETER 'MAXZ_VELOUT' in WMODELF"
c$$$            nz_velout = MAXZ_VELOUT
c$$$         end if
c$$$
c$$$         if (nx_velout .eq. 0) nx_velout = MAXX_VELOUT/2
c$$$         if (nz_velout .eq. 0) nz_velout = MAXZ_VELOUT/2
c$$$         dx = (right_border - left_border) / float(nx_velout -1)
c$$$         do i = 1, nx_velout
c$$$            x_output(i) = x_int(1,1) + float(i-1) * dx
c$$$            do actint = 1, n_int
c$$$               call get_z(x_output(i),x_int(1,actint),z_int(1,actint),d_int(1,actint),c_int(1,actint),b_int(1,actint),
c$$$     &              npoints_int(actint),
c$$$     &              zint(i,actint),actp,
c$$$     &              MAXPOINTS_INT)
c$$$            end do
c$$$         end do
c$$$         call get_min(x_int(1,1),z_int(1,1),d_int(1,1),c_int(1,1),b_int(1,1),npoints_int(1),
c$$$     &        zmin,
c$$$     &        MAXPOINTS_INT)
c$$$         call get_max(x_int(1,n_int),z_int(1,n_int),d_int(1,n_int),c_int(1,n_int),b_int(1,n_int),npoints_int(n_int),
c$$$     &        zmax,
c$$$     &        MAXPOINTS_INT)
c$$$         dz = (zmax - zmin) / float(nz_velout - 1)
c$$$         do j = 1, nz_velout
c$$$            z_output(j) = zmin + float(j-1) * dz
c$$$         end do
c$$$
c$$$         open(f_velout,file=fname_velout,form='unformatted',err=910)
c$$$         ilay = 1
c$$$         do i = 1, nx_velout
c$$$            do j = 1, nz_velout
c$$$               found = 0
c$$$               do while (found .eq. 0)
c$$$                  if (zint(i,ilay) .gt. z_output(j)) then
c$$$                     if (ilay .gt. 1) then
c$$$                        ilay = ilay - 1
c$$$                     else
c$$$                        found = -1
c$$$                     end if
c$$$                  else if (zint(i,ilay+1) .lt. z_output(j)) then
c$$$                     if (ilay .lt. (n_int-1)) then
c$$$                        ilay = ilay + 1
c$$$                     else
c$$$                        found = -1
c$$$                     end if
c$$$                  else
c$$$                     found = 1
c$$$                  end if
c$$$               end do
c$$$               call velonly(x_output(i),z_output(j),nx_grid(ilay),nz_grid(ilay),x_grid(1,ilay),z_grid(1,ilay),
c$$$     &              v_grid(1,1,1,ilay),1,vpvs_ratio(ilay),ilay,
c$$$     &              v_output(j),
c$$$     &              MAXP_XGRID, MAXP_ZGRID)
c$$$            end do
c$$$            write(f_velout) (v_output(j),j = 1, nz_velout)
c$$$         end do
c$$$         close(f_velout)
c$$$      end if
c$$$ 910  write(*,*) "MODEL: Error opening file "
c$$$      goto 995

      goto 999

 1100 format(" Changed x_int(1,",i2,") to ",f9.5)
 1110 format(" Changed x_int(",i2,",",i2,") to",f9.5)
 1300 format(/" Interface ",i2," may not have fictious parts"/,
     &     "  or an interface index (III) greater than 0"/)
 1310 format(/" Interface no. ",i2," has edge point(s) (III != 0).",/,
     &     " Defining the layer above or below as a layer of isovelocity interfaces (INDEX_GRID=1)",/,
     &     " will lead to wrong velocity distribution (velocity derivations not smooth).",/, 
     &     " Input a velocity grid (INDEX_GRID = 0)"/)
 1320 format(/" Point ", i2," of interface ",i2,": (",f10.5,",",f10.5,") is coinciding with",/,
     &     " point ", i2," of interface ",i2,": (",f10.5,",",f10.5,").",/,
     &     " But x- or/and z-coordinates are not the same! Change x- or/and z-coordinate!"/)
 1330 format(/" X-coordinates of interface no.",i2," not in ascending order! (point ",i2,")"/)
 1340 format(" WARNING: Changed index III of interface ",i2," point ",i2," to -1")


 920  write(*,1300) actint
      goto 995
 925  write(*,1310) actint
      goto 995
 930  write(*,1320) actp,actint,x_int(actp,actint),z_int(actp,actint),
     &     actp2,actint2,x_int(actp2,actint2),z_int(actp2,actint2)
      goto 995
 935  write(*,1320) actp+1,actint,x_int(actp+1,actint),z_int(actp,actint),
     &     actp2+1,actint2,x_int(actp2+1,actint2),z_int(actp2+1,actint2)
      goto 995
 940  write(*,1330) actint,actp+1
      goto 995
 980  write(*,*) "MODEL: Error MAX_INT_LOCAL!"
      goto 995
 990  write(*,*) "MODEL: Error while writing in file no.",f_out
 995  error = 1

 999  end
      
c*********************************************************************
c Subroutine for both interfaces and boreholes
c
      subroutine model_common(n_int,npoints_int,x_int,z_int,iii,
     &     b_int,c_int,d_int,error,
     &     MAXPOINTS_INT)

      implicit none

      integer MAXPOINTS_INT, n_int, error
      integer iii(MAXPOINTS_INT,n_int), npoints_int(n_int)
      real x_int(MAXPOINTS_INT,n_int),z_int(MAXPOINTS_INT,n_int)
      real b_int(MAXPOINTS_INT,n_int),c_int(MAXPOINTS_int,n_int)
      real d_int(MAXPOINTS_INT,n_int)

      integer MAXPOINTS_LOCAL
      parameter(MAXPOINTS_LOCAL = 200)
      real fx(MAXPOINTS_LOCAL), dx
      integer i, j, point0, point1, point2

      do 1 i = 1, n_int
         if (npoints_int(i) .gt. MAXPOINTS_LOCAL) then
            write(*,*)
            write(*,*) "MODEL_COMMON: Too many interface points!"
            write(*,*) "Enlarge local PARAMETER 'MAXPOINTS_LOCAL' in MODEL_COMMON"
            write(*,*) "at least to ",npoints_int(i)
            error = 1
            goto 999
         end if
         do j = 1, npoints_int(i)
            fx(j) = z_int(j,i)
         end do
         point1 = 1
         point2 = 1
         do while (point1 .lt. npoints_int(i))
            do while (iii(point2,i) .eq. 0 .and. point2 .lt. npoints_int(i))
               point2 = point2 + 1
            end do
            if (point1 .ne. point2) then
               if ((point2-point1) .eq. 1) then
                  d_int(point1,i) = 0.0
                  c_int(point1,i) = 0.0
                  b_int(point1,i) = (z_int(point2,i)-z_int(point1,i)) / (x_int(point2,i)-x_int(point1,i))
               else
                  fx(point1) = z_int(point1,i)
                  call splin(x_int(1,i),fx,point1,point2,error,MAXPOINTS_INT)
                  if (error .ne. 0) goto 999
                  do j = point1, point2-1
                     dx = x_int(j+1,i) - x_int(j,i)
                     d_int(j,i) = (fx(j+1) - fx(j)) / (3.*dx)
                     c_int(j,i) = fx(j)
                     b_int(j,i) = (z_int(j+1,i) - z_int(j,i)) / dx - dx * (fx(j+1) + 2.*fx(j)) / 3.0
                  end do
               end if
               point1 = point2
            else if (iii(point2,i) .gt. 0) then
               point0 = iii(point2,i)
               b_int(point2,i) = b_int(point0,i-1)
               c_int(point2,i) = c_int(point0,i-1)
               d_int(point2,i) = d_int(point0,i-1)
               point1 = point1 + 1
               point2 = point1
            else
               point2 = point2 + 1
            end if
         end do
 1    continue

 999  end
      
c*********************************************************************
c Unterprogramm SPLIN
c ...stammt aus dem Programmpaket SEIS88
c
c     Cubic spline interpolation with zero curvatures at the end points
c
      subroutine splin(x,fx,nmin,nmax,error,MAXPOINTS_INT)
      
      implicit none
      
      integer MAXPOINTS_INT, nmin, nmax, error
      real x(MAXPOINTS_INT), fx(MAXPOINTS_INT)

      integer MAXPOINTS_LOCAL
      parameter (MAXPOINTS_LOCAL = 200)
      real a(MAXPOINTS_LOCAL), b(MAXPOINTS_LOCAL), h(MAXPOINTS_LOCAL), f(MAXPOINTS_LOCAL)
      real d1,d2,xpom
      integer nmin1,nmax1,i,j
      
      if (nmax .gt. MAXPOINTS_LOCAL) then
         write(*,*)
         write(*,*) "SPLIN: Too many interface points!"
         write(*,*) "Enlarge local PARAMETER 'MAXPOINTS_LOCAL' in SPLIN"
         write(*,*) "at least to ",nmax
         error = 1
         goto 999
      end if
    
      if ((nmax-nmin) .eq. 1) goto 4
      nmin1 = nmin + 1
      nmax1 = nmax - 1
      h(nmin1) = x(nmin1) - x(nmin)
      d2 = (fx(nmin1) - fx(nmin)) / h(nmin1)
      do 1 i = nmin1, nmax1
        h(i+1) = x(i+1) - x(i)
        d1 = d2
        d2 = (fx(i+1) - fx(i)) / h(i+1)
        b(i) = h(i) + h(i+1)
        fx(i) = 3. * (d2-d1) / b(i)
        a(i) = h(i) / b(i)
        b(i) = h(i+1) / b(i)
    1 continue
    4 fx(nmin) = 0.
      fx(nmax) = 0.
      if ((nmax-nmin) .eq. 1) goto 999
      h(nmin) = 0.
      f(nmin) = 0.
      do 2 i = nmin1, nmax1
        xpom = 2. + a(i) * h(i-1)
        h(i) = -b(i) / xpom
        f(i) = (fx(i) - a(i)*f(i-1)) / xpom
    2 continue
      do 3 i = nmin, nmax1
        j = nmax1 - (i-nmin)
        fx(j) = h(j)*fx(j+1) + f(j)
    3 continue
      
  999 end

c*****************************************************************
c Unterprogramm BIAP / Hauptprogramm WFRONT
c Bestimmt die Splineparameter V_GRID(2-4,,,) fuer die Geschwindigkeitsverteilung
c anhand von bikubischer Splineinterpolation zwischen Stuetzpunkten mit Geschwindigkeiten V_GRID(1,,,)
c Die Routine stammt aus SEIS88, mit veraenderten Bezeichnungen
c V_GRID(2,,,) : A20
c V_GRID(3,,,) : A02
c V_GRID(4,,,) : A22
c
      subroutine biap(nx_grid,nz_grid,x_grid,z_grid,
     &        v_grid,error,
     &        MAXP_XGRID, MAXP_ZGRID)

      implicit none

      integer MAXP_XGRID, MAXP_ZGRID
      integer nx_grid, nz_grid, error
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)

      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      integer MAXPOINTS_LOCAL
      parameter (MAXPOINTS_LOCAL = 200)
      real fx(MAXPOINTS_LOCAL)
      integer i,j

C     THIS ROUTINE PERFORMS THE DETERMINATION OF COEFFICIENTS
C     OF THE BICUBIC SPLINE INTERPOLATION
C
      if (nx_grid .gt. MAXPOINTS_LOCAL .or. nz_grid .gt. MAXPOINTS_LOCAL) then
         write(*,*)
         write(*,*) "BIAP: Too many interface points!"
         write(*,*) "Enlarge local PARAMETER 'MAXPOINTS_LOCAL' in BIAP"
         if (nx_grid .gt. nz_grid) then
            write(*,*) "at least to ",nx_grid
         else
            write(*,*) "at least to ",nz_grid
         end if
         error = 1
         goto 999
      end if
    
      do 3 i = 1, nz_grid
         do j = 1, nx_grid
            fx(j) = v_grid(1,i,j)
         end do
         call splin(x_grid,fx,1,nx_grid,error,MAXP_XGRID)
         if (error .ne. 0) goto 999
         do j = 1, nx_grid
            v_grid(2,i,j) = fx(j)
         end do
 3    continue

      do 6 j = 1, nx_grid
         do i = 1, nz_grid
            fx(i) = v_grid(1,i,j)
         end do
         call splin(z_grid,fx,1,nz_grid,error,MAXP_ZGRID)
         if (error .ne. 0) goto 999
         do i = 1, nz_grid
            v_grid(3,i,j) = fx(i)
         end do
 6    continue

      do 9 i = 1, nz_grid
         do j = 1, nx_grid
            fx(j) = v_grid(3,i,j)
         end do
         call splin(x_grid,fx,1,nx_grid,error,MAXP_XGRID)
         if (error .ne. 0) goto 999
         do j = 1, nx_grid
            v_grid(4,i,j) = fx(j)
         end do
 9    continue

 999  end


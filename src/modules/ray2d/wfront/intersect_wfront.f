c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm INTERSECT_WFRONT / Hauptprogramm WFRONT
c Wird aufgerufen von WFRONT,GET_BOUNDRAY
c Aufruf der Unterprogramme INTERW2, FIND_EDGE, PARAXIAL_P, VELOCITY, GET_Z, GET_Z1, GET_Z1_AT_INT
c Berechnet den Schnittpunkt der 'Wellenfront' mit einer Grenzflaeche.
c Genauer: Schnittpunkt einer Geraden durch RAY1 und RAY2 mit der Grenzflaeche.
c RAY1 liegt in der neuen Schicht (vielleicht! ist naemlich evt. bereits wieder herausgetreten)
c RAY2 liegt garantiert in der alten Schicht. (bei Aufruf von GET_BOUNDRAY auch anders moeglich)
c Ausgabe sind:
c x-,z,px-,pz-Komponente, dz/dx der Grenzflaeche am Schnittpunkt
c Vorgehensweise:
c Aufruf von INTERW2. Ausgabe sind der Schnittpunkt (RAYSECT(1),RAYSECT(2)) und dz/dx.
c Slowness (RAYSECT(3),RAYSECT(4)) wird anhand von PARAXIAL_P paraxial extrapoliert.
c
      subroutine intersect_wfront(ray1,ray2,interf,cross,
     &     npoints_int,x_int,z_int,b_int,c_int,d_int,iii,iii2,
     &     veltype, ptos,actlayer,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     dzdx,raysect,abbruch,
     &     MAX_INT,MAXPOINTS_INT,N_PARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer MAXPOINTS_INT,MAX_INT,N_PARAM,MAXP_XGRID,MAXP_ZGRID
      integer npoints_int(MAX_INT),iii(MAXPOINTS_INT,MAX_INT),iii2(MAXPOINTS_INT,MAX_INT)
      integer interf,veltype,actlayer,cross
      real ray1(N_PARAM),ray2(N_PARAM)
      real x_int(MAXPOINTS_INT,MAX_INT),z_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT),c_int(MAXPOINTS_INT,MAX_INT),b_int(MAXPOINTS_INT,MAX_INT)
      real ptos
      integer nx_grid,nz_grid,abbruch
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real dzdx,raysect(4)

      integer actp,i,actp1,actp2,flagnew, lastactp, iactp, firstactp, flag
      real wfront1(4), wfront2(4), raytmp(4)
      real z2tmp,dz2dxtmp,z1tmp,dz1dxtmp,vel(6),s1,s2
      real dzdxtmp,stmp,epsilon2,smin,dz1,dz2

      epsilon2 = 0.00005
      abbruch = 0
      actp2 = 1
      dzdx = 0.0

      call get_z1(ray2(1),x_int(1,interf),z_int(1,interf),d_int(1,interf),c_int(1,interf),b_int(1,interf),npoints_int(interf),
     &     z2tmp,dz2dxtmp,actp2,
     &     MAXPOINTS_INT)

      actp = actp2
      actp1 = actp2

      call get_z1(ray1(1),x_int(1,interf),z_int(1,interf),d_int(1,interf),c_int(1,interf),b_int(1,interf),npoints_int(interf),
     &     z1tmp,dz1dxtmp,actp1,
     &     MAXPOINTS_INT)

c Pruefe auf direkten Schnitt:
      dz1 = ray1(2) - z1tmp
      dz2 = ray2(2) - z2tmp
      if (abs(dz1) .lt. epsilon2) then
         raysect(1) = ray1(1)
         raysect(2) = ray1(2)
         goto 100
      else if (abs(dz2) .lt. epsilon2) then
         raysect(1) = ray2(1)
         raysect(2) = ray2(2)
         goto 100
      end if

c Pruefe, in welcher Schicht RAY1 und RAY2 liegen:
c FLAG =  0: in verschiedenen Schichten
c FLAG = -1: beide in alter Schicht
c sonst: beide in neuer Schicht; abbrechen, aber nicht mit ABBRUCH, sondern: Strahl schneidet zweimal die Grenzflaeche!
      if (cross .lt. 0) then
         if (dz1 .gt. 0.0) then
            flag = -1
         else if (dz2 .lt. 0.0) then
c Trick! Mit DZDX = 0.0 wird das Ende vom aufrufenden Programm GET_BOUNDRAY simuliert. Dieses bricht danach
c sowieso ab, weil der Strahl zweimal die Grenzschicht ueberschritten hat. :
            dzdx = 0.0
            raysect(1) = ray2(1)
            call get_z(raysect(1),x_int(1,interf),z_int(1,interf),d_int(1,interf),c_int(1,interf),b_int(1,interf),
     &           npoints_int(interf),
     &           raysect(2),actp,
     &           MAXPOINTS_INT)
            goto 100
         else
            flag = 0
         end if
      else
         if (dz1 .lt. 0.0) then
            flag = -1
         else if (dz2 .gt. 0.0) then
c Trick! Mit DZDX = 0.0 wird das Ende vom aufrufenden Programm GET_BOUNDRAY simuliert. Dieses bricht danach
c sowieso ab, weil der Strahl zweimal die Grenzschicht ueberschritten hat. :
            dzdx = 0.0
            raysect(1) = ray2(1)
            call get_z(raysect(1),x_int(1,interf),z_int(1,interf),d_int(1,interf),c_int(1,interf),b_int(1,interf),
     &           npoints_int(interf),
     &           raysect(2),actp,
     &           MAXPOINTS_INT)
            goto 100
         else
            flag = 0
         end if
      end if

      do i = 1, 2
         wfront1(i) = ray1(i)
         wfront2(i) = ray2(i)
      end do

      wfront1(3) = -ray1(4)
      wfront1(4) = ray1(3)
      wfront2(3) = -ray2(4)
      wfront2(4) = ray2(3)

c***********************************
      flagnew = 0
      if (actp1 .ne. actp2) then
         if (actp1 .gt. actp2) then
            iactp = actp2
            firstactp = actp2
            lastactp = actp1
c            istep = -1
         else
            iactp = actp1
            firstactp = actp1
            lastactp = actp2
c            istep = 1
         end if
         do while (iactp .ne. lastactp .and. flagnew .eq. 0)
            iactp = iactp + 1
            call find_edge(iactp,interf,cross,iii,iii2,flagnew,
     &           MAX_INT,MAXPOINTS_INT)
         end do
      end if
c-----------------------      
c flagnew .ne. 0: Mindestens ein Eckpunkt liegt zwischen RAY1 und RAY2
c
      if (flagnew .ne. 0) then
c         write(*,*) "INTERSECT_WFRONT: Eckpunkt dazwischen!"
         smin = 1.0E+22

         do iactp = firstactp, lastactp
            call get_z1_at_actp(ray1(1),x_int(iactp,interf),z_int(iactp,interf),d_int(iactp,interf),
     &           c_int(iactp,interf),b_int(iactp,interf),
     &           z1tmp,dz1dxtmp)
            call get_z1_at_actp(ray2(1),x_int(iactp,interf),z_int(iactp,interf),d_int(iactp,interf),
     &           c_int(iactp,interf),b_int(iactp,interf),
     &           z2tmp,dz2dxtmp)
c RAY2 MUSS ausserhalb der neuen Schicht sein:
            if ((ray2(2)-z2tmp)*float(cross) .le. 0.0) then
               call interw2(wfront1,wfront2,interf,z1tmp,dz1dxtmp,z2tmp,dz2dxtmp,
     &              x_int,z_int,d_int,c_int,b_int,npoints_int,
     &              dzdxtmp,raytmp,actp,
     &              MAX_INT,MAXPOINTS_INT)
               if (flag .eq. 0) then
                  if (iactp .eq. actp) then
                     stmp = (raytmp(1)-ray2(1))**2 + (raytmp(2)-ray2(2))**2
                     if (stmp .lt. smin) then
                        do i = 1, 4
                           raysect(i) = raytmp(i)
                        end do
                        smin = stmp
                        dzdx = dzdxtmp
                     end if
                  end if
               else if (dz2*(ray2(2)-raytmp(2)) .gt. 0.0) then
                  stmp = (raytmp(1)-ray2(1))**2 + (raytmp(2)-ray2(2))**2
                  if (stmp .lt. smin) then
                     do i = 1, 4
                        raysect(i) = raytmp(i)
                     end do
                     smin = stmp
                     dzdx = dzdxtmp
                  end if
               end if
            end if
         end do
c------
c Leider reichte diese Berechnung nicht zur Bestimmung des Schnittpunkts aus, daher eine etwas unsichere Art:
c
         if (smin .gt. 1.0E+20) then
            do iactp = firstactp, lastactp
               call get_z1_at_actp(ray1(1),x_int(iactp,interf),z_int(iactp,interf),d_int(iactp,interf),
     &              c_int(iactp,interf),b_int(iactp,interf),
     &              z1tmp,dz1dxtmp)
               call get_z1_at_actp(ray2(1),x_int(iactp,interf),z_int(iactp,interf),d_int(iactp,interf),
     &              c_int(iactp,interf),b_int(iactp,interf),
     &              z2tmp,dz2dxtmp)
               call interw2(wfront1,wfront2,interf,z1tmp,dz1dxtmp,z2tmp,dz2dxtmp,
     &              x_int,z_int,d_int,c_int,b_int,npoints_int,
     &              dzdxtmp,raytmp,actp,
     &              MAX_INT,MAXPOINTS_INT)
               if (flag .eq. 0) then
                  stmp = (raytmp(1)-ray2(1))**2 + (raytmp(2)-ray2(2))**2
                  if (stmp .lt. smin) then
                     do i = 1, 4
                        raysect(i) = raytmp(i)
                     end do
                     smin = stmp
                     dzdx = dzdxtmp
                  end if
               else if (dz2*(ray2(2)-raytmp(2)) .gt. 0.0) then
                  stmp = (raytmp(1)-ray2(1))**2 + (raytmp(2)-ray2(2))**2
                  if (stmp .lt. smin) then
                     do i = 1, 4
                        raysect(i) = raytmp(i)
                     end do
                     smin = stmp
                     dzdx = dzdxtmp
                  end if
               end if
            end do

            if (smin .gt. 1.0E+20) then
               write(*,*) "ERROR in INTERSECT_WFRONT"
            end if
         end if

c-------------
c Kein Eckpunkt dazwischen
c Normal
c
      else
         call interw2(wfront1,wfront2,interf,z1tmp,dz1dxtmp,z2tmp,dz2dxtmp,
     &        x_int,z_int,d_int,c_int,b_int,npoints_int,
     &        dzdx,raysect,actp,
     &        MAX_INT,MAXPOINTS_INT)
      end if
c NEW Ende
c****************************************************

 100  s1 = (raysect(1)-ray1(1))**2 + (raysect(2)-ray1(2))**2
      s2 = (raysect(1)-ray2(1))**2 + (raysect(2)-ray2(2))**2

c Aendern: dieselbe Genauigkeit, wie in get_rayonint:
      if (s1 .lt. s2) then
         if (abs(ray1(5)) .lt. 0.00001) then
            write(*,*) 'INTERSECT_W: CAUSTIC!'
            abbruch = 1
            goto 999
         end if
         call velocity(ray1(1), ray1(2), nx_grid, nz_grid,
     &        x_grid, z_grid, v_grid,
     &        veltype, ptos,actlayer,
     &        vel,
     &        MAXP_XGRID, MAXP_ZGRID)
         call paraxial_p(ray1,vel,
     &        raysect,N_PARAM)
      else
         if (abs(ray2(5)) .lt. 0.00001) then
            write(*,*) 'INTERSECT_W: CAUSTIC!'
            abbruch = 1
            goto 999
         end if
         call velocity(ray2(1), ray2(2), nx_grid, nz_grid,
     &        x_grid, z_grid, v_grid,
     &        veltype, ptos,actlayer,
     &        vel,
     &        MAXP_XGRID, MAXP_ZGRID)
         call paraxial_p(ray2,vel,
     &        raysect,N_PARAM)
      end if

 999  end

c************************************************************************************
c Unterprogramm INTERW2 / Hauptprogramm WFRONT
c Wird aufgerufen von INTERSECT_WFRONT
c Aufruf der Unterprogramme GET_ACTP, GET_Z1
c Berechnet den Schnittpunkt eines Strahls und einer Grenzflaeche
c Ausgabe ist:
c dz/dx auf Grenzflaeche am Schnittpunkt, Schnittpunkt RAYSECT
c
      subroutine interw2(ray1,ray2,interf,z1int,dz1dx,z2int,dz2dx,
     &     x_int,z_int,d_int,c_int,b_int,npoints_int,
     &     dzdx,raysect,actp,
     &     MAX_INT,MAXPOINTS_INT)

      implicit none

      real ray1(4),ray2(4),z2int,dz2dx,z1int,dz1dx
      integer MAXPOINTS_INT,MAX_INT,npoints_int(MAX_INT),interf
      real x_int(MAXPOINTS_INT,MAX_INT),z_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT),c_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT)
     
      real dzdx,raysect(2)
      integer actp

      double precision dx,dz,xleft,xright,xmid,zmid,dzleft,dzright,dzmid
      double precision bray,z0ray
      double precision c,b,dint,bint,cint,z0int,dboth,cboth,bboth,z0both
      double precision angleray, angleint, angletransf, cosa, sina
      double precision yint(2,2),yray(2,2),gradint(2),gradray(2)
      double precision yinttrans(2,2),yraytrans(2,2),gradinttrans(2)
      double precision epsilon1, epsilon2, pihalf

      integer flagray,count,MAX_COUNT,i
      parameter (MAX_COUNT = 20)

      pihalf = 1.5707963267
      epsilon1 = 0.001
      epsilon2 = 0.00005

      if (abs(ray2(2) - z2int) .lt. epsilon2) then
         raysect(1) = ray2(1)
c         raysect(2) = ray2(2)
         raysect(2) = z2int
         dzdx = dz2dx
         call get_actp(raysect(1),interf,x_int,npoints_int,
     &        actp,
     &        MAX_INT,MAXPOINTS_INT)
         goto 999
      else if (abs(ray1(2) - z1int) .lt. epsilon2) then
         raysect(1) = ray1(1)
c         raysect(2) = ray1(2)
         raysect(2) = z1int
         dzdx = dz1dx
         call get_actp(raysect(1),interf,x_int,npoints_int,
     &        actp,
     &        MAX_INT,MAXPOINTS_INT)
         goto 999
      end if

      yint(2,2) = z2int
      gradint(2) = dz2dx
      yint(2,1) = z1int
      gradint(1) = dz1dx

      do i = 1, 2
         yray(i,1) = ray1(i)
         yray(i,2) = ray2(i)
      end do

      if (abs(ray1(3)) .lt. epsilon1 .or. abs(ray2(3)) .lt. epsilon1) then
         flagray = 1
         gradray(1) = ray1(3)/ray1(4)
         gradray(2) = ray2(3)/ray2(4)
         angleray = dsign(1.D0,gradray(2))*pihalf - atan(gradray(2))
      else
         flagray = 0
         gradray(1) = ray1(4)/ray1(3)
         gradray(2) = ray2(4)/ray2(3)
         angleray = atan(gradray(2))         
      end if

      angleint = atan(gradint(2))
      angletransf = (angleray + angleint)*0.5
      if (angleint*angleray .lt. 0.0) then
         if (abs(angleint) + abs(angleray) .gt. pihalf) then
            angletransf = angletransf - dsign(1.D0,angletransf)*pihalf
         end if
      end if

      yint(1,2) = yray(1,2)
      dx = abs(yray(1,2) - yray(1,1))
      if (dx .lt. epsilon2) then
         yint(1,1) = yint(1,2) - epsilon2
         if (yint(1,1) .lt. x_int(1,interf)) then
            yint(1,1) = yint(1,2) + epsilon2
         end if
      else
         yint(1,1) = yray(1,1)
      end if

      cosa = cos(angletransf)
      sina = sin(angletransf)

      do i = 1, 2
         yraytrans(1,i) = yray(1,i)*cosa + yray(2,i)*sina
         yraytrans(2,i) = yray(2,i)*cosa - yray(1,i)*sina
         yinttrans(1,i) = yint(1,i)*cosa + yint(2,i)*sina
         yinttrans(2,i) = yint(2,i)*cosa - yint(1,i)*sina
         gradinttrans(i) = (gradint(i)*cosa - sina)/(cosa + gradint(i)*sina)
      end do

c------------------------------------------------
c Der Strahl in der Form zray(x) = bray * x + z0ray
c
      dx = yraytrans(1,2) - yraytrans(1,1)
      dz = yraytrans(2,2) - yraytrans(2,1)

      bray = dz/dx
      z0ray = yraytrans(2,1) - bray*yraytrans(1,1)

c------------------------------------------------
c Die Grenzflaeche in der Form zint(x) = dint * x**3 + cint * x**2 + bint * x + zint0
c
      dx = yinttrans(1,2) - yinttrans(1,1)
      dz = yinttrans(2,2) - yinttrans(2,1)

      b = gradinttrans(1)
      c = (3.*dz/dx - (gradinttrans(2) + 2.*gradinttrans(1)))/dx
      dint = ((gradinttrans(2) - gradinttrans(1))/dx - 2.*c)/(3.*dx)
      cint = c - 3.*dint*yinttrans(1,1)
      bint = b - (cint + c)*yinttrans(1,1)
      z0int = ((c - dint*yinttrans(1,1))*yinttrans(1,1) - b)*yinttrans(1,1) + yinttrans(2,1)
c----------------------------------------
      dboth = -dint
      cboth = -cint
      bboth = bray - bint
      z0both = z0ray - z0int
c----------------------------------------

c      dx1 = (yraytrans(1,2)-yraytrans(1,1))/20.
c      do i = 1, 20
c         xx = yraytrans(1,1) + dx1*float(i-1)
c         zz = bray*xx + z0ray
c         write(72,*) xx,zz
c      end do
c      write(72,*)

c      dx1 = (yinttrans(1,2)-yinttrans(1,1))/20.
c      do i = 1, 20
c         xx = yinttrans(1,1) + dx1*float(i-1)
c         zz = ((dint*xx + cint)*xx + bint)*xx + z0int
c         write(73,*) xx,zz
c      end do
c      write(73,*)

      xleft = yraytrans(1,1)
      dzleft = ((dboth*xleft + cboth)*xleft + bboth)*xleft + z0both
      if (abs(dzleft) .lt. epsilon2) then
         xmid = xleft
         goto 20
      end if
      xright = yraytrans(1,2)
      dzright = ((dboth*xright + cboth)*xright + bboth)*xright + z0both
      if (abs(dzright) .lt. epsilon2) then
         xmid = xright
         goto 20
      end if
c Untersuchungsgebiet wird vergroessert:
c
      if (dzright*dzleft .gt. 0.0) then
         dx = abs(xleft-xright)
         if (xright .lt. xleft) then
            do while (dzright*dzleft .gt. 0.0)
               xright = xright - dx
               dzright = ((dboth*xright + cboth)*xright + bboth)*xright + z0both
               if (dzright*dzleft .gt. 0.0) then
                  xleft = xleft + dx
                  dzleft = ((dboth*xleft + cboth)*xleft + bboth)*xleft + z0both
               end if
            end do
         else
            do while (dzright*dzleft .gt. 0.0)
               xright = xright + dx
               dzright = ((dboth*xright + cboth)*xright + bboth)*xright + z0both
               if (dzright*dzleft .gt. 0.0) then
                  xleft = xleft - dx
                  dzleft = ((dboth*xleft + cboth)*xleft + bboth)*xleft + z0both
               end if
            end do
         end if
      end if
            
      xmid = 0.5*(yraytrans(1,1) + yraytrans(1,2))
      dzmid = ((dboth*xmid + cboth)*xmid + bboth)*xmid + z0both

      count = 0
      do while (abs(dzmid) .gt. epsilon2)
         count = count + 1
         if (count .gt. MAX_COUNT) then
            write(*,*) "Interw2: max_count iterations!"
            goto 20
         end if
         if (dzmid*dzright .lt. 0.0) then
            xleft = xmid
            dzleft = dzmid
         else
            xright = xmid
            dzright = dzmid
         end if
         xmid = 0.5*(xleft + xright)
         dzmid = ((dboth*xmid + cboth)*xmid + bboth)*xmid + z0both
      end do
 20   continue
      zmid = bray*xmid + z0ray

c-----------------------
c Ruecktrafo:
c
      raysect(1) = xmid*cosa - zmid*sina

c      xmidinf = xmid + delta
c      zmidinf = bray*xmidinf + z0ray
c      zsect = zmid*cosa + xmid*sina
c      xsectinf = xmidinf*cosa - zmidinf*sina = (xmid+delta)*cosa - (bray*(xmid+delta) +z0ray)*sina
c      = xmid*cosa - zmid*sina + delta*cosa - bray*delta*sina
c      = xsect + delta*cosa - bray*delta*sina
c      zsectinf = zmidinf*cosa + xmidinf*sina
c      xreal = xmid*cosa - zmid*sina

      call get_z1(raysect(1),x_int(1,interf),z_int(1,interf),d_int(1,interf),c_int(1,interf),b_int(1,interf),
     &     npoints_int(interf),
     &     raysect(2),dzdx,actp,
     &     MAXPOINTS_INT)

 999  end



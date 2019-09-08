c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm INTERSECTION / Hauptprogramm WFRONT
c Wird aufgerufen von INTERSECT_RAY
c Aufruf des Unterprogramms GET_ACTP
c Berechnet den (ungefaehren) Schnittpunkt eines Strahls und einer Grenzflaeche.
c Ausgabe sind:
c ACTP, der Index des Grenzflaechenstuecks des Schnittpunktes, Weglaenge S vom Start des Strahls (RAY1)
c zum Schnittpunkt
c
      subroutine intersection(ray1,ray2,interface,z1int,dz1dx,z2int,dz2dx,
     &     x_int,npoints_int,
     &     s,actp,
     &     MAX_INT,MAXPOINTS_INT)

      implicit none

      real ray1(4),ray2(4),z2int,dz2dx,z1int,dz1dx,epsilon
      integer MAXPOINTS_INT,MAX_INT,npoints_int(MAX_INT),interface,actp      
      real x_int(MAXPOINTS_INT,MAX_INT)
     
      double precision xsect,zsect
      real s

      double precision dx,dz,xleft,xright,xmid,zmid,dzleft,dzright,dzmid
      double precision cray1,bray1,z0ray1,cray2,bray2,z0ray2,cray,bray,z0ray
      double precision c,b,dint,bint,cint,z0int,dboth,cboth,bboth,z0both
      double precision angleray, angleint, angletransf, cosa, sina
      double precision yint(2,2),yray(2,2),gradint(2),gradray(2)
      double precision yinttrans(2,2),yraytrans(2,2),gradinttrans(2),gradraytrans(2)
      double precision epsilon2, pihalf

      real xreal

      integer flagray,count,MAX_COUNT,i
      parameter (MAX_COUNT = 50)

      pihalf = 1.5707963267
      epsilon = 0.00005
      epsilon2 = 0.001

      if (abs(ray2(2) - z2int) .lt. epsilon) then
c         dzdx = dz2dx
c         write(*,*) "intersection: Warning!"
         dx = ray2(1) - ray1(1)
         dz = ray2(2) - ray1(2)
         s = sqrt(dx*dx + dz*dz)
         xreal = ray2(1)
         goto 999
      else if (abs(ray1(2) - z1int) .lt. epsilon) then
c         dzdx = dz1dx
         xreal = ray1(1)
         s = 0.0
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

      if (abs(ray1(3)) .lt. epsilon2 .or. abs(ray2(3)) .lt. epsilon2) then
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
      if (dx .lt. epsilon) then
         yint(1,1) = yint(1,2) - epsilon
         if (yint(1,1) .lt. x_int(1,interface)) then
            yint(1,1) = yint(1,2) + epsilon
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

      if (flagray .gt. 0) then
         do i = 1, 2
            gradraytrans(i) = (cosa - gradray(i)*sina)/(gradray(i)*cosa + sina)
         end do
      else
         do i = 1, 2
            gradraytrans(i) = (gradray(i)*cosa - sina)/(cosa + gradray(i)*sina)
         end do
      end if            

      angleray = atan(gradraytrans(1))
      angleint = atan(gradinttrans(2))

c------------------------------------------------
c Der Strahl in der Form zray(x) = cray * x**2 + bray * x + zray0
c
      dx = yraytrans(1,2) - yraytrans(1,1)
      dz = yraytrans(2,2) - yraytrans(2,1)

      cray1 = (dz/dx - gradraytrans(1))/dx
      bray1 = gradraytrans(1) - 2.*cray1*yraytrans(1,1)
      z0ray1 = yraytrans(2,1) + (cray1*yraytrans(1,1) - gradraytrans(1))*yraytrans(1,1)

c denke dran: hier negativ, weil -dx,-dz
      cray2 = (gradraytrans(2) - dz/dx)/dx
      bray2 = gradraytrans(2) - 2.*cray2*yraytrans(1,2)
      z0ray2 = yraytrans(2,2) + (cray2*yraytrans(1,2) - gradraytrans(2))*yraytrans(1,2)

      cray = 0.5*(cray1 + cray2)
      bray = 0.5*(bray1 + bray2)
      z0ray = 0.5*(z0ray1 + z0ray2)

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
c Differenz 'ray - interface': Abstand DZ bei gegebenem X
c
      dboth = -dint
      cboth = cray - cint
      bboth = bray - bint
      z0both = z0ray - z0int
c----------------------------------------

c      dx1 = (yraytrans(1,2)-yraytrans(1,1))/20.
c      if (abs(dx1) .gt. 0.1) then
c      do i = 1, 20
c         xx = yraytrans(1,1) + dx1*float(i-1)
c         zz = (cray*xx + bray)*xx + z0ray
c         write(72,*) xx,zz
c      end do       
c      write(72,*)
c      end if

c      dx1 = (yinttrans(1,2)-yinttrans(1,1))/20.
c      if (abs(dx1) .gt. 0.1) then
c         do i = 1, 20
c            xx = yinttrans(1,1) + dx1*float(i-1)
c            zz = ((dint*xx + cint)*xx + bint)*xx + z0int
c            write(73,*) xx,zz
c         end do
c         write(73,*)
c      end if

      xleft = yraytrans(1,1)
      dzleft = ((dboth*xleft + cboth)*xleft + bboth)*xleft + z0both
      if (abs(dzleft) .lt. epsilon) then
         xmid = xleft
         goto 20
      end if
      xright = yraytrans(1,2)
      dzright = ((dboth*xright + cboth)*xright + bboth)*xright + z0both
      if (abs(dzright) .lt. epsilon) then
         xmid = xright
         goto 20
      end if
c-----------------
c Fall: Durch die Transformation sind sowohl Anfangs- wie Endpunkt des Strahlstuecks
c beide innerhalb eine Schicht gerutscht: Das Untersuchungsfeld muss vergroessert werden
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
      do while (abs(dzmid) .gt. epsilon)
         count = count + 1
         if (count .gt. MAX_COUNT) then
            write(*,*) "Intersection: max_count iterations!"
            goto 20
         end if
         if (dzmid*dzleft .lt. 0.0) then
            xright = xmid
            dzright = dzmid
         else
            xleft = xmid
            dzleft = dzmid
         end if
         xmid = 0.5*(xleft + xright)
         dzmid = ((dboth*xmid + cboth)*xmid + bboth)*xmid + z0both
      end do
 20   continue
      zmid = (cray*xmid + bray)*xmid + z0ray

c-----------------------
c Ruecktrafo:
c
      xsect = xmid*cosa - zmid*sina
      zsect = zmid*cosa + xmid*sina

      dx = xsect - yray(1,1)
      dz = zsect - yray(2,1)

      s = sqrt(dx*dx + dz*dz)

      xreal = xsect

 999  call get_actp(xreal,interface,x_int,npoints_int,
     &     actp,
     &     MAX_INT,MAXPOINTS_INT)

      end



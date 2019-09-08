c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm INTERSECT_RAY / Hauptprogramm WFRONT
c Wird aufgerufen von GET_BOUNDRAY, PROPAGATE, DO_INTERSECTION
c Aufruf der Unterprogramme INTERSECTION, FIND_EDGE, RUNKUTTA, RUNKUTTA2, DERY, VELOCITY, GET_Z, GET_Z1, GET_Z1_AT_INT
c Berechnet den genauen Schnittpunkt des Strahls mit einer Grenzflaeche.
c Extra: dz/dx und d2z/dx2 am Schnittpunkt werden berechnet
c CROSS: <0: Strahl kommt von unten; >0: Strahl kommt von oben.
c FLAG: =-1: Strahl ist gerade von Grenzflaeche gestartet; =0: sonst.
c XIND: = 2: Schnitt mit Bohrloch; =1: sonst.
c FLAGRAY: 0: normal; !=0: spezielle Informationen sind bereits vorhanden (DTNEW)
c T1: aktueller Zeitpunkt
c Vorgehensweise:
c Aufruf von INTERSECTION. Rueckgabewert S ist die ungefaehre Weglaenge vom Start des Strahls (RAY1)
c zum Schnittpunkt. Die Zeit, die der Strahl fuer die Weglaenge S braucht, wird bestimmt;
c anhand von RUNKUTTA wird der Strahl RAY1 um diese Zeit weiterpropagiert.
c Ist die Genauigkeit noch nicht erreicht, beginnt ein iteratives Verfahren, dass
c mit RUNKUTTA den Schnittpunkt genauer bestimmt.
c Ausgabe sind:
c RAYNEW : Strahlenparameter am Schnittpunkt
c TNEW   : Zeit am Schnittpunkt
c DZDX,D2ZDX2 : Erste und zweite Ableitung der Grenzflaeche am Schnittpunkt
c
c RAYSAVE,DTSAVE : Wird INTERSECT_RAY nach GET_FLAG aufgerufen, ist FLAGRAY != 0.
c                  Dann wird der Endpunkt des Strahls (RAY2) neu gesetzt (notwendigerweise).
c                  Alte Werte werden in RAYSAVE und DTSAVE gespeichert und am Ende wieder zurueckgesetzt
c
      subroutine intersect_ray(ray1,ray2,actp2,cross,z2int,t1,flag,dtray,xind,dtnew,
     &     x_int,d_int,c_int,b_int,z_int,npoints_int,iii,iii2,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     actlayer,ptos,veltype,flagray,
     &     raynew,tnew,fictnew,dzdx,d2zdx2,
     &     N_PARAM,MAX_INT,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer MAX_INT,MAXPOINTS_INT,N_PARAM,MAXP_XGRID,MAXP_ZGRID
      integer actp2, veltype, actlayer, cross, flag, xind, flagray
      real ray1(N_PARAM),ray2(N_PARAM),z2int,t1,dtray
      integer npoints_int(MAX_INT), iii(MAXPOINTS_INT,MAX_INT),  iii2(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT),c_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT),z_int(MAXPOINTS_INT,MAX_INT)
      real x_int(MAXPOINTS_INT,MAX_INT),ptos
      integer nx_grid,nz_grid
      real x_grid(MAXP_XGRID),z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real raynew(N_PARAM), tnew, dzdx, d2zdx2
      integer fictnew

      real s,p1,dx,dtnew,dtnewhalf,dtleft,dtright,divide,z1int,dz1dx,kh1(N_PARAM)
      real dzleft,dzright,dzmid,zintnew,dz2dx
      real epsilon,distance,xtmp,smin,dz1dxtmp,dz2dxtmp,z1tmp,z2tmp,signum,chkbound,dividesave
      real factor,dzdxnew
      integer interface, zind, i, actp1, actp, flagnew, istep, iactp, firstactp, lastactp

      integer count, MAX_COUNTS
      parameter (MAX_COUNTS = 40)
      real dtsave, raysave(N_PARAM)

      epsilon = 0.000001
      interface = abs(cross)
      signum = float(cross/interface)

c-------------
c Bohrloch: xind = 2; sonst xind = 1
      if (xind .eq. 1) then
         zind = 2
      else
         zind = 1
      end if
c-------------
      dx = ray2(xind) - x_int(actp2,interface)
      dz2dx = (3.*d_int(actp2,interface)*dx + 2.*c_int(actp2,interface))*dx + b_int(actp2,interface)
      
      actp = actp2
      actp1 = actp2
      call get_z1(ray1(xind),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &     b_int(1,interface),npoints_int(interface),
     &     z1int,dz1dx,actp1,
     &     MAXPOINTS_INT)

c-------------
c Es sind bereits Informationen vorhanden...
c Aender: pruefen!
c RAY2 wird hier neu gesetzt, am Ende aber wieder auf Originalwert gesetzt
      if (flagray .ne. 0) then
         do i = 1, N_PARAM
            raysave(i) = ray2(i)
         end do
         dtsave = dtray
         if (dtnew .lt. epsilon) then
            dtnew = epsilon
         end if
         dtray = dtnew
         dtnewhalf = 0.5*dtnew
         call runkutta(ray1,dtnew,dtnewhalf,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        ray2,
     &        N_PARAM,MAXP_XGRID,MAXP_ZGRID)
      end if
c***********************************
c NEW
      flagnew = 0
      if (actp1 .ne. actp2) then
         if (actp1 .gt. actp2) then
            iactp = actp2
            firstactp = actp2
            lastactp = actp1
         else
            iactp = actp1
            firstactp = actp1
            lastactp = actp2
         end if
         do while (iactp .ne. lastactp .and. flagnew .eq. 0)
            iactp = iactp + 1
            call find_edge(iactp,interface,cross,iii,iii2,flagnew,
     &           MAX_INT,MAXPOINTS_INT)
         end do
      end if
c-----------------------      
c flagnew .ne. 0: Ein Eckpunkt liegt zwischen RAY1 und RAY2
c
      if (flagnew .ne. 0 .and. xind .eq. 1) then
c     write(*,*) "INTERSECT_RAY: Eckpunkt dazwischen!"
         if (actp2 .lt. actp1) then
            istep = -1
         else
            istep = 1
         end if
         smin = 1.0E+22
         do iactp = firstactp, lastactp
            call get_z1_at_actp(ray1(xind),x_int(iactp,interface),z_int(iactp,interface),d_int(iactp,interface),
     &           c_int(iactp,interface),b_int(iactp,interface),
     &           z1tmp,dz1dxtmp)
            call get_z1_at_actp(ray2(xind),x_int(iactp,interface),z_int(iactp,interface),d_int(iactp,interface),
     &           c_int(iactp,interface),b_int(iactp,interface),
     &           z2tmp,dz2dxtmp)
            if ((ray1(zind)-z1tmp)*(ray2(zind)-z2tmp) .le. 0.0) then
               call intersection(ray1,ray2,interface,z1tmp,dz1dxtmp,z2tmp,dz2dxtmp,
     &              x_int,npoints_int,
     &              s,actp,
     &              MAX_INT,MAXPOINTS_INT)
               if (actp .eq. iactp .and. s .lt. smin) smin = s
            end if
         end do
         s = smin
c was ist, wenn s = 1.0E+22 ?
c-------------
c Kein Eckpunkt dazwischen
      else
c Schnittpunkt gefunden:
         if (abs(ray2(zind)-z2int) .lt. epsilon) then
            do i = 1, N_PARAM
               raynew(i) = ray2(i)
            end do
            zintnew = z2int
            dtnew = dtray
            count = 0
            goto 101
         else
c-------------
c Normal
c--------------
c Bohrloch
            if (xind .eq. 2) then
               xtmp = ray1(2)
               ray1(2) = ray1(1)
               ray1(1) = xtmp
               xtmp = ray1(4)
               ray1(4) = ray1(3)
               ray1(3) = xtmp
               xtmp = ray2(2)
               ray2(2) = ray2(1)
               ray2(1) = xtmp
               xtmp = ray2(4)
               ray2(4) = ray2(3)
               ray2(3) = xtmp
            end if
            call intersection(ray1,ray2,interface,z1int,dz1dx,z2int,dz2dx,
     &           x_int,npoints_int,
     &           s,actp,
     &           MAX_INT,MAXPOINTS_INT)
         end if
      end if
c NEW Ende
c****************************************************
c--------------
c Bohrloch
      if (xind .eq. 2) then
         xtmp = ray1(2)
         ray1(2) = ray1(1)
         ray1(1) = xtmp
         xtmp = ray1(4)
         ray1(4) = ray1(3)
         ray1(3) = xtmp
         xtmp = ray2(2)
         ray2(2) = ray2(1)
         ray2(1) = xtmp
         xtmp = ray2(4)
         ray2(4) = ray2(3)
         ray2(3) = xtmp
      end if
c--------------
      distance = sqrt((ray2(zind)-ray1(zind))**2 + (ray2(xind)-ray1(xind))**2)
      dtnew = dtray*s/distance

c slowness am Start des Strahls (RAY1)
      p1 = sqrt(ray1(3)*ray1(3) + ray1(4)*ray1(4))

c=============================================================================================
c Bestimmung von dtleft, dzleft, dtright, dzright:
c
c--------------------------------------------------------------------------------
      call dery(ray1,nx_grid,nz_grid,x_grid,z_grid,
     &     v_grid,ptos,veltype,actlayer,
     &     kh1,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)
c-------------------------
c Wenn Anfangsstrahl (von Grenzflaeche loslaufender Strahl) sofort eine Grenzflaeche schneidet...
c
      if (dtnew .lt. epsilon .and. flag .eq. -1) then
         chkbound = signum
         dividesave = 0.99
         factor = 0.0

         do while (chkbound*signum .gt. 0.0 .and. (factor+0.1) .lt. 0.99)
            divide = dividesave
            factor = factor + 0.1
            dtnew = dtray * divide
            dtnewhalf = 0.5*dtnew
            call runkutta2(ray1,dtnew,dtnewhalf,kh1,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           ptos,veltype,actlayer,
     &           raynew,
     &           4,MAXP_XGRID,MAXP_ZGRID)
            call get_z1(raynew(xind),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &           b_int(1,interface),npoints_int(interface),
     &           zintnew,dzdxnew,actp,
     &           MAXPOINTS_INT)
            dzleft = raynew(zind) - zintnew
            dzright = ray2(zind) - z2int
            do while (dzleft*dzright .gt. 0.0)
               dividesave = divide
               divide = divide*factor
               dtnew = dtray*divide
               dtnewhalf = 0.5*dtnew
               call runkutta2(ray1,dtnew,dtnewhalf,kh1,
     &              nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &              ptos,veltype,actlayer,
     &              raynew,
     &              4,MAXP_XGRID,MAXP_ZGRID)
               call get_z1(raynew(xind),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &              b_int(1,interface),npoints_int(interface),
     &              zintnew,dzdxnew,actp,
     &              MAXPOINTS_INT)
               dzleft = raynew(zind) - zintnew
            end do
            chkbound = dzdxnew*raynew(3) - raynew(4)
         end do

         if (abs(dzleft) .lt. epsilon .and. dzleft*signum .le. 0.0) goto 100

         dtleft = dtnew
         dtright = dtray
c----------
c
      else if (dtnew .gt. dtray) then
c         write(*,*) "interr:Warning!___________"
         dzright = ray2(zind) - z2int
         dzleft = ray1(zind) - z1int
         dtleft = 0.0
         dtright = dtray
c----------
c Normal
c
      else
         dtnewhalf = dtnew*0.5
         call runkutta2(ray1,dtnew,dtnewhalf,kh1,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        raynew,
     &        4,MAXP_XGRID,MAXP_ZGRID)
         call get_z(raynew(xind),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &        b_int(1,interface),npoints_int(interface),
     &        zintnew,actp,
     &        MAXPOINTS_INT)
      
         dzleft = raynew(zind) - zintnew

         if (abs(dzleft) .lt. epsilon .and. dzleft*signum .le. 0.0) goto 100

         dtleft = dtnew
         dtnew = dtnew - 2. * signum * dzleft * p1
         if (dtnew .gt. dtray) then
            dtnew = dtray
            dzright = ray2(zind) - z2int
         else if (dtnew .lt. 0.0) then
            dtnew = 0.0
            dzright = ray1(zind) - z1int
         else
            dtnewhalf = 0.5 * dtnew      
            call runkutta2(ray1,dtnew,dtnewhalf,kh1,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           ptos,veltype,actlayer,
     &           raynew,
     &           4,MAXP_XGRID,MAXP_ZGRID)
            call get_z(raynew(xind),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &           b_int(1,interface),npoints_int(interface),
     &           zintnew,actp,
     &           MAXPOINTS_INT)
         end if

         dzright = raynew(zind) - zintnew
         if (abs(dzright) .lt. epsilon .and. signum*dzright .le. 0.0) goto 100
         count = 0

         do while(dzleft*dzright .gt. 0.0 .and. count .lt. MAX_COUNTS)
            count = count + 1
            dtnew = dtnew - signum * dzleft * p1
            if (dtnew .gt. dtray) then
               dtnew = dtray
               dzright = ray2(zind) - z2int
            else if (dtnew .lt. 0.0) then
               dtnew = 0.0
               dzright = ray1(zind) - z1int
            else
               dtnewhalf = 0.5 * dtnew
               call runkutta2(ray1,dtnew,dtnewhalf,kh1,
     &              nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &              ptos,veltype,actlayer,
     &              raynew,
     &              4,MAXP_XGRID,MAXP_ZGRID)
               call get_z(raynew(xind),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &              b_int(1,interface),npoints_int(interface),
     &              zintnew,actp,
     &              MAXPOINTS_INT)
               dzright = raynew(zind) - zintnew
            end if
         end do
         if (count .eq. MAX_COUNTS) then
            dtleft = 0.0
            dtright = dtray
            dzleft = ray1(zind) - z1int
            dzright = ray2(zind) - z2int
         end if
         if (abs(dzright) .lt. epsilon .and. signum*dzright .le. 0.0) goto 100
         dtright = dtnew
      end if

c==========================================================================================
c Iteration: finde genauen Schnittpunkt, fuer den gilt: | RAYNEW(2) - ZINTNEW | < ESPILON
c (RAYNEW(2): z-Koord. des Schnittpunkts, ZINTNEW: z-Koord. des eigentlichen Schnittpunktes mit
c  der Grenzflaeche bei RAYNEW(1))
c

      dtnew = 0.5*(dtleft + dtright)
      dtnewhalf = 0.5*dtnew
      
      call runkutta2(ray1,dtnew,dtnewhalf,kh1,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     raynew,
     &     4,MAXP_XGRID,MAXP_ZGRID)
      call get_z(raynew(xind),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &     b_int(1,interface),npoints_int(interface),
     &     zintnew,actp,
     &     MAXPOINTS_INT)
      
      dzmid = raynew(zind) - zintnew

      count = 0
c Abbruch, falls  | DZMID |  kleiner EPSILON ist und RAYNEW innerhalb der alten Schicht liegt (DZMID*SIGNUM > 0.0)
c Extra-Abbruch, falls MAX_COUNTS ueberschritten wird
      do while((abs(dzmid) .ge. epsilon .or. dzmid*signum .gt. 0.0) .and. count .lt. MAX_COUNTS)
         count = count + 1
         if (dzmid*dzright .gt. 0.0) then
            dzright = dzmid
            dtright = dtnew
         else
            dzleft = dzmid
            dtleft = dtnew
         end if
         dtnew = 0.5*(dtleft + dtright)
         dtnewhalf = 0.5*dtnew
         call runkutta2(ray1,dtnew,dtnewhalf,kh1,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &        ptos,veltype,actlayer,
     &        raynew,
     &        4,MAXP_XGRID,MAXP_ZGRID)
         call get_z(raynew(xind),x_int(1,interface),z_int(1,interface),d_int(1,interface),c_int(1,interface),
     &        b_int(1,interface),npoints_int(interface),
     &        zintnew,actp,
     &        MAXPOINTS_INT)
         dzmid = raynew(zind) - zintnew
      end do

 100  call runkutta2(ray1,dtnew,dtnewhalf,kh1,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     raynew,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

 101  tnew = t1 + dtnew
      raynew(zind) = zintnew
      dx = raynew(xind) - x_int(actp,interface)
      dzdx = (3.*d_int(actp,interface)*dx + 2.*c_int(actp,interface))*dx + b_int(actp,interface)
      d2zdx2 = 6.*d_int(actp,interface)*dx + 2.*c_int(actp,interface)
c-------------------------------------------------
c Setzen von FICTNEW (nicht bei Bohrloch (XIND = 2))

      if (xind .eq. 1) then
         fictnew = actp
         if (actp .eq. 0) then
            write(*,*) 'intersect_ray: actp=0!'
         end if
      end if

c Zuruecksetzen auf Originalwerte:
      if (flagray .ne. 0) then
         do i = 1, N_PARAM
            ray2(i) = raysave(i)
         end do
         dtray = dtsave
      end if

      end
      


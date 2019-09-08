c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

      SUBROUTINE COEF8(P,VP1,VS1,RO1,VP2,VS2,RO2,NCODE,ND,RMOD,RPH)
C
C     THE ROUTINE COEF8 IS DESIGNED FOR THE COMPUTATION OF REFLECTION
C     AND TRANSMISSION COEFFICIENTS AT A PLANE INTERFACE BETWEEN TWO
C     HOMOGENEOUS SOLID HALFSPACES OR AT A FREE SURFACE OF A HOMOGENEOUS
C     SOLID HALFSPACE.
C
      implicit none

      real p,vp1,vs1,ro1,vp2,vs2,ro2,rmod,rph
      integer ncode,nd

      real a1,a2,a3,a4,a5,a6,a7,a8,a9
      real e1,e2,e3,e4
      real g1,g2,g3,g4,g5,g6
      real s1,s2,s3,s4,s5,s6
      real z2
      real r,pp,x,qp,s,sb,q,y,z,z3,sc

      COMPLEX B(4),RR,C1,C2,C3,C4,H1,H2,H3,H4,H5,H6,H,HB,HC
      real PRMT(4),D(4),DD(4)
      integer i

      IF(NCODE.GE.9)GO TO 300
c---------------
c Trick fuer Reflexion an Erdoberflaeche:
c
      IF(RO2.LT.0.000001)GO TO 150

      PRMT(1)=VP1
      PRMT(2)=VS1
      PRMT(3)=VP2
      PRMT(4)=VS2
      A1=VP1*VS1
      A2=VP2*VS2
      A3=VP1*RO1
      A4=VP2*RO2
      A5=VS1*RO1
      A6=VS2*RO2
      Q=2.*(A6*VS2-A5*VS1)
      PP=P*P
      QP=Q*PP
      X=RO2-QP
      Y=RO1+QP
      Z=RO2-RO1-QP
      G1=A1*A2*PP*Z*Z
      G2=A2*X*X
      G3=A1*Y*Y
      G4=A4*A5
      G5=A3*A6
      G6=Q*Q*PP
      DO 21 I=1,4
      DD(I)=P*PRMT(I)
   21 D(I)=SQRT(ABS(1.-DD(I)*DD(I)))
      IF(DD(1).LE.1..AND.DD(2).LE.1..AND.DD(3).LE.1..AND.DD(4).LE.1.)
     1GO TO 100
C
C     COMPLEX COEFFICIENTS
      DO 22 I=1,4
      IF(DD(I).GT.1.)GO TO 23
      B(I)=CMPLX(D(I),0.)
      GO TO 22
   23 B(I)= CMPLX(0.,D(I))
   22 CONTINUE
      C1=B(1)*B(2)
      C2=B(3)*B(4)
      C3=B(1)*B(4)
      C4=B(2)*B(3)
      H1=G1
      H2=G2*C1
      H3=G3*C2
      H4=G4*C3
      H5=G5*C4
      H6=G6*C1*C2
      H=1./(H1+H2+H3+H4+H5+H6)
      HB=2.*H
      HC=HB*P
      GO TO (1,2,3,4,5,6,7,8),NCODE
    1 RR=H*(H2+H4+H6-H1-H3-H5)
      GO TO 26
    2 RR=VP1*B(1)*HC*(Q*Y*C2+A2*X*Z)
      IF(ND.NE.0)RR=-RR
      GO TO 26
    3 RR=A3*B(1)*HB*(VS2*B(2)*X+VS1*B(4)*Y)
      GO TO 26
    4 RR=-A3*B(1)*HC*(Q*C4-VS1*VP2*Z)
      IF(ND.NE.0)RR=-RR
      GO TO 26
    5 RR=-VS1*B(2)*HC*(Q*Y*C2+A2*X*Z)
      IF(ND.NE.0)RR=-RR
      GO TO 26
    6 RR=H*(H2+H5+H6-H1-H3-H4)
      GO TO 26
    7 RR=A5*B(2)*HC*(Q*C3-VP1*VS2*Z)
      IF(ND.NE.0)RR=-RR
      GO TO 26
    8 RR=A5*B(2)*HB*(VP1*B(3)*Y+VP2*B(1)*X)
      GO TO 26
C     REAL COEFFICIENTS
  100 E1=D(1)*D(2)
      E2=D(3)*D(4)
      E3=D(1)*D(4)
      E4=D(2)*D(3)
      S1=G1
      S2=G2*E1
      S3=G3*E2
      S4=G4*E3
      S5=G5*E4
      S6=G6*E1*E2
      S=1./(S1+S2+S3+S4+S5+S6)
      SB=2.*S
      SC=SB*P
      GO TO (101,102,103,104,105,106,107,108),NCODE
  101 R=S*(S2+S4+S6-S1-S3-S5)
      GO TO 250
  102 R=VP1*D(1)*SC*(Q*Y*E2+A2*X*Z)
      IF(ND.NE.0)R=-R
      GO TO 250
  103 R=A3*D(1)*SB*(VS2*D(2)*X+VS1*D(4)*Y)
      GO TO 250
  104 R=-A3*D(1)*SC*(Q*E4-VS1*VP2*Z)
      IF(ND.NE.0)R=-R
      GO TO 250
  105 R=-VS1*D(2)*SC*(Q*Y*E2+A2*X*Z)
      IF(ND.NE.0)R=-R
      GO TO 250
  106 R=S*(S2+S5+S6-S1-S3-S4)
      GO TO 250
  107 R=A5*D(2)*SC*(Q*E3-VP1*VS2*Z)
      IF(ND.NE.0)R=-R
      GO TO 250
  108 R=A5*D(2)*SB*(VP1*D(3)*Y+VP2*D(1)*X)
      GO TO 250
C
C     EARTHS SURFACE,COMPLEX COEFFICIENTS AND COEFFICIENTS OF CONVERSION
  150 A1=VS1*P
      A2=A1*A1
      A3=2.*A2
      A4=2.*A1
      A5=A4+A4
      A6=1.-A3
      A7=2.*A6
      A8=2.*A3*VS1/VP1
      A9=A6*A6
      DD(1)=P*VP1
      DD(2)=P*VS1
      DO 151 I=1,2
  151 D(I)=SQRT(ABS(1.-DD(I)*DD(I)))

      IF(DD(1).LE.1..AND.DD(2).LE.1.)GO TO 200

      DO 154 I=1,2
         IF(DD(I).GT.1.)GO TO 155
         B(I)=CMPLX(D(I),0.)
         GO TO 154
 155     B(I)= CMPLX(0.,D(I))
 154  CONTINUE

      H1=B(1)*B(2)
      H2=H1*A8
      H=1./(A9+H2)
      GO TO (161,162,163,164,165,166,167,168),NCODE

  161 RR=(-A9+H2)*H
      GO TO 26

  162 RR=-A5*B(1)*H*A6
      IF(ND.NE.0)RR=-RR
      GO TO 26

  163 RR=A5*B(2)*H*A6*VS1/VP1
      IF(ND.NE.0)RR=-RR
      GO TO 26

  164 RR=-(A9-H2)*H
      GO TO 26

  165 RR=A5*H1*H
      IF(ND.NE.0)RR=-RR
      GO TO 26

  166 RR=-A7*B(1)*H
      GO TO 26

  167 RR=A7*B(2)*H
      GO TO 26

  168 RR=A5*VS1*H1*H/VP1
      IF(ND.NE.0)RR=-RR

   26 Z2=REAL(RR)
      Z3=AIMAG(RR)
      IF(Z2.EQ.0..AND.Z3.EQ.0.)GO TO 157
      RMOD=SQRT(Z2*Z2+Z3*Z3)
      RPH=ATAN2(Z3,Z2)
      RETURN

  157 RMOD=0.
      RPH=0.
      RETURN
c================================================================================
c     EARTHS SURFACE,REAL COEFFICIENTS AND COEFFICIENTS OF CONVERSION
  200 S1=D(1)*D(2)
      S2=A8*S1
      S=1./(A9+S2)
      GO TO (201,202,203,204,205,206,207,208),NCODE

  201 R=(-A9+S2)*S
      GO TO 250

  202 R=-A5*D(1)*S*A6
      IF(ND.NE.0)R=-R
      GO TO 250

  203 R=A5*D(2)*S*A6*VS1/VP1
      IF(ND.NE.0)R=-R
      GO TO 250

  204 R=(S2-A9)*S
      GO TO 250

  205 R=A5*S1*S
      IF(ND.NE.0)R=-R
      GO TO 250

  206 R=-A7*D(1)*S
      GO TO 250

  207 R=A7*D(2)*S
      GO TO 250

  208 R=A5*VS1*S1*S/VP1
      IF(ND.NE.0)R=-R
  250 IF(R.LT.0.)GO TO 251
      RMOD=R
      RPH=0.
      RETURN

  251 RMOD=-R
      RPH=-3.14159
      RETURN
C
C     SH COEFFICIENTS OF REFLECTION, TRANSMISSION AND CONVERSION
C
  300 if(ro2.lt..000001)go to 302
      A1=P*VS1
      A2=P*VS2
      A3=RO1*VS1
      A4=RO2*VS2
      A5=SQRT(ABS(1.-A1*A1))
      A6=SQRT(ABS(1.-A2*A2))
      C1=A5
      IF(A2.LE.1.)C2=A6
      IF(A2.GT.1.)C2=CMPLX(0.,A6)
      C1=C1*A3
      C2=C2*A4
      H=1./(C1+C2)
      IF(NCODE.EQ.10)GO TO 301
      RR=H*(C1-C2)
      GO TO 26
  301 RR=2.*C1*H
      GO TO 26
  302 if(ncode.eq.10)go to 303
      RMOD=1.
      RPH=0.
      RETURN
  303 RMOD=2.
      RPH=0.
      RETURN
      END


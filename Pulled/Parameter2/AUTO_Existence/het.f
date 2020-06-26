C----------------------------------------------------------------------
C---------------------------------------------------------------------- 
C  het: Finding heteroclinic connections for trigger fronts in CGL
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
C 
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
C     ---------- ---- 
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION U(NDIM),PAR(*),F(NDIM),DFDU(NDIM,NDIM),DFDP(NDIM,*)
      DOUBLE PRECISION PI,LE,W,G,C,K,RM,RP,MU1R,MU1I,MU2,UW,UC, A
      DOUBLE PRECISION UG,MM,MP,UUWM,UUCM,UUG,UUWP,UUCP, AMP, KY, KY2
C 
       PI=2.*acos(0.0d0)
C      COEFFICIENT OF LINEAR TERM ON RIGHT SIDE< MAKES numerics closer to absolute spectrum
       AMP = 1.0d0
C     Variables are in unscaled numerical form (except 5 - 8)
       W=PAR(1)
       C= PAR(2)
       G = PAR(3)
       A = PAR(4)
       MU1R = PAR(5)
       MU1I = PAR(6)
       MU2 = PAR(7)
       LE = PAR(8)
       KY = PAR(9)
       KY2 = KY**2
CCCCCCCCCCC
C     First set of scalings (tilde to hat)
C       UW = PAR(1)+PAR(4) - 0.5d0*PAR(4)*PAR(2)**2/(1.d0+PAR(4)**2)
C       UC = PAR(2)/sqrt(1.d0+PAR(4)**2)
C       UG = PAR(3)
C

C
C
C      ODE: 1-3 left side, 4 - 6 right side, 1&4 real part, 2&4 Im part, 3&6 radius
       F(1)=-(U(1)**2-U(2)**2)-C*(U(1)+A*U(2))/(1+A**2)
       F(1)=F(1)+(A*(W+A*KY2)-1.d0+KY2+(1.d0+A*G)*U(3))/(1.d0+A**2)
       F(2)=-2.d0*U(1)*U(2)
       F(2)=F(2)+(A+W-C*(U(2)-A*U(1))+(G-A)*U(3))/(1.d0+A**2)
       F(3)=2.d0*U(3)*U(1)
C
       F(4)=-(U(4)**2-U(5)**2)-C*(U(4)+A*U(5))/(1.d0+A**2)
       F(4)=F(4)+(A*(W+A*KY2)+1.d0*AMP+KY2)/(1.d0+A**2)
       F(4)=F(4)+((1.d0+A*G)*U(6))/(1.d0+A**2)
       F(5)=-2.d0*U(4)*U(5)
       F(5)=F(5)+(W-A*AMP-C*(U(5)-A*U(4)))/(1.d0+A**2)
       F(5)=F(5)+(G-A)*U(6)/(1.d0+A**2)
       F(6)=2.d0*U(6)*U(4)
C
       F(1)=F(1)*LE
       F(2)=F(2)*LE
       F(3)= F(3)*LE
       F(4)= F(4)*LE
       F(5)= F(5)*LE
       F(6)= F(6)*LE
C 
C      IF(IJAC.EQ.0)RETURN 
C 
C       DFDU(1,1)=0.0 
C       DFDU(1,2)=1 
C       DFDU(2,1)=-PAR(1)*E 
C       DFDU(2,2)=0.0 
C 
C      IF(IJAC.EQ.1)RETURN 
C 
C      *Parameter derivatives
C       DFDP(1,1)=0.0 
C       DFDP(2,1)=-E 
C 
      RETURN 
      END 
C 
      SUBROUTINE STPNT(NDIM,U,PAR,T) 
C     ---------- ----- 
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION U(NDIM),PAR(*) 
      DOUBLE PRECISION PI,LE,K,SC1,SC2,UW,UC,UG,SQQ,SQQ1,ZPR,ZPI,ZMR,ZMI
      DOUBLE PRECISION MM,MP,UUCM,UUWM,UUCP,UUWP,UUG,PS1,PS2,PS11,GMA
      COMPLEX ZP, PS, SQP
C 
       PI=2*ACOS(0.0d0)
C      Variables in unscaled form. Be sure to imput the numerical form of the variable.

       PAR(1)=0.0d0
       PAR(2)=   1.0d0
       PAR(3)= -0.9d0
       PAR(4)= -0.1d0
       PAR(9) = 0.0d0



C     Periodic Orbit in Write-up Variables
       GMA = PAR(3)-PAR(4)
       K =-PAR(2)
       K =K+sqrt(PAR(2)**2+4.d0*(PAR(3)+PAR(1)-GMA*PAR(9)**2)*GMA)
       K =K/(2.d0*GMA)
C 
C
       ZMR = 0.d0
       ZMI = K
C
C       PS1= PAR(2)**2*(1.d0-PAR(4)**2)
C       PS1=PS1+ 4.d0*(1.d0+PAR(4)**2)*(1.d0+PAR(4)*PAR(1))
C       PS1=-4.d0-PAR(2)**2+4.d0*PAR(1)*PAR(4)
C
C       PS2=-2.d0*PAR(4)*PAR(2)**2
C       PS2=PS2+4.d0*(1+PAR(4)**2)*(-PAR(4)+PAR(1))
C       PS2=4.d0*(-PAR(4)-PAR(1))
C
       PS11=-1.d0-PAR(9)**2+PAR(4)*(PAR(4)*PAR(9)**2+PAR(1))
       PS1=PAR(2)**2-4.d0*PS11
       PS2=4.d0*(PAR(4)*(+1.d0+2.d0*PAR(9)**2)+PAR(1))
       PS = COMPLEX(PS1,PS2)
       SQP = sqrt(PS)
       ZP = COMPLEX(-PAR(2),0.d0) + SQP
       IF (REAL(SQP) .GT. 0) ZP= ZP-2.d0*SQP
       ZP = ZP/(2.d0*COMPLEX(1.d0,PAR(4)))
C
       ZPR = REAL(ZP)
       ZPI = AIMAG(ZP)
C      ZPR = -(PAR(2)/2.d0 + SQQ)
C      ZPI = -SQQ1
       PAR(5)=ZMR - ZPR
       PAR(6)=ZMI - ZPI
       PAR(7)=(1 - K**2-PAR(9)**2)
       PAR(8) = 1.d0

C      Initial Conditions
       U(1)=ZMR
       U(2)=ZMI
       U(3)=(1-K**2-PAR(9)**2)
       U(4)=ZPR
       U(5)=ZPI
       U(6)=0.d0
 
      RETURN 
      END 
C 
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
C     ---------- ---- 
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION PAR(*),ICP(*),U0(NDIM),U1(NDIM),FB(NBC),DBC(NBC,*)
      DOUBLE PRECISION PI,LE,K,SC1,SC2,UW,UC,UG,SQQ1,PS1,PS2
      DOUBLE PRECISION SQQ,ZPI,ZPR,ZMI,ZMR,UUW,UUC,UUG,MM,MP
      DOUBLE PRECISION PS11,GMA
      COMPLEX ZP, PS, SQP
C
C       Boundary Value parameters
C      Periodic Orbit in Write-up Variables
C       K =-PAR(2)
C       K =K+sqrt(PAR(2)**2+4*(PAR(3)+PAR(1))*(PAR(3) - PAR(4)))
C     K =K/(2.d0*(PAR(3)-PAR(4)))
       GMA = PAR(3)-PAR(4)
       K =-PAR(2)
       K =K+sqrt(PAR(2)**2+4.d0*(PAR(3)+PAR(1)-GMA*PAR(9)**2)*GMA)
       K =K/(2.d0*GMA)
C
       ZMR = 0.d0
       ZMI = (K)
C Make sure in bar variables!!!!
C For plus equations use P, for minus use M!!!
C       PS1= PAR(2)**2*(1.d0-PAR(4)**2)
C       PS1=PS1+ 4.d0*(1.d0+PAR(4)**2)*(1.d0+PAR(4)*PAR(1))
C       PS1=-4.d0-PAR(2)**2+4.d0*PAR(1)*PAR(4)
C
C       PS2=-2.d0*PAR(4)*PAR(2)**2
C       PS2=PS2+4.d0*(1+PAR(4)**2)*(-PAR(4)+PAR(1))
C       PS2=-4.d0*(PAR(4)+PAR(1))
C
C       PS = CMPLX(PS1,PS2)
C     ZP = CMPLX(0.d0,PAR(2)) + sqrt(PS)
C       PS1=
C       PS2=
C       PS=CMPLX(PS1,PS2)
C       ZP=CMPLX(PAR(2))+sqrt(PS)
C     IF (PAR(1) .GE. -PAR(4)) ZP= ZP-2.d0*sqrt(PS)
C       IF (REAL(sqrt(PS))   .LE. 0.d0) ZP=ZP-2.d0*sqrt(PS)
C     ZP = ZP/(2.d0*CMPLX(PAR(4),-1.d0))
       PS11=1.d0+PAR(9)**2-PAR(4)*(PAR(4)*PAR(9)**2+PAR(1))
       PS1=PAR(2)**2+4.d0*PS11
       PS2=4.d0*(PAR(4)*(+1.d0+2.d0*PAR(9)**2)+PAR(1))
       PS = COMPLEX(PS1,PS2)
       SQP = sqrt(PS)
       ZP = COMPLEX(-PAR(2),0.d0) + SQP
       IF (REAL(SQP) .GT. 0) ZP=ZP-2.d0*SQP
       ZP = ZP/(2.d0*COMPLEX(1.d0,PAR(4)))
C
       ZPR = REAL(ZP)
       ZPI = AIMAG(ZP)
C      Boundary conditions
       FB(1)=U0(1) - ZMR 
       FB(2)=U0(2) - ZMI
       FB(3)=U1(4) - ZPR
       FB(4)=U1(5) - ZPI
       FB(5)=U1(1) - U0(4) - PAR(5)
       FB(6)=U1(2) - U0(5) - PAR(6)
       FB(7)=U1(3) - U0(6) - PAR(7)
C 
      RETURN 
      END 
C 
      SUBROUTINE ICND
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      DIMENSION U(NDIM),UOLD(NDIM),UDOT(NDIM),UPOLD(NDIM)
C      DIMENSION FI(NINT),DINT(NINT,*),ICP(*),PAR(*)

C
      RETURN 
      END 
C 
      SUBROUTINE FOPT 
      RETURN 
      END 
C 
      SUBROUTINE PVLS
      RETURN 
      END 

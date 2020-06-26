C----------------------------------------------------------------------
C---------------------------------------------------------------------- 
C  het: Finding heteroclinic connections for trigger fronts in CGL
C---------------------------------------------------------------------- 
C---------------------------------------------------------------------- 
C DEFINE A HEAVISIDE FUNCTION
C      FUNCTION HEAVISIDE(X) 
C      USE NUMERICAL
C      IMPLICIT NONE
C      REAL (KIND=NPRECISION) :: HEAVISIDE
C      REAL (KIND=NPRECISION) :: X 
!
C      IF(X .LT. 0.d0) THEN
C          HEAVISIDE=0.d0
C        ELSE IF(X .GT. 0.d0) THEN
C          HEAVISIDE=1.d0
C        ELSE
C          HEAVISIDE=HALF
C      END IF
!
C      END FUNCTION HEAVISIDE
C 
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
C     ---------- ---- 
C CCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION U(NDIM),PAR(*),F(NDIM),DFDU(NDIM,NDIM),DFDP(NDIM,*)
      DOUBLE PRECISION PI,LE,W,G,C,K,RM,RP,MU1R,MU1I,MU2,UW,UC
      DOUBLE PRECISION UG,MM,MP,UUWM,UUCM,UUG,UUWP,UUCP, AMP, CHI
      DOUBLE PRECISION MU TAU RHO X SH KK RR
C 
       PI=2.*ACOS(0.0d0)

C      COEFFICIENT OF LINEAR TERM ON RIGHT SIDE< MAKES numerics closer to absolute spectrum
       AMP = 1.0d0
C     Variables are in unscaled numerical form:
C     MU is linear onset parameter
C     TAU is homotope parameter to set up heaviside function
C     LE is domain length (scaled in)
       W=PAR(1)
       C= PAR(2)
       G = PAR(3)
       A = PAR(4)
       B=PAR(5)
       MU = PAR(6)
       TAU = PAR(7)
       LE = PAR(8)
       RHO = PAR(9)
       SH = PAR(10)
       KK = PAR(11)
       RR = PAR(12)
CCCCCCCCCCC
CCCCCCCCCCC
C     First set of scalings (tilde to hat)
C       UW = PAR(1)+PAR(4) - 0.5d0*PAR(4)*PAR(2)**2/(1.d0+PAR(4)**2)
C       UC = PAR(2)/sqrt(1.d0+PAR(4)**2)
C       UG = PAR(3)
C

C
C
       X = U(4)
C       CHI=1.d0
       CHI = 1.d0-5.5d0*TAU*(TANH(10000.d0*(X-0.05d0))+1.d0)
C      ODE: 1-3 left side, 1 real part, 2 Im part, 3 radius, 4 is time (deal with non-autonomy)
       F(1)=-(U(1)**2 - U(2)**2)-C*(U(1)+A*U(2))/(1+A**2.d0)
       F(1)=F(1)+(A*W-CHI-(MU+A*G)*U(3))/(1.d0+A**2.d0)
       F(1)=F(1)+((1.d0+A*B)*U(3)**2.d0)/(1.d0+A**2.d0)
       F(2)=-2.d0*U(1)*U(2)
       F(2)=F(2)+(A*CHI+W-C*(U(2)-A*U(1))+(-G+A*MU)*U(3))/(1.d0+A**2.0)
       F(2)=F(2)+((B-A)*U(3)**2.d0)/(1.d0+A**2.d0)
       F(3)=2.d0*U(3)*U(1)
       F(4) =1.d0
C       F(4)=-(U(4)**2-U(5)**2)- C*(U(4)+A*U(5))/(1.d0+A**2)
C       F(4)=F(4)+(A*W + 1.d0*AMP +(1.d0+A*G)*U(6))/(1.d0+A**2)
C       F(5)=-2.d0*U(4)*U(5)
C       F(5)=F(5)+(W-A*AMP-C*(U(5)-A*U(4))+(G-A)*U(6))/(1.d0+A**2)
C       F(6)=2.d0*U(6)*U(4)
C
       F(1) = F(1)*LE
       F(2) = F(2)*LE
       F(3)= F(3)*LE
       F(4)=F(4)*LE
C       F(4) = F(4)*LE
C       F(4)= F(4)*LE
C       F(5)= F(5)*LE
C       F(6)= F(6)*LE
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
      DOUBLE PRECISION PS1,PS2,PHIS,EXPS,TANN,BOT
      COMPLEX ZP, PS
C 
       PI=2*ACOS(0.0d0)
C      Variables in unscaled form. Be sure to imput the numerical form of the variable.
C Start with all complex parameters equal zero.  Want to solve real equation & find pushed front.
       PAR(1)=0.0d0
       PAR(3)= 0.0d0
       PAR(4)= 0.0d0
       PAR(5) = 0.0d0
       PAR(6) = 4.0d0
       PAR(2)=(-PAR(6)+2.d0*SQRT(PAR(6)**2.d0+4.d0))/SQRT(3.d0)
       PAR(7)=0.d0
       PAR(8) = 200.d0
C  Have to set RHO to something (?)
       PAR(9) = -0.05d0
C     Shift point xi=0 for initial profile (in AUTO time, not actual (U(4) time)
       PAR(10)=0.5
       PAR(11)=0.d0
       PAR(12)=0.5d0*(PAR(6)+SQRT(PAR(6)**2.d0+4.d0))
       



       SQQ=SQRT(3.d0)

C      Initial Conditions: Pushed solution from Van Saarloos
       PHIS =0.5d0*(PAR(6)+SQRT(PAR(6)**2.d0+4.d0))
       EXPS= EXP(2.d0*PHIS*(T-PAR(10))*PAR(8)/SQRT(3.d0))
       BOT=SQRT(3.d0)*(EXPS+1/PHIS)
       U(1)=-EXPS*PHIS/BOT
C       U(1)=-PHIS/(1.d0+EXPS)*EXPS/SQQ
C       TANN = TANH(PHIS*T*PAR(8)/SQRT(3.d0))
C       U(1)=-PHIS*(1.d0+TANH(PHIS*(T-0.01)*PAR(8)/SQQ))/SQRT(12.d0)
       U(2)=0.0d0
C       U(3)=PHIS/(1+EXPS)
       U(3)=SQRT(3.d0)/BOT
       U(4)=(T-PAR(10))*PAR(8)

 
      RETURN 
      END 
C  
      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
C     ---------- ---- 
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION PAR(*),ICP(*),U0(NDIM),U1(NDIM),FB(NBC),DBC(NBC,*)
      DOUBLE PRECISION PI,LE,K,SQQ1,PS1,PS2
      DOUBLE PRECISION SQQ,ZPI,ZPR,ZMI,ZMR,CHI
      COMPLEX*16 ZP
      COMPLEX*16 PS 
      COMPLEX*16 RHSS
C
C       Boundary Value parameters
C      Periodic Orbit in Write-up Variables
C
C  NOT RIGHT!! disp. rel. for cubic not quintic!!!!>>!
C       IF (PAR(3) .EQ. PAR(4)) THEN
C          K = -PAR(4)/SQRT(1+PAR(4)**2.d0)
C       ELSE
C          K =-PAR(2)
C          K =K+SQRT(PAR(2)**2 + 4*(PAR(3)+PAR(1))*(PAR(3) - PAR(4)))
C          K =K/(2.d0*(PAR(3) - PAR(4)))
C       ENDIF

       ZMR = 0.d0
       ZMI= PAR(11)
C Make sure in bar variables!!!!
C For plus equations use P, for minus use M!!!
C       PS1= PAR(2)**2*(1.d0-PAR(4)**2)
C       PS1=PS1+ 4.d0*(1.d0+PAR(4)**2)*(1.d0+PAR(4)*PAR(1))
       CHI = 1.d0-11.d0*PAR(7)
       PS1=PAR(2)**2.d0-4.d0*(PAR(4)*PAR(1)+CHI)
C       PS1 = PAR(2)**2.d0-4.d0*CHI
C
C       PS2=-2.d0*PAR(4)*PAR(2)**2
C       PS2=PS2+4.d0*(1+PAR(4)**2)*(-PAR(4)+PAR(1))
       PS2=-4.d0*(-PAR(1)+PAR(4)*CHI)
C       PS2 = 4.d0*PAR(1)
C
       PS = CMPLX(PS1,PS2)
       ZP = CMPLX(-PAR(2),0.d0)-CDSQRT(PS)
C       IF (PAR(1) .LT. -PAR(4)) ZP= ZP+2.d0*sqrt(PS)
       ZP = ZP/CMPLX(2.d0,2.d0*PAR(4))
C
       ZPR = REAL(ZP)
C       ZPR=-(PAR(6)+SQRT(PAR(6)**2.d0 +4.d0))/SQRT(12.d0)
C       ZPR=(-PAR(2)-SQRT(PAR(2)**2-4*CHI))/2.D0
       ZPI = IMAG(ZP)
C       ZPI=0.d0
C
       RHSS=PAR(4)*PAR(11)**2-PAR(3)*PAR(12)+PAR(5)*PAR(12)**2
C      Boundary conditions
       FB(1)=U0(1) - ZMR 
       FB(2)=U0(2) -ZMI 
       FB(3)=U1(1) - ZPR
       FB(4)=U1(2) - ZPI
C       FB(5)=U0(4) + PAR(8)*(-0.05)
       FB(5)=U1(4) - PAR(8)*(1-PAR(10))
       FB(6)=PAR(11)**2-1.d0-PAR(6)*PAR(12)+PAR(12)**2
       FB(7)=RHSS-(-PAR(1)+PAR(2)*PAR(11))
C       FB(5)=U0(4) - 0.d0
C       FB(6)=U1(4) - 1.d0

C No more matching coditions in the middle (Replaced by a phase condition, and heaviside function)
C       FB(5)=U1(1) - U0(4) - PAR(5)
C       FB(6)=U1(2) - U0(5) - PAR(6)
C       FB(7)=U1(3) - U0(6) - PAR(7)
C 
      RETURN 
      END 
C 
      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT) 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(NDIM),UOLD(NDIM),UDOT(NDIM),UPOLD(NDIM)
      DIMENSION FI(NINT),DINT(NINT,*),ICP(*),PAR(*)
      DOUBLE PRECISION SR RRS PIS
      PIS = 2.d0*ACOS(0.0d0)
C   Integral condition for phase condition
      SR = 0.6d0
      RRS = 1/(SR*SQRT(2.d0*PIS))
C      FI(1)=U(3)*RRS*EXP(-(U(4)+0.05d0)**2.d0/(2.d0*SR**2.d0))-PAR(9)
      FI(1)=U(3)*RRS*EXP(-(U(4)-PAR(9))**2.d0/(2.d0*SR**2))-0.0076857
C   -0.0076857 computed from mathematica for PAR(8) = 200
C     -0.0051238 computed from mathemmatica for PAR(8) = 300
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

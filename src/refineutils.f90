!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
!
! Version: :VERSION:
! Date:    :DATE:
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build:   :BUILD:
! Author:  :AUTHOR:
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of felixsim.
!
!  felixsim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixsim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REAL(RKIND) FUNCTION PhaseCorrelate(RImageSim,RImageExpiDummy,IErr,IXsizeIn,IYSizeIn)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  USE MPI
  USE MyMPI
  USE MyFFTW

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,IXsizeIn,IYSizeIn

  REAL(RKIND),DIMENSION(IXSizeIn,IYSizeIn) :: &
       RImageExpiDummy,RImageSim
  
  type(C_PTR) :: &
       Iplan

  real(C_DOUBLE), pointer :: &
       RImageSimDummy(:,:)
  type(C_PTR) :: & 
       p1,p2,p3,p4
  complex(C_DOUBLE_COMPLEX), pointer :: &
       CDummy1(:,:),CDummy2(:,:),CCorrelatedImage(:,:)
  integer(C_INT) :: IX,IY
  
  IX = IXSizeIn
  IY = IYSizein

  !PRINT*,"IX, IY =",IX,IY

  p1 = fftw_alloc_real(INT(IXSizeIn*IYSizeIn, C_SIZE_T))
  p2 = fftw_alloc_complex(INT(IXSizeIn*IYSizeIn, C_SIZE_T))
  p3 = fftw_alloc_complex(INT(IXSizeIn*IYSizeIn, C_SIZE_T))
  p4 = fftw_alloc_complex(INT(IXSizeIn*IYSizeIn, C_SIZE_T))
  call c_f_pointer(p1, RImageSimDummy, [IXSizeIn,IYSizeIn])
  call c_f_pointer(p2, CDummy1, [IXSizeIn,IYSizeIn])
  call c_f_pointer(p3, CDummy2, [IXSizeIn,IYSizeIn])
  call c_f_pointer(p4, CCorrelatedImage, [IXSizeIn,IYSizeIn])
  !...use arr and arr(i,j) as usual...
  
  ! Set the dummy array to the input simulated data

  RImageSimDummy = RImageSim


  !PRINT*,RImageSimDummy(:2,:2)

  ! Plan and Execute the fft of the Simulated Data 
  
  Iplan = FFTW_PLAN_DFT_r2c_2D(IX,IY,RImageSimDummy,CDummy1,FFTW_ESTIMATE)
  CALL FFTW_EXECUTE_DFT_R2C(Iplan,RImageSimDummy,CDummy1)
  CALL FFTW_DESTROY_PLAN(Iplan)

  ! Set the dummy array to the input experimental data

  RImageSimDummy = RImageExpiDummy 


  !PRINT*,RImageSimDummy(:2,:2)
 
  ! Plan and Execute the fft of the Experimental Data 

  Iplan = FFTW_PLAN_DFT_R2C_2D(IX,IY,RImageSimDummy,CDummy2,FFTW_ESTIMATE)

  CALL FFTW_EXECUTE_DFT_R2C(Iplan,RImageSimDummy,CDummy2)
  CALL FFTW_DESTROY_PLAN(Iplan)

  !Calculate the Phase Correlation

  WHERE(ABS(CDummy1*CONJG(CDummy2)).NE.ZERO)
     CCorrelatedImage = (CDummy1*CONJG(CDummy2))/&
          (&
          ABS(CDummy1*CONJG(CDummy2)))
  ELSEWHERE 
     CCorrelatedImage = CZERO
  END WHERE
  
  ! Plan and Execute the inverse fft of the phase correlation

  Iplan = FFTW_PLAN_DFT_C2R_2D(IX,IY,CCorrelatedImage,RImageSimDummy,FFTW_ESTIMATE)

  CALL FFTW_EXECUTE_DFT_C2R(Iplan,CCorrelatedImage,RImageSimDummy)

  CALL FFTW_DESTROY_PLAN(Iplan)

  
!!$  RCrossCorrelation = MAXVAL(RImageSimDummy)/(IX*IY)
  PhaseCorrelate = MAXVAL(RImageSimDummy)/(IX*IY)
  IOffset = MAXLOC(RImageSimDummy)
  
  !PRINT*,RImageSimDummy(:2,:2)

  call fftw_free(p1)
  call fftw_free(p2)
  call fftw_free(p3)
  call fftw_free(p4)
  
END FUNCTION  PhaseCorrelate

SUBROUTINE ReSortUgs( ISymmetryIntegers,CUgs, N )
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER (IKIND) N,IDummy,ISymmetryIntegers(N,2)
  REAL(RKIND) RhklarraySearch(THREEDIM), RhklarrayCompare(THREEDIM)
  COMPLEX(CKIND) CUgSearch,CUgCompare,CUgs(N)
  REAL(KIND=RKIND) ALN2I, LocalTINY
  PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
  
  INTEGER (IKIND) NN,M,L,K,J,I,LOGNB2, index
  COMPLEX(CKIND) Cdummy

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ReSort()"
  END IF
  
  LOGNB2=INT(LOG(REAL(N))*ALN2I+LocalTINY)
  M=N
  DO 12 NN=1,LOGNB2
     M=M/2
     K=N-M
     DO 11 J=1,K
        I=J
3       CONTINUE
        L=I+M
        
        CUgSearch = CUgs(L)
        CUgCompare = CUgs(I)
        IF( &
             (REAL(CUgSearch**2)) .GT. &
             (REAL(CUgCompare**2))) THEN
 !          DO 100
              !IF(my_rank.eq.0) THEN
              !   PRINT*,I
              !END IF
              Cdummy = CUgs(I)
              CUgs(I)= CUgs(L)
              Cugs(L)= Cdummy
              Idummy = ISymmetryIntegers(I,2)
              ISymmetryIntegers(I,2)= ISymmetryIntegers(L,2)
              ISymmetryIntegers(L,2)= Idummy
!100        ENDDO
           
           I=I-M
           IF(I.GE.1) GOTO 3
        ENDIF
11   ENDDO
12 ENDDO
  
  !PRINT*,"Finishing ResortHKL"

  !	PRINT*,"array0(1),array0(N)",array0(1),array0(N)
  RETURN

END SUBROUTINE ReSortUgs

REAL(RKIND) FUNCTION ResidualSumofSquares(RImage1,RImage2,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       IErr
  REAL(RKIND),DIMENSION(2*IPixelCount,2*IPixelCount) :: &
       RImage1,RImage2
!!$  REAL(RKIND) :: &
!!$       ResidualSumofSquares


!!$  RImage1 = RImage1*(2.0**16.0)
  RImage2 =  RImage2/(2.0**16.0)

  PRINT*,"Residual Sum of Squares Before",ResidualSumofSquares,MAXVAL(RImage1),MAXVAL(RImage2)

  ResidualSumofSquares = SUM((RImage2-RImage1)**2)

  PRINT*,"Residual Sum of Squares After",ResidualSumofSquares

END FUNCTION ResidualSumofSquares

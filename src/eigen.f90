!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
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
!  This file is part of FelixSim.
!
!  FelixSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  FelixSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with FelixSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: eigen.f90,v 1.10 2014/03/25 15:35:34 phsht Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE EigenSpectrum(IMatrixDimension, MatrixToBeDiagonalised, EigenValues, EigenVectors, IErr)

  USE WriteToScreen
  USE MyNumbers

  USE IConst
  USE IPara

  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IMatrixDimension, IErr
  COMPLEX(RKIND) :: &
       MatrixToBeDiagonalised(IMatrixDimension,IMatrixDimension), &
       EigenValues(IMatrixDimension), EigenVectors(IMatrixDimension,IMatrixDimension)
  INTEGER(IKIND) :: &
       WorkSpaceDimension
  COMPLEX(CKIND),DIMENSION(:), ALLOCATABLE :: &
       CWorkSpace
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: &
       WorkSpace
  EXTERNAL ZGEEV



  IF (my_rank.EQ.0.AND.IMinWeakBeams.NE.0) THEN
     DO WHILE (IMessageCounter .LT.8)
        CALL Message("EigenSpectrum",IMust,IErr)
        IMessageCounter = IMessageCounter +1
     END DO
  ELSE IF (my_rank.EQ.0) THEN
     DO WHILE (IMessageCounter .LT.6)
        CALL Message("EigenSpectrum",IMust,IErr)
        IMessageCounter = IMessageCounter +1
     END DO
  END IF

  ! ------------------------------------------------
  ! find optimum size of arrays
  ! ------------------------------------------------

  WorkSpaceDimension=1

  ALLOCATE(&
       CWorkSpace(WorkSpaceDimension), &
       STAT = IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (query stage)"
     RETURN
  ENDIF

  ALLOCATE(&
       WorkSpace(2*IMatrixDimension), &
       STAT = IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (query stage)"
     RETURN
  ENDIF

  WorkSpaceDimension=-1

  CALL ZGEEV('N','V', IMatrixDimension, MatrixToBeDiagonalised, IMatrixDimension,&
       EigenValues, 0,1, EigenVectors,IMatrixDimension, &
       CWorkSpace, WorkSpaceDimension, WorkSpace, IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ZGEEV determining work arrays"
     RETURN
  ENDIF

  WorkSpaceDimension = INT(CWorkSpace(1))

  ! ------------------------------------------------
  ! ALLOCATE necessary memory
  ! ------------------------------------------------

  DEALLOCATE(&
       CWorkSpace,&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (final stage)"
     RETURN
  ENDIF

  ALLOCATE(&
       CWorkSpace(WorkSpaceDimension),&
       STAT = IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (final stage)"
     RETURN
  ENDIF

  ! ------------------------------------------------
  ! do the actual call to get the spectrum
  ! ------------------------------------------------

  CALL ZGEEV('N','V', IMatrixDimension, MatrixToBeDiagonalised, IMatrixDimension,&
       EigenValues, 0,1, EigenVectors,IMatrixDimension, &
       CWorkSpace, WorkSpaceDimension, WorkSpace, IErr )
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error ", IErr, " in ZGEEV"
     RETURN
  ENDIF

  DEALLOCATE(&
       CWorkSpace,&
       STAT = IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (final stage)"
     RETURN
  ENDIF

  DEALLOCATE(&
       WorkSpace,&
       STAT = IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"EigenSpectrum: error in ALLOCATE() for work arrays (final stage)"
     RETURN
  ENDIF

  RETURN

END SUBROUTINE EigenSpectrum

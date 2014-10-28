!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all rights reserved
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

SUBROUTINE GMatrixInitialisation (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr,IUniqueKey,knd,IFound

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"GMatrixInitialisation()"
  END IF

  DO ind=1,nReflections
     DO jnd=1,nReflections
        RgMatMat(ind,jnd,:)= RgVecMatT(ind,:)-RgVecMatT(jnd,:)
        RgMatMag(ind,jnd)= SQRT(DOT_PRODUCT(RgMatMat(ind,jnd,:),RgMatMat(ind,jnd,:)))
     ENDDO
  ENDDO
   
  RgMatMag = RgMatMag/TWOPI
  
END SUBROUTINE GMatrixInitialisation

SUBROUTINE DetermineSymmetryRelatedUgs (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI
  USE CPara
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr,knd,Iuid


  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DetermineSymmetryRelatedUgs(",my_rank,")"
  END IF

  !Immediately set all the zeros to Relation 1
  
  ISymmetryRelations = 0
  
  Iuid = 0
  
  Iuid = Iuid + 1

  WHERE (ABS(CUgMat).LE.RTolerance)
     ISymmetryRelations = Iuid
  END WHERE
  
  DO ind = 1,nReflections
     DO jnd = 1,ind
        IF(ISymmetryRelations(ind,jnd).NE.0) THEN
           CYCLE
        ELSE
           Iuid = Iuid + 1
           WHERE (ABS(ABS(CUgMat)-ABS(CUgMat(ind,jnd))).LE.RTolerance)
              ISymmetryRelations = Iuid
           END WHERE
        END IF
     END DO
  END DO

  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"Unique Ugs = ",Iuid
  END IF

  ALLOCATE(&
       ISymmetryStrengthKey(Iuid,2),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"UgCalculation(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE(&
       CSymmetryStrengthKey(Iuid),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"UgCalculation(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
END SUBROUTINE DetermineSymmetryRelatedUgs

!---------------------------------------------------------------------
SUBROUTINE UgCalculation (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,ierr,imaxj,IFound,ICount,currentatom
  INTEGER(IKIND),DIMENSION(2) :: &
       IPos
  COMPLEX(CKIND) CVgij
  REAL(RKIND) :: &
       RMeanInnerPotentialVolts,RAtomicFormFactor  

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"UgCalculation()"
  END IF
  
  DO ind=1,nReflections
     imaxj = ind
     DO jnd=1,imaxj
        
        CVgij= 0.0D0
        
        DO IAtom=1, INAtomsUnitCell
           ICurrentAtom = IAtoms(IAtom)
           ! calculate f_e(q) as in Eq. (C.15) of Kirkland, "Advanced Computing in EM"
           
           SELECT CASE (IScatterFactorMethodFLAG)
              
           CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians 
              
              RAtomicFormFactor = &
                   ! 3 Lorentzians
                   RScattFactors(ICurrentAtom,1) / (RgMatMag(ind,jnd)**2 + RScattFactors(ICurrentAtom,2)) + &
                   RScattFactors(ICurrentAtom,3) / (RgMatMag(ind,jnd)**2 + RScattFactors(ICurrentAtom,4)) + &
                   RScattFactors(ICurrentAtom,5) / (RgMatMag(ind,jnd)**2 + RScattFactors(ICurrentAtom,6)) + &
                   ! 3 Gaussians
                   RScattFactors(ICurrentAtom,7) * EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(ICurrentAtom,8)) + &
                   RScattFactors(ICurrentAtom,9) * EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(ICurrentAtom,10)) + &
                   RScattFactors(ICurrentAtom,11) * EXP(-RgMatMag(ind,jnd)**2 * RScattFactors(ICurrentAtom,12))
              
           CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
              RAtomicFormFactor = &
                   RScattFactors(ICurrentAtom,1) * EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,5)) + &
                   RScattFactors(ICurrentAtom,2) * EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
                   RScattFactors(ICurrentAtom,3) * EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,7)) + &
                   RScattFactors(ICurrentAtom,4) * EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,8))

           CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)

              RAtomicFormFactor = &
                   RScattFactors(ICurrentAtom,1) * EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,2)) + &
                   RScattFactors(ICurrentAtom,3) * EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,4)) + &
                   RScattFactors(ICurrentAtom,5) * EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
                   RScattFactors(ICurrentAtom,7) * EXP(-(RgMatMag(ind,jnd)**2)/4 * RScattFactors(ICurrentAtom,8))

           END SELECT
              
          ! initialize potential as in Eq. (6.10) of Kirkland

           RAtomicFormFactor = RAtomicFormFactor*ROcc(iAtom)

           IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
              
              IF(RDWF(IAtom).GT.10.0D0.OR.RDWF(iAtom).LE.0.0D0.OR.RDWF(iAtom).LE.TINY) THEN
                 RDWF(IAtom) = RDebyeWallerConstant
              END IF
              
              SELECT CASE (IScatterFactorMethodFLAG)

              CASE (0)

              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-((RgMatMag(ind,jnd)/2.D0)**2)*RDWF(iAtom))
                               
              CASE(1)

              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-((RgMatMag(ind,jnd)/2.D0)**2)*RDWF(iAtom))
                          
                               
              CASE(2)

              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-((RgMatMag(ind,jnd)/2.D0)**2)*RDWF(iAtom))
                               
              END SELECT
              
              !RAtomicFormFactor = RAtomicFormFactor * &
              !     EXP(-RgMatMag(ind,jnd)**2*RDWF(iAtom))
              
           ELSE
              
              RAtomicFormFactor = RAtomicFormFactor * &
                   EXP(-TWOPI*DOT_PRODUCT(RgMatMat(ind,jnd,:), &
                   MATMUL( RAnisotropicDebyeWallerFactorTensor( &
                   IAnisoDWFT(iAtom),:,:), &
                   RgMatMat(ind,jnd,:))))
              
           END IF
           
           CVgij = CVgij + &
                RAtomicFormFactor * &
                EXP(-CIMAGONE* &
                DOT_PRODUCT(RgMatMat(ind,jnd,:), RrVecMat(iAtom,:)) &
                )
        ENDDO

        CUgMat(ind,jnd)=((((TWOPI**2)* RRelativisticCorrection) / &
             (PI * RVolume)) * CVgij)
!!$        CUgMat(ind,jnd) = RRelativisticCorrection*&
!!$             (((RPlanckConstantBar**2)*FOUR*PI)/(2*RElectronMass*RVolume*RElectronCharge))*&
!!$             CVgij
              
     ENDDO
  ENDDO

  RMeanInnerCrystalPotential= REAL(CUgMat(1,1))
  RMeanInnerPotentialVolts = RMeanInnerCrystalPotential*((RPlanckConstantBar**2)/ &
       (TWO*RElectronMass*RElectronCharge))
       

  IF((IWriteFLAG.GE.2.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"UgCalculation(",my_rank,") RMeanInnerCrystalPotential = ",RMeanInnerCrystalPotential,RMeanInnerPotentialVolts
  END IF
  
  DO ind=1,nReflections
     CUgMat(ind,ind)=CUgMat(ind,ind)-RMeanInnerCrystalPotential
  ENDDO

  CUgMat = CUgMat + CONJG(TRANSPOSE(CUgMat))

END SUBROUTINE UgCalculation

SUBROUTINE UgAddAbsorption(IErr)         


  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  INTEGER(IKIND) IErr,ind,jnd,knd
  REAL(RKIND) :: &
       RIntegrationParameterGMagPrime,&
       RAbsorpativeAtomicFormFactor,&
       RAbsorpativeAtomicFormFactorInterval,&
       RVgij
  INTEGER(IKIND),DIMENSION(2) :: &
       IPos
  COMPLEX(CKIND) CVgij

  IF((my_rank.EQ.0.AND.IWriteFLAG.GE.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"UgAddAbsorption(",my_rank,")"
  END IF

  CUgMatPrime = CZERO
    
  SELECT CASE (IAbsorbFLAG)

  CASE(1)

!!$     THE PROPORTIONAL MODEL OF ABSORPTION
     
     CUgMatPrime = CUgMatPrime+(REAL(CUgMat)*(RAbsorptionPercentage/100_RKIND)*CIMAGONE)
     
     DO ind = 1,SIZE(CUgMat,DIM=1)
        CUgMatPrime(ind,ind) = REAL(RMeanInnerCrystalPotential)*(RAbsorptionPercentage/100_RKIND)*CIMAGONE
     END DO
    
  CASE Default 
     RUniqueUgPrimeValues = ZERO

!!$     The Einstein TDS model 
!!$     Equations come from Bird and King 1990 Acta Cryst A46, 202-208 

!!$     Calculate the Mean Absorption
     
!!$     DO IAtom=1, INAtomsUnitCell
!!$        RAbsorpativeAtomicFormFactor = ZERO
!!$        ROuterIntegralLowerBound = 0.0D0
!!$        ROuterIntegralUpperBound = 30.0D0
!!$        CALL IntegrateForMeanAbsorptionValue(RAbsorpativeAtomicFormFactor,IErr)
!!$        CVgij = CVgij + RAbsorpativeAtomicFormFactor
!!$     END DO
     
!!$           Calculate Vg'
     
!!$     PRINT*,REAL((FOUR*TWOPI**2*RPlanckConstantBar**2*RPlanckConstant)/(&
!!$             (TWO*RElectronMass**2*((RElectronVelocity/RSpeedOfLight)*&
!!$             RSpeedofLight * RVolume)))*CVgij,RKIND)/RElectronCharge
     
!!$           Convert Vg' into Ug'
     
!!$     RUniqueUgPrimeValues(1) = REAL(RUniqueUgPrimeValues(1)*&
!!$          (TWO*RElectronMass)/&
!!$          (RPlanckConstantBar**2))
     
!!$           Correct for relativistic effects
     
!!$     RUniqueUgPrimeValues(1) = REAL(RUniqueUgPrimeValues(1)*&
!!$          RRelativisticCorrection,RKIND)
                    
!!$     Calculate the Anomolous Absorption

     DO knd = 1,SIZE(ISymmetryStrengthKey,DIM=1)
        ind = ISymmetryStrengthKey(knd,1)
        jnd = ISymmetryStrengthKey(knd,2)
        RVgij = ZERO
        RSVector = RgMatMat(ind,jnd,:)/(FOUR*PI)
        RSVectorMagnitude = RgMatMag(ind,jnd)/TWO
        
        DO IAtom=1, INAtomsUnitCell
           ICurrentAtom = IAtoms(IAtom)
           
           RAbsorpativeAtomicFormFactor = ZERO
           
           ROuterIntegralLowerBound = 0.0D0
           ROuterIntegralUpperBound = 100.0D0
           
           CALL RIntegrateForAnomolousAbsorption(RAbsorpativeAtomicFormFactor,IErr)
           
           RAbsorpativeAtomicFormFactor=&
                RAbsorpativeAtomicFormFactor*ROcc(IAtom)
           
           RVgij = RVgij + &
                RAbsorpativeAtomicFormFactor * &
                EXP(-CIMAGONE* &
                DOT_PRODUCT(RgMatMat(ind,jnd,:), RrVecMat(IAtom,:)) &
                )*EXP(-RDWF(IAtom)*(RSVectorMagnitude)**2)
        ENDDO
        
!!$           Calculate Vg'
        
        RUniqueUgPrimeValues(knd) = ((TWOPI*FOUR*RPlanckConstantBar**2*RPlanckConstant)/((&
             (TWO*RElectronMass**2*((RElectronVelocity)/RSpeedOfLight)*&
             RSpeedofLight * RVolume))))*RVgij/RElectronCharge
!!$
!!$           RUniqueUgPrimeValues(knd) =REAL(CVgij*&
!!$                ((RPlanckConstantBar**2)/&
!!$                (TWOPI*RElectronMass*RElectronCharge*RVolume)),RKIND)
        
!!$           PRINT*,RUniqueUgPrimeValues(knd)
        
!!$           Convert Vg' into Ug'
        
!!$           RUniqueUgPrimeValues(knd) = REAL(RUniqueUgPrimeValues(knd)*&
!!$                (TWO*RElectronMass)/&
!!$                (RPlanckConstantBar**2),RKIND)

!!$           Correct for relativistic effects

!!$           RUniqueUgPrimeValues(knd) = REAL(RUniqueUgPrimeValues(knd)*&
!!$                RRelativisticCorrection,RKIND)
        
     ENDDO
     
  END SELECT
  
  IErr = 0
  
END SUBROUTINE UgAddAbsorption

REAL(RKIND) FUNCTION RSPrimeIntegration(RThetaIntegralParameterIn)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 

  REAL(RKIND) :: &
       RAbsoluteError
  
  REAL(RKIND),INTENT(IN) :: &
       RThetaIntegralParameterIn
  INTEGER(IKIND) :: &
       IErr,IIntegrationSteps
   
  EXTERNAL RAnomolousAbsorpativeIntegrand

  RThetaIntegralParameter = RThetaIntegralParameterIn

  CALL DQNG(RAnomolousAbsorpativeIntegrand,ROuterIntegralLowerBound,ROuterIntegralUpperBound,&
       0.0D0,1.0D-3,RSPrimeIntegration,RAbsoluteError,IIntegrationSteps,IErr)

END FUNCTION RSPrimeIntegration

SUBROUTINE RIntegrateForAnomolousAbsorption(RAbsorpativeAtomicFormFactor,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 

  INTEGER(IKIND) :: &
       IErr,IIntegrationSteps
  REAL(RKIND),INTENT(OUT) :: &       
       RAbsorpativeAtomicFormFactor
  REAL(RKIND) :: &
       RAbsoluteError
      
  EXTERNAL RSPrimeIntegration

  RAbsorpativeAtomicFormFactor = ZERO

  CALL DQNG(RSPrimeIntegration,0.0D0,TWOPI,&
       0.0D0,1.0D-4,RAbsorpativeAtomicFormFactor,RAbsoluteError,IIntegrationSteps,IErr)

END SUBROUTINE RIntegrateForAnomolousAbsorption

SUBROUTINE IntegrateForMeanAbsorptionValue(RAbsorpativeAtomicFormFactor,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 

  INTEGER(IKIND) :: &
       IErr,IIntegrationSteps
  REAL(RKIND),INTENT(OUT) :: &
       RAbsorpativeAtomicFormFactor
  REAL(RKIND) :: &
       RAbsoluteError
      
  EXTERNAL RMeanAbsorpativeIntegrand
  
  CALL DQNG(RMeanAbsorpativeIntegrand,ROuterIntegralLowerBound,ROuterIntegralUpperBound,&
       0.0D0,1.0D-3,RAbsorpativeAtomicFormFactor,RAbsoluteError,IIntegrationSteps,IErr)
END SUBROUTINE IntegrateForMeanAbsorptionValue

REAL(RKIND) FUNCTION  RMeanAbsorpativeIntegrand(RSPrimeMagnitude)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  REAL(RKIND) :: &
       X

  REAL(RKIND) :: &
       RAtomicFormFactorGMagPrime,&
       RAtomicFormFactorGMagMinusGMagPrime,&
       RAbsorpativeIntegrand,&
       RIntegrationParameterGMagPrime,&
       OneDIntegral,RMeanAbsorpativeIntegrand

  REAL(RKIND),INTENT(IN) :: &
       RSPrimeMagnitude


  SELECT CASE (IScatterFactorMethodFLAG)
     
  CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians 
     RMeanAbsorpativeIntegrand = (&
          RSPrimeMagnitude * &
          ( &
          ! 3 Lorentzians
          RScattFactors(ICurrentAtom,1) / (RSPrimeMagnitude**2 + RScattFactors(ICurrentAtom,2)) + &
          RScattFactors(ICurrentAtom,3) / (RSPrimeMagnitude**2 + RScattFactors(ICurrentAtom,4)) + &
          RScattFactors(ICurrentAtom,5) / (RSPrimeMagnitude**2 + RScattFactors(ICurrentAtom,6)) + &
          ! 3 Gaussians
          RScattFactors(ICurrentAtom,7) * EXP(-(RSPrimeMagnitude**2)* RScattFactors(ICurrentAtom,8)) + &
          RScattFactors(ICurrentAtom,9) * EXP(-(RSPrimeMagnitude**2)* RScattFactors(ICurrentAtom,10)) + &
          RScattFactors(ICurrentAtom,11) * EXP(-(RSPrimeMagnitude**2)* RScattFactors(ICurrentAtom,12)))* &
          ( &
          ! 3 Lorentzians
          RScattFactors(ICurrentAtom,1) / (RSPrimeMagnitude**2 + RScattFactors(ICurrentAtom,2)) + &
          RScattFactors(ICurrentAtom,3) / (RSPrimeMagnitude**2 + RScattFactors(ICurrentAtom,4)) + &
          RScattFactors(ICurrentAtom,5) / (RSPrimeMagnitude**2 + RScattFactors(ICurrentAtom,6)) + &
          ! 3 Gaussians
          RScattFactors(ICurrentAtom,7) * EXP(-(RSPrimeMagnitude**2) * RScattFactors(ICurrentAtom,8)) + &
          RScattFactors(ICurrentAtom,9) * EXP(-(RSPrimeMagnitude**2) * RScattFactors(ICurrentAtom,10)) + &
          RScattFactors(ICurrentAtom,11) * EXP(-(RSPrimeMagnitude**2) * RScattFactors(ICurrentAtom,12))))
        
  CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
     RAtomicFormFactorGMagPrime = &
          RScattFactors(ICurrentAtom,1) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,5)) + &
          RScattFactors(ICurrentAtom,2) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,6)) + &
          RScattFactors(ICurrentAtom,3) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,7)) + &
          RScattFactors(ICurrentAtom,4) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,8))
     
     RAtomicFormFactorGMagMinusGMagPrime = &
          RScattFactors(ICurrentAtom,1) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,5)) + &
          RScattFactors(ICurrentAtom,2) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
          RScattFactors(ICurrentAtom,3) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,7)) + &
          RScattFactors(ICurrentAtom,4) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,8))
     
  CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)
     
     RAtomicFormFactorGMagPrime = &
          RScattFactors(ICurrentAtom,1) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,2)) + &
          RScattFactors(ICurrentAtom,3) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,4)) + &
          RScattFactors(ICurrentAtom,5) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,6)) + &
          RScattFactors(ICurrentAtom,7) * &
          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,8))
     
     RAtomicFormFactorGMagMinusGMagPrime = &
          RScattFactors(ICurrentAtom,1) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,2)) + &
          RScattFactors(ICurrentAtom,3) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,4)) + &
          RScattFactors(ICurrentAtom,5) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
          RScattFactors(ICurrentAtom,7) * &
          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,8))
     
  END SELECT
  
  ! initialize potential as in Eq. (6.10) of Kirkland
  
  IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
     
     RMeanAbsorpativeIntegrand = RMeanAbsorpativeIntegrand * &
          (1.0D0-EXP(-TWO*RDWF(IAtom)*RSPrimeMagnitude**2))
  ELSE
     
     RAbsorpativeIntegrand = RAbsorpativeIntegrand * &
          EXP(-TWOPI*DOT_PRODUCT(RGVector, &
          MATMUL( RAnisotropicDebyeWallerFactorTensor( &
          IAnisoDWFT(IAtom),:,:), &
          RGVector)))
       
  END IF

END FUNCTION RMeanAbsorpativeIntegrand

REAL(RKIND) FUNCTION  RAnomolousAbsorpativeIntegrand(RSPrimeMagnitude)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  REAL(RKIND) :: &
       X

  REAL(RKIND) :: &
       RAbsorpativeIntegrand,&
       RSPrimeMagnitude,RSMinusMagnitude,RSPlusMagnitude
  REAL(RKIND) :: &
       RAnomolousAbsorpativeIntegrand
  REAL(RKIND),DIMENSION(THREEDIM) :: &
       RSPrimeVector

  RSMinusMagnitude = SQRT(&
       ((RSVector(1)/2.0D0)-RSPrimeMagnitude*DCOS(RThetaIntegralParameter))**2+&
       ((RSVector(2)/2.0D0)-RSPrimeMagnitude*DSIN(RThetaIntegralParameter))**2)
  RSPlusMagnitude = SQRT(&
       ((RSVector(1)/2.0D0)+RSPrimeMagnitude*DCOS(RThetaIntegralParameter))**2+&
       ((RSVector(2)/2.0D0)+RSPrimeMagnitude*DSIN(RThetaIntegralParameter))**2) 

  RAbsorpativeIntegrand = ZERO
  
  SELECT CASE (IScatterFactorMethodFLAG)
     
  CASE(0) ! Kirkland Method using 3 Gaussians and 3 Lorentzians 
     
     RAbsorpativeIntegrand = &
          ! 3 Lorentzians
          RSPrimeMagnitude * &
          ( &
          RScattFactors(ICurrentAtom,1) / (RSMinusMagnitude**2 + RScattFactors(ICurrentAtom,2)) + &
          RScattFactors(ICurrentAtom,3) / (RSMinusMagnitude**2 + RScattFactors(ICurrentAtom,4)) + &
          RScattFactors(ICurrentAtom,5) / (RSMinusMagnitude**2 + RScattFactors(ICurrentAtom,6)) + &
          ! 3 Gaussians
          RScattFactors(ICurrentAtom,7) * EXP(-RSMinusMagnitude**2* RScattFactors(ICurrentAtom,8)) + &
          RScattFactors(ICurrentAtom,9) * EXP(-RSMinusMagnitude**2* RScattFactors(ICurrentAtom,10)) + &
          RScattFactors(ICurrentAtom,11) * EXP(-RSMinusMagnitude**2* RScattFactors(ICurrentAtom,12)))* &
          ( &
          ! 3 Lorentzians
          RScattFactors(ICurrentAtom,1) / (RSPlusMagnitude**2 + RScattFactors(ICurrentAtom,2)) + &
          RScattFactors(ICurrentAtom,3) / (RSPlusMagnitude**2 + RScattFactors(ICurrentAtom,4)) + &
          RScattFactors(ICurrentAtom,5) / (RSPlusMagnitude**2 + RScattFactors(ICurrentAtom,6)) + &
          ! 3 Gaussians 
          RScattFactors(ICurrentAtom,7) * EXP(-RSPlusMagnitude**2 * RScattFactors(ICurrentAtom,8)) + &
          RScattFactors(ICurrentAtom,9) * EXP(-RSPlusMagnitude**2 * RScattFactors(ICurrentAtom,10)) + &
          RScattFactors(ICurrentAtom,11) * EXP(-RSPlusMagnitude**2 * RScattFactors(ICurrentAtom,12)))
     
!!$  CASE(1) ! 8 Parameter Method with Scattering Parameters from Peng et al 1996 
!!$     RAtomicFormFactorGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,5)) + &
!!$          RScattFactors(ICurrentAtom,2) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,7)) + &
!!$          RScattFactors(ICurrentAtom,4) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$     RAtomicFormFactorGMagMinusGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,5)) + &
!!$          RScattFactors(ICurrentAtom,2) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,7)) + &
!!$          RScattFactors(ICurrentAtom,4) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$  CASE(2) ! 8 Parameter Method with Scattering Parameters from Doyle and Turner Method (1968)
!!$     
!!$     RAtomicFormFactorGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-(X**2)/4 * RScattFactors(ICurrentAtom,8))
!!$     
!!$     RAtomicFormFactorGMagMinusGMagPrime = &
!!$          RScattFactors(ICurrentAtom,1) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,2)) + &
!!$          RScattFactors(ICurrentAtom,3) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,4)) + &
!!$          RScattFactors(ICurrentAtom,5) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,6)) + &
!!$          RScattFactors(ICurrentAtom,7) * &
!!$          EXP(-((RGVectorMagnitude-X)**2)/4 * RScattFactors(ICurrentAtom,8))
     
  END SELECT
  
!!$  RAbsorpativeIntegrand = RSPrimeMagnitude*RAtomicFormFactorGMagPrime*RAtomicFormFactorGMagMinusGMagPrime
  
  ! initialize potential as in Eq. (6.10) of Kirkland
  
  IF (IAnisoDebyeWallerFactorFlag.EQ.0) THEN
     
!!$     IF(RDWF(IAtom).GT.10.OR.RDWF(IAtom).LE.0.OR.RDWF(iAtom).LE.TINY) THEN
!!$        RDWF(IAtom) = RDebyeWallerConstant
!!$     END IF
     
     RAbsorpativeIntegrand = RAbsorpativeIntegrand * &
          (1.0D0-EXP(-TWO*RDWF(IAtom)*(RSPrimeMagnitude**2-((RSVectorMagnitude)**2))/(TWOPI**2)))
!!$     RAbsorpativeIntegrand = RAbsorpativeIntegrand * (&
!!$          (1.0d0-EXP(-TWO*RDWF(IAtom)*(RSPrimeMagnitude**2-DOT_PRODUCT(RGVector,RSPrimeVector)))))

  ELSE
     
     RAbsorpativeIntegrand = RAbsorpativeIntegrand * &
          EXP(-TWOPI*DOT_PRODUCT(RGVector, &
          MATMUL( RAnisotropicDebyeWallerFactorTensor( &
          IAnisoDWFT(IAtom),:,:), &
          RGVector)))
       
  END IF
  
  RAnomolousAbsorpativeIntegrand = RAbsorpativeIntegrand

END FUNCTION RAnomolousAbsorpativeIntegrand



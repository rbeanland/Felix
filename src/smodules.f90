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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: smodules.f90,v 1.63 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!--------------------------------------------------------------------
MODULE CConst
  CHARACTER*30, PARAMETER :: RStr= "Version :VERSION:, Build :BUILD:"
  CHARACTER*30, PARAMETER :: DStr= "Date :DATE:, :TIME:"
  CHARACTER*30, PARAMETER :: AStr= "Status :RLSTATUS:, Author :AUTHOR:"

  CHARACTER*8 CSpaceGrp(230)

  DATA CSpaceGrp/"P1","P-1","P2","P21","C2","Pm","Pc","Cm",&
       "Cc","P2/m","P21/m","C2/m","P2/c","P21/c","C2/c", &
       "P222","P2221","P21212","P212121","C2221","C222","F222", &
       "I222","I212121","Pmm2","Pmc21","Pcc2","Pma2","Pca21", &
       "Pnc2","Pmm21","Pba2","Pna21","Pnn2","Cmm2","Cmc21","Ccc2", &
       "Amm2","Aem2","Ama2","Aea2","Fmm2","Fdd2","Imm2","Iba2", &
       "Ima2","Pmmm","Pnnn","Pccm","Pban","Pmma","Pnna","Pmna","Pcca", &
       "Pbam","Pccn","Pbcm","Pnnm","Pmmn","Pbcn","Pbca","Pnma","Cmcm", &
       "Cmce","Cmmm","Cccm","Cmme","Ccce","Fmmm","Fddd","Immm","Ibam", &
       "Ibca","Imma","P4","P41","P42","P43","I4","I41","P-4","I-4", &
       "P4/m","P42/m","P4/n","P42/n","I4/m","I41/a","P422","P4212", &
       "P4122","P41212","P4222","P42212","P4322","P43212","I422", &
       "I4122","P4mm","P4bm","P42cm","P42nm","P4cc","P4nc","P42mc", &
       "P42bc","I4mm","I4cm","I41md","I41cd","P-42m","P-42c","P-421m", &
       "P-421c","P-4m2","P-4c2","P-4b2","P-4n2","I-4m2","I-4c2", &
       "I-42m","I-42d","P4/mmm","P4/mcc","P4/nbm","P4/nnc","P4/mbm", &
       "P4/mnc","P4/nmm","P4/ncc","P42/mmc","P42/mcm","P42/nbc","P42/nnm", &
       "P42/mbc","P42/mnm","P42/mnc","P42/ncm","I4/mmm","I4/mcm","I41/amd", &
       "I41/acd","P3","P31","P32","R3","P-3","R-3","P312","P321","P3112", &
       "P3121","P3212","P3221","R32","P3m1","P31m","P3c1","P31c","R3m", &
       "R3c","P-31m","P-31c","P-3m1","P-3c1","R-3m","R-3c","P6","P61", &
       "P65","P62","P64","P63","P-6","P6/m","P63/m","P622","P6122", &
       "P6522","P6222","P6422","P6322","P6mm","P6cc","P63cm","P63mc", &
       "P-6m2","P-6c2","P-62m","P-62c","P6/mmm","P6/mcc","P63/mcm", &
       "P63/mmc","P23","F23","I23","P213","I213","Pm-3","Pn-3","Fm-3", &
       "Fd-3","Im-3","Pa-3","Ia-3","P432","P4232","F432","F4132","I432", &
       "P4332","P4132","I4132","P-43m","F-43m","I-43m","P-43n","F-43c", &
       "I-43d","Pm-3m","Pn-3n","Pm-3n","Pn-3m","Fm-3m","Fm-3c","Fd-3m", &
       "Fd-3c","Im-3m","Ia-3d"/
       

END MODULE CConst

!--------------------------------------------------------------------
MODULE IConst
  USE MyNumbers
  INTEGER(IKIND), PARAMETER :: &
       MAXWriteFLAG= 10, &
       THREEDIM= 3, &
       ADD_OUT_INFO=6, &
       IParallelFLAG=0,&
       IRandomFLAG = 0, &
       IRefinementVariableTypes = 8,&
       NElements=103

  !PriorityFLAG values - to match to the WriteFLAG - will change eventually,
  !hence why the silent & Must are both 0, no Silent option yet.
  !IInfo is now IWriteFLAG = 1, IAllInfo is IWriteFLAG = 10
  INTEGER(IKIND), PARAMETER :: &
       ISilent = 0 ,&
       IMust = 1 ,&
       IInfo = 2 ,&
       IMoreInfo = 3, &
       IAllInfo = 10, &
       IDebug = 100, &
       IWarning = 4, &
       IPotError = 5, &
       ICritError = 6

END MODULE IConst

!--------------------------------------------------------------------
MODULE RConst
  USE MyNumbers
  
  REAL(RKIND), PARAMETER :: &
       RSpeedOfLight=REAL(2.99762458D+8,RKIND), &
       RElectronMass=REAL(9.10938291D-31,RKIND), &
       RElectronMassMeV=REAL(0.510998928,RKIND), &
       RPlanckConstant=REAL(6.62606957D-34,RKIND), &
       RElectronCharge=REAL(1.602176565D-19,RKIND), &
       RAngstromConversion=REAL(1.D10,RKIND)
  REAL(RKIND), PARAMETER :: &
       RTolerance =REAL( 1E-5,RKIND)
    
END MODULE RConst

!--------------------------------------------------------------------
MODULE IPara
  USE MyNumbers
  USE IConst 
  
  !Write Out
  
  INTEGER(IKIND) :: &
       IMAXRBuffer,  IMAXCBuffer     
  
  !Input Flags

  INTEGER(IKIND) :: &
       IWriteFLAG,IDebugFLAG,IScatterFactorMethodFLAG, &
       ICentralBeamFLAG, IMaskFLAG, IVolumeFLAG, &
       IZolzFLAG,IAbsorbFLAG, IAnisoDebyeWallerFactorFlag, &
       IImageFLAG,IOutputFLAG,IBeamConvergenceFLAG,  &
       IPseudoCubicFLAG,IXDirectionFLAG,IBinorTextFLAG, IDevFLAG, &
       IRefineModeFLAG,ISoftwareMode,IHKLSelectFLAG,IPrint,IRefineSwitch,&
       IWeightingFLAG,IContinueFLAG

  !Minimum Reflections etc
  INTEGER(IKIND) :: &
       IMinReflectionPool,IMinStrongBeams,IMinWeakBeams

  !OtherFLAGS

  INTEGER(IKIND) :: &
       IDiffractionFLAG=0

  !Disk Radius

  INTEGER(IKIND) :: &
       IPixelCount
  
  !Crystal Settings

  INTEGER(IKIND) :: &
       ITotalAtoms

  ! Name2Atom index
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IAtomNumber,IAtoms

  !Microscope Settings

  INTEGER(IKIND) :: &
       IIncidentBeamDirectionX, IIncidentBeamDirectionY, &
       IIncidentBeamDirectionZ, &
       IXDirectionX, IXDirectionY, IXDirectionZ, &
       INormalDirectionX,INormalDirectionY,INormalDirectionZ

  !Iterative Ug

  INTEGER(IKIND) :: &
       INoofUgs

  !LACBED Input

  INTEGER(IKIND) :: &
       IReflectOut

  !Beams from selection criteria

  INTEGER(IKIND) :: &  
       nReflections,nStrongBeams,nWeakBeams,nBeams,IHKLMAXValue
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IAnisotropicDWFTensor, IAnisoDWFT

  !Main

  INTEGER(IKIND) :: &
       IPixelTotal, INAtomsUnitCell,IPixelComputed

  INTEGER, DIMENSION(2) :: & 
       IImageSizeXY

  !Refinement FLAGS
  INTEGER(IKIND) :: &
       IImageOutputFLAG
  
  !LACBED

  INTEGER(IKIND),DIMENSION(:,:), ALLOCATABLE :: &
       ILACBEDStrongBeamList, IPixelLocation, ISymmetryRelations,ISymmetryStrengthKey
  INTEGER(IKIND),DIMENSION(:), ALLOCATABLE :: &
       InBeams,IStrongBeamList,IOutputReflections

  !inpcif

  INTEGER(IKIND) :: &
      ISymCount
  
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       IFullAtomNumber, IFullAnisotropicDWFTensor

  INTEGER(IKIND) :: &
       IPixelCountTotal

  !LACBED Writing

  INTEGER(IKIND) :: &
       ISeperateFolderFlag
  ! Thickness loop Variables

  INTEGER(IKIND) :: &
       IThicknessCount

  INTEGER(IKIND),DIMENSION(:,:),ALLOCATABLE :: &
       IPixelLocations

  !Refine Parameters

  INTEGER(IKIND) :: &
       IFluxIterationSteps,IElements

  INTEGER(IKIND), DIMENSION(2) :: &
       IOffset

  INTEGER(IKIND), DIMENSION(:),ALLOCATABLE :: &
       IElementList

  !Ug Calculation

   INTEGER(IKIND) :: &
        ICurrentAtom,IAtom

   !Refine Mode Binary Selection
   
   INTEGER(IKIND),DIMENSION(IRefinementVariableTypes) :: &
        IRefineModeSelectionArray

   !List of Atomic Sites for Refinement

   INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: &
        IAtomicSitesToRefine

   !Simplex Variables

   INTEGER(IKIND) :: &
        IIndependentVariables

   !Iterative Ids

   INTEGER(IKIND),DIMENSION(:,:),ALLOCATABLE :: &
        IIterativeVariableUniqueIDs

   ! Refinement Vectors

   INTEGER(IKIND) :: &
        IAllowedVectors

   INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: &
        IAllowedVectorIDs

   INTEGER(IKIND) :: & 
        IFelixCount,IPreviousPrintedIteration,IStandardDeviationCalls
   
   !Message Counter (Avoid subroutines printing out 
   !entering message more than once)

   INTEGER(IKIND) :: &
        IMessageCounter=0

END MODULE IPara

!--------------------------------------------------------------------
MODULE RPara
  USE MyNumbers
  USE RConst
  USE IConst

  !INPUT Section
  
  !Beam Selection Criteria
  
  REAL(RKIND) :: &
       RBSMaxDeviationPara, RBSMaxGVecAmp, RBSBethePara, &
       RConvergenceTolerance,RBSBmax, RBSPMax
  
  !Crystal Settings
  
  REAL(RKIND) :: &
       RLengthX, RLengthY, RLengthZ, RVolume, &
       RAlpha, RBeta, RGamma, &
       RDebyeWallerConstant,RAbsorptionPercentage
  
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: &
       RIsotropicDebyeWallerFactors, RAtomicSitePartialOccupancy, RDWF, ROcc
  
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RSymVec,RAtomSiteFracCoordVec, MNP,&
       RUniqueKey
  
  REAL(RKIND), DIMENSION(:,:,:), ALLOCATABLE :: &
       RSymMat

  !Microscope Parameters

  REAL(RKIND) :: &
       RConvergenceAngle, RAcceleratingVoltage,RInnerConvergenceAngle

  !LACBED Input

  REAL(RKIND) :: &
       RInitialThickness, &
       RFinalThickness, &
       RDeltaThickness, &
       RInitialDebyeWallerFactor, &
       RFinalDebyeWallerFactor,&
       RDeltaDebyeWallerFactor

  !Iterative Ugs

  REAL(RKIND) :: &
       RPercentageUgChange

  !Debye Waller Factor not sure if we use this 
  REAL(RKIND) :: & 
       RMeanSquaredDisplacement


  !HKL indices 
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RHKL 
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RInputHKLs

  ! scattering factors
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
      RScattFactors 

  ! Microscopy Settings
  REAL(RKIND) :: &
       RElectronVelocity, RElectronWaveLength, &
       RElectronWaveVectorMagnitude, RRelativisticCorrection, &
       RRelativisticMass, RBraggCentral

  ! Crystallography 
  ! Real Space and Reciprocal Lattice Vectors in Orthogonal and Microscope
  ! reference framce
  REAL(RKIND), DIMENSION(THREEDIM) :: &
       RXDirM, RYDirM, RZDirM,& 
       RaVecO, RbVecO, RcVecO, &
       RaVecM, RbVecM, RcVecM, &
       RarVecO, RbrVecO, RcrVecO, &
       RarVecM, RbrVecM, RcrVecM, &
       RXDirC, RZDirC, &
       RNormDirC,RNormDirM
  
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RrVecMat

  REAL(RKIND) :: &
       RBaseVec(THREEDIM,THREEDIM), &
       RInvBaseVec(THREEDIM,THREEDIM)
  
  REAL(RKIND), DIMENSION(:,:,:), ALLOCATABLE :: &
       RAnisotropicDebyeWallerFactorTensor
  
  !Diffraction Pattern Definitions
  
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: &
       RgVecMag, RSg
  
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RgVecMat, RgVecMatT
  
  REAL(RKIND), DIMENSION(THREEDIM,THREEDIM) :: &
       RTMat

  REAL(RKIND) :: &
       RDeltaK, RMinimumGMag,RGVectorMagnitude
  
  REAL(RKIND),DIMENSION(THREEDIM) :: &
       RGVector

  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RGn

  !Image Initialisation
  
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: & 
       Rhklpositions
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RFinalMontageImage

  !Main Program
  
  REAL(RKIND) :: &
       RMeanInnerCrystalPotential
  
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: & 
       RMask, RgMatMag
  
  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: & 
       RgMatMat

  REAL(RKIND) :: &
       ROuterIntegralLowerBound,ROuterIntegralUpperBound,&
       RInnerIntegralLowerBound,RInnerIntegralUpperBound,&
       RInnerIntegrationParameterGMagPrime,&
       ROuterIntegrationParameterGMagPrime
  
  !LACBED Program
  
  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: &
       RFullAtomicFracCoordVec
  REAL(RKIND), DIMENSION(:), ALLOCATABLE :: &
       RFullPartialOccupancy, RFullIsotropicDebyeWallerFactor

  !WaveFunction Arrays

  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RWaveIntensity,RFullWaveIntensity
    
  REAL(RKIND), DIMENSION(:,:,:), ALLOCATABLE :: &
       RIndividualReflections

  !Refinement Variables

  REAL(RKIND), DIMENSION(:,:),ALLOCATABLE :: &
       RImageIn

  REAL(RKIND) :: &
       RCrossCorrelation,RDeltaUgChange,RlowerBoundUgChange,RUpperBoundUgChange

  !Ug' Unique Values
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RUniqueUgPrimeValues

  
  ! Experimental Images for felixrefine

  REAL(RKIND),DIMENSION(:,:,:),ALLOCATABLE :: &
       RImageExpi  

  ! Independent Variable Values
!!$  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
!!$       RIndependentVariableValues

  !Iterative Variable Value
  REAL(RKIND) :: &
       RValue

  !Refinement Vectors

  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RAllowedVectors
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RAllowedVectorMagnitudes
  REAL(RKIND) :: &
       RSimplexLengthScale,RSimplexStandardDeviation,RSimplexMean

  !Refinement Initial Coordinates

  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: &
       RInitialAtomSiteFracCoordVec

  !Weighting Coefficients for figure of merit combination
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: &
       RWeightingCoefficients

  

END MODULE RPara

MODULE CPara

  USE MyNumbers

  COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: &
       CUgMat,CUgMatPrime, CUgMatEffective,CEigenValuesChunk
  COMPLEX(CKIND), DIMENSION(:,:,:), ALLOCATABLE :: &
       CEigenVectorsChunk
  COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: &
       CAlphaWeightingCoefficients, CPsi0,CSymmetryStrengthKey
  COMPLEX(CKIND),DIMENSION(:,:), ALLOCATABLE :: &
       CEigenValueDependentTerms,CInvertedEigenVectors, &
       CBeamProjectionMatrix,CDummyBeamMatrix
  COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: &
       CEigenValues,CGammaValues, CWaveFunctions,CFullWaveFunctions
  COMPLEX(CKIND),DIMENSION(:,:),ALLOCATABLE :: &
       CEigenVectors
  COMPLEX(CKIND), DIMENSION(:,:,:), ALLOCATABLE :: &
       CAmplitudeandPhase

END MODULE CPara

!--------------------------------------------------------------------
MODULE SPara
  USE MyNumbers
  
  CHARACTER*1 :: &
       SSpaceGroupName
  CHARACTER*10 :: &
       SSpaceGrp
  CHARACTER*2, DIMENSION(:), ALLOCATABLE :: &
       SFullAtomicNameVec
  
  CHARACTER*2, DIMENSION(:), ALLOCATABLE :: &
       SAtomName, SMNP

  CHARACTER*1,DIMENSION(:),ALLOCATABLE :: &
       SWyckoffSymbols
  
END MODULE SPara

!--------------------------------------------------------------------
! Input- and Outputchannels
MODULE IChannels
  INTEGER, PARAMETER :: &
       IChInp= 40, &
       IChOutWF= 41, IChOutWI= 42, &
       IChOutEV= 43, IChOutEX= 44, &
       IChOutUM= 45, IChOut=46, &
       IChInImage = 51
  INTEGER :: &
       IChOutWF_MPI, IChOutWI_MPI, &
       IChOutES_MPI, IChOutUM_MPI, &
       IChOut_MPI 
  INTEGER, PARAMETER :: &
       IChOutWFImageReal= 47, IChOutWFImagePhase= 48, &
       IChOutWIImage= 49, MontageOut = 50,IChOutSimplex = 52
END MODULE IChannels

MODULE BlochPara
  
  USE MyNumbers
  USE IConst
  USE MPI
  USE MyMPI
  
  !--------------------------------------------------------------------
  ! eigen problem variables
  !--------------------------------------------------------------------
  
  REAL(RKIND) RBigK
  
  INTEGER(IKIND) IStrongBeamIndex, IWeakBeamIndex
  INTEGER(IKIND),DIMENSION(:), ALLOCATABLE :: &
       IWeakBeamList
  REAL(RKIND),DIMENSION(:), ALLOCATABLE :: &
       RDevPara
  REAL(RKIND), DIMENSION(THREEDIM) :: &
       RTiltedK
  REAL(8), DIMENSION(:), ALLOCATABLE :: &
       RROutArray, RIOutArray
  COMPLEX(CKIND),DIMENSION(:,:), ALLOCATABLE :: &
       CEigenSaveTemp
END MODULE BlochPara

MODULE Refinement

USE MyNumbers

REAL(RKIND),PARAMETER :: &
  RExitCondition = -10000.0,&
  RStayCondition = 10000.0

END MODULE 

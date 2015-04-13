#!/bin/bash

binSIM=felixsim
binREFINE=felixrefine
#binDRAW=felixdraw

inpfile=felix.inp
scafile=felix.sca
hklfile=felix.hkl

#Name of Job
#-----------
Job_Name= #job name here (FS... or FR...)


#Directory Input
#---------------
#input file directory - point to where all input file(s) are located (at the very least the .sca file)
# NEEDS TO EXIST
inputfiledir=$HOME #rest here

# point this to where the felixsim/refine executable is
# NEEDS TO EXIST
binarydir=$HOME #rest here

# settings for executable
exe=${1:felixsim}

# submission directory - point to what you want the PBS working directory to be 
# Doesn't need to exist
submitdir=$HOME/SUBMIT/ #rest here

# tmp directory to run felix in
# Doesn't need to exist
bashtmpdir=$HOME/TEMP/ #rest here


#Input Arguments
#---------------
# settings for parallel submission
cores=${2:-1}
JobId=${3:-1}
wtime=${4:-00:10:00}
queue=${5:devel}
InputImages=${6:-1}

let ranks=${cores}
let nodes=${cores}
let mppwidth=${cores}


if [ ${ranks} -lt 1 ]
then
    ranks=1
fi

if [ ${nodes} -lt 1 ]
then
    nodes=1
fi

if [ ${ranks} -lt 4 ]
then
    taskspernode=${ranks}
else
    taskspernode=4
fi

if [ ${mppwidth} -lt 4 ]
then
    mppwidth=4
fi


#loop starts here for looped input file parameters or CIFs#
#---------------------------------------------------------#

echo "Submitting jobs for ${cores} cores, ${ranks} MPI ranks, with maximum wall time ${wtime}, submitted to queue ${queue}"

[ -d ${submitdir} ] || mkdir -p ${submitdir}
cd ${submitdir}
job_dir=${submitdir}


job_file=`printf "$Job_Name}_${JobId}.sh"` 

echo "Job File location:" ${job_dir}/${job_file}

[ -d ${job_dir} ] || mkdir -p ${job_dir}

cat > ${job_dir}/${job_file} << EOD
#!/bin/bash --login
#PBS -l pvmem=2000mb
#PBS -M $USER@warwick.ac.uk
#PBS -m a
#PBS -r y
#PBS -V
#PBS -k oe
#PBS -j oe

#       The jobname
#PBS -N ${exe}_${JobId}

#       The total number of parallel tasks for your job.
#PBS -l nodes=${nodes}:ppn=1

#       Specify the wall clock time required for your job.
#PBS -l walltime=${wtime}

#       Specify which budget account that your job will be charged to.
##PBS -A midpluswarwick
  
# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=\$(readlink -f \$PBS_O_WORKDIR)

# The base directory is the directory that the job was submitted from.
# All simulations are in subdirectories of this directory.
basedir=\$PBS_O_WORKDIR
echo "basedir=" \${basedir}

# The binary directory is the directory where the scripts are
echo "binarydir=" ${binarydir}

# do the run in a TMPdir for safekeeping
tmpdir=$HOME/RUNS/`basename ${job_file} .sh`

[ -d \${tmpdir} ] || mkdir -p \${tmpdir}

# copy the needed binaries
cp ${binarydir}/${exe} \${tmpdir}/
cp \${basedir}/${ciffile} \${tmpdir}/
cp ${scadir}/${scafile} \${tmpdir}/
cp ${scadir}/${hklfile} \${tmpdir}/

# copy any input images if present; only needed for felixrefine
echo $InputImages
for imgno in {1..$InputImages}
do
 echo \${imgno}
 imgfile=felix.00\${imgno}.img
 cp ${scadir}/\${imgfile} \${tmpdir}/
done

# construct the input files

cd \${tmpdir}
echo $HOSTNAME
pwd

rm -rf $inpfile

touch $inpfile

echo "# Input file for felixrefine version :VERSION: Build :BUILD:"            >> $inpfile
echo "# ------------------------------------"                                  >> $inpfile
echo ""                                                                        >> $inpfile
echo "# ------------------------------------"                                  >> $inpfile
echo "# felixrefine input"                                                     >> $inpfile
echo ""                                                                        >> $inpfile
echo "# control flags"                                                         >> $inpfile
echo "IWriteFLAG                = 1"                                           >> $inpfile
echo "IImageFLAG                = 01"                                          >> $inpfile
echo "IOutputFLAG               = 0"                                           >> $inpfile
echo "IBinorTextFLAG            = 0"                                           >> $inpfile 
echo "IScatterFactorMethodFLAG  = 0"                                           >> $inpfile
echo "ICentralBeamFLAG          = 1"                                           >> $inpfile
echo "IMaskFLAG                 = 0"                                           >> $inpfile
echo "IZolzFLAG                 = 1"                                           >> $inpfile
echo "IAbsorbFLAG               = 1"                                           >> $inpfile
echo "IAnisoDebyeWallerFlag     = 0"                                           >> $inpfile
echo "IBeamConvergenceFLAG      = 1"                                           >> $inpfile
echo "IPseudoCubicFLAG          = 0"                                           >> $inpfile
echo "IXDirectionFLAG           = 1"                                           >> $inpfile
echo ""                                                                        >> $inpfile
echo "# radius of the beam in pixels"                                          >> $inpfile
echo "IPixelCount               = 64"                                          >> $inpfile
echo ""                                                                        >> $inpfile
echo "# beam selection criteria"                                               >> $inpfile
echo "IMinReflectionPool        = 300"                                         >> $inpfile
echo "IMinStrongBeams           = 50"                                          >> $inpfile
echo "IMinWeakBeams             = 5"                                           >> $inpfile
echo "RBSBMax                   = 0.1"                                         >> $inpfile
echo "RBSPMax                   = 0.1"                                         >> $inpfile
echo "RConvergenceTolerance (%) = 1.0"                                         >> $inpfile
echo ""                                                                        >> $inpfile
echo "# crystal settings"                                                      >> $inpfile
echo "RDebyeWallerConstant      = 0.4668"                                      >> $inpfile
echo "RAbsorptionPer            = 2.9"                                         >> $inpfile
echo ""                                                                        >> $inpfile
echo "# microscope settings"                                                   >> $inpfile
echo "ROuterConvergenceAngle    = 3.0"                                         >> $inpfile
echo "RInnerConvergenceAngle    = 0.0"                                         >> $inpfile
echo "IIncidentBeamDirectionX   = 1"                                           >> $inpfile
echo "IIncidentBeamDirectionY   = 1"                                           >> $inpfile
echo "IIncidentBeamDirectionZ   = 1"                                           >> $inpfile
echo "IXDirectionX              = 0"                                           >> $inpfile
echo "IXDirectionY              = -1"                                          >> $inpfile
echo "IXDirectionZ              = 1"                                           >> $inpfile
echo "INormalDirectionX         = 1"                                           >> $inpfile
echo "INormalDirectionY         = 1"                                           >> $inpfile
echo "INormalDirectionZ         = 1"                                           >> $inpfile
echo "RAcceleratingVoltage (kV) = 200.0"                                       >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Image Output Options"                                                  >> $inpfile
echo ""                                                                        >> $inpfile
echo "RInitialThickness        = 1000.0"                                       >> $inpfile
echo "RFinalThickness          = 1000.0"                                       >> $inpfile
echo "RDeltaThickness          = 10.0"                                         >> $inpfile
echo "IReflectOut              = 7"                                            >> $inpfile
echo ""                                                                        >> $inpfile
echo "# felixrefine Input"                                                     >> $inpfile
echo ""                                                                        >> $inpfile
echo "#Refinement Specific Flags"                                              >> $inpfile
echo "IRefineModeFLAG          = 0"                                            >> $inpfile
echo "IWeightingFLAG           = 0"                                            >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Debye Waller Factor Iteration"                                         >> $inpfile
echo ""                                                                        >> $inpfile
echo "RInitialDebyeWallerFactor = 0.1"                                         >> $inpfile
echo "RFinalDebyeWallerFactor = 1.0"                                           >> $inpfile
echo "RDeltaDebyeWallerFactor = 0.1"                                           >> $inpfile
echo "IElementsforDWFchange = {0}"                                             >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Ug Iteration"                                                          >> $inpfile
echo  ""                                                                       >> $inpfile 
echo "INoofUgs                  = 10"                                          >> $inpfile
echo "RLowerBoundUgChange       = 50.0"                                        >> $inpfile
echo "RUpperBoundUgChange       = 50.0"                                        >> $inpfile
echo "RDeltaUgChange            = 50.0"                                        >> $inpfile
echo  ""                                                                       >> $inpfile  
echo "# Structural Refinement"                                                 >> $inpfile
echo  ""                                                                       >> $inpfile  
echo "IAtomicsSites             = (1,2,3,4,5,6,7)"                             >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Refinement Output"                                                     >> $inpfile
echo  ""                                                                       >> $inpfile 
echo  "IPrint                   = 10"                                          >> $inpfile
echo  ""                                                                       >> $inpfile   

cat $inpfile
ls -al \${tmpdir}

# Make sure ranks are numbered appropriately
export MPICH_RANK_REORDER_METHOD=1 
export MPICH_PTL_MATCH_OFF=1
export MPICH_FAST_MEMCPY=1
export IPATH_NO_CPUAFFINITY=1
  
# Launch the parallel job, ${ranks} MPI ranks in ${taskspernode} tasks per node with ${ompthreads} threads
#export OMP_NUM_THREADS=${ompthreads}

echo "--- starting the SIM run"
mpirun ./$exe

echo "--- starting the IMAGE post-processing run"
pwd
ls -l

#convert -size 539x539 -depth 64 -define quantum:format=floating-point -define quantum:scale=65535.0 -endian lsb GRAY:F-WI-M-T01000-P00539-sampleoutput-highq.bin F-WI-M-T01000-P00539-sampleoutput-highq.tif

echo "--- starting the image part"

# copy the result of the run in the detination for safekeeping

targetdir=${submitdir}/`basename ${job_file} .sh`

[ -d \${targetdir} ] || mkdir -p \${targetdir}

cp -vr * \${targetdir}

echo "copying output files from" 
pwd
echo "to \${targetdir}" 

wait
#exit 0
EOD

#job submission arguments, uncomment queue required

chmod 755 ${job_dir}/${job_file}
#(cd ${job_dir} ; qsub -q serial ./${job_file})
#(cd ${job_dir} ; qsub -q parallel ./${job_file})
(cd ${job_dir} ; qsub -q $queue ./${job_file})
#(cd ${job_dir} ; qsub -q taskfarm ./${job_file})
#(cd ${job_dir} ; qsub ./${job_file})

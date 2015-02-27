#!/bin/bash

binSIM=felixsim
binREFINE=felixrefine
#binDRAW=felixdraw

#pytIMAGE=txt2png.py

inpfile=felix.inp
scafile=felix.sca
hklfile=felix.hkl

# point this to where felixsim and felixdraw are
binarydir=`pwd`
#binarydir=$HOME/D-LACBED/EXE

# point to where felix.sca is located
scadir=${binarydir}

submitdir=`pwd`
tmpdir=`pwd`

# settings for parallel submission
cores=${1:-1}
JobId=${2:-1}
InputImages=${3:-1}
wtime=${4:-00:10:00}

let ranks=${cores}
let nodes=${cores}
let mppwidth=${cores}
curr_dir=`pwd`

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

for ((Strongbeamsind=50;Strongbeamsind<=300;Strongbeamsind+25))
{

echo "Submitting jobs for ${cores} cores, ${ranks} MPI ranks, OpenMP threads with maximum wall time ${wtime}"

[ -d ${submitdir} ] || mkdir ${submitdir}
job_dir=${submitdir}

ciffile=felix.cif


#unsure what this does? I'm confused about what job_file is saved as,
# just the string FS_...etc, or with the command printf? - also can you submit 
#a shell script rather than a .pbs file?
job_file=`printf "FS_StrongBeams_${Strongbeamsind}.sh"` 

echo ${job_dir}/${job_file}

[ -d ${job_dir} ] || mkdir ${job_dir}

cat > ${job_dir}/${job_file} << EOD
#!/bin/bash --login
#PBS -l pvmem=400mb
##PBS -M keith.evans@warwick.ac.uk
#PBS -m a
#PBS -r y
#PBS -V
##PBS -k oe
#PBS -j oe

#       The jobname
#PBS -N FelixSim_${Strongbeamsind}

#       The total number of parallel tasks for your job.
#       This is the sum of the number of parallel tasks required by each
#       of the aprun commands you are using. In this example we have
#       ${mppwidth} tasks
##PBS -l mppwidth=${mppwidth}
#PBS -l nodes=${nodes}:ppn=1

#       Specify how many processes per node.
##PBS -l mppnppn=32

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

[ -d \${tmpdir} ] || mkdir \${tmpdir}

#cp ${binarydir}/${binSIM} \${tmpdir}/
#cp ${binarydir}/${binDRAW} \${tmpdir}/
#cp ${binarydir}/${binREFINE} \${tmpdir}/
cp \${basedir}/${ciffile} \${tmpdir}/
cp ${scadir}/${scafile} \${tmpdir}/
#cp ${scadir}/${hklfile} \${tmpdir}/


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
echo "IXDirectionFLAG           = 0"                                           >> $inpfile
echo ""                                                                        >> $inpfile
echo "# radius of the beam in pixels"                                          >> $inpfile
echo "IPixelCount               = 16"                                          >> $inpfile
echo ""                                                                        >> $inpfile
echo "# beam selection criteria"                                               >> $inpfile
echo "IMinReflectionPool        = 1000"                                        >> $inpfile
echo "IMinStrongBeams           = ${Strongbeamsind}"                           >> $inpfile
echo "IMinWeakBeams             = 0"                                           >> $inpfile
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
echo "RInitialThickness        = 800"                                          >> $inpfile
echo "RFinalThickness          = 1200"                                         >> $inpfile
echo "RDeltaThickness          = 10.0"                                         >> $inpfile
echo "IReflectOut              = 15"                                           >> $inpfile
echo ""                                                                        >> $inpfile


cat $inpfile
ls -al \${tmpdir}

# Make sure ranks are numbered appropriately
export MPICH_RANK_REORDER_METHOD=1 
export MPICH_PTL_MATCH_OFF=1
export MPICH_FAST_MEMCPY=1
export IPATH_NO_CPUAFFINITY=1
  
# Launch the parallel job, ${ranks} MPI ranks in ${taskspernode} tasks per node with ${ompthreads} threads
#export OMP_NUM_THREADS=${ompthreads}

#does this need -n etc?
echo "--- starting the SIM run"
mpirun ./felixsim

#echo "--- starting the DRAW run"
#mpirun -n ${cores} \${tmpdir}/${binDRAW}

#echo "--- starting the IMAGE post-processing run"
#pwd
#ls

#echo "--- starting the image part"
# copy the result of the run in the detination for safekeeping

targetdir=$HOME/RUNS/`basename ${job_file} .sh`
#targetdir=\${basedir}/`basename ${job_file} .sh`

[ -d \${targetdir} ] || mkdir \${targetdir}

cp -vr * \${targetdir}

wait
#exit 0
EOD

chmod 755 ${job_dir}/${job_file}
#(cd ${job_dir} ; qsub -q serial ./${job_file})
#(cd ${job_dir} ; qsub -q parallel ./${job_file})
(cd ${job_dir} ; qsub -q devel ./${job_file})     #what does the devel mean - different types of queue? 
#(cd ${job_dir} ; qsub -q taskfarm ./${job_file})
#(cd ${job_dir} ; qsub ./${job_file})
}
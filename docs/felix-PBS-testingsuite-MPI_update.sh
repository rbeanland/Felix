#!/bin/bash

#binSIM=felixsim
#binDRAW=felixdraw

#pytIMAGE=txt2png.py

inpfile=felix.inp
scafile=felix.sca
hklfile=felix.hkl

#cif directories
cifdir1=$HOME/Felix/samples/Si
cifdir2=$HOME/Felix/samples/Si_min
cifdir3=$HOME/Felix/samples/Si_max

# point this to where felixsim and felixdraw are - is this correct
binarydir=$HOME/project/Felix/src

# point to where felix.sca is located
scadir=${binarydir}

submitdir=`pwd`
tmpdir=/tmp/

# settings for parallel submission
cores=${1:-1}
wtime=${2:-00:10:00}

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

if [ ${ompthreads} -lt 1 ]
then
    ompthreads=1
fi

echo "Submitting jobs for ${cores} cores, ${ranks} MPI ranks, ${ompthreads} OpenMP threads with maximum wall time ${wtime}"

[ -d ${submitdir} ] || mkdir ${submitdir}
job_dir=${submitdir}


#loop over cases that will test a variety of parameters




for casenumber in 1:14
do
    
    ciffile=${cifdir1}/felix.cif

#default input
    writeind=2
    imageind=012
    outputind=0
    binorind=0
    scatterind=0
    centralbeamind=1
    maskind=0
    zolzind=1
    absorbind=1
    anisodebyeind=0
    beamconvergeind=1
    pseudoind=0
    xdirectionflagind=1
    pixelind=64
    minreflectionind=200
    strongbeamsind=50
    weakbeamsind=0
    rbsbind=0.1
    rbspind=0.1
    convergetolind=1.0
    debyewallerind=0.4668
    absorptionind=2.9   
    outerind=3.0
    innerind=0.0
    incidentxind=1
    incidentyind=1
    incidentzind=1
    xdirectionxind=0
    xdirectionyind=-1
    xdirectionzind=1
    normalxind=1
    normalyind=1
    normalzind=1
    acceleratingind=200.0
    initialthickind=1000.0
    finalthickind=1000.0
    deltathickind=10.0
    reflectind=15

    case "$casenumber" in
#normal test (as the normal user would use - same as default) 
	1)
	    ;;
#information test(check slightly obscure features - still should run) - need HOLZ stuff
	2) 
	    writeind=104
	    scatterind=1
	    zolzind=0
	    absorbind=0
	    anisodebyeind=1
	    xdirectionflagind=0
	    minreflectionind=1000
	    strongbeamsind=200
	    weakbeamsind=130
	    outerind=5.0
	    innerind=1.0
	    incidentxind=0
	    incidentyind=-1
	    incidentzind=-1
	    normalxind=0
	    normalyind=-1
	    normalzind=-1
	    acceleratingind=180.0
	    initialthickind=500.0
	    finalthickind=1500.0
	    reflectind=100
	    ;;
#Strong and weak beam mismatch
	3)
            #remove hkl file from previous iteration hkl test
	    rm -rf ${hklfile}
	    weakbeamsind=51
	    minreflectionind=52
	    ;;
#Weak beams higher than strong beams
	4)
	    weakbemasind=51
	    \\
#inner convergence outside outer convergence
	    5)
	    innerind=4.0
	    ;;
#test unfeatured flags/rogue values
	6)
	    outerind=0123
	    binorind=1
	    centralbeamind=1
	    maskind=1
	    absorbind=2
	    anisodebyeind=1
	    beamconvergeind=1
	    pseudoind=1
	    ;;
#strange debye waller constant &Absorptionper
	7)
	    debyewallerind=15.369
	    absorptionind=127.5
	    ;;
#random coordinates (see if felix can handle any direction)
	8)
	    incidentxind=$RANDOM
	    incidentyind=$RANDOM
	    incidentzind=$RANDOM
	    normalxind=$RANDOM
	    normalyind=$RANDOM
	    normalzind=$RANDOM
	    ;;
#strange rbsb & rpsp numbers
	9)
	    rbsbind=0.00000001
	    rbspind=10000000.0
	    ;;    
#crazy test
	10)
	    absorbind=3
	    outerind=6.0
	    innerind=1.0
	    writeind=104
	    scatterind=2
	    zolzind=0
	    minreflectionind=50
	    strongbeamsind=15
	    weakbeamsind=7
	    imageind=0{$RANDOM}
	    binorind=1
	    outputind=123
	    pseudoind=1
	    ;;
#mismatch hkl file
	11)
	    xdirectionflagind=0
	    incidentxind=5
	    incidentyind=-17
	    incidentzind=8
	    normalxind=5
	    normalyind=-17
	    normalzind=8
	    ;;
#cif file minimum
	12)
	    ciffile=${cifdir2}/felix.cif
	    ;;
#cif file maximum  
	13)
	    ciffile=${cifdir3}/felix.cif 
	    ;;
#run sample input file
	14)
	    inpfile=${cifdir3}/felix.inp
	    ;;
    esac

    job_file=`printf "FS_Test_Case_%s.sh" "${casenumber}"` 

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
#PBS -N Felix_${cores}_${casenumber}

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
tmpdir=/tmp/`basename ${job_file} .sh`

[ -d \${tmpdir} ] || mkdir \${tmpdir}

#cp ${binarydir}/${binSIM} \${tmpdir}/
#cp ${binarydir}/${binDRAW} \${tmpdir}/
cp ${ciffile} \${tmpdir}/
cp ${scadir}/${scafile} \${tmpdir}/

# construct the input files

cd \${tmpdir}
echo $HOSTNAME
pwd

if [$casenumber -ne 14]
rm -rf $inpfile

touch $inpfile

echo "# Input file for FelixSIM version :VERSION: Build :BUILD:"               >> $inpfile
echo "# ------------------------------------"                                  >> $inpfile
echo ""                                                                        >> $inpfile
echo "# ------------------------------------"                                  >> $inpfile
echo "# felixsim input"                                                        >> $inpfile
echo ""                                                                        >> $inpfile
echo "# control flags"                                                         >> $inpfile
echo "IWriteFLAG                = ${writeind}"                                 >> $inpfile
echo "IImageFLAG                = ${imageind}"                                 >> $inpfile
echo "IOutputFLAG               = ${outputind}"                                >> $inpfile
echo "IBinorTextFLAG            = ${binorind}"                                 >> $inpfile 
echo "IScatterFactorMethodFLAG  = ${scatterind}"                               >> $inpfile
echo "ICentralBeamFLAG          = ${centralbeamind}"                           >> $inpfile
echo "IMaskFLAG                 = ${maskind}"                                  >> $inpfile
echo "IZolzFLAG                 = ${zolzind}"                                  >> $inpfile
echo "IAbsorbFLAG               = ${absorbind}"                                >> $inpfile
echo "IAnisoDebyeWallerFlag     = ${anisodebyeind}"                            >> $inpfile
echo "IBeamConvergenceFLAG      = ${beamconvergeind}"                          >> $inpfile
echo "IPseudoCubicFLAG          = ${pseudoind}"                                >> $inpfile
echo "IXDirectionFLAG           = ${xdirectionflagind}"                        >> $inpfile
echo ""                                                                        >> $inpfile
echo "# radius of the beam in pixels"                                          >> $inpfile
echo "IPixelCount               = ${pixelind}"                                 >> $inpfile
echo ""                                                                        >> $inpfile
echo "# beam selection criteria"                                               >> $inpfile
echo "IMinReflectionPool        = ${minreflectionind}"                         >> $inpfile
echo "IMinStrongBeams           = ${strongbeamsind}"                           >> $inpfile
echo "IMinWeakBeams             = ${weakbeamsind}"                             >> $inpfile
echo "RBSBMax                   = ${rbsbind}"                                  >> $inpfile
echo "RBSPMax                   = ${rbspind}"                                  >> $inpfile
echo "RConvergenceTolerance (%) = ${convergtolind}"                            >> $inpfile
echo ""                                                                        >> $inpfile
echo "# crystal settings"                                                      >> $inpfile
echo "RDebyeWallerConstant      = ${debyewallerind}"                           >> $inpfile
echo "RAbsorptionPer            = ${absorptionind}"                            >> $inpfile
echo ""                                                                        >> $inpfile
echo "# microscope settings"                                                   >> $inpfile
echo "ROuterConvergenceAngle    = ${outerind}"                                 >> $inpfile
echo "RInnerConvergenceAngle    = ${innerind}"                                 >> $inpfile
echo "IIncidentBeamDirectionX   = ${incidentxind}"                             >> $inpfile
echo "IIncidentBeamDirectionY   = ${incidentyind}"                             >> $inpfile
echo "IIncidentBeamDirectionZ   = ${incidentzind}"                             >> $inpfile
echo "IXDirectionX              = ${xdirectionxind}"                           >> $inpfile
echo "IXDirectionY              = ${xdirectionyind}"                           >> $inpfile
echo "IXDirectionZ              = ${xdirectionzind}"                           >> $inpfile
echo "INormalDirectionX         = ${normalxind}"                               >> $inpfile
echo "INormalDirectionY         = ${normalyind}"                               >> $inpfile
echo "INormalDirectionZ         = ${normalzind}"                               >> $inpfile
echo "RAcceleratingVoltage (kV) = ${acceleratingind}"                          >> $inpfile
echo ""                                                                        >> $inpfile
echo "# Image Output Options"                                                  >> $inpfile
echo ""                                                                        >> $inpfile
echo "RInitialThickness        = ${initialthickind}"                           >> $inpfile
echo "RFinalThickness          = ${finalthickind}"                             >> $inpfile
echo "RDeltaThickness          = ${deltathickind}"                             >> $inpfile
echo "IReflectOut              = ${reflectind}"                                >> $inpfile
echo ""                                                                        >> $inpfile 



#create hkl file for failed test  
if  [ ${casenumber} -eq 11 ]
#hkl file

touch $hklfile

echo "[1,0,0]"                                                                 >> $hklfile

cat $hklfile
fi

#create hkl file for reflection test
if  [ ${casenumber} -eq 2 ]
#hkl file

touch $hklfile

echo "[0,-4,-4]"                                                                >> $hklfile
echo "[3,1,1]"                                                                  >> $hklfile
echo "[2,0,0]"                                                                  >> $hklfile

cat $hklfile
fi

cat $inpfile
ls -al \${tmpdir}

fi #disregard only in case 14

# Make sure ranks are numbered appropriately
export MPICH_RANK_REORDER_METHOD=1 
export MPICH_PTL_MATCH_OFF=1
export MPICH_FAST_MEMCPY=1
export IPATH_NO_CPUAFFINITY=1
  
# Launch the parallel job, ${ranks} MPI ranks in ${taskspernode} tasks per node with ${ompthreads} threads
export OMP_NUM_THREADS=${ompthreads}

echo "--- starting the SIM run"
mpirun -n ${cores} \${tmpdir}/${binSIM} 

#echo "--- starting the DRAW run"
#mpirun -n ${cores} \${tmpdir}/${binDRAW} 

echo "--- starting the IMAGE post-processing run"
pwd
ls

echo "--- starting the image part"
# copy the result of the run in the detination for safekeeping
targetdir=\${basedir}/`basename ${job_file} .sh`

[ -d \${targetdir} ] || mkdir \${targetdir}

cp -vr * \${targetdir}

wait
#exit 0
EOD

chmod 755 ${job_dir}/${job_file}
#(cd ${job_dir} ; qsub -q serial ./${job_file})
#(cd ${job_dir} ; qsub -q parallel ./${job_file})
#(cd ${job_dir} ; qsub -q devel ./${job_file})
#(cd ${job_dir} ; qsub -q taskfarm ./${job_file})
(cd ${job_dir} ; qsub ./${job_file})

done



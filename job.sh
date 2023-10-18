#! bin/bash

# For scripts required see: https://gist.github.com/matteoferla/e0496d5766c12a0ae1738b943b41a536
# A few things don't work in CentOS 7 due to GNU lib C (glibc) 2.17, so it has to be run in a Singularity container

#: << USAGE #########
#export JOB_SCRIPT=/data/xchem-fragalysis/shared/singularity.sh;
#export APPTAINER_CONTAINER=/data/xchem-fragalysis/shared/singularity/rockyplus.sif; # CHANGE?
#export JOB_INNER_SCRIPT=/data/xchem-fragalysis/kfieseler/repos/elaborate/job.sh; # CHANGE TO THIS PATH
#condor_submit /data/xchem-fragalysis/shared/target_script.condor
#USAGE #########

#: << INSTALLATION #############
#export PIP_NO_USER=1; # gets ignored.
#export PYTHONUSERBASE=$CONDA_PREFIX; # double tap
#conda install -y -n base -c conda-forge openbabel plip;
#pip install -q fragmenstein
#pip install -q pyrosetta-help
#PYROSETTA_USERNAME=levinthal PYROSETTA_PASSWORD=paradox install_pyrosetta
#INSTALLATION ############

export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
source /etc/os-release;
echo "Running script ${0} as $USER in $HOST which runs $PRETTY_NAME"
# ---------------------------------------------------------------

source /data/xchem-fragalysis/kfieseler/.bashrc
conda activate fragmenstein

cd /data/xchem-fragalysis/kfieseler/retrievesynthesizable

pwd;
export TEMPLATE_DIR="/data/xchem-fragalysis/kfieseler/D68EV3CPROA/apo_desolv"
export HOME_DIR="/data/xchem-fragalysis/kfieseler/D68EV3CPROA/elabs/1_step"
export INPUT_SDF="/data/xchem-fragalysis/kfieseler/D68EV3CPROA/fragalysis/D68EV3CPROA_combined.sdf"

N_CORES=$(cat /proc/cpuinfo | grep processor | wc -l)
nice -19 python fragmenstein_batch.py \
-d HOME_DIR \
-t TEMPLATE_DIR \
-i INPUT_SDF \
-p "D68EV3CPROA-";

echo 'COMPLETE'

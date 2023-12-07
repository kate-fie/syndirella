#! bin/bash

# For scripts required see: https://gist.github.com/matteoferla/e0496d5766c12a0ae1738b943b41a536
# A few things don't work in CentOS 7 due to GNU lib C (glibc) 2.17, so it has to be run in a Singularity container

#: << USAGE #########
# source run_job.env
# condor_submit /data/xchem-fragalysis/shared/target_script.condor
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

cd /data/xchem-fragalysis/kfieseler/repos/elaborate

pwd;
export TEMPLATE_DIR="/data/xchem-fragalysis/kfieseler/A71EV2A/apo_desolv/";
export HOME_DIR="/data/xchem-fragalysis/kfieseler/A71EV2A/elabs/1_step_dec7_batched/batch_8";
export INPUT_SDF="/data/xchem-fragalysis/kfieseler/A71EV2A/clean_hits.sdf";
export N_CORES=$(cat /proc/cpuinfo | grep processor | wc -l)
export STEP="1_of_1"
export BATCH_RANGE="0-3"
export EXACT_TEMPLATE="A71EV2A-x0310_0A"

echo "Running fragmenstein_batch_A71EV2A.py"
echo "HOME_DIR: $HOME_DIR"
echo "TEMPLATE_DIR: $TEMPLATE_DIR"
echo "INPUT_SDF: $INPUT_SDF"
echo "N_CORES: $N_CORES"
echo "STEP: $STEP"
echo "BATCH_RANGE: $BATCH_RANGE"
echo "EXACT_TEMPLATE: $EXACT_TEMPLATE"

time nice -19 python fragmenstein_batch_A71EV2A.py \
-d $HOME_DIR \
-t $TEMPLATE_DIR \
-i $INPUT_SDF \
--step $STEP \
--n_cores $(($N_CORES - 1)) \
--cutoff \
--wictor \
--batch_range $BATCH_RANGE \
--exact_template $EXACT_TEMPLATE;

echo 'COMPLETE'


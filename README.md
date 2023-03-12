# bamboo-testanalysis
- Analysis uses [Bamboo RDataFrame](https://bamboo-hep.readthedocs.io/en/latest/index.html)

## Bamboo installation(1st time):
```bash
mkdir bamboodev
cd bamboodev

# make a virtualenv
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate

# clone and install bamboo
git clone -o upstream https://gitlab.cern.ch/cp3-cms/bamboo.git
pip install ./bamboo

# clone and install plotIt
git clone -o upstream https://github.com/cp3-llbb/plotIt.git
mkdir build-plotit
cd build-plotit
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
make -j2 install
cd -

# These last two cmd are needed everytime you upgrade your LCG working version!
# To use scalefactors and weights in the new CMS JSON format, the correctionlib package should be installed with
# you can ignore torch and sphinx pip errors !
pip install --no-binary=correctionlib correctionlib

# To use the calculators modules for jet and MET corrections and systematic variations
pip install git+https://gitlab.cern.ch/cp3-cms/CMSJMECalculators.git
```
- Let's make things more simpler, in your ``~/.bashrc`` you can add:
```bash
function cms_env() {
    module --force purge
    module load cp3
    module load cms/cmssw
    module load grid/grid_environment_sl7
    module load crab/crab3
    module load slurm/slurm_utils
}
alias bamboo_env="source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh"
alias bambooenv="source $HOME/bamboodev/bamboovenv/bin/activate"

function bamboo_start(){
    cms_env
    voms-proxy-init --voms cms --valid 168:00
    bamboo_env
    bambooenv
}
```
- And, in your ``~/.config/bamboorc`` add:
```
[batch]
backend = slurm

[slurm]
sbatch_qos = cp3
sbatch_partition = cp3
sbatch_additionalOptions = --licenses=cms_storage:3
sbatch_time = 6:59:00
sbatch_memPerCPU = 7000

[das]
sitename = T2_BE_UCL
storageroot = /storage/data/cms
checklocalfiles = True
xrootdredirector = xrootd-cms.infn.it
```
## Environment Setup (always):
- Every time you want to setup your bamboo enviroment, just do:
```bash
bamboo_start
```
## To test your bamboo instalation:
```bash
bambooRun -m /path/to/your/bambooclone/examples/nanozmumu.py:NanoZMuMu /path/to/your/bambooclone/examples/test1.yml -o test1
```
## To test your code, move to the ```python``` directory and run:
```bash
bambooRun -m controlPlotter.py ../config/analysis.yml -o ../test/myPlots --samples ../config/samples_template.yml --test
```

## Notes
- this code is based on [ttbbRun2Bamboo](https://gitlab.cern.ch/swertz/ttbbRun2Bamboo/)

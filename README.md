# Analysis/Tools

**Core Analysis Framework of the DESY CMS Higgs -> bb group**

See also the code documentation [Doxygen](https://www.desy.de/~walsh/docs/analysis-framework/doxygen/latest) page


:warning: This package is supposed to contain only general codes to be used for analysis.
Codes for specific analysis must be developed in dedicated packages/repositories, e.g.
for the MSSM Hbb analyses one can use for developments the package
[Analysis/MssmHbb](https://github.com/desy-cms/analysis-mssmhbb),
which is currently under construction.


- [Installation](#installation)
- [Calibrations](#calibrations)
- [Ntuples](#ntuples)
- [Example](#example)
- [NAF Submission](#naf-submission)
- [Example Detailed Description](#example-detailed-description)
- [Luminosity calculations on the NAF](#luminosity-calculations-on-the-naf)



## Installation

The codes as well as the ntuples are independent of CMSSW. However, in order to compile it uses `scram`.
So the latest version in any architecture should be fine.

```bash
cmsrel CMSSW_X_Y_Z
cd CMSSW_X_Y_Z/src
cmsenv
git clone https://github.com/desy-cms/analysis-tools.git Analysis/Tools
scram b -j4 USER_CXXFLAGS="-Wno-misleading-indentation"
```

:zzz: The `USER_CXXFLAGS="-Wno-misleading-indentation"` prevents a large number of warnings
from misleading identation in modules of the boost library. User can also define the
environment variable below in `.bashrc` or every time after the command `cmsenv`
```bash
export USER_CXXFLAGS="-Wno-misleading-indentation"
```

## Calibrations

Scale factors, efficiencies etc can be obtained from the [analysis-calibrations](https://github.com/desy-cms/analysis-calibrations) repository (see also the README.md in each run period directory). It is recommended to install it in the Analysis/Tools/data directory, for the framework loads the calibrations files from that path

```bash
cd $CMSSW_BASE/src
git clone https://github.com/desy-cms/analysis-calibrations.git Analysis/Tools/data/calibrations
```

## Ntuples

The lists of ntuples files can be obtained from the [analysis-ntuples](https://github.com/desy-cms/analysis-ntuples.git) repository (see also the README.md in each run period directory). The repository can be installed in a directory of your convenience, e.g.

```bash
cd $CMSSW_BASE/src
git clone https://github.com/desy-cms/analysis-ntuples.git Analysis/Tools/data/ntuples
```

The ntuples path name is the global one for general grid access. If you are running your jobs on the NAF, it should
be faster to use the `/pnfs/...` path name. To make the conversion simply issue the following command once for each
installation of the analysis-ntuples as described above

```bash
ntuples_at_pnfs.sh
```

The script will modify all the files in `$CMSSW_BASE/src/Analysis/Tools/data/ntuples`.


## Example

A simple example macro can be found in [Analysis/Tools/bin/AnalyserSimpleExample.cc](bin/AnalyserSimpleExample.cc).
The example is a simple offline and online dijet selection (+ muon selection in the semileptonic case)
using signal MC samples and triggers from 2017 and 2018 periods. The macro uses a configuration file as
input. Configuration files are avaialble for both the all-hadronic and semileptonic cases for both 2017
and 2018 periods:
* [analyser_example_allhad_2017.cfg](test/analyser_example_allhad_2017.cfg)
* [analyser_example_allhad_2018.cfg](test/analyser_example_allhad_2018.cfg)
* [analyser_example_semilep_2017.cfg](test/analyser_example_semilep_2017.cfg)
* [analyser_example_semilep_2018.cfg](test/analyser_example_semilep_2018.cfg)

To execute an example
```bash
cd Analysis/Tools/test
AnalyserSimpleExample -c analyser_example_semilep_2018.cfg
```

## NAF Submission

A python script to submit to NAF condor queue, `naf_submit.py`, is available. 

**N.B.: So far the script does not make a single submission of multiple jobs. So be careful not to make too many submissions.**

```
usage: naf_submit.py [-h] [--exe EXE] [--config CONFIG] [--ntuples NTUPLES] [--nfiles NFILES] [--json JSON] [--label LABEL] [--events EVENTS_MAX] [--test NJOBS] [--dir DIR] [--status] [--resubmit] [--expert]

Prepare, submit and check jobs to NAF HTCondor batch system

optional arguments:
      -h, --help                     show this help message and exit

submission:
      prepare and submit jobs

      --exe EXE, -e EXE              Executable (REQUIRED)
      --config CONFIG, -c CONFIG     Configuration file (REQUIRED)
      --ntuples NTUPLES, -n NTUPLES  List of ntuples file
      --nfiles NFILES, -x NFILES     Number of ntuple files per job
      --json JSON, -j JSON           JSON file with certified data
      --label LABEL, -l LABEL        user label for the submission
      --events EVENTS_MAX            override eventsMax in the config file (default = -1)
      --test NJOBS                   *** expert only ***:produce njobs, no automatic submission

status:
      show and modify status

      --dir DIR                      an existing condor directory (REQUIRED)
      --status                       -> returns the status of the jobs in --dir
      --resubmit                     -> resubmits aborted and finished-with-error jobs in --dir
      --expert                       -> *** expert mode ***
```

**If you provide a configuration file with the NTUPLES and JSON parameters, you do not need to parse them, the script will read out that information from the configuration file.**

Using the example above to be submitted to the naf using [HTCondor](https://confluence.desy.de/display/ITPublic/HTCondor%3A+Job+Submission)
```bash
naf_submit.py --exe AnalyserSimpleExample --config analyser_example_semilep_2018.cfg -nfiles 2
```

After the jobs were submitted there will be a directory called `Condor_AnalyserSimpleExample_analyser_example_semilep_2018` containing several subdirectories called `job_xxxx`. In each of this subdirectory there will be several files. Files to take note will be:
* `<output>.root` - the output file containing histograms etc of your analysis
* `finished.txt` - a file that indicates whether the executable ran until the end (only if you are using the framework)
* `job.submit` - HTCondor submit configuration file. If a job is not finished you can resubmit that job manually by issuing the command
```bash
condor_submit job.submit
```

Once the directory `Condor_AnalyserSimpleExample_analyser_example_semilep_2018`
is created, one can use the status options, e.g.

```bash
naf_submit.py --dir Condor_AnalyserSimpleExample_analyser_example_semilep_2018 --status
```

will display a table like this


```
 
                          ***  STATUS OF JOBS  ***

   Condor_AnalyserSimpleExample_analyser_example_semilep_2018
 
  ----------------------------------------------------------------------------------------------------------
     job        finished       running       submitted       aborted       error       condor_id (latest)
  ----------------------------------------------------------------------------------------------------------
   job_0000        v                                                         v             15962854.0
   job_0001        v                                                         v             15962854.1
   job_0002                       v                                                        15962854.2
   job_0003                                       v                                        15968568.0
   job_0004                                                      X                         15962854.4
   job_0005        v                                                         X             15962854.5
  ----------------------------------------------------------------------------------------------------------

  N.B.: Good *finished* jobs will no longer appear in future "--status" calls

```

where `v` will appear as a green tick mark and `x` as a red X.


If jobs were aborted or finished with error, you can resubmit all of those jobs with the command 

```bash
naf_submit.py --dir Condor_AnalyserSimpleExample_analyser_example_semilep_2018 --resubmit
```

which, in the example above, will resubmit the jobs `job_0004` and `job_0005`.


## Example Detailed Description

### Creating a macro

#### Main

In the macro the first thing done is the main `Analyser` instantiation

```cpp
#include "Analysis/Tools/interface/Analyser.h"
using namespace std;
using namespace analysis;
using namespace analysis::tools;

int main(int argc, char ** argv)
{
   TH1::SetDefaultSumw2();
   Analyser analyser(argc,argv);
...   
```

The `[Info]` block in the `.cfg` file
```
[Info]
process = MssmHbb
eventsMax = -1
ntuplesList = rootFileList.txt
isMC = true
output = histograms.root
```
passes general information to initialise the `Analyser` in the macro. The `process` parameter is the name of the parent `TDirectoryFile` in the ntuple. To find it
```bash
root -l ntuple.root
# in the root prompt type
.ls
# it will show
KEY: TDirectoryFile	MssmHbb;1	MssmHbb
```
Where `ntuple.root` is any input ntuple file, such as the ones in your `roorFileList.txt` or in the [analysis-ntuples](https://github.com/desy-cms/analysis-ntuples).

Hopefully the other names of the config parameters are self-explanatory.

In **data** one must specify the json file with certified lumis, e.g.
```
[Info]
...
isMC = false
json = certified_json.txt
```

#### Event loop

The event loop must start like this
```cpp
for ( int i = 0 ; i < analyser.nEvents() ; ++i )
{
   if ( ! analyser.event(i) )   continue;
   // do actions and selections
}
```
where the `analyser.event(i)` reads the event from the ntuple and performs some actions, such as applying generator event weights in MC and JSON certified lumis selection in data.

#### Selections

The `Analyser` has several predefined selection methods that reads parameters from the configuration file and apply to the event. The selections must be within the event loop.

##### Jets

For example, if the analysis involves jets, one must define which collection to be used, the minimum number of jets, the jet identification, the jet pT etc. In the configuration there is a block for jets with the relevant parameters
```
[Jets]
jets = updatedPatJets
nMin = 2
id = tight
puId = loose
ptMin = 60
ptMin = 50
etaMax = 2.2
etaMax = 2.2
extendedFlavour = true
```
:warning: Parameters such as ptMin or etaMax are vectors, so the order in which they are put in the configuration makes a difference, so the first entry corresponds to the leading object, the second entry to the second leading object and so on.

In the macro, the selections are performed within the event loop calling the methods, which automatically reads the parameters from the configuration, e.g.
```cpp
      // jet identification selection
      if ( ! analyser.selectionJetId()          )   continue;
      if ( ! analyser.selectionJetPileupId()    )   continue;
      if ( ! analyser.selectionNJets()          )   continue;
```
This will prepare a list of jets containing only jets passing the required identification criteria and with a certain number of jets defined in the configuration.

```cpp
      //  1st and 2nd jet kinematic selection, pt and eta
      if ( ! analyser.selectionJet(1)          )   continue;
      if ( ! analyser.selectionJet(2)          )   continue;

```
In the `Analyser::selectionJet` method the argument is the jet rank, i.e., `1` refers to the leading jets, `2` refers to the second leading jet. the method will select the jet according to its pt and eta criteria defined in the configuration.

##### b-tag

For the b-tagging there is a block in the configuration file, where one defines the algorithm, working points etc.

```
[BTag]
nMin  = 2
wp = medium
wp = medium
algorithm = deepflavour
loose  = 0.0494
medium = 0.2770
tight  = 0.7264
```
With this configuration, one intends to select events with at least the two leading jets tagged with medium working point using the deepflavour algorithm. The thresholds of the working points myst be specified.

To perform the selection in the event loop:

```cpp
      if ( ! analyser.selectionBJet(1)         )   continue;
      if ( ! analyser.selectionBJet(2)         )   continue;

```
where the argument of the `Analyser::selectionBJet` is the jet rank as in the jet kinematic selection above.

##### Muons

Muon selection has its own configuration block and the procedure works in a similar way as jets

```
[Muons]
muons = slimmedMuons
nMin = 1
id = tight
ptMin = 13.
etaMax = 2.2
```

The selection code can be for example

```cpp
      // muon identification selection
      if ( ! analyser.selectionMuonId()         )   continue;
      if ( ! analyser.selectionNMuons()         )   continue;
      // muon kinematic selection
      if ( ! analyser.selectionMuons()          )   continue;
```
:warning: The  method `Analyser::selectionMuons` for muon kinematic selection differs from the one used for jet kinematic selection, for it makes a list of muons passing the required kinematic selection and the event is good if there is at list the minimum number of jets required. This method is useful for the muon-jet association in analyses requiring jets containing a muon.

If one wants to select muons depending on their ranking, then one can use, like in the jets case,
```cpp
      // leading muon kinematic selection
      if ( ! analyser.selectionMuon(1)          )   continue;
```


#### Histograms

Histograms can be predefined in the `Analyser`, e.g. for jets one can use the method `jetHistograms`, which receives the arguments `number of jets`, which must be at most the minimum number of jets required, and a name of a directory. One can create as many directories as needed and fill the histograms at different steps of the analysis workflow, e.g. before the event loop:
```cpp
   analyser.jetHistograms(2,"initial");
   analyser.jetHistograms(2,"final");
```
then within the event loop:
```cpp
...
for ( int i = 0 ; i < analyser.nEvents() ; ++i )
{
   if ( ! analyser.event(i) )   continue;
   // do something
   analyser.fillJetHistograms("initial");
   // do something else
   analyser.fillJetHistograms("final");
...
```
## Luminosity calculations on the NAF

In order to use `brilcalc` on the NAF, based on the [BRIL Work Suite documentation](https://cms-service-lumi.web.cern.ch/cms-service-lumi/) and using python2 (python3 did not work(?)), do the following:

```bash
wget https://cern.ch/cmslumisw/installers/linux-64/Brilconda-1.1.7-Linux-x86_64.sh
bash Brilconda-1.1.7-Linux-x86_64.sh -b -p /nfs/dust/cms/user/<your_naf_username>/brilconda
```
Replace `<your_naf_username>` by your NAF username, also below.

Add the following to your shell environment, e.g. bash
```bash
export PATH=/nfs/dust/cms/user/<your_naf_username>/brilconda/bin:$PATH
```
Open a new terminal and test the command
```bash
brilcalc -h
```
If you see the help for the command, then the installation is ok and you can follow the instructions for luminosity calculations in the [BRIL Work Suite documentation](https://cms-service-lumi.web.cern.ch/cms-service-lumi/).

### Scripts using brilcalc

To facilitate the life of the users, some scripts are available, e.g. to calculate the luminostities of HLT paths given a certified JSON file.

:warning: While analyses are not fully migrated to CMSSW releases using `python3`, namely Run2 analyses, the scripts will reside in a different branch of `analysis-tools` that may not be up-to-date with the `master`. Therefore, use it only for the luminosity calculations. Once Run2 analyses are finished, the scripts should be part of the `master` branch.

In any case, the scripts should work in any machine, lxplus or naf, as long as brilcalc is installed.

Usually, running brilcalc takes a while. The scripts use multithreading to help speeding up the computation.

The installation, e.g. on the NAF `el8`, is done like this:
```bash
export SCRAM_ARCH=el8_amd64_gcc11
cmsrel CMSSW_13_2_4
cd CMSSW_13_2_4/src
cmsenv

git clone https://github.com/robervalwalsh/analysis-tools.git Analysis/Tools
cd Analysis/Tools
git checkout brilcalc
git clone https://github.com/desy-cms/analysis-calibrations.git data/calibrations
cd $CMSSW_BASE/src
scram b -j4
```

The scripts are available in `Analysis/Tools/scripts`

#### HLT paths luminosities

The script that calculates the luminosity for HLT paths given a certified JSON is `hlt_lumi.py`
```bash
$ hlt_lumi.py --help
usage: hlt_lumi.py [-h] [--json JSON] [--triggers TRIGGERS] [--normtag NORMTAG] [--unit UNIT] [--threads THREADS] [--decimals DECIMALS]
                   [--output OUTPUT]

Obtain HLT paths luminosities (brilcal)

optional arguments:
  -h, --help           show this help message and exit
  --json JSON          Path to the Golden JSON file
  --triggers TRIGGERS  List of triggers (comma-separated or in a text file)
  --normtag NORMTAG    Path to the normtag file (default: /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json)
  --unit UNIT          Unit default: /pb)
  --threads THREADS    Number of threads (default: 10)
  --decimals DECIMALS  Number of decimals in results (default: 4)
  --output OUTPUT      Output file

```
Notice that the list of triggers can be given in the command line, separated by commas, from a text file, where the list can be given one trigger per line, or comma separated. The obligatory parameters are the `--json` and `--triggers` parameters
For example, the lumit 

```bash
cd $CMSSW_BASE/src/Analysis/Tools/test
hlt_lumi.py \
--json=../data/calibrations/2017/certified/Cert_Run2017CDEF_13TeV_UL2017_Collisions17_GoldenJSON.txt \
--triggers=HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v*,HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33_v*
```
where the output is a markdown table that can be redirected to an `--output` file.
```
| HLT Path | Recorded Luminosity [/pb] |
| --- | ---: |
| HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33 | 36263.6748 |
| HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagCSV_p33 | 36674.5111 |
```

#### Active L1 seeds

In the situation where the L1 seed is an OR of triggers, one may obtain the certified JSON file containing the lumi sections where a given L1 trigger is active/inactive. The script scan all the runs in the JSON file and find the lumi sections where the L1 trigger has a prescale different from zero, i.e., the trigger is active. Prescale zero means trigger is inactive. The list of active lumi sections is compared to the JSON file and all LS where the trigger is active are kept in the final JSON.
```bash
$ active_l1seed_json.py --help
usage: active_l1seed_json.py [-h] [--json JSON] [--hlt HLT] [--l1 L1] [--threads THREADS]

Obtain certified LS when the L1 seed is active/inactive (brilcal)

optional arguments:
  -h, --help         show this help message and exit
  --json JSON        Path to the Golden JSON file
  --hlt HLT          HLT Path
  --l1 L1            L1 Seed
  --threads THREADS  Number of threads (default: 20)
```

For example, the MSSM Hbb full hadronic trigger in 2017
```bash
active_l1seed_json.py \
--json ../data/calibrations/2017/certified/Cert_Run2017CDEF_13TeV_UL2017_Collisions17_GoldenJSON.txt \
--hlt HLT_DoublePFJets100MaxDeta1p6_DoubleCaloBTagCSV_p33_v* \
--l1 L1_DoubleJet100er2p3_dEta_Max1p6
```
The output are two files, with the certified LS where the L1 trigger was active(inactive), whose names are built from the original JSON file name:
- `Cert_Run2017CDEF_13TeV_UL2017_Collisions17_GoldenJSON_L1Active.txt`
- `Cert_Run2017CDEF_13TeV_UL2017_Collisions17_GoldenJSON_L1Inactive.txt`


#!/usr/bin/env python
import os
import glob
import re
from argparse import ArgumentParser
from argparse import HelpFormatter
from shutil import copyfile, rmtree, move
from time import sleep

# --- functions ---
def getConfigParameter( config, parameter ):
   p = None
   exist = False
   with open(config) as f:
      for line in f:
         line = line.replace(" ","").strip()
         if len(line) == 0:
            continue
         if line[0] != '#' and  line.split("=")[0] == parameter:
            par = line.split("=")[1]
            exist = True
            p = [parameter,par]
            if parameter=="ntuplesList" and "tools:" in p[1]:
               ntp_path = os.getenv('CMSSW_BASE')
               ntp_path += "/src/Analysis/Tools/data/ntuples/"
               p[1] = p[1].replace("tools:",ntp_path)
            break
            
   return p

def createConfigParameter( config, parameter ):
   exist = False
   with open(config,"r") as f:
      lines = f.readlines()
      
   for line in lines:
      line = line.replace(" ","").strip()
      if len(line) == 0:
         continue
      if line[0] != '#' and  line.split("=")[0] == parameter:
         exist = True
         break
         
   if not exist:
      with open(config, "w") as f:
         f.write(parameter+" = \n")
         for line in lines:
            f.write(line)
   

def replaceConfigParameter( config, parameter, newpar ):
   par = None
   with open(config, "r") as f:
      lines = f.readlines()
   with open(config, "w") as f:
      for line in lines:
         l = line.replace(" ","").strip()
         if len(l) < 1:
            f.write(line)
            continue
         if l[0] != '#' and l.split('=')[0] == parameter:
            if l.split('=')[1] == "" :
               par    = "="
               newpar = " = " + newpar
            else:
               par = l.split('=')[1]
            f.write(re.sub(par, newpar, line))
         else:
            f.write(line)

def removeConfigParameter( config, parameter ):
   with open(config, "r") as f:
      lines = f.readlines()
   with open(config, "w") as f:
      for line in lines:
         l = line.replace(" ","").strip()
         if len(l) < 1:
            f.write(line)
            continue
         if l[0] != '#' and l.split('=')[0] == parameter:
            continue;
         f.write(line)

   
def basenameConfigParameter( config, name ):
   with open(config, "r") as f:
      lines = f.readlines()
   with open(config, "w") as f:
      for line in lines:
         f.write(re.sub(name, os.path.basename(name), line))



# --- main code ---

# parsing arguments
parser = ArgumentParser(prog='naf_submit.py', formatter_class=lambda prog: HelpFormatter(prog,indent_increment=6,max_help_position=80,width=280), description='Prepare and submit jobs to NAF HTCondor batch system')
parser.add_argument("-e", "--exe", dest="exe", help="Executable")
parser.add_argument("-n", "--ntuples", dest="ntuples", help="List of ntuples file")
parser.add_argument("-x", "--nfiles", dest="nfiles", type=int, default=1, help="Number of ntuple files per job")
parser.add_argument("-c", "--config", dest="config", help="Configuration file")
parser.add_argument("-j", "--json", dest="json", help="JSON file with certified data")
parser.add_argument("-l", "--label", dest="label", help="user label for the submission")
parser.add_argument("--events", dest="events_max", help="override eventsMax in the config file")
parser.add_argument("--test", dest="njobs", help="produce njobs, no automatic submission")
args = parser.parse_args()
if not args.exe:
   print "nothing to be done" 
   quit()
   
ntuples = args.ntuples
json = args.json
config = args.config
events_max = args.events_max
test = args.njobs

configMC = getConfigParameter( config, "isMC" )
isMC = configMC[1] == 'true'

if test:
   print('TEST MODE:', test, 'jobs')

configNtuples = None
# get parameter from configuration 
if config:
   if not os.path.isfile(config):
      print "Configuration file does not exist"
      quit()
   configNtuples = getConfigParameter( config, "ntuplesList" )
   if not ntuples:
      if configNtuples:
         ntuples = configNtuples[1]
      if not ntuples:
         print "*error* You must define the parameter ntuplesList in your configuration."
         quit()
   if not isMC:
      configJson    = getConfigParameter( config, "json" )
      if not json:
         if configJson:
            json = configJson[1]
   else:
      json = None
      
# checking if require files exist
if ntuples:
   if not os.path.isfile(ntuples):      
      print "Ntuples list file does not exist"
      quit()
if json:
   if (not 'tools:' in json) and (not os.path.isfile(json)):      
      print "Json  file does not exist"
      quit()
      
# directory where the jobs will be stored
maindir = "Condor_"+os.path.basename(args.exe)
if config:
   maindir = maindir+"_"+ os.path.splitext(os.path.basename(config))[0]
if args.label:
   maindir += '_'+args.label
cwd = os.getcwd()
if os.path.exists(cwd+"/"+maindir):
   print maindir + "already exists. Rename or remove it and then resubmit"
   quit()
os.mkdir(maindir)

# splitting the file list
if ntuples:
   pid = os.getpid()
   tmpdir = ".tmp_" + str(pid)
   os.mkdir(tmpdir)
   copyfile(ntuples, tmpdir+"/"+os.path.basename(ntuples))
   if config:
      copyfile(config, tmpdir+"/"+os.path.basename(config))
   os.chdir(tmpdir)
   if config:
      # replace json and ntuples in the local exe config by their basenames
      if json:
         createConfigParameter(os.path.basename(config),'json')
         replaceConfigParameter(os.path.basename(config), 'json', os.path.basename(json))
      else:
         removeConfigParameter(os.path.basename(config),'json')
      # ntuples list
      createConfigParameter(os.path.basename(config),'ntuplesList')
      replaceConfigParameter(os.path.basename(config), 'ntuplesList', os.path.basename(ntuples))
      if events_max:
         replaceConfigParameter(os.path.basename(config), 'eventsMax', events_max)

   
   splitcmd = "split.csh" + " " + str(args.nfiles) + " " + os.path.basename(ntuples)
   os.system(splitcmd)
   os.remove(os.path.basename(ntuples))  # not needed anymore after the splitting
   files = glob.glob('.*_x????.txt')
   files.sort()
   os.chdir(cwd)

   # loop over the splitted files, each will correspond to a job on the NAF
   
   for i,f in enumerate(files):
      if test:
         if i >= int(test) and int(test)>=0:
            os.chdir(cwd)
            break
      jobnum = os.path.splitext(f)[0][-4:]
      jobid = "job_"+jobnum
      exedir = maindir+"/"+jobid
      os.mkdir(exedir)
      # moving stuff to the proper directories
      move(tmpdir+"/"+f,exedir+"/"+os.path.basename(ntuples))
      if json:
         if not 'tools:' in json:
            copyfile(json, exedir+"/"+os.path.basename(json))
      if config:
         copyfile(tmpdir+"/"+os.path.basename(config),exedir+"/"+os.path.basename(config))      
         condorcmd = "condor_scripts.csh" + " " + jobid + " " + args.exe + " " + os.path.basename(config)
      else:
         condorcmd = "condor_scripts.csh" + " " + jobid + " " + args.exe
      # make the submissions
      os.chdir(exedir)
      jobf = open('./seed.txt', 'w+')
      print >> jobf, int(jobnum)+1
      jobf.close()
      print "Creating ",jobid,"..."
      os.system(condorcmd)
      if not test:
         os.system('condor_submit job.submit')
      sleep(0.2)
      # back to original directory
      os.chdir(cwd)

else:
   exedir = maindir+"/job_0000"
   os.mkdir(exedir)
   if os.path.isfile(args.exe):
      copyfile(args.exe, exedir+"/"+os.path.basename(args.exe))
   os.chdir(exedir)
   jobf = open('./seed.txt', 'w+')
   print >> jobf, 1
   jobf.close()
   condorcmd = "condor_scripts.csh job_0000" + " " + os.path.basename(args.exe)
   os.system(condorcmd)
   if not test:
      os.system('condor_submit job.submit')
   os.chdir(cwd)
         
# remove the temporary directory
os.chdir(cwd)
if ntuples: 
   os.remove(tmpdir+"/"+os.path.basename(config))
   rmtree(tmpdir)



#!/usr/bin/env python3
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


def createCondorDirsSubmit(maindir):
   libpath = os.environ['LD_LIBRARY_PATH']
   condor_submit = """
Requirements = ( OpSysAndVer == "CentOS7" )
getenv      = True
executable  = $(jobdir)/job.sh
log         = $(jobdir)/xjob_$(Cluster)_$(Process).log
output      = $(jobdir)/xjob_$(Cluster)_$(Process).out
error       = $(jobdir)/xjob_$(Cluster)_$(Process).err
environment = "LD_LIBRARY_PATH_STORED=$ENV(LD_LIBRARY_PATH)"

queue jobdir matching dirs job_*
   """
   with open(f'{maindir}/jobs.submit', 'w') as condor_file:
      print(condor_submit,file=condor_file)

def createCondorSubmit(jobdir):
   libpath = os.environ['LD_LIBRARY_PATH']
   condor_submit = """
Requirements = ( OpSysAndVer == "CentOS7" )
getenv      = True
executable  = job.sh
log         = xjob_$(Cluster)_$(Process).log
output      = xjob_$(Cluster)_$(Process).out
error       = xjob_$(Cluster)_$(Process).err
environment = "LD_LIBRARY_PATH_STORED=$ENV(LD_LIBRARY_PATH)"

queue
   """
   with open(f'{jobdir}/job.submit', 'w') as condor_file:
      print(condor_submit,file=condor_file)

def createJobScript(exedir,jobid, exe, cfg):
   js = """
if [[ ! -e "{}" ]]; then
   cd {}
fi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_STORED
{} -c {}
   """.format(cfg,jobid,exe,cfg)
   js_name  = f'{exedir}/job.sh'
   with open(js_name, 'w') as job_script:
      print(js,file=job_script)
   os.chmod(js_name, 0o744)


def submission():
   
   ntuples = args.ntuples
   json = args.json
   config = args.config
   events_max = -1
   if args.events_max:
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
         print ("Configuration file does not exist")
         quit()
      configNtuples = getConfigParameter( config, "ntuplesList" )
      if not ntuples:
         if configNtuples:
            ntuples = configNtuples[1]
         if not ntuples:
            print ("*error* You must define the parameter ntuplesList in your configuration.")
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
         print ("Ntuples list file does not exist")
         quit()
   if json:
      if (not 'tools:' in json) and (not os.path.isfile(json)):      
         print ("Json  file does not exist")
         quit()
         
   # directory where the jobs will be stored
   maindir = "Condor_"+os.path.basename(args.exe)
   if config:
      maindir = maindir+"_"+ os.path.splitext(os.path.basename(config))[0]
   if args.label:
      maindir += '_'+args.label
   cwd = os.getcwd()
   if os.path.exists(cwd+"/"+maindir):
      print (maindir,"already exists. Rename or remove it and then resubmit")
      quit()
   os.mkdir(maindir)
   os.mkdir(maindir+'/finished_jobs')
   
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
            createCondorSubmit(cwd+'/'+exedir)
            createJobScript(exedir,jobid,args.exe,os.path.basename(config))
         else:
            condorcmd = "condor_scripts.csh" + " " + jobid + " " + args.exe
         # make the submissions
   #      os.chdir(exedir)
         jobf = open(f'{exedir}/seed.txt', 'w+')
         jobf.write(str(int(jobnum)+1))
   #      print >> jobf, int(jobnum)+1
         jobf.close()
         print ("Creating ",jobid,"...")
   #      os.system(condorcmd)
   #      if not test:
   #         os.chdir(exedir)
   #         os.system('condor_submit job.submit')
   #         sleep(0.2)
         # back to original directory
   #      os.chdir(cwd)
         
      createCondorDirsSubmit(cwd+'/'+maindir)
      if not test:
         os.chdir(cwd+'/'+maindir)
         os.system('condor_submit jobs.submit')
   
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


# given the directory of a job_xxx return the most recent log file
def get_job_log(jobdir):
   files = os.listdir(jobdir)
   log = [ jobdir+'/'+x for x in files if x.endswith('.log') ]
   if len(log) > 1:
      return max(log, key=os.path.getctime)
   elif len(log) == 1:
      return log[0]
   else:
      return None

# given the directory of a job_xxx check if there is anything in the most recent err file
def job_err(jobdir):
   files = os.listdir(jobdir)
   errs = [ jobdir+'/'+x for x in files if x.endswith('.err') ]
   if len(errs) > 1:
      err = max(errs, key=os.path.getctime)
   elif len(errs) == 1:
      err = errs[0]
   else:
      # if err file does not exist, there is no error
      return False
      
   return (os.path.getsize(err) > 0)
   


# given the directory of a job_xxxx check the status in the most recent log file
def get_job_status(jobdir):
   job_status  = {}
   job_status['submission'] = False
   job_status['termination'] = False
   job_status['execution'] = False
   job_status['abortion'] = False
   job_status['error'] = True
   job_status['jobid'] = ' '
   log = get_job_log(jobdir)
   if not log:
      return job_status
   logname = os.path.basename(log)
   job_status['jobid'] = logname.split('.')[0].split('_')[1] + '.' + logname.split('.')[0].split('_')[2]
   f = open(log,'r')
   for l in f:
      if re.search('Job submitted from host',l):
         job_status['submission'] = True
      if re.search('Job executing on host',l):
         job_status['execution'] = True
      if re.search('Job terminated.',l):
         job_status['termination'] = True         
      if re.search('Job was evicted',l) or re.search('Job was aborted',l):
         job_status['abortion'] = True
      if re.search('return value 0',l):
         job_status['error'] = False
      
   return job_status  
         
       
   

def status(submission_dir, failed_only=False):
   failed_jobs = 0
#   submission_dir = args.status
   finished_dir = submission_dir+'/finished_jobs'
   jobs_dir = os.listdir(submission_dir)
   jobs_dir = [ x for x in jobs_dir if 'job_' in x ]
   jobs_dir.sort()
   print(' ')
   header = '                          ***  Status of jobs  ***'
   if failed_only:
      header = '                          ***  FAILED JOBS  ***'
   print(header)
   print('\n   '+submission_dir)
   print(' ')
   dash= '  ----------------------------------------------------------------------------------------------------------'
   print(dash)
   print('     job        finished       running       submitted       aborted       error       condor_id (latest)')
   print(dash)
   for jj in jobs_dir:
      j = submission_dir+'/'+jj
      js = get_job_status(j)
      je = job_err(j)
      finished = ' '
      running = ' '
      submitted = ' '
      aborted = ' '
      error = ' '  # check
      jobid = js['jobid']
      if js['termination']:
         finished = '\u2705'  # check
         error = '\u2705 '  # check
#         if js['error']:
#            finished = '\u2705'
         if js['error'] and not js['abortion']:
            error = '\u274C' # red x
         if not js['abortion'] and not js['error']:
            # search for finished.txt before moving job dir, guarantee that all everything finished
            if os.path.isfile(j+'/finished.txt'):
               sleep(1)
               move(j,finished_dir)
            else:
               finished = '\u26A0\uFE0F ' # warning 
#      elif js['execution'] and not js['abortion'] and not js['error']:
      elif js['execution'] and not js['abortion']:
         running = '\u2705'
      elif js['submission'] and not js['abortion']:
         submitted = '\u2705'
      if js['abortion']:
         aborted = '\u274C '  # red cross
         
      if not ( js['submission'] or js['termination'] or js['execution'] or js['abortion'] ):
         submitted = '\u26A0\uFE0F ' # warning
      if js['abortion'] or ( js['termination'] and js['error']):
         failed_jobs +=1
      if failed_only and (not (js['abortion'] or ( js['termination'] and js['error']))):
            continue
      print('   '+jj+'        '+finished+'              '+running+'              '+submitted+'              '+aborted+'           '+error+'            '+jobid)
   if len(jobs_dir) == 0:
      print('  No jobs to be checked!')
   print(dash)
   if not failed_only:
      print('\n  N.B.: Good *finished* jobs will no longer appear in future "--status" calls')
      
   return failed_jobs
   
# resubmit jobs according to their status
def resubmit(submission_dir):
   failed_jobs = status(submission_dir,True)
   if failed_jobs == 0:
      return
   jobs_dir = os.listdir(submission_dir)
   jobs_dir = [ x for x in jobs_dir if 'job_' in x ]
   jobs_dir.sort()
   cwd = os.getcwd()
   print(' ')
   print('   ***  Resubmit failed jobs  ***')
   print(dash)
   confirmed = confirm()
   if not confirmed:
      print(dash)
      print(' ')
      return
   for jj in jobs_dir:
      j = submission_dir+'/'+jj
      js = get_job_status(j)
      if js['abortion'] or ( js['termination'] and js['error'] ):
         os.chdir(j)
         os.system('condor_submit job.submit')
         os.chdir(cwd)
   print(dash)
   print(' ')
   os.chdir(cwd)
   sleep(2)
   status(submission_dir)
   
   
   
# resubmit jobs expert mode
def resubmit_expert(submission_dir):
   cwd = os.getcwd()
#   submission_dir = args.resubmit_expert
   os.chdir(submission_dir)
   print(' ')
   print('                       ***  Resubmit jobs (EXPERT)  ***')
   print('\n  '+submission_dir)
   print(' ')
   print(dash)
   confirmed = confirm()
   if not confirmed:
      print(dash)
      print(' ')
      return
   os.system('condor_submit jobs.submit')
   print(dash)
   print(' ')
   os.chdir(cwd)
   sleep(2)
   status(submission_dir)
   

def confirm():
    answer = " "
    while answer not in ["y", "n", ""]:
        answer = input("   Please confirm this action (Y/[N])? ").lower()
    return (answer == "y" or answer != "")


# --- main code ---

# parsing arguments
parser = ArgumentParser(prog='naf_submit.py', formatter_class=lambda prog: HelpFormatter(prog,indent_increment=6,max_help_position=80,width=280), description='Prepare and submit jobs to NAF HTCondor batch system')
parser.add_argument("-e", "--exe", dest="exe", help="Executable")
parser.add_argument("-n", "--ntuples", dest="ntuples", help="List of ntuples file")
parser.add_argument("-x", "--nfiles", dest="nfiles", type=int, default=1, help="Number of ntuple files per job")
parser.add_argument("-c", "--config", dest="config", help="Configuration file")
parser.add_argument("-j", "--json", dest="json", help="JSON file with certified data")
parser.add_argument("-l", "--label", dest="label", help="user label for the submission")
parser.add_argument("--events", dest="events_max", default="-1", help="override eventsMax in the config file (default = -1)")
parser.add_argument("--test", dest="njobs", help="produce njobs, no automatic submission")
parser.add_argument("--status", dest="status", help="status of a given submission")
parser.add_argument("--resubmit", dest="resubmit", help="resubmit aborted and finished-with-error jobs")
parser.add_argument("--resubmit_expert", dest="resubmit_expert", help="*** expert only ***: resubmit available jobs")
args = parser.parse_args()

dash = '  ----------------------------------------------------------------------------------------------------------'
   
if args.status:
   status(args.status)
elif args.resubmit:
   resubmit(args.resubmit)
elif args.resubmit_expert:
   resubmit_expert(args.resubmit_expert)
else:
   if not args.exe and not args.config:
      print("nothing to be done")
      quit()
   submission()
   

#!/bin/csh -f

########################################################################
#                                                                      #
# Find root histograms inside some directories and add them with hadd  #
#                                                                      #
# This script is to be used mostly on the output of the naf_submit.py  #
#                                                                      #
# R. Walsh, Jan 2019                                                   #
#                                                                      #
########################################################################

if ( $#argv < 1 ) then
   echo "Usage: hadd2.csh directory_keyword [expression]"
   exit
endif

set dirk     = $1
set basedir  = `pwd`
set finished = "finished!"
set expression = ""
if ( $#argv == 2 ) then
   set expression = $2
endif   

set list_dirs = `ls -1d */ | grep $dirk`

foreach dir ($list_dirs)
#   cd $dir
   @ njobs = `ls -1 $dir/ | grep "job_" | wc -l`
   @ nend = `grep -r $finished $dir/job*/*.out | wc -l`
   if ( $nend < $njobs ) then
      echo "*** jobs not finished! ***"
   else
      if ( $expression == "" ) then
         set root_name = `find $dir/ -name "*.root" | head -n 1`
         set root_basename = `basename $root_name`
         hadd -j 4 $root_basename `find $dir/ -name "*.root"`
      else
         set root_basename = "$expression.root"
         hadd -j 4 $root_basename `find $dir/ -name "*$expression*.root"`
      endif
#      mv $root_basename $basedir
   endif
#   cd -
end

exit

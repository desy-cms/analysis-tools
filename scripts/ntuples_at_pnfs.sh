#!/usr/bin/csh -f

echo "Changing ntuples path from root://dcache-cms-xrootd.desy.de/ to /pnfs/desy.de/cms/tier2 for NAF usage..."

set mydir = $CMSSW_BASE/src/Analysis/Tools/data/ntuples/

set files = `find $mydir -name "*_rootFileList.txt"`

foreach f ( $files )
   sed -i 's~root://dcache-cms-xrootd.desy.de/~/pnfs/desy.de/cms/tier2~g' $f
end


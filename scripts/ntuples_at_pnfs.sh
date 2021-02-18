#!/usr/bin/csh -f

echo "Changing ntuples path from root://cms-xrd-global.cern.ch/ to /pnfs/desy.de/cms/tier2 for NAF usage..."

set mydir = $CMSSW_BASE/src/Analysis/Tools/data/ntuples/

set files = `find $mydir -name "*_rootFileList.txt"`

foreach f ( $files )
   sed -i 's~root://cms-xrd-global.cern.ch/~/pnfs/desy.de/cms/tier2~g' $f
end


#include "TRandom3.h"
#include "Analysis/Tools/interface/Analyser.h"

// this produces histogram with the true pileup

using namespace std;
using namespace analysis;
using namespace analysis::tools;

int main(int argc, char ** argv)
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   Analyser analyser(argc,argv);
   
//   analyser.pileupHistogram();
   
   TRandom3 * rnd = new TRandom3();
   
//   float ps = 1./analyser.config()->prescale();
   
   int n = analyser.config()->n();
   int evts_max = analyser.nEvents();
   int evts = evts_max;
//   float ps = float(n)/evts;
   
   if ( n < evts ) evts = n;
   
   // events loop
//   for ( int i = 0 ; i < analyser.nEvents() ; ++i )
   std::vector<int> v_evts;
   for ( int i = 0 ; i < evts ; ++i )
   {
      // avoiding bias for some samples
      double x; 
      if ( n < evts )
          x = rnd->Rndm()*(evts_max-1);
      else
         x = double(i);
//      if ( x > ps ) continue;
      int ii = int(x);
      if ( v_evts.size() > 0 )
      {
         if (std::find(v_evts.begin(), v_evts.end(),ii)!=v_evts.end())
         {
            std::cout << "repeated event" << std::endl;
            continue;
         }
      }
      v_evts.push_back(ii);
      
      analyser.event(ii);
      analyser.fillPileupHistogram();
      
   }  //end events loop
   
} // end main
      

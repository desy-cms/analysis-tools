#include "Analysis/Tools/interface/Analyser.h"

// this example has no selection
// for MC it will apply the generator weights

using namespace std;
using namespace analysis;
using namespace analysis::tools;

int main(int argc, char ** argv)
{
   Analyser analyser(argc,argv);
   auto config = analyser.config();
   
// HISTOGRAMS definitions  
   analyser.jetHistograms("nobtag_selection");
   analyser.jetHistograms("btagjet1_selection");
   analyser.jetHistograms("btagjet2_selection");
   
   for ( int i = 0 ; i < analyser.nEvents() ; ++i )
   {
      if ( ! analyser.event(i)                  )   continue;
      
   // JETS
      if ( ! analyser.selectionJetId()          )   continue;  // selection  : jet identification 
      if ( ! analyser.selectionJetPileupId()    )   continue;  // selection  : jet Pileup identification 
      if ( ! analyser.selectionNJets()          )   continue;  // selection  : number of jets 
      analyser.actionApplyJER();                               // correction : jet energy resolution smearing
      if ( ! analyser.selectionJet(1)           )   continue;  // selection  : jet1 pt and eta 
      if ( ! analyser.selectionJet(2)           )   continue;  // selection  : jet2 pt and eta 
      analyser.fillJetHistograms("nobtag_selection");               // histograms : jets fill
      if ( config->btagEfficiencies(1) == "" )
      {
         if ( ! analyser.selectionBJet(1)          )   continue;  // apply btag selection jet 1
         analyser.actionApplyBtagSF(1);
         analyser.fillJetHistograms("btagjet1_selection");               // histograms : jets fill
         if ( ! analyser.selectionBJet(2)          )   continue;  // apply btag selection jet 2
         analyser.actionApplyBtagSF(2);
         analyser.fillJetHistograms("btagjet2_selection");               // histograms : jets fill
      }
      else
      {
         analyser.actionApplyBtagEfficiency(1);  // weight according to btag efficiency
         analyser.fillJetHistograms("btagjet1_selection");               // histograms : jets fill
         analyser.actionApplyBtagEfficiency(2);  // weight according to btag efficiency
         analyser.fillJetHistograms("btagjet2_selection");               // histograms : jets fill
      }
      
   }  //end event loop
   

} // end main
      

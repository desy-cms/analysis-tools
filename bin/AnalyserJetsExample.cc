#include <iomanip>

#include "Analysis/Tools/interface/Analyser.h"

// this example has no selection
// for MC it will apply the generator weights

using namespace analysis;
using namespace analysis::tools;

int main(int argc, char ** argv)
{
   Analyser analyser(argc,argv);
   auto config = analyser.config();
   
   
// HISTOGRAMS definitions  
   // create some predefined jet histograms
   // if not defined, the number of jets is nJetMin from the configurations
   analyser.jetHistograms("initial");
   analyser.jetHistograms("final");
   
   auto output = analyser.output();
   output->cd();
   std::shared_ptr<TH1F>  hist = std::make_shared<TH1F>("event_weight","event weight", 200 ,-10,10);
   analyser.histogram("event_weight",hist);
   
   
   for ( int i = 0 ; i < analyser.nEvents() ; ++i )
   {
      if ( ! analyser.event(i)                  )   continue;
      output->cd();
      hist->Fill(analyser.weight());
      if ( ! analyser.selectionL1 ()            )   continue;    // L1  trigger (selection)
      if ( ! analyser.selectionHLT ()           )   continue;    // HLT trigger
      if ( ! analyser.selectionJetId()          )   continue;    // selection  : jet identification 
      if ( ! analyser.selectionJetPileupId()    )   continue;    // selection  : jet Pileup identification 
      if ( ! analyser.selectionNJets()          )   continue;    // selection  : number of jets 
      // Jet corrections
      analyser.actionApplyJEC();                                 // Jet energy scale (systematics)
      analyser.actionApplyJER();                                 // Jet energy resolution smearing (nominal+systematics)
      analyser.actionApplyBjetRegression();                      // Jet energy regression
      // Kinematic selections
      if ( ! analyser.selectionJet(1)           )   continue;    // selection  : jet1 pt and eta 
      if ( ! analyser.selectionJet(2)           )   continue;    // selection  : jet2 pt and eta 
      if ( ! analyser.selectionJetDphi(1,2)     )   continue;    // selection  : delta_phi_jets (1,2) [or  MIN(neg): analyser.selectionJetDphi(1,2,-2.0) / MAX(pos): analyser.selectionJetDphi(1,2,+2.0)]
      if ( ! analyser.onlineJetMatching(1)      )   continue;
      analyser.fillJetHistograms("initial");
      if ( ! analyser.onlineJetMatching(2)      )   continue;
      analyser.fillJetHistograms("final");
   }  //end event loop
   

} // end main
      

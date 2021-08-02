#ifndef Analysis_Tools_TriggerAnalyser_h
#define Analysis_Tools_TriggerAnalyser_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      TriggerAnalyser
// 
/**\class xxx TriggerAnalyser.h Analysis/Tools/interface/TriggerAnalyser.h

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh Bastos Rangel
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
//
//

// system include files
#include <memory>
#include <vector>
#include <string>
// 
#include "Analysis/Tools/interface/BaseAnalyser.h"
//
// class declaration
//

namespace analysis {
   namespace tools {

      class TriggerAnalyser  : virtual public analysis::tools::BaseAnalyser {
         
         public:
            /// constructors
            TriggerAnalyser();
            TriggerAnalyser(int argc, char * argv[]);
            /// desctructor
           ~TriggerAnalyser();
           
         
            // ----------member data ---------------------------
         protected:
            bool triggeranalysis_;
            bool l1tjetsanalysis_;
            bool l1tmuonsanalysis_;
            
            int   n_hl1tjets_;
         
         private:
               
         public:
            
            // Actions
            virtual bool analysisWithTrigger();
            virtual bool selectionTrigger();
            virtual bool selectionHLT();
            virtual bool selectionL1();
            
            bool l1tJetsAnalysis() const;
            bool l1tMuonsAnalysis() const;
            
            virtual bool selectionTriggerEmulated(const bool & , const bool &, const std::string& , const int &, const float &, const float &);
            
            /// Creates pre-defined histograms in directory 'label' for analysis with 'n' jets
            virtual void l1tjetHistograms(const int & n, const std::string & label = "x");
            /// Fill the pre-defined histograms created by the l1tjetHistograms() method
            virtual void fillL1TJetHistograms(const std::string & label = "x");
            
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsL1Jets();
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsCaloJets();
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsPFJets();
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsL1Muons();
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsL3Muons();

      };
   }
}

#endif  // Analysis_Tools_TriggerAnalyser_h

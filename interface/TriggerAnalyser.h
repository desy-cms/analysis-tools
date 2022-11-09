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
            std::vector< std::shared_ptr<L1TJet> > l1tjets_;
            std::vector< std::shared_ptr<L1TJet> > selected_l1tjets_;
            std::vector< std::shared_ptr<L1TMuon> > l1tmuons_;
            std::vector< std::shared_ptr<L1TMuon> > selected_l1tmuons_;
            
            std::vector<float> l1tjets_etabins_;
            std::vector<float> l1tjets_phibins_;
            std::vector<float> l1tmuons_etabins_;
            std::vector<float> l1tmuons_phibins_;
         
         private:
               
         public:
            
            // Actions
            virtual bool analysisWithTrigger();
            virtual bool selectionTrigger();
            virtual bool selectionHLT();
            virtual bool selectionL1();
            
            virtual bool analysisWithL1TJets();
            virtual bool analysisWithL1TMuons();
            
            virtual bool selectionNL1TJets(const int & nmin);
            virtual bool selectionL1TJet(const float & ptmin, const float & etamax);
            virtual bool selectionL1TDijet(const float & pt1min, const float & eta1max, const float & pt2min = -1, const float & eta2max = -1);
            virtual bool selectionL1TDijetDeta(const float & detamax);
            
            virtual bool selectionNL1TMuons(const int & nmin);
            virtual bool selectionL1TMuon(const float & ptmin, const float & etamax);
            virtual bool selectionL1TMuonQuality(const int & qual);
            
            virtual bool selectionL1TMuonJet(const float & drmax);
            
            bool l1tJetsAnalysis() const;
            bool l1tMuonsAnalysis() const;
            
            virtual bool selectionTriggerEmulated(const bool & , const bool &, const std::string& , const int &, const float &, const float &);
            
            /// Creates pre-defined histograms in directory 'label' for analysis with 'n' jets
            virtual void l1tjetHistograms(const std::string & label = "x");
            /// Fill the pre-defined histograms created by the l1tjetHistograms() method
            virtual void fillL1TJetHistograms(const std::string & label);
            /// Fill the pre-defined histograms created by the l1tjetHistograms() method
            virtual void fillL1TJetHistograms(const std::string & label, std::vector<std::shared_ptr<L1TJet> > sel_l1tjets);
            
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsL1Jets();
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsCaloJets();
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsPFJets();
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsL1Muons();
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsL3Muons();
            std::vector< std::shared_ptr<TriggerObject> > triggerObjectsBJets();
            
            std::vector< std::shared_ptr<L1TJet> > l1tJets();
            std::vector< std::shared_ptr<L1TJet> > selectedL1TJets();
            
            std::vector<float> l1tJetsEtaBins();
            std::vector<float> l1tJetsPhiBins();
            
            std::vector<float> l1tMuonsEtaBins();
            std::vector<float> l1tMuonsPhiBins();
            
      };
   }
}

#endif  // Analysis_Tools_TriggerAnalyser_h

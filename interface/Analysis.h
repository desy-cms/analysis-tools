#ifndef Analysis_Tools_Analysis_h
#define Analysis_Tools_Analysis_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      Analysis
//
/**\class xxx Analysis.h Analysis/Tools/interface/Analysis.h

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
#include <typeinfo>
#include <any>
#include <boost/core/demangle.hpp>
#include <boost/algorithm/string.hpp>
#include "stdlib.h"
//
// user include files

#include "TChain.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TH2.h"

#include "Analysis/Tools/interface/Utils.h"

#include "Analysis/Tools/interface/PhysicsObjectTree.h"
#include "Analysis/Tools/interface/Collection.h"
#include "Analysis/Tools/interface/BTagCalibrationStandalone.h"
#include "Analysis/Tools/interface/PileupWeight.h"
#include "Analysis/Tools/interface/MuonIdWeight.h"
#include "Analysis/Tools/interface/JetOnlineTriggerEfficiencies.h"
#include "Analysis/Tools/interface/OnlineBtagTriggerEfficiencies.h"
#include "Analysis/Tools/interface/MuonOnlineTriggerEfficiencies.h"

#include "Analysis/Tools/interface/Config.h"

//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//
// class declaration
//

using namespace JME;

namespace analysis {
   namespace tools {

      class Analysis
      {
         public:
            Analysis(const std::string & inputFilelist, const std::string & evtinfo = "MssmHbb/Events/EventInfo");
            Analysis(std::shared_ptr<Config> cfg);
           ~Analysis();
           
            // Info
            void tag(const std::string &);
            std::string tag();
            
            /// Configuration
            void config(std::shared_ptr<Config> cfg);
            
            /// seed for random number generator read from a txt file given as a parameter
            int seed(const std::string &);

            // Event
            int  numberEvents();
            int  size();
            void event(const int & event, const bool & addCollections = true);
            int  event();
            int  run();
            int  lumiSection();
            bool isMC();
            std::string fileName();
            std::string fileFullName();
            
            // PileupInfo
            int   nPileup();
            float nTruePileup();
            
            float lumiPileup();
            float instantLumi();
            
            // GenEventInfo
            double genWeight();
            double genScale();
            PDF    pdf();
            
            /// fixedGridRhoAll
            double rho();

            // Trees
            template<class Object>
            std::shared_ptr< PhysicsObjectTree<Object> > addTree(const std::string & unique_name, const std::string & path );
            template<class Object>
            std::shared_ptr< PhysicsObjectTree<Object> > tree(const std::string & unique_name);
            
            // Collections
            template<class Object>
            std::shared_ptr< Collection<Object> > addCollection(const std::string & unique_name);
            template<class Object>
            std::shared_ptr< Collection<Object> > addCollection(const Collection<Object> & collection);
            template<class Object>
            std::shared_ptr< Collection<Object> > addCollection(const std::vector<Object> & objects, const std::string & unique_name );
            template<class Object>
            std::shared_ptr< Collection<Object> > collection(const std::string & unique_name);
            
            template<class Object>
            void defaultCollection(const std::string & unique_name);
            template<class Object>
            std::string defaultCollection();
            
            // Cross sections
            int    crossSections(const std::string & path);
            double crossSection();
            double crossSection(const std::string & title);
            void   listCrossSections();
            
            // Luminosity
            double luminosity();
            double luminosity(const std::string & title);

            // Trigger results
            /// Open the TriggerResults tree and make it Friend, loop over all triggers, HLT and L1, retrieving if fired and prescale
            bool triggerResults(const std::string & path);
            /// Given a trigger name, HLT or L1 returns if fired or not
            bool triggerResult(const std::string & trig);
            /// Given a trigger name, HLT or L1 returns the prescale
            int triggerPrescale(const std::string & trig);
            /// Given a list of trigger names, HLT and L1, returns a map containing the prescales
            std::map<std::string,int> triggerPrescale(const std::vector<std::string> & trigs);

            // Generator Filter
            FilterResults generatorFilter(const std::string & path);
            void listGeneratorFilter();
            
            // Event Filter
            FilterResults eventFilter(const std::string & path);
            
            // Matching to trigger objects
            template <class Object1, class Object2>
            void match(const std::string & collection, const std::string & match_collection, const float & deltaR = 0.5);
            template <class Object1, class Object2>
            void match(const std::string & collection, const std::vector<std::string> & match_collections, const float & deltaR = 0.5);
            
            // good Json files
            int processJsonFile(const std::string & fileName = "goodJson.txt");
            bool selectJson();
            
            // btag efficiencies
            void addBtagEfficiencies(const std::string & );
            float btagEfficiency(const analysis::tools::Jet &, const int & rank = 0);
            void  btagEfficienciesAlgo(const std::string & );
            void  btagEfficienciesFlavour(const std::string & );
            
            std::shared_ptr<JetResolutionInfo> jetResolutionInfo(const std::string &, const std::string & );
            
            std::shared_ptr<BTagCalibrationReader> btagCalibration(const std::string & tagger,
                                 const std::string & filename,
                                 const std::string & wp,
                                 const std::string & sysType="central",
                                 const std::vector<std::string> & otherSysTypes={"up","down","up_correlated","down_correlated","up_uncorrelated","down_uncorrelated"});
            
            
            std::shared_ptr<BTagCalibrationReader> btagCalibration();
            
            float scaleLuminosity(const float & lumi);  // in pb-1

            std::shared_ptr<PileupWeight> pileupWeights(const std::string & , const std::string &, const bool &);
            std::shared_ptr<MuonIdWeight> muonIDWeights(const std::vector<std::string> &);

            /// L1 prefiring weights
            double prefiringWeight(const int & var = 0);

            // ----------member data ---------------------------
         protected:

            TFileCollection * fileCollection_;
            TCollection * fileList_;
            std::string inputFilelist_;
            
            // Configuration
            std::shared_ptr<Config>   config_;
            
            // Info
            std::string tag_;
            
            // btagging efficiencies
            TFile * fileBtagEff_;
            std::map<std::string,TH2F *> h2_btageff_;
            std::string btageff_flavour_;
            std::string btageff_algo_;
            
            std::shared_ptr<BTagCalibration> btagcalib_;
            std::shared_ptr<BTagCalibrationReader> btagcalibread_;
            
            std::shared_ptr<JetResolutionInfo> jerinfo_;
            
            // pileup weight
            std::shared_ptr<PileupWeight> puweights_;
	    
	         // muonID weight
	         std::shared_ptr<MuonIdWeight> muonIDweights_;

            std::map<std::string, double> xsections_;
            std::map<std::string, bool> triggerResults_;
            std::map<std::string, int> triggerResultsPS_;
            std::map<int,std::vector<std::string> > goodLumi_;
            FilterResults genfilter_;
            FilterResults evtfilter_;
            
            // default collections
            std::string defaultGenParticle_;

            int event_;
            int run_;
            int lumi_;
            bool is_mc_;
            
            int n_pu_;
            float n_true_pu_;
            
            float lumi_pu_;
            float inst_lumi_;
            
            double genWeight_;
            double genScale_;
            PDF    pdf_;
            
            double rho_;


            // L1 prefiring weights
            double prefw_;
            double prefw_up_;
            double prefw_down_;

            int nevents_;

            // TREES
            void treeInit_(const std::string & unique_name, const std::string & path);
            TChain * t_xsection_;
            TChain * t_genfilter_;
            TChain * t_evtfilter_;
            TChain * t_event_;
            TChain * t_triggerResults_;

         // Physics objects
            // root trees
            std::map<std::string, TChain*> tree_;

            // Framework trees and types
            std::map<std::string, std::any  > t_any_;
            std::map<std::string, std::string > t_type_;
            
            // Collections
            std::map<std::string, std::any > c_any_;
            
            // Luminosity
            float mylumi_;
            
            // json
            std::map<int,std::vector<int> > json_;
            
         private:


      }; // END OF CLASS DECLARATIONS!!!

// ========================================================
//                         IMPLEMENTATIONS!
// ========================================================

// +++++++++++++++++++++++ IMPORTANT ++++++++++++++++++++++
// Need to put the implementations in the header file when using template!!!
// otherwise would have to specialize, like the other functions.
// This explains the problems I have been getting.
// -------------------------------------------------------
      // TREES
      template <class Object>
      std::shared_ptr< PhysicsObjectTree<Object> >  Analysis::addTree(const std::string & unique_name, const std::string & path)
      {
         if ( path == "" || unique_name == "" ) return nullptr;
         this->treeInit_(unique_name,path);
         t_any_[unique_name] = std::shared_ptr< PhysicsObjectTree<Object> > ( new PhysicsObjectTree<Object>(tree_[unique_name], unique_name) );
         std::string type = boost::core::demangle(typeid(Object).name());
         std::vector<std::string> tmp;
         boost::split( tmp, type, boost::is_any_of("::"));
         t_type_[unique_name] = tmp.back();
         return std::any_cast< std::shared_ptr< PhysicsObjectTree<Object> > > (t_any_[unique_name]);
      }
      // --
      template <class Object>
      std::shared_ptr< PhysicsObjectTree<Object> >  Analysis::tree(const std::string & unique_name)
      {
         // If tree does not exist, return NULL
         std::map<std::string, std::any >::iterator it = t_any_.find(unique_name);
         if ( it == t_any_.end() )
            return nullptr;
         return std::any_cast< std::shared_ptr< PhysicsObjectTree<Object> > > (t_any_[unique_name]);
      }
// -------------------------------------------------------
      // COLLECTIONS
      template <class Object>
      std::shared_ptr< Collection<Object> >  Analysis::addCollection(const std::string & unique_name)
      {
         // Still need to see how to deal with collections not originating from the ntuple,
         // e.g. a selected jets collection from the ntuple jets collection.
         
         // If tree does not exist, return NULL
         std::map<std::string, std::any >::iterator it = t_any_.find(unique_name);
         if ( it == t_any_.end() )
            return nullptr;
         
         auto tree = std::any_cast< std::shared_ptr< PhysicsObjectTree<Object> > > (t_any_[unique_name]);
         c_any_[unique_name] = std::shared_ptr< Collection<Object> > ( new Collection<Object>(tree -> collection()));
         std::shared_ptr< Collection<Object> > ret = std::any_cast< std::shared_ptr< Collection<Object> > > (c_any_[unique_name]);
         
         return ret;
      }
      
      template <class Object>
      std::shared_ptr< Collection<Object> >  Analysis::addCollection(const Collection<Object> & collection)
      {
         std::string unique_name = collection.name();
         t_any_[unique_name] = nullptr;
         c_any_[unique_name] = std::shared_ptr< Collection<Object> > ( new Collection<Object>(collection) );
         std::shared_ptr< Collection<Object> > ret = std::any_cast< std::shared_ptr< Collection<Object> > > (c_any_[unique_name]);
         
         return ret;
      }
      
      template <class Object>
      std::shared_ptr< Collection<Object> >  Analysis::addCollection(const std::vector<Object> & objects , const std::string & unique_name )
      {
         Collection<Object> collection(objects,unique_name);
         t_any_[unique_name] = nullptr;
         c_any_[unique_name] = std::shared_ptr< Collection<Object> > ( new Collection<Object>(collection) );
         std::shared_ptr< Collection<Object> > ret = std::any_cast< std::shared_ptr< Collection<Object> > > (c_any_[unique_name]);
         return ret;
      }
      
      template <class Object>
      std::shared_ptr< Collection<Object> >  Analysis::collection(const std::string & unique_name)
      {
         std::shared_ptr< Collection<Object> > ret = std::any_cast< std::shared_ptr< Collection<Object> > > (c_any_[unique_name]);
         return ret;
      }
      //--
      template <class Object1, class Object2>
      void Analysis::match(const std::string & collection, const std::string & match_collection, const float & deltaR)
      {
         if ( match_collection == "" ) return;
         auto o1 = std::any_cast< std::shared_ptr< Collection<Object1> > > (c_any_[collection]);
         auto o2 = std::any_cast< std::shared_ptr< Collection<Object2> > > (c_any_[match_collection]);
         o1->matchTo(o2->vectorCandidates(),o2->name(), deltaR);
      }
      //--
      template <class Object1, class Object2>
      void Analysis::match(const std::string & collection, const std::vector<std::string> & match_collections, const float & deltaR)
      {
         for ( auto & mc : match_collections )
            this->match<Object1,Object2>(collection, mc, deltaR);
      }

      template<class Object> void Analysis::defaultCollection(const std::string & unique_name)
      { 
         if ( std::is_same<Object,GenParticle>::value ) defaultGenParticle_ = unique_name; 
      }
      
      template<class Object> std::string Analysis::defaultCollection()
      { 
         std::string ret;
         if ( std::is_same<Object,GenParticle>::value ) ret = defaultGenParticle_ ;
         return ret; 
      }
      
// ========================================================
      
   }
}

#endif  // Analysis_Tools_Analysis_h

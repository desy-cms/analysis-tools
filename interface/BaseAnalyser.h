#ifndef Analysis_Tools_BaseAnalyser_h
#define Analysis_Tools_BaseAnalyser_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      Analysis
//
/**

 Description: A base class for analyser

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh
//         Created:  6 Sep 2018
//
//

// system include files
#include <memory>
#include <vector>
#include <string>
//
// user include files
#include "boost/program_options.hpp"

#include "Analysis/Tools/interface/Analysis.h"
#include "Analysis/Tools/interface/Config.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"

namespace po = boost::program_options;


//
// class declaration
//

namespace analysis {
   namespace tools {

      class BaseAnalyser {

         public:
            /// constructor
            BaseAnalyser();
            /// constructor
            BaseAnalyser(int argc, char * argv[]);
            /// desctructor
           virtual ~BaseAnalyser();


            // ----------member data ---------------------------
         protected:
            /// Analysis objects
            std::shared_ptr<Analysis> analysis_;
            /// Config objects
            std::shared_ptr<Config>   config_;

            /// Cutflow index
            int cutflow_;

            /// seed value
            int seed_;

            /// event weight
            float weight_;

            /// cross section
            float xsection_;

            /// output root file
            std::shared_ptr<TFile> hout_;
            /// 1D histograms mapping
            std::map<std::string, std::shared_ptr<TH1F> > h1_;
            /// 2D histograms mapping
            std::map<std::string, std::shared_ptr<TH2F> > h2_;

            std::map<std::string, std::shared_ptr<TGraphAsymmErrors> > btageff_[3];

            bool genpartsanalysis_;
            bool genjetsanalysis_;

            /// primary vertex
            bool primaryvtxanalysis_;

            /// pileup weight
            std::shared_ptr<PileupWeight> puweights_;
            bool puweights_ok_;
            
            /// pileup weight label
            std::string puw_label_;
            
            /// true pileup
            bool istruepu_;

            /// overall scaling
            float scale_;

            /// flag for MC sample
            bool isMC_;
            /// flag for DATA sample
            bool isData_;

            /// output tree
            std::shared_ptr<TTree> analyser_tree_;
            
            /// emulated triggers
            std::map<std::string,bool> trg_emul_;

            /// external scale data
            std::map<std::string, std::vector<float> > scale_data_;
            float scale_correction_;

            /// is muon analysis
            bool is_mouns_analysis_;

         private :

            /// name of the executable
            std::string exe_;

         public:
            // Gets
            /// returns pointer to Analysis object
            std::shared_ptr<Analysis> analysis();
            /// returns pointer to Config object
            std::shared_ptr<Config>   config();

            /// number of events to be processed
            int nEvents();


            /// returns all 1D histograms
            std::map<std::string, std::shared_ptr<TH1F> > histograms();
            /// returns a given 1D histogram
            std::shared_ptr<TH1F> histogram(const std::string &);
            /// 
            void histogram(const std::string &, std::shared_ptr<TH1F>);

            /// get cutflow index
            int  cutflow();
            /// set cutflow index
            void cutflow(const int &);
            /// create and update cutflow entry in the cutflow histogram 
            void cutflow(const std::string & label, const bool & ok = true);

            /// get is muons analysis
            bool isMuonsAnalysis();
            /// get is muons analysis
            void isMuonsAnalysis(const bool &);

            // Actions
            /// event entry to be readout and processed
            virtual bool event(const int &);
            /// create n histograms of a given type
            virtual void histograms(const std::string &);

            /// returns a seed for random number generator
            int  seed();
            /// reads a file containing a seed and returns the seed or -1 if fails
            int  seed(const std::string &);
            /// sets a seed for random number generators
            void seed(const int &);

            /// returns event weight
            float weight();
            /// sets event weight
            void weight(const float &);
            
            /// generator weight
            void generatorWeight();

            /// returns cross section
            float crossSection() const;

            /// returns pointer to output root file
            std::shared_ptr<TFile> output();

            /// returns whether analysis of gen particles can be done
            bool genParticlesAnalysis() const;
            /// returns whether analysis of gen jets can be done
            bool genJetsAnalysis() const;
            /// returns whether analysis of primary vertex can be done
            bool primaryVertexAnalysis() const;

            /// returns pointer to pileup weights (MC-only)
            std::shared_ptr<PileupWeight> pileupWeights() const;
            /// returns pileup weight given the true pileup and uncertainty variation in values of sigma
            float pileupWeight(const float & truepu, const int & var) const;
            /// returns true number of interactions
            float trueInteractions() const;
            /// check conditions for pileup weights is fulfilled
            bool pileupWeightsApplicable() const;

            /// creates pileup histogram
            void pileupHistogram();
            /// fills pileup histogram
            void fillPileupHistogram();

            /// sets a scale
            void scale(const float &);

            /// returns the basename of a path
            std::string basename(const std::string &);

            /// apply pileup weight given a systematic variation
            void actionApplyPileupWeight(const int & var);
            void actionApplyPileupWeight();

            /**
            btag efficiencies

            Given as TGraphAsymmErrors for each flavour
            */
            std::map<std::string, std::shared_ptr<TGraphAsymmErrors> > btagEfficiencies(const int & model) const;

            /// root tree
            void analyserTree();
            /// fill root tree
            void fillAnalyserTree();

            /// emulate l1 trigger
            virtual bool triggerEmulation(const std::string &, const int &, const float &, const float &, const std::string & );
            
            /// return status of all emulated l1 triggers
            std::map<std::string, bool> triggersEmulated();
            
            /// return status of an emulated l1 trigger 
            bool triggerEmulated(const std::string &);

            /// apply prefiring weight given a systematic variation
            void actionApplyPrefiringWeight(const int & var);
            void actionApplyPrefiringWeight();

            /// add histogram
            virtual void add1DHistogram(const std::string & label, const std::string & name, const std::string & title, const int & nbins, const float & min, const float & max);
            /// fill added histogram 
            virtual void fill1DHistogram(const std::string & label, const std::string & name, const float & value, const float & weight = 1.);

            /// apply scale correction
            void actionApplyScaleCorrection(const std::string & title="");


      };
   }
}

#endif  // Analysis_Tools_BaseAnalyser_h

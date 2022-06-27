#ifndef Analysis_Tools_Config_h
#define Analysis_Tools_Config_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      xxx
//
/**\class xxx Config.cc Analysis/Tools/src/Config.cc

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
// user include files
#include "boost/program_options.hpp"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"


namespace po = boost::program_options;
using TH1s = std::map<std::string, TH1F*>;
using TH2s = std::map<std::string, TH2F*>;

//
// class declaration
//

namespace analysis {
   namespace tools {

      struct Options
      {
         po::options_description cmd;
         po::options_description cfg;
      };

      class Config {
         public:
            /// constructors
            Config();
            Config(int argc, char ** argv);
            /// desctructor
           ~Config();

            po::options_description & optionsCMD();
            po::options_description & optionsCFG();

            void loadOptions();

            // configuration variables - some basic ones out in the private section and make set/get methods

         // gets (to replace the public variables)
            std::string configFile() const;
         // analysis control
            std::string ntuplesList() const;
            std::string sampleName() const;
            std::string eventInfo() const;
            std::string crossSectionTree() const;
            std::string crossSectionType() const;
            float crossSection() const;
            float luminosity() const;
            int  nEventsMax() const;
            bool isMC() const;
            bool signalRegion() const;
            bool blind() const;
            bool override() const;
            bool nlo() const;
            bool fullGenWeight() const;
            int  workflow() const;
            int  index() const;
            
            std::string process() const;
            std::string eventsDir() const;

            std::string seedFile() const;
            int seed() const;

            bool pythia8() const;

            float scale() const;

            std::vector<float> erasLumi() const;
            std::vector<std::string> eras() const;

            std::string pileupWeights() const;
            int         pileupWeightSystematics() const;

         // jets
            std::string jetsCollection() const;
            int nJetsMin() const;
            int nJetsMax() const;
            int nJets() const;
            std::vector<float> jetsPtMin() const;
            std::vector<float> jetsPtMax() const;
            std::vector<float> jetsEtaMax() const;
            std::string jetsId() const;
            std::string jetsPuId() const;
            std::string jerPtRes() const;
            std::string jerSF() const;
            int         jerSystematics() const;
            int         jecSystematics() const;
            
            std::string l1tJetsCollection() const;

            std::string btagAlgorithm() const;
            std::string btagScaleFactors() const;
            std::vector<std::string> jetsBtagWP() const;
            std::vector<float> jetsBtagProbB() const;
            std::vector<float> jetsBtagProbBB() const;
            std::vector<float> jetsBtagProbLepB() const;
            std::vector<float> jetsBtagProbC() const;
            std::vector<float> jetsBtagProbG() const;
            std::vector<float> jetsBtagProbLight() const;
            bool bRegression() const;
            std::string revBtagWP()  const;
            int   revBtagJet() const;
            bool useJetsExtendedFlavour() const;
            bool doDijet() const;
            bool doDijetFlavour() const;
            int nBJetsMin() const;
            bool prefiringWeight() const;
            int  prefiringWeightSystematics() const;

            std::vector<float>  jetsQGmin() const;
            std::vector<float>  jetsQGmax() const;



         // jet-jet
            float jetsDetaMax()          const;
            float jetsDetaMin()          const;
            float jetsDphiMax()          const;
            float jetsDphiMin()          const;
            float jetsDrMax()            const;
            float jetsDrMin()            const;
            float jetsPtImbalanceMax()   const;
            float jetsPtImbalanceMin()   const;
            bool  jetsWithMuons()        const;
            float jetsMuonsDRMax()       const;

         // muons
            std::string muonsCollection() const;
            int nMuonsMin() const;
            int nMuonsMax() const;
            std::vector<float> muonsPtMin() const;
            std::vector<float> muonsPtMax() const;
            std::vector<float> muonsEtaMax() const;
            std::string muonsId() const;
            bool muonsVeto() const;
            std::string l1tMuonsCollection() const;

         // muon-muon
            float muonsDrMax()            const;
            float muonsDrMin()            const;

         // trigger
            std::string triggerResults()         const;
            std::string triggerObjectsDir()      const;
            std::string triggerObjectsL1Muons()  const;
            std::string triggerObjectsL3Muons()  const;
            std::string triggerObjectsBJets()    const;
            std::string triggerObjectsL1Jets()   const;
            std::string triggerObjectsCaloJets() const;
            std::string triggerObjectsPFJets()   const;
            int         triggerObjectsNJets()    const;
            int         triggerObjectsNBJets()   const;
            int         triggerObjectsNMuons()   const;

            /// L1 muon trigger emulation
            std::string triggerEmulateL1Muons()       const;
            int         triggerEmulateL1MuonsNMin()   const;
            float       triggerEmulateL1MuonsPtMin()  const;
            float       triggerEmulateL1MuonsEtaMax() const;

            /// L3 muon trigger emulation
            std::string triggerEmulateL3Muons()       const;
            int         triggerEmulateL3MuonsNMin()   const;
            float       triggerEmulateL3MuonsPtMin()  const;
            float       triggerEmulateL3MuonsEtaMax() const;


            /// L1 jet trigger emulation
            std::string triggerEmulateL1Jets()       const;
            int         triggerEmulateL1JetsNMin()   const;
            float       triggerEmulateL1JetsPtMin()  const;
            float       triggerEmulateL1JetsEtaMax() const;
            
            /// Calo jet trigger emulation
            std::string triggerEmulateCaloJets()       const;
            int         triggerEmulateCaloJetsNMin()   const;
            float       triggerEmulateCaloJetsPtMin()  const;
            float       triggerEmulateCaloJetsEtaMax() const;
            
            
            /// PF jet trigger emulation
            std::string triggerEmulatePFJets()       const;
            int         triggerEmulatePFJetsNMin()   const;
            float       triggerEmulatePFJetsPtMin()  const;
            float       triggerEmulatePFJetsEtaMax() const;
            
            
         // matching trigger objects
            float triggerMatchL1MuonsDrMax()       const;
            float triggerMatchL3MuonsDrMax()       const;
            float triggerMatchL1JetsDrMax()        const;
            float triggerMatchCaloJetsDrMax()      const;
            float triggerMatchPFJetsDrMax()        const;
            float triggerMatchCaloBJetsDrMax()     const;

         // generator level
            std::string genJetsCollection() const;
            std::string genParticlesCollection() const;

         // btag
            float btagWP(const std::string &) const;

         // general
            float massMin() const;
            float massMax() const;

         // AI
            std::vector<std::string> variablesAI(const std::string & t = "F") const;
            std::string directoryAI() const;
            std::string methodAI() const;
            float discriminatorMaxAI() const;
            float discriminatorMinAI() const;
            float efficiencyMinAI() const;

         // output tree
            bool doTree() const;

         // generic options
            int prescale() const;
            int n() const;
            float min() const;
            float max() const;

         // histograms
            bool histogramJetsRegionSplit() const;
            bool histogramJetsPerFlavour() const;

         // btag
            std::string btagEfficiencies(const int & model) const;

         // ========================

         // analysis control
            
            int nlumis_;
            int runmin_;
            int runmax_;
            std::string outputRoot() const;
            std::string json() const;


            //
            float trgmatchdrmax_;


            // btag SF csv file
            std::vector<float> jetsbtagmin_;

            // additional cuts of unidentified objects or for extra selections
            int nmin_;
            int nmax_;
            std::vector<float> ptmin_;
            std::vector<float> ptmax_;
            std::vector<float> etamax_;



            float drmin_;
            float drmax_;
            float detamax_;
            float detamin_;
            float dphimin_;
            float dphimax_;

            std::string hltPath_;
            std::string l1Seed_;
            
            // User vectors
            std::vector<float> vectorFloat() const;
            std::vector<int>   vectorInt() const;
            

            // ----------member data ---------------------------
         protected:

            int argc_;
            char ** argv_;

            std::string cfg_; // config file
            
            bool cmdl_mc_;
            bool cmdl_data_;
            bool cmdl_bweight_;
            std::string cmdl_inputlist_;
            int cmdl_evtmax_;
            int cmdl_jer_;
            int cmdl_jec_;

            po::options_description opt_cmd_;
            po::options_description opt_cfg_;

         // analysis control
            std::string inputlist_;
            std::string samplename_;
            std::string eventinfo_;
            std::string process_;
            std::string eventsdir_;
            std::string xsectiontree_;
            std::string xsectiontype_;
            float xsection_;
            float lumi_;
            int nevtmax_;
            bool isMC_;
            bool signalregion_;
            bool blind_;
            bool override_;
            bool nlo_;
            bool fullgenweight_;
            int workflow_;
            int index_;

            bool apply_correct_;

            int seed_;
            std::string seedfile_;

            bool pythia8_;

            float scale_;

            std::vector<std::string> eras_;
            std::vector<float> eraslumi_;

            std::string puweight_;
            int puweightsyst_;

            int prefwsyst_;

         // generic options
            int prescale_;
            int n_;
            float min_;
            float max_;


         // jets
            std::string jetsCol_;
            int njetsmin_;
            int njetsmax_;
            int njets_;
            std::vector<float> jetsptmin_;
            std::vector<float> jetsptmax_;
            std::vector<float> jetsetamax_;
            std::string jetsid_;
            std::string jetspuid_;

            std::vector<float>  qgmin_;
            std::vector<float>  qgmax_;

            float jetsdetamax_;
            float jetsdetamin_;
            float jetsdphimin_;
            float jetsdphimax_;
            float jetsdrmin_;
            float jetsdrmax_;
            float jetsptimbalmax_;
            float jetsptimbalmin_;
            bool  jetswithmuons_;
            float jetsmuonsdrmax_;

            // muon-muons
            float muonsdrmin_;
            float muonsdrmax_;

            // JER resolution and scale factors from txt file
            std::string jerptres_;
            std::string jersf_;
            int         jersyst_;
            // JEC systematics
            int         jecsyst_;
            //
            std::string l1tjetsCol_;
            //
            //
            std::string btagalgo_;
            std::string btagsf_;


            std::vector<std::string> jetsbtagwp_;
            std::string revbtagwp_;
            int revbtagjet_;
            // cuts on the probabilities
            std::vector<float> jetsbtagprobb_;
            std::vector<float> jetsbtagprobbb_;
            std::vector<float> jetsbtagproblepb_;
            std::vector<float> jetsbtagprobc_;
            std::vector<float> jetsbtagprobg_;
            std::vector<float> jetsbtagproblight_;

            //
            bool usejetsextflv_;
            bool dodijet_;

            std::vector<int> dijet_ranks_;

            bool bregression_;

            bool prefw_;


         // muons
            std::string muonsCol_;
            int nmuonsmin_;
            int nmuonsmax_;
            std::vector<float> muonsptmin_;
            std::vector<float> muonsptmax_;
            std::vector<float> muonsetamax_;
            std::string muonsid_;
            std::string l1tmuonsCol_;
            bool muonsveto_;

         // trigger
            std::string triggerCol_;
            std::string triggerObjDir_;
            std::string trgObjsL1Muons_;
            std::string trgObjsL3Muons_;
            std::string trgObjsBJets_;
            std::string trgObjsL1Jets_;
            std::string trgObjsCaloJets_;
            std::string trgObjsPFJets_;
            int         trgObjsNJets_;
            int         trgObjsNBJets_;
            int         trgObjsNMuons_;

         /// L1 muon trigger emulation
            std::string l1muonemul_;
            int         l1muonemulnmin_;
            float       l1muonemulptmin_;
            float       l1muonemuletamax_;

         /// L3 muon trigger emulation
            std::string l3muonemul_;
            int         l3muonemulnmin_;
            float       l3muonemulptmin_;
            float       l3muonemuletamax_;

         /// L1 jet trigger emulation
            std::string l1jetemul_;
            int         l1jetemulnmin_;
            float       l1jetemulptmin_;
            float       l1jetemuletamax_;

         /// Calo jet trigger emulation
            std::string calojetemul_;
            int         calojetemulnmin_;
            float       calojetemulptmin_;
            float       calojetemuletamax_;

         /// PF jet trigger emulation
            std::string pfjetemul_;
            int         pfjetemulnmin_;
            float       pfjetemulptmin_;
            float       pfjetemuletamax_;

            float matchTrgL1MuonsDrMax_;
            float matchTrgL3MuonsDrMax_;
            float matchTrgL1JetsDrMax_;
            float matchTrgCaloJetsDrMax_;
            float matchTrgPFJetsDrMax_;
            float matchTrgCaloBJetsDrMax_;

         // generator level collections
            std::string genjetsCol_;
            std::string genpartsCol_;



         // btag
            float btagwploose_;
            float btagwpmedium_;
            float btagwptight_;
            float btagwpxxx_;
            int nbjetsmin_;
            // btag efficiency
            std::string btageff_[3];


         // general
            float massmin_;
            float massmax_;

         // AI
            std::vector<std::string> varsf_ai_;
            std::vector<std::string> varsi_ai_;
            std::string dir_ai_;
            std::string method_ai_;
            float disc_max_ai_;
            float disc_min_ai_;
            float eff_min_ai_;
         // output tree

            bool do_tree_;

         // histograms
            bool histjets_rsplit_;
            bool histjets_flavour_;

         // analysis control
            std::string outputRoot_;
            std::string json_;
            
         // user vector
            std::vector<float> vfloat_;
            std::vector<int>   vint_;

         private:


      };
   }
}

#endif  // Analysis_Tools_Config_h

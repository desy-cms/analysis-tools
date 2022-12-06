// system include files
#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
// 
// user include files
#include "Analysis/Tools/interface/Analyser.h"

//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;

Analyser::Analyser()
{
}


Analyser::Analyser(int argc, char * argv[]) : BaseAnalyser(argc,argv),  // not sure the BaseAnalyser should be called here
                                              TriggerAnalyser(argc,argv),
                                              JetAnalyser(argc,argv),
                                              MuonAnalyser(argc,argv)
{
}

Analyser::~Analyser()
{
}


//
// member functions
//
// ------------ method called for each event  ------------

bool Analyser::event(const int & i)
{
   bool ok = true;
   analysis_->event(i);
   cutflow_ = -1;
   weight_ = 1.;  // reset weight at the beginning of the event analysis
  
   if ( config_->sampleName() != "" ) 
      cutflow(config_->sampleName());
   else
      cutflow("Total events read");
   
   // Generator weight
   if ( isMC_ )
   {
      this -> generatorWeight();
      std::string genweight_type = "sign of weights";
      if ( config_->fullGenWeight() ) genweight_type = "full weights";
      
      cutflow(Form("Generated weighted events (%s)",genweight_type.c_str()));
   }
   
   if ( config_->runmin_ > 0 && analysis_->run() < config_->runmin_ ) return false;
   if ( config_->runmax_ > 0 && analysis_->run() > config_->runmax_ ) return false;

   if (! config_->isMC() && config_->json() != "" ) 
   {
       auto json = basename(config_->json());
       ok = analysis_->selectJson();
       cutflow(Form("Certified data: %s",json.c_str()),ok);
       if ( ! ok ) return false;
       
   }
   
   if ( this->genParticlesAnalysis() )
      cutflow(Form("Using GenParticles collection: %s",(config_->genParticlesCollection()).c_str()));
   
   if ( this->genJetsAnalysis() )
      cutflow(Form("Using GenJets collection: %s",(config_->genJetsCollection()).c_str()));

   if ( this->primaryVertexAnalysis())
      cutflow(Form("Using Primary Vertex Collection: %s",(config_->primaryVertexCollection()).c_str()));
      
   if ( this -> l1tJetsAnalysis() )
      cutflow(Form("Using L1TJet collection: %s", (config_->l1tJetsCollection()).c_str()));

   if ( this -> l1tMuonsAnalysis() )
      cutflow(Form("Using L1TMuon collection: %s", (config_->l1tMuonsCollection()).c_str()));
      
   analysisWithJets();
   analysisWithMuons();
   analysisWithL1TJets();
   analysisWithL1TMuons();
      
   // PILEUP RE-WEIGHT
   this->actionApplyPileupWeight(); 

   return ok;
   
}

bool Analyser::muonJet(const int & r)
{
   
   if ( ! muonsanalysis_ ) return true;  // will skip this selection

   int j = r-1;
   auto jet = selectedJets_[j];
   jet -> addMuon(selectedMuons_,config_->jetsMuonsDRMax());
   bool isMuonJet = (jet -> muon() != nullptr);
   cutflow(Form("Jet %d: Jet-muon association (delta_R < %4.2f)",r, config_->jetsMuonsDRMax()),isMuonJet);
   
   return isMuonJet;
   
}

bool Analyser::preselection()
{
// IDENTIFICATIONS
      if ( ! this->selectionMuonId()         )   return false;
      if ( ! this->selectionJetId()          )   return false;
      if ( ! this->selectionJetPileupId()    )   return false;
      return true;
      
}

bool Analyser::selectionPrimaryVertex()
{
   if ( ! this->primaryVertexAnalysis() ) return true;

   auto pvs = analysis_->collection<Vertex>("PrimaryVertex");
   std::shared_ptr<Vertex> pv = std::make_shared<Vertex>(pvs->at(0));
   auto fake = pv->fake();
   auto rho  = pv->rho();
   auto ndof = pv->ndof();
   auto absz = fabs(pv->z());

   // a bit confusing, maybe ok
   std::string notfake = "!fake";
   if ( ! config_->primaryVertexNotFake() ) notfake = "fake";
   std::string label = Form("Primary vertex selection: %s && ndof>%3.1f && |z|<%3.1f && rho<%3.1f", notfake.c_str(),config_->primaryVertexNdofMin(),config_->primaryVertexAbsZMax(),config_->primaryVertexRhoMax());

   bool goodPV = (           (!fake && config_->primaryVertexNotFake()) );
   goodPV =      ( goodPV && (ndof   > config_->primaryVertexNdofMin()) );
   goodPV =      ( goodPV && (absz   < config_->primaryVertexAbsZMax()) ); 
   goodPV =      ( goodPV && (rho    < config_->primaryVertexRhoMax() ) );

   cutflow(label, goodPV);

   return goodPV;
}

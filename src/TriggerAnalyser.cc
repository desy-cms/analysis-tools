#include "Analysis/Tools/interface/TriggerAnalyser.h"
// system include files
#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
// 
// user include files
#include "TString.h" 

//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;

TriggerAnalyser::TriggerAnalyser()
{
}


TriggerAnalyser::TriggerAnalyser(int argc, char * argv[]) : BaseAnalyser(argc,argv)
{

   triggeranalysis_ = false;
   
   if ( config_->triggerResults() != "" )
      triggeranalysis_  = analysis_->triggerResults(config_->triggerResults());
   
   if ( config_->triggerObjectsDir() != "" )
   {
      // online jets
      if ( config_->triggerObjectsL1Jets() != "l1tJets")
         analysis_->addTree<TriggerObject> (config_->triggerObjectsL1Jets()  ,Form("%s/%s", config_->triggerObjectsDir().c_str(),config_->triggerObjectsL1Jets().c_str()));
      analysis_->addTree<TriggerObject> (config_->triggerObjectsCaloJets(),Form("%s/%s", config_->triggerObjectsDir().c_str(),config_->triggerObjectsCaloJets().c_str()));
      analysis_->addTree<TriggerObject> (config_->triggerObjectsPFJets()  ,Form("%s/%s", config_->triggerObjectsDir().c_str(),config_->triggerObjectsPFJets().c_str()));
      // online b jets
      analysis_->addTree<TriggerObject> (config_->triggerObjectsBJets(),Form("%s/%s", config_->triggerObjectsDir().c_str(),config_->triggerObjectsBJets().c_str()));
      // online muons
      analysis_->addTree<TriggerObject> (config_->triggerObjectsL1Muons(),Form("%s/%s",config_->triggerObjectsDir().c_str(),config_->triggerObjectsL1Muons().c_str()));
      analysis_->addTree<TriggerObject> (config_->triggerObjectsL3Muons(),Form("%s/%s",config_->triggerObjectsDir().c_str(),config_->triggerObjectsL3Muons().c_str()));
   }
   if ( config_ -> l1tJetsCollection() != "")
   {
      l1tjetsanalysis_ = ( analysis_ -> addTree<L1TJet> ("l1tJets",config_ -> l1tJetsCollection()) != nullptr );
   }
   if ( config_ -> l1tMuonsCollection() != "")
   {
      l1tmuonsanalysis_ = ( analysis_ -> addTree<L1TMuon> ("l1tMuons",config_ -> l1tMuonsCollection()) != nullptr );
   }

}

TriggerAnalyser::~TriggerAnalyser()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
// ------------ method called for each event  ------------


bool TriggerAnalyser::selectionTrigger() // Maybe not use this, use selectionHLT and selectionL1
{
   bool hlt = selectionHLT();
   bool l1  = selectionL1();
   
   /// Emulated triggers
   // L1 muon trigger
   bool l1muon = true;
   if ( config_->triggerEmulateL1Muons() != "" &&  config_->triggerEmulateL1MuonsNMin() > 0 && config_->triggerObjectsL1Muons() != "" )
   {
      int nmin = config_->triggerEmulateL1MuonsNMin();
      float ptmin = config_->triggerEmulateL1MuonsPtMin();
      float etamax = config_->triggerEmulateL1MuonsEtaMax();
      l1muon = selectionTriggerEmulated(l1,hlt,config_->triggerEmulateL1Muons(),nmin,ptmin,etamax);
   }
   // L3 muon trigger
   bool l3muon = true;
   if ( config_->triggerEmulateL3Muons() != "" &&  config_->triggerEmulateL3MuonsNMin() > 0  && config_->triggerObjectsL3Muons() != "" )
   {
      int nmin = config_->triggerEmulateL3MuonsNMin();
      float ptmin = config_->triggerEmulateL3MuonsPtMin();
      float etamax = config_->triggerEmulateL3MuonsEtaMax();
      l3muon = selectionTriggerEmulated(l1,hlt,config_->triggerEmulateL3Muons(),nmin,ptmin,etamax);
   }
      
   // L1 jet trigger
   bool l1jet = true;
   if ( config_->triggerEmulateL1Jets() != "" &&  config_->triggerEmulateL1JetsNMin() > 0 && config_->triggerObjectsL1Jets() != "" )
   {
      int nmin = config_->triggerEmulateL1JetsNMin();
      float ptmin = config_->triggerEmulateL1JetsPtMin();
      float etamax = config_->triggerEmulateL1JetsEtaMax();
      l1jet = selectionTriggerEmulated(l1,hlt,config_->triggerEmulateL1Jets(),nmin,ptmin,etamax);
   }
   
   // Calo jet trigger
   bool calojet = true;
   if ( config_->triggerEmulateCaloJets() != "" &&  config_->triggerEmulateCaloJetsNMin() > 0 && config_->triggerObjectsCaloJets() != "" )
   {
      int nmin = config_->triggerEmulateCaloJetsNMin();
      float ptmin = config_->triggerEmulateCaloJetsPtMin();
      float etamax = config_->triggerEmulateCaloJetsEtaMax();
      calojet = selectionTriggerEmulated(l1,hlt,config_->triggerEmulateCaloJets(),nmin,ptmin,etamax);
   }
   
   
   // PF jet trigger
   bool pfjet = true;
   if ( config_->triggerEmulatePFJets() != "" &&  config_->triggerEmulatePFJetsNMin() > 0 && config_->triggerObjectsPFJets() != "" )
   {
      int nmin = config_->triggerEmulatePFJetsNMin();
      float ptmin = config_->triggerEmulatePFJetsPtMin();
      float etamax = config_->triggerEmulatePFJetsEtaMax();
      pfjet = selectionTriggerEmulated(l1,hlt,config_->triggerEmulatePFJets(),nmin,ptmin,etamax);
   }
   
   
   bool emul = l1muon && l3muon && l1jet && calojet && pfjet;
   
   return (hlt && l1 && emul);
   
}

bool TriggerAnalyser::selectionHLT()
{
   if ( config_->hltPath_ == "" ) return true;
   
   ++cutflow_;
   if ( ! analysis_->triggerResult(config_->hltPath_) ) return false;
   
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
      h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,(config_->hltPath_).c_str());
   
   h1_["cutflow"] -> Fill(cutflow_,weight_);

   return true;
}

bool TriggerAnalyser::selectionL1()
{
   if ( config_->l1Seed_ == "" ) return true;
   
   ++cutflow_;
   if ( ! analysis_->triggerResult(config_->l1Seed_)  ) return false;
   
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
      h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,(config_->l1Seed_).c_str());
   
   h1_["cutflow"] -> Fill(cutflow_,weight_);

   return true;
}

bool TriggerAnalyser::selectionTriggerEmulated(const bool & l1, const bool & hlt, const std::string & name, const int & nmin, const float & ptmin, const float & etamax)
{
   
   ++cutflow_;
   
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
      h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,Form("Emulated: %s (n >= %d, pT >= %4.1f GeV, |eta| <= %4.1f)",name.c_str(),nmin,ptmin,etamax));
   
   
   if ( ! ( l1 && hlt ) ) return false;
   if ( ! triggerEmulated(name) ) return false;
   
   h1_["cutflow"] -> Fill(cutflow_,weight_);

   
   return true;
}




bool TriggerAnalyser::analysisWithTrigger()
{
   return triggeranalysis_;
}


std::vector< std::shared_ptr<TriggerObject> > TriggerAnalyser::triggerObjectsL1Jets()
{
   auto collection = analysis_->collection<TriggerObject>(config_->triggerObjectsL1Jets());
   std::vector< std::shared_ptr<TriggerObject> > objects;
   for ( int j = 0 ; j < collection->size() ; ++j )
      objects.push_back(std::make_shared<TriggerObject>(collection->at(j)));
   return objects;
}

std::vector< std::shared_ptr<TriggerObject> > TriggerAnalyser::triggerObjectsCaloJets()
{
   auto collection = analysis_->collection<TriggerObject>(config_->triggerObjectsCaloJets());
   std::vector< std::shared_ptr<TriggerObject> > objects;
   for ( int j = 0 ; j < collection->size() ; ++j )
      objects.push_back(std::make_shared<TriggerObject>(collection->at(j)));
   return objects;
}

std::vector< std::shared_ptr<TriggerObject> > TriggerAnalyser::triggerObjectsPFJets()
{
   auto collection = analysis_->collection<TriggerObject>(config_->triggerObjectsPFJets());
   std::vector< std::shared_ptr<TriggerObject> > objects;
   for ( int j = 0 ; j < collection->size() ; ++j )
      objects.push_back(std::make_shared<TriggerObject>(collection->at(j)));
   return objects;
}

std::vector< std::shared_ptr<TriggerObject> > TriggerAnalyser::triggerObjectsL1Muons()
{
   auto collection = analysis_->collection<TriggerObject>(config_->triggerObjectsL1Muons());
   std::vector< std::shared_ptr<TriggerObject> > objects;
   for ( int j = 0 ; j < collection->size() ; ++j )
      objects.push_back(std::make_shared<TriggerObject>(collection->at(j)));
   return objects;
}

std::vector< std::shared_ptr<TriggerObject> > TriggerAnalyser::triggerObjectsL3Muons()
{
   auto collection = analysis_->collection<TriggerObject>(config_->triggerObjectsL3Muons());
   std::vector< std::shared_ptr<TriggerObject> > objects;
   for ( int j = 0 ; j < collection->size() ; ++j )
      objects.push_back(std::make_shared<TriggerObject>(collection->at(j)));
   return objects;
}

bool  TriggerAnalyser::l1tJetsAnalysis() const
{
   return l1tjetsanalysis_;
}

bool  TriggerAnalyser::l1tMuonsAnalysis() const
{
   return l1tmuonsanalysis_;
}

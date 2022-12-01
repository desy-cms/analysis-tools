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
#include "Analysis/Tools/interface/Utils.h"

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
   l1tjetsanalysis_ = false;
   l1tmuonsanalysis_ = false;
   
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
   l1tjets_etabins_ = utilsL1TJetsEtaBins();
   l1tjets_phibins_ = utilsL1TJetsPhiBins();
   l1tmuons_etabins_ = utilsL1TMuonsEtaBins();
   l1tmuons_phibins_ = utilsL1TMuonsPhiBins();
   
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
   bool l1  = selectionL1();
   bool hlt = selectionHLT();
   
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

std::vector< std::shared_ptr<TriggerObject> > TriggerAnalyser::triggerObjectsBJets()
{
   auto collection = analysis_->collection<TriggerObject>(config_->triggerObjectsBJets());
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


void TriggerAnalyser::l1tjetHistograms(const std::string & label )
{
   this->output()->cd();
   if ( ! this->output()->FindObjectAny(label.c_str()) )
   {
      this->output()->mkdir(label.c_str());
      this->output()->cd(label.c_str());
   }
   else
   {
      if ( h1_.find(Form("n_l1tjet_%s"  , label.c_str())) != h1_.end() ) // the jet histograms already exist
      {
         return;
      }
   }
   
   n_hl1tjets_ = 13;
   
   h1_[Form("n_l1tjet_%s"  , label.c_str())]  = std::make_shared<TH1F>("n_l1tjet" , Form("n_l1tjet_%s" ,label.c_str()) ,12 , 0   , 12  );
   
   for ( int j = 0; j < n_hl1tjets_; ++j ) // loop over jets
   {
      // 1D histograms
      h1_[Form("pt_l1tjet%d_%s"  , j+1,label.c_str())]  = std::make_shared<TH1F>(Form("pt_l1tjet%d"  , j+1) , Form("pt_l1tjet%d_%s"  , j+1,label.c_str()) ,1000 , 0   , 1000  );
      h1_[Form("eta_l1tjet%d_%s" , j+1,label.c_str())]  = std::make_shared<TH1F>(Form("eta_l1tjet%d" , j+1) , Form("eta_l1tjet%d_%s" , j+1,label.c_str()) , 500 , -5, 5 );
      h1_[Form("phi_l1tjet%d_%s" , j+1,label.c_str())]  = std::make_shared<TH1F>(Form("phi_l1tjet%d" , j+1) , Form("phi_l1tjet%d_%s" , j+1,label.c_str()) , 360 , -180, 180 );
   }
   
   this->output()->cd();
}

void TriggerAnalyser::fillL1TJetHistograms(const std::string & label, std::vector<std::shared_ptr<L1TJet> > sel_l1tjets)
{
   this->output()->cd();
   this->output()->cd(label.c_str());
   
   int n = n_hl1tjets_;
   
   
   h1_[Form("n_l1tjet_%s"  , label.c_str())] -> Fill(sel_l1tjets.size(), weight_);
   
   for ( size_t j = 0; j < sel_l1tjets.size(); ++j )
   {
      int r = j+1;
      if ( r > n ) r = n;
      h1_[Form("pt_l1tjet%d_%s"   , r,label.c_str())] -> Fill(sel_l1tjets[j]->pt(),weight_);
      h1_[Form("eta_l1tjet%d_%s"  , r,label.c_str())] -> Fill(sel_l1tjets[j]->eta(),weight_);
      h1_[Form("phi_l1tjet%d_%s"  , r,label.c_str())] -> Fill(sel_l1tjets[j]->phi(),weight_);
   }
   this->output()->cd();
   
   cutflow(Form("*** Filling jets histograms - %s",label.c_str()));
   
}

void TriggerAnalyser::fillL1TJetHistograms(const std::string & label)
{
   fillL1TJetHistograms(label,selected_l1tjets_);
}

bool TriggerAnalyser::analysisWithL1TJets()
{
   l1tjets_.clear();
   selected_l1tjets_.clear();
   
   if ( ! l1tjetsanalysis_ ) return false;
   
   auto l1tjets = analysis_->collection<L1TJet>("l1tJets");
   for ( int j = 0 ; j < l1tjets->size() ; ++j ) 
   {
      l1tjets_.push_back(std::make_shared<L1TJet>(l1tjets->at(j)));
      selected_l1tjets_.push_back(l1tjets_.back());
   }
   
   return true;
}

std::vector< std::shared_ptr<L1TJet> > TriggerAnalyser::l1tJets()
{
   return l1tjets_;
}

std::vector< std::shared_ptr<L1TJet> > TriggerAnalyser::selectedL1TJets()
{
   return selected_l1tjets_;
}


bool TriggerAnalyser::selectionL1TJet(const float & ptmin, const float & etamax)
{
   bool isgood = true;
   std::string label = Form("L1TJet pt>=%4.2f, |eta|<=%4.2f", ptmin, etamax);;
   
   auto jets = selected_l1tjets_;
   selected_l1tjets_.clear();
   
   for ( auto & j : jets )
   {
      if ( j->pt() >= ptmin && fabs(j->eta()) <= etamax )
         selected_l1tjets_.push_back(j);
   }
   isgood = ( selected_l1tjets_.size() >= 1 );
   
   cutflow(label,isgood);
   return isgood;
   
}

bool TriggerAnalyser::selectionL1TDijet(const float & pt1min, const float & eta1max, const float & pt2min, const float & eta2max)
{
   bool isgood = true;
   std::string label;
   
   float pt1 = pt1min;
   float eta1 = eta1max;
   float pt2 = pt2min;
   float eta2 = eta2max;
   
   if ( pt2 < 0 ) pt2 = pt1;
   if ( eta2 < 0 ) eta2 = eta1;
   
   // Just to be flexible with inputs
   if ( pt1min < pt2min )
   {
      pt1 = pt2min;
      pt2 = pt1min;
      eta1 = eta2max;
      eta2 = eta1max;
   }
   
   if ( pt1 == pt2 && eta1 == eta2 ) // same thresholds
   {
      label = Form("L1TDijet pt>=%4.2f, |eta|<=%4.2f", pt1, eta1);
   }
   else
   {
      label = Form("L1TJet1: pt>=%4.2f, |eta|<=%4.2f; L1TJet2: pt>=%4.2f, |eta|<=%4.2f", pt1, eta1, pt2, eta2);
   }
   
   auto jets = selected_l1tjets_;
   selected_l1tjets_.clear();
   
   for ( auto & j : jets )
   {
      if ( j->pt() >= pt2 && fabs(j->eta()) <= eta2 )
         selected_l1tjets_.push_back(j);
   }
   isgood = ( selected_l1tjets_.size() >= 2 );
   isgood = ( isgood && ( selected_l1tjets_[0]->pt()>=pt1 && fabs(selected_l1tjets_[0]->eta()) <= eta1 ) );
   
   cutflow(label,isgood);
   return isgood;
   
}

bool TriggerAnalyser::selectionNL1TJets(const int & nmin)
{
   if ( nmin < 0 ) return true;

   std::string label = Form("L1TJets N >= %d",nmin);
   bool isgood = ( (int)selected_l1tjets_.size() >= nmin );
   
   cutflow(label,isgood);
   
   return isgood;
   
}

bool TriggerAnalyser::selectionL1TDijetDeta(const float & detamax)
{
   if ( detamax == 0 ) return true;

   std::string label = Form("L1TDijet delta_eta <= %4.2f",detamax);
   bool isgood = false;
   
   for ( auto & j1 : selected_l1tjets_ )
   {
      for ( auto & j2 : selected_l1tjets_ )
      {
         if ( j1 == j2 ) continue;
         
         if ( j1->deltaEta(*j2) <= detamax )   // a few percent difference wrt to the precise method below
//         if ( utilsL1TJetsDeta(j1->eta(),j2->eta()) <= detamax )
         {
            isgood = true;
            break;
         }
      }
      if ( isgood ) break;
   }
   cutflow(label,isgood);
   
   return isgood;
   
}

bool TriggerAnalyser::analysisWithL1TMuons()
{
   l1tmuons_.clear();
   selected_l1tmuons_.clear();
   
   if ( ! l1tmuonsanalysis_ ) return false;
   
   auto l1tmuons = analysis_->collection<L1TMuon>("l1tMuons");
   for ( int m = 0 ; m < l1tmuons->size() ; ++m ) 
   {
      l1tmuons_.push_back(std::make_shared<L1TMuon>(l1tmuons->at(m)));
      selected_l1tmuons_.push_back(l1tmuons_.back());
   }
   
   return true;
}

bool TriggerAnalyser::selectionNL1TMuons(const int & nmin)
{
   if ( nmin < 0 ) return true;

   std::string label = Form("L1TMuons N >= %d",nmin);
   bool isgood = ( (int)selected_l1tmuons_.size() >= nmin );
   
   cutflow(label,isgood);
   
   return isgood;
   
}

bool TriggerAnalyser::selectionL1TMuon(const float & ptmin, const float & etamax)
{
   bool isgood = true;
   std::string label = Form("L1TMuon pt>=%4.2f, |eta|<=%4.2f", ptmin, etamax);
   
   auto muons = selected_l1tmuons_;
   selected_l1tmuons_.clear();
   
   for ( auto & m : muons )
   {
      if ( m->pt() >= ptmin && fabs(m->eta()) <= etamax )
         selected_l1tmuons_.push_back(m);
   }
   isgood = ( selected_l1tmuons_.size() >= 1 );
   cutflow(label,isgood);
   
   return isgood;
   
}

bool TriggerAnalyser::selectionL1TMuonQuality(const int & qual)
{
   bool isgood = true;
   std::string label = Form("L1TMuon quality >= %d", qual);
   
   auto muons = selected_l1tmuons_;
   selected_l1tmuons_.clear();
   
   for ( auto & m : muons )
   {
      if ( m->hwQual() >= qual)
         selected_l1tmuons_.push_back(m);
   }
   isgood = ( selected_l1tmuons_.size() >= 1 );
   cutflow(label,isgood);
   
   return isgood;
   
}
bool TriggerAnalyser::selectionL1TMuonJet(const float & drmax)
{
   bool isgood = true;
   std::string label = Form("L1TMuonJet drmax <= %4.2f", drmax);
   
   auto muons = selected_l1tmuons_;
   auto jets = selected_l1tjets_;
   selected_l1tmuons_.clear();

   for ( auto & m : muons )
   {
      for ( auto & j : jets )
      {
         if ( m->deltaR(*j) <= drmax )
//         float dr = utilsL1TMuonJetDr(m->eta(),m->phi(),j->eta(),j->phi());
//         if ( dr <= drmax )
         {
            selected_l1tmuons_.push_back(m);
         }
      }
   }
   isgood = ( selected_l1tmuons_.size() >= 1 );
   cutflow(label,isgood);
   
   return isgood;
   
}


std::vector<float> TriggerAnalyser::l1tJetsEtaBins()
{
   return l1tjets_etabins_;
}
std::vector<float> TriggerAnalyser::l1tJetsPhiBins()
{
   return l1tjets_phibins_;
}
std::vector<float> TriggerAnalyser::l1tMuonsEtaBins()
{
   return l1tmuons_etabins_;
}
std::vector<float> TriggerAnalyser::l1tMuonsPhiBins()
{
   return l1tmuons_phibins_;
}

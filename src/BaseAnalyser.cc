// system include files
#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "TString.h"
// 
// user include files
#include "Analysis/Tools/interface/BaseAnalyser.h"

//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;

//
// constructors and destructor
//
BaseAnalyser::BaseAnalyser()
{
}


BaseAnalyser::BaseAnalyser(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();
   
   exe_ = std::string(argv[0]);
   
   // some inits
   scale_    = -1.;
   weight_   = 1.;
   xsection_ = -1.;
   genpartsanalysis_ = false;
   genjetsanalysis_  = false;
   primaryvtxanalysis_ = false;
   puweights_ok_ = false;
    
   // the heavy stuff
   config_   = std::make_shared<Config>(argc,argv);
//   analysis_ = std::make_shared<Analysis>(config_->ntuplesList(),config_->eventInfo());
   analysis_ = std::make_shared<Analysis>(config_);
   
   // output file
   if ( config_->outputRoot() != "" )
   {
      hout_= std::make_shared<TFile>(config_->outputRoot().c_str(),"recreate",Form("%s %s %s",argv[0],argv[1],argv[2]));
      hout_ -> cd();
   }
   
   seed_ = analysis_ ->seed(config_->seedFile());
   
   // Workflow
   h1_["cutflow"] = std::make_shared<TH1F>("workflow",Form("Workflow #%d: %s",config_->workflow(),config_->workflowTitle().c_str()), 100,0,100);
      
   
   isMC_ = analysis_->isMC();
   isData_ = !isMC_;
   
   // Some MC-only stuff
   if ( isMC_ )
   {
      // Cross sections
      if ( analysis_ -> crossSections(config_->crossSectionTree()) == 0 )
         xsection_ = analysis_->crossSection(config_->crossSectionType());

      // gen part analysis
      genpartsanalysis_  = ( analysis_->addTree<GenParticle> ("GenParticles",config_->genParticlesCollection()) != nullptr );
      // gen jets analysis
      genjetsanalysis_  = ( analysis_->addTree<GenJet> ("GenJets",config_->genJetsCollection()) != nullptr );
      
   }
   puweights_ok_ = (isMC_ || (isData_ && config_->pileupData() != "")) && (config_->pileupWeights() != "");

   // Pileup weights
   if ( puweights_ok_ )
   {
      puweights_ = analysis_->pileupWeights(config_->pileupWeights(),config_->pileupData(), isMC_);
      puw_label_ = basename(config_->pileupWeights());
   }
   else
   {
      if ( isMC_ )  puw_label_ = "*** missing *** assuming MC pileup weight = 1";
   }

   
   // primary vertex analysis
   primaryvtxanalysis_ = ( analysis_->addTree<Vertex> ("PrimaryVertex",config_->primaryVertexCollection()) != nullptr );
   
   // JSON for data   
   if( isData_ && config_->json() != "" ) analysis_->processJsonFile(config_->json());
   
   // btag efficiencies
   for ( int mb = 0; mb < 3; ++mb )
   {
      if ( config_->btagEfficiencies(mb+1) != "" )
      {
         TFile f(config_->btagEfficiencies(mb+1).c_str(),"old");
         auto list = f.GetListOfKeys();
         for ( int i = 0; i < list -> GetSize(); ++i)
         {
            TString item(list -> At(i) -> GetName());
            if ( ! item.BeginsWith("eff_")) continue;
            item.Remove(0,4);
            btageff_[mb][item.Data()] = std::shared_ptr<TGraphAsymmErrors>((TGraphAsymmErrors*)f.Get(("eff_"+item).Data()));
         }
         f.Close();
      }
   }
   this -> pileupHistogram();

   scale_correction_ = 1.;
   if ( config_->scaleFilename() != "" && config_->scaleParameter() != "" )
   {
      scale_data_ = readParameterDataCSVFile(config_->scaleFilename());
      if ( ! scale_data_.empty() )
         scale_correction_ = scale_data_[config_->scaleParameter()][0];
   }
}

BaseAnalyser::~BaseAnalyser()
{
   std::cout << std::endl;
   // get last bin
   int lastbin = 0;
   for ( int i = 1; i <= h1_["cutflow"] ->GetNbinsX(); ++i )
   {
      std::string label = std::string(h1_["cutflow"]->GetXaxis()->GetBinLabel(i));
      if ( label == "" )
      {
         lastbin = i-1;
         break;
      }
   }
   float fevts =  h1_["cutflow"] -> GetBinContent(lastbin);
   // overall scale
   float scale = 1.; 
   bool doscale = false;
   if ( scale_ > 0. )  // superseeds the scale from config
   {
      doscale = true;
      scale = scale_;
   }
   else
   {
      // scale from config
      if ( config_ -> scale() > 0. )
      {
         doscale = true;
         scale = config_ -> scale();
      }
   }
   if ( doscale )
   {
      if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(lastbin+1)) == "" )
      {
         h1_["cutflow"] -> GetXaxis()-> SetBinLabel(lastbin+1,Form("Number of events after scaling to %10.5f",scale));
      }
      h1_["cutflow"] -> Fill(lastbin,fevts*scale);
   }
    
   
   for ( auto h : h1_ )
   {
      if ( h.first == "cutflow" || h.first == "pileup" || h.first == "pileup_w" )    continue;
      if ( doscale ) h.second -> Scale(scale);
      bool is_empty =  ( h.second -> GetEntries() != 0 || h.second -> GetSumOfWeights() != 0 );
      if ( is_empty )
         continue;
      
   }
   for ( auto h : h2_ )
   {
      if ( doscale ) h.second -> Scale(scale);
      bool is_empty =  ( h.second -> GetEntries() != 0 || h.second -> GetSumOfWeights() != 0 );
      if ( is_empty ) continue;
   }
   
   if ( hout_ )
   {
//      std::cout << std::endl;
//      std::cout << "output root file: " << config_->outputRoot() << std::endl;
      hout_ -> cd();
      hout_ -> Write();
      hout_ -> Close();
      // print workflow using the Workflow macro
      try
      {
         system(Form("Workflow %s",config_->outputRoot().c_str()));
      }
      catch(...)
      {
         std::cout << "Problems with Workflow macro or the output file, no summary printed" << std::endl;
      }
   }   
   
   std::cout << std::endl;
   std::cout << exe_ << " finished!" << std::endl;
   printf("%s\n", std::string(100,'_').c_str());
   std::cout << std::endl;
   
   std::ofstream finished;
   finished.open("finished.txt");
   finished << exe_ << "\n";
   finished.close();
   
}


//
// member functions
//
// ------------ method called for each event  ------------

bool BaseAnalyser::event(const int & i) { return true; }
void BaseAnalyser::histograms(const std::string & s) { }

std::shared_ptr<Analysis> BaseAnalyser::analysis()
{
   return analysis_;
}

std::shared_ptr<Config> BaseAnalyser::config()
{
   return config_;
}

int BaseAnalyser::nEvents()
{
   int maxevt = config_->nEventsMax();
   if ( maxevt < 0 || maxevt > analysis_->size() ) maxevt = analysis_->size();
   return maxevt;
}

std::shared_ptr<TH1F> BaseAnalyser::histogram(const std::string & hname)
{
   if ( h1_.find(hname) == h1_.end() ) 
   {
      std::cout << "-e- BaseAnalyser::H1F(const string & hname) -> no histogram with hname = " << hname << std::endl;
      return nullptr;
   }
   
   return h1_[hname];
}

void BaseAnalyser::histogram(const std::string & hname, std::shared_ptr<TH1F> h1)
{
   h1_[hname] = h1;
}

std::map<std::string, std::shared_ptr<TH1F> > BaseAnalyser::histograms()
{
   return h1_;
}


int  BaseAnalyser::seed()
{
   return seed_;
}

int  BaseAnalyser::seed(const std::string & f)
{
   seed_ = analysis_ ->seed(f);
   return seed_;
}

void BaseAnalyser::seed(const int & seed)
{
   seed_ = seed;
}

void  BaseAnalyser::weight(const float & w)
{
   weight_ = w;
}

float  BaseAnalyser::weight()
{
   return weight_;
}

std::shared_ptr<TFile> BaseAnalyser::output()
{
   return hout_;
}

bool  BaseAnalyser::genParticlesAnalysis() const
{
   return genpartsanalysis_;
}

bool  BaseAnalyser::genJetsAnalysis() const
{
   return genjetsanalysis_;
}

bool  BaseAnalyser::primaryVertexAnalysis() const
{
   return primaryvtxanalysis_;
}

float BaseAnalyser::crossSection() const
{
   return xsection_;
}

std::shared_ptr<PileupWeight>  BaseAnalyser::pileupWeights() const
{
   return puweights_;
}

float BaseAnalyser::pileupWeight(const float & truepu, const int & var) const
{
   if ( ! puweights_ ) return 1.;   
   return puweights_->weight(truepu,var);
}

float BaseAnalyser::trueInteractions() const
{
   if ( ! config_->isMC() ) return -1;
    
   return float(analysis_->nTruePileup());
}

void BaseAnalyser::actionApplyPileupWeight(const int & var)
{
   if ( ! puweights_ok_ ) return;
      
   if ( puweights_ )
   {
      float truepu = -1;
      if ( isMC_ )
      {
         truepu = analysis_->nTruePileup();
         if ( truepu < 0 )
         {
            weight_ *= 1;
            std::cout << "-w- BaseAnalyser::actionApplyPileupWeight: pileup negative!? ";
            std::cout << "Please check! Assuming pileup weight = 1. " << std::endl;
         }
         else
         {
            weight_ *= this->pileupWeight(truepu,var);
         }
      }
      else
      {
         truepu = puweights_->getPileupFromData(analysis_->run(),analysis_->lumiSection());
         if ( truepu < 0 )
         {
            weight_ *= 1;
            // std::cout << "-w- BaseAnalyser::actionApplyPileupWeight: pileup negative ";
            // std::cout << "(run = " << analysis_->run() << ", ls = " << analysis_->lumiSection();
            // std::cout << ")!? Please check! Assuming pileup weight = 1. " << std::endl;
         }
         else
         {
            weight_ *= this->pileupWeight(truepu,var);
         }
      }
   }
   else
   {
      weight_ *= 1;
   }
   
   if ( var != 0 )
   {
      puw_label_ += Form(", syst: %+d sig",var);
   }
   
   cutflow(puw_label_);
   
   this -> fillPileupHistogram();

}

void BaseAnalyser::actionApplyPileupWeight()
{
   actionApplyPileupWeight(config_->pileupWeightSystematics());
}

void BaseAnalyser::pileupHistogram()
{
   this->output()->cd();
   if ( config_->min() > 0 && config_->max() > 0 )
   {
      h1_["pileup"] = std::make_shared<TH1F>("pileup" , "pileup" , config_->n() , config_->min() , config_->max() );
      h1_["pileup_w"] = std::make_shared<TH1F>("pileup_w" , "weighted pileup" , config_->n() , config_->min() , config_->max() );
   }
   else
   {
      h1_["pileup"] = std::make_shared<TH1F>("pileup" , "pileup" , 100 , 0 , 100 );
      h1_["pileup_w"] = std::make_shared<TH1F>("pileup_w" , "weighted pileup" , 100 , 0 , 100 );
   }
   
}
void BaseAnalyser::fillPileupHistogram()
{
   
   h1_["pileup"] -> Fill(analysis_->nTruePileup());
   h1_["pileup_w"] -> Fill(analysis_->nTruePileup(),this->pileupWeight(analysis_->nTruePileup(),0));
}

bool BaseAnalyser::pileupWeightsApplicable() const
{
   return puweights_ok_;
}

int BaseAnalyser::cutflow()
{
   return cutflow_;
}

void BaseAnalyser::cutflow(const int & c)
{
   cutflow_ = c;
}

void BaseAnalyser::cutflow(const std::string & label, const bool & ok)
{
   ++cutflow_;
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
   {
      h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,label.c_str());
   }
   if ( ok ) h1_["cutflow"] -> Fill(cutflow_,weight_);
   
}

void BaseAnalyser::scale(const float & scale)
{
   scale_ = scale;
}

std::string BaseAnalyser::basename(const std::string & name)
{
   std::string bn = "";
   std::vector<std::string> paths;
   if ( name != "" )
   {
      boost::split(paths, name, boost::is_any_of("/"));
      bn = paths.back();
   }
   return bn;

   
}

std::map<std::string, std::shared_ptr<TGraphAsymmErrors> > BaseAnalyser::btagEfficiencies(const int& model) const
{
   return btageff_[model-1];
}


void BaseAnalyser::generatorWeight()
{
   if ( ! config_->isMC() ) return;
   
   float weight = analysis_->genWeight();
   if ( config_->fullGenWeight() )
   {
      weight_ *= weight;
   }
   else
   {
      float sign =  (weight > 0) ? 1 : ((weight < 0) ? -1 : 0);
      weight_ *= sign;
   }
      
}

bool BaseAnalyser::triggerEmulation(const std::string & name, const int & nmin, const float & ptmin, const float & etamax, const std::string & newname)
{
   trg_emul_[newname] = true;
   std::vector<TriggerObject> new_objects;
      
   
   if ( name != "l1tJets" && name != "l1tMuons" )
   {
      std::shared_ptr< Collection<TriggerObject> > objects = analysis_->collection<TriggerObject>(name);
      
      for ( int i = 0 ; i < objects->size() ; ++i )
      {
         TriggerObject obj = objects->at(i);
         if ( obj.pt() >= ptmin && fabs(obj.eta()) <= etamax ) 
         {
            new_objects.push_back(obj);
         }
      }
      
   }
   
   if ( name == "l1tJets" )
   {
      std::shared_ptr< Collection<L1TJet> > l1tjets = analysis_->collection<L1TJet>(name);

      for ( int i = 0 ; i < l1tjets->size() ; ++i )
      {
         L1TJet l1tjet = l1tjets -> at(i);
         float pt = l1tjet.pt();
         float eta = l1tjet.eta();
         float phi = l1tjet.phi();
         float e = l1tjet.e();
         TriggerObject obj(pt,eta,phi,e);
         if ( obj.pt() >= ptmin && fabs(obj.eta()) <= etamax ) 
         {
            new_objects.push_back(obj);
         }
      }
      
   }
   if ( name == "l1tMuons" )
   {
      std::shared_ptr< Collection<L1TMuon> > l1tmuons = analysis_->collection<L1TMuon>(name);
      
      for ( int i = 0 ; i < l1tmuons->size() ; ++i )
      {
         L1TMuon l1tmuon = l1tmuons -> at(i);
         float pt = l1tmuon.pt();
         float eta = l1tmuon.eta();
         float phi = l1tmuon.phi();
         float e = l1tmuon.e();
         TriggerObject obj(pt,eta,phi,e);
         if ( obj.pt() >= ptmin && fabs(obj.eta()) <= etamax ) 
         {
            new_objects.push_back(obj);
         }
      }
      
   }
   trg_emul_[newname] = ( (int)new_objects.size() >= nmin );
      
   Collection<TriggerObject> new_collection(new_objects,newname);
   analysis_->addCollection(new_collection);
   
   return trg_emul_[newname];
   
   
   
}


std::map<std::string,bool> BaseAnalyser::triggersEmulated()
{
   return trg_emul_;
}

bool BaseAnalyser::triggerEmulated(const std::string & name)
{
   return trg_emul_[name];
}

void BaseAnalyser::actionApplyPrefiringWeight(const int & var)
{
   if ( ! config_->isMC() ) return;
   if ( ! config_->prefiringWeight() ) return;

   std::string prefw_label = "Prefiring weight";
 
   weight_ *= analysis_->prefiringWeight(var);
   
   if ( var != 0 )  prefw_label += Form(", syst: %+d sig",var);
   
   cutflow(prefw_label);
   
}
void BaseAnalyser::actionApplyPrefiringWeight()
{
   actionApplyPrefiringWeight(config_->prefiringWeightSystematics());
}

void BaseAnalyser::add1DHistogram(const std::string & label, const std::string & name, const std::string & title, const int & nbins, const float & min, const float & max)
{
   this->output()->cd();
   if (!this->output()->FindObjectAny(label.c_str()))
   {
      std::cout << "-warning- BaseAnalyser::add1DHistogram - directory with label " << label << " does not exist!" << std::endl;
      return;
   }
   this->output()->cd(label.c_str());

   std::string hist_tag = Form("%s_%s", name.c_str(), label.c_str());
   std::string hist_title = title;
   if ( title == "") hist_title = hist_tag;
   std::string hist_name = name;

   h1_[hist_tag] = std::make_shared<TH1F>(hist_name.c_str(), hist_title.c_str(), nbins, min, max);

   this->output()->cd();
}

void BaseAnalyser::fill1DHistogram(const std::string & label, const std::string & name, const float & value, const float & weight)
{
   std::string hist_tag = Form("%s_%s", name.c_str(), label.c_str());
   this->output()->cd(label.c_str());
   h1_[hist_tag]->Fill(value, weight);
   this->output()->cd();
}

void BaseAnalyser::actionApplyScaleCorrection(const std::string &title)
{
   if ( config_->scaleParameter() == "" || config_->scaleFilename() == "" || scale_data_.empty() ) return;
   weight_ *= scale_correction_;

   std::string label_type = "Scale correction";
   std::string label = Form("%6.4f for parameter %s (%s)", scale_correction_, config_->scaleParameter().c_str(), basename(config_->scaleFilename()).c_str());
   if ( title != "" ) 
   {
      label_type = "scale correction";
      label = Form("%s %s: %s", title.c_str(),label_type.c_str(), label.c_str());
   }
   cutflow(label);
}

/// get is muons analysis
bool BaseAnalyser::isMuonsAnalysis()
{
   return is_mouns_analysis_;
}
/// get is muons analysis
void BaseAnalyser::isMuonsAnalysis(const bool & is_muon_analysis)
{
   is_mouns_analysis_ = is_muon_analysis;
}

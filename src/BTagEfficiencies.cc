#include <iostream>
#include "Analysis/Tools/interface/BTagEfficiencies.h"
#include "TFile.h"

using namespace analysis::tools;
//
// constructors and destructor
//

BTagEfficiencies::BTagEfficiencies()
{
   filename_ = "";
   readonly_ = "";
   eff_ = 1.;
   etabins_ = {"0.0","2.5"};
}

BTagEfficiencies::BTagEfficiencies(const std::string & filename, const bool & readonly):
      filename_(filename),
      readonly_(readonly)
{
   eff_=1.;
   std::string opt = readonly ? "OLD" : "NEW";
   // file extension
   std::string filext = "";
   int dotpos = filename_.find_last_of(".");
   if ( dotpos > 0 ) filext = filename_.substr(dotpos+1);
   
   etabins_.clear();
   etabins_.push_back("0.0");
   
   if ( filext == "root" )
   {
      TFile * f = new TFile(filename.c_str(),"OLD");
      auto keys = f->GetListOfKeys();
      std::string prefix = "btag_eff_";
      for ( int i = 0; i < keys->GetSize(); ++i )
      {
         std::string name = keys->At(i)->GetName();
         if ( name.find_first_of(prefix) != 0 ) continue;
         std::string flavour = name;
         flavour.erase(0,prefix.length());
         flavour.erase(flavour.find_first_of("_"));
         std::string etabin = name;
         etabin.erase(0,etabin.find_last_of("_")+1);
         geff_[flavour][etabin] = (TGraphAsymmErrors*) f -> Get(name.c_str());
         auto etabin_i = etabin;
         auto etabin_f = etabin;
         etabin_i.erase(etabin_i.find_first_of("-"));
         etabin_f.erase(0,etabin_f.find_last_of("-")+1);
         if ( flavour == "b" )
            etabins_.push_back(etabin_f);
      }
      // open the root file
      
      // fill the Tgraphs geff_, e.g. geff_["udsg"]["0.0-0.5"] = ...
      //                                            ^         ^
      //                                       flavour  eta_range
      
      // maybe perform smoothing to a TSpline3?
      
      // close root file
      
      f->Close();
   }
}

BTagEfficiencies::~BTagEfficiencies()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


// methods
float BTagEfficiencies::efficiency(const std::string & flavour, const float & pt, const float & eta)
{
   // e.g. eff_ = TGraph / TSpline3 -> Eval() given the flavour pt and eta of the jet
   float eff = 1.;
   std::string etabin = "";
   
   for ( size_t i = 0; i < etabins_.size()-1; ++i )
   {
      if ( std::stof(etabins_[i]) <= fabs(eta) && fabs(eta) < std::stof(etabins_[i+1]) ) etabin = etabins_[i]+"-"+etabins_[i+1];
   }
   
   if ( etabin == "" )
   {
      std::cout << "*** warning *** eta = " << eta << " not in the valid range (|eta| < 2.2), returning eff=1" << std::endl;
      return eff;
   }
   
   eff = geff_[flavour][etabin]->Eval(pt);
   
   return eff;
}

float BTagEfficiencies::efficiency(const int & flavour, const float & pt, const float & eta)
{
   // e.g. eff_ = TGraph / TSpline3 -> Eval() given the flavour pt and eta of the jet
   return eff_;
}


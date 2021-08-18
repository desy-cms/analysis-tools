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
   // need to improve this to not hard code
   std::string etabin = "";
   if ( 0.0 <= fabs(eta) && fabs(eta) < 0.5 ) etabin = "0.0-0.5";
   if ( 0.5 <= fabs(eta) && fabs(eta) < 1.0 ) etabin = "0.5-1.0";
   if ( 1.0 <= fabs(eta) && fabs(eta) < 1.4 ) etabin = "1.0-1.4";
   if ( 1.4 <= fabs(eta) && fabs(eta) < 2.2 ) etabin = "1.4-2.2";
   
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


#include "Analysis/Tools/interface/BTagEfficiencies.h"

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
   // file extension
   std::string filext = "";
   int dotpos = filename_.find_last_of(".");
   if ( dotpos > 0 ) filext = filename_.substr(dotpos+1);
   
   if ( filext == "root" )
   {
      // open the root file
      
      // fill the Tgraphs geff_, e.g. geff_["udsg"]["0.0-0.5"] = ...
      //                                            ^         ^
      //                                       flavour  eta_range
      
      // maybe perform smoothing to a TSpline3?
      
      // close root file
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
   return eff_;
}

float BTagEfficiencies::efficiency(const int & flavour, const float & pt, const float & eta)
{
   // e.g. eff_ = TGraph / TSpline3 -> Eval() given the flavour pt and eta of the jet
   return eff_;
}


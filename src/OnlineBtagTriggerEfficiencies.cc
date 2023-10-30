#include <iostream>
#include <map>
#include "Analysis/Tools/interface/OnlineBtagTriggerEfficiencies.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1.h"


using namespace analysis::tools;
//
// constructors and destructor
//

BtagTriggerEfficiencies::BtagTriggerEfficiencies()
{
}


BtagTriggerEfficiencies::BtagTriggerEfficiencies(const std::string & filename)
{
   sf_name_ = "online_btag_scale_factors";
   sigma_variations_ = {0,+1,-1,+2,-2};
   sigma_variations_str_[0]  = "fit";
   sigma_variations_str_[+1] = "fit_plus_1sigma";
   sigma_variations_str_[-1] = "fit_minus_1sigma";
   sigma_variations_str_[+2] = "fit_plus_2sigma";
   sigma_variations_str_[-2] = "fit_minus_sigma";

   TFile f(filename.c_str(),"READ");
   for(auto & var : sigma_variations_)
   {
      std::string funct_name = Form("%s_%s",sf_name_.c_str(),sigma_variations_str_[var].c_str());
      scale_factors_[var] = (TGraph*)f.Get(funct_name.c_str());
   }
   f.Close();
   float bin_size = scale_factors_[0]->GetX()[1] - scale_factors_[0]->GetX()[0];
   pt_boundaries_.first = scale_factors_[0]->GetX()[0] - bin_size/2.;
   pt_boundaries_.second = scale_factors_[0]->GetX()[scale_factors_[0]->GetN()-1] + bin_size/2.;

}

BtagTriggerEfficiencies::~BtagTriggerEfficiencies()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


double BtagTriggerEfficiencies::findSF(const float & pT, const int & sigma)
{
   // scale factors must exist
   if ( !scale_factors_[sigma] || !scale_factors_[0] ) return 1.;

   // if pt within boundaries for any sigma variation
   if ( ( pT >= pt_boundaries_.first &&  pT <= pt_boundaries_.second ) ) 
      return scale_factors_[sigma]->Eval(pT);

   float pt_boundary = pt_boundaries_.first;
   if ( pT > pt_boundaries_.second ) pt_boundary = pt_boundaries_.second;

   // in case of uncertainty variations outside the boundaries
   float nominal_sf = scale_factors_[0]->Eval(pt_boundary);
   float uncertainty_sf = scale_factors_[sigma]->Eval(pt_boundary)-nominal_sf; // preserve sign; in case nominal uncertainty is zero
   return nominal_sf + 2.*uncertainty_sf;

}

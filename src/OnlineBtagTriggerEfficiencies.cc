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

//s   TH1::AddDirectory(kFALSE);

   f = new TFile(filename.c_str(),"READ");
   
   std::string var[5] = {"nominal","+1sigma","-1sigma","+2sigma","-2sigma"};
   std::string funct_name;
   

   for(int i = 0; i<3; i++)
   {
      for(int j = 0; j<5; j++)
      {
         funct_name = Form("SF_%s",var[j].c_str());
         functions[var[j].c_str()] = *((TF1*)f -> Get(funct_name.c_str()));
      }
   }

   functions["nominal"].GetRange(xmin,xmax);
   f->Close();
}

BtagTriggerEfficiencies::~BtagTriggerEfficiencies()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


double BtagTriggerEfficiencies::findSF(const float & pT, const int & sigma)
{
   TF1 sf_function, var_funct, nominal_funct;
   double sf = 1;
   
   std::string var = "";

   if (sigma == 0)
      var = "nominal";
   else
   {
      if (sigma > 0)
      var = Form("+%dsigma",int(fabs(sigma)));
      if (sigma < 0)
      var = Form("-%dsigma",int(fabs(sigma)));
   }
   
   nominal_funct = functions["nominal"];
      
   if (pT < xmax)
   {
      sf_function = functions[var.c_str()];
      sf = sf_function(pT);
   }
   else
   {
      var_funct = functions[var.c_str()];
      sf = 2* var_funct(xmax) - nominal_funct(xmax);// for syst take twice the variation from pT > upper range of the function
   }

   return sf;
}

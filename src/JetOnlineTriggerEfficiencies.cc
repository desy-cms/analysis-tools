#include <iostream>
#include <map>
#include "Analysis/Tools/interface/JetOnlineTriggerEfficiencies.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1.h"


using namespace analysis::tools;
//
// constructors and destructor
//

JetTriggerEfficiencies::JetTriggerEfficiencies()
{
}


JetTriggerEfficiencies::JetTriggerEfficiencies(const std::string & filename)
{

//s   TH1::AddDirectory(kFALSE);

   f = new TFile(filename.c_str(),"READ");
   
   std::string etabins[3] = {"0p0_to_1p0","1p0_to_1p4","1p4_to_2p2"};
   std::string var[5] = {"nominal","1s_up","1s_down","2s_up","2s_down"};
   std::string funct_name;
   

//                etabin                 var   function
//   std::map<std::string, std::map<std::string, TF1>> functions;

   for(int i = 0; i<3; i++)
   {
      for(int j = 0; j<5; j++)
      {
         funct_name = Form("jetOnlineTriggerScaleFactor_etabin_%s_%s",etabins[i].c_str(),var[j].c_str());
         functions[etabins[i].c_str()][var[j].c_str()] = *((TF1*)f -> Get(funct_name.c_str()));
       }
   }

   functions["0p0_to_1p0"]["nominal"].GetRange(xmin,xmax);
   f->Close();
}

JetTriggerEfficiencies::~JetTriggerEfficiencies()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


double JetTriggerEfficiencies::findSF(const float & eta, const float & pT, const int & sigma) 
{
   TF1 sf_function, var_funct, nominal_funct;
   double sf = 1;
   double etalimits[4] = {0.0, 1.0, 1.4, 2.2} ; 
   std::string etabin = "";
   std::string var = "";
   std::string etabins[3] = {"0p0_to_1p0","1p0_to_1p4","1p4_to_2p2"};

   if (fabs(eta) >= etalimits[2])// jets with eta > 2.2 will not pass offline selection in further steps of the analysis, included them in last eta bin
   etabin = etabins[2];
   else
   {
      for (int i = 0; i < 3; i++)
      {
         if(fabs(eta) >= etalimits[i] && fabs(eta) < etalimits[i+1] )
         etabin = etabins[i];
      }
   }

 
   if (sigma == 0)
      var = "nominal";
   else
   {
      if (sigma > 0)
      var = Form("%ds_up",int(fabs(sigma)));
      if (sigma < 0)
      var = Form("%ds_down",int(fabs(sigma)));
   }
   
   nominal_funct = functions[etabin.c_str()]["nominal"];
      
   if (pT < xmax)
   {
      sf_function = functions[etabin.c_str()][var.c_str()];
      sf = sf_function(pT);
   }
   else
   {
      var_funct = functions[etabin.c_str()][var.c_str()];
      sf = 2* var_funct(xmax) - nominal_funct(xmax);// for syst take twice the variation from pT > upper range of the function
   }

   return sf;
}


#include <iostream>
#include <map>
#include "Analysis/Tools/interface/MuonOnlineTriggerEfficiencies.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1.h"


using namespace analysis::tools;
//
// constructors and destructor
//

MuonTriggerEfficiencies::MuonTriggerEfficiencies()
{
}


MuonTriggerEfficiencies::MuonTriggerEfficiencies(const std::string & filename)
{

   f = new TFile(filename.c_str(),"READ");
   
   std::string var[5] = {"nominal","1s_up","1s_down","2s_up","2s_down"};
   std::string graph_name;
   

   for(int i = 0; i<5; i++)
   {
      graph_name = Form("muonOnlineTriggerScaleFactor_%s",var[i].c_str());
      graphs[var[i].c_str()] = *((TGraph*)f -> Get(graph_name.c_str()));
   }
   

   f->Close();
}

MuonTriggerEfficiencies::~MuonTriggerEfficiencies()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


double MuonTriggerEfficiencies::findSF(const float & pT, const int & sigma) 
{
   TGraph sf_graph, var_graph, nominal_graph;

   pTmax = 30; //max pT [GeV] for which the SF is defind

   double sf = 1;
   std::string var = "";

   if (sigma == 0)
      var = "nominal";
   else
   {
      if (sigma > 0)
      var = Form("%ds_up",int(fabs(sigma)));
      if (sigma < 0)
      var = Form("%ds_down",int(fabs(sigma)));
   }
   
   nominal_graph = graphs["nominal"];
      
   if (pT < pTmax)
   {
      sf_graph = graphs[var.c_str()];
      sf = sf_graph.Eval(pT);
   }
   else
   {
      var_graph = graphs[var.c_str()];
      sf = 2 * var_graph.Eval(pTmax) - nominal_graph.Eval(pTmax);// for syst take twice the variation from pT > upper range of the graph
   }

   return sf;
}


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


float MuonTriggerEfficiencies::findSF(const float & pT, const int & sigma) 
{
   TGraph sf_graph, var_graph, nominal_graph;

   std::vector <float> pTranges = {11.5, 12.5, 13.5, 18.5, 30.};
   int pTbin = 0;
   float sf = 1;
   std::string var = "";

   float pTmax = pTranges[pTranges.size() - 1 ]; 

   if (pT >= pTmax)// muons with pT > pTmax included in last bin
   pTbin = pTranges.size() - 2;

   else
   {
      for (int i = 0; i < int(pTranges.size()); i++)
      {
         if(pT > pTranges[i] && pT < pTranges[i+1] )
         pTbin = i;
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
   
   nominal_graph = graphs["nominal"];
   double x, y, xn, yn;

   if (pT <= pTmax)
   {
      sf_graph = graphs[var.c_str()];
      sf_graph.GetPoint(pTbin, x, y);
      sf = y;
      //sf = sf_graph.GetPoint(pTbin);
   }
   else
   {
      var_graph = graphs[var.c_str()];
      var_graph.GetPoint(pTbin, x, y);
      var_graph.GetPoint(pTbin, xn, yn);
      sf = 2 * y - yn;// for syst take twice the variation from pT > upper range of the graph
   }

   return sf;
}


/*
float MuonTriggerEfficiencies::findSF(const TGraph & graph, const int & sigma) 
{
   float
}*/

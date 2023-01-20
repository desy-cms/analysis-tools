#ifndef Analysis_Tools_MuonTriggerEfficiencies_h
#define Analysis_Tools_MuonTriggerEfficiencies_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class MuonOnlineTriggerEfficiencies
//

#include <map>
#include "TGraphAsymmErrors.h"
#include "TFile.h"

namespace analysis {
   namespace tools {
      class MuonTriggerEfficiencies
      {
         public:
            /// Constructor
            MuonTriggerEfficiencies();
            MuonTriggerEfficiencies(const std::string & filename);
            //MuonTriggerEfficiencies(const std::string & filename, const bool & readonly = true);
           ~MuonTriggerEfficiencies();
           
            std::map<std::string, TGraph> graphs;
            float findSF(const float & pT, const int & sigma);
            float pTmax=0;
           
            /// estimate the weight from a file
            //float efficiency(const std::string & flavour, const float & pt, const float & eta);
            //float efficiency(const int & flavour, const float & pt, const float & eta);
            
         protected:
            /// filename
            std::string filename_;
            /// read only
            bool readonly_;
         
            TFile * f;
     };
  } 
}          
            



#endif  // Analysis_Tools_BTagEfficiencies_h

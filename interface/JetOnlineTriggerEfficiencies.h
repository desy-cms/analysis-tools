#ifndef Analysis_Tools_JetTriggerEfficiencies_h
#define Analysis_Tools_JetTriggerEfficiencies_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class JetOnlineTriggerEfficiencies
//

#include <map>
#include "TGraphAsymmErrors.h"
#include "TFile.h"

namespace analysis {
   namespace tools {
      class JetTriggerEfficiencies
      {
         public:
            /// Constructor
            JetTriggerEfficiencies();
            JetTriggerEfficiencies(const std::string & filename);
            //JetTriggerEfficiencies(const std::string & filename, const bool & readonly = true);
           ~JetTriggerEfficiencies();
            double findSF(const float & eta, const float & pT, const std::string & filename, const int & sigma);
            TF1 findfunction(const std::string & filename, const int & sigma, const float & eta);
            double findSF(const float & eta, const float & pT, const int & sigma);
            TF1 findfunction(const int & sigma, const float & eta);
            std::map<std::string, std::map<std::string, TF1>> functions;
            double xmin=0,xmax=0;
           
            /// estimate the weight from a file
            //float efficiency(const std::string & flavour, const float & pt, const float & eta);
            //float efficiency(const int & flavour, const float & pt, const float & eta);
            
         protected:
            /// filename
            std::string filename_;
            /// read only
            bool readonly_;
            /// btag efficiency
            float eff_;
            /// btag efficiency graphs: mapping to flavour and eta range
            std::map<std::string, std::map<std::string, TGraphAsymmErrors *> > geff_;
            /// eta bins
            std::vector<std::string> etabins_;
            TFile * f;
     };
  } 
}          
            



#endif  // Analysis_Tools_BTagEfficiencies_h

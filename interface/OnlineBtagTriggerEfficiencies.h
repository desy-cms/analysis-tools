#ifndef Analysis_Tools_BtagTriggerEfficiencies_h
#define Analysis_Tools_BtagTriggerEfficiencies_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class BtagTriggerEfficiencies
//

#include <map>
#include "TGraphAsymmErrors.h"
#include "TFile.h"

namespace analysis {
   namespace tools {
      class BtagTriggerEfficiencies
      {
         public:
            /// Constructor
            BtagTriggerEfficiencies();
            BtagTriggerEfficiencies(const std::string & filename);
            //BtagTriggerEfficiencies(const std::string & filename, const bool & readonly = true);
           ~BtagTriggerEfficiencies();
            double findSF(const float & pT, const std::string & filename, const int & sigma);
            TF1 findfunction(const std::string & filename, const int & sigma);
            double findSF(const float & pT, const int & sigma);
            TF1 findfunction(const int & sigma);
            std::map<std::string, TF1> functions;
            double xmin=0,xmax=0;
           
            
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

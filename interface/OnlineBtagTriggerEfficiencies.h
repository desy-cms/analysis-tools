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
            // TF1 findfunction(const std::string & filename, const int & sigma);
            double findSF(const float & pT, const int & sigma);
            // TF1 findfunction(const int & sigma);
            // std::map<std::string, TF1> functions;
            double xmin=0,xmax=0;
           
            
         protected:
            /// filename
            std::string filename_;
            ///
            std::map<int, TGraph* > scale_factors_;
            /// read only
            bool readonly_;
            /// scale factor name
            std::string sf_name_;
            ///
            std::vector<int> sigma_variations_;
            ///
            std::map<int,std::string> sigma_variations_str_;
            ///
            std::pair<float,float> pt_boundaries_;
     };
  } 
}          
            



#endif  // Analysis_Tools_BTagEfficiencies_h

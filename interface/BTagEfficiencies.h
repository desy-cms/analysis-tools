#ifndef Analysis_Tools_BTagEfficiencies_h
#define Analysis_Tools_BTagEfficiencies_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      BTagEfficiencies
// 
/**\class xxx BTagEfficiencies.cc Analysis/Tools/src/BTagEfficiencies.cc

 Description: btag efficiencies from a given file

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh
//         Created:  Wed, 02 Jun 2021 18:45:00 CEST
//
//

#include <map>
#include "TGraphAsymmErrors.h"

namespace analysis {
   namespace tools {
      class BTagEfficiencies
      {
         public:
            /// Constructor
            BTagEfficiencies();
            BTagEfficiencies(const std::string & filename, const bool & readonly = true);
           ~BTagEfficiencies();
           
            /// estimate the weight from a file
            float efficiency(const std::string & flavour, const float & pt, const float & eta);
            float efficiency(const int & flavour, const float & pt, const float & eta);
            
         protected:
            /// filename
            std::string filename_;
            /// read only
            bool readonly_;
            /// btag efficiency
            float eff_;
            /// btag efficiency graphs: mapping to flavour and eta range
            std::map<std::string, std::map<std::string, TGraphAsymmErrors *> > geff_;
     };
  } 
}          
            



#endif  // Analysis_Tools_BTagEfficiencies_h

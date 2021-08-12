#include "Analysis/Tools/interface/Utils.h"
#include "TMath.h"
#include <iostream>

namespace analysis {
   namespace tools {

      
/* L1 trigger eta and phi bins (scales)
   http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=DN2014_029_v2.pdf
   https://globaltrigger.web.cern.ch/globaltrigger/upgrade/ugt/downloads
*/
                  
// L1TJets eta bins
      std::vector<float> utilsL1TJetsEtaBins()
      {
         int eta_bins = 115;
         float eta_max = 5.;
         float eta_step = 0.010875;
         float eta_0 = 0.;
         std::vector<float> eta(eta_bins*2+1);
         for ( int i = 0 ; i < eta_bins; ++i )
         {
            float min = eta_0+i*eta_step;
            float max = eta_0+(i+1)*eta_step;
            if ( max > eta_max ) max = eta_max;
            eta[i+eta_bins] = min;
            eta[eta_bins-1-i] = -max;
         }
         eta[eta_bins*2] = 5.;
         return eta;
         
      }

// L1TJets phi bins
      std::vector<float> utilsL1TJetsPhiBins()
      {
         float pi = TMath::Pi();
         int phi_bins = 144;
         std::vector<float> phi(phi_bins+1);
         float phi_step = 2*pi/phi_bins;
         for ( int i = 0 ; i < phi_bins; ++i )
         {
            float min = i*phi_step;
            phi[i] = min;
         }
         phi[phi_bins] = phi[phi_bins-1]+phi_step;
         return phi;
      }

// L1TMuons eta bins
      std::vector<float> utilsL1TMuonsEtaBins()
      {
         int eta_bins = 225;
         float eta_max = 2.45;
         float eta_step = 0.010875;
         float eta_0 = -eta_step/2.;
         std::vector<float> eta(eta_bins*2+1);
         for ( int i = 0 ; i < eta_bins; ++i )
         {
            float min = eta_0+i*eta_step;
            float max = eta_0+(i+1)*eta_step;
            if ( max > eta_max ) max = eta_max;
            eta[i+eta_bins] = min;
            eta[eta_bins-1-i] = -max;
         }
         eta[eta_bins*2] = eta_max;
         return eta;
         
      }

// L1TMuons phi bins
      std::vector<float> utilsL1TMuonsPhiBins()
      {
         float pi = TMath::Pi();
         int phi_bins = 576;
         float phi_step = 2*pi/phi_bins;
         std::vector<float> phi(phi_bins+1);
         for ( int i = 0 ; i < phi_bins; ++i )
         {
            float min = i*phi_step;
            phi[i] = min;
         }
         phi[phi_bins] = phi[phi_bins-1]+phi_step;
         return phi;
      }

   }
}

#include "Analysis/Tools/interface/Utils.h"
#include "TMath.h"
#include "TVector2.h"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
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
         float eta_step = 0.0435;
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

      float utilsL1TJetsDeta(const float & eta1, const float & eta2)
      {
         auto etabins = utilsL1TJetsEtaBins();
         float e1 = eta1;
         float e2 = eta2;
         if ( e1 > e2 )
         {
            e1 = eta2;
            e2 = eta1;
         }
         int eb1 = -1;
         int eb2 = -1;
         for ( size_t i = 0; i < etabins.size()-1; ++i )
         {
            if ( e1 >= etabins[i] && e1 < etabins[i+1] ) eb1 = i;
            if ( e2 >= etabins[i] && e2 < etabins[i+1] ) eb2 = i;
            if ( eb1 > 0  && eb2 > 0 ) break;
         }
         if ( eb1 == eb2 ) return 0.;
         
         return fabs(etabins[eb2]-etabins[eb1+1]);
         
      }
      float utilsL1TJetsDphi(const float & phi1, const float & phi2)
      {
         auto phibins = utilsL1TJetsPhiBins();
         float ff1 = TVector2::Phi_0_2pi(phi1);
         float ff2 = TVector2::Phi_0_2pi(phi2);
         float f1 = ff1;
         float f2 = ff2;
         if ( ff1 > ff2 )
         {
            f1 = ff2;
            f2 = ff1;
         }
         int fb1 = -1;
         int fb2 = -1;
         for ( size_t i = 0; i < phibins.size()-1; ++i )
         {
            if ( f1 >= phibins[i] && f1 < phibins[i+1] ) fb1 = i;
            if ( f2 >= phibins[i] && f2 < phibins[i+1] ) fb2 = i;
            if ( fb1 > 0  && fb2 > 0 ) break;
         }
         if ( fb1 == fb2 ) return 0.;
         
         float dphi = TVector2::Phi_mpi_pi(phibins[fb2]-phibins[fb1+1]);
         return fabs(dphi);
         
      }
      

      float utilsL1TMuonJetDeta(const float & m_eta, const float & j_eta)
      {
         auto j_etabins = utilsL1TJetsEtaBins();
         auto m_etabins = utilsL1TMuonsEtaBins();
         
         int m_b = -1;
         int j_b = -1;
         for ( size_t i = 0; i < m_etabins.size()-1; ++i )
         {
            if ( m_eta >= m_etabins[i] && m_eta < m_etabins[i+1] ) m_b = i;
            if ( m_b > 0 ) break;
         }
         for ( size_t i = 0; i < j_etabins.size()-1; ++i )
         {
            if ( j_eta > j_etabins[i] && j_eta < j_etabins[i+1] ) j_b = i;

            if ( j_b > 0 ) break;
         }
         
         float deta;
         
         if ( j_eta > m_eta ) deta = fabs(j_etabins[j_b]-m_etabins[m_b+1]);
         else                 deta = fabs(m_etabins[m_b]-j_etabins[j_b+1]);
         
         return deta;
         
      }
      float utilsL1TMuonJetDphi(const float & m_phi, const float & j_phi)
      {
         auto j_phibins = utilsL1TJetsPhiBins();
         auto m_phibins = utilsL1TMuonsPhiBins();
         
         float m_f = TVector2::Phi_0_2pi(m_phi);
         float j_f = TVector2::Phi_0_2pi(j_phi);

         int m_b = -1;
         int j_b = -1;
         for ( size_t i = 0; i < m_phibins.size()-1; ++i )
         {
            if ( m_f >= m_phibins[i] && m_f < m_phibins[i+1] ) m_b = i;
            if ( m_b > 0 ) break;
         }
         for ( size_t i = 0; i < j_phibins.size()-1; ++i )
         {
            if ( j_f >= j_phibins[i] && j_f < j_phibins[i+1] ) j_b = i;
            if ( j_b > 0 ) break;
         }
         
         
         float dphi;
         if ( j_phi > m_phi ) dphi = fabs(j_phibins[j_b]-m_phibins[m_b+1]);
         else                 dphi = fabs(m_phibins[m_b]-j_phibins[j_b+1]);
         
         // FIXME (print eta and deta above)
         
         std::cout << "PHI(MU)  " << m_phi << "   " << m_f << "    " << m_phibins[m_b] << "   " << m_phibins[m_b+1] << std::endl;
         std::cout << "PHI(JE)  " << j_phi << "   " << j_f << "    " << j_phibins[j_b] << "   " << j_phibins[j_b+1] << std::endl;
         std::cout << "DPHI     " << dphi << std::endl;
         
         return dphi;
         
      }
      
      
      float utilsL1TMuonJetDr(const float & eta1, const float & phi1, const float & eta2, const float & phi2)
      {
         float deta = utilsL1TMuonJetDeta(eta1,eta2);
         float dphi = utilsL1TMuonJetDphi(phi1,phi2);
         float dr = sqrt(deta*deta+dphi*dphi);
         
         return dr;
         
      }

      std::map<std::string, std::vector<float> > readParameterDataCSVFile(const std::string &filename)
      {
         std::map<std::string, std::vector<float> > csv_data;
         std::ifstream file(filename);
         if (file.is_open())
         {
            std::string line;
            std::getline(file, line);

            std::istringstream header(line);
            // std::string title;
            // while (std::getline(header, title, ','))
            // {
            //    // don't need to do anything
            // }

            while (std::getline(file, line))
            {
               std::istringstream row(line);
               std::string value;
               std::string parameter_name;
               int i = 0;
               while (std::getline(row>>std::ws, value, ','))
               {
                  if ( i == 0) // first column is the parameter name
                  {
                     parameter_name = value;
                     csv_data[parameter_name] = std::vector<float>();
                     ++i;
                     continue;
                  }
                  csv_data[parameter_name].push_back(std::stod(value));

                  // csv_data[it->first].push_back(std::stod(value));
                  // csv_data[csv_data.begin()->first].push_back(std::stod(value));
                  i++;
               }
            }
            file.close();
         }

         return csv_data;
      }
   }
}

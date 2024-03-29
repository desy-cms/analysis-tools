#ifndef Analysis_Tools_Utils_h
#define Analysis_Tools_Utils_h 1

#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp> // Include for boost::split
#include "Analysis/Tools/interface/JetResolution.h"

namespace analysis {
   namespace tools {

      struct FilterResults
      {
         int   total;
         int   filtered;
         float efficiency;
      };
      
      struct PDF
      {
         std::pair<int,int> id;
         std::pair<double,double> x;
      };

      template <typename T> 
      struct less_than_pt
      {
          inline bool operator() (const T& gp1, const T& gp2)
          {
              return (gp1.pt() < gp2.pt());
          }
      };
      
      template <typename T> 
      struct greater_than_pt
      {
          inline bool operator() (const T& gp1, const T& gp2)
          {
              return (gp1.pt() > gp2.pt());
          }
      };
      
      
      inline std::vector<std::string> getWords(const std::string & line)
      {
         std::vector<std::string> words;
         boost::split(words, line, boost::is_any_of(" "), boost::token_compress_on);
         return words;
      }
      
      struct jerDefinitions
      {
      };
      
      struct JetResolutionInfo
      {
//         JetResolutionInfo(JME::JetResolution r, JME::JetResolutionScaleFactor s) :
//               resolution(r), scalefactor(s)
//         {}
         JME::JetResolution resolution;
         JME::JetResolutionScaleFactor scalefactor;
      };      

      struct ScaleFactors
      {
         float nominal;
         float up;
         float down;
      };
      
      std::vector<float> utilsL1TJetsEtaBins();
      std::vector<float> utilsL1TJetsPhiBins();
      std::vector<float> utilsL1TMuonsEtaBins();
      std::vector<float> utilsL1TMuonsPhiBins();
      
      float utilsL1TJetsDeta(const float & eta1, const float & eta2);
      float utilsL1TJetsDphi(const float & phi1, const float & phi2);
      float utilsL1TMuonJetDeta(const float & m_eta, const float & j_eta);
      float utilsL1TMuonJetDphi(const float & phi1, const float & phi2);
      float utilsL1TMuonJetDr(const float & eta1, const float & phi1, const float & eta2, const float & phi2);

      std::map<std::string, std::vector<float> > readParameterDataCSVFile(const std::string & filename);
      }
}




#endif  // Analysis_Tools_Utils_h

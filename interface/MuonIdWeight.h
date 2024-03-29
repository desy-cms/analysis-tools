#ifndef Analysis_Tools_MuonIdWeight_h
#define Analysis_Tools_MuonIdWeight_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      xxx
// 
/**\class xxx MuonIdWeight.cc Analysis/Tools/src/MuonIdWeight.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Antonio Vagnerini
//
//

// system include files
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <iostream>
// 
// user include files

#include "TH2.h"


//
// class declaration
//

namespace analysis {
   namespace tools {

      class MuonIdWeight {
         public:
            /// constructors
            MuonIdWeight();
            MuonIdWeight(const std::vector<std::string> & fnames );
            //MuonIdWeight(const std::string &); //TO DO, read vector of files
            /// desctructor
           ~MuonIdWeight();
           
            // ----------member data ---------------------------
         protected:
            std::map<int,std::shared_ptr<TH2D> > histos_;
                        
         private:
               
         public:
            float weight(const float & pt, const float & eta, const int & var = 0);
            float findSF(const float & pt, const float & eta , const int & var = 0);
            float sf_pTmin = 100000, sf_pTmax = -1;
            float fabs_etamax = -1;
            std::string findfile(const float & pT);
            std::vector <std::pair<std::pair<double,double>, std::string> > pTranges_file; //pTmin, pTmax, filename
            std::map<std::string, TH2F> filename_sf_map;
            std::map<std::string, TH2F> filename_unc_map;
            //std::vector<std::pair<std::string, TH2F>> filename_sf_pairs; //filename, histograme
            //std::vector<std::pair<std::string, TH2F>> filename_unc_pairs; //filename, histogram
            //std::vector<TH2F>  sf_hists;
            //std::shared_ptr<TH2D> histogram(const int & var = 0);
            
      };
   }
}

#endif  // Analysis_Tools_MuonIdWeight_h

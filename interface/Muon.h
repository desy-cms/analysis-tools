#ifndef Analysis_Tools_Muon_h
#define Analysis_Tools_Muon_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      Muon
// 
/**\class xxx Muon.cc Analysis/Tools/src/Muon.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh Bastos Rangel
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
//
//

// system include files
#include <memory>
// 
// user include files
#include "Analysis/Tools/interface/Candidate.h"
//
// class declaration
//

namespace analysis {
   namespace tools {

      class Muon : public Candidate {
         public:
            /// default constructor
            Muon();
            /// constructor from 4-momentum information
            Muon(const float & pt, const float & eta, const float & phi, const float & e, const float & q);
            /// destructor
           ~Muon();
           // Gets
           
           /// is particle flow muon
           bool  isPFMuon()        const ;
           /// is global muon
           bool  isGlobalMuon()    const ; 
           /// is tracker muon
           bool  isTrackerMuon()   const ; 
           /// is loose id muon
           bool  isLooseMuon()     const ; 
           /// is medium id muon
           bool  isMediumMuon()    const ;
           /// is tight id muon
           bool  isTightMuon()     const ;
           /// muon id
           bool  id(const std::string & ) const;
           

           float validFraction()          const ;
           float segmentCompatibility()   const ; 
           float trkKink()                const ; 
           float chi2LocalPos()           const ;

           float normChi2()               const ;

           // Sets
           /// set flag for particle flow muon
           void isPFMuon     (const bool & isPFMuon)      ;    
           /// set flag for global muon
           void isGlobalMuon (const bool & isGlobalMuon)  ; 
           /// set flag for tracker muon
           void isTrackerMuon(const bool & isTrackerMuon) ;  
           /// set flag for loose muon id
           void isLooseMuon  (const bool & isLooseMuon)   ; 
           /// set flag for medium muon id
           void isMediumMuon (const bool & isMediumMuon)  ;
           /// set flag for tight muon id
           void isTightMuon  (const bool & isTightMuon)  ;

           void validFraction( const float & validFraction) ;     
           void segmentCompatibility(const float & segmentCompatibility) ; 
           void trkKink      ( const float & trkKink     )  ;
           void chi2LocalPos ( const float & chi2LocalPos)  ;                     

           void normChi2     ( const float & normChi2 )    ;           

       protected:
       // ----------member data ---------------------------
           bool isPFMuon_;             
           bool isGlobalMuon_;         
           bool isTrackerMuon_;        
           bool isLooseMuon_;          
           bool isMediumMuon_;         
           bool isTightMuon_;

           float validFraction_;        
           float segmentCompatibility_; 
           float trkKink_;              
           float chi2LocalPos_;         

           float normChi2_;
           /// muon id

       private:
            // ----------member data ---------------------------
            
            // 
      };
   }
}

#endif  // Analysis_Tools_Muon_h

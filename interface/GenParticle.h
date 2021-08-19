#ifndef Analysis_Tools_GenParticle_h
#define Analysis_Tools_GenParticle_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      GenParticle
// 
/**\class xxx GenParticle.cc Analysis/Tools/src/GenParticle.cc

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

      class GenParticle : public Candidate {
         public:
            /// default contructor
            GenParticle();
            /// constructor from particle kinematics
            GenParticle(const float & pt, const float & eta, const float & phi, const float & e, const float & q);
            /// destructor
           ~GenParticle();
           
            /// set the pdgId of the particle 
            void pdgId(const int & pdgId);
            /// set the generator status of the particle 
            void status(const int & status);
            /// set whether the particle is a Higgs daughter 
            void higgsDaughter(const bool & higgs_dau);
            /// set the index of the particle from the generator particle list 
            void index(const int &);
            /// set the indices of the particle mother (for some reason there are two, usually they are the same)
            void mother(const int &, const int &);
            /// set the indices of the particle daughter
            void daughter(const int &, const int &);
            
            /// returns the pdgID of the particle
            int pdgId();
            /// returns the generator status of the particle
            int status();
            /// returns whether the particle is a Higgs daughter
            bool higgsDaughter();
            /// returns the index of the particle in the generator particle list
            int index();
            /// returns the index of the 2 mothers, where n=1 for first mother and n=2 for second mother
            int mother(const int & n = 1);
            /// returns the index of the 2 daughters, where n=1 for first daughter and n=2 for second daughter
            int daughter(const int & n = 1);
      
         private:
            // ----------member data ---------------------------
               
            /// the pdgId of the particle 
            int   pdgid_;
            /// the generator status of the particle 
            int   status_;
            /// is the particle a Higgs daughter
            bool  higgs_dau_;
            /// index of the particle in the generator particle list
            int   indx_;
            /// array of mothers
            int   mo_[2];
            /// array of daughters
            int   da_[2];
            // 
      };
      // ===============================================
      // INLINE IMPLEMENTATIONS
         
      inline int   GenParticle::pdgId()                    { return pdgid_    ; }                   
      inline int   GenParticle::status()                   { return status_   ; }                   
      inline bool  GenParticle::higgsDaughter()            { return higgs_dau_; }                   
      inline int   GenParticle::index()                    { return indx_     ; }
      inline int   GenParticle::mother(const int & n)
      {
         if ( n < 1 || n > 2 ) { std::cout << "*w* GenParticle has two possible mothers"   << std::endl; return -1; }
         return mo_[n-1]    ;
      }
      inline int   GenParticle::daughter(const int & n)
      {
         if ( n < 1 || n > 2 ) { std::cout << "*w* GenParticle has two possible daughters" << std::endl; return -1; }
         return da_[n-1]    ;
      }

      inline void GenParticle::pdgId  (const int & pdgId)             { pdgid_  = pdgId  ; } 
      inline void GenParticle::status (const int & status)            { status_ = status ; } 
      inline void GenParticle::higgsDaughter (const bool & higgs_dau) { higgs_dau_ = higgs_dau; } 
      inline void GenParticle::index(const int & indx)                { indx_ = indx    ; } 
      inline void GenParticle::mother(const int & n, const int & mo)
      {
         if ( n < 1 || n > 2 ) { std::cout << "*w* GenParticle has two possible mothers" << std::endl; }
         else  mo_[n-1] = mo  ; 
      }
      inline void GenParticle::daughter(const int & n, const int & da)
      {
         if ( n < 1 || n > 2 ) { std::cout << "*w* GenParticle has two possible daughters" << std::endl; }
         else  da_[n-1] = da  ; 
      }
         
   }
}

#endif  // Analysis_Tools_GenParticle_h

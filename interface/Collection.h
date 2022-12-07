#ifndef Analysis_Tools_Collection_h
#define Analysis_Tools_Collection_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      Collection
// 
/**\class xxx Collection.h Analysis/Tools/interface/Collection.h

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
#include <vector>
// 
// user include files
#include "Analysis/Tools/interface/Candidate.h"
#include "Analysis/Tools/interface/Jet.h"
#include "Analysis/Tools/interface/TriggerObject.h"
#include "Analysis/Tools/interface/GenParticle.h"
#include "Analysis/Tools/interface/GenJet.h"
#include "Analysis/Tools/interface/Vertex.h"

//
// class declaration
//

namespace analysis {
   namespace tools {
      
      template <class Object>
      class Collection {
         
         typedef std::vector<Object> Objects;

         public:
            /// default constructor
            Collection();
            /// constructor
            Collection(const Objects & objects, const std::string & name_ = "");
            /// destructor
           ~Collection();

            /// name of the collection
            std::string name() const;
            /// size of the collection
            int size();
            /// set the size of the collection
            void setSize(const int & size);
            /// retrieve the object at the specified position/index
            Object & at(const int & index);
            /// add object to the collection
            void add(const Object & object);
            
            /// deltaR match a vector of candidates to the objects of the collection, must give a name for the candidates
            void matchTo( const std::vector<Candidate>* vectorcandidates, const std::string & name , const float & deltaR = 0.5 );
            /// deltaR match a collections of candidates to the objects of the collection
            void matchTo( const Collection<Candidate> & collection, const float & deltaR = 0.5 );
            /// deltaR and deltaPT match a collection of candidates to the objects of the collection
            void matchTo( const Collection<Candidate> & collection, const float & delta_pT, const float & deltaR);
            /// deltaR and deltaPT match a collection of jets to the objects of the collection
            void matchTo( const Collection<Jet> & collection, const float & delta_pT, const float & deltaR);
            /// deltaR match a collection of trigger objects to the objects of the collection
            void matchTo( const Collection<TriggerObject> & collection, const float & deltaR = 0.5 );
            /// deltaR match a collection of trigger objects (shared_ptr) to the objects of the collection
            void matchTo( const std::shared_ptr<Collection<TriggerObject> > collection, const float & deltaR = 0.5 );

            /// returns a vector of the collection objects shared pointers
            std::vector< std::shared_ptr<Object> > vector();
            /// returns a vector of candidates objects
            std::vector<Candidate>* vectorCandidates() const;
            
            /// detalR match of partons to the collection objects
            void associatePartons(const std::shared_ptr<Collection<GenParticle> > & , const float & deltaR = 0.4, const float & ptMin = 1., const bool & pythia8 = true);

            void addGenJets(const std::shared_ptr<Collection<GenJet> > &);
            void btagAlgo(const std::string &);           
            // ----------member data ---------------------------
         protected:
               
         private:
            Objects objects_;
            mutable std::vector<Candidate> candidates_; // maybe not the best idea but need to make code work
            int size_;
            std::string name_;

      };
      // Implementation of non-specialised functions for any template object must be in the header file
      // ===============================================
      // Gets
      template <class Object> int         Collection<Object>::size()                { return size_; }
      template <class Object> Object  &   Collection<Object>::at(const int & index) { return objects_.at(index); }
      template <class Object> std::string Collection<Object>::name() const          { return name_; }
      // Sets
      template <class Object> void   Collection<Object>::add(const Object & object) { objects_.push_back(object); ++size_;  }
      template <class Object> void   Collection<Object>::setSize(const int & size) { size_ = size; }


      // member functions specialization declarations - needed to be declared in the same namespace as the class
      template <> Collection<Vertex>::Collection(const Objects & objects, const std::string & name);
      template <> std::vector<Candidate> * Collection<Vertex>::vectorCandidates() const;
      template <> void Collection<Jet>::matchTo( const Collection<Jet> & collection, const float & deltaR, const float & delta_pt);
      template <> void Collection<Jet>::associatePartons( const std::shared_ptr<Collection<GenParticle> > & particles, const float & deltaR, const float & ptMin, const bool & pythia8  );
      template <> void Collection<Jet>::btagAlgo( const std::string & algo );
      template <> void Collection<Jet>::addGenJets( const std::shared_ptr<Collection<GenJet> > & cands );
      /// definitions for Vertex candidates, which has no matching purposes
      template <> void Collection<Vertex>::matchTo( const std::vector<Candidate>* vectorcandidates, const std::string & name, const float & deltaR );
      template <> void Collection<Vertex>::matchTo( const Collection<Candidate> & collection, const float & deltaR );
      template <> void Collection<Vertex>::matchTo( const Collection<Candidate> & collection, const float & deltaR, const float & delta_pt);
      template <> void Collection<Vertex>::matchTo( const Collection<TriggerObject> & collection, const float & deltaR );
      template <> void Collection<Vertex>::matchTo( const std::shared_ptr<Collection<TriggerObject> > collection, const float & deltaR );


   }
}

#endif  // Analysis_Tools_Tree_h

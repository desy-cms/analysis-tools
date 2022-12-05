// system include files
#include <iostream>
#include <set>
//
#include "TRandom2.h"
// user include files
#include "Analysis/Tools/interface/Candidate.h"
#include "Analysis/Tools/interface/Jet.h"
#include "Analysis/Tools/interface/MET.h"
#include "Analysis/Tools/interface/Muon.h"
#include "Analysis/Tools/interface/Vertex.h"
#include "Analysis/Tools/interface/GenParticle.h"
#include "Analysis/Tools/interface/GenJet.h"
#include "Analysis/Tools/interface/JetTag.h"
#include "Analysis/Tools/interface/L1TMuon.h"
#include "Analysis/Tools/interface/L1TJet.h"
#include "Analysis/Tools/interface/RecoMuon.h"
#include "Analysis/Tools/interface/RecoTrack.h"

#include "Analysis/Tools/interface/Collection.h"

//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;

//
// constructors and destructor
//
template <class Object>
Collection<Object>::Collection()
{
	size_ = 0;
	candidates_.clear();
}
template <class Object>
Collection<Object>::Collection(const Objects & objects, const std::string & name)
{
   objects_ = objects;
   size_ = (int) objects_.size();
   name_ = name;
   candidates_.clear();
   for ( int i = 0; i < size_ ; ++i ) candidates_.push_back(objects_[i]);
}
template <>
Collection<Vertex>::Collection(const Objects & objects, const std::string & name)
{
   objects_ = objects;
   size_ = (int) objects_.size();
   name_ = name;
}
template <class Object>
Collection<Object>::~Collection()
{

}


//
// member functions
//

template <class Object>
std::vector< std::shared_ptr<Object> >  Collection<Object>::vector()
{
   std::vector< std::shared_ptr<Object> > objects;
   for ( auto & obj : objects_ )
   {
      std::shared_ptr<Object> sp_obj = std::shared_ptr<Object> ( new Object(obj));
      objects.push_back(sp_obj);
   }
   return objects;
}
template <class Object>
void Collection<Object>::btagAlgo(const std::string & algo  )
{
}
template <>
void Collection<Jet>::btagAlgo(const std::string & algo  )
{
   for ( auto & jet : objects_ )
      jet.btagAlgo(algo);
}
template <class Object>
void Collection<Object>::addGenJets(const std::shared_ptr<Collection<GenJet> > & cands)
{
}
template <>
void Collection<Jet>::addGenJets(const std::shared_ptr<Collection<GenJet> > & cands)
{
   if ( objects_.size() < 1 ) return;
   std::vector< std::shared_ptr<GenJet> > vec = cands->vector();
   for ( auto & jet : objects_ )
      jet.genJets(vec);
   
}
template <class Object>
void Collection<Object>::associatePartons(const std::shared_ptr<Collection<GenParticle> > & particles, const float & deltaR, const float & ptMin, const bool & pythia8  )
{
}
template <>
void Collection<Jet>::associatePartons(const std::shared_ptr<Collection<GenParticle> > & particles, const float & deltaR, const float & ptMin, const bool & pythia8  )
{
   if ( objects_.size() < 1 ) return;
   auto vec = particles->vector();
   
   for ( auto & jet : objects_ )
      jet.associatePartons(vec,deltaR,ptMin,pythia8);
   
   // resolving ambiguities
   // if a parton belongs to more than one jet, than remove it 
   // from the jet with largest deltaR wrt the parton.
   // std::map<Jet *, std::vector<int> > removeFromJet;
   std::map<Jet *, std::set<int> > removeFromJet;
   
   for ( size_t j1 = 0 ; j1 < objects_.size()-1 ; ++j1 )
   {
      Jet * jet1 = &(objects_.at(j1));
      auto partons1 = jet1->partons();
      for ( auto it1 = partons1.begin(); it1 != partons1.end(); it1++ )
      {
         auto parton1 = *it1;
         auto p1 = std::distance(partons1.begin(), it1 );
         for ( size_t j2 = j1+1 ; j2 < objects_.size() ; ++j2 )
         {
            Jet * jet2 = &(objects_.at(j2));
            auto partons2 = jet2->partons();
            for ( auto it2 = partons2.begin(); it2 != partons2.end(); ++it2 )
            {
               auto parton2 = *it2;
               int p2 = std::distance(partons2.begin(), it2 );
               if ( parton1 == parton2 )
               {
                  float dR1 = jet1->p4().DeltaR(parton1->p4());
                  float dR2 = jet2->p4().DeltaR(parton2->p4());
                  if ( dR1 < dR2 ) { removeFromJet[jet2].insert(p2) ;}
                  else             { removeFromJet[jet1].insert(p1) ;}
               }
            }
         }
         
      }
   }
   
   // Remove ambiguous partons...
   for ( auto & jet : removeFromJet )
   {
      int counts = 0;
      for ( auto & p : jet.second )
      {
         jet.first -> removeParton(p-counts);
         ++counts;
      }
   }
   
}
template <class Object>
void Collection<Object>::matchTo( const std::vector<Candidate>* vectorcandidates, const std::string & name, const float & deltaR )
{
   for ( auto & obj : objects_ )
      obj.matchTo(vectorcandidates,name,deltaR);
}
template <class Object>
void Collection<Object>::matchTo( const Collection<Candidate> & collection, const float & deltaR )
{
   for ( auto & obj : objects_ )
      obj.matchTo(collection.vectorCandidates(),collection.name(), deltaR);
}
template <class Object>
void Collection<Object>::matchTo( const Collection<Candidate> & collection, const float & delta_pT, const float & deltaR)
{
	for (auto & obj : objects_)
   {
		obj.matchTo(collection.vectorCandidates(),collection.name(),delta_pT,deltaR);
	}
}
template <>
void Collection<Jet>::matchTo( const Collection<Jet> & collection, const float & delta_pT, const float & deltaR)
{
	for (auto & obj : objects_)
   {
		obj.matchTo(collection.vectorCandidates(),collection.name(),delta_pT*obj.jerPtResolution()*obj.pt(),deltaR);
	}
}
template <class Object>
void Collection<Object>::matchTo( const Collection<TriggerObject> & collection, const float & deltaR )
{
   for ( auto & obj : objects_ )
      obj.matchTo(collection.vectorCandidates(),collection.name(),deltaR);
}
template <class Object>
void Collection<Object>::matchTo( const std::shared_ptr<Collection<TriggerObject> > collection, const float & deltaR )
{
   this->matchTo(*collection, deltaR);
}
template <class Object>
std::vector<Candidate>* Collection<Object>::vectorCandidates() const
{
   return &candidates_;
}
// try to find how the enable_if works to avoid this specialization
// typename std::enable_if<std::is_base_of<Candidate, Object>::value, std::vector<Candidate>* >::type
// std::is_base_of<Foo, Bar>::value
template <>
std::vector<Candidate> * Collection<Vertex>::vectorCandidates() const
{
   return nullptr;
}
// ------------ method called for each event  ------------


// ------------ declare the template classes  ------------
template class Collection<Candidate>;
template class Collection<Jet>;
template class Collection<MET>;
template class Collection<Muon>;
template class Collection<Vertex>;
template class Collection<TriggerObject>;
template class Collection<GenParticle>;
template class Collection<GenJet>;
template class Collection<JetTag>;
template class Collection<L1TMuon>;
template class Collection<L1TJet>;
template class Collection<RecoMuon>;
template class Collection<RecoTrack>;

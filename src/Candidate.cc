// system include files
//
#include <iostream>
#include <iomanip>
// user include files
#include "Analysis/Tools/interface/Candidate.h"


//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;

//
// constructors and destructor
//
Candidate::Candidate()
{
   q_ = 0;
   p4_.SetPtEtaPhiE(0.,0.,0.,0.);
   components_.clear();
}

Candidate::Candidate(const float & pt, const float & eta, const float & phi, const float & e, const float & q)
{
   q_ = q;
   p4_.SetPtEtaPhiE(pt,eta,phi,e);
   components_.clear();
}

Candidate::Candidate(const float & px, const float & py, const float & pz)
{
   q_ = 0;
   p4_.SetXYZM(px,py,pz,0.);
   components_.clear();
}

Candidate::Candidate(const float & px, const float & py, const float & pz, const float & q)
{
   q_ = q;
   p4_.SetXYZM(px,py,pz,0.);
   components_.clear();
}

Candidate::Candidate(const TLorentzVector & p4, const float & q)
{
   q_ = q;
   p4_ = p4;
   components_.clear();
}

Candidate::~Candidate()
{
}

// operators
Candidate Candidate::operator + (const Candidate & cand) const
{
    auto p4 = p4_+cand.p4();
    auto q = q_+cand.q();
    Candidate sum(p4,q);
    sum.addComponent(*this);
    sum.addComponent(cand);
    return sum;
}

// operators
Candidate & Candidate::operator += (const Candidate & cand) 
{
    p4_+=cand.p4();
    q_+=cand.q();
    this->addComponent(cand);
    return *this;
}


//
// member functions
//
void Candidate::addComponent(const Candidate & cand)
{
   components_.push_back(&cand);
}

std::vector<const Candidate*> Candidate::components() const
{
   return components_;
}

bool Candidate::matchTo(const std::vector<Candidate> * cands, const std::string & name, const float & deltar_max)
{
   bool status = false;


   if ( ! cands )
   {
      this -> matched_[name] = nullptr;
      return status;
   }

   const Candidate * cand = nullptr;
   const Candidate * nearcand = nullptr;
   float minDeltaR = 100.;
   for ( size_t i = 0; i < cands->size() ; ++i )
   {
      cand = &(cands->at(i));
      if(this->deltaR(*cand) < minDeltaR)
      {
         minDeltaR = this->deltaR(*cand);
         nearcand = cand;
      }
   }

   if(minDeltaR < deltar_max)
   {
     this->matched_[name]=nearcand;
     status = true;
   }

   else {
     this -> matched_[name] = nullptr;
   }

   return status;
}

bool Candidate::matchTo(const std::vector<Candidate> * cands, const std::string & name, const float & deltar_max, const float & deltaptrel_max)
{
   bool status = false;


   if ( ! cands )
   {
      this -> matched_[name] = nullptr;
      return status;
   }

   const Candidate * cand = nullptr;
   const Candidate * nearcand = nullptr;
   float minDeltaR = deltar_max + 1; 		// Assign more real value;
   float dpT = 0.;
   float dpTmin = deltaptrel_max + 1;
   for ( size_t i = 0; i < cands->size() ; ++i )
   {
      cand = &(cands->at(i));
      dpT = std::abs(this->pt() - cand->pt());
      if(this->deltaR(*cand) < minDeltaR && dpT < dpTmin)
      {
         minDeltaR = this->deltaR(*cand);
         dpTmin    = dpT;
         nearcand = cand;
      }
   }

   if(minDeltaR < deltar_max && dpTmin < deltaptrel_max)
   {
     this->matched_[name]=nearcand;
     status = true;
   }

   else {
     this -> matched_[name] = nullptr;
   }

   return status;
}

// Gets
float Candidate::px()   const { return p4_.Px() ; }
float Candidate::py()   const { return p4_.Py() ; }
float Candidate::pz()   const { return p4_.Pz() ; }
float Candidate::pt()   const { return p4_.Pt() ; }
float Candidate::eta()  const { return p4_.Eta(); }
float Candidate::phi()  const { return p4_.Phi(); }
float Candidate::e()    const { return p4_.E()  ; }
float Candidate::m()    const { return p4_.M()  ; }
float Candidate::mass() const { return p4_.M()  ; }
int   Candidate::q()    const { return q_;   }
float Candidate::deltaR(const Candidate &cand) const { return p4_.DeltaR(cand.p4()) ;}
float Candidate::deltaPhi(const Candidate &cand) const { return p4_.DeltaPhi(cand.p4()) ;}
float Candidate::deltaEta(const Candidate &cand) const { return fabs(this->eta()-cand.eta());}

TLorentzVector Candidate::p4() const { return p4_; }
TVector3       Candidate::p3() const { return p4_.Vect(); }


const Candidate * Candidate::matched(const std::string & name) { return matched_[name]; }
const Candidate * Candidate::matched(const std::string & name) const { return matched_.find(name) != matched_.end() ? matched_.find(name)->second : 0; }

void Candidate::unmatch(const std::string & name)
{
    matched_[name] = nullptr;
}

float Candidate::matchedDeltaR(const std::string &name) const
{
   return this->deltaR(*(this->matched(name))) ;
}

// Sets
void  Candidate::p4(const TLorentzVector & p4) { p4_ = p4; }
void  Candidate::px(const float & px) { p4_.SetPx(px); }
void  Candidate::py(const float & py) { p4_.SetPy(py); }
void  Candidate::pz(const float & pz) { p4_.SetPy(pz); }
void  Candidate::e (const float & e ) { p4_.SetE(e);   }
void  Candidate::q (const float & q)  { q_ = q; }

// print 
void Candidate::printInfo(const std::string & type) const
{
   std::string out = type == "" ? "Candidate :" : type+" :"; 
   std::cout << std::fixed << std::setprecision(5);
   std::cout << out << this << ": "; 
   std::cout << "pt  = " << std::setw(10) << this->pt()  << " | ";
   std::cout << "eta = " << std::setw(10) << this->eta() << " | ";
   std::cout << "phi = " << std::setw(10) << this->phi() << " \n ";
   
}

void Candidate::printMatchedInfo(const std::string & name, const std::string & type) const
{
   std::string out = (type == "") ? "Candidate" : type;
   out += " matched to "+name+" :";
   auto match = this->matched(name);
   std::cout << std::fixed << std::setprecision(5);
   if ( ! match )
   {
      std::cout << out << "  NO MATCHING" << std::endl;
   }
   else
   {
      std::cout << out << match << ": "; 
      std::cout << "dR  = " << std::setw(10) << this->matchedDeltaR(name)  << " | ";
      std::cout << "pt  = " << std::setw(10) << match->pt()  << " | ";
      std::cout << "eta = " << std::setw(10) << match->eta() << " | ";
      std::cout << "phi = " << std::setw(10) << match->phi() << " \n ";
   }
}


void Candidate::listMatchedNames() const
{
   std::cout << "List of names of matched candidates\n";
   for ( auto const& match: matched_ )
   {
      std::cout << match.first << std::endl;
   }
}

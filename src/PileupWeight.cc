// system include files
// user include files
#include <iostream>
#include <fstream>
#include "Analysis/Tools/interface/PileupWeight.h"
#include "TFile.h"
#include "ROOT/RCsvDS.hxx"
#include <boost/tokenizer.hpp>



//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;

//
// constructors and destructor
PileupWeight::PileupWeight()
{
}

// Main constructor
PileupWeight::PileupWeight(const std::string & puweight_name, const std::string & pudata_name, const bool & is_mc )
{
   puweight_name_ = puweight_name;
   pudata_name_ = pudata_name;
   is_mc_ = is_mc;

   TH1D * tmp;
   std::shared_ptr<TFile> f = std::make_shared<TFile>(puweight_name.c_str(),"old");
   tmp = (TH1D*) f->Get("weight_2down");
   if ( tmp )  puweight_histos_[-2] = std::make_shared<TH1D>(*tmp);
   
   tmp = (TH1D*) f->Get("weight_1down");
   if ( tmp )  puweight_histos_[-1] = std::make_shared<TH1D>(*tmp);
   
   tmp = (TH1D*) f->Get("weight");
   if ( tmp )  puweight_histos_[0] = std::make_shared<TH1D>(*tmp);
   
   tmp = (TH1D*) f->Get("weight_1up");
   if ( tmp )  puweight_histos_[1] = std::make_shared<TH1D>(*tmp);
   
   tmp = (TH1D*) f->Get("weight_2up");
   if ( tmp )  puweight_histos_[2] = std::make_shared<TH1D>(*tmp);

   if ( ! puweight_histos_[-2] ) 
      std::cout << "WARNING - PileupWeight::PileupWeight | Histogram weight_2down not found. Weight = 1" << std::endl;
   if ( ! puweight_histos_[-1] ) 
      std::cout << "WARNING - PileupWeight::PileupWeight | Histogram weight_1down not found. Weight = 1" << std::endl;
   if ( ! puweight_histos_[0] ) 
      std::cout << "WARNING - PileupWeight::PileupWeight | Histogram weight not found. Weight = 1" << std::endl;
   if ( ! puweight_histos_[1] ) 
      std::cout << "WARNING - PileupWeight::PileupWeight | Histogram weight_1up not found. Weight = 1" << std::endl;
   if ( ! puweight_histos_[2] ) 
      std::cout << "WARNING - PileupWeight::PileupWeight | Histogram weight_2up not found. Weight = 1" << std::endl;


   if ( ! is_mc_ )
   {
      if ( pudata_name_ != "" )
      {         
         //open csv with RDataFrame (not workin on old root verions)
         // df_pudata_ = std::make_shared<ROOT::RDataFrame>(ROOT::RDF::MakeCsvDataFrame(pudata_name_));

         // open CSV and store on map         
         std::ifstream file(pudata_name_);
         // Ignore the first row
         std::string first_row;
         std::getline(file, first_row);
         
         std::string line;
         while (std::getline(file, line)) {
            // Use Boost tokenizer to split the line into tokens
            typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;
            Tokenizer tokens(line);

            auto it = tokens.begin();
            if (it != tokens.end()) {
                  int run = std::stoi(*it++);
                  if (it != tokens.end()) {
                     int ls = std::stoi(*it++);
                     if (it != tokens.end()) {
                        float avgpu = std::stof(*it++);
                        m_pudata_[run][ls] = avgpu;
                     }
                  }        
            } else {
                  std::cerr << "Failed to parse line: " << line << std::endl;
            }
         }

         file.close();

      }
   }
}

PileupWeight::~PileupWeight()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
// ------------ method called for each event  ------------

// output tree
float PileupWeight::weight(const float & truepu, const int & var)
{
   float weight = 1.;
   if ( ! puweight_histos_[0] ) return weight;
   if ( ! puweight_histos_[var] ) return weight;
   if ( ! is_mc_ ) // For Data
   {
      weight = puweight_histos_[var]->Interpolate(truepu);
      return weight;
   }
   // For MC
   int bin = puweight_histos_[var] -> FindBin(truepu);
   weight = puweight_histos_[var] -> GetBinContent(bin);
   return weight;

   return weight;
}

float PileupWeight::getPileupFromData(const int & myrun, const int & myls)
{

   float pileup = -1;
   if ( pudata_name_ == "" ) return pileup;

   pileup = m_pudata_[myrun][myls];

   return pileup;

}

// float PileupWeight::getPileupFromDataFrame(const int & myrun, const int & myls)
// {

//    float pileup = -1;
//    if ( pudata_name_ == "" ) return pileup;

//    auto selectedData = df_pudata_->Filter(
//       [=](Long64_t run, Long64_t lumi_section) {
//          return run == myrun && lumi_section == myls;
//       },
//       {"run", "lumi_section"}
//    );
//    std::vector<double> avgpuValues = *selectedData.Take<double>("avgpu");
//    if ( avgpuValues.size() == 0 ) return -1;
//    pileup = avgpuValues[0];
//    return pileup;

// }
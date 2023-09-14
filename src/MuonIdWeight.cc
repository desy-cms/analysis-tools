// system include files
// user include files
#include "Analysis/Tools/interface/MuonIdWeight.h"
#include "TFile.h"

//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;

//
// constructors and destructor
MuonIdWeight::MuonIdWeight()
{
}

// Main constructor
MuonIdWeight::MuonIdWeight(const std::vector<std::string> & fnames ) 
{
   TFile * file;
   TH2F * hist_sf;
   TH2F * hist_unc;
   float range_pTmin = -1, range_pTmax = -1;

   
   for(unsigned int f = 0; f < fnames.size(); f++)
   {
      file = new TFile(fnames[f].c_str(),"READ");
      hist_sf = (TH2F*)((TH2F*)file->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt"))->Clone();
      hist_unc = (TH2F*)((TH2F*)file->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_unc"))->Clone();
      if(! hist_sf){std::cout<<"WARNING no muonID scale factor histogram on the file provided"; continue;}
      if(! hist_unc){std::cout<<"WARNING no muonID scale factor uncertainty histogram on the file provided"; continue;}

      range_pTmin = hist_sf -> GetYaxis() -> GetBinLowEdge(1);
      range_pTmax = hist_sf -> GetYaxis() -> GetBinUpEdge(hist_sf -> GetYaxis() -> GetNbins());
      fabs_etamax = hist_sf -> GetXaxis() -> GetBinUpEdge(hist_sf -> GetXaxis() -> GetNbins());
      pTranges_file.push_back(std::make_pair(std::make_pair(range_pTmin, range_pTmax), fnames[f]));

      if (sf_pTmin > range_pTmin) {sf_pTmin = range_pTmin;}
      if (sf_pTmax < range_pTmax) {sf_pTmax = range_pTmax;}

      //sf_hists.push_back(hist_sf); //do you need this?
      //filename_sf_pairs.push_back(std::make_pair(fnames[f], * hist_sf));
      //filename_unc_pairs.push_back(std::make_pair(fnames[f], * hist_unc));
      filename_sf_map[fnames[f]] = *hist_sf;
      filename_unc_map[fnames[f]] = *hist_unc;

      delete hist_sf;
      delete hist_unc;
      file -> Close();
   }

}



MuonIdWeight::~MuonIdWeight()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//
// ------------ method called for each event  ------------


float MuonIdWeight::findSF(const float & pt, const float & eta , const int & var) 
{
   float eta_aux = eta, pt_aux = pt;
   float sf = 1.;
   float nominal_sf = 1., sf_unc = 1.;
   bool twice_unc = false;
   std::string filename = findfile(pt); //use findfile to know which histogram of sf_files_hists to access
   TH2F scalefactor_hist;
   TH2F uncert_hist;

   if (filename_sf_map.find(filename) != filename_sf_map.end())   { scalefactor_hist = filename_sf_map[filename];}  else {std::cout << "WARNING File name not found." << std::endl; }
   if (filename_unc_map.find(filename) != filename_unc_map.end()) { uncert_hist      = filename_unc_map[filename];} else {std::cout << "WARNING File name not found." << std::endl; }

   //some minor checks/corrections   
   if (fabs(eta) >= fabs_etamax){eta_aux = fabs_etamax - 0.00001;}        // muons with eta > 2.4 (very unlikely), take max eta to avoid errors
   if (pt <= sf_pTmin) {pt_aux = sf_pTmin + 0.00001;}                     // muons with pT < lowest pT for which the SF is defined (sf_pTmin), take sf of sf_pTmin, also rather unlikely 
   if (pt >= sf_pTmax) {pt_aux = sf_pTmax - 0.00001; twice_unc = true;}   // muons with pT > highest pT for which the SF is defined (sf_pTmax), take twice the uncertainty

   if (var == 0) // get nominal scale factor
   {
      // Find the bin indices for eta_aux and pt_aux
      Int_t binX = scalefactor_hist.GetXaxis()->FindBin(fabs(eta_aux));//symmetric in eta
      Int_t binY = scalefactor_hist.GetYaxis()->FindBin(pt_aux);

      // Get the scale factor value from the histogram
      nominal_sf = scalefactor_hist.GetBinContent(binX, binY);
      sf = nominal_sf;

   }    
   else //get systematic uncertainties
   {
      Int_t binX = scalefactor_hist.GetXaxis()->FindBin(fabs(eta_aux));
      Int_t binY = scalefactor_hist.GetYaxis()->FindBin(pt_aux);
      
      nominal_sf = scalefactor_hist.GetBinContent(binX, binY);
      sf_unc     = uncert_hist.GetBinContent(binX, binY);
      if (!twice_unc)
      sf = nominal_sf + var * sf_unc;
      else
      sf = nominal_sf + 2. * var * sf_unc;
   }

   return sf;
}


std::string MuonIdWeight::findfile(const float & pT)
{
   std::string name = "";

   for (const auto& entry : pTranges_file)
   {
      const std::pair<double, double>& pTrange = entry.first;
      const std::string& filename = entry.second;
      if (pTrange.first <= pT && pT <= pTrange.second) { name = filename; break;}
   }
   

   if(name == "" && pT > sf_pTmax) { name = findfile(sf_pTmax); } //do you need to do this here?

   return name;
}


// unused functions
/*float MuonIdWeight::weight(const float & pt, const float & eta , const int & var) 
{
   float weight = 1.;
   if ( ! histos_[0] ) return weight;
   if ( ! histos_[var] ) return weight;
   int binx = histos_[var] -> GetXaxis()->FindBin(pt);
   int biny = histos_[var] -> GetYaxis()->FindBin(eta);
   weight = histos_[var] -> GetBinContent(binx, biny);
   return weight;
}*/


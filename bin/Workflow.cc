#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>

#include "TFile.h" 
#include "TH1.h" 

#define ANSI_COLOR_RED     "\033[0;31m"
#define ANSI_COLOR_RESET   "\033[0m"

// =============================================================================================   
int main(int argc, char * argv[])
{

   std::string color_red = "\033[0;31m";
   std::string color_reset  = "\033[0m";
   if ( argc <= 1 ) return -1;
   const char * rootfile = argv[1];
   std::string base_filename = std::string(rootfile).substr(std::string(rootfile).find_last_of("/") + 1);

   std::string opt = "";
   if  ( argc == 3 )  opt = std::string(argv[2]);
   bool do_csv = false;
   std::string csv_filename = "";
   std::ofstream csvfile;
   if ( opt == "csv" )
   {
      do_csv = true;
      csv_filename=base_filename.substr(0,base_filename.find_last_of('.'))+".csv";
      csv_filename="workflow_"+csv_filename;
   }
   
   TFile f(rootfile,"old");
   
   TH1F * h = (TH1F*) f.Get("workflow");
   
   printf("+%s%s%s+\n", color_red.c_str(), std::string(170,'-').c_str(),color_reset.c_str());
   printf("| %-108s |    %10s |   %16s |   %16s |\n",h->GetTitle(),"n events","ratio wrt first","ratio wrt previous");
   printf("+%s%s%s+\n", color_red.c_str(), std::string(170,'-').c_str(),color_reset.c_str());
   if ( do_csv )
   {
      csvfile.open(csv_filename.c_str());
      csvfile << h->GetTitle()<<";label;n events;ratio wrt first;ratio wrt previous" << std::endl;
   }
   int firstbin = 2;
   for ( int i = 1; i <= h ->GetNbinsX(); ++i )
   {
      std::string label = std::string(h->GetXaxis()->GetBinLabel(i));
      if ( label == "" ) continue;
//      if ( firstbin < 0 ) firstbin = i;
      float n = h->GetBinContent(i);
      float rn1 = h->GetBinContent(i)/h->GetBinContent(firstbin);
      float rni = 0;
      if ( i == 1 )
      {
         printf("| %2d - %-103s |    %10.1f |   %16s |  %19s |\n",i-1,label.c_str(),n,"-","-");
         csvfile << i-1 << ";" << label << ";" << n << ";" << ";" << ";" << std::endl;
      }
      else if ( i == 2 )
      {
         printf("| %2d - %-103s |    %10.1f |   %16.6f |  %19s |\n",i-1,label.c_str(),n,rn1,"-");
         csvfile << i-1 << ";" << label << ";" << n << ";" << rn1 << ";" << std::endl;
      }
      else
      {
         rni = h-> GetBinContent(i)/h->GetBinContent(i-1);
         printf("| %2d - %-103s |    %10.1f |   %16.6f |     %16.6f |\n",i-1,label.c_str(),n,rn1,rni);
         csvfile << i-1 << ";" << label << ";" << n << ";" << rn1 << ";" << rni << std::endl;
      }
      
   }
   printf("+%s%s%s+\n", color_red.c_str(), std::string(170,'-').c_str(),color_reset.c_str());

   if ( do_csv )
      csvfile.close();

   return 0;
   
}


#include "TFile.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <utility>
#include "TSpectrum.h"
#include <fstream>
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include "TString.h"
#include <sstream>


// MACRO TO combine the efficiencies and put into a format for BAT-hl
// Toby Dixon: toby.dixon@universite-paris-saclay.fr 23/5/2023    
// Possibility to convolve the systematics (??)


TH1D* Combine(TString path_P,TString path_M)
{



  // Get the histograms
  //------------------------------------------------------------------

  TFile *fP = new TFile(path_P);
  TH1D *hP = (TH1D*)fP->Get("h1_fit_observable_fQ");
  

  TFile *fM = new TFile(path_M);
  TH1D *hM = (TH1D*)fM->Get("h1_fit_observable_fQ");


  TH1D *heff= new TH1D("heff","heff",1000,0,1);
  

  for (int i=0;i<1000000;i++)
    {


      double effP=hP->GetRandom();
      double effM=hM->GetRandom();

      heff->Fill(effP*effM);

    }
  heff->Scale(1/(heff->GetBinWidth(1)*heff->Integral()));
  
  return heff;
}


int main()
{

  std::vector<int>DS{3601,3602,3603,3604,3605,3606,3607,3607,3608,3609};

  TString path = "output/CUORE_eff";
  TH1D *heff;
  for (auto &D:DS)
    {
      heff = (TH1D*)Combine(Form("%s_P/ds%i/output_eff.root",path.Data(),D),
			    Form("%s_M/ds%i/output_eff.root",path.Data(),D));

      heff->SaveAs(Form("%s/efficiency_ds%i.root",path.Data(),D));

    }
}

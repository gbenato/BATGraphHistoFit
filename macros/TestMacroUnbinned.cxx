
#include "src/BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"

int main(int argc,char**argv)
{


  // A SET OF TESTS



  // **************** HISTOGRAM FIT **********************************
  // *****************************************************************

  TF1 *fGauss = new TF1("fGauss","[0]+[1]*TMath::Gaus(x,[2],[3])/(sqrt(2*TMath::Pi())*[3])",2500,2700);

  fGauss->SetParameter(0,1);
  fGauss->SetParameter(1,150);
  fGauss->SetParameter(2,2615);
  fGauss->SetParameter(3,3);

  TH1D *fHist = new TH1D("fHist","fHist",200,2500,2700);
  fHist->SetTitle(" ; Energy [keV]; Counts ;  ");
  double N = fGauss->Integral(2500,2700);
  std::vector<double>energy;
  for (int i=0;i<round(N);i++)
    {
      double E=fGauss->GetRandom();
      fHist->Fill(E);
      energy.push_back(E);
    }
  fGauss->SetParameter(0,5);
  fGauss->SetParLimits(0,0,20);
  fGauss->SetParameter(1,60);
  fGauss->SetParLimits(1,20,400);
  fGauss->SetParameter(2,2600);
  fGauss->SetParLimits(2,2550,2650);
  fGauss->SetParameter(3,5);
  fGauss->SetParLimits(3,0,10);



  //Binned fit
  BatGraphFitter *fitter2 = new BatGraphFitter(fHist,fGauss);
  fitter2->SetPrecison(3);
  TCanvas *c2= new TCanvas("c2","c2");
  c2->cd();
  fitter2->Fit();
  c2->Draw();
  c2->Print("fit.pdf(","pdf");
  fitter2->fModel->PrintAllMarginalized("out_hist_plots.pdf");
  fitter2->fModel->WriteMarginalizedDistributions("output_hist_dist.root", "RECREATE");
  fitter2->fModel->WriteMarkovChain("output_hist_mcmc.root", "RECREATE");

  BatGraphFitter *fitter3 = new BatGraphFitter(fHist,fGauss,energy);
  fitter3->SetPrecison(3);
  c2->cd();
  fitter3->Fit();
  c2->Draw();
  c2->Print("fit.pdf)","pdf");
  fitter3->fModel->PrintAllMarginalized("out_unbinned_hist_plots.pdf");
  fitter3->fModel->WriteMarginalizedDistributions("output_unbinned_hist_dist.root", "RECREATE");
  fitter3->fModel->WriteMarkovChain("output_unbinned_hist_mcmc.root", "RECREATE");






  

}

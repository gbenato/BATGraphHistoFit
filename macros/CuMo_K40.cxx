
#include "BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include "TH3D.h"
int main(int argc,char**argv)
{


  // A SET OF TESTS

  //**************** GRAPH FIT ******************
  //*********************************************

  
  TApplication *app  =new TApplication("app",&argc,argv);
  TCanvas *c =new TCanvas();

  TFile *file = new TFile("~/Downloads/output.root");

  file->ls();
  TH1D *h1=(TH1D*)file->Get("hfirst");
  TH1D *h2=(TH1D*)file->Get("hmiddle");
  TH1D *h3=(TH1D*)file->Get("hlast");


  h1->GetXaxis()->SetRangeUser(1420,1500);
  gStyle->SetOptStat(0);
  BatGraphFitter *fitter1 = new BatGraphFitter(h1);
  fitter1->SetPrecison(3);
  BatGraphFitter *fitter2 = new BatGraphFitter(h2);
  fitter1->SetPrecison(3);
  BatGraphFitter *fitter3 = new BatGraphFitter(h3);
  fitter3->SetPrecison(3);

  TF1 * f = new TF1("f","[0]+[1]*(x-1460)+[4]*TMath::Gaus(x,[2],[3])/(sqrt(2*3.14)*[3])",1440,1480);

  f->SetParLimits(0,100,200);
  f->SetParLimits(1,-2,2);
  f->SetParLimits(2,1455,1470);
  f->SetParLimits(3,1,5);
  f->SetParLimits(4,0,400);

  TCanvas *c1=new TCanvas();
  fitter1->Fit(f,"R","",1440,1480);
  c1->Draw();
  c1->Print("fit.pdf(","pdf");
  fitter1->fModel->WriteMarkovChain("output_mcmc1.root", "RECREATE");

  TCanvas *c2=new TCanvas();
  c2->cd();
  fitter2->Fit(f,"R","",1440,1480);
  c2->Draw();
  c2->Print("fit.pdf","pdf");
  fitter2->fModel->WriteMarkovChain("output_mcmc2.root", "RECREATE");

  TCanvas *c3=new TCanvas();
  c3->cd();
  fitter3->Fit(f,"R","",1440,1480);
  c->Draw();
  c->Print("fit.pdf)","pdf");
  fitter3->fModel->WriteMarkovChain("output_mcmc3.root", "RECREATE");

  app->Run();


}


#include "BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>

#include "TF1.h"
#include "TApplication.h"
int main(int argc,char**argv)
{

  
  TCanvas *c =new TCanvas();
  TFile *f = new TFile("../Cd113Shape/data_CWOnat.root");

  TF1 *func =(TF1*)f->Get("funeff");
  func->SetParameter(0,0.01);
  func->SetParLimits(0,0,0.2);
  func->SetParameter(1,50);
  func->SetParLimits(1,0,100);
  func->SetParameter(2,0.01);
  func->SetParLimits(2,-.1,0.2);
  TGraphAsymmErrors * g =(TGraphAsymmErrors*)f->Get("Graph_Efficiency");
  
  BatGraphFitter *fitter = new BatGraphFitter(g);
  std::cout<<"First a ROOT Fit"<<std::endl;
 g->Fit(func);
   func->SetParameter(0,0.01);
  func->SetParLimits(0,0,0.2);
  func->SetParameter(1,50);
  func->SetParLimits(1,0,100);
  func->SetParameter(2,0.01);
  func->SetParLimits(2,-.1,0.2);

  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"NOW FIT with bat"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  
  fitter->Fit(func);
  c->Draw();
  c->Print("fit.pdf");
  
  fitter->fModel->PrintAllMarginalized("out_plots.pdf");
  fitter->fModel->WriteMarkovChain("output_mcmc.root", "RECREATE");


  



}


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

  TFile *file = new TFile("resolution_output.root");
  file->ls();
  TGraphErrors *g = (TGraphErrors*)file->Get("graph_energybias");
  BatGraphFitter *fitter = new BatGraphFitter(g);
  fitter->SetPrecison(4);
  fitter->SetGraphMaxMin(5,-5);

  TF1 * f = new TF1("f","pol2",0,3000);
  std::cout<<"First a ROOT Fit"<<std::endl;
  
  g->Fit(f);
  f->SetParLimits(0,-2,2);
  f->SetParLimits(1,-2/1000.,1/1000.);
  f->SetParLimits(2,-1e-6,1e-6);

  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"NOW FIT with bat"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  
  fitter->Fit(f,"C");
  
  c->Draw();
  c->Print("fit.pdf","pdf");
  
  fitter->fModel->PrintAllMarginalized("out_plots.pdf");
  fitter->fModel->WriteMarkovChain("output_mcmc.root", "RECREATE");
  TFile *fout = new TFile("output_marg.root","recreate");
    
  TTree *T = (TTree*)fitter->fModel->GetMarkovChainTree();

  TH3D *h_bias = new TH3D("h_bias","h_bias",200,-5e-7,5e-7,200,-2/1000.,1/1000.,200,-2,2);

  T->Draw("p0:p1:p2>>h_bias","phase==1");
  
  h_bias->Write();
  

  app->Run();


}

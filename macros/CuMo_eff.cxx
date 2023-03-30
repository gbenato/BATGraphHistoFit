
#include "BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include "TH3D.h"
#include "TH2D.h"

int main(int argc,char**argv)
{


  // A SET OF TESTS

  //**************** GRAPH FIT ******************
  //*********************************************

  std::string name=argv[1];
  std::string type=argv[2];
  std::cout<<argv[0]<<" "<<argv[1]<<" "<<argv[2]<<" "<<std::endl;

  TApplication *app  =new TApplication("app",&argc,argv);
  TCanvas *c =new TCanvas();

  TFile *file = new TFile((TString)name);
  file->ls();
  TGraphAsymmErrors *g = (TGraphAsymmErrors*)file->Get("graph");
  BatGraphFitter *fitter = new BatGraphFitter(g);
  fitter->SetQbb(3034.4);
  fitter->SetPrecison(4);
  fitter->SetGraphMaxMin(1,0);

  TF1 * f = new TF1("f","pol1",0,3000);
  std::cout<<"First a ROOT Fit"<<std::endl;
  
  g->Fit(f);
  f->SetParLimits(0,0,1.2);
  f->SetParLimits(1,-0.0001,0.0001);

  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"NOW FIT with bat"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  
  fitter->Fit(f,"C");
  
  c->Draw();
  c->Print("fit.pdf","pdf");
  
  fitter->fModel->PrintAllMarginalized(type+"_out_plots.pdf");
  fitter->fModel->WriteMarkovChain(type+"_output_mcmc.root", "RECREATE");
  TFile *fout = new TFile(TString(type+"_output_marg.root"),"recreate");
    
  TTree *T = (TTree*)fitter->fModel->GetMarkovChainTree();

  TH2D *h_bias = new TH2D("h_bias","h_bias",2000,-0.0001,0.0001,2000,0.8,1.2);

  T->Draw("p0:p1>>h_bias","phase==1");
  
  h_bias->Write();
  

  app->Run();


}

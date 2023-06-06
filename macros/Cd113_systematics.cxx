
#include "src/BatGraphFitter.h"
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

  TCanvas *c =new TCanvas();

  TFile *file = new TFile("inputs/cross/data_CWOnat.root","UPDATE");
  TGraphAsymmErrors *g = (TGraphAsymmErrors*)file->Get("Graph_Efficiency");
  TF1 * f =(TF1*)file->Get("funeff");
  f->SetParameter(0,0.1);
  f->SetParameter(1,40);
  f->SetParameter(2,0.1);
  g->Fit(f);
  f->SetParLimits(0,f->GetParameters()[0]-5*f->GetParErrors()[0],f->GetParameters()[0]+5*f->GetParErrors()[0]);
  f->SetParLimits(1,0,f->GetParameters()[1]+5*f->GetParErrors()[1]);
  f->SetParLimits(2,0,f->GetParameters()[2]*5);

  BatGraphFitter *fitter = new BatGraphFitter(g,f);
  fitter->fModel->WriteMarkovChain("output/cross/eff_output_mcmc.root", "RECREATE");

  fitter->SetPrecison(5);
  fitter->SetGraphMaxMin(1,0);  

  fitter->Fit("C");
  
  c->Draw();
  c->Print("fit.pdf(","pdf");

  TGraphAsymmErrors *gbias = (TGraphAsymmErrors*)file->Get("Graph_Energy_Bias_Residual_corrected");

  TF1 *fbias = new TF1("fbias","[0]",0,2000);
  fbias->SetParLimits(0,-1,1);
  
  BatGraphFitter *fitterbias = new BatGraphFitter(gbias,fbias);
  fitterbias->fModel->WriteMarkovChain("output/cross/bias_output_mcmc.root", "RECREATE");

  fitterbias->SetPrecison(5);
  fitterbias->SetGraphMaxMin(10,-10);
  fitterbias->Fit("C");
  c->Draw();
  c->Print("fit.pdf)","pdf");
  file->cd();
  f->Write();
  fbias->Write();
  file->Write();

  
  fitter->fModel->PrintAllMarginalized("output/cross/eff_out_plots.pdf");
  fitterbias->fModel->PrintAllMarginalized("output/cross/bias_out_plots.pdf");
  

  // TFile *fout = new TFile(TString("output/cross/eff_output_marg.root"),"recreate");



  // NOW FIT THE BIAS ----------------------------------------------------------
  // ***************************************************************************

  

}

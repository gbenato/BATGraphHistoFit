#include "../src/BatGraphFitter.h"
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
  TString path=argv[1];

  TCanvas *c =new TCanvas();


  TFile *file = new TFile(Form("%s/data_CWOnat.root",path.Data()),"UPDATE");
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
  fitter->fModel->WriteMarkovChain(Form("%s/eff_output_mcmc.root",path.Data()), "RECREATE");

  fitter->fModel->SetNIterationsRun(50000000);
  fitter->SetPrecison(5);
  fitter->SetGraphMaxMin(1,0);

  fitter->Fit("C");

  c->Draw();
  c->Print(Form("%s/fit_eff.C",path.Data()));
  TGraphAsymmErrors *gbias = (TGraphAsymmErrors*)file->Get("Graph_Energy_Bias_Residual_corrected");

  TF1 *fbias = new TF1("fbias","[0]",0,2000);
  fbias->SetParLimits(0,-1,1);

  BatGraphFitter *fitterbias = new BatGraphFitter(gbias,fbias);
  fitterbias->fModel->WriteMarkovChain(Form("%s/bias_output_mcmc.root",path.Data()), "RECREATE");
  
  fitterbias->SetPrecison(5);
  fitterbias->SetGraphMaxMin(10,-10);
  fitterbias->Fit("C");
  c->Draw();
  c->Print(Form("%s/fit_bias.C",path.Data()));	

  file->cd();
  f->Write();
  fbias->Write();


  fitter->fModel->PrintAllMarginalized(Form("%s/eff_out_plots.pdf",path.Data()));
  fitterbias->fModel->PrintAllMarginalized(Form("%s/bias_out_plots.pdf",path.Data()));




  // Fit the reso                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
  // *****************************************                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
  TGraphAsymmErrors *greso = (TGraphAsymmErrors*)file->Get("Graph_Resolution");

  TF1 *freso = new TF1("freso","sqrt([0]*[0]+[1]*[1]*x+[2]*[2]*x*x)",0,2000);
  freso->SetParameter(0,5);
  freso->SetParameter(1,sqrt(5/2000.));
  freso->SetParameter(2,sqrt(5/(2000.*2000.)));

  greso->Fit(freso);
  freso->SetParLimits(0,0,freso->GetParameters()[0]+5*freso->GetParErrors()[0]);
  freso->SetParLimits(1,0,freso->GetParameters()[1]+5*freso->GetParErrors()[1]);
  freso->SetParLimits(2,0,freso->GetParameters()[2]*5);


  BatGraphFitter *fitterreso = new BatGraphFitter(greso,freso);
  fitterreso->fModel->WriteMarkovChain(Form("%s/reso_output_mcmc.root",path.Data()), "RECREATE");

  fitterreso->SetPrecison(5);
  fitterreso->SetGraphMaxMin(10,0);
  fitterreso->Fit("C");
  TGraphAsymmErrors * grint = (TGraphAsymmErrors*)fitterreso->fGrint;

  TGraph *gLow=new TGraph();
  TGraph *gMode = new TGraph();
  TGraph *gHigh=new TGraph();
  for (int i=0;i<grint->GetN();i++)
    {
      double x= grint->GetX()[i];
      double y=grint->GetY()[i];
      double ey_low= grint->GetErrorYlow(i);
      double ey_high=grint->GetErrorYhigh(i);

      gLow->SetPoint(i,x,y-ey_low);
      gMode->SetPoint(i,x,y);
      gHigh->SetPoint(i,x,y+ey_high);

    }
  


  c->Draw();
  c->Print(Form("%s/fit_reso.C",path.Data()));	

  file->cd();
  gLow->SetName("low_reso");
  gMode->SetName("best_fit_reso");
  gHigh->SetName("high_reso");
  gLow->Write();
  gHigh->Write();
  gMode->Write();
  f->Write();
  freso->Write();
  

  fitterreso->fModel->PrintAllMarginalized(Form("%s/reso_out_plots.pdf",path.Data()));

  file->Close();
					   



}



#include "BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"

int main(int argc,char**argv)
{


  // A SET OF TESTS

  //**************** GRAPH FIT ******************
  //*********************************************

  
  TApplication *app  =new TApplication("app",&argc,argv);
  TF1 * f = new TF1("f","[0]*log(x+[1])+[2]",0,4000);

  f->SetParameter(0,0.1);
  f->SetParameter(1,40);
  f->SetParameter(2,0.1);

  TRandom3 *rand = new TRandom3();

  std::vector<double>E{1,5,10,50,100,150,200,500,1000,2000,2500,3000};
  TGraphAsymmErrors *g = new TGraphAsymmErrors();
  for (int i=0;i<E.size();i++)
    {
      double error =abs(rand->Gaus(0.02,0.01));
      double y=rand->Gaus(f->Eval(E[i]),error);
      g->SetPoint(i,E[i],y);
      g->SetPointEYhigh(i,error);
      g->SetPointEYlow(i,error);
    }
  
  g->SetTitle("Test data for the Graph ; Energy [keV]; y ; ");
  // ******* GAUSSIAN GRAPH ********* //
  TCanvas *c =new TCanvas();
  f->SetParameter(0,0.01);
  f->SetParLimits(0,0,0.2);
  f->SetParameter(1,50);
  f->SetParLimits(1,0,120);
  f->SetParameter(2,0.01);
  f->SetParLimits(2,-.1,0.3);
 
  BatGraphFitter *fitter = new BatGraphFitter(g);
  fitter->SetPrecison(3);
  fitter->SetGraphMaxMin(1,0);

  std::cout<<"First a ROOT Fit"<<std::endl;
  g->Fit(f);
  
  f->SetParameter(0,0.01);
  f->SetParLimits(0,0,0.3);
  f->SetParameter(1,50);
  f->SetParLimits(1,0,120);
  f->SetParameter(2,0.01);
  f->SetParLimits(2,-.2,0.3);

  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"NOW FIT with bat"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  
  fitter->Fit(f);
  
  c->Draw();
  c->Print("fit.pdf(","pdf");
  
  fitter->fModel->PrintAllMarginalized("out_plots.pdf");
  fitter->fModel->WriteMarkovChain("output_mcmc.root", "RECREATE");


  // **************** HISTOGRAM FIT **********************************
  // *****************************************************************

  TF1 *fGauss = new TF1("fGauss","[0]+[1]*TMath::Gaus(x,[2],[3])/(sqrt(2*TMath::Pi())*[3])",2500,2700);

  fGauss->SetParameter(0,10);
  fGauss->SetParameter(1,150);
  fGauss->SetParameter(2,2615);
  fGauss->SetParameter(3,3);

  TH1D *fHist = new TH1D("fHist","fHist",200,2500,2700);
  fHist->SetTitle(" ; Energy [keV]; Counts ;  ");
  double N = fGauss->Integral(2500,2700);
  for (int i=0;i<round(N);i++)
    {
      fHist->Fill(fGauss->GetRandom());
    }
  fGauss->SetParameter(0,5);
  fGauss->SetParLimits(0,0,20);
  fGauss->SetParameter(1,60);
  fGauss->SetParLimits(1,20,400);
  fGauss->SetParameter(2,2600);
  fGauss->SetParLimits(2,2550,2650);
  fGauss->SetParameter(3,5);
  fGauss->SetParLimits(3,0,10);
  
  BatGraphFitter *fitter2 = new BatGraphFitter(fHist);
  fitter2->SetPrecison(3);
  TCanvas *c2= new TCanvas("c2","c2");
  c2->cd();
  fitter2->Fit(fGauss);
  c2->Draw();
  c2->Print("fit.pdf","pdf");
  fitter2->fModel->PrintAllMarginalized("out_hist_plots.pdf");
  fitter2->fModel->WriteMarginalizedDistributions("output_hist_dist.root", "RECREATE");

  fitter2->fModel->WriteMarkovChain("output_hist_mcmc.root", "RECREATE");
  


  // lets make a counting analysis


  std::vector<std::pair<double,double>> range;
  range.push_back(std::make_pair(2580,2600));
  range.push_back(std::make_pair(2605,2625));
  range.push_back(std::make_pair(2630,2650));

  std::vector<double>probs{0,1,0};
  BatGraphFitter *fitter_count = new BatGraphFitter(fHist,range,probs,400);

  fitter_count->SetPrecison(3);
  TF1 *fInt=new TF1("fInt","[0]*x",2500,2700);
  fInt->SetParLimits(0,0,20);
  TCanvas *c3=new TCanvas("c3","c3");
  c3->cd();
  fitter_count->Fit(fInt);

  c3->Draw();
  c3->Print("fit.pdf)","pdf");
  fitter_count->fModel->PrintAllMarginalized("out_count_plots.pdf");
  fitter_count->fModel->WriteMarginalizedDistributions("output_count_dist.root", "RECREATE");
    
  fitter_count->fModel->WriteMarkovChain("output_count_mcmc.root", "RECREATE");




  app->Run();


}

#include "BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include "TRandom3.h"


void GenToy(TH1D*&h,TRandom3 *&rand,double Nb,double Elow,double Ehigh,double Ns,double Q,double dE)
{
  // generate a toy given a certain number of background counts (on avg) in range [Nb,Ns] and a certain number of signal counts at energy Q with resolution dE (sigma)

  int N=rand->Poisson(Nb);
  for (int i=1;i<h->GetNbinsX()+1;i++)
    h->SetBinContent(i,0);
  std::cout<<"Nb = "<<Nb<<" Ns = "<<Ns<<" N = "<<N<<std::endl;
  for (int i=0;i<N;i++)
    {
      h->Fill(rand->Uniform(100)+Elow);
    }
  for (int i=0;i<rand->Poisson(Ns);i++)
    {
      h->Fill(rand->Gaus(Q,dE));
    }
}


double scale_2_rate(double eff,double m,double eta,double W)
{
  // convert half-life to a number of events per year

  return log(2)*m*eff*eta*6.066e26/(W);
}

		  

double T12_2_mbb(double T,double M,double G)
{
  return 511*1e6*sqrt(1/(G*T))/(M*1.27*1.27);
}


int main(int argc,char**argv)
{
  // CUPID toys


  // SET THE PARAMETERS

  double Qbb=3034.4;
  double eff=0.9*0.79;
  double dE = 5/2.355;
  double T=10;
  double m = 472;
  double eta=0.95;
  double W=177.7;
  double b=1e-4;

  TString mode=argv[1];
  
  TRandom3 *rand=new TRandom3(0);
    
  TCanvas *can = new TCanvas();
  if (mode=="E")
    {
      can->Print("CUPID_sens_exclusion/toys_0_sig.pdf(","pdf");
    }
  else
    {
      can->Print("CUPID_sens_discovery/toys_0_sig.pdf(","pdf");
    }
  gStyle->SetOptStat(0);
  TH1D * h = new TH1D("h","h",100,2984,3084);
  double scale = T*log(2)*m*eff*eta*6.066e26/(W);

  double Nb=100*b*m*T;

  TF1 *fit = new TF1("f",Form("[0]*%f+(%f*[1])*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",m*T,scale,Qbb,dE,dE),2984,3084);
  TF1 *fit_bk_only = new TF1("f",Form("[0]*%f",m*T),2984,3084);

  fit->SetParLimits(0,0,1e-3);
  fit->SetParLimits(1,0,3e-26);
  fit_bk_only->SetParLimits(0,0,1e-3);
  fit->SetParNames("b","T_{1/2}");
  fit_bk_only->SetParNames("b");

  double G=15.92*pow(10,-15);
  double Nlow=3.90;
  double Nhigh=6.588;
  
  double mlow=18.4;
  double mmax=50;

  TH1D *hlimits= new TH1D("hlimits","hlimits",1000,0,pow(10,28));
  TH1D *hmodes= new TH1D("hmodes","hmodes",1000,0,pow(10,28));

  TH1D *hmlow= new TH1D("hmlow","hmlow",1000,0,100);
  TH1D *hmhigh= new TH1D("hmhigh","hmhigh",1000,0,100);
  TH1D *margdistro;
  TH1D *margdistrob;
  TH1D*margdistrobonly;
  TH1D *hbkg= new TH1D("hbkg","hbkg",100,0,5e-4);;
  TH1D *hbkg_only= new TH1D("hbkg_only","hbkg_only",100,0,5e-4);
  
  TH1D *hprob_SB= new TH1D("hprob_SB","hprob_SB",10000,0,100);
  TH1D *hprob_B= new TH1D("hprob_B","hprob_B",10000,0,100);

  TH2D *hprob_SB_2D= new TH2D("hprob_SB_2D","hprob_SB_2D",1000,0,3e-27,2000,0,100);
  TH2D *hprob_B_2D= new TH2D("hprob_B_2D","hprob_B_2D",1000,0,3e-27,2000,0,100);
  TH2D *hlimits_2D= new TH2D("hlimits_2D","hlimits_2D",1000,0,3e-27,2000,0,1e-26);
  TH2D *hmodes_2D= new TH2D("hmodes_2D","hmodes_2D",1000,0,3e-27,2000,0,1e-26);

  if (mode=="E")
    can->Print("CUPID_sens_exclusion/toys_0_sig.pdf(","pdf");
  else
    can->Print("CUPID_sens_discovery/toys_0_sig.pdf(","pdf");

  GenToy(h,rand,Nb,2984,3084,0,3034.4,dE);

  // make the fitter
  BatGraphFitter *fitter2= new BatGraphFitter(h,fit);
  BatGraphFitter *fitterbkg= new BatGraphFitter(h,fit_bk_only);



  // Loop over the toys
  //-----------------------------------------------------------------------------------------
  for (int i=0;i<10000;i++)
    {
      double Ns;
      double inV;
      std::cout<<" "<<std::endl;
      std::cout<<" "<<std::endl;
      std::cout<<" "<<std::endl;
      std::cout<<"Fitting toy "<<i<<std::endl;
      std::cout<<"-------------------------------------"<<std::endl;
      
      // Create the toy
      if (mode=="E")
	GenToy(h,rand,Nb,2984,3084,0,3034.4,dE);
      else
	{
	  inV=rand->Uniform(0,3e-27);
	  Ns=scale*inV;
	  GenToy(h,rand,Nb,2984,3084,Ns,3034.4,dE);
	  
	}




      
      if (mode=="E")
	h->SetTitle(Form("toy %i ; Energy [keV] ; counts/0.1 keV ; ",i));
      else
	h->SetTitle(Form("toy %i T_{1/2}=%f ; Energy [keV] ; counts/0.1 keV ; ",i,1/inV));

      h->Draw("E");

      fitterbkg->SetTH1(h);
      fitterbkg->fModel->SetNChains(2);
      fitterbkg->SetPrecison(1);
      fitterbkg->Fit();


      
      // set some parameters
      fitter2->SetTH1(h);
      fitter2->fModel->SetNChains(2);
      fitter2->SetPrecison(1);
      fitter2->Fit();
      
      
      
	  
      // draw the outputs
      can->Draw();
      if (mode=="E")
	{
	  can->Print("CUPID_sens_exclusion/toys_0_sig.pdf","pdf");
	}
      else
	{
	  can->Print("CUPID_sens_discovery/toys_0_sig.pdf","pdf");
	}


      margdistro = (TH1D*)fitter2->fModel->GetMarginalizedHistogram("T_{1/2}" );
      margdistro->SetTitle(Form("Posterior on T_{1/2}^{-1} for toy %i ; T_{1/2}^{-1} ; Probability [arb. units] ; ",i));
      margdistro->Draw();
      can->Draw();


      // get the limits                                                                                                                                                                                
      double x,q;
      q=0.9;
      margdistro->GetQuantiles(1,&x,&q);

      double mode_in=margdistro->GetBinCenter(margdistro->GetMaximumBin());
      
      double mlow =T12_2_mbb(1/x,Nlow,G);
      double mhigh=T12_2_mbb(1/x,Nhigh,G);
      
      std::cout<<"Limit = "<<1/x<<" years"<<std::endl;
      std::cout<<"Range = "<<mlow<<" "<<mhigh<<std::endl;


      
      if (mode=="E")
	{
	  can->Print("CUPID_sens_exclusion/toys_0_sig.pdf","pdf");


	  // save them
	  hlimits->Fill(1/x);
	  hmlow->Fill(mlow);
	  hmhigh->Fill(mhigh);
	      
	  margdistrob = (TH1D*)fitter2->fModel->GetMarginalizedHistogram("b" );
	  margdistrobonly = (TH1D*)fitterbkg->fModel->GetMarginalizedHistogram("b" );
	  
	  double b_est = margdistrob->GetBinCenter(margdistrob->GetMaximumBin());
	  double bonly_est = margdistrob->GetBinCenter(margdistrobonly->GetMaximumBin());
	  
	  hbkg->Fill(b_est);
	  hbkg_only->Fill(b_est);

	}
      else
	{
	  hlimits_2D->Fill(inV,x);
	  hmodes_2D->Fill(inV,mode_in);
	  can->Print("CUPID_sens_discovery/toys_0_sig.pdf","pdf");
	}

      // Get the evidences
      fitter2->fModel->SetRelativePrecision(1e-4);
      fitterbkg->fModel->SetRelativePrecision(1e-4);

      double evidence_SB = fitter2->fModel->Normalize();
      double evidence_B = fitterbkg->fModel->Normalize();
      if (mode=="E")
	{
	  hprob_SB->Fill(100*evidence_SB/(evidence_SB+evidence_B));
	  hprob_B->Fill(100*evidence_B/(evidence_SB+evidence_B));
	}
      else
	{
	  hprob_SB_2D->Fill(inV,100*evidence_SB/(evidence_SB+evidence_B));
	  hprob_B_2D->Fill(inV,100*evidence_B/(evidence_SB+evidence_B));
	}
      std::cout<<"evidences S+B, B = "<<evidence_SB<<" , "<<evidence_B<<std::endl;
      fitter2->fModel->ResetResults();
      fitterbkg->fModel->ResetResults();
	  
    }
      


  if (mode=="E")
    {

      // make some plots
      hlimits->SetTitle("Distribution of limits for CUPID baseline ; T_{1/2} [90 % c.i.] ; Number of toys ; ");
      hlimits->Draw();

	      
      
      can->Draw();
      can->Print("CUPID_sens_exclusion/toys_0_sig.pdf)","pdf");
      can->Print("CUPID_sens_exclusion/sensitivity.pdf");
      can->Print("CUPID_sens_exclusion/sens.C");

      
      hbkg->SetTitle("Mode background index ; b [counts/keV/kg/yr] ; Probability [arb. units] ;");
      TLegend *lb = new TLegend(0.7,0.7,0.9,0.9);
      lb->AddEntry(hbkg,"Signal + background");
      hbkg_only->SetLineColor(2);
      lb->AddEntry(hbkg_only,"Background only");
      hbkg->Draw();
      hbkg_only->Draw();
      lb->Draw();
      
      can->Draw();
      can->Print("CUPID_sens_exclusion/bkg.pdf");
      
      hmlow->SetTitle("Distribution of exclusion of m_{#beta#beta}; m_{#beta#beta}[ 90 % c.i.] ; Number of toys ; ");
      hmhigh->SetTitle("Distribution of exclusion of m_{#beta#beta} with largest NME; m_{#beta#beta}[ 90 % c.i.] ; Number of toys ; ");

      double maxi=1.2*hmhigh->GetMaximum();
      
      hmhigh->GetYaxis()->SetRangeUser(0,maxi);
	      
      TLegend * l = new TLegend(0.7,0.7,0.9,0.9);
      l->AddEntry(hmlow,"Smallest NMEs");
      l->AddEntry(hmhigh,"Largest NMEs");
      hmhigh->Draw();
      hmlow->SetLineColor(2);
      hmlow->Draw("same");
      TBox *b = new TBox(mlow,0,mmax,maxi);
      b->SetFillColorAlpha(9,0.3);
      l->AddEntry(b,"IO region");
      
      b->Draw();
      l->Draw();
      can->Draw();
      
      can->Print("CUPID_sens_exclusion/mhigh.pdf");
      can->Print("CUPID_sens_exclusion/mhigh.C");
      

      hprob_SB->SetTitle("Probability of models ; Probability(model|Data) [%] ; Number of toys ; ");
      hprob_B->SetLineColor(2);

      
      TLegend *le =new TLegend(0.7,0.7,0.9,0.9);
      le->AddEntry(hprob_SB,"Signal + Background");
      le->AddEntry(hprob_B,"Background");
      
      hprob_SB->Draw();
      hprob_B->Draw("same");
      can->Draw();
      le->Draw();
      can->Draw();
      can->Print("CUPID_sens_exclusion/evidences.C");
      
    }
  else
    {
      
      // discovery plots

      hlimits_2D->SetTitle("90 % limits on (T_{1/2})^{-1} as a function of injected signal ; Injected (T_{1/2}^{-1}) [yr^{-1}] ; T_{1/2}^{-1} 90 % [yr^{-1}]; ");
      hlimits_2D->Draw("colz");
      TF1 *fline = new TF1("fline","x",0,1e-26);
      fline->Draw("same");
      can->Draw();
      can->Print("CUPID_sens_discovery/limits.C");
      

      hmodes_2D->SetTitle("Mode of (T_{1/2})^{-1} as a function of injected signal ; Injected (T_{1/2})^{-1} [yr^{-1}] ; T_{1/2}^{-1} [yr^{-1}]; ");
      hmodes_2D->Draw("colz");
      //TF1 *fline = new TF1("fline","x",0,1e-26);
      fline->Draw("same");
      can->Draw();
      can->Print("CUPID_sens_discovery/modes.C");

      hprob_SB_2D->SetTitle("Probability of S+B model ; Injected Signal [yr^{-1}] ; Probability(model|Data) [%] ; ");
      hprob_B_2D->SetTitle("Probability of B model ; Injected Signal [yr^{-1}] ; Probability(model|Data) [%] ; ");

      hprob_SB_2D->Draw("colz");
      can->Draw();
      can->Print("CUPID_sens_discovery/SigBkg_prob.C");
      
      hprob_B_2D->Draw("colz");
      
      
      TF1 *fline2 = new TF1("fline2","0.27",0,1e-26);
      fline2->Draw("same");
      can->Draw();
      can->Print("CUPID_sens_discovery/Bkg_prob.C");
      

	    

      can->Print("CUPID_sens_discovery/toys_0_sig.pdf)","pdf");


    }
      
  
}


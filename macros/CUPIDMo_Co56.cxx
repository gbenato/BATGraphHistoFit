#include "TCut.h"
#include "../src/BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include "TH1D.h"

void CUPIDMo_Co56(TString path)
{


  TFile *fdata= new TFile(path);

  TTree * T = (TTree*)fdata->Get("signaltree");

  TCut cut_M1 = "RiseTimeCut==1&&abs(NormTVR)<10&&abs(NormTVL)<10&&abs(NormBaseline)<15&&abs(PCA1)<9&&Multiplicity==1&&Channel1!=3&&abs(LightDist)<4";
  TCut cut_M2 = "abs(PCA1)<23 &&Multiplicity==2 &&Channel1!=3 &&Channel2!=3 &&NormLightCut";



  TH1D * hM1 = new TH1D("hM1","hM1",4000/2.,0,4000);
  TH1D *hM1_fine = new TH1D("hM1_fine","hM1_fine",8000,0,4000);
  TH2D *hM2 = new TH2D("hM2","hM2",4000/2.,0,4000,4000/2.,0,4000);
  TH1D *hE1_fine = new TH1D("hE1_fine","hE1_fine",8000,0,4000);
  



  T->Draw("Ener1>>hM1",cut_M1);
  T->Draw("Ener1:Ener2>>hM2",cut_M2);
  T->Draw("Ener1>>hM1_fine",cut_M1);
  T->Draw("Ener1>>hE1_fine",cut_M2);

  hM1->SetTitle(" ; Energy [keV] ; counts/2 keV ; ");
  hM2->SetTitle(" ; Energy 1 [keV] ; Energy 2 [keV] ;");

  hM2->SetContour(1000);

  hM1->Draw();
  TPad *subpad = new TPad("subpad","",0.15,0.1,0.5,0.45);
  subpad->Draw();
  subpad->cd();
  TH1D* hM1_zoom = (TH1D*)hM1->Clone("hM1_zoom");
  hM1_zoom->GetXaxis()->SetRangeUser(2800,3600);
  hM1_zoom->Draw();
  gStyle->SetOptStat(0);
  TCanvas::MakeDefCanvas();
  gStyle->SetOptFit(1);
  hM2->Draw("colz");

  hM2->GetXaxis()->SetRangeUser(0,1500);
  hM2->GetYaxis()->SetRangeUser(0,3500);


  int counter=0;
  std::vector<double>peaks{846.8,1037.8,1238.3,1771.4,2598.5,3202.0,3253.5,3451.2};
  std::vector<bool> is_two{0,0,0,0,0,0,1,0};
  double peak2= 3273.1;
  TCanvas * can = new TCanvas("can","can");
  can->Print("peaks.pdf(","pdf");
  TGraphErrors * g_sigma = new TGraphErrors();
  TGraphErrors *g_bias = new TGraphErrors();
  for (int i=0;i<peaks.size();i++)
    {

      
      hM1_fine->GetXaxis()->SetRangeUser(peaks[i]-50,peaks[i]+50);
      TF1 * fit = new TF1("fit","[0]+[1]*(x-[2])+[3]*TMath::Gaus(x,[2],[4]/2.355)",peaks[i]-30,peaks[i]+30);


      fit->SetParameter(0,hM1_fine->GetBinContent(hM1_fine->FindBin(peaks[i]-30)));
      fit->SetParLimits(0,0,2*hM1_fine->GetBinContent(hM1_fine->FindBin(peaks[i]-30))+10);
      fit->SetParameter(1,0);
      fit->SetParLimits(1,-(hM1_fine->GetBinContent(hM1_fine->FindBin(peaks[i]-30))+10)/50.,(hM1_fine->GetBinContent(hM1_fine->FindBin(peaks[i]-30))+10)/50.);
      fit->SetParameter(2,peaks[i]);
      fit->SetParLimits(2,peaks[i]-4,peaks[i]+4);
      fit->SetParameter(3,hM1_fine->GetMaximum());
      fit->SetParLimits(3,0,2*hM1_fine->GetMaximum());
      fit->SetParameter(4,5);
      fit->SetParLimits(4,2,12);
      fit->SetNpx(10000);


      
      if (is_two[i]==1)
	{
	  fit = new TF1("fit","[0]+[1]*(x-[2])+[3]*TMath::Gaus(x,[2],[4]/2.355)+[5]*TMath::Gaus(x,[6],[7]/2.355)",peaks[i]-30,peaks[i]+50);

	  fit->SetParameter(0,hM1_fine->GetBinContent(hM1_fine->FindBin(peaks[i]-30))+2);
	  fit->SetParLimits(0,0,2*hM1_fine->GetBinContent(hM1_fine->FindBin(peaks[i]-30))+10);
	  fit->SetParameter(1,0);
	  fit->SetParLimits(1,-(hM1_fine->GetBinContent(hM1_fine->FindBin(peaks[i]-30))+10)/50.,(hM1_fine->GetBinContent(hM1_fine->FindBin(peaks[i]-30))+10)/50.);

	  fit->SetParameter(2,peaks[i]);
	 fit->SetParLimits(2,peaks[i]-4,peaks[i]+4);

	  fit->SetParameter(3,hM1_fine->GetMaximum());
	  fit->SetParLimits(3,0,2*hM1_fine->GetMaximum());

	  fit->SetParameter(4,5);
	  fit->SetParLimits(4,2,12);

	  fit->SetParameter(7,5);
	  fit->SetParLimits(7,2,12);

	  fit->SetParameter(6,peaks[i]+20);
	  fit->SetParLimits(6,peaks[i]+10,peaks[i]+30);
	  fit->SetParameter(5,hM1_fine->GetMaximum()/4.);
	  fit->SetParLimits(5,0,hM1_fine->GetMaximum());

	  fit->SetNpx(1000);
	}

      
      BatGraphFitter *fitter= new BatGraphFitter(hM1_fine,fit);
      //fitter->SetPrecison(3);
      fitter->Fit("R","",peaks[i]-40,peaks[i]+40);
      
      hM1_fine->Draw();
      can->Draw();
      fit->Draw("Csame");
      //  fitter->SetPrecison(3);
   
      can->Print("peaks.pdf","pdf");

      can->Print(Form("output/CUPIDMo/fit_%i.pdf",int(peaks[i])));
      can->Print(Form("output/CUPIDMo/fit_%i.C",int(peaks[i])));
     
      fitter->fModel->PrintAllMarginalized(Form("output/CUPIDMo/out_hist_plots_%i.pdf",int(peaks[i])));
      can->cd();
      TH1D *dist_sigma = (TH1D*)fitter->fModel->GetMarginalizedHistogram(4);
      TH1D *dist_mu = (TH1D*)fitter->fModel->GetMarginalizedHistogram(2);

      TF1 * fgauss_sigma = new TF1("fgauss_sigma","gaus",2,12);
      TF1 * fgauss_bias = new TF1("fgauss_bias","gaus",peaks[i]-4,peaks[i]+4);

      fgauss_sigma->SetNpx(1000);
      fgauss_bias->SetNpx(1000);


      dist_sigma->Draw();
      dist_sigma->Fit(fgauss_sigma);
      can->Draw();
      can->Print("peaks.pdf","pdf");

      dist_mu->Draw();
      dist_mu->Fit(fgauss_bias);
      can->Draw();
      can->Print("peaks.pdf","pdf");

      g_sigma->SetPoint(counter,peaks[i],fgauss_sigma->GetParameters()[1]);
      g_sigma->SetPointError(counter,0,fgauss_sigma->GetParameters()[2]);

      g_bias->SetPoint(counter,peaks[i],fgauss_bias->GetParameters()[1]-peaks[i]);
      g_bias->SetPointError(counter,0,fgauss_bias->GetParameters()[2]);

      counter++;
      if (is_two[i]==1)
	{
	  dist_sigma = (TH1D*)fitter->fModel->GetMarginalizedHistogram(7);
	  dist_mu = (TH1D*)fitter->fModel->GetMarginalizedHistogram(6);

	  fgauss_sigma = new TF1("fgauss_sigma","gaus",2,12);
	  fgauss_bias = new TF1("fgauss_bias","gaus",peak2-4,peak2+4);
	
	  fgauss_sigma->SetNpx(1000);
	  fgauss_bias->SetNpx(1000);
	  
	  
	  dist_sigma->Draw();
	  dist_sigma->Fit(fgauss_sigma);
	  can->Draw();
	  can->Print("peaks.pdf","pdf");
	  
	  dist_mu->Draw();
	  dist_mu->Fit(fgauss_bias);
	  can->Draw();
	  can->Print("peaks.pdf","pdf");
	  
	  g_sigma->SetPoint(counter,peak2,fgauss_sigma->GetParameters()[1]);
	  g_sigma->SetPointError(counter,0,fgauss_sigma->GetParameters()[2]);
	  
	  g_bias->SetPoint(counter,peak2,fgauss_bias->GetParameters()[1]-peak2);
	  g_bias->SetPointError(counter,0,fgauss_bias->GetParameters()[2]);
	  counter++;
	}
      
    }
  g_sigma->SetTitle(" ; Energy [keV]; #Delta E FWHM [keV] ; ");
  g_bias->SetTitle(" ; Energy [keV]; E_{fit} - E_{nom} [keV]; ");

  TF1 *f_sigma = new TF1("f_sigma","sqrt([0]*[0]+[1]*[1]*x)",0,4000);
  f_sigma->SetParLimits(0,0,5);
  f_sigma->SetParLimits(1,0,0.4);

  BatGraphFitter *fitter_sigma = new BatGraphFitter(g_sigma,f_sigma);
  fitter_sigma->SetQbb(3034.4);
  fitter_sigma->SetPrecison(3);
  fitter_sigma->SetGraphMaxMin(20,0);
  fitter_sigma->fModel->WriteMarkovChain(Form("output/CUPIDMo/sigma_mcmc.root",path.Data()), "RECREATE");	\
  can->cd();
  fitter_sigma->Fit("C");
  can->SaveAs("output/CUPIDMo/reso.C");

  can->Draw();
  
  can->Print("peaks.pdf","pdf");

  TF1 *f_bias= new TF1("f_bias","pol2",0,4000);

  g_bias->Fit(f_bias);
  
  f_bias->SetParLimits(0,-fabs(f_bias->GetParameters()[0])*2,2*fabs(f_bias->GetParameters()[0]));
  f_bias->SetParLimits(1,-fabs(f_bias->GetParameters()[1])*2,2*fabs(f_bias->GetParameters()[1]));
  f_bias->SetParLimits(2,-fabs(f_bias->GetParameters()[2])*2,2*fabs(f_bias->GetParameters()[2]));

  BatGraphFitter *fitter_bias = new BatGraphFitter(g_bias,f_bias);
  fitter_bias->SetPrecison(3);
  fitter_bias->SetQbb(3034.4);

  fitter_bias->SetGraphMaxMin(10,-10);
  fitter_bias->fModel->WriteMarkovChain(Form("output/CUPIDMo/bias_mcmc.root",path.Data()), "RECREATE");
  can->cd();

  fitter_bias->Fit("C");
  
  can->Draw();
  can->Print("peaks.pdf)","pdf");
  can->SaveAs("output/CUPIDMo/bias.C");



  // Fit the DE peaks

  double DE_1 =3253.5-1022;
  double DE_2 = 2598.5-1022;
  std::vector<double>DE{2598.5-1022,3253.5-1022};
  for (int i=0;i<DE.size();i++)
    {
      hE1_fine->GetXaxis()->SetRangeUser(DE[i]-50,DE[i]+50);
      TF1 * fit = new TF1("fit","[0]+[1]*(x-[2])+[3]*TMath::Gaus(x,[2],[4]/2.355)",DE[i]-30,DE[i]+30);


      fit->SetParameter(0,hE1_fine->GetBinContent(hE1_fine->FindBin(DE[i]-30)));
      fit->SetParLimits(0,0,2*hE1_fine->GetBinContent(hE1_fine->FindBin(DE[i]-30))+10);
      fit->SetParameter(1,0);
      fit->SetParLimits(1,-(hE1_fine->GetBinContent(hE1_fine->FindBin(DE[i]-30))+10)/50.,(hE1_fine->GetBinContent(hE1_fine->FindBin(DE[i]-30))+10)/50.);
      fit->SetParameter(2,DE[i]);
      fit->SetParLimits(2,DE[i]-4,DE[i]+4);
      fit->SetParameter(3,hE1_fine->GetMaximum());
      fit->SetParLimits(3,0,2*hE1_fine->GetMaximum());
      fit->SetParameter(4,5);
      fit->SetParLimits(4,2,12);
      fit->SetNpx(10000);


      BatGraphFitter *fitter= new BatGraphFitter(hE1_fine,fit);
      //fitter->SetPrecison(3);                                                                                                                                                                                                                                       
      fitter->Fit("R","",DE[i]-40,DE[i]+40);

      hE1_fine->Draw();
      can->Draw();
      fit->Draw("Csame");
      //  fitter->SetPrecison(3);                                                                                                                                                                                                                                     

      can->Print("peaks.pdf","pdf");

      can->Print(Form("output/CUPIDMo/fit_DE_%i.pdf",int(DE[i])));
      can->Print(Form("output/CUPIDMo/fit_DE_%i.C",int(DE[i])));

      fitter->fModel->PrintAllMarginalized(Form("output/CUPIDMo/out_hist_plots_%i.pdf",int(DE[i])));
      can->cd();
      TH1D *dist_sigma = (TH1D*)fitter->fModel->GetMarginalizedHistogram(4);
      TH1D *dist_mu = (TH1D*)fitter->fModel->GetMarginalizedHistogram(2);

      TF1 * fgauss_sigma = new TF1("fgauss_sigma","gaus",2,12);
      TF1 * fgauss_bias = new TF1("fgauss_bias","gaus",DE[i]-4,DE[i]+4);

      fgauss_sigma->SetNpx(1000);
      fgauss_bias->SetNpx(1000);


      dist_sigma->Draw();
      dist_sigma->Fit(fgauss_sigma);
      dist_sigma->SaveAs(Form("output/CUPIDMo/DE_%i_sigma.root",int(DE[i])));
      
      can->Draw();
      can->Print("peaks.pdf","pdf");

      dist_mu->Draw();
      dist_mu->SaveAs(Form("output/CUPIDMo/DE_%i_mean.root",int(DE[i])));

      dist_mu->Fit(fgauss_bias);
      can->Draw();
      can->Print("peaks.pdf","pdf");


  
    }






}
  


int main()
{
  CUPIDMo_Co56("inputs/cupid_mo/Calibration56Co_10ms_Reduced_ds8249.root");
}

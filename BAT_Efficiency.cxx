#include "TString.h"

#include "TFile.h"
#include "TH1D.h"
#include "BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"



std::vector<TH1D*>GetHistograms(TH1D* &h,TString title,std::vector<double>energies={1173,1333,1461,2615},double dE_center=10,double dE_side=30)
{
  std::vector<TH1D*>out;
  TCanvas *ci= new TCanvas();
  ci->Print(Form("%s.pdf(",title.Data()),"pdf");
  for (auto & E : energies)
    {

      /* CREATE THE CANVAS
	 -----------------------------------------------------
      */
      
      ci= new TCanvas(Form("c_%i",(int)E),Form("c_%i",(int)E));
      ci->cd();
      
      // PARAMETERS FOR THE FIT
      // ------------------------------------------------------
      std::vector<std::pair<double,double>> range;
      range.push_back(std::make_pair(E-dE_side,E-dE_center));
      range.push_back(std::make_pair(E-dE_center,E+dE_center));
      range.push_back(std::make_pair(E+dE_center,E+dE_side));
      std::vector<double>probs{0,1,0};


      h->GetXaxis()->SetRangeUser(E-dE_side-5,E+dE_side+5);


      // INITIAL PARS
      // ------------------------------------------------------
      double max=h->GetMaximum();
      int bin1 = h->FindBin(E-dE_side);
      int bin2 = h->FindBin(E-dE_center);
      double binwidth=h->GetBinWidth(bin1);

      double B = h->Integral(bin1,bin2)/(bin2-bin1);
      double Smax = max-B;
      double Sest = Smax*(sqrt(2*TMath::Pi()*5))/binwidth;

      // MAKE THE FITTER
      // ------------------------------------------------------
      BatGraphFitter *fitter_count_pass = new BatGraphFitter(h,range,probs,3*Sest);
      fitter_count_pass->SetPrecison(3);

      // MAKE FIT FUNCTIOn
      //-------------------------------------------------------
      TF1 *fInt=new TF1("fInt",Form("[0]*x+[1]*pow(x-%f,2)/2.",E),E-dE_side-5,E+dE_side+5);
      fInt->SetParLimits(0,0,3*B);
      fInt->SetParameter(1,0);
      fInt->SetParLimits(1,-B*0.5/30.,+B*0.5/30.);


      // RUN THE FIT AND DRAW
      //-------------------------------------------------------
      fitter_count_pass->Fit(fInt);
      ci->Draw();
      ci->Print(Form("%s.pdf",title.Data()),"pdf");
      

      // GET THE POSTERIOR
      // -------------------------------------------------------
      TH1D* margdistro = (TH1D*)fitter_count_pass->fModel->GetMarginalizedHistogram("S" );
      margdistro->SetTitle(Form("Peak at %f keV ; ; ; ",E));
      margdistro->Draw();
      
      ci->Draw();
      ci->Print(Form("%s.pdf",title.Data()),"pdf");
      fitter_count_pass->fModel->PrintAllMarginalized(Form("%s_%i_plots.pdf",title.Data(),(int)E));
      out.push_back(margdistro);
      
    }
  ci->Print(Form("%s.pdf)",title.Data()),"pdf");
  return out;
}

TH1D * toy_combine(TH1D*hpass,TH1D*hfail,double &mu,double &sigma,double E,int Ntoys=1e5)
{
  TH1D * hout = new TH1D("hout","hout",1200,0,1.2);
  for (int i=0;i<Ntoys;i++)
    {
      double Np=hpass->GetRandom();
      double Nf=hfail->GetRandom();

      hout->Fill(Np/(Np+Nf));

    }
  hout->GetXaxis()->SetRangeUser(0.8,1.1);
  hout->SetTitle(Form("Efficiency posterior for energy %i keV ; Efficiency ; Probability ; ",(int)E));
  TF1 *fGauss= new TF1("fGauss","gaus",0,1);
  hout->Fit(fGauss);

  mu=fGauss->GetParameters()[1];
  sigma=fGauss->GetParameters()[2];

  return hout;
}
      
  
int main(int argc, char **argv)
{

  // calculation of efficiency by peaks
  
  // please replace with a cfg file
  //std::map<TString,PeakInfo> peak_map;
  TApplication *app  =new TApplication("app",&argc,argv);

  int ds=3601;
  TFile *f = new TFile(Form("/home/tdixon/Downloads/m1_eff_histo_withtimecut_PCACut_directsum_ds%i.root",ds));

  

  TH1D *h = (TH1D*)f->Get(Form("h_ds%i",ds));
  TH1D *hb = (TH1D*)f->Get(Form("hb_ds%i",ds));
  TH1D *hm = (TH1D*)f->Get(Form("hm_ds%i",ds));
  TH1D *hmb = (TH1D*)f->Get(Form("hmbds%i",ds));
  TH1D *hpca = (TH1D*)f->Get(Form("hpca_ds%i",ds));
  TH1D *hm_fail = (TH1D*)f->Get(Form("hmbfds%i",ds));
  TH1D *hpca_fail = (TH1D*)f->Get(Form("hpcafds%i",ds));
  hpca->SetLineColor(1);
  
  hpca->Rebin(1/hpca->GetBinWidth(2));
  hpca_fail->Rebin(1/hpca_fail->GetBinWidth(2));
  std::vector<double> energies={1173,1333,1461,2615};

  std::vector<TH1D*>hpass = GetHistograms(hpca,"output/out_pca_pass");
  std::vector<TH1D*>hfail = GetHistograms(hpca_fail,"output/out_pca_fail");

  
  TGraphErrors *gerror = new TGraphErrors(); 
  TCanvas *ce = new TCanvas();
  ce->Print("output/eff.pdf(","pdf");
  for (int i=0;i<hpass.size();i++)
    {
      double mu,sigma;
      TH1D * hout=toy_combine(hpass[i],hfail[i],mu,sigma,energies[i]);
      gerror->SetPoint(i,energies[i],mu);
      gerror->SetPointError(i,0,sigma);
      hout->Draw();
      ce->Print("output/eff.pdf","pdf");
    }

  gerror->SetTitle(Form("PCA efficiency for ds %i ; Energy [keV] ; Efficiency ",ds));
  TF1 *fEff = new TF1("fEff","[0]+[1]*x/4000",0,4000);
  fEff->SetParameter(0,0.8);
  fEff->SetParLimits(0,0,1.2);
  fEff->SetParameter(1,0);
  fEff->SetParLimits(1,-0.5,0.5);
  BatGraphFitter *fitter = new BatGraphFitter(gerror);
  fitter->SetPrecison(3);
  fitter->SetQbb(2527);
  fitter->SetGraphMaxMin(1,0);

  fitter->Fit(fEff,"C");

  ce->Draw();
  ce->Print("output/eff.pdf)","pdf");

  ce->SaveAs("test.C");
  fitter->fModel->PrintAllMarginalized("output/eff_fit.pdf");
  fitter->fModel->WriteMarginalizedDistributions("output_eff.root", "RECREATE");


      
  int done;
  std::cin>>done;
  
}
      

#include "../src/BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TLine.h"
#include "TMath.h"
#include "TF1.h"
#include "TApplication.h"
#include "TRandom3.h"
#include "TLatex.h"
#include <getopt.h>
#include <chrono>
// Method to get the integral assuming                                                                                                                                                                    //-------------------------------------------------------                                                                                                                                                  

double GetIntegral(double low,double high,TH1D*&h)
{


  int lowbin = h->FindBin(low);
  int highbin=h->FindBin(high);

  double integral = (h->Integral(lowbin+1,highbin)+h->GetBinContent(lowbin)*(double)rand()/RAND_MAX)/h->Integral();

  return integral;





}
double p_2_z(double p)
{
  // convert p-value to significance
  
  return TMath::ErfcInverse(2*p)*sqrt(2);
}

void GetLimit(double &sens,TH1D *&hsens,TH1D* &hinv_sens,double t,bool & foundsens,TH1D *&test_stat,double p_old,int j,double Nsignals,double maxr)
{

  double p =GetIntegral(t,150,test_stat);
  if (!foundsens&& p<0.1)
    {
      double diff=p_old-p;
      double frac= (p_old-0.10)/diff;
      sens = ((j-1)/(double)Nsignals)*(maxr)+(maxr)*frac/((double)Nsignals);
 
      hsens->Fill((1e-27)*1/sens);
      hinv_sens->Fill(1e27*sens);
      foundsens=1;
    }









}

// METHOD TO MAKE A PROFILED LIKELIHOOD
//-----------------------------------------------------------------------------------------------------------
TGraphErrors *profile(TH2D* &h,TString axis,double logLL)
{

  TH1D *hout;
  int counter=0;
  double Nmax =h->GetMaximum();
  TGraphErrors *g = new TGraphErrors();
  if (axis=="X")
    {
      hout->Clear();
      hout->Reset();
      
      TH1D *h1;
      TH1D *h2;
      for (int i=1;i<h->GetNbinsX();i++)
	{
	  
	  h1=(TH1D*)h->ProjectionY("h1",i,i+1);
	  h2=(TH1D*)h->ProjectionX("h2",i,i+1);

	  if (h1->GetMaximum()>0)
	    {

	      g->SetPoint(counter,h2->GetBinCenter(h2->FindBin(i)),-log(h1->GetMaximum())+log(Nmax)-logLL);
	      g->SetPointError(counter,0,1/sqrt(h1->GetMaximum()));
	      counter++;
	      
	    }
	}
    }
  
  else
    {

      TH1D *h1;
      TH1D *h2;
      for (int i=1;i<h->GetNbinsY();i++)
        {
          h1=(TH1D*)h->ProjectionX("h1",i,i+1);
	  h2=(TH1D*)h->ProjectionY("h2",i,i+1);

	  if (h1->GetMaximum()>0)
	    {

	      g->SetPoint(counter,h2->GetBinCenter(h2->FindBin(i)),-log(h1->GetMaximum())+log(Nmax)-logLL);
	      g->SetPointError(counter,0,1/sqrt(h1->GetMaximum()));
	      counter++;

            }
	  
	  
	}
      
    }
  return g;


  
}	  




double  GenToy(TH1D*&h,std::vector<double>&vec,TRandom3 *&rand,double Nb,double Elow,double Ehigh,double Ns,double Q,double dE)
{
  // generate a toy given a certain number of background counts (on avg) in range [Nb,Ns] and a certain number of signal counts at energy Q with resolution dE (sigma)
  
  int N=rand->Poisson(Nb);
  double Na=rand->Poisson(Ns);
  for (int i=1;i<h->GetNbinsX()+1;i++)
    h->SetBinContent(i,0);
  vec.clear();
  for (int i=0;i<N;i++)
    {
      double E=rand->Uniform(Ehigh-Elow)+Elow;
  
      h->Fill(E);
      vec.push_back(E);
   }
  for (int i=0;i<Na;i++)
    {
      double E=rand->Gaus(Q,dE);
      h->Fill(E);
      vec.push_back(E);
    }
  return Na;
}
double GenFancyToy(TH1D*&h,std::vector<double>&vec,TRandom3*&rand,TF1 *&fb,double Elow,double Ehigh,double Ns,double Q,double dE)
{
  int N = rand->Poisson(fb->Integral(Elow,Ehigh));

  double Na=rand->Poisson(Ns);
  vec.clear();
  for (int i=1;i<h->GetNbinsX()+1;i++)
    h->SetBinContent(i,0);
  
  for (int i=0;i<N;i++)
    {
      double E=fb->GetRandom(Elow,Ehigh);
      h->Fill(E);
      vec.push_back(E);
    }
  for (int i=0;i<Na;i++)
    {
      double E=rand->Gaus(Q,dE);
      h->Fill(E);
      vec.push_back(E);
    }
  return Na;
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



void CreatePValuePlot(std::vector<TH1D*>p_value,std::vector<TH1D*>test_stat,std::vector<TH1D*>test_stat_bkg,int Nsignals,double cut,TString label,TCanvas *&can,TString name,double maxr,TGraphErrors *&prob_plot)
{
  // fill graphs and compute exclusion sensitivity                                                                                                                                                     
  //----------------------------------------------------------------------------------------------
  
  TGraphAsymmErrors* p1=new TGraphAsymmErrors();
  TGraphAsymmErrors* p2=new TGraphAsymmErrors();
  TGraph *gmed=new TGraph();


  prob_plot=new TGraphErrors();
  std::cout<<"p - value has size "<<p_value.size()<<std::endl;
  std::cout<<"test stat / test stat bkg has size"<<test_stat.size()<<" "<<test_stat_bkg.size()<<" Nsig = "<<Nsignals<<std::endl;

  for (int j=0;j<Nsignals;j++)
    {
      double invT12=(j/(double)Nsignals)*maxr;



      p_value[j]->SetTitle(Form("Distribution of p_S for T_{1/2}^{-1}  =  %E yr^{-1} ; p ; counts",invT12));
      p_value[j]->Draw();
      can->Draw();
      can->Print(Form("output/CUPID_sens/%s/%s_p_value_%i.pdf",name.Data(),label.Data(),j));

      double Np=p_value[j]->Integral(0,p_value[j]->FindBin(cut));
      double Nt= p_value[j]->Integral();

      double p = Np/Nt;
      double q=1-p;
      double ep=sqrt(p*q/Nt);

      prob_plot->SetPoint(j,invT12,p);
      prob_plot->SetPointError(j,0,ep);
      
      

      q;
      q=(1-0.9544)/2.;


      double m2sig;


      p_value[j]->GetQuantiles(1,&m2sig,&q);

      q=(1-0.682)/2.;
      double m1sig;
      p_value[j]->GetQuantiles(1,&m1sig,&q);

      double med;
      q=0.5;
      p_value[j]->GetQuantiles(1,&med,&q);


      double p1sig;
      q=1-((1-0.682)/2.);
      p_value[j]->GetQuantiles(1,&p1sig,&q);

      double p2sig;
      q=1-((1-0.9544)/2.);

      p_value[j]->GetQuantiles(1,&p2sig,&q);

      
      if (med<cut)
        {
          std::cout<<"Sensitivity = "<<j<<" or "<<Form("%E",1/invT12)<<std::endl;
        }
      std::cout<<"quantiles = "<<m2sig<<" "<<m1sig<<" "<<med<<" "<<p1sig<<" "<<p2sig<<std::endl;

      p1->SetPoint(j,invT12,med);
      p1->SetPointError(j,0,0,med-m1sig,p1sig-med);

      p2->SetPoint(j,invT12,med);
      p2->SetPointError(j,0,0,med-m2sig,p2sig-med);

      gmed->SetPoint(j,invT12,med);

      test_stat[j]->Draw();
      test_stat_bkg[j]->SetLineColor(2);
      test_stat_bkg[j]->Draw("same");
      can->Draw();
      can->Print(Form("output/CUPID_sens/%s/%s_test_stat_comp_%i.C",name.Data(),label.Data(),j));

    }


  // create some more plots                           
  // ------------------------------------------------------------------------------                                                                                                                                                                                                                                                                                                                                                                                       

  p1->SetLineColorAlpha(3,0.2);
  p2->SetLineColorAlpha(kOrange-3,0.2);
  p1->SetFillColorAlpha(3,0.7);
  p2->SetFillColorAlpha(kOrange-3,0.7);

  gmed->SetTitle(" ; (T_{1/2})^{-1} [yr^{-1}] ; p-value ; ");

  TLine * l =new TLine(0,cut,maxr,cut);
  l->SetLineColor(2);
  
  gmed->Draw("AL");
  p2->Draw("3same");
  p1->Draw("3same");
  gmed->Draw("Lsame");
  l->Draw("same");
  can->Draw();
  can->Print(Form("output/CUPID_sens/%s/%s_exc.C",name.Data(),label.Data()));
  can->Print(Form("output/CUPID_sens/%s/%s_toys_bkg_0_sig.pdf)",name.Data(),label.Data()),"pdf");




}

void TestStatDist(double bkg,TString name,bool fancy,double maxR,int Nsignals,int Ntoys,int idx)
{
  // SET THE PARAMETERS
  // -----------------------------------------------------------------------------------------------------------------------------------------
  
  double Qbb=3034.4;
  double eff=0.9*0.79;
  double dE = 5/2.355;
  double T=10;
  double m = 472;
  double eta=0.95;
  double W=177.7;
  double b=bkg;
  double slope=-6e-7;
  double Elow=2960;
  double Ehigh=3100;
  double minR=-0.5*maxR;
  
  std::vector<double>peaks{2978.9,3000.0,3053.9,3081.8};
  std::vector<double>norm{21e-5,9.4e-5,27.4e-5,6.1e-5};

  // Build the background prediction
  //------------------------------------------------------------------------------------------------------------------------------------------

  TF1 *fb = new TF1("fb","[0]*([1]+(x-3034.4)*[2]+[3]*TMath::Gaus(x,[4],[5]/2.355)+[6]*TMath::Gaus(x,[7],[5]/2.355)+[8]*TMath::Gaus(x,[9],[5]/2.355)+[10]*TMath::Gaus(x,[11],[5]/2.355))",Elow,Ehigh);

  fb->SetParameter(0,m*T);
  fb->SetParameter(1,b);
  fb->SetParameter(2,slope);
  fb->SetParameter(5,5);
  fb->SetParameter(3,norm[0]);
  fb->SetParameter(4,peaks[0]);
  fb->SetParameter(6,norm[1]);
  fb->SetParameter(7,peaks[1]);
  fb->SetParameter(8,norm[2]);
  fb->SetParameter(9,peaks[2]);
  fb->SetParameter(10,norm[3]);
  fb->SetParameter(11,peaks[3]);

  
  fb->SetNpx(1000);

  // Create some objects
  //------------------------------------------------------------------------------------------------------------------------------------------
  TRandom3 *rand=new TRandom3(0);
  TCanvas *can = new TCanvas();
  
  can->Print(Form("output/CUPID_sens/%s/toys_%i_sig.pdf(",name.Data(),idx),"pdf");
  gStyle->SetOptStat(0);
  TH1D * h = new TH1D("h","h",Ehigh-Elow,Elow,Ehigh);

  double scale = T*log(2)*m*eff*eta*6.066e26/(W);
  std::vector<double>vec;
  double Nb=(Ehigh-Elow)*b*m*T;
  double invT12=(idx/(double)Nsignals)*maxR;


  // CREATE FIT FUNCTIONS
  //-------------------------------------------------------------------------------------------------------------------------------------------
  TF1 *fit;
  TF1 *fit_fix;
  TF1 *fit_bk_only;

  TF1 *fit_int;
  TF1 *fit_fix_int;
  TF1 *fit_bk_only_int;

  
  int Rpar;
  if (!fancy)
    {
      fit = new TF1("f",Form("[0]*%f+(%f*[1])*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",m*T,scale,Qbb,dE,dE),Elow,Ehigh);
      fit_int =new TF1("f_int",Form("[0]*%f+%f*[1]",m*T*(Ehigh-Elow),scale),Elow,Ehigh);
    
      fit_fix = new TF1("fit_fix",Form("[0]*%f+%f*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",m*T,0.,Qbb,dE,dE),Elow,Ehigh);
      fit_fix_int = new TF1("fit_fix_int",Form("[0]*%f+%f",m*T*(Ehigh-Elow),0.),Elow,Ehigh);
      
      
      fit_bk_only = new TF1("fbk",Form("[0]*%f",m*T),Elow,Ehigh);
      fit_bk_only_int = new TF1("fbk_int",Form("[0]*%f",m*T*(Ehigh-Elow)),Elow,Ehigh);

      fit->SetParLimits(0,0,1e-3);
      fit->SetParLimits(1,minR,2*maxR);
      fit_fix->SetParLimits(0,0,2e-3);
      fit->SetParNames("b","T_{1/2}^{-1}");
      fit_fix->SetParNames("b");
      fit_bk_only->SetParNames("b");
      fit_bk_only->SetParLimits(0,0,2e-3);
      Rpar=1;
    }
  else
    {

      
      fit = new TF1("f",Form("%f*([0]+(x-3034.4)*[1]+[2]*TMath::Gaus(x,%f,%f)+[3]*TMath::Gaus(x,%f,%f)+[4]*TMath::Gaus(x,%f,%f)+[5]*TMath::Gaus(x,%f,%f))+(%f*[6])*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",
			     m*T,peaks[0],dE,peaks[1],dE,peaks[2],dE,peaks[3],dE,scale,Qbb,dE,dE),Elow,Ehigh);

      fit_int = new TF1("fit_int",Form("%f*([0]*%f+[2]*sqrt(2*TMath::Pi())*%f+[3]*sqrt(2*TMath::Pi())*%f+[4]*sqrt(2*TMath::Pi())*%f+[5]*sqrt(2*TMath::Pi())*%f)+(%f*[6])",
				       m*T,Ehigh-Elow,dE,dE,dE,dE,scale),Elow,Ehigh);

      
      fit_fix = new TF1("fit_fix",Form("%f*([0]+(x-3034.4)*[1]+[2]*TMath::Gaus(x,%f,%f)+[3]*TMath::Gaus(x,%f,%f)+[4]*TMath::Gaus(x,%f,%f)+[5]*TMath::Gaus(x,%f,%f))+%f*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",
				       m*T,peaks[0],dE,peaks[1],dE,peaks[2],dE,peaks[3],dE,0.,Qbb,dE,dE),
			Elow,Ehigh);

      fit_fix_int = new TF1("fit_fix_int",Form("%f*([0]*%f+[2]*sqrt(2*TMath::Pi())*%f+[3]*sqrt(2*TMath::Pi())*%f+[4]*sqrt(2*TMath::Pi())*%f+[5]*sqrt(2*TMath::Pi())*%f)+%f",
					       m*T,Ehigh-Elow,dE,dE,dE,dE,0.),
                        Elow,Ehigh);



      
      fit_bk_only = new TF1("fbk",Form("%f*([0]+(x-3034.4)*[1]+[2]*TMath::Gaus(x,%f,%f)+[3]*TMath::Gaus(x,%f,%f)+[4]*TMath::Gaus(x,%f,%f)+[5]*TMath::Gaus(x,%f,%f))",
				       m*T,peaks[0],dE,peaks[1],dE,peaks[2],dE,peaks[3],dE),
			    Elow,Ehigh);
      
       fit_bk_only_int = new TF1("fbk_int",Form("%f*([0]*%f+[2]*sqrt(2*TMath::Pi())*%f+[3]*sqrt(2*TMath::Pi())*%f+[4]*sqrt(2*TMath::Pi())*%f+[5]*sqrt(2*TMath::Pi())*%f)",
						m*T,Ehigh-Elow,dE,dE,dE,dE),
                            Elow,Ehigh);

      fit->SetParLimits(0,0,1e-3);
      fit->SetParLimits(1,-1e-5,1e-5);
      fit->SetParLimits(2,0,10e-4);
      fit->SetParLimits(3,0,10e-4);
      fit->SetParLimits(4,0,10e-4);
      fit->SetParLimits(5,0,10e-4);
      fit->SetParLimits(6,minR,2*maxR);
      
      fit_fix->SetParLimits(0,0,1e-3);
      fit_fix->SetParLimits(1,-1e-5,1e-5);
      fit_fix->SetParLimits(2,0,10e-4);
      fit_fix->SetParLimits(3,0,10e-4);
      fit_fix->SetParLimits(4,0,10e-4);
      fit_fix->SetParLimits(5,0,10e-4);
      fit->SetParNames("b","s","n1","n2","n3","n4","T_{1/2}^{-1}");
      fit_fix->SetParNames("b","s","n1","n2","n3","n4");
      fit_bk_only->SetParNames("b","s","n1","n2","n3","n4");

      fit_bk_only->SetParLimits(0,0,1e-3);
      fit_bk_only->SetParLimits(1,-1e-5,1e-5);
      fit_bk_only->SetParLimits(2,0,5e-4);
      fit_bk_only->SetParLimits(3,0,5e-4);
      fit_bk_only->SetParLimits(4,0,5e-4);
      fit_bk_only->SetParLimits(5,0,5e-4);
      Rpar=6;
    }


  // Some parameters
  //----------------------------------------------------------------------

  can->Print(Form("output/CUPID_sens/%s/toys_%i_sig.pdf(",name.Data(),idx),"pdf");
  double Ns,Nin;
  Ns=scale*invT12;



  if (!fancy)
    GenToy(h,vec,rand,Nb,Elow,Ehigh,0,3034.4,dE);
  else
    GenFancyToy(h,vec,rand,fb,Elow,Ehigh,0,3034.4,dE);

  h->Draw();
  fb->Draw("Csame");
  can->Draw();
  
  

  gErrorIgnoreLevel=kFatal;

  // make the fitter                                                                                                                                                                                        //------------------------------------------------------------
  BatGraphFitter *fitter2= new BatGraphFitter(h,fit,fit_int,vec);
  BatGraphFitter *fitterbkg= new BatGraphFitter(h,fit_fix,fit_fix_int,vec);
  BatGraphFitter *fitterbkg2= new BatGraphFitter(h,fit_bk_only,fit_bk_only_int,vec);


  // CONTAINERS FOR OUTPUT
  //-----------------------------------------------------------------------------------------
  TFile *file_out = new TFile(Form("output/CUPID_sens/%s/test_stat/test_stat_%i_%.2E.root",name.Data(),idx,invT12),"RECREATE");
  file_out->cd();
  TTree *Tout = new TTree("test_stat","test_stat");

  double fix_rate;
  double fix_bkg;
  Tout->Branch("fix_rate",&fix_rate,"fix_rate/D");
  Tout->Branch("fix_bkg",&fix_bkg,"fix_bkg/D");

  double fix_peak1;
  double fix_peak2;
  double fix_peak3;
  double fix_peak4;
  double fix_slope;

  double float_peak1;
  double float_peak2;
  double float_peak3;
  double float_peak4;
  double float_slope;

  
  double bkg_only_peak1;
  double bkg_only_peak2;
  double bkg_only_peak3;
  double bkg_only_peak4;
  double bkg_only_slope;

  if (fancy)
    {
      Tout->Branch("fix_peak1",&fix_peak1,"fix_peak1/D");
      Tout->Branch("fix_peak2",&fix_peak2,"fix_peak2/D");
      Tout->Branch("fix_peak3",&fix_peak3,"fix_peak3/D");
      Tout->Branch("fix_peak4",&fix_peak4,"fix_peak4/D");
      Tout->Branch("fix_slope",&fix_slope,"fix_slope/D");

    }

  double float_rate;
  double float_bkg;
  Tout->Branch("float_rate",&float_rate,"float_rate/D");
  Tout->Branch("float_bkg",&float_bkg,"float_bkg/D");

  if (fancy)
    {
     
      Tout->Branch("float_peak1",&float_peak1,"float_peak1/D");
      Tout->Branch("float_peak2",&float_peak2,"float_peak2/D");
      Tout->Branch("float_peak3",&float_peak3,"float_peak3/D");
      Tout->Branch("float_peak4",&float_peak4,"float_peak4/D");
      Tout->Branch("float_slope",&float_slope,"float_slope/D");

    }
  

  double bkg_only_rate;
  double bkg_only_bkg;
  Tout->Branch("bkg_only_rate",&bkg_only_rate,"bkg_only_rate/D");
  Tout->Branch("bkg_only_bkg",&bkg_only_bkg,"bkg_only_bkg/D");


  if (fancy)
    {
     
      Tout->Branch("bkg_only_peak1",&bkg_only_peak1,"bkg_only_peak1/D");
      Tout->Branch("bkg_only_peak2",&bkg_only_peak2,"bkg_only_peak2/D");
      Tout->Branch("bkg_only_peak3",&bkg_only_peak3,"bkg_only_peak3/D");
      Tout->Branch("bkg_only_peak4",&bkg_only_peak4,"bkg_only_peak4/D");
      Tout->Branch("bkg_only_slope",&bkg_only_slope,"bkg_only_slope/D");

    }

  double test;
  double test_zero;
  double rate=invT12;
  int index=idx;

  Tout->Branch("test_stat",&test,"test_stat/D");
  Tout->Branch("test_stat_zero",&test_zero,"test_zero/D");
  Tout->Branch("rate",&rate,"rate/D");
  Tout->Branch("idx",&index,"idx/I");
  
  
  



  TF1 *fit2_fix=nullptr;
  TF1 *fit2_fix_int=nullptr;



  // Loop over the toys                                                                                                                                                                                   
  //-----------------------------------------------------------------------------------------

  // Loop over signal strength      
  for (int i=0;i<Ntoys;i++)
    {
      
      bool quiet;
      if (i%100000==0)
	quiet=0;
      else
	quiet=1;
      
      if (!quiet)
	{
	  
	  std::cout<<"Making fit "<<i<<" with "<<invT12<<std::endl;
	  std::cout<<"**************************************"<<std::endl;
	}

	  

      // generate toy
      //----------------------------------
      
      double Ns,Nin;
      Ns=scale*invT12;
	  
	  
      if (!fancy)
	GenToy(h,vec,rand,Nb,Elow,Ehigh,Ns,3034.4,dE);
      else
	    GenFancyToy(h,vec,rand,fb,Elow,Ehigh,Ns,3034.4,dE);
      if (!quiet)
	{
	  std::cout<<"N events = "<<vec.size()<<" histo int = "<<h->Integral()<<std::endl;
	}
      h->SetTitle(Form("toy %i T_{1/2}=%E ; Energy [keV] ; counts/0.1 keV ; ",i,invT12));
      h->Draw("E");



      // Fit with S+B model
      //-----------------------------------
      fitter2->SetTH1(h);
      fitter2->SetVector(vec);
      fitter2->fModel->SetNChains(2);
      fitter2->SetPrecison(1);
      
      if (!quiet)
	    {
	      std::cout<<" "<<std::endl;
	      std::cout<<"making fit to toy with "<<Form("%E",double(1/invT12))<<" with S+B model"<<std::endl;
	    }
      fitter2->Fit(" "," ",Elow,Ehigh,0,quiet);

      std::vector<double>modes =fitter2->fModel->GetBestFitParameters();

      if (!fancy)
	{
	  float_bkg=modes[0];
	  float_rate=modes[1];
	  float_peak1=-1;
	  float_peak2=-1;
	  float_peak3=-1;
	  float_peak4=-1;
	  float_slope=-1;

	}
      else
	{
	       

	  float_bkg=modes[0];
          float_rate=modes[6];
          float_peak1=modes[2];
          float_peak2=modes[3];
          float_peak3=modes[4];
          float_peak4=modes[5];
          float_slope=modes[1];
	}
	  


      double logLL=fitter2->fModel->LogLikelihood(modes);

      // fit with 0 signal                                                                                                                                                                                      //---------------------------------                                                                                                                                                              
      
      fitterbkg2->SetTH1(h);
      fitterbkg2->SetVector(vec);
      
      fitterbkg2->fModel->SetNChains(2);
      fitterbkg2->SetPrecison(1);
      if (!quiet)
	{
	  std::cout<<" "<<std::endl;
	  std::cout<<"making fit to toy with "<<Form("%E",1/invT12)<<" with B model"<<std::endl;
	}

      fitterbkg2->Fit(" "," ",Elow,Ehigh,0,quiet);
      std::vector<double>modes_fix2 =fitterbkg2->fModel->GetBestFitParameters();
      double logLL_0=fitterbkg2->fModel->LogLikelihood(modes_fix2);

      if (!fancy)
        {
          bkg_only_bkg=modes_fix2[0];
          bkg_only_rate=0;
          bkg_only_peak1=-1;
          bkg_only_peak2=-1;
          bkg_only_peak3=-1;
          bkg_only_peak4=-1;
          bkg_only_slope=-1;

        }
      else
        {


          bkg_only_bkg=modes_fix2[0];
          bkg_only_rate=0;
          bkg_only_peak1=modes_fix2[2];
          bkg_only_peak2=modes_fix2[3];
          bkg_only_peak3=modes_fix2[4];
          bkg_only_peak4=modes_fix2[5];
          bkg_only_slope=modes_fix2[1];
        }


      if (modes[Rpar]<0)
	{
	  logLL=logLL_0;
	}

      
      // Fit fixing signal
      //-----------------------------------
      
	  
      if (fancy==0)
	{
	  fit2_fix = new TF1("fit2_fix",Form("[0]*%f+(%f)*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",m*T,scale*invT12,Qbb,dE,dE),Elow,Ehigh);
	  fit2_fix_int = new TF1("fit2_fix_int",Form("[0]*%f+(%f)",m*T*(Ehigh-Elow),scale*invT12),Elow,Ehigh);
	  
	  fit2_fix->SetParLimits(0,0,1e-3);
	  fit2_fix->SetParNames("b");
	}
      else
	{
	  fit2_fix = new TF1("fit2_fix",Form("%f*([0]+(x-3034.4)*[1]+[2]*TMath::Gaus(x,%f,%f)+[3]*TMath::Gaus(x,%f,%f)+[4]*TMath::Gaus(x,%f,%f)+[5]*TMath::Gaus(x,%f,%f))+%f*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",
						m*T,peaks[0],dE,peaks[1],dE,peaks[2],dE,peaks[3],dE,scale*invT12,Qbb,dE,dE),

			     Elow,Ehigh);

	  fit2_fix_int = new TF1("fit2_fix_int",Form("%f*(%f*[0]+[2]*(sqrt(2*TMath::Pi())*%f)+[3]*(sqrt(2*TMath::Pi())*%f)+[4]*(sqrt(2*TMath::Pi())*%f)+[5]*(sqrt(2*TMath::Pi())*%f))+%f",
							 m*T,Ehigh-Elow,dE,dE,dE,dE,scale*invT12),

                                 Elow,Ehigh);
	  
	  fit2_fix->SetParLimits(0,0,1e-3);
	  fit2_fix->SetParLimits(1,-1e-5,1e-5);
	  fit2_fix->SetParLimits(2,0,10e-4);
	  fit2_fix->SetParLimits(3,0,10e-4);
	  fit2_fix->SetParLimits(4,0,10e-4);
	  fit2_fix->SetParLimits(5,0,10e-4);
	  fit2_fix->SetParNames("b","s","n1","n2","n3","n4");
	}




      fitterbkg->ResetTF1(fit2_fix);
      fitterbkg->ResetTF1Int(fit2_fix_int);
      fitterbkg->SetTH1(h);
      fitterbkg->SetVector(vec);
      fitterbkg->fModel->SetNChains(2);
      fitterbkg->SetPrecison(1);
      if (!quiet)
	{
	  std::cout<<" "<<std::endl;
	  std::cout<<"making fit to toy with "<<Form("%E",1/invT12)<<" with S fixed model"<<std::endl;
	}
      
      fitterbkg->Fit(" "," ",Elow,Ehigh,0,quiet);
      std::vector<double>modes_fix =fitterbkg->fModel->GetBestFitParameters();

      if (!fancy)
        {
          fix_bkg=modes_fix[0];
          fix_rate=invT12;
          fix_peak1=-1;
          fix_peak2=-1;
          fix_peak3=-1;
          fix_peak4=-1;
          fix_slope=-1;

        }
      else
        {


          fix_bkg=modes_fix[0];
          fix_rate=invT12;
          fix_peak1=modes_fix[2];
          fix_peak2=modes_fix[3];
          fix_peak3=modes_fix[4];
          fix_peak4=modes_fix[5];
          fix_slope=modes_fix[1];
        }





      // get output - save test stats
      // -----------------------------------
      double logLL_fix=fitterbkg->fModel->LogLikelihood(modes_fix);
      
      double t =-2*(logLL_fix-logLL);
	  


      test=t;


      // Svae the t0
      double t0=-2*(logLL_0-logLL);
      test_zero=t0;
     
	  
	  
      if (!quiet)
	{
	  std::cout<<"maxlogL = "<<logLL<<" maxlogL(S= "<<invT12<<" ) = "<<logLL_fix<<" maxlogL(S= 0 ) = "<<logLL_0<<std::endl;
	  std::cout<<"test stat (S) = "<<-2*(logLL_fix-logLL)<<std::endl;
	  std::cout<<"t0            = "<<t0<<std::endl;
	  
	  std::cout<<" "<<std::endl;
	  std::cout<<" "<<std::endl;
	  std::cout<<" "<<std::endl;
	}
      
      
	  
      // reset results
      fitter2->fModel->ResetResults();
      fitterbkg->fModel->ResetResults();
      fitterbkg2->fModel->ResetResults();
      


	  

      //save plots
      //-----------------------------------
      if (!quiet)
	{
	  h->Draw();
	  fit->Draw("Csame");
	  can->Draw();
	  can->SaveAs(Form("output/CUPID_sens/%s/toys_%i_%i_sig.C",name.Data(),i,idx));
	  fit2_fix->SetLineColor(3);
	  fit2_fix->Draw("Csame");
	  fit_bk_only->SetLineColor(4);
	  fit_bk_only->Draw("Csame");
	  can->Draw();
	      
	  can->Print(Form("output/CUPID_sens/%s/toys_%i_sig.pdf",name.Data(),idx),"pdf");
	}
      
      Tout->Fill();
      delete fit2_fix_int;
      delete fit2_fix;
    }
  
  Tout->Write();
  file_out->Close();
  can->Print(Form("output/CUPID_sens/%s/toys_%i_sig.pdf)",name.Data(),idx),"pdf");
  

}
  /*
  
  
  // Now generate a new array of bkg only toys
  //---------------------------------------------------------------------------
  
  can->Print(Form("output/CUPID_sens/%s/toys_bkg_0_sig.pdf(",name.Data()),"pdf");

  TH1D *hinvsens = new TH1D("hinvsens","hinvsens",1000,0,40);

  TH1D *hsens = new TH1D("hsens","hsens",1000,0,40);

  // Loop over the toys
  //---------------------------

  for (int i=0;i<Ntoys;i++)
    {

      bool quiet=1;
      if (i%100000==0)
	{
	  quiet=0;
	}

      // Generate the toy
      //---------------------------------------------------------------
      
      double Ns,Nin;
      Ns=0;

      if (!fancy)
	Nin=GenToy(h,vec,rand,Nb,Elow,Ehigh,Ns,3034.4,dE);
      else
	Nin=GenFancyToy(h,vec,rand,fb,Elow,Ehigh,Ns,3034.4,dE);

      h->SetTitle(Form("toy %i  ; Energy [keV] ; counts/0.1 keV ; ",i));
      h->Draw("E");


      //Fit with S+B model
      // ------------------------------------------------------------
      
      fitter2->SetTH1(h);
      fitter2->SetVector(vec);
      fitter2->fModel->SetNChains(2);
      fitter2->SetPrecison(1);
      if (!quiet)
	{
	  std::cout<<"making fit to toy with 0 signal to S+B model"<<std::endl;
	}

      fitter2->Fit(" "," ",Elow,Ehigh,0,quiet);
      std::vector<double>modes =fitter2->fModel->GetBestFitParameters();
      double logLL=fitter2->fModel->LogLikelihood(modes);
      bool foundsens=0;
      double p=0;

      TF1 *fit3_fix=nullptr;
      TF1 *fit3_fix_int=nullptr;


      // fit with S+B model fixing S
      // --------------------------------------------------------------
      for (int j=0;j<Nsignals;j++)
	{

	  double invT12=(j/(double)Nsignals)*maxR;


	  if (!quiet)
	    {
	      std::cout<<"Making fit "<<i<<" with "<<invT12<<std::endl;
	      std::cout<<"**************************************"<<std::endl;
	    }

	  // create fit function + fit
	  // ------------------------------------------------------------


	  if (fancy==0)
            {
              fit3_fix = new TF1("fit3_fix",Form("[0]*%f+(%f)*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",m*T,scale*invT12,Qbb,dE,dE),Elow,Ehigh);
	      fit3_fix_int = new TF1("fit3_fix_int",Form("[0]*%f+(%f)",m*T*(Ehigh-Elow),scale*invT12),Elow,Ehigh);

	      fit3_fix->SetParLimits(0,0,1e-3);
              fit3_fix->SetParNames("b");
            }
          else
            {
              fit3_fix = new TF1("fit3_fix",Form("%f*([0]+(x-3034.4)*[1]+[2]*TMath::Gaus(x,%f,%f)+[3]*TMath::Gaus(x,%f,%f)+[4]*TMath::Gaus(x,%f,%f)+[5]*TMath::Gaus(x,%f,%f))+%f*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",
                                                m*T,peaks[0],dE,peaks[1],dE,peaks[2],dE,peaks[3],dE,scale*invT12,Qbb,dE,dE),
                                 Elow,Ehigh);
	      fit3_fix_int = new TF1("fit3_fix_int",Form("%f*(%f*[0]*[2]*(sqrt(2*TMath::Pi())*%f)+[3]*(sqrt(2*TMath::Pi())*%f)+[4]*(sqrt(2*TMath::Pi())*%f)+[5]*(sqrt(2*TMath::Pi())*%f))+%f",
							 m*T,Ehigh-Elow,dE,dE,dE,dE,scale*invT12),
                                 Elow,Ehigh);
	      
              fit3_fix->SetParLimits(0,0,1e-3);
              fit3_fix->SetParLimits(1,-1e-5,1e-5);
              fit3_fix->SetParLimits(2,0,10e-4);
              fit3_fix->SetParLimits(3,0,10e-4);
              fit3_fix->SetParLimits(4,0,10e-4);
              fit3_fix->SetParLimits(5,0,10e-4);
              fit3_fix->SetParNames("b","s","n1","n2","n3","n4");
      	    }

	  fitterbkg->ResetTF1(fit3_fix);
	  fitterbkg->ResetTF1Int(fit3_fix_int);
	  fitterbkg->SetTH1(h);
	  fitterbkg->SetVector(vec);
          fitterbkg->fModel->SetNChains(2);
          fitterbkg->SetPrecison(1);
	  if (!quiet)
            {
              std::cout<<"making fit to toy with 0 signal to model with "<<Form("%E",1/invT12)<<std::endl;
            }

          fitterbkg->Fit(" "," ",Elow,Ehigh,0,quiet);



	  // get test stat
	  //----------------------------------------------------------------
	  
          std::vector<double>modes_fix =fitterbkg->fModel->GetBestFitParameters();
	  
	  double logLL_fix=fitterbkg->fModel->LogLikelihood(modes_fix);

	   if (j==0 &&modes[Rpar]<0)
            {
	      logLL=logLL_fix;
	    }
				

	  double t=-2*(logLL_fix-logLL);
	  test_stat_bkg[j]->Fill(t);
	 

	  // get the p-value
	  //--------------------------------------------------------------------------
	  double sens;
	  GetLimit(sens,hsens,hinvsens,t,foundsens,test_stat[j],p,j,(double)Nsignals,maxR);
	  
	  
	 
	  p = GetIntegral(t,150,test_stat[j]); 
	  p_value[j]->Fill(p);
	 


	  if (!quiet)
	    {
	      std::cout<<" "<<std::endl;
	      std::cout<<" "<<std::endl;
	    }
	 fitter2->fModel->ResetResults();
	 fitterbkg->fModel->ResetResults();


	 

	 // make plots
	 //-----------------------------------------------------------------------
	 if (!quiet)
	   {
	     h->SetTitle(Form("toy %i T_{1/2}=%E ; Energy [keV] ; counts/0.1 keV ; ",i,invT12));
	     h->Draw();

	     for (int m;m<modes.size();m++)
	       {
		 fit->SetParameter(m,modes[m]);
	       }
	     fit3_fix->SetLineColor(3);

	     for (int m=0;m<modes_fix.size();m++)
	       {
		 fit3_fix->SetParameter(m,modes_fix[m]);
	       }
	     fit3_fix->Draw("Csame");
	     fit->Draw("Csame");
	     can->Draw();
	     can->Print(Form("output/CUPID_sens/%s/toys_bkg_0_sig.pdf",name.Data()),"pdf");
	   }
	 delete fit3_fix_int;
	 delete fit3_fix;
	}
    }
    hsens->SaveAs(Form("output/CUPID_sens/%s/sens.C",name.Data()));
    hinvsens->SaveAs(Form("output/CUPID_sens/%s/invsens.C",name.Data()));

    TGraphErrors *gex;
    TGraphErrors *gd;
    CreatePValuePlot(p_value,test_stat,test_stat_bkg,Nsignals,0.1,"exclude",can,name,maxR,gex);
    CreatePValuePlot(p_b,test_stat_0,test_stat_fit_0,Nsignals,0.14e-2,"discover",can,name,maxR,gd);

    gex->SetTitle("Probability of exclusion ; (T_{1/2})^{-1} ; Prob  ");
    gd->SetTitle("Probability of discovery ; (T_{1/2})^{-1} ; Prob  ");
	
    gex->SaveAs(Form("output/CUPID_sens/%s/ex_prob.C",name.Data()));
    gd->SaveAs(Form("output/CUPID_sens/%s/disc_prob.C",name.Data()));


  

    }*/

      
      









int BayesianLimits(int Ntoys,double b,TString name,double maxR,bool full_bkg,int index_low,int index_high,int group_index)
{
  // CUPID toys

  gErrorIgnoreLevel=kFatal;

  // SET THE PARAMETERS
  // ---------------------------------------------------------------------------
  
  double Qbb=3034.4;
  double eff=0.9*0.79;
  double dE = 5/2.355;
  double T=10;
  double m = 472;
  double eta=0.95;
  double W=177.7;
  //double b=bkg;
  double slope=-6e-7;
  double G=15.92*pow(10,-15);
  double Nlow=3.90;
  double Nhigh=6.588;
  double mlow=18.4;
  double mmax=50;
  bool fancy=full_bkg;

  double Elow=2960;
  double Ehigh=3100;

  std::vector<double>peaks{2978.9,3000.0,3053.9,3081.8};
  std::vector<double>norm{21e-5,9.4e-5,27.4e-5,6.1e-5};
  TTree *Tt  = new TTree("output","output");
  double limit;
  double mode;
  double bkg;
  double limit_mlow;
  double limit_mhigh;
  Tt->Branch("limit",&limit,"limit/D");
  Tt->Branch("bkg",&bkg,"bkg/D");
  Tt->Branch("limit_mlow",&limit_mlow,"limit_mlow/D");
  Tt->Branch("limit_mhigh",&limit_mhigh,"limit_mhigh/D");
  Tt->Branch("mode",&mode,"mode/D");
  
  // Build the background prediction                                                                                                                                                                      
  //------------------------------------------------------------------------------------------------------------------------------------------                                                            

  TF1 *fb = new TF1("fb","[0]*([1]+(x-3034.4)*[2]+[3]*TMath::Gaus(x,[4],[5])+[6]*TMath::Gaus(x,[7],[5])+[8]*TMath::Gaus(x,[9],[5])+[10]*TMath::Gaus(x,[11],[5]))",Elow,Ehigh);

  fb->SetParameter(0,m*T);
  fb->SetParameter(1,b);
  fb->SetParameter(2,slope);
  fb->SetParameter(5,5);
  fb->SetParameter(3,norm[0]);
  fb->SetParameter(4,peaks[0]);
  fb->SetParameter(6,norm[1]);
  fb->SetParameter(7,peaks[1]);
  fb->SetParameter(8,norm[2]);
  fb->SetParameter(9,peaks[2]);
  fb->SetParameter(10,norm[3]);
  fb->SetParameter(11,peaks[3]);

  

  fb->SetNpx(1000);




  TRandom3 *rand=new TRandom3(0);
  TCanvas *can = new TCanvas();
  can->Print(Form("output/CUPID_sens/%s/toys_0_sig.pdf(",name.Data()),"pdf");
  gStyle->SetOptStat(0);

  std::vector<double>vec;
  TH1D * h = new TH1D("h","h",Ehigh-Elow,Elow,Ehigh);
  double scale = T*log(2)*m*eff*eta*6.066e26/(W);
  double Nb=(Ehigh-Elow)*b*m*T;

  

  // Create the fits
  // ---------------------------------------------------------------------------
  TF1 *fit;
  TF1 *fit_bk_only;
  TF1 *fit_int;
  TF1 *fit_bk_only_int;
  if (!fancy)
    {
      fit = new TF1("f",Form("[0]*%f+(%f*[1])*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",m*T,scale,Qbb,dE,dE),Elow,Ehigh);
      fit_int = new TF1("f_int",Form("[0]*%f+(%f*[1])",m*T*(Ehigh-Elow),scale),Elow,Ehigh);

      fit_bk_only = new TF1("fbk",Form("[0]*%f",m*T),Elow,Ehigh);
      fit_bk_only_int=new TF1("fbk_int",Form("[0]*%f",m*T*(Ehigh-Elow)),Elow,Ehigh);
      fit->SetParLimits(0,0,1e-3);
      fit->SetParLimits(1,0,2*maxR);
      fit->SetParNames("b","T_{1/2}^{-1}");
      fit_bk_only->SetParNames("b");
      fit_bk_only->SetParLimits(0,0,2e-3);

    }
  else
    {
      fit = new TF1("f",Form("%f*([0]+(x-3034.4)*[1]+[2]*TMath::Gaus(x,%f,%f)+[3]*TMath::Gaus(x,%f,%f)+[4]*TMath::Gaus(x,%f,%f)+[5]*TMath::Gaus(x,%f,%f))+(%f*[6])*TMath::Gaus(x,%f,%f)/(sqrt(2*TMath::Pi())*%f)",
                             m*T,peaks[0],dE,peaks[1],dE,peaks[2],dE,peaks[3],dE,scale,Qbb,dE,dE),Elow,Ehigh);

      fit_int = new TF1("f_int",Form("%f*(%f*[0]+[2]*(sqrt(2*TMath::Pi())*%f)+[3]*(sqrt(2*TMath::Pi())*%f)+[4]*(sqrt(2*TMath::Pi())*%f)+[5]*(sqrt(2*TMath::Pi())*%f))+(%f*[6])",
				 m*T,Ehigh-Elow,dE,dE,dE,dE,scale),Elow,Ehigh);

      
      fit_bk_only = new TF1("fbk",Form("%f*([0]+(x-3034.4)*[1]+[2]*TMath::Gaus(x,%f,%f)+[3]*TMath::Gaus(x,%f,%f)+[4]*TMath::Gaus(x,%f,%f)+[5]*TMath::Gaus(x,%f,%f))",
                                       m*T,peaks[0],dE,peaks[1],dE,peaks[2],dE,peaks[3],dE),
                            Elow,Ehigh);
      fit_bk_only_int = new TF1("fbk_int",Form("%f*(%f*[0]+[2]*(sqrt(2*TMath::Pi())*%f)+[3]*(sqrt(2*TMath::Pi())*%f)+[4]*(sqrt(2*TMath::Pi())*%f)+[5]*(sqrt(2*TMath::Pi())*%f))",
					       m*T,Ehigh-Elow,dE,dE,dE,dE),
				Elow,Ehigh);
      
      
      fit->SetParLimits(0,0,1e-3);
      fit->SetParLimits(1,-20e-6,20e-6);
      fit->SetParLimits(2,0,10e-4);
      fit->SetParLimits(3,0,10e-4);
      fit->SetParLimits(4,0,10e-4);
      fit->SetParLimits(5,0,10e-4);
      fit->SetParLimits(6,0,2*maxR);

      fit->SetParNames("b","s","n1","n2","n3","n4","T_{1/2}^{-1}");
      fit_bk_only->SetParNames("b","s","n1","n2","n3","n4");

      fit_bk_only->SetParLimits(0,0,1e-3);
      fit_bk_only->SetParLimits(1,-20e-6,20e-6);
      fit_bk_only->SetParLimits(2,0,10e-4);
      fit_bk_only->SetParLimits(3,0,10e-4);
      fit_bk_only->SetParLimits(4,0,10e-4);
      fit_bk_only->SetParLimits(5,0,10e-4);

    }





  // Create all the objects
  //-----------------------------------------------------------------------------
  
  can->Print(Form("output/CUPID_sens/%s/toys_%i_sig.pdf(",name.Data(),group_index),"pdf");

  
  if (!fancy)
    GenToy(h,vec,rand,Nb,Elow,Ehigh,0,3034.4,dE);
  else
    GenFancyToy(h,vec,rand,fb,Elow,Ehigh,0,3034.4,dE);

  // make the fitter
  // -----------------------------------------------------------------------------------
  BatGraphFitter *fitter2= new BatGraphFitter(h,fit,fit_int,vec);
  BatGraphFitter *fitterbkg= new BatGraphFitter(h,fit_bk_only,fit_bk_only_int,vec);

  TH1D *margdistro;
  TH1D *margdistrob;

  TLatex *tlat =new TLatex();
	
  // Loop over the toys
  //-----------------------------------------------------------------------------------------
  for (int i=index_low;i<index_high;i++)
    {
      double Ns;
      double inV;
      if (i%100==0)
	{
	  std::cout<<" "<<std::endl;
	  std::cout<<" "<<std::endl;
	  std::cout<<" "<<std::endl;
	  std::cout<<"Fitting toy "<<i<<std::endl;
	  std::cout<<"-------------------------------------"<<std::endl;
	}
      double Nin;
      double Nout;

      bool zero_sig=0;
      
      // Create the toy
      inV=0;
	  

      Ns=scale*inV;
      
      
      if (!fancy)
	GenToy(h,vec,rand,Nb,Elow,Ehigh,Ns,3034.4,dE);
      else
	GenFancyToy(h,vec,rand,fb,Elow,Ehigh,Ns,3034.4,dE);


      
      h->SetTitle(Form("toy %i ; Energy [keV] ; counts/0.1 keV ; ",i));
      
      h->Draw("E");
      bool quiet=!(i%100==0);

      /*
      fitterbkg->SetTH1(h);
      fitterbkg->SetVector(vec);
      fitterbkg->fModel->SetNChains(2);
      fitterbkg->SetPrecison(2);
      fitterbkg->Fit(" "," ",Elow,Ehigh,1,quiet);
      
      */

  
      fitter2->SetTH1(h);
      fitter2->SetVector(vec);
      fitter2->fModel->SetNChains(2);
      fitter2->SetPrecison(2);
      fitter2->Fit(" "," ",Elow,Ehigh,1,quiet);

      h->Draw("HIST");
      
     fitter2->fModel->fTF1->Draw("Csame");
          // draw the outputs
      can->Draw();
      if (i%100==0)
	{
	  can->Print(Form("output/CUPID_sens/%s/toys_%i_sig.pdf",name.Data(),group_index),"pdf");
	  can->SaveAs(Form("output/CUPID_sens/%s/data_%i.C",name.Data(),i));

	}

      margdistro = (TH1D*)fitter2->fModel->GetMarginalizedHistogram("T_{1/2}^{-1}" );
      margdistro->SetTitle(Form("Posterior on T_{1/2}^{-1} for toy %i ; T_{1/2}^{-1} ; Probability [arb. units] ; ",i));
      
      margdistro->Draw();
      if (i%100==0)
	margdistro->SaveAs(Form("output/CUPID_sens/%s/limit_%i.root",name.Data(),i));

      // get the limits                                                                                                                                                                                    
      double x,q;
      q=0.9;
      margdistro->GetQuantiles(1,&x,&q);
      
      mode=margdistro->GetBinCenter(margdistro->GetMaximumBin());
      q=0.5;
      
      limit=1/x;
      
      limit_mlow =T12_2_mbb(1/x,Nlow,G);
      limit_mhigh=T12_2_mbb(1/x,Nhigh,G);

      margdistrob = (TH1D*)fitter2->fModel->GetMarginalizedHistogram("b" );

      bkg=margdistro->GetBinCenter(margdistrob->GetMaximumBin());

      Tt->Fill();

      can->Draw();
      


      fitter2->fModel->ResetResults();
      // fitterbkg->fModel->ResetResults();
	  
    }

  Tt->SaveAs(Form("output/CUPID_sens/%s/limits/limits_%i.root",name.Data(),group_index));

  can->Print(Form("output/CUPID_sens/%s/toys_%i_sig.pdf)",name.Data(),group_index),"pdf");

	     

      
  return 1;
}

void Usage()
{
  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  if( std::getenv("USER") != NULL )
    std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
  else
    std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
  std::cout<<"./CUPID_Sens"<<std::endl;
  std::cout<<"options "<<std::endl;
  std::cout<<"------------------------------------------------------------"<<std::endl;
  std::cout<<"-f (or --fit-type)        [either B for bayesian or F for frequentst] (default: B)"<<std::endl;
  std::cout<<"-n (or --number-toys)     [number of toys]                            (default: 10000)] "<<std::endl;
  std::cout<<"Next 2 options for bayesian analysis"<<std::endl;
  std::cout<<"-F (or --first-index)     [first toy]                            (default: 0)] "<<std::endl;
  std::cout<<"-L (or --last-index)     [last toy]                            (default: 10000)] "<<std::endl;
  std::cout<<"-G (or --group-index)     [index for sets]                            (default: 0)] "<<std::endl;
  std::cout<<"-t (or --test-stat-idx)   [which index for signal strength to generate - freq analysis] (default 0)] "<<std::endl;
  std::cout<<"-s (or --number-signals)  [number of signals for profile likelihood - freq analysis]  (default 250)] "<<std::endl;
  std::cout<<"-b (or --bkg-index)       [background index]                          (default: 1e-4)"<<std::endl;
  std::cout<<"-f (or --full-bkg)        [ful background? 0 or 1]                    (default: 0)"<<std::endl;
  std::cout<<"-n (or --name)            [name for output]                           (default bayesian_baseline_simple_bkg)"<<std::endl;
  std::cout<<"-r (or --max-rate)        [maximum rate:                              (default 1e-26 yr-1) "<<std::endl;
  std::cout<<"-h (or --help)              [this help]"<<std::endl;

}

int main(int argc,char **argv)
 {

   bool bayesian_fit=1;
   int Ntoys=10000;
   int Nsignals=250;
   int test_stat_idx=0;
   double bkg_index=1e-4;
   bool full_bkg=0;
   TString name="bayesian_baseline_simple_bkg";
   double max_rate=1e-26;
   int first,last,group;
   first=0;
   last=Ntoys;
   group=0;

   {
     static struct option long_options[] = {
					    { "fit-type",        required_argument,  nullptr, 'm' },
					    { "number-toys",       required_argument,  nullptr,'n'},
					    { "number-signals",required_argument,nullptr,'s'},
					    {"first-index",required_argument,nullptr,'F'},
					    {"last-index",required_argument,nullptr,'L'},
					    {"group-index",required_argument,nullptr,'G'},
					    {"test_stat_idx",required_argument,nullptr,'t'},
					    { "bkg-index",required_argument, nullptr,'b'},
					    { "full-bkg",      required_argument,nullptr,'f'  },
					    { "name",         required_argument,  nullptr,'e'},
					    { "max-rate",    required_argument, nullptr,'r'},
					    { "help",              no_argument,        nullptr,'h'},
					    {nullptr, 0, nullptr, 0}
  };
     const char* const short_options = "m:n:b:f:e:r:s:F:G:L:t:h";
     
     
     int c;
  
     while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 )
       {
	 switch (c)
	   {
	   case 'm':
	     {
	       std::cout<<optarg<<std::endl;
	       if (std::string(optarg)=="B")
		 bayesian_fit=1;
	       else if (std::string(optarg)=="F")
		 bayesian_fit=0;
	       else
		 {
		   std::cout<<"Error fit type must be either 'B' or 'F'"<<std::endl;
		   return -1;
		 }
	       break;
	     }
	   case 'n':
	     {
	       Ntoys = atoi(optarg);
	       break;
	     }
	   case 'F':
             {
               first = atoi(optarg);
               break;
             }
	   case 't':{
	     test_stat_idx=atoi(optarg);
	     break;
	   }
	   case 'L':
             {
               last = atoi(optarg);
               break;
             }
	   case 'G':
             {
               group= atoi(optarg);
               break;
             }

	   case 's':
	     {
	       Nsignals=atoi(optarg);
	       break;
	     }
	   case 'b':
	     {
	       bkg_index=atof(optarg);
	       break;
	     }
	     
	   case 'f':
	     {
	       full_bkg=atoi(optarg);
	       break;
	     }
	   case 'e':
	     {
	       name=optarg;
	       break;
	     }
      
	   case 'r':
	     {
	       max_rate=atof(optarg);
	       break;
	       
	     }

	   case'h':
	     {
	       Usage();
	       return 0;
	     }
	   default: {
	     exit(1);
	   }
	     
	   }
    

    
       }
   }

   std::cout<<"Running fit with"<<std::endl;
   std::cout<<"Ntoys = "<<Ntoys<<std::endl;
   std::cout<<"bkg index = "<<bkg_index<<std::endl;
   std::cout<<"name = "<<name<<std::endl;
   std::cout<<"test_stat_idx = "<<test_stat_idx<<std::endl;
   std::cout<<"full bkg = "<<full_bkg<<std::endl;
   if (bayesian_fit==1)
     BayesianLimits(Ntoys,bkg_index,name,max_rate,full_bkg,first,last,group);
   else
     TestStatDist(bkg_index,name,full_bkg,max_rate,Nsignals,Ntoys,test_stat_idx);
   
       

   return 1;
 }



// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "BAT_GraphFit.h"

#include <BAT/BCMath.h>
#include <exception>
#include "TLegend.h"

void BAT_GraphFit::PlotCI(TGraphAsymmErrors *g1,TGraphAsymmErrors *g2,TGraphAsymmErrors *g3)
{
  g3->SetTitle("Test fit ; Energy [keV] ; Efficiency ; ");
  g3->Draw("3A");
  g2->Draw("3same");
  g1->Draw("3same");

  fGraph->Draw("PEsame");

  fTF1->Draw("Csame");
  TLegend * l = new TLegend(0.7,0.7,0.9,0.9);
  l->AddEntry(fGraph,"Data","PE");
  l->AddEntry(fTF1,"Fit","L");
  l->AddEntry(g1,"1 #sigma","L");
  l->AddEntry(g2,"2 #sigma","L");
  l->AddEntry(g3,"3 #sigma","L");

  l->Draw();
}
void BAT_GraphFit::Plot(TString goption)
{
  if (goption=="" &&(fMode=="B"||fMode=="G"))
    {
      goption="APE";
    }
  else if (goption=="" &&(fMode=="P"||fMode=="C"))
    {
      goption="HIST";
    }

  // for the counting

  TH1D *hSide ;
  TH1D *hSignal ;

  
  if (fMode=="C")
    {
      hSide= (TH1D*)fTH1->Clone("hSide");
      hSignal =      (TH1D*)fTH1->Clone("hSignal");

      for (int i=1;i<fTH1->GetNbinsX()+1;i++)
	{
	  double E=fTH1->GetBinLowEdge(i);
	  
	  if (E<fCountingRange[0].first)
	    {
	      hSide->SetBinContent(i,0);
	      hSignal->SetBinContent(i,0);
	    }
	  else if(E<fCountingRange[0].second)
	    {
	      hSignal->SetBinContent(i,0);
	    }
	  else if (E<fCountingRange[1].first)
	    {
	      hSide->SetBinContent(i,0);
	      hSignal->SetBinContent(i,0);
	    }
	  else if (E<fCountingRange[1].second)
	    {
	      hSide->SetBinContent(i,0);
	    }
	  else if (E<fCountingRange[2].first)
	    {
	      hSignal->SetBinContent(i,0);
	      hSide->SetBinContent(i,0);
	    }
	  else if (E<fCountingRange[2].second)
	    {
	      hSignal->SetBinContent(i,0);
	    }
	  else
	    {
	      hSignal->SetBinContent(i,0);
	      hSide->SetBinContent(i,0);
	    }

	  TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
	  l->AddEntry(hSignal,"Signal region","F");
	  l->AddEntry(hSide,"Sideband region","F");
	  l->AddEntry(fTF1,"Background model");
	  hSignal->SetFillColor(8);
	  hSide->SetFillColor(9);
	  gStyle->SetOptStat(0);
	  fTH1->Draw("HIST");
	  hSignal->Draw("HISTsame");
	  hSide->Draw("HISTsame");
	  l->Draw();
      
	}
    }
	  
  if (fMode=="G")
    fGraph->Draw(goption);
  if (fMode=="P")
    fTH1->Draw(goption);
  if (fMode=="B")
    fGraph->Draw(goption);


 
  // draw derivative for CDF mode and function of PDF, nothing for counting
  if (fMode=="C"||(fType=="C" &&fMode=="P"))
    fTF1->DrawDerivative("Csame");
  else if (fMode!="C")
    fTF1->Draw("Csame");
}

void BAT_GraphFit::SetGraph(TGraphAsymmErrors *&g,double max,double min)
{
  fGraphMaximum=max;
  fGraphMinimum=min;
  fGraph=g;
  AddObservable("fQ",fMinObs,fMaxObs,"f(Q_{#beta#beta})","[]");

  GetObservables().SetNBins(1000);
  
}
void BAT_GraphFit::SetHisto(TH1D*&h,TString type)
{
  fType=type;
  fTH1=h;

 
}



void BAT_GraphFit::SetCountingPars(  std::vector<std::pair<double,double>>range,std::vector<double> prob)
{
  fCountingRange=range;
  fCountingProb=prob;
  fCountingN.resize(3);
  // lets enusre the range match the bin edges - we make a choice

  // add a check the 3 intervals dont overlap
  for (int i=0;i<3;i++)
    {
      
      double E1= fCountingRange.at(i).first;
      double E2=fCountingRange.at(i).second;
      
      int bin1=fTH1->FindBin(E1);

      int bin2=fTH1->FindBin(E2);

      double E1h=fTH1->GetBinLowEdge(bin1);
      double E2h=fTH1->GetBinLowEdge(bin2);

      if (E1!=E1h)
	{
	  // edit the edges
	  fCountingRange.at(i).first=E1h;
	}
      if (E2!=E2h)
	{
	  fCountingRange.at(i).second=E2h;
	}

      fCountingN.at(i)=fTH1->Integral(bin1,bin2-1);
      
    }

  fFitLow=fCountingRange[0].first;
  fFitHigh=fCountingRange[2].second;
  
  
}

	  

	  
  
  



void BAT_GraphFit::SetBinomGraph(TGraph *&gTrial,TGraph *&gSuccess)
{
  fGraphBinomTrial =gTrial;
  fGraphBinomSuccess = gSuccess;

  // some checks
  for (int i=0;i<gTrial->GetN();i++)
    {
      double x1= gTrial->GetX()[i];
      double x2=gSuccess->GetX()[i];
      double x3=fGraph->GetX()[i];
      if (x1!=x2 or x1!=x3)
	{
	  std::cout<<"ERROR for Binomal fittting the X-axis on the 3 graphs arent the same"<<std::endl;
	}
    }

}
// ---------------------------------------------------------
BAT_GraphFit::BAT_GraphFit(const std::string& name,TF1 *f,bool verbose,TString mode,double Sm,double Qbb)
    : BCModel(name)
{


  
  // get the information on the parameters from the TF1
  
  fCountingSmax=Sm;
  fVerbose=verbose;
  fMode=mode;
  // Default range from the TF1
  fQbb=Qbb;
  this->SetTF1(f);
}

void BAT_GraphFit::AddParameters()
{
  
  // for the counting instead the fit range comes from fCountingRanges

  if (fMode!="C")
    {
    

      if (fVerbose)
	{
	  std::cout<<"RUNNING FIT WITH RANGE: "<<fFitLow<<" ' "<<fFitHigh<<std::endl;
	  std::cout<<"***************************************************"<<std::endl;
	}
    }
  for (int N=0;N<fNpar;N++)
    {

      // get the name
      TString name=fTF1->GetParName(N);

      // get the range
      double low, high;
      fTF1->GetParLimits(N,low,high);

      if (low==0 &&high==0)
	{
	  std::cout<<Form("WARNING - Parameter %i %s has low = high = 0, set par limits !!!",N,name.Data())<<std::endl;
	  throw;
	}
      AddParameter(name.Data(),low,high,name.Data(),"[]");
    }
  // add signal for counting analysis
  if (fMode=="C")
    {
      AddParameter("S",0,fCountingSmax,"S","[]");
    }

  // how to set priors
  GetParameters().SetPriorConstantAll();
  GetParameters().SetNBins(300);
  
}






void BAT_GraphFit::CalculateObservables(const std::vector<double>&pars)
{
  if (fMode=="G")
    {
      for (int n=0;n<fNpar;n++)
	{
	  fTF1->SetParameter(n,pars[n]);
	}
      
      GetObservable(0)=fTF1->Eval(fQbb);
    }
}

// ---------------------------------------------------------
double BAT_GraphFit::LogLikelihood(const std::vector<double>& pars)
{
  // implement the gaussian log lilihood

  double logL=0;

  // gaussian is g(x)=1/(sigma*sqrt(2*pi))*TMath::Exp(-(x-mu)^2/2sigma^2)
  // hmm i think we should use the split gaussian
  // g(x) = N*TMath::Exp(-(x-mu)^2/2*sigma_l^2) x<mu //
  //      = N*TMath::Exp(-(x-mu)^2/(2*sigma_r^2) x>mu //
  // now we ave to find the N that is
  // N = sqrt(2pi)*sigma_l/2+sqrt(2pi)*sigma_r/2

  for (int n=0;n<fNpar;n++)
    {
      fTF1->SetParameter(n,pars[n]);
    }
  
  /// ************ GAUSSIAN *************************
  /// ***********************************************
  if (fMode=="G")
    {
      for (int i=0;i<fGraph->GetN();i++)
	{
	  // check if x is in the fit range
	  
	  double x =fGraph->GetX()[i];
	  if (x<fFitLow || x>fFitHigh)
	    {
	      continue;
	    }
	  double y=fGraph->GetY()[i];
	  
	  double mod_y = fTF1->Eval(x);
	  
	  double elow= fGraph->GetErrorYlow(i);
	  double ehigh = fGraph->GetErrorYhigh(i);
	  double exp;
	  double N = sqrt(2*TMath::Pi())*(elow+ehigh)/2.;
	  if (mod_y<y)
	    {
	      exp=-pow(mod_y-y,2)/(2*elow*elow);
	    }
	  else
	    {
	      exp=-pow(mod_y-y,2)/(2*ehigh*ehigh);
	    }
	  if (mod_y>=fGraphMinimum &&mod_y<fGraphMaximum)
	    {
	      logL+=log(N);
	      logL+=exp;
	    }
	  else
	    {
	      logL-=pow(10,9);
	    }
	}
    }

  //  ************ POISSON  *************************
  /// ***********************************************                                                                                                                                                     

  else if (fMode=="P")
    {

      double last=-1;
      for (int i=1;i<fTH1->GetNbinsX();i++)
	{
	  double E = fTH1->GetBinCenter(i);
	  if (E<fFitLow ||E>fFitHigh)
	    {
	      continue;
	    }
	  double width=fTH1->GetBinWidth(i);
	  double Elow = E-width/2.;
	  double Ehigh=E+width/2.;
	  double N = fTH1->GetBinContent(i);


	  // set pars of TF1
	  double mu;
	  if (last==-1)
	    {
	      last =fTF1->Eval(Elow);
	    }


	  // depends if the TF1 is the PDF or the CDF
	  if (fType=="C")
	    {
	      mu = fTF1->Eval(Ehigh)-last;
	      last = fTF1->Eval(Ehigh);
	    }
	  else
	    {
	      mu=fTF1->Eval(E);
	    }
	  if (mu>=0)
	    logL+=N*log(mu)-mu-BCMath::LogFact(N);
	  else
	    {
	      // penalty for negative model
	      logL-=pow(10,8);
	    }


	}
    }
  // **************COUNTING ANALYSIS*****************
  // ************************************************

  else if (fMode=="C")
    {
      
      // two cases - constant background or a TF1

      // calculate the values in the 3 bins, if the bkg is TF1 it must be integrate

      double B1 = fTF1->Eval(fCountingRange[0].second)-fTF1->Eval(fCountingRange[0].first);
      double B2 = fTF1->Eval(fCountingRange[1].second)-fTF1->Eval(fCountingRange[1].first);
      double B3 = fTF1->Eval(fCountingRange[2].second)-fTF1->Eval(fCountingRange[2].first);
      
      double S=pars[fNpar];
      
      double S1=S*fCountingProb[0];
      double S2=S*fCountingProb[1];
      double S3=S*fCountingProb[2];

      double N1=fCountingN[0];
      double N2=fCountingN[1];
      double N3=fCountingN[2];
      
      double mu1=B1+S1;
      double mu2=B2+S2;
      double mu3=B3+S3;
      
	  
      if (mu1>=0)
	logL+=N1*log(mu1)-mu1-BCMath::LogFact(N1);
      else
	{
	  // penalty for negative model                                                                                                                                                               
	  logL-=pow(10,8);
	}
      if (mu2>=0)
	logL+=N2*log(mu2)-mu2-BCMath::LogFact(N2);
      else
	{
	  // penalty for negative model                                                                                                                                                               
	  logL-=pow(10,8);
	}
      
      if (mu3>=0)
	logL+=N3*log(mu3)-mu3-BCMath::LogFact(N3);
      else
	{
	  // penalty for negative model                                                                                                                                                               
	  logL-=pow(10,8);
	}

	  

      
	
    }
  /// ************ BINOMIAL  ************************
  /// ***********************************************                                                                                                                                                     
  else if (fMode=="B")
    {
        for (int i=0;i<fGraphBinomTrial->GetN();i++)
        {
          // check if x is in the fit range                                                                                                                                                                

          double x =fGraphBinomTrial->GetX()[i];
          if (x<fFitLow || x>fFitHigh)
            {
              continue;
            }
	  double Ntrial = fGraphBinomTrial->GetY()[i];
	  double Nsucc = fGraphBinomSuccess->GetY()[i];
	  double p = fTF1->Eval(x);
	  double q=1-p;

	  if (p<1 &&p>0)
	    {
	      logL+=BCMath::LogBinomFactor(Ntrial,Nsucc);
	      logL+=Nsucc*log(p)+(Ntrial-Nsucc)*log(1-p);
	    }
	  else
	    {
	      // dont allow eps > 1 or <0 
	      logL+=pow(10,7);

	    }
      

	}
    }
  return logL;
  
  

  
}

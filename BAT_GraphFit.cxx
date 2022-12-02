// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "BAT_GraphFit.h"

 #include <BAT/BCMath.h>



void BAT_GraphFit::Plot(TString goption)
{
  if (goption=="" &&(fMode=="B"||fMode=="G"))
    {
      goption="APE";
    }
  else if (goption=="" &&fMode=="P")
    {
      goption="HIST";
    }
  if (fMode=="G")
    fGraph->Draw(goption);
  if (fMode=="P")
    fTH1->Draw(goption);
  if (fMode=="B")
    fGraph->Draw(goption);
  fTF1->Draw("Csame");
}

void BAT_GraphFit::SetGraph(TGraphAsymmErrors *&g)
{
  fGraph=g;
}
void BAT_GraphFit::SetHisto(TH1D*&h,TString type)
{
  fType=type;
  fTH1=h;
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
BAT_GraphFit::BAT_GraphFit(const std::string& name,TF1 *&f,bool verbose,TString mode)
    : BCModel(name)
{
  
  // get the information on the parameters from the TF1

  fTF1 = f;
  //fGraph =g;
  fNpar=fTF1->GetNpar();
  fVerbose=verbose;
  fMode=mode;
  // Default range from the TF1

  fTF1->GetRange(fFitLow,fFitHigh);
  if (fVerbose)
    {
      std::cout<<"RUNNING FIT WITH RANGE: "<<fFitLow<<" ' "<<fFitHigh<<std::endl;
      std::cout<<"***************************************************"<<std::endl;
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
	  std::cout<<"WARNING - Parameter "<<N<<" "<<name<<" has low = high = 0, set par limits !!!"<<std::endl;
	}
      AddParameter(name.Data(),low,high,name.Data(),"[]");
    }

  // how to set priors
  GetParameters().SetPriorConstantAll();
  GetParameters().SetNBins(300);
  
}




// ---------------------------------------------------------
BAT_GraphFit::~BAT_GraphFit()
{
    // destructor
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
	  logL+=log(N);
	  logL+=exp;
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
	      logL+=pow(10,8);
	    }


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
	      logL+Nsucc*log(p)+(Ntrial-Nsucc)*log(1-p);
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

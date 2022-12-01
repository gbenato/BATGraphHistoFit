// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "BAT_GraphFit.h"

// #include <BAT/BCMath.h>



void BAT_GraphFit::Plot(TString goption)
{
  if (goption=="")
    {
      goption="APE";
    }
  fGraph->Draw(goption);
  fTF1->Draw("Csame");
}

void BAT_GraphFit::SetGraph(TGraphAsymmErrors *&g)
{
  fGraph=g;
}
void BAT_GraphFit::SetHisto(TH1D*&h)
{
  fTH1=h;
}
void BAT_GraphFit::SetBinomGraph(TGraph *&gTrial,TGraph *&gSuccess)
{
  fGraphBinomTrial =gTrial;
  fGraphBinomSuccess = gSuccess;
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
	  
	  for (int n=0;n<fNpar;n++)
	    {
	      fTF1->SetParameter(n,pars[n]);
	    }
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
    }

  /// ************ BINOMIAL  ************************
  /// ***********************************************                                                                                                                                                     
  else if (fMode=="B")
    {

    }
  return logL;
  
  

  
}

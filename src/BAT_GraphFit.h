// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__BAT_GRAPHFIT__H
#define __BAT__BAT_GRAPHFIT__H

#include <BAT/BCModel.h>
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <string>
#include <vector>
#include "TH1D.h"
#include <tuple>
#include <array>
#include "TStyle.h"
//#include <pair>
// This is a BAT_GraphFit header file.
// Model source code is located in file BAT_GraphFit/BAT_GraphFit.cxx

// ---------------------------------------------------------
class BAT_GraphFit : public BCModel
{

public:

  TF1 *fTF1=nullptr;
  TGraphAsymmErrors *fGraph=nullptr;
  TGraph *fGraphBinomTrial;
  TGraph *fGraphBinomSuccess;
  TH1D *fTH1=nullptr;
  int fNpar;
  double fFitLow,fFitHigh;
  double fGraphMinimum,fGraphMaximum;
  TString fType;
  bool fVerbose;
  // Constructor
  void SetTF1(TF1*&f,bool addpars=1)
  {

    fTF1=f;
    fNpar=fTF1->GetNpar();
    fTF1->GetRange(fFitLow,fFitHigh);
    if (addpars)
      this->AddParameters();
    
  };

  void AddParameters();
  
  BAT_GraphFit(const std::string& name,TF1 *f,bool verbose,TString mode,double Sm=1000,double Qbb=2527);
  void CalculateObservables(const std::vector<double>&pars);

  void SetGraph(TGraphAsymmErrors *&g,double max=pow(10,5),double min=-pow(10,5));
  void SetHisto(TH1D*&h,TString type="P");
  void SetBinomGraph(TGraph *&gTrial,TGraph *&gSuccess);
  void PlotCI(TGraphAsymmErrors *g1,TGraphAsymmErrors *g2,TGraphAsymmErrors *g3);
  void SetCountingPars(  std::vector<std::pair<double,double>>range,std::vector<double> prob);

  void Plot(TString goption="APE");
  TString fMode;
  // Destructor
  ~BAT_GraphFit()
  {
  };

  double fQbb;
  // Overload LogLikelihood to implement model
  double LogLikelihood(const std::vector<double>& pars);
  
  std::vector<double> fCountingN;
  std::vector<double> fCountingProb;
  double fCountingSmax;
  std::vector<std::pair<double,double>>fCountingRange;
};
// ---------------------------------------------------------

#endif

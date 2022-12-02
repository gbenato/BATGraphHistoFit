
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "BAT_GraphFit.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
class BatGraphFitter
{
 public:

  // 5 constructors
  // 3 for Gaussian graph fit
  // 1 for Binomial fit
  BatGraphFitter(TGraphAsymmErrors *&g);
  BatGraphFitter(TH1D*&h);
  BatGraphFitter(TGraph*&g,double error=1);
  BatGraphFitter(TGraphErrors*&g);
  BatGraphFitter(TGraph *&gTrial,TGraph *&gSuccess);

  void Fit(TF1 *&f,TString option="",TString goption="",double low=0,double high=0);
  void SetPrecison(int prec){fPrec=prec;};

  int fPrec;
  ~BatGraphFitter(){};

  BAT_GraphFit *fModel;
  TGraphAsymmErrors *fGraph;
  TF1  *fTF1;
  TString fMode;
  TH1D *fHisto;
  TGraph *fGraphBinomTrial;
  TGraph *fGraphBinomSuccess;
  // outputs
  std::map<int,TH1D*> f1DPosteriors;
  std::map<int,std::map<int,TH2D*>>f2DPosteriors;
  
};

  

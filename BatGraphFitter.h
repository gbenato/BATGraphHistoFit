
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "BAT_GraphFit.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include <tuple>
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include <vector>
class BatGraphFitter
{
 public:

  // 5 constructors
  // 3 for Gaussian graph fit
  // 1 for Binomial fit
  BatGraphFitter(TGraphAsymmErrors *&g); // graph fit with errors
  BatGraphFitter(TH1D*&h);  // histo fit
  BatGraphFitter(TGraph*&g,double error=1); //graph fit without errors 
  BatGraphFitter(TGraphErrors*&g); // graph fit with symmetric errors
  BatGraphFitter(TGraph *g,TGraph *&gTrial,TGraph *&gSuccess); //binomial fit
  BatGraphFitter(TH1D *&h,std::vector<std::pair<double,double>> ranges,std::vector<double>probs,double Sm=100); // counting analysis

  void SetCountingPar(std::vector<std::pair<double,double>>range, std::vector<double>prob){fCountingRange=range;fCountingProb=prob;};
  void GetCredibleInterval(TGraphAsymmErrors * &grint1,TGraphAsymmErrors *&grint2,TGraphAsymmErrors *&grint3,int n,TTree *T_markov,double ylow,double yhigh);

  void Fit(TF1 *&f,TString option="",TString goption="",double low=0,double high=0);
  void SetPrecison(int prec){fPrec=prec;};
  void SetGraphMaxMin(double max,double min){fMax=max; fMin=min;};
    
  int fPrec;
  ~BatGraphFitter(){};
  double fMax,fMin;
  double fQbb=2527;
  void SetQbb(double Q){fQbb=Q;};
  BAT_GraphFit *fModel;
  TGraphAsymmErrors *fGraph;
  double fSmax;
  TF1  *fTF1;
  TString fMode;
  TH1D *fHisto;
  TGraph *fGraphBinomTrial;
  TGraph *fGraphBinomSuccess;
  double fCountingLow,fCountingHigh;
  std::vector<std::pair<double,double>> fCountingRange;
  std::vector<double> fCountingProb;
  // outputs
  std::map<int,TH1D*> f1DPosteriors;
  std::map<int,std::map<int,TH2D*>>f2DPosteriors;
  
};

  


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
#include <chrono>
class BatGraphFitter
{
 public:

  // 5 constructors
  // 3 for Gaussian graph fit
  // 1 for Binomial fit
  BatGraphFitter(TGraphAsymmErrors *&g,TF1 *&f); // graph fit with errors
  BatGraphFitter(TH1D*&h,TF1 *&f);  // histo fit
  BatGraphFitter(TGraph*&g,TF1*&f,double error=1); //graph fit without errors 
  BatGraphFitter(TGraphErrors*&g,TF1 *&f); // graph fit with symmetric errors
  BatGraphFitter(TGraph *g,TGraph *&gTrial,TGraph *&gSuccess,TF1 *&f); //binomial fit
  BatGraphFitter(TH1D *&h,TF1 *&f,std::vector<std::pair<double,double>> ranges,std::vector<double>probs,double Sm=100); // counting analysis
  BatGraphFitter(TH1D*&h,TF1*&f,TF1*&f_int,std::vector<double>energy);

  void SetTH1(TH1D*&h){fHisto=h;};
  void SetVector(std::vector<double>&vec){fEnergyVector=vec;};
  
  void SetCountingPar(std::vector<std::pair<double,double>>range, std::vector<double>prob){fCountingRange=range;fCountingProb=prob;};
  void GetCredibleInterval(TGraphAsymmErrors * &grint1,TGraphAsymmErrors *&grint2,TGraphAsymmErrors *&grint3,int n,TTree *&T_markov,double ylow,double yhigh);

  void Fit(TString option="",TString goption="",double low=0,double high=0,bool marg=1,bool quiet=0);
  void SetPrecison(int prec){fPrec=prec;};
  void SetGraphMaxMin(double max,double min){fMax=max; fMin=min;};
  void ResetTF1(TF1*&f);
  void ResetTF1Int(TF1*&fi);

  int fPrec;
  ~BatGraphFitter()
  {
    if (fGraph!=nullptr)delete fGraph;
    if (fGraph!=nullptr)delete fGraph;
  };
  TGraphAsymmErrors *fGrint;
  double fMax,fMin;
  double fQbb;
  void SetQbb(double Q){fQbb=Q;};
  BAT_GraphFit *fModel=nullptr;
  TGraphAsymmErrors *fGraph=nullptr;
  double fSmax;
  TF1  *fTF1=nullptr;
  TF1  *fTF1Int=nullptr;

  TString fMode;
  TH1D *fHisto=nullptr;
  TGraph *fGraphBinomTrial;
  TGraph *fGraphBinomSuccess;
  double fCountingLow,fCountingHigh;
  std::vector<std::pair<double,double>> fCountingRange;
  std::vector<double> fCountingProb;
  // outputs
  std::map<int,TH1D*> f1DPosteriors;
  std::map<int,std::map<int,TH2D*>>f2DPosteriors;
  std::vector<double>fEnergyVector;
  
};

  

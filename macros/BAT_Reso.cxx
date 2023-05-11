#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "src/BatGraphFitter.h"
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include <sstream>
#include "TGraphErrors.h"
// MACRO TO SCALE THE RESOLUTION USING BAT
// Toby Dixon: toby.dixon@universite-paris-saclay.fr 11/5/2023

void GetCI(TH1D *&h)
{

  double pe=h->GetBinCenter(h->GetMaximumBin());

  std::cout<<"point estimate = "<<pe<<std::endl;
  double q;
  double x=pe;


  q=h->Integral(0,h->FindBin(pe))/h->Integral();
  std::cout<<q<<std::endl;
  double qdown=q-0.683/2.;
  double qup=q+0.683/2.;
  if (qdown<0)
    {
      q=0.9;
      h->GetQuantiles(1,&x,&q);
      std::cout<<"no measurment can only set limit  < "<<x<<std::endl;
    }
  else if (qup>1)
    {
      q=0.1;
      h->GetQuantiles(1,&x,&q);
      std::cout<<"no measurment can only set limit  > "<<x<<std::endl;
    }
  else
    {

      double low,high;
      h->GetQuantiles(1,&low,&qdown);
      h->GetQuantiles(1,&high,&qup);

      double e_low = pe-low;
      double e_high=high-pe;

      std::cout<<" Measurment: "<<pe<<" +/- "<<e_high<<" / "<<e_low<<std::endl;
    }
}


struct LineShapeResults
{
  double EnergyNom;
  double EnergyFit;
  double ErrorEnergyFit;
  double Scaling;
  double ErrorScaling;
};
  
std::vector<LineShapeResults> ReadInputText2Map(TString path)
{
  // this method reads the input txt file and saves it in c++ objects

  std::vector<LineShapeResults> input;
  std::fstream file;
  file.open(path);

  while (file.is_open() &&!file.eof())
    {
      std::stringstream ss;
      std::string line;
      while (getline(file, line))
	{
	  double En,Ef,ErrEf,S,ErrS;
	  std::istringstream iss(line);
	  
	  iss >> En;
	  iss >> Ef;
	  iss >> ErrEf;
	  iss >> S;
	  iss >>ErrS;

	  LineShapeResults tmp;
	  tmp.EnergyNom = En;
	  tmp.EnergyFit=Ef;
	  tmp.ErrorEnergyFit=ErrEf;
	  tmp.Scaling=S;
	  tmp.ErrorScaling=ErrS;
	  input.push_back(tmp);
	}

    }

  std::cout<<"Reading the input files "<<path<<std::endl;
  std::cout<<"------------------------"<<std::endl;
  for (int i=0;i<input.size();i++)
    {
      std::cout<<"EnergyNom, EnergyFit, Scaling = "<<input[i].EnergyNom<<" "<<input[i].EnergyFit<<" +/- "<<input[i].ErrorEnergyFit<<
	" , "<<input[i].Scaling<<" +/- "<<input[i].ErrorScaling<<std::endl;
    }

  return input;
}
  

  
int main(int argc, char **argv)
{

  
  // please replace with a cfg file
  //std::map<TString,PeakInfo> peak_map;
  TApplication *app  =new TApplication("app",&argc,argv);
  gROOT->SetBatch(1);

  TString reso_type="sqrt(pol1)";

  TString function_bias="[0]+[1]*x+[2]*x*x";
  TString function_reso="[0]+[1]*x+[2]*x*x";
  int NparReso=3;
  if (reso_type=="pol1")
    {
      function_reso="[0]+[1]*x";
      NparReso=2;
    }
  else if (reso_type=="sqrt(pol1)")
    {
      function_reso="sqrt([0]*[0]+[1]*x)";
      NparReso=2;
    }
  else if (reso_type=="sqrt(pol2)")
    {
      function_reso="sqrt([0]*[0]+[1]*x+[2]*x*x)";
    }

  std::cout<<function_reso<<std::endl;
  // READ THE DATA
  // ----------------------------------------------------------
  
  // hardcoded should change !!!!
  int ds=atoi(argv[1]);
  double Qbb=2527;
  TString path="inputs/reso/fitresults";
  std::vector<LineShapeResults> input = ReadInputText2Map(Form("%s_ds%i.dat",path.Data(),ds));


  // create graphs

  TGraphErrors * gerror_bias = new TGraphErrors();
  TGraphErrors * gerror_reso = new TGraphErrors();

  for (int i=0;i<input.size();i++)
    {
      gerror_bias->SetPoint(i,input[i].EnergyNom,input[i].EnergyFit-input[i].EnergyNom);
      gerror_bias->SetPointError(i,0,input[i].ErrorEnergyFit);
      gerror_reso->SetPoint(i,input[i].EnergyNom,input[i].Scaling);
      gerror_reso->SetPointError(i,0,input[i].ErrorScaling);
      
    }
  gerror_bias->SetTitle(Form("Fit to energy bias for ds %i ; Energy Nominal [keV] ; Energy Fit - Energy Nominal [keV] ;",ds));
  gerror_reso->SetTitle(Form("Fit to energy resolution scaling for ds %i ; Energy Nominal [keV] ; Scaling [] ;",ds));
  
  path = Form("output/CUORE_reso/ds%i/",ds);

  TCanvas *ce = new TCanvas();


  
  // CREATE A FUNCTION TO FIT THE BIAS
  //------------------------------------------
  TF1 *fBias = new TF1("fBias","[0]+[1]*x+[2]*x*x",0,4000);

  gerror_bias->Fit(fBias);
  fBias->SetParLimits(0,fBias->GetParameters()[0]-5*fBias->GetParErrors()[0],fBias->GetParameters()[0]+5*fBias->GetParErrors()[0]);
  fBias->SetParLimits(1,fBias->GetParameters()[1]-5*fBias->GetParErrors()[1],fBias->GetParameters()[1]+5*fBias->GetParErrors()[1]);
  fBias->SetParLimits(2,fBias->GetParameters()[2]-5*fBias->GetParErrors()[2],fBias->GetParameters()[2]+5*fBias->GetParErrors()[2]);

  
  TF1 *fReso = new TF1("fReso",function_reso,0,4000);
  gerror_reso->Fit(fReso);
  
  fReso->SetParLimits(0,fReso->GetParameters()[0]-5*fReso->GetParErrors()[0],fReso->GetParameters()[0]+5*fReso->GetParErrors()[0]);
  fReso->SetParLimits(1,fReso->GetParameters()[1]-5*fReso->GetParErrors()[1],fReso->GetParameters()[1]+5*fReso->GetParErrors()[1]);
  fReso->SetParLimits(2,fReso->GetParameters()[2]-5*fReso->GetParErrors()[2],fReso->GetParameters()[2]+5*fReso->GetParErrors()[2]);

  // CREATE THE BAT fitter (bias)
  // -----------------------------------------------
  
  BatGraphFitter *fitter_bias = new BatGraphFitter(gerror_bias,fBias);
  fitter_bias->fModel->SetObsRange(-2,2);
  fitter_bias->SetPrecison(4);
  fitter_bias->SetQbb(Qbb);
  fitter_bias->SetGraphMaxMin(20,-20);  
  fitter_bias->Fit();

  // important - 3rd parameter is hardcoded
  TH1D* margdistro_bias = (TH1D*)fitter_bias->fModel->GetMarginalizedHistogram(3);
  GetCI(margdistro_bias);

  ce->Draw();
  ce->Print(Form("%s/bias.pdf(",path.Data()),"pdf");
  margdistro_bias->Draw();
  ce->Draw();
  ce->Print(Form("%s/bias.pdf)",path.Data()),"pdf");
  fitter_bias->fModel->PrintAllMarginalized(Form("%s/bias_fit.pdf",path.Data()));
  fitter_bias->fModel->WriteMarginalizedDistributions(Form("%s/output_bias.root",path.Data()), "RECREATE");


  // SAME FOR RESO
  // --------------------------------------------------------

  
  BatGraphFitter *fitter_reso = new BatGraphFitter(gerror_reso,fReso);
  fitter_bias->fModel->SetObsRange(0.8,1.2);

  fitter_reso->SetPrecison(4);
  fitter_reso->SetQbb(Qbb);
  fitter_reso->SetGraphMaxMin(2,0);
  fitter_reso->Fit();

  TH1D* margdistro_reso = (TH1D*)fitter_reso->fModel->GetMarginalizedHistogram(NparReso);
  GetCI(margdistro_reso);

  // save
  ce->Draw();
  ce->Print(Form("%s/reso.pdf(",path.Data()),"pdf");
  margdistro_reso->Draw();
  ce->Draw();
  ce->Print(Form("%s/reso.pdf)",path.Data()),"pdf");
  fitter_reso->fModel->PrintAllMarginalized(Form("%s/reso_fit.pdf",path.Data()));
  fitter_reso->fModel->WriteMarginalizedDistributions(Form("%s/output_reso.root",path.Data()), "RECREATE");

  

  
  return 1;
  
}
      

#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <utility>
#include "TSpectrum.h"
#include <fstream>
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
  std::cout<<qdown<<" "<<qup<<std::endl;
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
  std::ifstream file;
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

void Usage()
{
  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  if( std::getenv("USER") != NULL )
    std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
  else
    std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
  std::cout<<"./BAT_Reso"<<std::endl;
  std::cout<<"options "<<std::endl;
  std::cout<<"------------------------------------------------------------"<<std::endl;
  std::cout<<"-i (or --input-path)        [directory for input txt] (default: inputs/reso/fitresults)"<<std::endl;
  std::cout<<"-o (or --output-path)       [directory for output file (default: output/CUORE_reso/] "<<std::endl;
  std::cout<<"-d (or --dataset)           [dataset number : default 3602] "<<std::endl;
  std::cout<<"-r (or --reso-function)     [resolution function: default sqrt([0]*[0]+[1]*x)]"<<std::endl;
  std::cout<<"-b (or --bias-function)     [bias function: default [0]*[0]+[1]*x+[2]*x*x]"<<std::endl;
  std::cout<<"-l (or --label)             [label for outputs default "" ]"<<std::endl;
  std::cout<<"-q (or --q-value)           [q-value to evaluate the function at default: 2527 keV]"<<std::endl;
  std::cout<<"-p (or --precision)         [precision for efficiency fit on single peaks. Default: 3 (kHigh)]"<<std::endl;
  std::cout<<"-h (or --help)              [this help]"<<std::endl;

}
  
int main(int argc, char **argv)
{

  
  int ds=3602;
  TString bias_function="pol2";
  TString reso_function="[0]+[1]*x";//"sqrt([0]*[0]+[1]*x)";
  TString path="inputs/reso/fitresults";
  TString out_path = Form("output/CUORE_reso/");
  TString label="";
  double Qbb=2527;
  int precision=3;
  {
  static struct option long_options[] = {
					 { "input-path",        required_argument,  nullptr, 'i' },
					 { "output-path",       required_argument,  nullptr,'o'},
					 { "dataset",required_argument, nullptr,'d'},
					 { "bias-function",      required_argument,nullptr,'b'  },
					 { "reso-function",         required_argument,  nullptr,'r'},
					 { "label",    required_argument, nullptr,'l'},
					 { "q-value",           required_argument,  nullptr, 'q' },
					 { "precision", required_argument, nullptr, 'p' },
					 { "help",              no_argument,        nullptr,'h'},
					 {nullptr, 0, nullptr, 0}
  };
  const char* const short_options = "i:o:d:r:b:l:q:p:h";
  
  int c;
  
  while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 ) {
    switch (c) {
    case 'i': {
      path = optarg;
      break;
    }

    case 'o': {
      out_path = optarg;
      break;
    }
    case 'l':{
      label=optarg;
      break;
	}
      
    case 'r': {
      reso_function = optarg;
      break;
    }

    case 'd':{
      ds=atoi(optarg);
      break;
    }
    case 'b': {
      bias_function = optarg;
      break;
    }
    case 'q': {
      Qbb = atof(optarg);
      break;
    }
    case 'p':{
	precision = atoi(optarg);
	break;
    }
    case'h': {
      Usage();
      return 0;
    }
    default: {
      exit(1);
    }
      
    }

    
      
   } 
  }

    gROOT->SetBatch(1);
    

    // READ THE DATA
    // ----------------------------------------------------------
    
    std::vector<LineShapeResults> input = ReadInputText2Map(Form("%s_ds%i.dat",path.Data(),ds));
    
    out_path+=Form("/ds%i",ds);
    std::system("mkdir -p " + out_path);
    // create graphs
    
    TGraphErrors * gerror_bias = new TGraphErrors();
    TGraphErrors * gerror_reso = new TGraphErrors();
    
    for (int i=0;i<input.size();i++)
      {
	gerror_bias->AddPoint(input[i].EnergyNom,input[i].EnergyFit-input[i].EnergyNom);
	gerror_bias->SetPointError(gerror_bias->GetN()-1,0,input[i].ErrorEnergyFit);
	if( fabs( input[i].EnergyNom - 511.   ) > 1. &&  // Skip annihilation and 208Tl SEP peaks from reso scaling
	    fabs( input[i].EnergyNom - 2103.5 ) > 1. )   // because they're affected by Doppler broadening.
	    {
		gerror_reso->AddPoint(input[i].EnergyNom,input[i].Scaling);
		gerror_reso->SetPointError(gerror_reso->GetN()-1,0,input[i].ErrorScaling);
	    }
	
      }
    gerror_bias->SetTitle(Form("Fit to energy bias for ds %i ; Energy Nominal [keV] ; Energy Fit - Energy Nominal [keV] ;",ds));
    gerror_reso->SetTitle(Form("Fit to energy resolution scaling for ds %i ; Energy Nominal [keV] ; Scaling [] ;",ds));
  
    
    TCanvas *ce = new TCanvas();


  
    // CREATE A FUNCTION TO FIT THE BIAS
    //------------------------------------------
    TF1 *fBias = new TF1("fBias",bias_function,0,4000);
    
    gerror_bias->Fit(fBias);
    fBias->SetParLimits(0,fBias->GetParameters()[0]-12.*fBias->GetParErrors()[0],fBias->GetParameters()[0]+12.*fBias->GetParErrors()[0]);
    fBias->SetParLimits(1,fBias->GetParameters()[1]-12.*fBias->GetParErrors()[1],fBias->GetParameters()[1]+12.*fBias->GetParErrors()[1]);
    fBias->SetParLimits(2,fBias->GetParameters()[2]-12.*fBias->GetParErrors()[2],fBias->GetParameters()[2]+12.*fBias->GetParErrors()[2]);

  
    TF1 *fReso = new TF1("fReso",reso_function,0,4000);
    gerror_reso->Fit(fReso);
    
    for( int p=0; p<fReso->GetNpar(); p++ )
	{
	    double parMin = fReso->GetParameter(p) - 12. * fReso->GetParError(p);
	    if( parMin < 0. ) parMin = 0.;
	    double parMax = fReso->GetParameter(p) + 12. * fReso->GetParError(p);
	    fReso->SetParLimits( p, parMin, parMax );
	}
    //fReso->SetParLimits(0,fReso->GetParameters()[0]-12.*fReso->GetParErrors()[0],fReso->GetParameters()[0]+12.*fReso->GetParErrors()[0]);
    //fReso->SetParLimits(1,fReso->GetParameters()[1]-12.*fReso->GetParErrors()[1],fReso->GetParameters()[1]+12.*fReso->GetParErrors()[1]);
    //fReso->SetParLimits(2,fReso->GetParameters()[2]-12.*fReso->GetParErrors()[2],fReso->GetParameters()[2]+12.*fReso->GetParErrors()[2]);

    // CREATE THE BAT fitter (bias)
    // -----------------------------------------------
  
    BatGraphFitter *fitter_bias = new BatGraphFitter(gerror_bias,fBias);
    fitter_bias->fModel->SetObsRange(-2,2);
    fitter_bias->SetPrecison(precision);
    fitter_bias->SetQbb(Qbb);
    fitter_bias->SetGraphMaxMin(20,-20);  
    fitter_bias->fModel->WriteMarkovChain(Form("%s/output_bias_%s_mcmc.root",out_path.Data(),label.Data()), "RECREATE");

    fitter_bias->Fit();
    gerror_bias->SetTitle(Form("Fit to energy bias for ds %i ; Energy Nominal [keV] ; Energy Fit - Energy Nominal [keV] ;",ds));


    // important - 3rd parameter is hardcoded
    TH1D* margdistro_bias = (TH1D*)fitter_bias->fModel->GetMarginalizedHistogram(fBias->GetNpar());
    GetCI(margdistro_bias);
    
    ce->Draw();
    ce->Print(Form("%s/bias_%s.pdf(",out_path.Data(),label.Data()),"pdf");
    ce->SaveAs(Form("%s/bias_%s.C",out_path.Data(),label.Data()));

    margdistro_bias->Draw();
    ce->Draw();
    ce->SaveAs(Form("%s/bias_%s_post.C",out_path.Data(),label.Data()));

    ce->Print(Form("%s/bias_%s.pdf)",out_path.Data(),label.Data()),"pdf");
    fitter_bias->fModel->PrintAllMarginalized(Form("%s/bias_fit_%s.pdf",out_path.Data(),label.Data()));
    fitter_bias->fModel->WriteMarginalizedDistributions(Form("%s/output_bias_%s.root",out_path.Data(),label.Data()), "RECREATE");

    
    // SAME FOR RESO
    // --------------------------------------------------------

  
    BatGraphFitter *fitter_reso = new BatGraphFitter(gerror_reso,fReso);
    fitter_reso->fModel->SetObsRange(0.3,1.8);

    fitter_reso->SetPrecison(precision);
    fitter_reso->SetQbb(Qbb);
    fitter_reso->SetGraphMaxMin(2,0);
    fitter_reso->fModel->WriteMarkovChain(Form("%s/output_reso_%s_mcmc.root",out_path.Data(),label.Data()), "RECREATE");

    fitter_reso->Fit();
    gerror_reso->SetTitle(Form("Fit to energy resolution scaling for ds %i ; Energy Nominal [keV] ; Scaling [] ;",ds));

    TH1D* margdistro_reso = (TH1D*)fitter_reso->fModel->GetMarginalizedHistogram(fReso->GetNpar());
    GetCI(margdistro_reso);

    // save
    ce->Draw();
    ce->Print(Form("%s/reso_%s.pdf(",out_path.Data(),label.Data()),"pdf");
    ce->SaveAs(Form("%s/reso_%s.C",out_path.Data(),label.Data()));
    
    margdistro_reso->Draw();
    ce->Draw();
    ce->SaveAs(Form("%s/reso_%s_post.C",out_path.Data(),label.Data()));

    ce->Print(Form("%s/reso_%s.pdf)",out_path.Data(),label.Data()),"pdf");
    fitter_reso->fModel->PrintAllMarginalized(Form("%s/reso_fit_%s.pdf",out_path.Data(),label.Data()));
    fitter_reso->fModel->WriteMarginalizedDistributions(Form("%s/output_reso_%s.root",out_path.Data(),label.Data()), "RECREATE");
    fitter_reso->fModel->WriteMarkovChain(Form("%s/output_reso_%s_mcmc.root",out_path.Data(),label.Data()), "RECREATE");

  

  
  return 1;
  
}
      

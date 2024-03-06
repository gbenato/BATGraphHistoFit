#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <utility>

#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include "TSpectrum.h"
#include "TCanvas.h"

#include <BAT/BCLog.h>

#include "src/BatGraphFitter.h"

// MACRO TO CALCULATE THE EFFICIENCY USING THE PEAKS IMPLEMENTED IN BAT
// Toby Dixon: toby.dixon@universite-paris-saclay.fr 11/5/2023
std::vector<TH1D*>GetHistogramsFrom2D( TH2D* &h,
				       TString title,
				       std::vector<double>energies={1173,1333,1461,2615},
				       std::vector<double> dE_center={10,10,10,10},
				       std::vector<double>dE_side={30,30,30,30},
				       std::vector<double>* single_energies=0,
				       int precision=3)
{
    std::vector<TH1D*>out;

    TCanvas *ci= new TCanvas();
    ci->Print(Form("%s.pdf(",title.Data()),"pdf");
    int counter=0;
    if (dE_side.size()!=energies.size()||dE_center.size()!=energies.size())
	{
	    std::cout<<"ERROR: the size of sideband not same as center energy, check the cfg file!!"<<std::endl;
	    throw;
	}

    double minX   = h->GetXaxis()->GetXmin();
    double maxX   = h->GetXaxis()->GetXmax();
    double dX     = h->GetXaxis()->GetBinWidth(1);
    int    nBinsX = h->GetNbinsX();
    double minY   = h->GetYaxis()->GetXmin();
    double maxY   = h->GetYaxis()->GetXmax();
    double dY     = h->GetYaxis()->GetBinWidth(1);
    int    nBinsY = h->GetNbinsY();
    std::cout << "Bin widths: " << dX << "\t" << dY << std::endl;    
    std::cout << "nBins: " << nBinsX << "\t" << nBinsY << std::endl;
    
    //TApplication* app = new TApplication("app",NULL,0);
    //TCanvas* can = new TCanvas();
    //can->cd();
    //can->Divide(2,1);
    for (auto & E : energies)
	{
	    int nSlides = 10;
	    double delta = E / nSlides;
	    for( int i=0; i<nSlides; i++ )
		{
		    double thisE = delta * ( 0.5 + i );

		    std::vector<std::pair<double,double>> range;
		    range.push_back(std::make_pair(E-dE_side[counter],E-dE_center[counter]));
		    range.push_back(std::make_pair(E-dE_center[counter],E+dE_center[counter]));
		    range.push_back(std::make_pair(E+dE_center[counter],E+dE_side[counter]));
		    std::vector<double>probs{0,1,0};
		    std::string thisHistoName = std::to_string(E)+"_"+std::to_string(i);
		    int nBins = range[2].second - range[0].first;
		    TH1D* thisH = new TH1D( thisHistoName.c_str(), thisHistoName.c_str(), nBins, range[0].first, range[2].second );

		    for( int xx=1; xx<=nBinsX; xx++ )
			for( int yy=1; yy<=nBinsY; yy++ )
			    {
				// Compute bin center
				double x = dX * ( 0.5 + xx );
				double y = dY * ( 0.5 + yy );

				// Cut on sum energy
				if( x+y <= range[0].first ||
				    x+y >  range[2].second )
				    continue;

				// Cut on energy of first event in the M2 duplet
				if( x < delta*i || x >= delta*(i+1) )
				    continue;

				// Fill histogram of sum energy
				int bin1d = thisH->FindBin(x+y);
				int bin2d = h->FindBin(x,y);
				thisH->AddBinContent(bin1d,h->GetBinContent(bin2d));				    
			    }
		    //thisH->Draw();
		    //app->Run(kTRUE);

		    // Initial parameters for fitting function
		    double max = thisH->GetMaximum();
		    int bin1   = thisH->FindBin(E-dE_side[counter]);
		    int bin2   = thisH->FindBin(E-dE_center[counter]);
		    int bin4   = thisH->FindBin(E+dE_side[counter]);
		    int bin3   = thisH->FindBin(E+dE_center[counter]);

		    double binwidth = thisH->GetBinWidth(bin1);
      
		    double Bl   = thisH->Integral(bin1,bin2)/(bin2-bin1);
		    double Br   = thisH->Integral(bin3,bin4)/(bin4-bin3);
		    double Bmax = 3.*0.5 * (Br+Bl);
		    if( Bmax < 10. ) Bmax = 10.;
		    double Smax = max-(Br+Bl)/2.;
		    if( Smax < 10. ) Smax = 10.;
		    double Sest = 50+fabs(thisH->Integral(bin2,bin3)-dE_center[counter]*(Br+Bl));
      
		    double slope_est = (Br-Bl)/(dE_center[counter]+dE_side[counter]);
		    double slope_max = (Br+2)/(dE_center[counter]+dE_side[counter]);

		    // Bkg fitting function (linear)
		    TF1 *fInt=new TF1("fInt",Form("[0]*x+[1]*pow(x-%f,2)/2.",E),E-dE_side[counter]-5,E+dE_side[counter]+5);
		    fInt->SetParLimits(0,0,Bmax);
		    fInt->SetParameter(0,(Br+Bl)/2.);
		    fInt->SetParameter(1,slope_est);
		    fInt->SetParLimits(1,-5*fabs(slope_max),+5*fabs(slope_max));

		    // Create fitter
		    std::cout<<"Sest = "<<Sest<<std::endl;
		    BatGraphFitter *fitter_count_pass = new BatGraphFitter(thisH,fInt,range,probs,5*Sest);
		    fitter_count_pass->SetPrecison(precision);

		    // Fit
		    fitter_count_pass->Fit();
		    ci->Print(Form("%s.pdf",title.Data()),"pdf");

		    // GET THE POSTERIOR
		    // -------------------------------------------------------
		    TH1D* margdistro = (TH1D*)fitter_count_pass->fModel->GetMarginalizedHistogram("S" );
		    margdistro->SetTitle(Form("Peak at %f keV ; ; ; ",E));
		    margdistro->Draw();
		    ci->Draw();
		    ci->Print(Form("%s.pdf",title.Data()),"pdf");
		    out.push_back(margdistro);
		    single_energies->push_back(thisE);
		}
	    counter ++;
	}
      ci->Print(Form("%s.pdf)",title.Data()),"pdf");
    
    return out;
}

std::vector<TH1D*>GetHistograms( TH1D* &h,
				 TString title,
				 std::vector<double>energies={1173,1333,1461,2615},
				 std::vector<double> dE_center={10,10,10,10},
				 std::vector<double>dE_side={30,30,30,30},
				 int precision=3)
{
  // ** THIS IS THE METHOD THAT EXTRACTS THE posterior of the number of counts for each peak
  // You should supply a vector of energies
  // 1) A TH1D* of the data
  // 2) Title is the name of a path to put summary plots
  // 3) energies of the peaks
  // 4) THe range for center is +/- dE_center, for left its -dE_side --> -dE_center and right is dE_center-- dE_side
  

  
  std::vector<TH1D*>out;
  TCanvas *ci= new TCanvas();
  ci->Print(Form("%s.pdf(",title.Data()),"pdf");
  int counter=0;
  if (dE_side.size()!=energies.size()||dE_center.size()!=energies.size())
    {
      std::cout<<"ERROR: the size of sideband not same as center energy, check the cfg file!!"<<std::endl;
      throw;
    }
  for (auto & E : energies)
    {

      /* CREATE THE CANVAS
	 -----------------------------------------------------
      */

	std::cout << "----------------------" << std::endl;
      std::cout<<"Running the fit of "<<E<<std::endl;
      std::cout<<"side = "<<dE_side[counter]<<" center = "<<dE_center[counter]<<std::endl;
      ci= new TCanvas(Form("c_%i",(int)E),Form("c_%i",(int)E));
      ci->cd();
      
      // PARAMETERS FOR THE FIT
      // ------------------------------------------------------
      std::vector<std::pair<double,double>> range;
      range.push_back(std::make_pair(E-dE_side[counter],E-dE_center[counter]));
      range.push_back(std::make_pair(E-dE_center[counter],E+dE_center[counter]));
      range.push_back(std::make_pair(E+dE_center[counter],E+dE_side[counter]));
      std::vector<double>probs{0,1,0};

      h->GetXaxis()->SetRangeUser(E-dE_side[counter]-5,E+dE_side[counter]+5);
      

      // INITIAL PARS
      // ------------------------------------------------------
      double max = h->GetMaximum();
      int bin1   = h->FindBin(E-dE_side[counter]);
      int bin2   = h->FindBin(E-dE_center[counter]);
      int bin4   = h->FindBin(E+dE_side[counter]);
      int bin3   = h->FindBin(E+dE_center[counter]);

      double binwidth = h->GetBinWidth(bin1);
      
      double Bl   = h->Integral(bin1,bin2)/(bin2-bin1);
      double Br   = h->Integral(bin3,bin4)/(bin4-bin3);
      double Smax = max-(Br+Bl)/2.;
      double Sest = 50+fabs(h->Integral(bin2,bin3)-dE_center[counter]*(Br+Bl));
      
      double slope_est = (Br-Bl)/(dE_center[counter]+dE_side[counter]);
      double slope_max = (Br+2)/(dE_center[counter]+dE_side[counter]);
      
      // MAKE BKG FIT FUNCTION
      //-------------------------------------------------------
      TF1 *fInt=new TF1("fInt",Form("[0]*x+[1]*pow(x-%f,2)/2.",E),E-dE_side[counter]-5,E+dE_side[counter]+5);
      fInt->SetParLimits(0,0,3*(Br+Bl)/2.);
      fInt->SetParameter(0,(Br+Bl)/2.);
      fInt->SetParameter(1,slope_est);
      fInt->SetParLimits(1,-5*fabs(slope_max),+5*fabs(slope_max));
     
      // MAKE THE FITTER
      // ------------------------------------------------------
      std::cout<<"Sest = "<<Sest<<std::endl;
      BatGraphFitter *fitter_count_pass = new BatGraphFitter(h,fInt,range,probs,5*Sest);
      fitter_count_pass->SetPrecison(precision);
      

      // RUN THE FIT AND DRAW
      //-------------------------------------------------------

      fitter_count_pass->Fit();
     
      ci->Print(Form("%s.pdf",title.Data()),"pdf");
      

      // GET THE POSTERIOR
      // -------------------------------------------------------
      TH1D* margdistro = (TH1D*)fitter_count_pass->fModel->GetMarginalizedHistogram("S" );
      margdistro->SetTitle(Form("Peak at %f keV ; ; ; ",E));
      margdistro->Draw();
      
      ci->Draw();
      ci->Print(Form("%s.pdf",title.Data()),"pdf");
      out.push_back(margdistro);
      counter++;

    }
  ci->Print(Form("%s.pdf)",title.Data()),"pdf");

  // TOby: Should delete all objects - I was lazy to do it since the method is only called 2 times
  return out;

  
}

// method to combine the two probability distributions of Npass and Nfail to one on efficiency
TH1D * toy_combine(TH1D*hpass,TH1D*hfail,double &mu,double &sigma,double E,int ds,TString mode="P")
{
    std::string name = "eff_" + std::to_string((int)E);
    TH1D * hout = new TH1D(name.c_str(),name.c_str(),6000,0,1.2);
    int Ntoys = hpass->GetEntries() / 10;
  for (int i=0;i<Ntoys;i++)
    {
      double Np=hpass->GetRandom();
      double Nf=hfail->GetRandom();

      hout->Fill(Np/(Np+Nf));

    }
  hout->GetXaxis()->SetRangeUser(0.5,1.1);
  hout->SetTitle(Form("DS %i - Efficiency posterior for energy %i keV ; Efficiency ; Probability ; ",ds,(int)E));

  if( mode == "P" )
      {
	  TF1 *fGauss= new TF1("fGauss","gaus",0,1);
	  hout->Fit(fGauss);

	  mu=fGauss->GetParameters()[1];
	  sigma=fGauss->GetParameters()[2];
      }
  else if( mode == "P2" )
      {
	  TF1 *fGauss= new TF1("fGauss","gaus",0,1);
	  hout->Fit(fGauss,"N");

	  mu=fGauss->GetParameters()[1];
	  sigma=fGauss->GetParameters()[2];
      }
  hout->Draw();
  
  return hout;
}


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
      double low,high;
      qdown=0;
      h->GetQuantiles(1,&low,&qdown);
      h->GetQuantiles(1,&high,&qup);
      double e_low = pe-low;
      double e_high=high-pe;

      std::cout<<" Measurment: "<<pe<<" +/- "<<e_high<<" / "<<e_low<<std::endl;


    }
  else if (qup>1)
    {
      q=0.1;
      h->GetQuantiles(1,&x,&q);
      std::cout<<"no measurment can only set limit  > "<<x<<std::endl;
      double low,high;
      qup=1;
      h->GetQuantiles(1,&low,&qdown);
      h->GetQuantiles(1,&high,&qup);
      double e_low = pe-low;
      double e_high=high-pe;

      std::cout<<" Measurment: "<<pe<<" +/- "<<e_high<<" / "<<e_low<<std::endl;

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

void ReadCfg(TString path,std::vector<double>&E,std::vector<double>&center,std::vector<double>&side,std::vector<bool>&on)
{
  // method to read the cfg options


  std::fstream file;
  file.open(path);
  std::cout<<"Reading file "<<path<<std::endl;
  if (!file.is_open())
    {
      std::cout<<"ERROR: cfg file doesnt exist"<<std::endl;
      throw;
    }
  while (file.is_open() &&!file.eof())
    {
      std::stringstream ss;
      std::string line;
      getline(file,line);
      while (getline(file, line))
        {

	  double E_tmp,Center_tmp,Side_tmp;
	  bool On_tmp;
          std::istringstream iss(line);

          iss >> E_tmp;
          iss >> Center_tmp;
          iss >> Side_tmp;
          iss >> On_tmp;
	  
	  if(On_tmp)
	      {
		  E.push_back(E_tmp);
		  center.push_back(Center_tmp);
		  side.push_back(Side_tmp);
		  on.push_back(On_tmp);
		  std::cout<<"Add peak at energy "<<E_tmp<<" and range "<<"center "<<Center_tmp<<" side "<<Side_tmp<<" which is on ? "<<On_tmp<<std::endl;
	      }
	}
    }
}


void Usage()
{
  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  if( std::getenv("USER") != NULL )
    std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
  else
    std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
  std::cout<<"./BAT_Efficiency"<<std::endl;
  std::cout<<"options "<<std::endl;
  std::cout<<"------------------------------------------------------------"<<std::endl;
  std::cout<<"-i (or --input-path)        [directory for input txt] (default: inputs/eff/m1_eff_histo_withtimecut_PCACut_directsum_)"<<std::endl;
  std::cout<<"-o (or --output-path)       [directory for output file (default: output/CUORE_eff/] "<<std::endl;
  std::cout<<"-c (or --cfg-path)          [directory for output file (default: cfs/EffPeaks.cfg] "<<std::endl;
  std::cout<<"-d (or --dataset)           [dataset number : default 3602] "<<std::endl;
  std::cout<<"-m (or --mode)              [mode either P for PCA, or M for multiplicity]"<<std::endl;
  std::cout<<"-f (or --function)          [eff function: default [0]+[1]*x/4000]"<<std::endl;
  std::cout<<"-l (or --label)             [label for outputs default "" ]"<<std::endl;
  std::cout<<"-q (or --q-value)           [q-value to evaluate the function at default: 2527 keV]"<<std::endl;
  std::cout<<"-p (or --precision)         [precision for efficiency fit on single peaks. Default: 3 (kHigh)]"<<std::endl;
  std::cout<<"-h (or --help)              [this help]"<<std::endl;

}


	  
  
int main(int argc, char **argv)
{

  int ds=3602;
  TString fit_function="[0]+[1]*x/4000.";
  TString path="inputs/eff/m1_eff_histo_withtimecut_PCACut_directsum_";
  TString out_path = Form("output/CUORE_eff");
  TString label="";
  double Qbb=2527;
  TString cfg_path="cfg/EffPeaks.cfg";
  TString mode="P";
  int precision=3;  
  {
    static struct option long_options[] = {
                                         { "input-path",        required_argument,  nullptr, 'i' },
                                         { "output-path",       required_argument,  nullptr,'o'},
					 { "cfg-path",       required_argument,  nullptr,'c'},
					 { "dataset",required_argument, nullptr,'d'},
                                         { "function",      required_argument,nullptr,'f'  },
                                         { "label",    required_argument, nullptr,'l'},
					 { "mode", required_argument,nullptr,'m'},
                                         { "q-value",           required_argument,  nullptr, 'q' },
					 { "precision", required_argument, nullptr, 'p' },
                                         { "help",              no_argument,        nullptr,'h'},
                                         {nullptr, 0, nullptr, 0}
  };

    const char* const short_options = "i:o:d:f:l:c:q:m:p:h";

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
      case 'd':{
	ds=atoi(optarg);
	break;
      }
      case 'c': {
	cfg_path = optarg;
        break;
      }
	
      case 'l':{
	label=optarg;
	break;
      }
	
      case 'f': {
	fit_function = optarg;
	break;
      }

	
      case 'q': {
	Qbb = atof(optarg);
	break;
      }
      case 'm':{
	mode=optarg;
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
	std::cout<<"def "<<optarg<<std::endl;
	exit(1);
      }
	
      }
      


    }
  }






  // READ THE DATA
  // ----------------------------------------------------------
  TFile *f = new TFile(Form("%sds%i.root",path.Data(),ds));

 
  //TH1D *h = (TH1D*)f->Get(Form("h_ds%i",ds));
  //TH1D *hb = (TH1D*)f->Get(Form("hb_ds%i",ds));
  //TH1D *hm = (TH1D*)f->Get(Form("hm_ds%i",ds));
  TH1D *hm          = (TH1D*)f->Get(Form("hmbds%i",ds));//TH1D *hmb = (TH1D*)f->Get(Form("hmbds%i",ds));
  TH1D *hpca        = (TH1D*)f->Get(Form("hpca_ds%i",ds));
  TH1D *hm_fail     = (TH1D*)f->Get(Form("hmbfds%i",ds));
  TH1D *hpca_fail   = (TH1D*)f->Get(Form("hpcafds%i",ds));
  TH2D *hpca2d      = (TH2D*)f->Get(Form("hm2pca_2d_ds%i",ds));
  TH2D *hpca2d_fail = (TH2D*)f->Get(Form("hm2pcaf_2d_ds%i",ds));

  hpca->SetLineColor(1);
  
  hpca->Rebin(1/hpca->GetBinWidth(2));
  hpca_fail->Rebin(1/hpca_fail->GetBinWidth(2));
  hm->Rebin(1/hm->GetBinWidth(2));
  hm_fail->Rebin(1/hm_fail->GetBinWidth(2));

  // --------------------------------------------------
  // SET THE ENERGIES - should change to cfg file

  std::vector<double>energies;
  std::vector<double>dE_center;
  std::vector<double>dE_side;
  std::vector<bool>is_on;
  
  ReadCfg(cfg_path,energies,dE_center,dE_side,is_on);

  
  // set the output path
  //---------------------------------------------------
  out_path+=TString("_")+mode;
  TString out_path_ds = Form("%s/ds%i%s",out_path.Data(),ds,label.Data());
  std::system("mkdir -p " + out_path_ds);
  
  // RUNT THE COUNTING ANALYSES
  // -------------------------------------------------
  std::vector<TH1D*>hpass;
  std::vector<TH1D*>hfail; 
  std::vector<double>* single_energies = new std::vector<double>();
  if (mode=="P")
    {
	hpass = GetHistograms(hpca,Form("%s/out_pca_pass",out_path_ds.Data()),energies,dE_center,dE_side,precision);
	hfail = GetHistograms(hpca_fail,Form("%s/out_pca_fail",out_path_ds.Data()),energies,dE_center,dE_side,precision);
    }
  else if ( mode=="P2")
      {
	  hpass = GetHistogramsFrom2D(hpca2d,Form("%s/out_pca_pass",out_path_ds.Data()),energies,dE_center,dE_side,single_energies,precision);
	  hfail = GetHistogramsFrom2D(hpca2d_fail,Form("%s/out_pca_fail",out_path_ds.Data()),energies,dE_center,dE_side,single_energies,precision);
      }
  else if (mode=="M")
    {
      // check this
      hpass = GetHistograms(hm,Form("%s/out_multi_pass",out_path_ds.Data()),energies,dE_center,dE_side,precision);
      hfail = GetHistograms(hm_fail,Form("%s/out_multi_fail",out_path_ds.Data()),energies,dE_center,dE_side,precision);

    }


  // FOR EACH PEAK GET A POSTERIOR ON EFFICINECY (TOY MC) and save them to graph
  // -------------------------------------------------
  TGraphErrors *gerror = new TGraphErrors(); 

  TCanvas *ce = new TCanvas();
  ce->Draw();
  ce->cd();
  ce->Print(Form("%s/eff.pdf(",out_path_ds.Data()),"pdf");

  TFile *peakeff_file = new TFile( Form("%s/eff.root",out_path_ds.Data()), "RECREATE");
  peakeff_file->cd();

  std::vector<TH1D*>* all_hout = new std::vector<TH1D*>();
  int counter=0;
  if ( mode=="P2")
      {
	  for (int i=0;i<hpass.size();i++)
	      {
		  double mu,sigma;
		  TH1D* hout = toy_combine(hpass[i],hfail[i],mu,sigma,single_energies->at(i),ds,mode);
		  
		  gerror->SetPoint(counter,single_energies->at(i),mu);
		  gerror->SetPointError(counter,0,sigma);
		  counter ++;
		  
		  ce->cd();
		  hout->Draw();
		  ce->Draw();
		  ce->Print(Form("%s/eff.pdf",out_path_ds.Data()),"pdf");
		  hout->Write();
		  all_hout->push_back( (TH1D*)hout->Clone() );
		  all_hout->back()->SetDirectory(gROOT);// We need this to keep all_hout in memory after we close this fucking ROOT file!

	      }
      }
  else
      {

	  for (int i=0;i<hpass.size();i++)
	      {
		  double mu,sigma;
		  TH1D * hout=toy_combine(hpass[i],hfail[i],mu,sigma,energies[i],ds,mode);
		  
		  if (is_on[i]==true)
		      {
			  gerror->SetPoint(counter,energies[i],mu);
			  gerror->SetPointError(counter,0,sigma);
			  counter++;
		      }
		  ce->cd();
		  hout->Draw();
		  ce->Draw();
		  ce->Print(Form("%s/eff.pdf",out_path_ds.Data()),"pdf");
		  hout->Write();
	      }
      }
    
  if( mode == "M" )
      {
	  ce->SaveAs(Form("%s/eff.pdf)",out_path_ds.Data()),"pdf");
	  return 0;
      }
  else if( mode == "P" )
      {
	  gerror->SetTitle(Form("PCA efficiency for ds %i ; Energy [keV] ; Efficiency ",ds));
	  gerror->SaveAs(Form("%s/graph.root",out_path_ds.Data()));
      }


  // CREATE A FUNCTION TO FIT THE EFFICIENCY
  //------------------------------------------
  TF1 *fEff = new TF1("fEff",fit_function,0,4000);
  fEff->SetParameter(0,0.8);
  fEff->SetParLimits(0,0,1.2);
  fEff->SetParameter(1,0);
  fEff->SetParLimits(1,-1,1);


  // CREATE THE BAT fitter
  // -----------------------------------------------
  
  BatGraphFitter *fitter;
  if ( mode=="P2")
      fitter = new BatGraphFitter( all_hout,fEff,*single_energies);
  else
      fitter = new BatGraphFitter(gerror,fEff);
  fitter->SetPrecison(precision);
  fitter->SetQbb(Qbb);
  fitter->SetGraphMaxMin(1,0);


  // RUN THE GRAPH FIT - adding a confidence interval calculation
  // -----------------------------------------------------------
  fitter->fModel->WriteMarkovChain(Form("%s/output_mcmc.root",out_path_ds.Data()), "RECREATE");
  fitter->Fit();


  
  TH1D* margdistro = (TH1D*)fitter->fModel->GetMarginalizedHistogram(fEff->GetNpar());

  GetCI(margdistro);
  
  // SAVE THE OUTPUT
  // ------------------------------------------------------------

  if( mode == "P" )
      {
	  ce->Draw();
	  ce->Print(Form("%s/eff.pdf",out_path_ds.Data()),"pdf");
      }
  else if( mode == "P2" )
      {
	  gerror->Draw("AP");
	  fitter->GetFittingFunction()->Draw("same");
	  ce->Draw();
	  ce->Print(Form("%s/eff.pdf",out_path_ds.Data()),"pdf");
      }
  margdistro->Draw();
  ce->Draw();
  ce->Print(Form("%s/eff.pdf)",out_path_ds.Data()),"pdf");

  fitter->fModel->PrintAllMarginalized(Form("%s/eff_fit.pdf",out_path_ds.Data()));
  fitter->fModel->WriteMarginalizedDistributions(Form("%s/output_eff.root",out_path_ds.Data()), "RECREATE");


  return 1;
  
}
      

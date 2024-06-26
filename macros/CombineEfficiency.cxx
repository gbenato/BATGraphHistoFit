#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <utility>
#include <fstream>
#include <iomanip>

#include "TFile.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include "TString.h"
#include "TCanvas.h"


// MACRO TO combine the efficiencies and put into a format for BAT-hl
// Toby Dixon: toby.dixon@universite-paris-saclay.fr 23/5/2023    
// Possibility to convolve the systematics (??)


TH1D* Combine(TString path_P,TString path_M)
{



  // Get the histograms
  //------------------------------------------------------------------

  TFile *fP = new TFile(path_P);
  TH1D *hP = (TH1D*)fP->Get("h1_fit_observable_fQ");
  

  TFile *fM = new TFile(path_M);
  TH1D *hM = (TH1D*)fM->Get("h1_fit_observable_fQ");


  TH1D *heff= new TH1D("heff","heff",1000,0,1);
  

  for (int i=0;i<1000000;i++)
    {


      double effP=hP->GetRandom();
      double effM=hM->GetRandom();

      heff->Fill(effP*effM);

    }
  heff->Scale(1/(heff->GetBinWidth(1)*heff->Integral()));
  
  return heff;
}

void Combine( std::vector<TH1D*>& inputhisto,
	      TH1D* outputhisto )
{
    int N = std::numeric_limits<int>::max();
    for( auto& h: inputhisto )
	if( h->GetEntries() < N )
	    N = h->GetEntries();
    //N /= 5;
    std::cout << "N: " << N << std::endl;
    std::cout << "nH: " << inputhisto.size() << std::endl;    
    for( int i=0; i<N; i++ )
	{
	    size_t nH = inputhisto.size();
	    std::vector<double> inputvalue(nH);
	    double outputvalue = 1.;
	    for( size_t v=0; v<nH; v++ )
		inputvalue[v] = inputhisto[v]->GetRandom();
	    for( auto& v: inputvalue )
		outputvalue *= v;
	    
	    outputhisto->Fill(outputvalue);
	}

    outputhisto->Scale( 1. / outputhisto->Integral(1,outputhisto->GetNbinsX()) );
    
    return;
}

void ComputeShortestInterval( int ds, TH1D* h)
{
    int    modeBin  = h->GetMaximumBin();
    double mode     = h->GetBinCenter(modeBin);
    int    binLeft  = modeBin;
    int    binRight = modeBin;
    double binWidth = h->GetBinWidth(1);
    
    double integral = h->GetBinContent(modeBin);
    double quantile = 0.6826895;// Corresponds to +-1 sigma for a Gaussian distribution
    while( integral < quantile )
	{
	    double intLeft  = h->GetBinContent( binLeft - 1 );
	    double intRight = h->GetBinContent( binRight + 1 );
	    if( intLeft > intRight )
		{
		    integral += intLeft;
		    binLeft --;
		}
	    else if( intLeft < intRight )
		{
		    integral += intRight;
		    binRight ++;
		}
	    else// It will probably never happen that intLeft == intRight
		{
		    integral += intLeft + intRight;
		    binLeft --;
		    binRight ++;
		}
	}
    double errShortLeft  = binWidth * ( modeBin - binLeft );
    double errShortRight = binWidth * ( binRight - modeBin );

    std::cout << "DS: " << ds << std::endl;
    std::cout << "Mode:  " << std::setprecision(10) << mode << std::endl;
    std::cout << "Low error for shortest interval method:  " << std::setprecision(10) << errShortLeft << std::endl;
    std::cout << "High error for shortest interval method: " << std::setprecision(10) << errShortRight << std::endl;
    std::cout << "Coverage for equal-areas method: "
	      << h->Integral( binLeft, binRight )
	      << std::endl << std::endl;
    return;
}

void Usage()
{
  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  if( std::getenv("USER") != NULL )
    std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
  else
    std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
  std::cout<<"./CombineEfficiency"<<std::endl;
  std::cout<<"options "<<std::endl;
  std::cout<<"------------------------------------------------------------"<<std::endl;
  std::cout<<"-i (or --input-path)        [directory for input] (default: none)"<<std::endl;
  std::cout<<"-n (or --input-file)        [input file name] (default: none)"<<std::endl;
  std::cout<<"-N (or --input-histo)       [input histo name] (default: none)"<<std::endl;
  std::cout<<"-a (or --input-ascii-file)  [input ascii file] (default: none)"<<std::endl;
  std::cout<<"-o (or --output-path)       [directory for output file (default: output/CUORE_eff/] "<<std::endl;
  std::cout<<"-d (or --dataset)           [dataset number : default none] "<<std::endl;
  std::cout<<"-f (or --first-dataset)     [dataset number : default 0] "<<std::endl;
  std::cout<<"-l (or --last-dataset)      [dataset number : default 0] "<<std::endl;
  std::cout<<"-h (or --help)              [this help]"<<std::endl;

}

int main(int argc, char **argv)
{

    std::vector<std::string> inputdirs;
    std::vector<std::string> inputname;
    std::vector<std::string> inputhistoname;
    std::vector<std::string> inputasciiname;
    std::string outpath;
    std::vector<int> ds;
    int firstds = 0;
    int lastds  = -1;
    static struct option long_options[] = {
	{ "input-path",        required_argument, nullptr,'i'},
	{ "input-file",        required_argument, nullptr,'n'},
	{ "input-histo",       required_argument, nullptr,'N'},
	{ "input-ascii-file",  required_argument, nullptr,'a' },
	{ "output-path",       required_argument, nullptr,'o'},
	{ "dataset",           required_argument, nullptr,'d'},
	{ "first-dataset",     required_argument, nullptr,'f'},
	{ "last-dataset",      required_argument, nullptr,'l'},
	{ "help",              no_argument,       nullptr,'h'},
	{nullptr, 0, nullptr, 0}
    };

    const char* const short_options = "i:n:N:a:o:d:f:l:h";
    int c;
    while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 )
	{
	    switch (c)
		{
		case 'i':
		    {
			inputdirs.emplace_back(optarg);
			break;
		    }
		case 'n':
		    {
			inputname.emplace_back(optarg);
			inputasciiname.emplace_back("");
			break;
		    }
		case 'N':
		    {
			inputhistoname.emplace_back(optarg);
			break;
		    }
		case 'a':
		    {
			inputname.emplace_back("");
			inputhistoname.emplace_back("");
			inputasciiname.emplace_back(optarg);
			break;
		    }
		case 'o':
		    {
			outpath = optarg;
			if( outpath.back() != '/' )
			    outpath.append("/");
			break;
		    }
		case 'd':
		    {
			ds.emplace_back(atoi(optarg));
			break;
		    }
		case 'f':
		    {
			firstds = atoi(optarg);
			break;
		    }
		case 'l':
		    {
			lastds = atoi(optarg);
			break;
		    }
		case'h':
		    {
			Usage();
			return 0;
		    }
		default:
		    {
			std::cout<<"def "<<optarg<<std::endl;
			exit(1);
		    }
		}
	}

    std::cout << "List of directories with efficiency files:" << std::endl;
    for( auto& d: inputdirs )
	{
	    if( d.back() != '/' )
		d.append("/");
	    std::cout << d << std::endl;
	}
    size_t nDirs = inputdirs.size();
    if( ds.size() == 0 )
	for( int d=firstds; d<=lastds; d++ )
	    ds.emplace_back(d);
    if( ds.size() == 0 )
	{
	    std::cout << "No dataset specified. Abort." << std::endl;
	    return 0;
	}
    
    std::cout << "List of datasets:" << std::endl;
    for( auto& d: ds )
	std::cout << d << "\t";
    std::cout << std::endl;
    
    size_t nDs = ds.size();
    int nBins = 1000;
    double minEff = 0.;
    double maxEff = 1.;
    std::vector<TH1D*> outputhisto;
    TF1* func = new TF1("func","gaus",0,1);
    func->SetNpx(10000);
    func->SetParameter(0,1);

    for( size_t d=0; d<nDs; d++ )
	{
	    std::string histoname = "CombinedEff_DS" + std::to_string(ds[d]);
	    outputhisto.emplace_back( new TH1D( histoname.c_str(), histoname.c_str(), nBins, minEff, maxEff ) );

	    std::vector<TFile*> rootfile;
	    std::vector<std::ifstream*> asciifile;
	    std::vector<TH1D*> inputhisto;
	    for( size_t dir=0; dir<nDirs; dir++ )
	    	{

		    if( inputname[dir].find(".root") != std::string::npos )// ROOT files
			{
			    std::string filename = inputdirs[dir] + "ds" + std::to_string(ds[d]) + "/" + inputname[dir];
			    std::ifstream testfile(filename.c_str());
			    if( !testfile.good() )
				{
				    std::cout << "File " << filename << " does not exist. Abort." << std::endl;
				    exit(0);
				}
			    
			    rootfile.emplace_back( new TFile(filename.c_str()) );
			    inputhisto.emplace_back( (TH1D*)(rootfile.back()->Get(inputhistoname[dir].c_str())) );
			    asciifile.emplace_back(new std::ifstream(""));
			}
		    else// ASCII files
			{
			    std::string filename = inputdirs[dir] + "ds" + std::to_string(ds[d])
				+ "/" + inputasciiname[dir] + "_ds" + std::to_string(ds[d]) + "_combined.txt";
			    
			    rootfile.emplace_back( new TFile() );
			    asciifile.emplace_back(new std::ifstream( filename.c_str() ));

			    std::string tmpS;			    
			    for( int i=0; i<3; i++ )
				std::getline( *(asciifile.back()), tmpS );
			    double tmpVal;
			    double eff, upper, lower;
			    for( int i=0; i<4; i++ )
				{
				    *(asciifile.back()) >> tmpS >> tmpVal;
				    if( tmpS == "fEff" )
					eff = tmpVal;
				    else if( tmpS == "fEffLower" )
					lower = tmpVal;
				    else if( tmpS == "fEffUpper" )
					upper = tmpVal;
				}

			    double errLower = eff - lower;
			    double errUpper = upper - eff;
			    double err = 0.5 * ( errLower + errUpper );

			    func->SetParameter(1,eff);
			    func->SetParameter(2,err);

			    std::string histoname = "PulserEff_ds" + std::to_string(ds[d]);
			    inputhisto.emplace_back( new TH1D( histoname.c_str(), histoname.c_str(), 30000, 0, 1 ) );
			    for( int i=0; i<1.e7; i++ )
				inputhisto.back()->Fill( func->GetRandom() );

			}
			/*
		    if( inputname[dir].find(".root") != std::string::npos )
			{
			    std::string filename = inputdirs[dir] + "ds" + std::to_string(ds[d]) + "/" + inputname[dir];
			    std::ifstream testfile(filename.c_str());
			    if( !testfile.good() )
				{
				    std::cout << "File " << filename << " does not exist. Abort." << std::endl;
				    exit(0);
				}
			    
			    rootfile.emplace_back( new TFile(filename.c_str()) );
			    inputhisto.emplace_back( (TH1D*)(rootfile.back()->Get(inputhistoname[dir].c_str())) );
			    asciifile.emplace_back(std::ifstream(""));
			}
		    else
			{
			    std::string filename = inputdirs[d] + "ds" + std::to_string(ds[d]) + "/" + inputasciiname[dir] + "_ds" + std::to_string(ds[d]) + "_combined.txt";
			    std::cout << "ASCII FILE: " << filename << std::endl;
			    rootfile.emplace_back( new TFile() );
			    asciifile.emplace_back(std::ifstream(filename.c_str()));
			    std::string tmpS;			    
			    for( int i=0; i<3; i++ )
				std::getline( asciifile.back(), tmpS );
			    double tmpVal;
			    double eff, upper, lower;
			    for( int i=0; i<4; i++ )
				{
				    asciifile.back() >> tmpS >> tmpVal;
				    std::cout << tmpS << "\t" << tmpVal;
				}
			}
		    */
	    	}


	    Combine( inputhisto, outputhisto[d] );

	    ComputeShortestInterval( ds[d], outputhisto[d] );
	    
	    for( auto& f: rootfile )
		f->Close();

	}

    std::string outputname = outpath + "CombinedEfficiencies.root";
    TFile* outputfile = new TFile( outputname.c_str(), "RECREATE" );
    outputfile->cd();
    for( auto& h: outputhisto )
	h->Write();
    outputfile->Close();

    
    /*
  std::vector<int>DS{3601,3602,3603,3604,3605,3606,3607,3607,3608,3609};

  TString path = "output/CUORE_eff";
  TH1D *heff;
  for (auto &D:DS)
    {
      heff = (TH1D*)Combine(Form("%s_P/ds%i/output_eff.root",path.Data(),D),
			    Form("%s_M/ds%i/output_eff.root",path.Data(),D));

      heff->SaveAs(Form("%s/efficiency_ds%i.root",path.Data(),D));

    }
    */
    return 0;
}

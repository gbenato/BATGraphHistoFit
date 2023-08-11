#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <utility>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "TFile.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include "TString.h"
#include "TCanvas.h"

void Usage()
{
  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  if( std::getenv("USER") != NULL )
    std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
  else
    std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
  std::cout<<"./CombineReso"<<std::endl;
  std::cout<<"options "<<std::endl;
  std::cout<<"------------------------------------------------------------"<<std::endl;
  std::cout<<"-i (or --input-path)        [directory for input] (default: none)"<<std::endl;
  std::cout<<"-o (or --output-path)       [directory for output file (default: output/CUORE_eff/] "<<std::endl;
  std::cout<<"-d (or --dataset)           [dataset number : default none] "<<std::endl;
  std::cout<<"-f (or --first-dataset)     [dataset number : default 0] "<<std::endl;
  std::cout<<"-l (or --last-dataset)      [dataset number : default 0] "<<std::endl;
  std::cout<<"-L (or --label)             [file label : default linearreso]"<<std::endl;
  std::cout<<"-h (or --help)              [this help]"<<std::endl;

}

int main(int argc, char **argv)
{
    std::string inputdir;
    std::string outpath;
    std::string label="linearreso";
    std::vector<int> ds;
    int firstds = 0;
    int lastds  = -1;
    static struct option long_options[] = {
	{ "input-path",        required_argument, nullptr,'i'},
	{ "output-path",       required_argument, nullptr,'o'},
	{ "dataset",           required_argument, nullptr,'d'},
	{ "first-dataset",     required_argument, nullptr,'f'},
	{ "last-dataset",      required_argument, nullptr,'l'},
	{ "label",             required_argument, nullptr,'L'},
	{ "help",              no_argument,       nullptr,'h'},
	{nullptr, 0, nullptr, 0}
    };

        const char* const short_options = "i:o:d:f:l:L:h";
    int c;
    while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 )
	{
	    switch (c)
		{
		case 'i':
		    {
			inputdir = optarg;
			if( inputdir.back() != '/' )
			    inputdir.append("/");
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
		case 'L':
		    {
			label = optarg;
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

    std::vector<TH1D*> biashisto;
    std::vector<TH1D*> resohisto;
    std::string histoname = "h1_fit_observable_fQ";
    for( size_t d=0; d<nDs; d++ )
	{
	    std::string biasfilename = inputdir + "ds" + std::to_string(ds[d]) + "/output_bias_" + label + ".root";
	    std::string resofilename = inputdir + "ds" + std::to_string(ds[d]) + "/output_reso_" + label + ".root";

	    TFile* biasfile = new TFile(biasfilename.c_str());
	    TFile* resofile = new TFile(resofilename.c_str());


	    biashisto.emplace_back( (TH1D*)biasfile->Get(histoname.c_str()) );
	    biashisto.back()->SetName( std::string( "bias_Qbb_ds" + std::to_string(ds[d]) ).c_str() );
	    biashisto.back()->SetTitle( std::string( "bias_Qbb_ds" + std::to_string(ds[d]) ).c_str() );
	    resohisto.emplace_back( (TH1D*)resofile->Get(histoname.c_str()) );
	    resohisto.back()->SetName( std::string( "reso_Qbb_ds" + std::to_string(ds[d]) ).c_str() );
	    resohisto.back()->SetTitle( std::string( "reso_Qbb_ds" + std::to_string(ds[d]) ).c_str() );
	    
	}

    std::string outputfilename = outpath + "CombinedReso.root";
    TFile* output = new TFile(outputfilename.c_str(),"RECREATE");
    output->cd();
    for( size_t d=0; d<nDs; d++ )
	{
	    biashisto[d]->Write();
	    resohisto[d]->Write();
	}
    output->Close();
    
    return 0;
}

#include <iostream>
#include <string>
#include <getopt.h>
#include <utility>
#include <vector>
#include <dirent.h>
#include <fstream>
#include <iomanip>

#include "TFile.h"
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"

void Usage()
{
    std::cout << std::endl << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    if( std::getenv("USER") != NULL )
	std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
    else
	std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
    std::cout<<"./PrintEfficiencyValues"<<std::endl;
    std::cout<<"options "<<std::endl;
    std::cout<<"------------------------------------------------------------"<<std::endl;
    std::cout<<"-i (or --input-dir)        [directory for input txt] (default: none)"<<std::endl;
    std::cout<<"-d (or --dataset)           [dataset number : default 0] "<<std::endl;
    std::cout<<"-n (or --histo-name)    [name of efficiency histograme] (default: eff_1461)"<<std::endl;
    std::cout<<"-h (or --help)              [this help]"<<std::endl;
    std::cout<<"------------------------------------------------------------"<<std::endl;
    
    return;
}

int main(int argc, char **argv)
{
    int ds=0;
    std::string inputdir;
    std::string histoname("eff_1461");
    static struct option long_options[] = {
	{ "input-dir",        required_argument,  nullptr, 'i' },
	{ "dataset",required_argument, nullptr,'d'},
	{ "histo-name", required_argument, nullptr, 'n' },
	{ "help",              no_argument,        nullptr,'h'},
	{nullptr, 0, nullptr, 0}
    };
    
    const char* const short_options = "i:d:n:h";
    int c;
    while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 )
	{
	    switch (c)
		{
		case 'i':
		    {
			inputdir = optarg;
			DIR *dir = opendir(inputdir.c_str());
			if( dir == nullptr )
			    {
				std::cout << "Directory " << inputdir << " does not exist. Abort." << std::endl;
				return 1;
			    }
			if( inputdir.back() != '/' )
			    inputdir.append("/");
			break;
		    }
		case 'd':
		    {
			ds=atoi(optarg);
			break;
		    }
		case 'n':
		    {
			histoname=optarg;
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

    //std::string filename = inputdir + "ds" + std::to_string(ds) + "/eff.root";
    std::string filename = inputdir + "ds" + std::to_string(ds) + "/output_eff.root";
    std::ifstream testfile(filename.c_str());
    if( !testfile.good() )
	{
	    std::cout << "File " << filename << " does not exist. Abort." << std::endl;
	    return 1;
	}

    TFile* file = new TFile(filename.c_str());
    file->cd();
    
    TH1D* heff = (TH1D*)file->Get(histoname.c_str());
    int nBins = heff->GetNbinsX();
    
    TH1D* shortest = new TH1D("Shortest","Shortest",nBins,heff->GetXaxis()->GetXmin(),heff->GetXaxis()->GetXmax());
    shortest->SetFillColor(kBlue-5);

    double quantile = 0.6826895;// Corresponds to +-1 sigma for a Gaussian distribution
    double binWidth = heff->GetBinWidth(1);
    int modeBin = heff->GetMaximumBin();
    double mode = heff->GetBinCenter(modeBin);
    int binLeft = modeBin;
    int binRight = modeBin;
    
    heff->Scale( 1./heff->Integral(1,nBins) );
    double integral = heff->GetBinContent(modeBin);
    while( integral < quantile )
	{
	    double intLeft  = heff->GetBinContent( binLeft - 1 );
	    double intRight = heff->GetBinContent( binRight + 1 );
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
    for( int b=binLeft; b<=binRight; b++ )
	shortest->SetBinContent( b, heff->GetBinContent(b) );
    
    double errShortLeft  = binWidth * ( modeBin - binLeft );
    double errShortRight = binWidth * ( binRight - modeBin );
    std::cout << "Mode:  " << std::setprecision(10) << mode << std::endl;
    std::cout << "Left:  " << std::setprecision(10) << binWidth*binLeft << std::endl;
    std::cout << "Right: " << std::setprecision(10) << binWidth*binRight << std::endl;
    std::cout << "Low error for shortest interval method:  " << std::setprecision(10) << errShortLeft << std::endl;
    std::cout << "High error for shortest interval method: " << std::setprecision(10) << errShortRight << std::endl;
    std::cout << "Coverage for equal-areas method: "
	      << heff->Integral( binLeft, binRight )
	      << std::endl << std::endl;

    TApplication *app = new TApplication("App", nullptr, nullptr);
    TCanvas* can = new TCanvas();
    can->cd();
    heff->Draw("hist");
    shortest->Draw("same");
    app->Run(kTRUE);

    
    return 0;
}

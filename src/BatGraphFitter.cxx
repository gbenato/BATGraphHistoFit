#include "BatGraphFitter.h"
#include "BAT_GraphFit.h"
#include "TFile.h"

#include "TCanvas.h"
#include <BAT/BCH1D.h>
#include "TH1D.h"
#include <BAT/BCParameter.h>
#include <BAT/BCModel.h>


BatGraphFitter::BatGraphFitter(TGraphAsymmErrors *&g,TF1 *&f)
{
  
    // constructor
    // in the constructor we need to create the
    //gaussian likelihood

    fMode="G";
    fTF1=f;
    fModel= new BAT_GraphFit("fit",f,1,fMode,fQbb);
    fGraph = g;
    fMin=-pow(10,9);
    fMax=pow(10,9);
}

BatGraphFitter::BatGraphFitter(TH1D*&h,TF1*&f)
{
    fMode="P";
    fHisto=h;
    fTF1=f;
    fModel= new BAT_GraphFit("fit",f,1,fMode,fQbb);
  
    this->SetPrecison(2);

}

BatGraphFitter::BatGraphFitter(TH1D*&h,TF1*&f,TF1 *&f_int,std::vector<double>energy)
{
    // Unbinned fit
    fMode="U";
    fHisto=h;
    fTF1=f;
    fTF1Int=f_int;
    fEnergyVector=energy;
    fModel= new BAT_GraphFit("fit",f,1,fMode,fQbb);
    fModel->SetTF1Int(f_int);
  
    this->SetPrecison(2);

}


BatGraphFitter::BatGraphFitter(TH1D *&h,TF1*&f,std::vector<std::pair<double,double>>ranges,std::vector<double>probs,double Sm)
{
    fMode="C";
    fHisto=h;
    fTF1=f;
    fModel= new BAT_GraphFit("fit",f,1,fMode,Sm);
	
    fSmax=Sm;
    this->SetCountingPar(ranges,probs);
    this->SetPrecison(2);

}

BatGraphFitter::BatGraphFitter(TGraph *&g,TF1*&f,double error)
{
    //constructor for graph
    fGraph=new TGraphAsymmErrors();
    fMode="G";
    std::cout<<"Warning : Attempt to fit a plain TGraph without errors - likelihood fit cannot be performed erros assumed constant (and 1):"<<std::endl;
    std::cout<<"You can set the value of the error to avoid numerical issies"<<std::endl;
    for (int i=0;i<g->GetN();i++)
	{
	    fGraph->SetPoint(i,g->GetX()[i],g->GetY()[i]);
	    fGraph->SetPointEYhigh(i,error);
	    fGraph->SetPointEYlow(i,error);
	}	

    fGraph->SetTitle(f->GetTitle());
    fMin=-pow(10,9);
    fMax=pow(10,9);
    fTF1=f;
    fModel= new BAT_GraphFit("fit",f,1,fMode,fQbb);

  
    this->SetPrecison(2);
}


BatGraphFitter::BatGraphFitter(TGraph *g,TGraph *&gTrial,TGraph *&gSuccess,TF1 *&f)
{
    fGraphBinomTrial = gTrial;
    fGraphBinomSuccess=gSuccess;
    if (gTrial->GetN()!=gSuccess->GetN())
	{
	    std::cout<<"ERROR number of points of trials not the same as successes"<<std::endl;
	}
    for (int i=0;i<gTrial->GetN();i++)
	{
	    fGraph->SetPoint(i,gTrial->GetX()[i],gSuccess->GetY()[i]/gTrial->GetY()[i]);
	    fGraph->SetPointEYhigh(i,1);
	    fGraph->SetPointEYlow(i,1);
	}
    fGraph->SetTitle(g->GetTitle());
    fMode="B";

    fModel= new BAT_GraphFit("fit",f,1,fMode,fQbb);

    this->SetPrecison(2);
    fTF1=f;
}



BatGraphFitter::BatGraphFitter(TGraphErrors *&g,TF1 *&f)
{
    // constructor for graph

    fMode="G";
    // in the constructor we need to create the
    fGraph=new TGraphAsymmErrors();
    for (int i=0;i<g->GetN();i++)
	{
	    fGraph->SetPoint(i,g->GetX()[i],g->GetY()[i]);
	    fGraph->SetPointEYhigh(i,g->GetEY()[i]);
	    fGraph->SetPointEYlow(i,g->GetEY()[i]);
	}
    fTF1=f;
    fGraph->SetTitle(g->GetTitle());
    fModel= new BAT_GraphFit("fit",f,1,fMode,fQbb);
    
    this->SetPrecison(2);

}

BatGraphFitter::~BatGraphFitter()
{
    if (fGraph!=nullptr)delete fGraph;
    if (fGraph!=nullptr)delete fGraph;
}




void BatGraphFitter::ResetTF1(TF1*&f)
{
    fModel->SetTF1(f,0);
}
void BatGraphFitter::ResetTF1Int(TF1*&fi)
{
    fModel->SetTF1Int(fi);
}

void BatGraphFitter::Fit( TString option,
			  TString goption,
			  double low,
			  double high,
			  bool marg,
			  bool quiet )
{

    auto start = std::chrono::high_resolution_clock::now();


    if (!quiet)
	std::cout<<"Mode is "<<fMode<<std::endl;
    bool CI=0;
    // parse the option - currently the only one to be implemented is R
    TString type;
    if (option.Contains("R")==1)
	{
	    fTF1->SetRange(low,high);
	}
  
    if (option.Contains("I")==1)
	{
	    type="C";
	}
    else
	{
	    type="P";
	}
    if (option.Contains("C")==1)
	{
	    // compute the confidence interval
	    CI = 1;
	}


    BCLog::SetLogLevelScreen(BCLog::error);

    // add some code to delete the old model (or reset it)

    if (fMode=="G")
	{
	    fModel->SetGraph(fGraph,fMax,fMin);
	}
    else if(fMode=="P"||fMode=="U")
	{
	    fHisto->GetXaxis()->SetRangeUser(low-10,high+10);
	    fModel->SetHisto(fHisto,type);
	    if (fMode=="U")
		{
		    fModel->SetEnergyVector(fEnergyVector);
		}
	}
    else if (fMode=="B")
	{

	    fModel->SetGraph(fGraph);
	    fModel->SetBinomGraph(fGraphBinomTrial,fGraphBinomSuccess);
	}
    else if (fMode=="C")
	{
	    fHisto->GetXaxis()->SetRangeUser(fCountingRange[0].first-10,fCountingRange[2].second+10);
	    fModel->SetHisto(fHisto,type);
	    fModel->SetCountingPars( fCountingRange,fCountingProb);
	}
  

    if (fPrec==1)
	{
	    fModel->SetPrecision(BCEngineMCMC::kLow);
	    fModel->SetNChains(2);
	}
    else if (fPrec==2)
	{
	    fModel->SetPrecision(BCEngineMCMC::kMedium);
	}
    else if (fPrec==3)
	{
	    fModel->SetPrecision(BCEngineMCMC::kHigh);
	}
    else if(fPrec>=4)
	{
	    fModel->SetPrecision(BCEngineMCMC::kVeryHigh);
	}
    
    if (marg==1)
	fModel->MarginalizeAll(BCIntegrate::kMargMetropolis);

    fModel->FindMode(fModel->GetBestFitParameters());


    if (!quiet)
	{
	    std::cout<<"                      "<<std::endl;
	    std::cout<<"Fit with BAT best fit value"<<std::endl;
	    std::cout<<"****************************************"<<std::endl;
      
	    for (int i=0;i<fModel->GetBestFitParameters().size(); i++)
		{
		    BCParameter  parameter = fModel->GetParameter(i);
		    fTF1->SetParameter(i,fModel->GetBestFitParameters()[i]);
		    std::cout<<"parameter "<<parameter.GetName().data()<<
			" "<<fModel->GetBestFitParameters()[i]<<std::endl;
	  
		}      std::cout<<" "<<std::endl;
    
      

	    auto stop = std::chrono::high_resolution_clock::now();
	    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	    std::cout<<"Fit took : "<<duration.count()<<" ms"<<std::endl;
	    std::cout<<" "<<std::endl;

	}
    // GET THE CREDIBLE INTERVAL
    if (CI)
	{
	    auto start = std::chrono::high_resolution_clock::now();

	    int n=500;
	    TGraphAsymmErrors *g1=new TGraphAsymmErrors(n);
	    TGraphAsymmErrors *g2=new	TGraphAsymmErrors(n);
	    TGraphAsymmErrors *g3=new	TGraphAsymmErrors(n);
	    TTree *T=(TTree*)fModel->GetMarkovChainTree();
	    this->GetCredibleInterval(g1,g2,g3,n,T,fMin,fMax);

	    for (int i=0;i<fModel->GetBestFitParameters().size(); i++)
		{
		    BCParameter  parameter = fModel->GetParameter(i);
		    fTF1->SetParameter(i,fModel->GetBestFitParameters()[i]);
		    fTF1->SetParError(i,fModel->GetMarginalized(parameter.GetName().data()).GetHistogram()->GetRMS());
		}
	    auto stop =std::chrono::high_resolution_clock::now();
	    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	    std::cout<<"Computing credible band took : "<<duration.count()<<" ms"<<std::endl;
	    std::cout<<" "<<std::endl;
	    delete T;
	    fModel->PlotCI(g1,g2,g3);
	    fGrint=g1;
	}
    else
	{
	    fModel->Plot(goption);
	}

  


}



void BatGraphFitter::GetCredibleInterval(TGraphAsymmErrors * &grint1,TGraphAsymmErrors *&grint2,TGraphAsymmErrors *&grint3,int n,TTree *&T_markov,double ylow,double yhigh)
{
    // method to get a Bayesian credible interval
    // we need to
    // 1) take the markov chain file and fill histograms
    // 2)


    // basic information we need
    // 1) fit low and high range

  
    double l,h;
    fTF1->GetRange(l,h);


    // get a range for the y-axis

    std::vector<double>pars;
    pars.resize(fTF1->GetNpar());
    int pass;

    for (int i=0;i<fTF1->GetNpar();i++)
	{
      
	    BCParameter  parameter = fModel->GetParameter(i);
	    T_markov->SetBranchAddress(parameter.GetName().data(),&pars[i]);
	}
    T_markov->SetBranchAddress("Phase",&pass);
  

    // create first a vector of TH1Ds

    std::vector<TH1D*> histos;
    for (int i=0;i<n;i++)histos.push_back(new TH1D(Form("h%i",i),Form("%i",i),10000,ylow,yhigh));

    for (long int i=0;i<T_markov->GetEntries()/10.;i++)
	{
	    T_markov->GetEntry(i);
	    if (i%1000000==0)
		{
		    std::cout<<100*double(i)/((double)T_markov->GetEntries())<<" %"<<std::endl;

		}
      
	    if (pass!=1){continue;}

	    // T_markov->Show(i);
	    for (int p=0;p<fTF1->GetNpar();p++)
		{
		    fTF1->SetParameter(p,pars[p]);

		}
	    double step = (h-l)/n;
	    for (int j=0;j<n;j++)
		{
		    // loop over histograms

		    double x = l+step/2.+step*j;

		    double y = fTF1->Eval(x);

		    histos[j]->Fill(y);
		    // std::cout<<"x = "<<x<<std::endl;
		    //  std::cout<<"y = "<<y<<std::endl;
		}
	}

    // now we have our n histograms - we need our c.i.
    double p0=0.683;
    double p1=0.955;
    double p2=0.9973;
    for (int j=0;j<n;j++)
	{

	    double step = (h-l)/n;

	    // for the c.i. we will do something a bit controversial use median and quantiles
	    double x = l+step/2.+step*j;

	    double q=0.5;
	    double point_est;
	    histos[j]->GetQuantiles(1,&point_est,&q);

	    grint1->SetPoint(j,x,point_est);
	    grint2->SetPoint(j,x,point_est);
	    grint3->SetPoint(j,x,point_est);


	    double up_error,down_error;

	    // ****** 1 sigma **********//
	    //**************************
	    q=0.5+p0/2;

	    histos[j]->GetQuantiles(1,&up_error,&q);
	    up_error-=point_est;

	    q=0.5-p0/2;
	    histos[j]->GetQuantiles(1,&down_error,&q);

	    //down_error+=point_est;
	    down_error=point_est-down_error;

	    grint1->SetPointError(j,0,0,down_error,up_error);

      
	    //***** 2 sigma ********** //
	    //*************************
	    q=0.5+p1/2;

	    histos[j]->GetQuantiles(1,&up_error,&q);
	    up_error-=point_est;

	    q=0.5-p1/2;
	    histos[j]->GetQuantiles(1,&down_error,&q);

	    down_error=point_est-down_error;

	    grint2->SetPointError(j,0,0,down_error,up_error);


	    // ******3 sigma***********
	    //***************************

	    q=0.5+p2/2;
	    histos[j]->GetQuantiles(1,&up_error,&q);
	    up_error-=point_est;
	    q=0.5-p2/2;
	    histos[j]->GetQuantiles(1,&down_error,&q);
	    down_error=point_est-down_error;
	    grint3->SetPointError(j,0,0,down_error,up_error);


      
	    histos[j]->GetXaxis()->SetRangeUser(point_est-down_error,point_est+up_error);
	    histos[j]->Draw();
     
	}
      

    grint1->SetLineColorAlpha(3,0.2);
    grint2->SetLineColorAlpha(kOrange-3,0.2);
    grint3->SetLineColorAlpha(kRed+1,0.2);
    grint1->SetFillColorAlpha(3,0.7);
    grint2->SetFillColorAlpha(kOrange-3,0.7);
    grint3->SetFillColorAlpha(kRed+1,0.7);


    for (int i=0;i<n;i++)delete histos[i];



}

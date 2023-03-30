#include <BAT/BCLog.h>
// ***************************************************************
// This file was created using the bat-project script
// for project BAT_GraphFit.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>

#include "BAT_GraphFit.h"
#include "TFile.h"
#include "TCanvas.h"
int main()
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);


    // create new BAT_GraphFit object

    // test it

    TFile *f = new TFile("../Cd113Shape/data_CWOnat.root");

    TF1 *func =(TF1*)f->Get("funeff");
    func->SetParameter(0,0.01);
    func->SetParLimits(0,0,0.2);
    func->SetParameter(1,50);
    func->SetParLimits(1,0,100);
    func->SetParameter(2,0.01);
    func->SetParLimits(2,-.1,0.2);
    TGraphAsymmErrors * g =(TGraphAsymmErrors*)f->Get("Graph_Efficiency");
    BAT_GraphFit m("Name_Me",func,g,1);
    
    // set precision
    m.SetPrecision(BCEngineMCMC::kVeryHigh);

    BCLog::OutSummary("Test model created");

    //////////////////////////////
    // perform your analysis here

    // Normalize the posterior by integrating it over the full parameter space
    // m.Normalize();

    // Write Markov Chain to a ROOT file as a TTree
    // m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // run mode finding; by default using Minuit
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

    // print summary plots
    // m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
    // m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
    // m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
    // m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

    // print results of the analysis into a text file
    m.PrintSummary();
    TCanvas *c = new TCanvas();
    c->cd();
    m.Plot();
    c->Print("test.pdf");
    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}

# BATGraphHistoFit
A tool to fit graphs and histograms with BAT


The idea of this code is to make it simple and easy to make Bayesian MCMC fits of Graphs and Histograms to functions.
The syntax follows closely root syntax for fitting TH1D and TGraph
3 modes are currently implemented with examples in TestMacro.cxx





To compile you can use:
g++ -std=c++0x `root-config  --cflags --glibs` macros/TestMacro.cxx src/BatGraphFitter.cxx src/BAT_GraphFit.cxx -o bin/TestMacro `bat-config --libs --cflags`

Or change the name of TestMacro.cxx to your macro.


To RUN the test do:
./bin/TestMacro


For the CUORE efficiencies compile and do:
./bin/CUORE_Efficiency




This will make a series of files in the output directory
out_pca_pass_${ENERGY}_plots.pdf - plots for the counting analysis for the peak at energy {ENERGY}
out_pca_fail_${ENERGY}_plots.pdf - same for events failing the cut

out_pca_pass.pdf/ out_pca_fail.pdf - plots of the fit reconstruction and posteriors on number of counts

eff.pdf                            - plots of the posteriors on efficiency and the graph fit
eff_fit.pdf                        - summary plots for the Graph fit

*most important*

output_eff.root                    - root file with all 1D and 2D posteriors for further analysis eg.


  KEY: TH1D	h1_fit_parameter_p0;1	
  KEY: TH2D	h2_fit_par_p1_vs_par_p0;1	
  KEY: TH2D	h2_fit_obs_fQ_vs_par_p0;1	
  KEY: TH1D	h1_fit_parameter_p1;1	
  KEY: TH2D	h2_fit_obs_fQ_vs_par_p1;1	
  KEY: TH1D	h1_fit_observable_fQ;1	




// ** BELOW IS AN EXPLANTION OF THE BAT CODE
// ---------------------------------------------------------------------------------------------------------



1) Fit a graph assuming gaussian errors
-----------------------------------------------------------------------------------------------

Syntax is something like:

  BatGraphFitter *fitter = new BatGraphFitter(g,f);
  fitter->SetPrecison(3);
  // sets the precion with 1 =kLow, 2 = kMedium, 3=kHigh, 4 = kVeryHigh (generally 2 or 3 are good choices)
  
  fitter->Fit("","",0,3000);

  ie the last Fit command has the same syntax as ROOT

  sometimes you have constraints eg the data are effieciencies or resoltution, to restrict to the physical region use
  fitter->SetGraphMaxMin(double max,double min);
  To limit the data to be within this range (multiply by this uniform distribution)






2) Fit a graph assuming a binomial error (for some efficiencies)
-----------------------------------------------------------------------------------------------

Constuct the fitter like

   BatGraphFitter(TGraph *g,TGraph *&gTrial,TGraph *&gSuccess,TF1 *&f);

WHere gTrial is the number of throws, gSuccess the number of catch and g the ratio (efficiency)





3) Fit a histogram assuming Poisson errors
-----------------------------------------------------------------------------------------------

Construct with
  BatGraphFitter *fitter = new BatGraphFitter(fHist,func);
  fitter->SetPrecison(3);
  fitter->Fit();




Currently accesing the outputs is in the "standard way" for BAT, I will try to add more methods soon

eg.
  // to print the marginalised posterior or write the markov chain do this
  fitter->fModel->PrintAllMarginalized("out_plots.pdf");
  fitter->fModel->WriteMarkovChain("output_mcmc.root", "RECREATE");

  // a plot of the best fit will automatically be produced and can be saves
  

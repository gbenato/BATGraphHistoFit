#include <iostream>

/*   2019.07.31 -- gfantini
 *   This version of the fitter is matched to the ROI BAT analysis
 *   for what concerns the parametrization of the polynomial scaling of the lineshape
 *   Tested successfully on the ds3021 lineshape data from early 2019
 *   -----------------------------------------------------------------------------------
 *   Contact: guido.fantini@gssi.it
 */

void BGResults(int ds=0,
	       const char* inputFolder = "fitresults",
	       const char* outputFolder= "" )
{
  if(ds == 0){
    cout << "Expect a valid dataset, will open file fitresults/fitresults_ds%d.dat" << endl;
    return;
  }
  ifstream in;
  vector<double> actualEnergy;
  vector<double> actualEnergyErr;
  vector<double> expectedEnergyForResiduals; // The Physics (true values)
  vector<double> expectedEnergyForResidualsErr; // 0. (or very very small)
  vector<double> expectedEnergyForSigma; // The Physics (true values)
  vector<double> expectedEnergyForSigmaErr; // 0. (or very very small)
  vector<double> residual;
  vector<double> residualErr;
  vector<double> widthFactor1;
  vector<double> widthFactor1Err;

  float ReadActualEnergy,ReadExpectedEnergy,ReadActualEnergyErr,ReadSigmaScale,ReadSigmaScaleErr;
  
  in.open(Form("%s/fitresults_ds%d.dat", inputFolder, ds));
  if( !in.good() ){
    cout << "ERROR: could not open file." << endl;
    return;
  }

  while (in >> ReadActualEnergy >> ReadExpectedEnergy >> ReadActualEnergyErr >> ReadSigmaScale >> ReadSigmaScaleErr){
    //peaks to be ignored removed directly from input files -- otherwise add conditions here 
    actualEnergy.push_back( ReadActualEnergy );
    actualEnergyErr.push_back( 0 ); // so that on the x axis there is no error ... this should change
    expectedEnergyForResiduals.push_back( ReadExpectedEnergy );
    expectedEnergyForResidualsErr.push_back(0.);
    residual.push_back( ReadExpectedEnergy - ReadActualEnergy );
    residualErr.push_back( ReadActualEnergyErr );

    expectedEnergyForSigma.push_back( ReadExpectedEnergy );
    expectedEnergyForSigmaErr.push_back(0.);
    widthFactor1.push_back( ReadSigmaScale );
    widthFactor1Err.push_back( ReadSigmaScaleErr );

    std::cout << "Loaded peak at:\t" << ReadExpectedEnergy << std::endl;
  }
  in.close();
  
  gStyle->SetStatX(0.47);
  gStyle->SetStatY(0.9);
  
  TVirtualFitter::SetDefaultFitter( "Minuit2" );


  // fitting the global scaling parameter
  TCanvas *c1 = new TCanvas();
  TGraphErrors *tge1 = new TGraphErrors( expectedEnergyForResiduals.size(),
					 &expectedEnergyForResiduals[0],
					 &residual[0],
					 &expectedEnergyForResidualsErr[0],
					 &residualErr[0] );
  tge1->SetMarkerColor(kBlack);
  tge1->SetMarkerStyle(22);
  TF1 *fun1 = new TF1("fun1","[0]+[1]*(x-[3])+[2]*(x-[3])*(x-[3])",0.,10.);
  //fun1->FixParameter( 3, 2527.518 );
  fun1->FixParameter( 3, 0. );
  TFitResultPtr FitResult1 = tge1->Fit("fun1", "EFS");
  ((TFitResult*)FitResult1.Get())->GetCorrelationMatrix().Print();
  // gfantini -- sanity check (TFitResultPtr is a smart pointer)
  if( (Int_t)FitResult1 != 0){
    std::cerr << "(PP) Fit did not converge!! Check code! Exiting . . ." << std::endl;
    std::cerr << "(PP) FitResult " << (Int_t)FitResult1 << std::endl;
    exit(0);
  }
  TMatrixDSym * TmpCovMat1 = &(((TFitResult*)FitResult1.Get())->GetCovarianceMatrix());
  TMatrixDSym CovMat1 = TmpCovMat1->GetSub(0,2,0,2);
  tge1->GetFunction("fun1")->SetLineColor(kBlack);
  double min = 835;
  min = 500;
  double max = 2614.511;
  TGraphErrors *tge1err = new TGraphErrors(101);
  for (int i = 0; i <= 100; i++)
    {
      double x    = min + i*(max - min)/100;
      double y    = fun1->Eval(x);
      double k    = fun1->GetParameter(3);
      double errx = 0.;
      double erry = sqrt( pow( fun1->GetParError(0), 2. ) +
			  pow( fun1->GetParError(1), 2. ) * pow( x-k, 2. ) +
			  pow( fun1->GetParError(2), 2. ) * pow( x-k, 4. ) +
			  2. * CovMat1( 0, 1 ) * ( x-k ) +
			  2. * CovMat1( 0, 2 ) * pow( x-k, 2. ) +
			  2. * CovMat1( 1, 2 ) * pow( x-k, 3. ) );
      tge1err->SetPoint(i, x, y);
      tge1err->SetPointError( i, errx, erry );
    }
  //(TVirtualFitter::GetFitter())->GetConfidenceIntervals(tge1err, 0.683);// --> This does NOT work!!
  tge1err->SetFillColor(kGray);
  tge1err->Draw("ape3");
  tge1->Draw("psame");
  tge1err->SetTitle( Form("Residual vs. energy, 1-sigma uncertainty (ds%d)",ds) );
  tge1err->GetYaxis()->SetRangeUser(-2,2);
  tge1err->GetYaxis()->SetTitle("Residual [keV]");
  tge1err->GetXaxis()->SetTitle("Energy [keV]");
  
  c1->SaveAs( Form("%s/residual_vs_energy_ds%d.pdf", outputFolder, ds) );
  c1->SaveAs( Form("%s/residual_vs_energy_ds%d.C", outputFolder, ds) );
  
  double Q = 2527.515;
  double errAtQ = sqrt( pow( fun1->GetParError(0), 2. ) +
			pow( fun1->GetParError(1), 2. ) * pow( Q-fun1->GetParameter(3), 2. ) +
			pow( fun1->GetParError(2), 2. ) * pow( Q-fun1->GetParameter(3), 4. ) +
			2. * CovMat1( 0, 1 ) * ( Q-fun1->GetParameter(3) ) +
			2. * CovMat1( 0, 2 ) * pow( Q-fun1->GetParameter(3), 2. ) +
			2. * CovMat1( 1, 2 ) * pow( Q-fun1->GetParameter(3), 3. ) );
  //(TVirtualFitter::GetFitter())->GetConfidenceIntervals (1, 1, Q, errAtQ, 0.683);// --> This does NOT work!!
  
  cout << "Dataset " << ds << endl;
  cout << "Shift at Q-value: " << fun1->Eval(Q) << " ± " << errAtQ << endl;

  // sanity check: is chi^2 computed OK?
  double chi2 = 0.;
  for( unsigned int ix = 0U; ix < expectedEnergyForResiduals.size(); ix++ )
    {
      // residual / residualErr / expectedEnergyForResiduals[0]
      double y = fun1->Eval( expectedEnergyForResiduals.at(ix) );
      double y0= residual.at(ix);
      double s = residualErr.at(ix);
      double X = (y-y0)*(y-y0)/(s*s);
      cout << "[Residuals] " << expectedEnergyForResiduals.at(ix) << "\t"
	   << (y-y0) << "\t" << s << "\t" << X << endl;
      chi2    += X;
    }
  cout << "[Residuals] chi2 manual computation: " << chi2 << endl;
  
  // fitting the sigma global scale factor with pol2 (Energy)
  TCanvas *c2 = new TCanvas();
  TGraphErrors *tge2 = new TGraphErrors( expectedEnergyForSigma.size(), 
					 &expectedEnergyForSigma[0],
					 &widthFactor1[0],
					 &expectedEnergyForSigmaErr[0],
					 &widthFactor1Err[0] );
  tge2->SetMarkerColor(kBlack);
  tge2->SetMarkerStyle(22);
  TF1 *fun2 = new TF1("fun2","[0]+[1]*x+[2]*x*x",0,10);
  fun2->FixParameter( 2, 0. );
  //fun2->SetParLimits(0,0.,2.);
  TFitResultPtr FitResult2 = tge2->Fit("fun2", "ES");
  tge2->GetFunction("fun2")->SetLineColor(kBlack);
  // gfantini -- sanity check (TFitResultPtr is a smart pointer)
  if( (Int_t)FitResult1 != 0){
    std::cerr << "(PP) Fit did not converge!! Check code! Exiting . . ." << std::endl;
    std::cerr << "(PP) FitResult " << (Int_t)FitResult1 << std::endl;
    exit(0);
  }
  TMatrixDSym CovMat2 = FitResult2->GetCovarianceMatrix();
  double min = 835;
  min = 500;
  double max = 2614.511;
  TGraphErrors *tge2err = new TGraphErrors(101);
  //for (int i = 0; i < 101; i++)
  //tge2err->SetPoint(i, min + i * (max - min) / 100, 0);

  for (int i = 0; i <= 100; i++)
    {
      double x    = min + i*(max - min)/100;
      double y    = fun2->Eval(x);
      double errx = 0.;
      double erry = sqrt( pow( fun2->GetParError(0), 2. ) +
			  pow( fun2->GetParError(1), 2. ) * pow( x, 2. ) +
			  pow( fun2->GetParError(2), 2. ) * pow( x, 4. ) +
			  2. * CovMat2( 0, 1 ) * x +
			  2. * CovMat2( 0, 2 ) * pow( x, 2. ) +
			  2. * CovMat2( 1, 2 ) * pow( x, 3. ) );
      tge2err->SetPoint(i, x, y);
      tge2err->SetPointError( i, errx, erry );
    }
  //(TVirtualFitter::GetFitter())->GetConfidenceIntervals(tge2err, 0.683);// --> This does NOT work!!
  tge2err->SetFillColor(kGray);
  tge2err->Draw("ape3");
  tge2->Draw("psame");
  tge2err->SetTitle( Form("Background resolution vs. energy, 1-sigma uncertainty (ds%d)",ds) );
  //tge2err->GetYaxis()->SetRangeUser(.1,2);
  tge2err->GetYaxis()->SetTitle("Effective resolution [keV]");
  tge2err->GetXaxis()->SetTitle("Energy [keV]");
  
  double errAtQ2 = sqrt( pow( fun2->GetParError(0), 2. ) +
			 pow( fun2->GetParError(1), 2. ) * pow( Q, 2. ) +
			 pow( fun2->GetParError(2), 2. ) * pow( Q, 4. ) +
			 2. * CovMat2( 0, 1 ) * Q +
			 2. * CovMat2( 0, 2 ) * pow( Q, 2. ) +
			 2. * CovMat2( 1, 2 ) * pow( Q, 3. ) );
  //(TVirtualFitter::GetFitter())->GetConfidenceIntervals (1, 1, Q, errAtQ2, 0.683);// --> This does NOT work!!
  
  cout << "Dataset " << ds << endl;
  cout << "Effective resolution at Q-value: " << fun2->Eval(2527.515) << " ± " << errAtQ2 << endl;

  // sanity check: is chi^2 computed OK?
  chi2 = 0.;
  for( unsigned int ix = 0U; ix < expectedEnergyForSigma.size(); ix++ )
    {
      // residual / residualErr / expectedEnergyForResiduals[0]
      double y = fun2->Eval( expectedEnergyForSigma.at(ix) );
      double y0= widthFactor1.at(ix);
      double s = widthFactor1Err.at(ix);
      double X = (y-y0)*(y-y0)/(s*s);
      cout << "[Sigma] " << expectedEnergyForSigma.at(ix) << "\t"
           << (y-y0) << "\t" << s << "\t" << X << endl;
      chi2    += X;
    }
  cout << "[Sigma] chi2 manual computation: " << chi2 << endl;
  
  c2->SaveAs( Form("%s/width_vs_energy_ds%d.pdf",outputFolder,ds) );
  c2->SaveAs( Form("%s/width_vs_energy_ds%d.C",outputFolder,ds) );
  // writing fit result and correlation matrix                                                                           
  string OutputFileName = Form("%s/residual_and_width_vs_energy_ds%d.root",outputFolder,ds);
  cout << "Plots, TF1s and covariance matrices will be stored to this file: "
       << OutputFileName.c_str() << endl;
  cout << "Use it to make the 0nbb BAT based fit with the scaling as NP." << endl;

  TFile* OutputFile = new TFile(OutputFileName.c_str(),"RECREATE");
  c1->Write("plot_Q");
  tge1->Write("residual_vs_energy");
  fun1->Write("fit_function_Q");
  FitResult1->Write("fit_result_Q");
  CovMat1.Write("covariance_matrix_Q");

  c2->Write("plot_sigma");
  tge2->Write("width_vs_energy");
  fun2->Write("fit_function_sigma");
  FitResult2->Write("fit_result_sigma");
  CovMat2.Write("covariance_matrix_sigma");
  OutputFile->Close();
  /*
    TGraphErrors *tge3 = new TGraphErrors(actualEnergy.size(), &actualEnergy[0], &widthFactor2[0], &actualEnergyErr[0], &widthFactor2Err[0]);
    tge3->SetMarkerColor(kBlack);
    tge3->SetMarkerStyle(22);
    TF1 *fun2 = new TF1("fun3","[0]+[1]*x+[2]*(2*x*x-1)",0,10);
    tge3->Fit("fun3");
    tge3->GetFunction("fun3")->SetLineColor(kBlack);
    double min = 835;
    double max = 2614.511;
    TGraphErrors *tge3err = new TGraphErrors(101);
    for (int i = 0; i < 101; i++)
      tge3err->SetPoint(i, min + i * (max - min) / 100, 0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(tge3err, 0.683);
    tge3err->SetFillColor(kGray);
    tge3err->Draw("ape3");
    tge3->Draw("psame");
    tge3err->SetTitle("Background resolution vs. energy, 1-#sigma fits (ds3518)");
    //tge2err->GetYaxis()->SetRangeUser(.1,2);
    tge3err->GetYaxis()->SetTitle("Effective resolution [keV]");
    tge3err->GetXaxis()->SetTitle("Energy [keV]");
    
    double errAtQ3[1] = {0};
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals (1, 1, Q, errAtQ3, 0.683);
    
    cout << "Effective resolution at Q-value: " << fun3->Eval(2527.515) << " ± " << *errAtQ3 << endl;
    
    c1->SaveAs("plots/width2_vs_energy.pdf");
    c1->SaveAs("plots/width2_vs_energy.C");
  */
}


double GetNoiseFWHM(const char* inputFolder,
		    int ds,
		    double noiseRangeKeV = 20.,
		    int noiseBins = 400,
		    const char* inputTree = "noise_pulser_tree" )
{
  // opening files for tower 1-19
  TChain* c = new TChain( inputTree );
  for( int tower=1; tower <= 19; tower++ )
    c->Add( Form("%s/ReducedNoisePulser_ds%d_T%03d.root",inputFolder,ds,tower) );
  std::cout << "Entries loaded: " << c->GetEntries() << std::endl;
  
  // producing histogram for noise events
  TH1D* hNoise = new TH1D("noise","noise events",noiseBins,-1.*noiseRangeKeV, noiseRangeKeV);
  c->Draw("Energy>>noise","IsNoise==1 && RBI && B4A");
  
  double hRMS = hNoise->GetRMS();

  // lorenzo moneta implements FWHM
  int bin1 = hNoise->FindFirstBinAbove(hNoise->GetMaximum()/2);
  int bin2 = hNoise->FindLastBinAbove(hNoise->GetMaximum()/2);
  double fwhm = hNoise->GetBinCenter(bin2) - hNoise->GetBinCenter(bin1);
  
  std::cout << "binned FWHM: " << fwhm << std::endl;
  return fwhm;
}

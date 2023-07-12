void Bias(TString path)
{

  // macro to plot the bias

  TChain *c= new TChain("test_stat");

  int N=  c->Add(path);
  std::cout<<"N = "<<N<<std::endl;
  TH2D * h_Gamma = new TH2D("h_Gamma","h_Gamma",500,0,1e-26,200,0,2e-26);

  c->Draw("(float_rate>0)*float_rate:rate>>h_Gamma");

  h_Gamma->SetTitle(" ; #Gamma_{inj} [yr^{-1}] ; #Gamma_{fit} [yr^{-1}]; ");
  h_Gamma->Fit("pol1");
  h_Gamma->Draw("colz");

  /// bias on bkg index

  TCanvas::MakeDefCanvas();
  TH2D * h_Bkg = new TH2D("h_Bkg","h_Bkg",500,0,1e-26,200,0,2e-4);

  c->Draw("(float_rate>0)*float_bkg+(float_rate<=0)*bkg_only_bkg:rate>>h_Bkg");

  h_Bkg->SetTitle(" ; #Gamma_{inj} [yr^{-1}] ; b_{fit} [counts/keV/kg/yr] ");
  TF1 *pol1_b = new TF1("pol1_b","pol1",0,1e-26);
  pol1_b->SetParameters(1e-4,0);
  h_Bkg->Fit(pol1_b);
  h_Bkg->Draw("colz");



}

  

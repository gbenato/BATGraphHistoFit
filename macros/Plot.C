TH1D* MedianError(TH1D*&h,double med)
{
  
  TH1D * hpost = new TH1D("hpost","hpost",1000,med*0.8,med*1.2);


  for (int i=0;i<10000;i++)
    {
      TH1D *htmp = new TH1D("htmp","htmp",200,0,5e27);

      for (int j=0;j<1e4;j++)
	{
	  htmp->Fill(h->GetRandom());
	}
      double q=0.5;
      double mode;
      htmp->GetQuantiles(1,&mode,&q);
      hpost->Fill(mode);
      delete htmp;
    }
  return hpost;
}
void mbb(TString name,bool use_eff)
{

  TChain * c1= new TChain("output");
  c1->Add(Form("output/CUPID_sens/bayesian_%s_simple_shape/limits/*",name.Data()));



  double M_shell_effective = 2.240;
  double M_shell_bare = 3.962;
  double M_IBM=4.22;
  double M_QRPA_1 = 3.9;
  double M_QRPA_2 = 5.402;
  double M_EDF=6.58;
  


  double G=15.92e-15;

  TH1D * h1= new TH1D("h1","h1",200,0,50);
  TH1D * h2= new TH1D("h2","h2",200,0,50);
  TH1D * h3= new TH1D("h3","h3",200,0,50);
  TH1D * h4= new TH1D("h4","h4",200,0,50);
  TH1D * h5= new TH1D("h5","h5",200,0,50);
  TH1D * h6= new TH1D("h6","h6",200,0,50);

  c1->Draw(Form("%f/sqrt(limit)>>h1",511*1e6*sqrt(1/G)/(M_shell_effective*1.27*1.27)));
  c1->Draw(Form("%f/sqrt(limit)>>h2",511*1e6*sqrt(1/G)/(M_shell_bare*1.27*1.27)));
  c1->Draw(Form("%f/sqrt(limit)>>h3",511*1e6*sqrt(1/G)/(M_IBM*1.27*1.27)));
  c1->Draw(Form("%f/sqrt(limit)>>h4",511*1e6*sqrt(1/G)/(M_QRPA_1*1.27*1.27)));
  c1->Draw(Form("%f/sqrt(limit)>>h5",511*1e6*sqrt(1/G)/(M_QRPA_2*1.27*1.27)));
  c1->Draw(Form("%f/sqrt(limit)>>h6",511*1e6*sqrt(1/G)/(M_EDF*1.27*1.27)));

  TLegend * l = new TLegend(0.5,0.5,0.9,0.9);
  l->AddEntry(h1,"Shell effective");
  l->AddEntry(h2,"Shell bare");
  l->AddEntry(h3,"IBM");
  l->AddEntry(h4,"QRPA 1");
  l->AddEntry(h5,"QRPA 2 ");
  l->AddEntry(h6,"EDF");


  h1->SetLineColor(kRed);
  h2->SetLineColor(kOrange-3);
  h3->SetLineColor(kGreen+1);
  h4->SetLineColor(kCyan+1);
  h5->SetLineColor(kBlue);
  h6->SetLineColor(kMagenta+1);

  h6->SetTitle(" ; m_{#beta#beta} ; probability [arb. units] ; ");

  
  TBox  *b  = new TBox(18.4,0,50,h6->GetMaximum());
  b->SetFillColorAlpha(9,0.3);
  
  h6->Draw();
  if (use_eff)
    h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");
  h5->Draw("same");
  h6->Draw("same");
  b->Draw();
  gStyle->SetOptStat(0);
  l->Draw();

    

    

  






}

void Plot(TString name,bool plot_full)
{

  TChain * c2= new TChain("output");
  c2->Add(Form("output/CUPID_sens/bayesian_%s_simple_shape/limits/*",name.Data()));

  TH1D * h2= new TH1D("h2","h2",200,0,5e27);
  c2->Draw("limit>>h2");


  gStyle->SetOptStat(0);


  double q=0.5;
  
  double mode2;

  h2->GetQuantiles(1,&mode2,&q);
  
  TLegend * l = new TLegend(0.5,0.5,0.9,0.9);
  //l->AddEntry(h1,"Full background shape");
  l->AddEntry(h2,"Distribution of limits");
  h2->SetTitle(" ; T_{1/2} 90% CI limit; counts ; ");
  h2->Draw("same");
  l->Draw();


  double maxi = 1.1*max(h2->GetMaximum(),h2->GetMaximum());

  h2->GetYaxis()->SetRangeUser(0,maxi);

  TLine *l2=new	TLine(mode2,0,mode2,maxi);
  l2->SetLineStyle(9);
  l2->SetLineColor(2);

  
  double low2,high2;
  q=0.16;
  h2->GetQuantiles(1,&low2,&q);
  l->AddEntry(l2,"Median sentivity");
  
  q=0.84;
  h2->GetQuantiles(1,&high2,&q);

  std::cout<<mode2<<" - "<<mode2-low2<<" + "<<high2-mode2<<std::endl;
  TBox * b = new TBox(low2,0,high2,maxi);
  b->SetFillColorAlpha(9,0.3);
  l->AddEntry(b,"68% of cases","F");

  h2->Draw();
  l2->Draw();
  b->Draw("same");
  l->Draw();
  
  TCanvas::MakeDefCanvas();

  TH1D* h2post=(TH1D*)MedianError(h2,mode2);

  h2post->SetTitle(" ; T_{1/2} median sensitivity [yrs] ; Probability [arb. units] ; ");
  h2post->Draw();
  mode2 = h2post->GetBinCenter(h2post->GetMaximumBin());

  double qlow=0.16;
  double low1,high1;
  double lowm2,highm2;
  double qhigh=0.84;
  
  h2post->GetQuantiles(1,&low2,&qlow);
  h2post->GetQuantiles(1,&high2,&qhigh);

  

  std::cout<<"Sens simple = "<<mode2<<" + "<<high2-mode2<<" - "<<mode2-low2<<std::endl;

}  

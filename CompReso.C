
void CompReso()
{

  std::vector<double>res_sqrt{1.04475,1.02225,1.01775,0.86925,0.95325};
  std::vector<double>ehigh_res_sqrt{0.0457493,0.0188829,0.0180573,0.0582138,0.0170038};
  std::vector<double>elow_res_sqrt{0.0327104,0.0170433,0.0154836,0.0491416,0.0437147};

  std::vector<double>res_pol1{1.06125,1.05675,1.20975,0.95625,1.06275};
  std::vector<double>ehigh_res_pol1{0.0636403,0.0545002,0.0501784,0.0845197,0.0479015};
  std::vector<double>elow_res_pol1{0.0625523,0.0488058,0.0477044,0.0782187,0.016812};

  std::vector<int>ds{3602,3603,3604,3605,3606,3608};
  

  TGraphAsymmErrors * g = new TGraphAsymmErrors();
  TGraphAsymmErrors * g2 = new TGraphAsymmErrors();

  g->SetTitle("Reso Scaling ; DS ; Efficiency ; ");
  for (int i=0;i<res_sqrt.size();i++)
    {
      g->SetPoint(i,ds[i],res_sqrt[i]);
      g->SetPointError(i,0,0,elow_res_sqrt[i],ehigh_res_sqrt[i]);
    }
  for (int i=0;i<res_pol1.size();i++)
    {
      g2->SetPoint(i,ds[i]+0.1,res_pol1[i]);
      g2->SetPointError(i,0,0,elow_res_pol1[i],ehigh_res_pol1[i]);
    }

  TLegend * l =new TLegend(0.7,0.7,.9,.9);
  l->AddEntry(g,"sqrt(pol1)","LP");
  l->AddEntry(g2,"pol1","LP");
  g2->SetLineColor(2);
  g2->SetMarkerColor(2);
  g->Draw("APE");
  g2->Draw("Psame");
  l->Draw();
}

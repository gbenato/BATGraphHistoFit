
void Comp()
{


  std::vector<double>pe_ridge{0.9419,0.9752,0.9596,0.9522,0.9609,0.9691,0.9646,0.9707,0.9472,0.9631};
  std::vector<double>low_ridge{ 0.0114, 0.0066 ,0.0058,0.0056 ,0.0136 ,0.0064, 0.0057,0.0048, 0.0056 ,0.0049 };
  std::vector<double>high_ridge{ 0.0114 , 0.0066,0.0058 ,0.0056,0.0136 ,0.0064, 0.0057,0.0048 , 0.0056 ,0.0049};
  std::vector<int>DS{3601,3602,3603,3604,3605,3606,3607,3608,3609,3610,3611,3612,3613,3614,3615};

  
  std::vector<double>pe_toby{0.88875,0.95325,0.97275,0.97125,0.99525,0.99525,0.96525,0.94275,0.96525,0.95325};
  std::vector<double>low_toby{0.029757,0.0236652,0.0147736,0.0153181,0.013254,0.011224,0.0158502,0.020563,0.0167204,0.0178818};
  std::vector<double>high_toby{0.031145,0.0248517,0.0180804,0.0171861,1-0.99525,1-0.99525,0.0194641,0.0215729,0.0185739,0.0200226};
  
  std::vector<double>pe_toby_flat{0.94875,0.97125,0.95025,0.95025,0.96975,0.96825,0.95625,0.96525,0.94575,0.95625};
  std::vector<double>low_toby_flat{0.00600697,0.00442919,0.00350722,0.00351003,0.00799782,0.00373897,0.0031744,0.00374553,0.00377446};
  std::vector<double>high_toby_flat{0.00870933,0.00642699,0.00643751,0.00647387,0.0105054,0.00706954,0.00674541,0.0077513,0.0061083,0.00641898};

  TGraphAsymmErrors * g = new TGraphAsymmErrors();
  TGraphAsymmErrors * g2 = new TGraphAsymmErrors();
  TGraphAsymmErrors * g3 = new TGraphAsymmErrors();

  g->SetTitle("AntiCoincidence efficiency ; DS ; Efficiency ; ");
  for (int i=0;i<pe_ridge.size();i++)
    {
      g->SetPoint(i,DS[i],pe_ridge[i]);
      g->SetPointError(i,0,0,low_ridge[i],high_ridge[i]);
    }
  for (int i=0;i<pe_toby.size();i++)
    {
      g2->SetPoint(i,DS[i]+0.1,pe_toby[i]);
      g2->SetPointError(i,0,0,low_toby[i],high_toby[i]);
    }

  for (int i=0;i<pe_toby_flat.size();i++)
    {
      g3->SetPoint(i,DS[i]+0.2,pe_toby_flat[i]);
      g3->SetPointError(i,0,0,low_toby_flat[i],high_toby_flat[i]);
    }

  TLegend * l =new TLegend(0.7,0.7,.9,.9);
  l->AddEntry(g,"Old method","LP");
  l->AddEntry(g2,"BAT fit pol1","LP");
  l->AddEntry(g3,"BAT fit pol0","LP");

  g3->SetLineColor(9);
  g3->SetMarkerColor(9);

  g2->SetLineColor(2);
  g2->SetMarkerColor(2);
  g->Draw("APE");
  g2->Draw("Psame");
  g3->Draw("Psame");
  l->Draw();
}


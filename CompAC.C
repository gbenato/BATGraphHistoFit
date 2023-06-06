
void Comp()
{


  std::vector<double>pe_ridge{0.9922,1.0032,1.0014,0.9917,0.9851,0.9937,0.9875,0.9992,0.9918,0.9883,0.9976,0.9931,0.9914,0.9876,0.9932};
  std::vector<double>low_ridge{0.0089,.0069,0.0050,0.0051,0.0129,0.0060,0.0061,0.0040 ,0.0046,0.0049,0.0048 ,0.0045,0.0040,0.0044,0.0038};
  std::vector<double>high_ridge{0.0087,.0072,0.0051,0.0050,0.0124,0.0060,0.0059,0.0040 ,0.0045,0.0048,0.0048 ,0.0044,0.0039,0.0043,0.0038};
  std::vector<int>DS{3601,3602,3603,3604,3605,3606,3607,3608,3609,3610,3611,3612,3613,3614,3615};

  
  std::vector<double>pe_toby{0.99225,0.99675,0.98025,0.98025,0.99675,0.98925,0.98775,0.99825,0.99075,0.98925};
  std::vector<double>low_toby{0.00616645,0.00365911,0.00470489,0.00470932,0.00572998,0.00573976,0.00455923,0.00175982,0.00414095,0.00449534};
  std::vector<double>high_toby{1-0.99225,1-0.99675,0.00841055,0.00851463,1-0.99675,0.00873683,0.00706297,1-0.99825,0.00631036,0.00636757};
  TGraphAsymmErrors * g = new TGraphAsymmErrors();
  TGraphAsymmErrors * g2 = new TGraphAsymmErrors();

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

  TLegend * l =new TLegend(0.7,0.7,.9,.9);
  l->AddEntry(g,"Old method","LP");
  l->AddEntry(g2,"BAT fit","LP");
  g2->SetLineColor(2);
  g2->SetMarkerColor(2);
  g->Draw("APE");
  g2->Draw("Psame");
  l->Draw();
}


void yield_simc(){
  string rootfile_dir="/u/group/c-csv/hdbhatt/yieldsimc2022/results/yield20binsptt_nov_widedelta";

  ifstream infile(Form("/u/group/c-csv/hdbhatt/my_analysis/da/rungroup_textfiles/rungrp_simc.txt"));
  if(!infile){
    cout<<"rungrp file not fount! Quiting!!"<<endl;
    exit(0);
  }
  infile.ignore(10000,'\n');

  double rungrp;
  string tgt, pol;
  vector<double>Rungrp;
  vector<string> Tgt;

  while(infile>>rungrp>>tgt){
    Rungrp.push_back(rungrp);
    Tgt.push_back(tgt);
  }
  int npt=Rungrp.size();
  cout<<npt<<endl;
  int nbin = 20;

  for(int ipt=0;ipt<npt;ipt++){
    TFile *fin=new TFile(rootfile_dir+"/yield_"+Tgt[ipt]+"_"+Rungrp[ipt]+"_simc.root");
    if(fin==NULL) {
      cout<<"rootfile doesn't exist"<<endl;
      continue;
    }
    if(Tgt[ipt]=="D2"){
      //define canvas
      TH1F *h_pos_inc_rad = (TH1F*)fin->Get("z_pos_inc_rad");
      TH1F *h_pos_inc_norad = (TH1F*)fin->Get("z_pos_inc_norad");
      TH1F *h_pos_exc = (TH1F*)fin->Get("z_pos_exc_rad");
      TH1F *h_pos_rho = (TH1F*)fin->Get("z_pos_rho");
      TH1F *h_pos_delta = (TH1F*)fin->Get("z_pos_delta");
      TH1F *h_pos_subtracted = (TH1F*)fin->Get("z_pos_inc_rad");
      h_pos_subtracted->Add(h_pos_rho,-1);
      h_pos_subtracted->Add(h_pos_delta,-1);
      h_pos_subtracted->Add(h_pos_exc,-1);
      h_pos_subtracted->Draw();
      TH1F *h_neg_inc_rad = (TH1F*)fin->Get("z_neg_inc_rad");
      TH1F *h_neg_inc_norad = (TH1F*)fin->Get("z_neg_inc_norad");
      TH1F *h_neg_exc = (TH1F*)fin->Get("z_neg_exc_rad");
      TH1F *h_neg_rho = (TH1F*)fin->Get("z_neg_rho");
      TH1F *h_neg_delta = (TH1F*)fin->Get("z_neg_delta");
      TH1F *h_neg_subtracted = (TH1F*)fin->Get("z_neg_inc_rad");
      h_neg_subtracted->Add(h_neg_rho,-1);
      h_neg_subtracted->Add(h_neg_delta,-1);
      h_neg_subtracted->Add(h_neg_exc,-1);
      h_neg_subtracted->Draw();
      double bin_content_pos, bin_content_neg, bin_error_pos, bin_error_neg;
      for(int ibin=1;ibin<=nbin;ibin++){
	bin_content_pos = h_pos_subtracted->GetBinContent(ibin);
	bin_content_neg = h_neg_subtracted->GetBinContent(ibin);
	bin_error_pos = h_pos_subtracted->GetBinError(ibin);
	bin_error_neg = h_neg_subtracted->GetBinError(ibin);
	cout<<Rungrp[ipt]<<"\t"<<Tgt[ipt]<<"\t"<<bin_content_pos<<"\t"<<bin_error_pos<<"\t"<<bin_content_neg<<"\t"<<bin_error_neg<<endl;
      }
      //canvas save
    }
    if(Tgt[ipt]=="H2"){
      //define canvas
      TH1F *h_pos_inc_rad = (TH1F*)fin->Get("z_pos_inc_rad");
      TH1F *h_pos_inc_norad = (TH1F*)fin->Get("z_pos_inc_norad");
      TH1F *h_pos_exc = (TH1F*)fin->Get("z_pos_exc_rad");
      TH1F *h_pos_rho = (TH1F*)fin->Get("z_pos_rho");
      TH1F *h_pos_delta = (TH1F*)fin->Get("z_pos_delta");
      TH1F *h_pos_subtracted = (TH1F*)fin->Get("z_pos_inc_rad");
      h_pos_subtracted->Add(h_pos_rho,-1);
      h_pos_subtracted->Add(h_pos_delta,-1);
      h_pos_subtracted->Add(h_pos_exc,-1);
      h_pos_subtracted->Draw();
      TH1F *h_neg_inc_rad = (TH1F*)fin->Get("z_neg_inc_rad");
      TH1F *h_neg_inc_norad = (TH1F*)fin->Get("z_neg_inc_norad");
      TH1F *h_neg_rho = (TH1F*)fin->Get("z_neg_rho");
      TH1F *h_neg_delta = (TH1F*)fin->Get("z_neg_delta");
      TH1F *h_neg_subtracted = (TH1F*)fin->Get("z_neg_inc_rad");
      h_neg_subtracted->Add(h_neg_rho,-1);
      h_neg_subtracted->Add(h_neg_delta,-1);
      h_neg_subtracted->Draw();
      double bin_content_pos, bin_content_neg, bin_error_pos, bin_error_neg;
      for(int ibin=1;ibin<=nbin;ibin++){
	bin_content_pos = h_pos_subtracted->GetBinContent(ibin);
	bin_content_neg = h_neg_subtracted->GetBinContent(ibin);
	bin_error_pos = h_pos_subtracted->GetBinError(ibin);
	bin_error_neg = h_neg_subtracted->GetBinError(ibin);
	cout<<Rungrp[ipt]<<"\t"<<Tgt[ipt]<<"\t"<<bin_content_pos<<"\t"<<bin_error_pos<<"\t"<<bin_content_neg<<"\t"<<bin_error_neg<<endl;
      }
      //canvas save
    }
  }
}

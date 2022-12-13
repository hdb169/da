void yield_data(double myrungrp, string mypol, string kingrp){
  TString textfile_dir = "/u/group/c-csv/hdbhatt/my_analysis/da";
  TString data_file_dir = "/u/group/c-csv/hdbhatt/yielddec21/yield_skim_nov_widedelta";
  TString simc_file_dir = "/u/group/c-csv/hdbhatt/yieldsimc2022/results/yield20binsptt_nov_widedelta";

  ifstream infile(Form("%s/%s.txt",textfile_dir.Data(),kingrp.c_str()));
  if(!infile){
    cout<<"Text file doesn't exist! Quiting!!"<<endl;
    exit(0);
  }
  infile.ignore(10000,'\n');
  double rungrp, q2, x, z, hmsp, shmsp, runnum;
  string tgt, pol;
  vector<double>Rungrp, Q2, X, Z, Hmsp, Shmsp, Runnum;
  vector<string>Tgt, Pol;
  while(infile>>rungrp>>tgt>>q2>>x>>z>>hmsp>>shmsp>>pol>>runnum){
    Rungrp.push_back(rungrp);
    Tgt.push_back(tgt);
    Q2.push_back(q2);
    X.push_back(x);
    Z.push_back(z);
    Hmsp.push_back(hmsp);
    Shmsp.push_back(shmsp);
    Pol.push_back(pol);
    Runnum.push_back(runnum);
  }
  int npt = Rungrp.size();
  cout<<npt<<endl;
  int count=0;
  int countAl=0;
  const int nbin=20;
  double wavg[nbin]={0};
  double werror[nbin]={0};
  double wi[nbin]={0};
  double ws[nbin]={0};
  double ww[nbin]={0};
  double wavgAl[nbin]={0};
  double werrorAl[nbin]={0};
  double wiAl[nbin]={0};
  double wsAl[nbin]={0};
  double wwAl[nbin]={0};
  double wavgCor[nbin]={0};
  double werrorCor[nbin]={0};
  double dummyFacD2=0.254;
  double dummyFacH2=0.262;
  double corrFac=0;
  
  for(int ipt=0;ipt<npt;ipt++){
    if(Rungrp[ipt]==myrungrp && Pol[ipt]==mypol){
      if(Tgt[ipt]!="Al"){
        count++;
}
      if(Tgt[ipt]=="Al"){
        countAl++;
      }
      double bincontent[count][nbin];
      double binerror[count][nbin];
      double bincontentAl[countAl][nbin];
      double binerrorAl[countAl][nbin];

      if(Tgt[ipt]!="Al"){
        TFile *dfile = new TFile(data_file_dir+"/yield_"+Runnum[ipt]+".root","READ");
        if(dfile->IsZombie()) continue;
        TH1D *hy = (TH1D*)dfile->Get("h_zhad_clean_mom_counts");
        for(int icount=0;icount<count;icount++){
          for(int ibin=1;ibin<=nbin;ibin++){
            bincontent[icount][ibin] = hy->GetBinContent(ibin);
            binerror[icount][ibin] = hy->GetBinError(ibin);
//            cout<<Runnum[ipt]<<"\t"<<ibin<<"\t"<<bincontent[icount][ibin]<<"\t"<<binerror[icount][ibin]<<endl;
            if(bincontent[icount][ibin]>0.0 && binerror[icount][ibin]<0.9*bincontent[icount][ibin]){//if binerror more than 90 %, ignore the bin.
               wi[ibin]=1.0/(binerror[icount][ibin]*binerror[icount][ibin]);
               ws[ibin]+=bincontent[icount][ibin]*wi[ibin];
               ww[ibin]+=1.0/(binerror[icount][ibin]*binerror[icount][ibin]);
            }
          }
        }
      }
      if(Tgt[ipt]=="Al"){
        TFile *dfile = new TFile(data_file_dir+"/yield_"+Runnum[ipt]+".root","READ");
        if(dfile->IsZombie()) continue;
        TH1D *hy = (TH1D*)dfile->Get("h_zhad_clean_mom_counts");
        for(int icount=0;icount<countAl;icount++){
          for(int ibin=1;ibin<=nbin;ibin++){
            bincontentAl[icount][ibin] = hy->GetBinContent(ibin);
            binerrorAl[icount][ibin] = hy->GetBinError(ibin);
//            cout<<Runnum[ipt]<<"\t"<<ibin<<"\t"<<bincontentAl[icount][ibin]<<"\t"<<binerrorAl[icount][ibin]<<endl;
            if(bincontentAl[icount][ibin]>0.0 && binerrorAl[icount][ibin]<0.2*bincontentAl[icount][ibin]){//if binerror more than 20 %, ignore the bin.
               wiAl[ibin]=1.0/(binerrorAl[icount][ibin]*binerrorAl[icount][ibin]);
               wsAl[ibin]+=bincontentAl[icount][ibin]*wiAl[ibin];
               wwAl[ibin]+=1.0/(binerrorAl[icount][ibin]*binerrorAl[icount][ibin]);
            }
          }
        }
      }
      count=0;
      countAl=0;
    }
  }
  for(int ipt=0;ipt<npt;ipt++){
    if(Rungrp[ipt]==myrungrp && Pol[ipt]==mypol){
      if(Tgt[ipt]!="Al"){
        if(Tgt[ipt]=="D2") corrFac=dummyFacD2;
        if(Tgt[ipt]=="H2") corrFac=dummyFacH2;
        for(int ibin=1;ibin<=nbin;ibin++){
          if(ww[ibin]>0.0 && wwAl[ibin]>0.0){
            wavg[ibin]=ws[ibin]/ww[ibin];
            werror[ibin]=sqrt(1.0/ww[ibin]);
            wavgAl[ibin]=wsAl[ibin]/wwAl[ibin];
            werrorAl[ibin]=sqrt(1.0/wwAl[ibin]);
            wavgCor[ibin]=wavg[ibin]-corrFac*wavgAl[ibin];
            werrorCor[ibin]=sqrt(werror[ibin]*werror[ibin]+corrFac*corrFac*werrorAl[ibin]*werrorAl[ibin]);
            cout<<Tgt[ipt]<<"\t"<<corrFac<<ibin<<"\t"<<wavg[ibin]<<"\t"<<werror[ibin]<<"\t"<<wavgAl[ibin]<<"\t"<<werrorAl[ibin]<<"\t"<<wavgCor[ibin]<<"\t"<<werrorCor[ibin]<<endl;
          }
        }
      }
    }
  }
}

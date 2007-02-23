
void plotall(){
  gROOT->Reset();
  gROOT->SetBatch(kTRUE);//don't plot to screen, save mungo time!

  TFile myFile("MonitorTrigT1BS.root");  
  myFile.cd();  
  
  myCanvas = new TCanvas("plots","Plots");
  myCanvas->SetFillColor(10);//white
  myPsFile = new TPostScript("AllPlots.ps",112);
  
  TDirectory *current_sourcedir = gDirectory;
  TIter nextkey(current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  TObject *obj;// = key->ReadObj(); 
  TH1 *h1;// = (TH1*)obj; 
  string nameof ;

  cout<<"loop"<<endl<<endl;

  for(int j=1; j<current_sourcedir->GetNkeys()+1 ;j+=4){
  myPsFile->NewPage();
  myCanvas->Divide(2,2);

    for(int i=1; i<5;i++){
      myCanvas->cd(i);//1->4
      key = (TKey*)nextkey(); 
      nameof = key->GetName();
      cout<<"Key Name = "<<nameof<<endl;
      obj = key->ReadObj(); 
      h1 = (TH1*)obj; 
      h1->Draw();
      gPad->Update();
    }
    myCanvas->Update(); 
    // break;
  }

  cout<<"eo loooop"<<endl<<endl;
  myPsFile->Close();

  /*
//to run in root do:
.L Plotter.cxx
plotall()
  */
}

void prettyplot(){
  gROOT->Reset();
  gROOT->SetBatch(kTRUE);//don't plot to screen, save mungo time!

  Dplots(); 
  CaloComparisonPlots();
  JEControlPlots();
  TTControlPlots();
  TSCTPlots();
  CalibrationPlots();


  /*
//to run in root do:
.L Plotter.cxx
prettyplot()
  */
}

void selectplots(){

  TFile myFile("MonitorTrigT1BS.root");  
  myFile.cd();  
  
  myCanvas = new TCanvas("plots","Plots");
  myCanvas->SetFillColor(10);//white
  myPsFile = new TPostScript("TalkPlots.ps",112);

  TDirectory *current_sourcedir = gDirectory;
  
  //  TLegend* Legend = new TLegend(0.1,0.1,0.1,0.1,NULL,"brNDC");
  gStyle->SetOptStat(0);//Get rid of legend 
  //Page 1
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  //  TT_Calo_EtTower->SetStats(kFALSE);
  TT_Calo_Em_EtTower->SetStats(kFALSE);
  TT_Calo_Em_EtTower->Draw();
  gPad->Update();
  myCanvas->cd(2);
  Calo_Em_Et->Draw("e");
  gPad->Update();
  
  myCanvas->Update(); 

  myPsFile->Close();



  /*
//to run in root do:
.L Plotter.cxx
selectplots()
  */
}

void Dplots(){
  //Plot the distributions for Towers where the Ratio of 
  //TT to TriggerStyleCaloTower Et (Had and Em seperately) 
  //is over or under a cutoff

  gROOT->SetBatch(kTRUE);//don't plot to screen, save mungo time!

  TFile myFile("MonitorTrigT1BS.root");  
  myFile.cd();  
  
  myCanvas = new TCanvas("plots","Plots");
  myCanvas->SetFillColor(10);//white
  myPsFile = new TPostScript("UnmatchedTowerPlots.ps",112);

  TDirectory *current_sourcedir = gDirectory;

  //Plot ratio distributions (TSCT/TT)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Ratio_D_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Ratio_D_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  // ==================================
  //  Em Et Divergence Plots
  // ==================================

  //Plot Et (Em and Had) of TT and TSCT 
  // TT[Em]<TSCT[Em] (i.e. under)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_DEm_under_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_DEm_under_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  Calo_DEm_under_TTEm_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(4);
  Calo_DEm_under_TTHad_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  // Plot TSCT eta and phi 
  // TT[Em]<TSCT[Em] (i.e. under)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_DEm_under_phi->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_DEm_under_eta->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Plot Et (Em and Had) of TT and TSCT 
  // TT[Em]>TSCT[Em] (i.e. over)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_DEm_over_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_DEm_over_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  Calo_DEm_over_TTEm_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(4);
  Calo_DEm_over_TTHad_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  // Plot TSCT eta and phi 
  // TT[Em]>TSCT[Em] (i.e. over)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_DEm_over_phi->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_DEm_over_eta->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //===================================
  // Had Et Divergence Plots
  //===================================

  //Plot Et (Em and Had) of TT and TSCT 
  // TT[Had]<TSCT[Had] (i.e. under)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_DHad_under_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_DHad_under_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  Calo_DHad_under_TTEm_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(4);
  Calo_DHad_under_TTHad_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  // Plot TSCT eta and phi 
  // TT[Had]<TSCT[Had] (i.e. under)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_DHad_under_phi->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_DHad_under_eta->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Plot Et (Em and Had) of TT and TSCT 
  // TT[Had]>TSCT[Had] (i.e. over)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_DHad_over_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_DHad_over_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  Calo_DHad_over_TTEm_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(4);
  Calo_DHad_over_TTHad_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  // Plot TSCT eta and phi 
  // TT[Had]>TSCT[Had] (i.e. over)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_DHad_over_phi->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_DHad_over_eta->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  myPsFile->Close();



  /*
//to run in root do:
.L Plotter.cxx
Dplots()
  */
}



void CaloComparisonPlots(){
  //Plot the distributions for Towers where the Ratio of 
  //TT to TriggerStyleCaloTower Et (Had and Em seperately) 
  //is over or under a cutoff

  gROOT->SetBatch(kTRUE);//don't plot to screen, save mungo time!

  TFile myFile("MonitorTrigT1BS.root");  
  myFile.cd();  
  
  myCanvas = new TCanvas("plots","Plots");
  myCanvas->SetFillColor(10);//white
  myPsFile = new TPostScript("CaloTowerControlPlots.ps",112);

  TDirectory *current_sourcedir = gDirectory;
 
  //Plot Et, eta and phi distributions (Et>0GeV)
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  CaloTower_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  CaloTower_phi->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  CaloTower_eta->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Plot Et, eta and phi distributions for Towers with Et>10GeV
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  CaloTower_Et10->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  CaloTower_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  CaloTower_eta_gt10->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Plot eta-phi distribution and hitmap
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  CaloTower_Eta_Phi->Draw();
  gPad->Update();
  myCanvas->cd(2);
  CaloTower_Eta_Phi_HM->Draw();
  gPad->Update();
  myCanvas->Update(); 

  myPsFile->Close();

  /*
//to run in root do:
.L Plotter.cxx
CaloComparisonPlots()
  */
}

void JEControlPlots(){
  //Plot JE distributions

  gROOT->SetBatch(kTRUE);//don't plot to screen, save mungo time!

  TFile myFile("MonitorTrigT1BS.root");  
  myFile.cd();  
  
  myCanvas = new TCanvas("plots","Plots");
  myCanvas->SetFillColor(10);//white
  myPsFile = new TPostScript("JEControlPlots.ps",112);

  TDirectory *current_sourcedir = gDirectory;

  //Plot EmEt HadEt, phi and eta
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  JE_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  JE_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  JE_eta->Draw("e");
  gPad->Update();
  myCanvas->cd(4);
  JE_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Plot Total Et (Had+Em), EmEt in the range 0-10, HadEt in the range 0-10 
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  JE_Tot_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  JE_EM10_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  JE_HAD10_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Calo region plots

  //Barrel EmEt, HadEt, phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Barrel_JE_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Barrel_JE_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  Barrel_JE_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Barrel Et in range 0-10
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Barrel10_JE_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Barrel10_JE_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //EC EmEt, HadEt, phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  EC_JE_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  EC_JE_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  EC_JE_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //EC Et in range 0-10
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  EC10_JE_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  EC10_JE_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //FCAL EmEt, HadEt, phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  FCAL_JE_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  FCAL_JE_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  FCAL_JE_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //FCAL Et in range 0-10
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  FCAL10_JE_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  FCAL10_JE_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //eta-phi distribution

  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Eta_Phi_JE->Draw();
  gPad->Update();
  myCanvas->cd(2);
  Eta_Phi_HM_JE->Draw();
  gPad->Update();
  myCanvas->Update(); 

  //For Et>10

  //eta and phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  JE_eta_gt10->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  JE_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //phi in calo regions

  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  JE_Barrel_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  JE_EC_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  JE_FCAL_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  myPsFile->Close();

  /*
//to run in root do:
.L Plotter.cxx
JEControlPlots()
  */
}

void TTControlPlots(){
  //Plot TT distributions

  gROOT->SetBatch(kTRUE);//don't plot to screen, save mungo time!

  TFile myFile("MonitorTrigT1BS.root");  
  myFile.cd();  
  
  myCanvas = new TCanvas("plots","Plots");
  myCanvas->SetFillColor(10);//white
  myPsFile = new TPostScript("TTControlPlots.ps",112);

  TDirectory *current_sourcedir = gDirectory;

  //Plot EmEt HadEt, phi and eta
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  TT_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  TT_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  TT_eta->Draw("e");
  gPad->Update();
  myCanvas->cd(4);
  TT_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Plot Total Et (Had+Em), EmEt in the range 0-10, HadEt in the range 0-10 
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  TT_Tot_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  TT_EM10_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  TT_HAD10_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Calo region plots

  //Barrel EmEt, HadEt, phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Barrel_TT_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Barrel_TT_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  Barrel_TT_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Barrel Et in range 0-10
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Barrel10_TT_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Barrel10_TT_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //EC EmEt, HadEt, phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  EC_TT_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  EC_TT_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  EC_TT_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //EC Et in range 0-10
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  EC10_TT_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  EC10_TT_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //FCAL EmEt, HadEt, phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  FCAL_TT_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  FCAL_TT_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  FCAL_TT_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //FCAL Et in range 0-10
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  FCAL10_TT_EM_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  FCAL10_TT_HAD_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //eta-phi distribution

  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Eta_Phi_TT->Draw();
  gPad->Update();
  myCanvas->cd(2);
  Eta_Phi_HM_TT->Draw();
  gPad->Update();
  myCanvas->Update(); 

  //For Et>10

  //eta and phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  TT_eta_gt10->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  TT_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //phi in calo regions

  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  TT_Barrel_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  TT_EC_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  TT_FCAL_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  myPsFile->Close();

  /*
//to run in root do:
.L Plotter.cxx
TTControlPlots()
  */
}

void TSCTPlots(){
  //Plot distributions for Trigger Style Calo Towers
  //i.e. Trigger Towers created directly from ESD level CaloCells
  //also plot comparisons between them and TTs.

  gROOT->SetBatch(kTRUE);//don't plot to screen, save mungo time!

  TFile myFile("MonitorTrigT1BS.root");  
  myFile.cd();  
  
  myCanvas = new TCanvas("plots","Plots");
  myCanvas->SetFillColor(10);//white
  myPsFile = new TPostScript("TSCTPlots.ps",112);

  //phi, eta, and phi and eta for EmEt>10GeV
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_phi->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_eta->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  Calo_phi_gt10->Draw("e");
  gPad->Update();
  myCanvas->cd(4);
  Calo_eta_gt10->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Et: Em, Had and in 0-10GeV range
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Calo_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  Calo_Em_Et10->Draw("e");
  gPad->Update();
  myCanvas->cd(4);
  Calo_Had_Et10->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Eta-Phi distributions
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Calo_Eta_Phi->Draw();
  gPad->Update();
  myCanvas->cd(2);
  Calo_Eta_Phi_HM->Draw();
  gPad->Update();
  myCanvas->Update(); 

  //calo region plots

  //Barrel Et Had, Em and phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Barrel_Calo_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Barrel_Calo_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  Barrel_Calo_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Barrel Ets in 0-10GeV range
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  Barrel10_Calo_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  Barrel10_Calo_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //EC Et Had, Em and phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  EC_Calo_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  EC_Calo_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  EC_Calo_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //EC Ets in 0-10GeV range
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  EC10_Calo_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  EC10_Calo_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //FCAL Et Had, Em and phi
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  FCAL_Calo_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  FCAL_Calo_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  FCAL_Calo_phi->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //FCAL Ets in 0-10GeV range
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  FCAL10_Calo_Em_Et->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  FCAL10_Calo_Had_Et->Draw("e");
  gPad->Update();
  myCanvas->Update(); 

  //Comparison plots TSCT with TT
  //Em and Had Ets, phi and eta
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  TT_Calo_Em_EtTower->Draw();
  gPad->Update();
  myCanvas->cd(2);
  TT_Calo_Had_EtTower->Draw();
  gPad->Update();
  myCanvas->cd(3);
  TT_Calo_EtaTower->Draw();
  gPad->Update();
  myCanvas->cd(4);
  TT_Calo_PhiTower->Draw();
  gPad->Update();
  myCanvas->Update(); 

  myPsFile->Close();

  /*
//to run in root do:
.L Plotter.cxx
TSCTPlots()
  */
}


void CalibrationPlots(){
  //Plot TT Calibration plots

  gROOT->SetBatch(kTRUE);//don't plot to screen, save mungo time!

  TFile myFile("MonitorTrigT1BS.root");  
  myFile.cd();  
  
  myCanvas = new TCanvas("plots","Plots");
  myCanvas->SetFillColor(10);//white
  myPsFile = new TPostScript("CalibrationPlots.ps",112);

  //Plotting Et vs Eta for TTEmEt TTHadEt and Trigger Style Calo Towers
  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  TTEM_EtEta->Draw();
  gPad->Update();
  myCanvas->cd(2);
  TTHAD_EtEta->Draw();
  gPad->Update();
  myCanvas->cd(3);
  CaloEM_EtEta->Draw();
  gPad->Update();
  myCanvas->cd(4);
  CaloHAD_EtEta->Draw();
  gPad->Update();
  myCanvas->Update(); 

  //Em vs eta for CaloTowers from ESD
  myPsFile->NewPage();
  myCanvas->cd(1);
  CaloT_EtEta->Draw();
  gPad->Update();
  myCanvas->Update(); 

  myPsFile->Close();

  /*
//to run in root do:
.L Plotter.cxx
CalibrationPlots()
  */
}





void Plots(){
  //Plot 

  gROOT->SetBatch(kTRUE);//don't plot to screen, save mungo time!

  TFile myFile("MonitorTrigT1BS.root");  
  myFile.cd();  
  
  myCanvas = new TCanvas("plots","Plots");
  myCanvas->SetFillColor(10);//white
  myPsFile = new TPostScript("Plots.ps",112);

  myPsFile->NewPage();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  ->Draw("e");
  gPad->Update();
  myCanvas->cd(2);
  ->Draw("e");
  gPad->Update();
  myCanvas->cd(3);
  ->Draw("e");
  gPad->Update();
  myCanvas->cd(4);
  ->Draw("e");
  gPad->Update();
  myCanvas->Update(); 


  myPsFile->Close();

  /*
//to run in root do:
.L Plotter.cxx
Plots()
  */
}


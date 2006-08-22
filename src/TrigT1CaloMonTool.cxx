// ********************************************************************
//
// NAME:     TrigT1CaloMonTool.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Ethan-Etienne Woehrling (eew@hep.ph.bham.ac.uk)
//           
//
// ********************************************************************

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ITHistSvc.h"

#include <TH1D.h>
#include <TH2D.h>

#include "TString.h"

#include "StoreGate/StoreGateSvc.h"

#include "CaloEvent/CaloTower.h"
#include "CaloEvent/CaloTowerContainer.h"

#include "TrigT1CaloMonitoring/TrigT1CaloMonTool.h"

#include "CLHEP/Units/SystemOfUnits.h"

//to open TT & JE Persistent objects
#include "TrigT1Calo/TriggerTowerCollection.h"
#include "TrigT1Calo/JetElementCollection.h"
#include "TrigT1Calo/TrigT1CaloDict.h"
#include "TrigT1Calo/TriggerTower_ClassDEF.h"

//For TrigT1Keys
#include "TrigT1Calo/JetElementKey.h"
#include "TrigT1Calo/TriggerTowerKey.h"

/*---------------------------------------------------------*/
TrigT1CaloMonTool::TrigT1CaloMonTool(const std::string & type, 
				 const std::string & name,
				 const IInterface* parent)
  : MonitorToolBase(type, name, parent)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 
  declareProperty("DataVector<TriggerTower>",  m_TriggerTowerContainerName = "LVL1TriggerTowers");
  declareProperty("JetElementContainer",  m_JetElementContainerName = "LVL1JetElements");
  
  //ROOT File directory
  declareProperty("histoPathBase",m_path = "/" );
//"/disk/f8a/home/eew/ATLAS/athena/analysis/TrigT1CaloMonitoring/run");
  /*
  declareProperty("energyThreshold",m_Threshold=50.); //Threshold in MeV
  */
  declareProperty("towersContainerName",m_towersContName="CombinedTower"); //SG Tower Container
}

/*---------------------------------------------------------*/
TrigT1CaloMonTool::~TrigT1CaloMonTool()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloMonTool:: initialize()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  
  StatusCode sc;

  m_stem=m_THistSvc_streamname+m_path;

  log << MSG::INFO << "ETHAN in initialize()" << endreq;
 

  sc = service( "StoreGateSvc", m_StoreGate);// m_eventStore or m_storeGate?
  if( sc.isFailure() ) {
    log << MSG::ERROR << name() << ": Unable to locate Service StoreGateSvc" << endreq;
    return sc;
  }

  SetBookStatus(false);

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloMonTool::bookHists()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  
  log << MSG::INFO << "ETHAN in bookHists()" << endreq;
  //  log << MSG::INFO << "Using base path " << m_stem << endreq;
    
  //  m_TilenTowers = book1D("tilenTowers","TileCal Tower number",100,0., 100.);
  //  m_TrigT1EtaPhi = book2D("tileTowerEtaPhi","TileCal Position of most en. Tower",80,-4., 4.,64,-3.15,3.15);
  //  m_TrigT1EtaPhi->GetXaxis()->SetTitle("#eta");
  //  m_TrigT1EtaPhi->GetYaxis()->SetTitle("#phi");

  // Bin TT & JE Histos

  m_h_TT_Em_Et = book1D("TT_EM_Et","TT EM Et monitor",255,0,255);
  m_h_TT_Had_Et = book1D("TT_HAD_Et","TT HAD Et monitor",255,0,255);
  m_h_TT_eta = book1D("TT_eta","Trigger Tower eta",100,-5,5);
  //otherwise consider plotting eta in calo regions .
  m_h_TT_phi = book1D("TT_phi","Trigger Tower phi ",64,0,2*M_PI);

  m_h_JE_Em_Et = book1D("JE_EM_Et","JE EM Et monitor",255,0,255);
  m_h_JE_Had_Et = book1D("JE_HAD_Et","JE HAD Et monitor",255,0,255);
  m_h_JE_eta = book1D("JE_eta","JE eta",100,-5,5);
  m_h_JE_phi = book1D("JE_phi","JE phi ",64,0,2*M_PI);

  m_h_TT_Tot_Et = book1D("TT_Tot_Et","TT EM+HAD Et monitor",10,0,10);
  m_h_JE_Tot_Et = book1D("JE_Tot_Et","JE EM+HAD Et monitor",10,0,10);

  m_h_TT_Em10_Et = book1D("TT_EM10_Et","TT EM Et magnified",11,-1,10);
  m_h_TT_Had10_Et = book1D("TT_HAD10_Et","TT HAD Et magnified",11,-1,10);
  m_h_JE_Em10_Et = book1D("JE_EM10_Et","JE EM Et magnified",11,-1,10); 
  m_h_JE_Had10_Et = book1D("JE_HAD10_Et","JE HAD Et magnified",11,-1,10); 

  m_h_TT_key = book1D("TT_Key","TT Key",100,0,8000);

  //TT & JE Histos in Calo Regions
    
  m_h_Barrel_TT_phi = book1D("Barrel_TT_phi","Barrel_TT phi",64,0,2*M_PI);
  m_h_Barrel_TT_Em_Et = book1D("Barrel_TT_EM_Et","Barrel_TT EM Et",255,0,255);
  m_h_Barrel_TT_Had_Et = book1D("Barrel_TT_HAD_Et","Barrel_TT HAD Et",255,0,255);
  m_h_Barrel_JE_phi = book1D("Barrel_JE_phi","Barrel_JE phi",64,0,2*M_PI);
  m_h_Barrel_JE_Em_Et = book1D("Barrel_JE_EM_Et","Barrel_JE EM Et",255,0,255);
  m_h_Barrel_JE_Had_Et = book1D("Barrel_JE_HAD_Et","Barrel_JE HAD Et",255,0,255);

  //note EndCap == HEC 
  //the 10 in the identifier denotes looking at the low energy, 0-10 GeV, range 
  m_h_EC10_TT_Em_Et = book1D("EC10_TT_EM_Et","EndCap_TT EM Et",10,0,10);
  m_h_EC10_TT_Had_Et = book1D("EC10_TT_HAD_Et","EndCap_TT HAD Et",10,0,10);
  m_h_EC10_JE_Em_Et = book1D("EC10_JE_EM_Et","EndCap_JE EM Et",10,0,10);
  m_h_EC10_JE_Had_Et = book1D("EC10_JE_HAD_Et","EndCap_JE HAD Et",10,0,10);

  m_h_Barrel10_TT_Em_Et = book1D("Barrel10_TT_EM_Et","Barrel10_TT EM Et",10,0,10);
  m_h_Barrel10_TT_Had_Et = book1D("Barrel10_TT_HAD_Et","Barrel10_TT HAD Et",10,0,10);
  m_h_Barrel10_JE_Em_Et = book1D("Barrel10_JE_EM_Et","Barrel10_JE EM Et",10,0,10);
  m_h_Barrel10_JE_Had_Et = book1D("Barrel10_JE_HAD_Et","Barrel10_JE HAD Et",10,0,10);

  m_h_FCAL10_TT_Em_Et = book1D("FCAL10_TT_EM_Et","FCAL10_TT EM Et",10,0,10);
  m_h_FCAL10_TT_Had_Et = book1D("FCAL10_TT_HAD_Et","FCAL10_TT HAD Et",10,0,10);
  m_h_FCAL10_JE_Em_Et = book1D("FCAL10_JE_EM_Et","FCAL10_JE EM Et",10,0,10);
  m_h_FCAL10_JE_Had_Et = book1D("FCAL10_JE_HAD_Et","FCAL10_JE HAD Et",10,0,10);

  m_h_EC_TT_phi = book1D("EC_TT_phi","EC_TT phi",64,0,2*M_PI);
  m_h_EC_TT_Em_Et = book1D("EC_TT_EM_Et","EC_TT EM Et",255,0,255);
  m_h_EC_TT_Had_Et = book1D("EC_TT_HAD_Et","EC_TT HAD Et",255,0,255); 
  m_h_EC_JE_phi = book1D("EC_JE_phi","EC_JE phi",64,0,2*M_PI);
  m_h_EC_JE_Em_Et = book1D("EC_JE_EM_Et","EC_JE EM Et",255,0,255);
  m_h_EC_JE_Had_Et = book1D("EC_JE_HAD_Et","EC_JE HAD Et",255,0,255);

  m_h_FCAL_TT_phi = book1D("FCAL_TT_phi","FCAL_TT phi",64,0,2*M_PI); 
  m_h_FCAL_TT_Em_Et = book1D("FCAL_TT_EM_Et","FCAL_TT EM Et",255,0,255);
  m_h_FCAL_TT_Had_Et = book1D("FCAL_TT_HAD_Et","FCAL_TT HAD Et",255,0,255);
  m_h_FCAL_JE_phi = book1D("FCAL_JE_phi","FCAL_JE phi",64,0,2*M_PI); 
  m_h_FCAL_JE_Em_Et = book1D("FCAL_JE_EM_Et","FCAL_JE EM Et",255,0,255);
  m_h_FCAL_JE_Had_Et = book1D("FCAL_JE_HAD_Et","FCAL_JE HAD Et",255,0,255);

  m_h_TT_etaphi = book2D("Eta_Phi_TT","TT Eta_Phi",100,-5,5, 64,0,2*M_PI);
  m_h_TT_etaphi_hitmap = book2D("Eta_Phi_HM_TT","TT Eta_Phi_HitMap",40,-5,5, 64,0,2*M_PI);

  m_h_JE_etaphi = book2D("Eta_Phi_JE","JE Eta_Phi",100,-5,5, 64,0,2*M_PI);
  m_h_JE_etaphi_hitmap = book2D("Eta_Phi_HM_JE","JE Eta_Phi_HitMap",100,-5,5, 64,0,2*M_PI);

  //JE & TT plots, Et>10GeV

  m_h_TT_eta_gt10 = book1D("TT_eta_gt10","10 Gev TT eta",100,-5,5);
  m_h_TT_phi_gt10 = book1D("TT_phi_gt10","10GeV TT phi",64,0,2*M_PI); 
  m_h_TT_EC_phi_gt10 = book1D("TT_EC_phi_gt10","10GeV TT EC phi",64,0,2*M_PI); 
  m_h_TT_Barrel_phi_gt10 = book1D("TT_Barrel_phi_gt10","10GeV TT Barrel phi",64,0,2*M_PI); 
  m_h_TT_FCAL_phi_gt10 = book1D("TT_FCAL_phi_gt10","10GeV TT FCAL phi",64,0,2*M_PI); 

  m_h_JE_eta_gt10 = book1D("JE_eta_gt10","10 Gev JE eta",100,-5,5);
  m_h_JE_phi_gt10 = book1D("JE_phi_gt10","10GeV JE phi",64,0,2*M_PI); 
  m_h_JE_EC_phi_gt10 = book1D("JE_EC_phi_gt10","10GeV JE EC phi",64,0,2*M_PI); 
  m_h_JE_Barrel_phi_gt10 = book1D("JE_Barrel_phi_gt10","10GeV JE Barrel phi",64,0,2*M_PI); 
  m_h_JE_FCAL_phi_gt10 = book1D("JE_FCAL_phi_gt10","10GeV JE FCAL phi",64,0,2*M_PI); 

  //Plots of the /Calorimeter information (combined LAr and Tile info)
  m_h_Calo_phi = book1D("Calo_phi","Calo phi",64,-M_PI,M_PI); 
  m_h_Calo_eta = book1D("Calo_eta","Calo eta",100,-5,5); 
  m_h_Calo_phi_gt10 = book1D("Calo_phi_gt10","Calo phi gt 10",64,-M_PI,M_PI); 
  m_h_Calo_eta_gt10 = book1D("Calo_eta_gt10","Calo eta gt 10",100,-5,5); 

  m_h_Calo_Et10 = book1D("Calo_Et10","Calo Et",100,0,10.*GeV); //Note, need to specify GeV 
  m_h_Calo_Et = book1D("Calo_Et","Calo Et",255,0,255.*GeV); 

  m_h_Calo_etaphi = book2D("Calo_Eta_Phi","Calo Eta_Phi",100,-5,5, 64, -M_PI, M_PI);
  m_h_Calo_etaphi_hitmap = book2D("Calo_Eta_Phi_HM","Calo Eta_Phi_HitMap",100,-5,5,64, -M_PI, M_PI);

  m_h_Calo_key = book1D("Calo_Key","Calo Key",100,0,8000);

  //Calo & TT ratio plots
  m_h_TT_Calo_eta = book1D("TT_Calo_eta","eta TT vs Calo",100, -5, 5);
  m_h_TT_Calo_Et = book1D("TT_Calo_Et","Et TT vs Calo", 255, 0, 255);

  m_h_TT_Calo_EtTower = book2D("TT_Calo_EtTower","Calo and TT Tower Et",255,0,255.*GeV, 255, 0, 255);
  m_h_TT_Calo_EtaTower = book2D("TT_Calo_EtaTower","Calo and TT Tower Eta",100,-5,5, 100, -5, 5);
  m_h_TT_Calo_PhiTower = book2D("TT_Calo_PhiTower","Calo and TT Tower Phi",64,0,2*M_PI,64,0,2*M_PI);

  //Discrepancy check plots

  m_h_TT_D_under_phi = book1D("TT_D_under_phi","TT phi Discrepant",64,0,2*M_PI);
  m_h_TT_D_under_eta = book1D("TT_D_under_eta","TT eta Discrepant",100, -5, 5);
  m_h_Calo_D_under_phi = book1D("Calo_D_under_phi","Calo phi Discrepant",64,-M_PI,M_PI);
  m_h_Calo_D_under_eta = book1D("Calo_D_under_eta","Calo eta Discrepant",100, -5, 5);

  m_h_TT_D_over_phi = book1D("TT_D_over_phi","TT phi Discrepant",64,0,2*M_PI);
  m_h_TT_D_over_eta = book1D("TT_D_over_eta","TT eta Discrepant",100, -5, 5);
  m_h_Calo_D_over_phi = book1D("Calo_D_over_phi","Calo phi Discrepant",64,-M_PI,M_PI);
  m_h_Calo_D_over_eta = book1D("Calo_D_over_eta","Calo eta Discrepant",100, -5, 5);

  m_h_TT_D_over_Et = book1D("TT_D_over_Et","Et TT Discrepant", 255, 0, 255);
  m_h_TT_D_under_Et = book1D("TT_D_under_Et","Et TT Discrepant", 255, 0, 255);
  m_h_TT_D_over_EmEt = book1D("TT_D_over_EmEt","Em Et TT Discrepant", 255, 0, 255);
  m_h_TT_D_under_EmEt = book1D("TT_D_under_EmEt","Em Et TT Discrepant", 255, 0, 255);
  m_h_TT_D_over_HadEt = book1D("TT_D_over_HadEt","Had Et TT Discrepant", 255, 0, 255);
  m_h_TT_D_under_HadEt = book1D("TT_D_under_HadEt","Had Et TT Discrepant", 255, 0, 255);

  m_h_Calo_D_over_Et = book1D("Calo_D_over_Et","Et Calo Discrepant", 255, 0, 255);
  m_h_Calo_D_under_Et = book1D("Calo_D_under_Et","Et Calo Discrepant", 255, 0, 255);


  SetBookStatus(true);
  
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloMonTool::fillHists()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  
  log << MSG::DEBUG << "Ethan in fillHists()" << endreq;

 //Retrieve Tower collection from SG
  //  StatusCode sc = m_eventStore->retrieve(tower_container, m_towersContName);
  /*
  StatusCode sc = service("StoreGateSvc", m_StoreGate);
  if (sc.isFailure()) {
    log << MSG::ERROR<< "Unable to retrieve pointer to StoreGateSvc"<< endreq;
  }
  */

  //Retrieve TriggerTowers from SG
  const TriggerTowerCollection* TriggerTowerTES = 0; 
  StatusCode sc = m_StoreGate->retrieve(TriggerTowerTES, "LVL1TriggerTowers"); 
  if( sc.isFailure()  ||  !TriggerTowerTES ) {
    log << MSG::ERROR<< "No ESD LVL1::TriggerTowerCollection container found"<< endreq; 
  }
  
  //Retrieve JetElements from SG
  const JetElementCollection* JetElementTES = 0;
  sc=m_StoreGate->retrieve( JetElementTES, m_JetElementContainerName);
  if( sc.isFailure()  ||  !JetElementTES ) {
    log << MSG::ERROR<< "No ESD LVL1::JetElementCollection container found"<< endreq; 
  }
  
  //Retrieve Tower collection from SG
  const CaloTowerContainer* CaloTowerTES = 0; 
  sc=m_StoreGate->retrieve(CaloTowerTES, "CombinedTower"); 
  if( sc.isFailure()  ||  !CaloTowerTES ) {
    log << MSG::ERROR<< "No CaloTowerContainer found"<< endreq; 
  }

  //Key Stuff for Trigger vs Calo Tower comparison
  TriggerTowerKey TTKey(0,0);
  JetElementKey* JEKey; 
  std::map<int, TriggerTower *>* m_ttContainer;
  m_ttContainer       =new std::map<int, TriggerTower*>;          //Create a map to hold the towers
  int key = 0;

  //TriggerTowers:

  TriggerTowerCollection::const_iterator TriggerTowerIterator    = TriggerTowerTES->begin(); 
  TriggerTowerCollection::const_iterator TriggerTowerIteratorEnd = TriggerTowerTES->end(); 

  //  int counter = 0;
  std::cout<<"before Trigger Loop"<<std::endl;
  log << MSG::ERROR << "before Trigger Loop"<<endreq;

  for (; TriggerTowerIterator != TriggerTowerIteratorEnd; ++TriggerTowerIterator) 
  {

    key = TTKey.ttKey((*TriggerTowerIterator)->phi(),(*TriggerTowerIterator)->eta());
    //fill m_ttcontainer
    std::map<int, TriggerTower *>::iterator test=m_ttContainer->find( key );
    if (test != m_ttContainer->end()){
         log << MSG::ERROR
           <<"key: "<<key<< "already used: TT exists!!!!!! (shouldn't happen - beat Ed up)"
           <<endreq ;
       }
       m_ttContainer->insert(std::map<int, TriggerTower*>::value_type(key,*TriggerTowerIterator));

    m_h_TT_Em_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_TT_Had_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 
    m_h_TT_eta->Fill( (*TriggerTowerIterator)->eta(), 1.); 
    m_h_TT_phi->Fill( (*TriggerTowerIterator)->phi(), 1.); 

    m_h_TT_etaphi_hitmap->Fill( (*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi(), 1.);
    m_h_TT_etaphi->Fill( (*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi(), (*TriggerTowerIterator)->emEnergy()+(*TriggerTowerIterator)->hadEnergy());

    //Fill Combined Energy plot
    m_h_TT_Tot_Et->Fill( (*TriggerTowerIterator)->hadEnergy()+(*TriggerTowerIterator)->emEnergy(), 1.); 
    //    m_h_TT_Tot_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 

    m_h_TT_Em10_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_TT_Had10_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 

    //   m_h_TT_key->Fill((*TriggerTowerIterator)->key(), 1.);

    double eta = (*TriggerTowerIterator)->eta();

    eta*=eta;
    eta=sqrt(eta);

    if (((*TriggerTowerIterator)->emEnergy()+(*TriggerTowerIterator)->hadEnergy())>10){
      m_h_TT_key->Fill((*TriggerTowerIterator)->key(), 1.);

      m_h_TT_eta_gt10->Fill( (*TriggerTowerIterator)->eta(),1.);      
      m_h_TT_phi_gt10->Fill( (*TriggerTowerIterator)->phi(), 1.); 
      if(1.5 < eta && eta< 3.2)m_h_TT_EC_phi_gt10->Fill( (*TriggerTowerIterator)->phi(), 1.); 
      if(0.1 < eta && eta < 1.5)m_h_TT_Barrel_phi_gt10->Fill( (*TriggerTowerIterator)->phi(), 1.); 
      if(3.2 < eta && eta< 4.9)m_h_TT_FCAL_phi_gt10->Fill( (*TriggerTowerIterator)->phi(), 1.); 
    }

    if(0.1 < eta && eta < 1.5){//Barrel 
    m_h_Barrel10_TT_Em_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_Barrel10_TT_Had_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 

    m_h_Barrel_TT_phi->Fill( (*TriggerTowerIterator)->phi(), 1.); 
    m_h_Barrel_TT_Em_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_Barrel_TT_Had_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 
    }

    if(1.5 < eta && eta< 3.2){//EC/HEC
    m_h_EC10_TT_Em_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_EC10_TT_Had_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 

    m_h_EC_TT_phi->Fill( (*TriggerTowerIterator)->phi(), 1.); 
    m_h_EC_TT_Em_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_EC_TT_Had_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 
    }

    if(3.2 < eta && eta< 4.9){//FCAL
    m_h_FCAL10_TT_Em_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_FCAL10_TT_Had_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 

    m_h_FCAL_TT_phi->Fill( (*TriggerTowerIterator)->phi(), 1.); 
    m_h_FCAL_TT_Em_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_FCAL_TT_Had_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 
    }
  }
  
  
  //JetElements

  JetElementCollection::const_iterator JetElementIterator    = JetElementTES->begin();
  JetElementCollection::const_iterator JetElementIteratorEnd = JetElementTES->end();

  for (; JetElementIterator != JetElementIteratorEnd; ++JetElementIterator)
  {
    m_h_JE_Em_Et->Fill( (*JetElementIterator)->emEnergy(), 1.);
    m_h_JE_Had_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.);
    m_h_JE_eta->Fill( (*JetElementIterator)->eta(), 1.);
    m_h_JE_phi->Fill( (*JetElementIterator)->phi(), 1.);

    //Fill combined energy
    m_h_JE_Tot_Et->Fill( (*JetElementIterator)->hadEnergy()/*+(*JetElementIterator)->emEnergy()*/, 1.);
    m_h_JE_Tot_Et->Fill( (*JetElementIterator)->emEnergy(), 1.);

    m_h_JE_Em10_Et->Fill( (*JetElementIterator)->emEnergy(), 1.); 
    m_h_JE_Had10_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.); 

    m_h_JE_etaphi->Fill( (*JetElementIterator)->eta(), (*JetElementIterator)->phi(),(*JetElementIterator)->emEnergy()+(*JetElementIterator)->hadEnergy());
    m_h_JE_etaphi_hitmap->Fill( (*JetElementIterator)->eta(), (*JetElementIterator)->phi(), 1.);

    double eta = (*JetElementIterator)->eta();

    eta*=eta;
    eta=sqrt(eta);

    if (((*JetElementIterator)->emEnergy()+(*JetElementIterator)->hadEnergy())>10){
      m_h_JE_eta_gt10->Fill( (*JetElementIterator)->eta(),1.);      
      m_h_JE_phi_gt10->Fill( (*JetElementIterator)->phi(), 1.); 
      if(1.5 < eta && eta< 3.2)m_h_JE_EC_phi_gt10->Fill( (*JetElementIterator)->phi(), 1.); 
      if(0.1 < eta && eta < 1.5)m_h_JE_Barrel_phi_gt10->Fill( (*JetElementIterator)->phi(), 1.); 
      if(3.2 < eta && eta< 4.9)m_h_JE_FCAL_phi_gt10->Fill( (*JetElementIterator)->phi(), 1.); 
    }

    if(0.1 < eta && eta< 1.5){//Barrel
      m_h_Barrel10_JE_Em_Et->Fill( (*JetElementIterator)->emEnergy(), 1.);
      m_h_Barrel10_JE_Had_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.);
      
      m_h_Barrel_JE_phi->Fill( (*JetElementIterator)->phi(), 1.);
      m_h_Barrel_JE_Em_Et->Fill( (*JetElementIterator)->emEnergy(), 1.);
      m_h_Barrel_JE_Had_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.);
    }

    if(1.5 < eta && eta< 3.2){//EC/EC
      m_h_EC10_JE_Em_Et->Fill( (*JetElementIterator)->emEnergy(), 1.);
      m_h_EC10_JE_Had_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.);
      
      m_h_EC_JE_phi->Fill( (*JetElementIterator)->phi(), 1.);
      m_h_EC_JE_Em_Et->Fill( (*JetElementIterator)->emEnergy(), 1.);
      m_h_EC_JE_Had_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.);
    }
    
    if(3.2 < eta && eta< 4.9){//FCAL
      m_h_FCAL10_JE_Em_Et->Fill( (*JetElementIterator)->emEnergy(), 1.);
      m_h_FCAL10_JE_Had_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.);
      
      m_h_FCAL_JE_phi->Fill( (*JetElementIterator)->phi(), 1.);
      m_h_FCAL_JE_Em_Et->Fill( (*JetElementIterator)->emEnergy(), 1.);
      m_h_FCAL_JE_Had_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.);
    }
   
  }

  //CaloTowers (combined LAr and Tile data)

  CaloTowerContainer::const_iterator CaloTowerIterator    = CaloTowerTES->begin();
  CaloTowerContainer::const_iterator CaloTowerIteratorEnd = CaloTowerTES->end();

  int TTkey; //change this to key and use previous int
  double TTtoCaloEnergyRatio;

 for (; CaloTowerIterator != CaloTowerIteratorEnd; ++CaloTowerIterator)
  {
    if((*CaloTowerIterator)->et()>0.*GeV){
      m_h_Calo_phi->Fill( (*CaloTowerIterator)->phi(), 1.);
      m_h_Calo_eta->Fill( (*CaloTowerIterator)->eta(), 1.);
      //need to convert phi to 0->2pi range:
   
      if((*CaloTowerIterator)->et()>10.*GeV){
	double calophi = (*CaloTowerIterator)->phi();
	if(calophi>M_PI)calophi-=2*M_PI;
	TTkey = TTKey.ttKey(calophi,(*CaloTowerIterator)->eta());
	m_h_Calo_key->Fill(TTkey,1);
	std::map<int, TriggerTower *>::iterator test=m_ttContainer->find( TTkey ); 

	//	TriggerTowerIterator->find()
	//	m_h_TT_Calo_EtTower->Fill(1,1,1);

	if (test != m_ttContainer->end()){
	  //	   double energyOfNeighbour=(test->second)->emEnergy();//energy();
	   m_h_TT_Calo_EtTower->Fill((*CaloTowerIterator)->et(),(test->second)->emEnergy()+(test->second)->hadEnergy(),1);
	   m_h_TT_Calo_PhiTower->Fill(calophi, (test->second)->phi(),1);
	   m_h_TT_Calo_EtaTower->Fill((*CaloTowerIterator)->eta(),(test->second)->eta(),1);
	   //log << MSG::ERROR << "TT of Key "<<TTkey<<" has emNRG "<<(test->second)->emEnergy()<<endreq;


	   //Find ratio of energy TT to Calo, and if off by 10% or more plot the eta phi distributions of those events
	   TTtoCaloEnergyRatio= ((test->second)->emEnergy()+(test->second)->hadEnergy())/((*CaloTowerIterator)->et());

	   if(TTtoCaloEnergyRatio >1.1){
	     m_h_TT_D_over_phi->Fill((test->second)->phi(),1); 
	     m_h_TT_D_over_eta->Fill((test->second)->eta(),1);
	     m_h_Calo_D_over_phi->Fill((*CaloTowerIterator)->phi(),1);
	     m_h_Calo_D_over_eta->Fill((*CaloTowerIterator)->eta(),1);

	     m_h_TT_D_over_Et->Fill((test->second)->emEnergy()+(test->second)->hadEnergy(),1);
	     m_h_TT_D_over_EmEt->Fill((test->second)->emEnergy(),1);
	     m_h_TT_D_over_HadEt->Fill((test->second)->hadEnergy(),1);

	     m_h_Calo_D_over_Et->Fill((*CaloTowerIterator)->et(),1);


	   }
	   if(TTtoCaloEnergyRatio <0.9){
	     m_h_TT_D_under_phi->Fill((test->second)->phi(),1); 
	     m_h_TT_D_under_eta->Fill((test->second)->eta(),1);
	     m_h_Calo_D_under_phi->Fill((*CaloTowerIterator)->phi(),1);
	     m_h_Calo_D_under_eta->Fill((*CaloTowerIterator)->eta(),1);

	     m_h_TT_D_under_Et->Fill((test->second)->emEnergy()+(test->second)->hadEnergy(),1);
	     m_h_TT_D_under_EmEt->Fill((test->second)->emEnergy(),1);
	     m_h_TT_D_under_HadEt->Fill((test->second)->hadEnergy(),1);

	     m_h_Calo_D_under_Et->Fill((*CaloTowerIterator)->et(),1);
	   }

	}else{
	   m_h_TT_Calo_EtTower->Fill((*CaloTowerIterator)->et(),0.0,1.0);
	   m_h_TT_Calo_PhiTower->Fill(calophi,0.,1.);
	   m_h_TT_Calo_EtaTower->Fill((*CaloTowerIterator)->eta(),-5.,1.);
	   log << MSG::ERROR << "TT of Key "<<TTkey<<" has emNRG "<<0<<endreq;
	}

	m_h_Calo_phi_gt10->Fill( (*CaloTowerIterator)->phi(), 1.);
	m_h_Calo_eta_gt10->Fill( (*CaloTowerIterator)->eta(), 1.);
      
      m_h_Calo_etaphi_hitmap->Fill( (*CaloTowerIterator)->eta(), (*CaloTowerIterator)->phi(), 1.);
      m_h_Calo_etaphi->Fill( (*CaloTowerIterator)->eta(), (*CaloTowerIterator)->phi(), (*CaloTowerIterator)->et());

      m_h_Calo_Et->Fill( (*CaloTowerIterator)->et(), 1.);
      }
    }

    m_h_Calo_Et10->Fill( (*CaloTowerIterator)->et(), 1.);
  }

  /*
  // Pointer to a Tile cell container
  const CaloTowerContainer* tower_container;
  //Retrieve Tower collection from SG
  StatusCode sc = m_eventStore->retrieve(tower_container, m_towersContName);
  
  if (sc.isFailure()){
    log << MSG::DEBUG << "TrigT1CaloMonTool: Retrieval of Tile towers from container " << m_towersContName << " failed" << endreq;
    return sc;
  } else {
    log << MSG::VERBOSE << "TrigT1CaloMonTool: Retrieval of Tile towers from container " << m_towersContName << " succeeded" << endreq;
  }

  CaloTowerContainer::const_iterator iTower = tower_container->begin();
  CaloTowerContainer::const_iterator lastTower  = tower_container->end();
  for( ; iTower != lastTower; ++iTower) {

    if ( (*iTower)->energy() > energy_most ) {
      const CaloTower* tower_ptr = *iTower;    
      energy_most = tower_ptr->energy();
      et_most     = tower_ptr->et();
      eta_most    = tower_ptr->eta();
      phi_most    = tower_ptr->phi();
      ncells_most = tower_ptr->getNumberOfCells();  
      set_most =true;
    }

  }
  
  iTower = tower_container->begin();
  for( ; iTower != lastTower; ++iTower) {

    const CaloTower* tower_ptr = *iTower;    // pointer to tower object
    double energy = tower_ptr->energy();
    double phi    = tower_ptr->phi();
    double eta    = tower_ptr->eta();

m_TileTowerEta->Fill(eta_most,1.);


  }
  */

  return StatusCode::SUCCESS;

}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloMonTool::finalHists()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "Ethan in finalHists()" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloMonTool::checkHists(bool /* fromFinalize */)
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "Ethan in checkHists()" << endreq;

  return StatusCode::SUCCESS;
}

TH1D* TrigT1CaloMonTool::book1D(std::string nam, std::string tit, int nx, double xmin, double xmax){
  TH1D *hist = new TH1D(TString(nam), TString(tit), nx, xmin, xmax);
  
  if(ToolRootHistSvc()->regHist(m_stem+"/"+nam, hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< m_stem+nam << endreq;
  }
  
  return hist;
}

TH2D* TrigT1CaloMonTool::book2D(std::string nam, std::string tit, int nx, double xmin, double xmax,  
	     int ny, double ymin, double ymax)
{		
  TH2D *hist = new TH2D(TString(nam), TString(tit), nx, xmin, xmax, ny, ymin, ymax);
  
  if(ToolRootHistSvc()->regHist(m_stem+"/"+nam, hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< m_stem+nam << endreq;
  }
  
  return hist;
}

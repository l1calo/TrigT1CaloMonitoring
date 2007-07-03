// ********************************************************************
//
// NAME:     TrigT1CaloMonTool.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Ethan-Etienne Woehrling (eew@hep.ph.bham.ac.uk)
//           Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
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

#include "CaloEvent/CaloCell.h"
#include "CaloEvent/CaloCellContainer.h"

#include "TrigT1CaloMonitoring/TrigT1CaloBSMonTool.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "CLHEP/Units/SystemOfUnits.h"

//to open TT & JE Persistent objects
#include "TrigT1Calo/TriggerTowerCollection.h"
#include "TrigT1Calo/JetElementCollection.h"
#include "TrigT1Calo/TrigT1CaloDict.h"
#include "TrigT1Calo/TriggerTower_ClassDEF.h"

//For TrigT1Keys
#include "TrigT1Calo/JetElementKey.h"
#include "TrigT1Calo/TriggerTowerKey.h"

//to form my own TT from calo cells
#include "TrigT1Calo/InternalTriggerTower.h"
#include "CaloEvent/CaloSampling.h"

/*---------------------------------------------------------*/
TrigT1CaloBSMonTool::TrigT1CaloBSMonTool(const std::string & type, const std::string & name,
					 const IInterface* parent)
  : ManagedMonitorToolBase ( type, name, parent )
/*---------------------------------------------------------*/
{
  //declareInterface<IMonitorToolBase>(this); 
  declareProperty("BS_TriggerTowerContainer",  m_TriggerTowerContainerName = "LVL1TriggerTowers");
  declareProperty("BS_JetElementContainer",  m_JetElementContainerName = "LVL1JetElements");
  
  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/CMM") ;
  declareProperty( "DataType", m_DataType="") ;

  declareProperty("towersContainerName",m_towersContName="CombinedTower"); //SG Tower Container
}

/*---------------------------------------------------------*/
TrigT1CaloBSMonTool::~TrigT1CaloBSMonTool()
/*---------------------------------------------------------*/
{
}


/*---------------------------------------------------------*/
StatusCode TrigT1CaloBSMonTool::bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "in TrigT1CaloBSMonTool::bookHistograms" << endreq;

  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR
	   << "Unable to retrieve pointer to StoreGateSvc"
	   << endreq;
      return sc;
    }
 

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  MonGroup Calo_TriggerTower ( this, (m_PathInRootFile+"/TriggerTower").c_str(), expert, eventsBlock );
  HistoBooker* TriggerTower_Booker = new HistoBooker(&Calo_TriggerTower, &log, m_DataType);

  MonGroup TT_ChannelEnergy( this, (m_PathInRootFile+"/ChannelEnergy").c_str(), expert, eventsBlock );
  HistoBooker* ChannelEnergy_Booker = new HistoBooker(&TT_ChannelEnergy, &log, m_DataType);

  MonGroup Calo_JetElements ( this, (m_PathInRootFile+"/JetElements").c_str(), expert, eventsBlock );
  HistoBooker* JetElements_Booker = new HistoBooker(&Calo_JetElements, &log, m_DataType);

  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
   //Energy distribution per Channel
      TH1F* help1;
      std::vector <TH1F*> help2;
      std::string name,title;
      std::stringstream buffer1, buffer2, etabuffer, phibuffer;

      double eta,phi;

   for (int i=0;i<4;i++)
     {
       buffer1.str("");
       //buffer1.clear();
       buffer1<<i;

       help2.clear();

       eta=-4.9+i*0.425;
       etabuffer.str("");
       etabuffer<<eta;
       log << MSG::DEBUG << "eta " << etabuffer.str()<< endreq;

       for (int j=0;j<16;j++)
	 {
	   buffer2.str("");
	   //buffer2.clear();
	   buffer2<<j;

	   phi=j*2*M_PI/16;
	   phibuffer.str("");
	   phibuffer<<phi;

	   name = "TTChannelNo_" + buffer1.str() + "_" + buffer2.str();
	   title = "TT Energy distribution: eta=" + etabuffer.str() + ", phi=" + phibuffer.str();
	   
	   log << MSG::DEBUG << "title" <<title<< endreq;

	   help1=ChannelEnergy_Booker->book1F(name,title,255,0,255,"Et [GeV]");
	   help2.push_back(help1);
	 }
       m_h_TT_channels.push_back(help2);
     }

  // Bin TT & JE Histos
  m_h_TT_Em_Et = TriggerTower_Booker->book1F("TT_EM_Et","TT EM Et monitor",255,0,255,"Et [GeV]");
  m_h_TT_Had_Et = TriggerTower_Booker->book1F("TT_HAD_Et","TT HAD Et monitor",255,0,255,"Et [GeV]");

  m_h_TT_eta = TriggerTower_Booker->book1F("TT_eta","Trigger Tower eta",100,-5,5,"#eta");
  m_h_TT_phi = TriggerTower_Booker->book1F("TT_phi","Trigger Tower phi ",64,0,2*M_PI,"#phi");
  
  m_h_TT_Em10_Et = TriggerTower_Booker->book1F("TT_EM10_Et","TT EM Et magnified",11,-1,10,"Et [GeV]");
  m_h_TT_Had10_Et = TriggerTower_Booker->book1F("TT_HAD10_Et","TT HAD Et magnified",11,-1,10,"Et [GeV]");
  m_h_TT_Tot_Et = TriggerTower_Booker->book1F("TT_Tot_Et","TT EM+HAD Et magnified",10,0,10,"Et [GeV]");

  //m_h_TT_key = TriggerTower_Booker->book1F("TT_Key","TT Key",100,0,8000,"Key");

  m_h_Barrel_TT_phi = TriggerTower_Booker->book1F("Barrel_TT_phi","Barrel_TT phi",64,0,2*M_PI, "#phi");
  m_h_Barrel_TT_Em_Et = TriggerTower_Booker->book1F("Barrel_TT_EM_Et","Barrel_TT EM Et",255,0,255,"Et [GeV]");
  m_h_Barrel_TT_Had_Et = TriggerTower_Booker->book1F("Barrel_TT_HAD_Et","Barrel_TT HAD Et",255,0,255,"Et [GeV]");
  m_h_Barrel10_TT_Em_Et = TriggerTower_Booker->book1F("Barrel10_TT_EM_Et","Barrel10_TT EM Et",10,0,10,"Et [GeV]");
  m_h_Barrel10_TT_Had_Et = TriggerTower_Booker->book1F("Barrel10_TT_HAD_Et","Barrel10_TT HAD Et",10,0,10,"Et [GeV]");

  m_h_EC_TT_phi = TriggerTower_Booker->book1F("EC_TT_phi","EC_TT phi",64,0,2*M_PI, "#phi");
  m_h_EC_TT_Em_Et = TriggerTower_Booker->book1F("EC_TT_EM_Et","EC_TT EM Et",255,0,255,"Et [GeV]");
  m_h_EC_TT_Had_Et = TriggerTower_Booker->book1F("EC_TT_HAD_Et","EC_TT HAD Et",255,0,255,"Et [GeV]"); 
  m_h_EC10_TT_Em_Et = TriggerTower_Booker->book1F("EC10_TT_EM_Et","EndCap_TT EM Et",10,0,10,"Et [GeV]");
  m_h_EC10_TT_Had_Et = TriggerTower_Booker->book1F("EC10_TT_HAD_Et","EndCap_TT HAD Et",10,0,10,"Et [GeV]");

  m_h_FCAL_TT_phi = TriggerTower_Booker->book1F("FCAL_TT_phi","FCAL_TT phi",64,0,2*M_PI, "#phi"); 
  m_h_FCAL_TT_Em_Et = TriggerTower_Booker->book1F("FCAL_TT_EM_Et","FCAL_TT EM Et",255,0,255,"Et [GeV]");
  m_h_FCAL_TT_Had_Et = TriggerTower_Booker->book1F("FCAL_TT_HAD_Et","FCAL_TT HAD Et",255,0,255,"Et [GeV]");
  m_h_FCAL10_TT_Em_Et = TriggerTower_Booker->book1F("FCAL10_TT_EM_Et","FCAL10_TT EM Et",10,0,10,"Et [GeV]");
  m_h_FCAL10_TT_Had_Et = TriggerTower_Booker->book1F("FCAL10_TT_HAD_Et","FCAL10_TT HAD Et",10,0,10,"Et [GeV]");

  m_h_TT_eta_gt10 = TriggerTower_Booker->book1F("TT_eta_gt10","10 Gev TT eta",100,-5,5, "#eta");
  m_h_TT_phi_gt10 = TriggerTower_Booker->book1F("TT_phi_gt10","10GeV TT phi",64,0,2*M_PI, "#phi"); 
  m_h_TT_Barrel_phi_gt10 = TriggerTower_Booker->book1F("TT_Barrel_phi_gt10","10GeV TT Barrel phi",64,0,2*M_PI, "#phi"); 
  m_h_TT_EC_phi_gt10 = TriggerTower_Booker->book1F("TT_EC_phi_gt10","10GeV TT EC phi",64,0,2*M_PI, "#phi"); 
  m_h_TT_FCAL_phi_gt10 = TriggerTower_Booker->book1F("TT_FCAL_phi_gt10","10GeV TT FCAL phi",64,0,2*M_PI, "#phi"); 

  m_h_TT_etaphi = TriggerTower_Booker->book2F("Eta_Phi_TT","TT Eta_Phi",100,-5,5, 64,0,2*M_PI, "#eta", "#phi");
  m_h_TT_etaphi_hitmap = TriggerTower_Booker->book2F("Eta_Phi_HM_TT","TT Eta_Phi_HitMap",40,-5,5, 64,0,2*M_PI, "#eta", "#phi");
  m_h_Calib_TTEM_EtEta = TriggerTower_Booker->book2F("TTEM_EtEta", "TT (EM) Et vs eta", 50, -5, 5, 10, 0, 10, "#eta", "Et [GeV]");
  m_h_Calib_TTHAD_EtEta = TriggerTower_Booker->book2F("TTHAD_EtEta", "TT (HAD) Et vs eta", 50, -5, 5, 10, 0, 10, "#eta", "Et [GeV]");



  m_h_JE_Em_Et = JetElements_Booker->book1F("JE_EM_Et","JE EM Et monitor",255,0,255,"Et [GeV]");
  m_h_JE_Had_Et = JetElements_Booker->book1F("JE_HAD_Et","JE HAD Et monitor",255,0,255,"Et [GeV]");
  m_h_JE_Tot_Et = JetElements_Booker->book1F("JE_Tot_Et","JE EM+HAD Et monitor",10,0,10,"Et [GeV]");

  m_h_JE_eta = JetElements_Booker->book1F("JE_eta","JE eta",100,-5,5,"#eta");
  m_h_JE_phi = JetElements_Booker->book1F("JE_phi","JE phi ",64,0,2*M_PI,"#phi");

  m_h_JE_Tot_Et = JetElements_Booker->book1F("JE_Tot_Et","JE EM+HAD Et monitor",10,0,10,"Et [GeV]");
  m_h_JE_Em10_Et = JetElements_Booker->book1F("JE_EM10_Et","JE EM Et magnified",11,-1,10,"Et [GeV]"); 
  m_h_JE_Had10_Et = JetElements_Booker->book1F("JE_HAD10_Et","JE HAD Et magnified",11,-1,10,"Et [GeV]"); 

  m_h_Barrel_JE_phi = JetElements_Booker->book1F("Barrel_JE_phi","Barrel_JE phi",64,0,2*M_PI, "#phi");
  m_h_Barrel_JE_Em_Et = JetElements_Booker->book1F("Barrel_JE_EM_Et","Barrel_JE EM Et",255,0,255,"Et [GeV]");
  m_h_Barrel_JE_Had_Et = JetElements_Booker->book1F("Barrel_JE_HAD_Et","Barrel_JE HAD Et",255,0,255,"Et [GeV]");

  m_h_EC10_JE_Em_Et = JetElements_Booker->book1F("EC10_JE_EM_Et","EndCap_JE EM Et",10,0,10,"Et [GeV]");
  m_h_EC10_JE_Had_Et = JetElements_Booker->book1F("EC10_JE_HAD_Et","EndCap_JE HAD Et",10,0,10,"Et [GeV]");

  m_h_Barrel10_JE_Em_Et = JetElements_Booker->book1F("Barrel10_JE_EM_Et","Barrel10_JE EM Et",10,0,10,"Et [GeV]");
  m_h_Barrel10_JE_Had_Et = JetElements_Booker->book1F("Barrel10_JE_HAD_Et","Barrel10_JE HAD Et",10,0,10,"Et [GeV]");

  m_h_FCAL10_JE_Em_Et = JetElements_Booker->book1F("FCAL10_JE_EM_Et","FCAL10_JE EM Et",10,0,10,"Et [GeV]");
  m_h_FCAL10_JE_Had_Et = JetElements_Booker->book1F("FCAL10_JE_HAD_Et","FCAL10_JE HAD Et",10,0,10,"Et [GeV]");

  m_h_EC_JE_phi = JetElements_Booker->book1F("EC_JE_phi","EC_JE phi",64,0,2*M_PI, "#phi");
  m_h_EC_JE_Em_Et = JetElements_Booker->book1F("EC_JE_EM_Et","EC_JE EM Et",255,0,255,"Et [GeV]");
  m_h_EC_JE_Had_Et = JetElements_Booker->book1F("EC_JE_HAD_Et","EC_JE HAD Et",255,0,255,"Et [GeV]");

  m_h_FCAL_JE_phi = JetElements_Booker->book1F("FCAL_JE_phi","FCAL_JE phi",64,0,2*M_PI, "#phi"); 
  m_h_FCAL_JE_Em_Et = JetElements_Booker->book1F("FCAL_JE_EM_Et","FCAL_JE EM Et",255,0,255,"Et [GeV]");
  m_h_FCAL_JE_Had_Et = JetElements_Booker->book1F("FCAL_JE_HAD_Et","FCAL_JE HAD Et",255,0,255,"Et [GeV]");

  m_h_JE_etaphi = JetElements_Booker->book2F("Eta_Phi_JE","JE Eta_Phi",100,-5,5, 64,0,2*M_PI, "#eta", "#phi");
  m_h_JE_etaphi_hitmap = JetElements_Booker->book2F("Eta_Phi_HM_JE","JE Eta_Phi_HitMap",100,-5,5, 64,0,2*M_PI, "#eta", "#phi");

  m_h_JE_eta_gt10 = JetElements_Booker->book1F("JE_eta_gt10","10Gev JE eta",100,-5,5, "#eta");
  m_h_JE_phi_gt10 = JetElements_Booker->book1F("JE_phi_gt10","10GeV JE phi",64,0,2*M_PI, "#phi"); 
  m_h_JE_EC_phi_gt10 = JetElements_Booker->book1F("JE_EC_phi_gt10","10GeV JE EC phi",64,0,2*M_PI, "#phi"); 
  m_h_JE_Barrel_phi_gt10 = JetElements_Booker->book1F("JE_Barrel_phi_gt10","10GeV JE Barrel phi",64,0,2*M_PI, "#phi"); 
  m_h_JE_FCAL_phi_gt10 = JetElements_Booker->book1F("JE_FCAL_phi_gt10","10GeV JE FCAL phi",64,0,2*M_PI, "#phi"); 
 
    }

  if( isNewRun ) { }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloBSMonTool::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  
  log << MSG::DEBUG << "in fillHistograms()" << endreq;

  //Retrieve TriggerTowers from SG
  const TriggerTowerCollection* TriggerTowerTES = 0; 
  StatusCode sc = m_storeGate->retrieve(TriggerTowerTES, m_TriggerTowerContainerName); 
  if( (sc==StatusCode::FAILURE) ) 
    {
      log << MSG::DEBUG
	   << "No TriggerTower found in TES at "
	   << m_TriggerTowerContainerName
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  //Retrieve JetElements from SG
  const JetElementCollection* JetElementTES = 0;
  sc=m_storeGate->retrieve( JetElementTES, m_JetElementContainerName);
  if( (sc==StatusCode::FAILURE) ) 
    {
      log << MSG::DEBUG << "No JetElements found in TES at "  << m_JetElementContainerName << endreq ;
      return StatusCode::SUCCESS;
    }

  //Key Stuff for Trigger vs Calo Tower comparison
  TriggerTowerKey TTKey(0,0);
  //JetElementKey* JEKey; 
  std::map<int, TriggerTower *>* m_ttContainer;
  m_ttContainer       =new std::map<int, TriggerTower*>;          //Create a map to hold the towers
  int key = 0;
  

  // =========================================================================================
  // ================================== TriggerTower Plots ===================================
  // =========================================================================================

  // TTs are only filled/stored if they contain energy, 
  // so looping over all TTs will only loop over those with energy

  TriggerTowerCollection::const_iterator TriggerTowerIterator    = TriggerTowerTES->begin(); 
  TriggerTowerCollection::const_iterator TriggerTowerIteratorEnd = TriggerTowerTES->end(); 

  //  int counter = 0;
  log << MSG::DEBUG << "Ethan:: before Trigger Loop"<<endreq;

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

    log << MSG::DEBUG << "Johanna:: eta = "<<(*TriggerTowerIterator)->eta()<<endreq;
    log << MSG::DEBUG << "Johanna:: phi = "<<(*TriggerTowerIterator)->phi()<<endreq;

    m_h_TT_channels[0][0]->Fill((*TriggerTowerIterator)->emEnergy(), 1.); 

    m_h_TT_Em_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_TT_Had_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 
    m_h_TT_eta->Fill( (*TriggerTowerIterator)->eta(), 1.); 
    m_h_TT_phi->Fill( (*TriggerTowerIterator)->phi(), 1.); 

    m_h_TT_etaphi_hitmap->Fill( (*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi(), 1.);
    m_h_TT_etaphi->Fill( (*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi(), (*TriggerTowerIterator)->emEnergy()+(*TriggerTowerIterator)->hadEnergy());

    //Calibration plots
    m_h_Calib_TTEM_EtEta->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->emEnergy(),1.);//Ethan et->phi
    m_h_Calib_TTHAD_EtEta->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->hadEnergy(),1.);

    //Fill Combined Energy plot
    m_h_TT_Tot_Et->Fill( (*TriggerTowerIterator)->hadEnergy()+(*TriggerTowerIterator)->emEnergy(), 1.); 

    m_h_TT_Em10_Et->Fill( (*TriggerTowerIterator)->emEnergy(), 1.); 
    m_h_TT_Had10_Et->Fill( (*TriggerTowerIterator)->hadEnergy(), 1.); 

    //   m_h_TT_key->Fill((*TriggerTowerIterator)->key(), 1.);

    double eta = (*TriggerTowerIterator)->eta();

    eta*=eta;
    eta=sqrt(eta);

    //Et (EM+Had) >10 :
    if (((*TriggerTowerIterator)->emEnergy()+(*TriggerTowerIterator)->hadEnergy())>10){
      //m_h_TT_key->Fill((*TriggerTowerIterator)->key(), 1.);

      m_h_TT_eta_gt10->Fill( (*TriggerTowerIterator)->eta(),1.);      
      m_h_TT_phi_gt10->Fill( (*TriggerTowerIterator)->phi(), 1.); 
      if(1.5 < eta && eta< 3.2)m_h_TT_EC_phi_gt10->Fill( (*TriggerTowerIterator)->phi(), 1.); 
      if(0.1 < eta && eta < 1.5)m_h_TT_Barrel_phi_gt10->Fill( (*TriggerTowerIterator)->phi(), 1.); 
      if(3.2 < eta && eta< 4.9)m_h_TT_FCAL_phi_gt10->Fill( (*TriggerTowerIterator)->phi(), 1.); 
    }

    //Calo Region Plots

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
  
  
  // ================================================================================
  // ============================= JetElements ======================================
  // ================================================================================
  
  // JEs are only filled/stored if they contain energy, 
  // so looping over all JEs will only loop over those with energy
  log << MSG::DEBUG << "Ethan:: before JE Loop"<<endreq;
  JetElementCollection::const_iterator JetElementIterator    = JetElementTES->begin();
  JetElementCollection::const_iterator JetElementIteratorEnd = JetElementTES->end();
  log << MSG::DEBUG << "Ethan:: before JE Trigger Loop"<<endreq;
  for (; JetElementIterator != JetElementIteratorEnd; ++JetElementIterator)
  {
    m_h_JE_Em_Et->Fill( (*JetElementIterator)->emEnergy(), 1.);
    m_h_JE_Had_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.);
    m_h_JE_eta->Fill( (*JetElementIterator)->eta(), 1.);
    m_h_JE_phi->Fill( (*JetElementIterator)->phi(), 1.);

    //Fill combined energy
    m_h_JE_Tot_Et->Fill( (*JetElementIterator)->hadEnergy()+(*JetElementIterator)->emEnergy(), 1.);

    m_h_JE_Em10_Et->Fill( (*JetElementIterator)->emEnergy(), 1.); 
    m_h_JE_Had10_Et->Fill( (*JetElementIterator)->hadEnergy(), 1.); 

    m_h_JE_etaphi->Fill( (*JetElementIterator)->eta(), (*JetElementIterator)->phi(),(*JetElementIterator)->emEnergy()+(*JetElementIterator)->hadEnergy());
    m_h_JE_etaphi_hitmap->Fill( (*JetElementIterator)->eta(), (*JetElementIterator)->phi(), 1.);

    double eta = (*JetElementIterator)->eta();

    eta*=eta;
    eta=sqrt(eta);

    //Et (Em+Had) > 10:

    if (((*JetElementIterator)->emEnergy()+(*JetElementIterator)->hadEnergy())>10){
      m_h_JE_eta_gt10->Fill( (*JetElementIterator)->eta(),1.);      
      m_h_JE_phi_gt10->Fill( (*JetElementIterator)->phi(), 1.); 
      if(1.5 < eta && eta< 3.2)m_h_JE_EC_phi_gt10->Fill( (*JetElementIterator)->phi(), 1.); 
      if(0.1 < eta && eta < 1.5)m_h_JE_Barrel_phi_gt10->Fill( (*JetElementIterator)->phi(), 1.); 
      if(3.2 < eta && eta< 4.9)m_h_JE_FCAL_phi_gt10->Fill( (*JetElementIterator)->phi(), 1.); 
    }

    //Calo Region Plots

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
  
  return StatusCode( StatusCode::SUCCESS );

}
//_______________________________ proc  Histograms ___________________________________________
StatusCode
TrigT1CaloBSMonTool::
procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
{
        if( isEndOfEventsBlock || isEndOfLumiBlock ) 
	  {

	}
	
	if( isEndOfRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}

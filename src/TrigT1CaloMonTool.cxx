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

#include "CaloEvent/CaloCell.h"
#include "CaloEvent/CaloCellContainer.h"

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

//to form my own TT from calo cells
#include "TrigT1Calo/InternalTriggerTower.h"
#include "CaloEvent/CaloSampling.h"

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
  
  //  log << MSG::INFO << "ETHAN in bookHists()" << endreq;
  //  log << MSG::INFO << "Using base path " << m_stem << endreq;
    
  // Bin TT & JE Histos

  m_h_TT_Em_Et = book1D("TT_EM_Et","TT EM Et",255,0,255);
  m_h_TT_Had_Et = book1D("TT_HAD_Et","TT HAD Et",255,0,255);
  m_h_TT_eta = book1D("TT_eta","Trigger Tower eta",100,-5,5);
  //otherwise consider plotting eta in calo regions .
  m_h_TT_phi = book1D("TT_phi","Trigger Tower phi ",63,0,2*M_PI);

  m_h_JE_Em_Et = book1D("JE_EM_Et","JE EM Et",255,0,255);
  m_h_JE_Had_Et = book1D("JE_HAD_Et","JE HAD Et",255,0,255);
  m_h_JE_eta = book1D("JE_eta","JE eta",100,-5,5);
  m_h_JE_phi = book1D("JE_phi","JE phi ",64,0,2*M_PI);

  m_h_TT_Tot_Et = book1D("TT_Tot_Et","TT EM+HAD Et 0-10",10,0,10);
  m_h_JE_Tot_Et = book1D("JE_Tot_Et","JE EM+HAD Et 0-10",10,0,10);

  m_h_TT_Em10_Et = book1D("TT_EM10_Et","TT EM Et 0-10",10,0,10);
  m_h_TT_Had10_Et = book1D("TT_HAD10_Et","TT HAD Et 0-10",10,0,10);

  m_h_JE_Em10_Et = book1D("JE_EM10_Et","JE EM Et 0-10",10,0,10); 
  m_h_JE_Had10_Et = book1D("JE_HAD10_Et","JE HAD Et 0-10",10,0,10); 

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

  m_h_Barrel10_TT_Em_Et = book1D("Barrel10_TT_EM_Et","Barrel TT EM Et 0-10",10,0,10);
  m_h_Barrel10_TT_Had_Et = book1D("Barrel10_TT_HAD_Et","Barrel TT HAD Et 0-10",10,0,10);

  m_h_Barrel10_JE_Em_Et = book1D("Barrel10_JE_EM_Et","Barrel JE EM Et 0-10",10,0,10);
  m_h_Barrel10_JE_Had_Et = book1D("Barrel10_JE_HAD_Et","Barrel JE HAD Et 0-10",10,0,10);

  m_h_FCAL10_TT_Em_Et = book1D("FCAL10_TT_EM_Et","FCAL TT EM Et 0-10",10,0,10);
  m_h_FCAL10_TT_Had_Et = book1D("FCAL10_TT_HAD_Et","FCAL TT HAD Et 0-10",10,0,10);

  m_h_FCAL10_JE_Em_Et = book1D("FCAL10_JE_EM_Et","FCAL JE EM Et 0-10",10,0,10);
  m_h_FCAL10_JE_Had_Et = book1D("FCAL10_JE_HAD_Et","FCAL JE HAD Et 0-10",10,0,10);

  m_h_EC_TT_phi = book1D("EC_TT_phi","EndCap_TT phi",64,0,2*M_PI);
  m_h_EC_TT_Em_Et = book1D("EC_TT_EM_Et","EndCap_TT EM Et",255,0,255);
  m_h_EC_TT_Had_Et = book1D("EC_TT_HAD_Et","EndCap_TT HAD Et",255,0,255); 

  m_h_EC_JE_phi = book1D("EC_JE_phi","EndCap_JE phi",64,0,2*M_PI);
  m_h_EC_JE_Em_Et = book1D("EC_JE_EM_Et","EndCap_JE EM Et",255,0,255);
  m_h_EC_JE_Had_Et = book1D("EC_JE_HAD_Et","EndCap_JE HAD Et",255,0,255);

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

  m_h_TT_eta_gt10 = book1D("TT_eta_gt10",">10 Gev TT eta",100,-5,5);
  m_h_TT_phi_gt10 = book1D("TT_phi_gt10",">10 GeV TT phi",64,0,2*M_PI); 
  m_h_TT_EC_phi_gt10 = book1D("TT_EC_phi_gt10",">10 GeV TT EndCap phi",64,0,2*M_PI); 
  m_h_TT_Barrel_phi_gt10 = book1D("TT_Barrel_phi_gt10",">10 GeV TT Barrel phi",64,0,2*M_PI); 
  m_h_TT_FCAL_phi_gt10 = book1D("TT_FCAL_phi_gt10",">10 GeV TT FCAL phi",64,0,2*M_PI); 

  m_h_JE_eta_gt10 = book1D("JE_eta_gt10",">10 GeV JE eta",100,-5,5);
  m_h_JE_phi_gt10 = book1D("JE_phi_gt10",">10 GeV JE phi",64,0,2*M_PI); 
  m_h_JE_EC_phi_gt10 = book1D("JE_EC_phi_gt10",">10 GeV JE EndCap phi",64,0,2*M_PI); 
  m_h_JE_Barrel_phi_gt10 = book1D("JE_Barrel_phi_gt10",">10 GeV JE Barrel phi",64,0,2*M_PI); 
  m_h_JE_FCAL_phi_gt10 = book1D("JE_FCAL_phi_gt10",">10 GeV JE FCAL phi",64,0,2*M_PI); 

  // ====================================================================================================
  // Plots of the /Calorimeter information (combined LAr and Tile info)
  // Using CaloTowers from the ESD
  // ====================================================================================================

  m_h_CaloT_phi = book1D("CaloTower_phi","CaloTower phi",64,-M_PI,M_PI); 
  m_h_CaloT_eta = book1D("CaloTower_eta","CaloTower eta",100,-5,5); 
  m_h_CaloT_phi_gt10 = book1D("CaloTower_phi_gt10","CaloTower phi gt 10",64,-M_PI,M_PI); 
  m_h_CaloT_eta_gt10 = book1D("CaloTower_eta_gt10","CaloTower eta gt 10",100,-5,5); 

  m_h_CaloT_Et10 = book1D("CaloTower_Et10","CaloTower Et (0-10 GeV)",100,0,10.*GeV); //Note, need to specify GeV 
  m_h_CaloT_Et = book1D("CaloTower_Et","CaloTower Et",255,0,255.*GeV); 

  m_h_CaloT_etaphi = book2D("CaloTower_Eta_Phi","CaloTower Eta_Phi",100,-5,5, 64, -M_PI, M_PI);
  m_h_CaloT_etaphi_hitmap = book2D("CaloTower_Eta_Phi_HM","CaloTower Eta_Phi_HitMap",100,-5,5,64, -M_PI, M_PI);

  m_h_CaloT_key = book1D("CaloTower_Key","CaloTower Key",100,0,8000);


  // =====================================================================================================
  // =================================== Trigger Style Calo Plots ========================================
  // =====================================================================================================

  m_h_Calo_phi = book1D("Calo_phi","TSCT phi",64,0,2*M_PI); 
  m_h_Calo_eta = book1D("Calo_eta","TSCT eta",100,-5,5); 
  m_h_Calo_phi_gt10 = book1D("Calo_phi_gt10","TSCT phi Et > 10 GeV",64,0,2*M_PI); 
  m_h_Calo_eta_gt10 = book1D("Calo_eta_gt10","TSCT eta Et > 10 GeV",100,-5,5); 

  m_h_Calo_Em_Et10 = book1D("Calo_Em_Et10","TSCT Em Et 0-10 GeV",100,0,10.*GeV); //Note, need to specify GeV 
  m_h_Calo_Em_Et = book1D("Calo_Em_Et","TSCT Em Et",255,0,255.*GeV); 
  m_h_Calo_Had_Et10 = book1D("Calo_Had_Et10","TSCT Had Et 0-10 GeV",100,0,10.*GeV); //Note, need to specify GeV 
  m_h_Calo_Had_Et = book1D("Calo_Had_Et","TSCT Had Et",255,0,255.*GeV); 

  m_h_Calo_etaphi = book2D("Calo_Eta_Phi","TSCT Eta_Phi",100,-5,5, 64, 0, 2*M_PI);
  m_h_Calo_etaphi_hitmap = book2D("Calo_Eta_Phi_HM","TSCT Eta_Phi_HitMap",100,-5,5,64, 0, 2*M_PI);

  m_h_Calo_key = book1D("Calo_Key","TSCT Key",100,0,8000);

  // region plots

  m_h_Barrel_Calo_Em_Et = book1D("Barrel_Calo_Em_Et", "TSCT Barrel Et (Em)",100,0,100.*GeV);
  m_h_Barrel_Calo_Had_Et = book1D("Barrel_Calo_Had_Et", "TSCT Barrel Et (Had)",100,0,100.*GeV);
  m_h_Barrel10_Calo_Em_Et = book1D("Barrel10_Calo_Em_Et", "TSCT Barrel Et (Em, 0-10 GeV)",10,0,10.*GeV);
  m_h_Barrel10_Calo_Had_Et = book1D("Barrel10_Calo_Had_Et", "TSCT Barrel Et (Had, 0-10 GeV)",10,0,10.*GeV);
  m_h_Barrel_Calo_phi = book1D("Barrel_Calo_phi", "TSCT Barrel phi",64,0,2*M_PI);
  
  m_h_EC_Calo_Em_Et = book1D("EC_Calo_Em_Et", "TSCT EndCap Et (Em)",100,0,100.*GeV);
  m_h_EC_Calo_Had_Et = book1D("EC_Calo_Had_Et", "TSCT EndCap Et (Had)",100,0,100.*GeV);
  m_h_EC10_Calo_Em_Et = book1D("EC10_Calo_Em_Et", "TSCT EndCap Et (Em, 0-10 GeV)",10,0,10.*GeV);
  m_h_EC10_Calo_Had_Et = book1D("EC10_Calo_Had_Et", "TSCT EndCap Et (Had, 0-10 GeV)",10,0,10.*GeV);
  m_h_EC_Calo_phi = book1D("EC_Calo_phi", "TSCT EndCap phi",64,0,2*M_PI);
  
  m_h_FCAL_Calo_Em_Et = book1D("FCAL_Calo_Em_Et", "TSCT FCAL Et (Em)",100,0,100.*GeV);
  m_h_FCAL_Calo_Had_Et = book1D("FCAL_Calo_Had_Et", "TSCT FCAL Et (Had)",100,0,100.*GeV);
  m_h_FCAL10_Calo_Em_Et = book1D("FCAL10_Calo_Em_Et", "TSCT FCAL Et (Em, 0-10 GeV)",10,0,10.*GeV);
  m_h_FCAL10_Calo_Had_Et = book1D("FCAL10_Calo_Had_Et", "TSCT FCAL Et (Had, 0-10 GeV)",10,0,10.*GeV);
  m_h_FCAL_Calo_phi = book1D("FCAL_Calo_phi", "TSCT FCAL phi",64,0,2*M_PI);

  //comparisons of TT:Calo

  //1D Ratio plots //not implemented
  m_h_TT_Calo_eta = book1D("TT_Calo_eta","eta TT vs TSCT",100, -5, 5);
  m_h_TT_Calo_Et = book1D("TT_Calo_Et","Et TT vs TSCT", 255, 0, 255);

  //2D plotted one against the other
  m_h_TT_Calo_Em_EtTower = book2D("TT_Calo_Em_EtTower","TSCT and TT Tower Et (Em)",100,0,100.*GeV, 100, 0, 100);
  m_h_TT_Calo_Had_EtTower = book2D("TT_Calo_Had_EtTower","TSCT and TT Tower Et (Had)",100,0,100.*GeV, 100, 0, 100);
  m_h_TT_Calo_EtaTower = book2D("TT_Calo_EtaTower","TSCT and TT Tower Eta",100,-5,5, 100, -5, 5);
  m_h_TT_Calo_PhiTower = book2D("TT_Calo_PhiTower","TSCT and TT Tower Phi",64,0,2*M_PI,64,0,2*M_PI);

  //TSCT Calibration plots - see calibration below

  //Discrepancy check plots
  
  m_h_Ratio_D_Em_Et = book1D("Ratio_D_Em_Et","TSCT Et / TT Et (Em)",30,0,3); 
  m_h_Ratio_D_Had_Et = book1D("Ratio_D_Had_Et","TSCT Et / TT Et (Had)",30,0,3); 

  m_h_Calo_DEm_under_phi = book1D("Calo_DEm_under_phi","TSCT phi (TT<TSCT [Em])",64,0,2*M_PI); 
  m_h_Calo_DEm_under_eta = book1D("Calo_DEm_under_eta","TSCT eta (TT<TSCT [Em])",100, -5, 5); 
  m_h_Calo_DEm_under_Em_Et = book1D("Calo_DEm_under_Em_Et","Em Et TSCT (TT<TSCT [Em])", 255, 0, 255.*GeV);
  m_h_Calo_DEm_under_Had_Et = book1D("Calo_DEm_under_Had_Et","Had Et TSCT (TT<TSCT [Em])", 255, 0, 255.*GeV);
  m_h_Calo_DEm_under_TTEm_Et = book1D("Calo_DEm_under_TTEm_Et","Em Et TT (TT<TSCT [Em])", 255, 0, 255.*GeV);
  m_h_Calo_DEm_under_TTHad_Et = book1D("Calo_DEm_under_TTHad_Et","Had Et TT (TT<TSCT [Em])", 255, 0, 255.*GeV);

  m_h_Calo_DEm_over_phi = book1D("Calo_DEm_over_phi","TSCT phi (TT>TSCT [Em])",64,0,2*M_PI); 
  m_h_Calo_DEm_over_eta = book1D("Calo_DEm_over_eta","TSCT eta (TT>TSCT [Em])",100, -5, 5); 
  m_h_Calo_DEm_over_Em_Et = book1D("Calo_DEm_over_Em_Et","Em Et TSCT (TT>TSCT [Em])", 255, 0, 255.*GeV);
  m_h_Calo_DEm_over_Had_Et = book1D("Calo_DEm_over_Had_Et","Had Et TSCT (TT>TSCT [Em])", 255, 0, 255.*GeV);
  m_h_Calo_DEm_over_TTEm_Et = book1D("Calo_DEm_over_TTEm_Et","Em Et TT (TT>TSCT [Em])", 255, 0, 255.*GeV);
  m_h_Calo_DEm_over_TTHad_Et = book1D("Calo_DEm_over_TTHad_Et","Had Et TT (TT>TSCT [Em])", 255, 0, 255.*GeV);

  m_h_Calo_DHad_under_phi = book1D("Calo_DHad_under_phi","TSCT phi (TT<TSCT [Had])",64,0,2*M_PI); 
  m_h_Calo_DHad_under_eta = book1D("Calo_DHad_under_eta","TSCT eta (TT<TSCT [Had])",100, -5, 5); 
  m_h_Calo_DHad_under_Em_Et = book1D("Calo_DHad_under_Em_Et","Em Et TSCT (TT<TSCT [Had])", 255, 0, 255.*GeV);
  m_h_Calo_DHad_under_Had_Et = book1D("Calo_DHad_under_Had_Et","Had Et TSCT (TT<TSCT [Had])", 255, 0, 255.*GeV);
  m_h_Calo_DHad_under_TTEm_Et = book1D("Calo_DHad_under_TTEm_Et","Em Et TT (TT<TSCT [Had])", 255, 0, 255.*GeV);
  m_h_Calo_DHad_under_TTHad_Et = book1D("Calo_DHad_under_TTHad_Et","Had Et TT (TT<TSCT [Had])", 255, 0, 255.*GeV);

  m_h_Calo_DHad_over_phi = book1D("Calo_DHad_over_phi","TSCT phi (TT>TSCT [Had])",64,0,2*M_PI); 
  m_h_Calo_DHad_over_eta = book1D("Calo_DHad_over_eta","TSCT eta (TT>TSCT [Had])",100, -5, 5); 
  m_h_Calo_DHad_over_Em_Et = book1D("Calo_DHad_over_Em_Et","Em Et TSCT (TT>TSCT [Had])", 255, 0, 255.*GeV);
  m_h_Calo_DHad_over_Had_Et = book1D("Calo_DHad_over_Had_Et","Had Et TSCT (TT>TSCT [Had])", 255, 0, 255.*GeV);
  m_h_Calo_DHad_over_TTEm_Et = book1D("Calo_DHad_over_TTEm_Et","Em Et TT (TT>TSCT [Had])", 255, 0, 255.*GeV);
  m_h_Calo_DHad_over_TTHad_Et = book1D("Calo_DHad_over_TTHad_Et","Had Et TT (TT>TSCT [Had])", 255, 0, 255.*GeV);


  // ====================================================================================
  // Calibration Plots:
  // ====================================================================================


  //TT:
  m_h_Calib_TTEM_EtEta = book2D("TTEM_EtEta", "TT (EM) Et vs eta", 100, -5, 5, 10, 0, 10);
  m_h_Calib_TTHAD_EtEta = book2D("TTHAD_EtEta", "TT (HAD) Et vs eta", 100, -5, 5, 10, 0, 10);
  //CaloTower
  m_h_Calib_CaloT_EtEta = book2D("CaloT_EtEta", "ESD CaloTower Et vs eta", 100, -5, 5, 100, 0, 10.*GeV);
  
  //need to add TSCT Em and Had
  m_h_Calib_CaloEM_EtEta = book2D("CaloEM_EtEta", "TSCT (EM) Et vs eta", 100, -5, 5, 100, 0, 10.*GeV);
  m_h_Calib_CaloHAD_EtEta = book2D("CaloHAD_EtEta", "TSCT (HAD) Et vs eta", 100, -5, 5, 100, 0, 10.*GeV);
  
  //Not implemented
  m_h_Calib_EMRatio_ETEta = book2D("CalibEM_EtEta", "Calib (EM) Et vs eta", 100, -5, 5, 2, 0, 2);
  m_h_Calib_HADRatio_ETEta = book2D("CalibHAD_EtEta", "Calib (HAD) Et vs eta", 100, -5, 5, 2, 0, 2);
  


  SetBookStatus(true);
  
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloMonTool::fillHists()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  
  log << MSG::DEBUG << "Ethan in fillHists()" << endreq;

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
  
  //Retrieve Calo Tower collection from SG
  const CaloTowerContainer* CaloTowerTES = 0; 
  sc=m_StoreGate->retrieve(CaloTowerTES, "CombinedTower"); 
  if( sc.isFailure()  ||  !CaloTowerTES ) {
    log << MSG::ERROR<< "No CaloTowerContainer found"<< endreq; 
  }

  //Retreive Calo Cell collection from SG
  const CaloCellContainer* CaloCellTES = 0;
  sc=m_StoreGate->retrieve(CaloCellTES, "AllCalo"); 
  if( sc.isFailure()  ||  !CaloCellTES ) {
    log << MSG::ERROR<< "No CaloCellContainer found"<< endreq; 
  }else{
    //    log << MSG::ERROR<< "ETHAN :: CaloCellContainer found"<< endreq; 
  }
  

  //Key Stuff for Trigger vs Calo Tower comparison
  TriggerTowerKey TTKey(0,0);
  JetElementKey* JEKey; 
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
  //  std::cout<<"before Trigger Loop"<<std::endl;
  //  log << MSG::ERROR << "before Trigger Loop"<<endreq;

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
      m_h_TT_key->Fill((*TriggerTowerIterator)->key(), 1.);

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

  JetElementCollection::const_iterator JetElementIterator    = JetElementTES->begin();
  JetElementCollection::const_iterator JetElementIteratorEnd = JetElementTES->end();

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

  // =============================================================================================
  // ================= CaloTowers (combined LAr and Tile data) From ESD ==========================
  // =============================================================================================

  // CaloTowers exist irrespective of how low the Energy is, 
  // so need an Et cut to select only those towers with a deposited energy

  CaloTowerContainer::const_iterator CaloTowerIterator    = CaloTowerTES->begin();
  CaloTowerContainer::const_iterator CaloTowerIteratorEnd = CaloTowerTES->end();

  int TTkey; //change this to key and use previous int
  double TTtoCaloEnergyRatio;

 for (; CaloTowerIterator != CaloTowerIteratorEnd; ++CaloTowerIterator)
  {

    //select only Towers with an energy deposit
    if((*CaloTowerIterator)->et()>0.*GeV){
      m_h_CaloT_phi->Fill( (*CaloTowerIterator)->phi(), 1.);
      m_h_CaloT_eta->Fill( (*CaloTowerIterator)->eta(), 1.);

      if((*CaloTowerIterator)->et()>1.*GeV){

	m_h_Calib_CaloT_EtEta->Fill( (*CaloTowerIterator)->eta(), (*CaloTowerIterator)->et(), 1.);


	// Everything Below is Calo-TT comparison, should be in the TSCT part
	//
	//determine phi in the 0->2Pi range so that it will work with TT functionality
	double calophi = (*CaloTowerIterator)->phi();
	if(calophi>M_PI)calophi-=2*M_PI; //Fixme, wrong way round
	TTkey = TTKey.ttKey(calophi,(*CaloTowerIterator)->eta());
	m_h_CaloT_key->Fill(TTkey,1);
	m_h_CaloT_phi_gt10->Fill( (*CaloTowerIterator)->phi(), 1.);
	m_h_CaloT_eta_gt10->Fill( (*CaloTowerIterator)->eta(), 1.);
	m_h_CaloT_etaphi_hitmap->Fill( (*CaloTowerIterator)->eta(), (*CaloTowerIterator)->phi(), 1.);
	m_h_CaloT_etaphi->Fill( (*CaloTowerIterator)->eta(), (*CaloTowerIterator)->phi(), (*CaloTowerIterator)->et());
	
	m_h_CaloT_Et->Fill( (*CaloTowerIterator)->et(), 1.);
      }//eo calotower et>Xgev
    }//eo calotower et>0gev

    m_h_CaloT_Et10->Fill( (*CaloTowerIterator)->et(), 1.);

  }
 
 // ==============================================================================================
 // ================================ Trigger-Style-Calo-Towers ===================================
 // ==============================================================================================


 //Start by creating Trigger-Style-Calo-Towers - use a TrigT1Calo objext, InternalTriggerTower
 
 InternalTriggerTower* TriggerStyleCaloTower = 0;
 //For each event create a map for TriggerStyleCaloTowers
 std::map<int, InternalTriggerTower*> * TriggerStyleCaloTowerContainer = new std::map<int, InternalTriggerTower*>;

 CaloCellContainer::const_iterator CaloCellIterator    = CaloCellTES->begin();
 CaloCellContainer::const_iterator CaloCellIteratorEnd = CaloCellTES->end();
 for (; CaloCellIterator != CaloCellIteratorEnd; ++CaloCellIterator)
   {
     const CaloDetDescrElement* cell_DDE = (*CaloCellIterator)->caloDDE(); // check acceptedCaloCell arguments
     cell_DDE = (*CaloCellIterator)->caloDDE(); // check acceptedCaloCell arguments

     //em is cell_DDE->getSampling()== CaloCell_ID::FCAL0 || LAREM
     //had is everything else

     //check if cell is in a tile gap
     CaloSampling::CaloSample  sampl = CaloSampling::getSampling(**CaloCellIterator);
     if ( cell_DDE->getSubCalo()!=CaloCell_ID::TILE || CaloSampling::TileGap3 != sampl ){
       //Either not Tile, OR if it is a Tile, it's not in the TileGap

       //Obtain TT Key from CaloCell

       double calocellphi = (*CaloCellIterator)->phi();
       //convert calocellphi from -pi->pi range to TT 0->2pi
       if(calocellphi<0)calocellphi+=2*M_PI;


       TTkey = TTKey.ttKey(calocellphi,(*CaloCellIterator)->eta());
       double tt_phi = TTKey.phi();
       double tt_eta = TTKey.eta();

       //check to see if TriggerStyleCaloTower exists or need to add it to the map of TriggerStyleCaloTowers
       std::map<int, InternalTriggerTower*>::iterator TriggerStyleCaloTowerIterator = TriggerStyleCaloTowerContainer->find( TTkey );
       if (TriggerStyleCaloTowerIterator == TriggerStyleCaloTowerContainer->end()){//i.e. no TriggerStyleCaloTower at that key
	 TriggerStyleCaloTower = new InternalTriggerTower(tt_phi,tt_eta, TTkey);
	 TriggerStyleCaloTowerContainer->insert(std::map<int, InternalTriggerTower*>::value_type(TTkey,TriggerStyleCaloTower));//Add TT to map
       }else{//TriggerStyleCaloTower already exists at that key
	 TriggerStyleCaloTower=(TriggerStyleCaloTowerIterator->second);
       }

       if ( cell_DDE->getSubCalo()==CaloCell_ID::LAREM ||  
	    cell_DDE->getSampling()==CaloCell_ID::FCAL0 ) {//Em 
	 TriggerStyleCaloTower->addEMPeak((*CaloCellIterator)->et());
	 std::vector<double> emamp = TriggerStyleCaloTower->EmAmps();
	 //log << MSG::ERROR << "emamp.at3 = "<<emamp.at(3)<<endreq ;
	 //	 if(emamp.at(3)<-0.2*GeV) log << MSG::ERROR <<// "-ve CellEt = "<<(*CaloCellIterator)->et()<<
	 //  " emamp.at3= "<<emamp.at(3)*GeV<<endreq ;
	   
       }else{//Had
	 TriggerStyleCaloTower->addHadPeak((*CaloCellIterator)->et());
	 //std::vector<double> hadamp = TriggerStyleCaloTower->HadAmps();
	 // log << MSG::ERROR << "hadamp.at(3)= "<<hadamp.at(3)<<endreq ;

       }
     }//EO gap check     
   }
 
 //Now loop over TT and get key, compare to internalTT - which are in fact calo-like towers.


 // ========================================================================================
 //  Fill Plots with TriggerTower and Trigger-Style-Calo-Tower Comparisons
 // ========================================================================================
 /*
 //put TriggerTowerIterator back to the start (previously used in TT plots)
 TriggerTowerIterator    = TriggerTowerTES->begin(); 

 //use already defined TriggerTowerIterator

 for (; TriggerTowerIterator != TriggerTowerIteratorEnd; ++TriggerTowerIterator) 
   {
     //Find Key for the Trigger Tower     
     key = TTKey.ttKey((*TriggerTowerIterator)->phi(),(*TriggerTowerIterator)->eta());

     //Get the TriggerStyleCaloTower for that Key
     std::map<int, InternalTriggerTower*>::iterator TriggerStyleCaloTowerIterator = TriggerStyleCaloTowerContainer->find(key);

     //Get EmAmps vector for the TriggerStyleCaloTower
     //Get 3rd (peak) element
     m_h_TT_Calo_over_EtTower->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),(*TriggerTowerIterator)->emEnergy(),1);

     //log << MSG::ERROR << "emamp "<<(TriggerStyleCaloTowerIterator->second->EmAmps()).at(3)<<" EMEt "<<(*TriggerTowerIterator)->emEnergy()<<endreq ;

   }
 */
 
 std::map<int, InternalTriggerTower*>::iterator TriggerStyleCaloTowerIterator = TriggerStyleCaloTowerContainer->begin();
 std::map<int, InternalTriggerTower*>::iterator TriggerStyleCaloTowerIteratorEnd = TriggerStyleCaloTowerContainer->end();

 for (; TriggerStyleCaloTowerIterator != TriggerStyleCaloTowerIteratorEnd; ++TriggerStyleCaloTowerIterator){

   double EmEt = (TriggerStyleCaloTowerIterator->second->EmAmps()).at(3);
   double HadEt = (TriggerStyleCaloTowerIterator->second->HadAmps()).at(3);
   double eta = TriggerStyleCaloTowerIterator->second->eta();
   
   if(EmEt>1.*GeV){
     //     log << MSG::ERROR << "emamp.at3 = "<<(TriggerStyleCaloTowerIterator->second->EmAmps()).at(3)<<endreq ;

   m_h_Calo_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1.);
   m_h_Calo_phi->Fill(TriggerStyleCaloTowerIterator->second->phi(),1.);
   m_h_Calo_eta->Fill(TriggerStyleCaloTowerIterator->second->eta(),1.);
   m_h_Calib_CaloEM_EtEta->Fill(TriggerStyleCaloTowerIterator->second->eta(),(TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1.);

   }
   if(HadEt>1.*GeV){
     m_h_Calo_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1.);
     m_h_Calib_CaloHAD_EtEta->Fill(TriggerStyleCaloTowerIterator->second->eta(),(TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1.);
   }

   //Et>10 GeV     
   if((HadEt+EmEt)>10.*GeV){
    
     m_h_Calo_phi_gt10->Fill(TriggerStyleCaloTowerIterator->second->phi(),1.);
     m_h_Calo_eta_gt10->Fill(TriggerStyleCaloTowerIterator->second->eta(),1.);

   }


   //Only want towers with energy, so arbitrarily call noise Et<1 GeV
   if((HadEt+EmEt)>1.*GeV){
     m_h_Calo_etaphi_hitmap->Fill( TriggerStyleCaloTowerIterator->second->eta(),TriggerStyleCaloTowerIterator->second->phi(), 1.);
     m_h_Calo_etaphi->Fill(TriggerStyleCaloTowerIterator->second->eta(),TriggerStyleCaloTowerIterator->second->phi(), HadEt+EmEt );

   m_h_Calo_Em_Et10->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1.);
   m_h_Calo_Had_Et10->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1.);


   //Region Plots

   if(1.5 < eta && eta< 3.2){//EndCap

     m_h_EC_Calo_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1.);
     m_h_EC_Calo_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1.);
     m_h_EC_Calo_phi->Fill(TriggerStyleCaloTowerIterator->second->phi(),1.);
     if((HadEt+EmEt)>10.*GeV){
       m_h_EC10_Calo_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1.);
       m_h_EC10_Calo_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1.);
     }
   }else if(0.1 < eta && eta < 1.5){//Barrel

     m_h_Barrel_Calo_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1.);
     m_h_Barrel_Calo_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1.);
     m_h_Barrel_Calo_phi->Fill(TriggerStyleCaloTowerIterator->second->phi(),1.);
     if((HadEt+EmEt)>10.*GeV){
       m_h_Barrel10_Calo_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1.);
       m_h_Barrel10_Calo_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1.);
     }

   }else if(3.2 < eta && eta< 4.9){//FCAL
  
     m_h_FCAL_Calo_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1.);
     m_h_FCAL_Calo_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1.);
     m_h_FCAL_Calo_phi->Fill(TriggerStyleCaloTowerIterator->second->phi(),1.);
     if((HadEt+EmEt)>10.*GeV){
       m_h_FCAL10_Calo_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1.);
       m_h_FCAL10_Calo_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1.);
     }
   }

   }//EO >1GeV noise cut

   //TT and TSCT comparison plots ============================================

   TTkey = TTKey.ttKey(TriggerStyleCaloTowerIterator->second->phi(),eta);
   m_h_Calo_key->Fill(TTkey,1);

   //loop over triggertowers and make tt:calo comparison plots
   std::map<int, TriggerTower *>::iterator test=m_ttContainer->find( TTkey ); 

   if (test != m_ttContainer->end()){//If TT exists

     m_h_TT_Calo_Em_EtTower->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),(test->second)->emEnergy(),1);
     m_h_TT_Calo_Had_EtTower->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),(test->second)->hadEnergy(),1);
     //Fixme, want one for Em and one for Had
     m_h_TT_Calo_PhiTower->Fill(TriggerStyleCaloTowerIterator->second->phi(), (test->second)->phi(),1);
     m_h_TT_Calo_EtaTower->Fill(eta,(test->second)->eta(),1);

     if((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3)>0.8*GeV && (test->second)->emEnergy()>0){
       double TSCToverTTEmEnergyRatio = ((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3)/1000)
	                                                                   /(test->second)->emEnergy();
       //log << MSG::ERROR << "Ethan TSCToverTTEmEnergyRatio = " << TSCToverTTEmEnergyRatio <<endreq;

       m_h_Ratio_D_Em_Et->Fill(TSCToverTTEmEnergyRatio,1);

       if(TSCToverTTEmEnergyRatio>1.2){//TSCT>TT -> "under"
	 m_h_Calo_DEm_under_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1);
	 m_h_Calo_DEm_under_phi->Fill(TriggerStyleCaloTowerIterator->second->phi(),1);
	 m_h_Calo_DEm_under_eta->Fill(eta,1);
	 m_h_Calo_DEm_under_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1); 
	 m_h_Calo_DEm_under_TTEm_Et->Fill((test->second)->emEnergy()); 
	 m_h_Calo_DEm_under_TTHad_Et->Fill((test->second)->hadEnergy()); 
       }
       else if (TSCToverTTEmEnergyRatio<0.8){//TSCT<TT -> "over"
       	 m_h_Calo_DEm_over_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1);
	 m_h_Calo_DEm_over_phi->Fill(TriggerStyleCaloTowerIterator->second->phi(),1);
	 m_h_Calo_DEm_over_eta->Fill(eta,1);
	 m_h_Calo_DEm_over_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1); 
	 m_h_Calo_DEm_over_TTEm_Et->Fill((test->second)->emEnergy()); 
	 m_h_Calo_DEm_over_TTHad_Et->Fill((test->second)->hadEnergy()); 
       }
     }

     if((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3)>0.8*GeV && (test->second)->hadEnergy()>0){
       double TSCToverTTHadEnergyRatio = ((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3)/1000)
	                                                                   /(test->second)->hadEnergy();

       m_h_Ratio_D_Had_Et->Fill(TSCToverTTHadEnergyRatio,1);

       if(TSCToverTTHadEnergyRatio>1.2){//TSCT>TT -> "under"
	 m_h_Calo_DHad_under_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1);
	 m_h_Calo_DHad_under_phi->Fill(TriggerStyleCaloTowerIterator->second->phi(),1);
	 m_h_Calo_DHad_under_eta->Fill(eta,1);
	 m_h_Calo_DHad_under_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1); 
	 m_h_Calo_DHad_under_TTEm_Et->Fill((test->second)->emEnergy()); 
	 m_h_Calo_DHad_under_TTHad_Et->Fill((test->second)->hadEnergy()); 
       }
       else if (TSCToverTTHadEnergyRatio<0.8){//TSCT<TT -> "over"
       	 m_h_Calo_DHad_over_Em_Et->Fill((TriggerStyleCaloTowerIterator->second->EmAmps()).at(3),1);
	 m_h_Calo_DHad_over_phi->Fill(TriggerStyleCaloTowerIterator->second->phi(),1);
	 m_h_Calo_DHad_over_eta->Fill(eta,1);
	 m_h_Calo_DHad_over_Had_Et->Fill((TriggerStyleCaloTowerIterator->second->HadAmps()).at(3),1); 
	 m_h_Calo_DHad_over_TTEm_Et->Fill((test->second)->emEnergy()); 
	 m_h_Calo_DHad_over_TTHad_Et->Fill((test->second)->hadEnergy()); 
       }

     }

   }

 }
 
  return StatusCode::SUCCESS;

}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloMonTool::finalHists()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  //   log << MSG::INFO << "Ethan in finalHists()" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TrigT1CaloMonTool::checkHists(bool /* fromFinalize */)
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  //  log << MSG::INFO << "Ethan in checkHists()" << endreq;

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

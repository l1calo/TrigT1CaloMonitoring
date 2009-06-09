// ********************************************************************
//
// NAME:        JEPTransPerfMon.cxx
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the transmission and performance in the JEP
//              further documentaion at https://twiki.cern.ch/twiki/bin/view/Atlas/TrigT1CaloMonitoring
//

// done checks:

// Transmission: only data path (no simulation)
// Functionality: comparison with simulation 

// ================= TriggerTowers -> JetElements ==============================================
// Check	Transmission & Functionality
//              JetElements(PPrDAQ.TriggerTowers) = JEMDAQ.JetElements

// ================= JetElements -> JEM Jet Hits ===============================================
// Check	Functionality
//              JEMJetHits(JEMDAQ.JetElements) = JEMDAQ.JEMJetHits

// ================= JetElements -> JEM Energy Sums ============================================
// Check	Functionality
//              JEMEnergySums(JEMDAQ.JetElements) = JEMDAQ.JEMEnergySums

// ================= JEM Jet Hits -> (Crate/System) CMM Jet Hits ===============================
// Check	Transmission
//              JEMDAQ.JEMJetHits = (Crate/System)CMMJetDAQ.JEMJetHits
// Check	Functionality
//              CMMJetHits(JEMDAQ.JEMJetHits) = (Crate/System)CMMJetDAQ.CrateCMMJetHits

// ================= JEMEnergySums -> (Crate/System) CMM Energy Sums ===========================
// Check	Transmission
//              JEMDAQ.JEMEnergySums = (Crate/System)CMMEnergySums.JEMEnergySums
// Check	Functionality
//              CMMEnergySums(JEMDAQ.JEMEnergySums) = (Crate/System)CMMEnergyDAQ.CrateCMMEnergySums

// ================= Crate CMM Jet Hits -> System CMM Jet Hits =================================
// Check	Transmission
//              CrateCMMJetDAQ.CrateCMMJetHits = SystemCMMJetDAQ.RemoteCMMJetHits
// Check	Functionality
//              TotalCMMJetHits(CrateCMMJetDAQ.CrateCMMJetHits+ SystemCMMJetDAQ.CrateCMMJetHits) = SystemCMMJetDAQ.SystemCMMJetHits

// ================= Crate CMM Energy Sums -> System CMM Energy Sums ===========================
// Check	Transmission
//              CrateCMMEnergyDAQ.CrateCMMEnergySums = SystemCMMEnergyDAQ.RemoteCMMEnergySums
// Check	Functionality
//              TotalCMMEnergySums(CrateCMMEnergyDAQ.CrateCMMEnergySums+ SystemCMMEnergyDAQ.CrateCMMEnergySums) = SystemCMMEnergyDAQ.SystemCMMEnergySums

// ================= RoI and Level-2 ===========================================================
// Check	Transmission
//              JEPRoI.TotalJetEtHits = SystemCMMJetDAQ.TotalJetEtHits
//              JEPRoI.SystemCMMEnergySums = SystemCMMEnergyDAQ.SystemCMMEnergySums
// Check	Functionality
//              JEMJetRoI(JEMDAQ.JetElements) = JEPRoI.JEMRoIJetHits

// ================= JEP -> CTP ================= ==============================================
// Check	Transmission
//              SystemCMMJetDAQ.SystemCMMJetHits = JEPCTP.SystemCMMJetHits
//              SystemCMMEnergyDAQ.TotalEtMissHits. = JEPCTP.TotalEtMissHits
//              SystemCMMEnergyDAQ.TotalEtSumHits. = JEPCTP.TotalEtSumHits
// =============================================================================================

// ********************************************************************

#include <sstream>

#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "StoreGate/StoreGateSvc.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "SGTools/StlVectorClids.h"

#include <TROOT.h>
#include <TColor.h> 
#include <TCanvas.h> 
#include "TH1.h"
#include "TStyle.h"
#include <algorithm>
#include <math.h>
#include <functional>
#include <iostream>
#include <cstdlib>
#include <utility>

#include "TrigT1CaloUtils/JetAlgorithm.h"
#include "TrigT1CaloUtils/CoordToHardware.h"
#include "TrigT1CaloEvent/CMMJetHits.h"
#include "TrigT1CaloEvent/JEMHits.h"
#include "TrigT1CaloEvent/JetElement.h"
#include "TrigT1CaloEvent/JEMRoI.h"
#include "TrigT1CaloEvent/CMMRoI.h"
#include "TrigT1CaloEvent/JEMEtSums.h"
#include "TrigT1CaloEvent/CMMEtSums.h"
#include "TrigT1CaloEvent/TriggerTower.h"
#include "TrigT1CaloToolInterfaces/IL1JEPHitsTools.h"
#include "TrigT1CaloToolInterfaces/IL1JetTools.h"
#include "TrigT1CaloToolInterfaces/IL1JetElementTools.h"
#include "TrigT1CaloToolInterfaces/IL1JEPEtSumsTools.h"
#include "TrigT1Interfaces/Coordinate.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"

#include "TrigT1CaloMonitoring/JEPTransPerfMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "Identifier/HWIdentifier.h"
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"



namespace LVL1 {
 class CMMRoI;
}

namespace LVL1 {
 class JEMRoI;
}

// *********************************************************************
// Public Methods
// *********************************************************************

/*---------------------------------------------------------*/
JEPTransPerfMon::JEPTransPerfMon( const std::string & type, const std::string & name,
		const IInterface* parent )
   : ManagedMonitorToolBase( type, name, parent ),
    m_jepHitsTool("LVL1::L1JEPHitsTools/L1JEPHitsTools"),
    m_jetTool("LVL1::L1JetTools/L1JetTools"),
    m_jetElementTool("LVL1::L1JetElementTools/L1JetElementTools"),
    m_etSumsTool("LVL1::L1JEPEtSumsTools/L1JEPEtSumsTools"),
    mLog(msgSvc(),name)
/*---------------------------------------------------------*/
{  
  
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "JEPHitsTool", m_jepHitsTool);
  declareProperty( "JetTool", m_jetTool);
  declareProperty( "JetElementTool", m_jetElementTool);
  declareProperty( "JEPEtSumsTool", m_etSumsTool);
  
  declareProperty( "BS_JetElementLocation", m_BS_JetElementLocation = LVL1::TrigT1CaloDefs::JetElementLocation); 
  declareProperty( "BS_TriggerTowerLocation", m_BS_TriggerTowerLocation = LVL1::TrigT1CaloDefs::TriggerTowerLocation); 
  declareProperty( "NoLUTSlices", m_NoLUTSlices = 5); 

  declareProperty( "BS_JEMHitsLocation", m_BS_JEMHitsLocation =  LVL1::TrigT1CaloDefs::JEMHitsLocation) ;
  declareProperty( "Sim_JEMHitsLocation", m_Sim_JEMHitsLocation =  LVL1::TrigT1CaloDefs::JEMHitsLocation) ;

  declareProperty( "BS_JEMEtSumsLocation", m_BS_JEMEtSumsLocation=   LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;
  declareProperty( "Sim_JEMEtSumsLocation", m_Sim_JEMEtSumsLocation=   LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;

  declareProperty( "BS_JEMRoILocation", m_BS_JEMRoILocation=   LVL1::TrigT1CaloDefs::JEMRoILocation) ;
  declareProperty( "Sim_JEMRoILocation", m_Sim_JEMRoILocation=   LVL1::TrigT1CaloDefs::JEMRoILocation) ;

  declareProperty( "BS_CMMJetHitsLocation", m_BS_CMMJetHitsLocation =  LVL1::TrigT1CaloDefs::CMMJetHitsLocation) ;
  declareProperty( "Sim_CMMJetHitsLocation", m_Sim_CMMJetHitsLocation =  LVL1::TrigT1CaloDefs::CMMJetHitsLocation) ;

  declareProperty( "BS_CMMEtSumsLocation", m_BS_CMMEtSumsLocation =  LVL1::TrigT1CaloDefs::CMMEtSumsLocation) ;  
  declareProperty( "Sim_CMMEtSumsLocation", m_Sim_CMMEtSumsLocation =  LVL1::TrigT1CaloDefs::CMMEtSumsLocation) ;  

  declareProperty( "BS_CMMRoILocation", m_BS_CMMRoILocation =  LVL1::TrigT1CaloDefs::CMMRoILocation) ;
  declareProperty( "Sim_CMMRoILocation", m_Sim_CMMRoILocation =  LVL1::TrigT1CaloDefs::CMMRoILocation) ;

  declareProperty( "Offline", m_Offline = 1) ;
  declareProperty( "CompareWithSimulation", m_CompareWithSimulation = 1) ;


  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/TransAndPerf") ;
  
}


/*---------------------------------------------------------*/
JEPTransPerfMon::~JEPTransPerfMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode JEPTransPerfMon:: initialize()
/*---------------------------------------------------------*/
{
    
  MsgStream mLog( msgSvc(), name() );
   
  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  
  if (sc.isFailure()) 
    {
      mLog << MSG::ERROR
	   << "Unable to retrieve pointer to StoreGateSvc"
	   << endreq;
      return sc;
    }

  // Get a pointer to DetectorStore services
  sc = service("DetectorStore", m_detStore);
  if (sc.isFailure()) 
    {
      mLog << MSG::ERROR << "Cannot access DetectorStore" << endreq;
      return StatusCode::FAILURE;;
    }

  // Retrieve the CaloIdManager from the detector store
  sc = m_detStore->retrieve(m_caloMgr);
  if (sc.isFailure()) 
    {
      mLog << MSG::ERROR << "Unable to retrieve CaloIdManager from DetectorStore" << endreq;
      return StatusCode::FAILURE;
    }
  
  // Use the CaloIdManager to get a pointer to an instance of the CaloLVL1_ID helper
  m_lvl1Helper = m_caloMgr->getLVL1_ID();
  if (!m_lvl1Helper) 
    {
      mLog << MSG::ERROR << "Could not access CaloLVL1_ID helper" << endreq;
      return StatusCode::FAILURE;
    }

  IToolSvc* toolSvc;
  StatusCode status = service( "ToolSvc",toolSvc  );

  if(status.isSuccess()) 
    {
      IAlgTool *algtool;
      sc = toolSvc->retrieveTool("L1CaloTTIdTools", algtool);
      mLog<<MSG::DEBUG<<"L1CaloTTIdTools retrieved"<<endreq;
      if (sc!=StatusCode::SUCCESS) 
	{
	  mLog << MSG::ERROR << " Cannot get L1CaloTTIdTools !" << endreq;
	  return sc;
	}
      m_l1CaloTTIdTools = dynamic_cast<L1CaloTTIdTools*> (algtool);
    } 
  else 
    {
      return StatusCode::FAILURE;
    }

  
  ISvcLocator* svcLoc = Gaudi::svcLocator( );
  toolSvc = 0; // Pointer to Tool Service
  sc = svcLoc->service( "ToolSvc",toolSvc  );
  if(sc.isSuccess()) 
    {
      sc = toolSvc->retrieveTool("CaloTriggerTowerService",m_ttSvc);
      if(sc.isFailure())
	{
	  mLog << MSG::ERROR << "Could not retrieve CaloTriggerTowerService Tool" << endreq;
	  return StatusCode::FAILURE;
	}
    } 
  else 
    {
      mLog << MSG::ERROR << "Could not retrieve ToolSvc" << endreq;
      return StatusCode::FAILURE;
    }


  // Use the CaloIdManager to get a pointer to an instance of the TTOnlineID helper
  m_l1ttonlineHelper = m_caloMgr->getTTOnlineID();
  if (!m_l1ttonlineHelper ) 
    {
      mLog << MSG::ERROR << "Could not access TTOnlineID helper" << endreq;
      return StatusCode::FAILURE;
    }
  
  sc = m_jetElementTool.retrieve();
  if( sc.isFailure() ) {
    mLog << MSG::ERROR << "Unable to locate Tool L1JetElementTools" << endreq;
    return sc;
  }

  sc = m_jetTool.retrieve();
  if( sc.isFailure() ) {
    mLog << MSG::ERROR << "Unable to locate Tool L1JetTools" << endreq;
    return sc;
  }

  sc = m_jepHitsTool.retrieve();
  if( sc.isFailure() ) {
    mLog << MSG::ERROR << "Unable to locate Tool L1JEPHitsTools" << endreq;
    return sc;
  }

  sc = m_etSumsTool.retrieve();
  if( sc.isFailure() ) {
    mLog << MSG::ERROR << "Unable to locate Tool L1JEPEtSumsTools" << endreq;
    return sc;
  }

   
  
  sc = service( "StoreGateSvc", m_eventStore);
  if( sc.isFailure() ) {
    mLog << MSG::FATAL << name() << ": Unable to locate Service EventStore" << endreq;
    return sc;
  }

  mLog << MSG::INFO << "JEPTransPerfMon initialised" << endreq;

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode JEPTransPerfMon::bookHistograms( bool isNewEventsBlock, 
				   bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{


  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in JEPTransPerfMon::bookHistograms" << endreq;
 
  ManagedMonitorToolBase::LevelOfDetail_t LevelOfDetail=shift;
  if (m_DataType=="Sim") LevelOfDetail = expert;

  MonGroup CMM_transmission ( this, (m_PathInRootFile ).c_str(), shift, run );
  HistoBooker transmission_Booker(&CMM_transmission, &mLog, "");
  
  MonGroup SimBSComparison_JEM ( this, (m_PathInRootFile).c_str(), shift, run );
  HistoBooker JEM_Booker(&SimBSComparison_JEM, &mLog, m_DataType);

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  if ( isNewEventsBlock|| isNewLumiBlock) { }
  
  if( isNewRun ) 
    {
    
      const EventInfo* ByteStreamEventInfo;  
      StatusCode sc = m_eventStore->retrieve(ByteStreamEventInfo);
      if (sc.isFailure())  { 
        mLog << MSG::ERROR << "no event info" << endreq ;
        m_evtNum    = 0;
        m_lumiBlock = 0;
        m_evtBCID   = 0;   
        m_runNum    = 0;
      }

      else {
      EventID* eventID = ByteStreamEventInfo->event_ID();
      if (eventID) {
        m_evtNum    = eventID->event_number();
        m_lumiBlock = eventID->lumi_block();
        m_evtBCID   = eventID->bunch_crossing_id();
        m_runNum    = eventID->run_number();
      }
      else {	
        m_evtNum    =0;
        m_lumiBlock =0;
        m_evtBCID   =0;
        m_runNum    =0;
      }


      } 
      
    	
      m_NoEvents=0;

      std::string name;
      std::stringstream buffer,buffer2;
      
      
     
      m_h_TransCheck_emJetElements=transmission_Booker.book2F("emJE_TransPerfCheck","em JE Transmission and Performance Check for (TT|JE)", (1*m_NoLUTSlices),0.5,(1*m_NoLUTSlices+.5), 35,0.5,35.5, "", "");
      
      m_h_TransCheck_emJetElements->SetStats(kFALSE);
      //      m_h_TransCheck_emJetElements->SetOption(colz);
      m_h_TransCheck_hadJetElements=transmission_Booker.book2F("hadJE_TransPerfCheck", "had JE Transmission and Performance Check for (TT|JE)",(1*m_NoLUTSlices),0.5,(1*m_NoLUTSlices+.5), 35,0.5,35.5, "", "");
      m_h_TransCheck_hadJetElements->SetStats(kFALSE);
      //      m_h_TransCheck_hadJetElements->SetOption(colz);
      
      for (int i = 0; i < 16; i++)
	{
	  buffer.str("");
	  buffer<<i;
	  
	  name = "JEM " + buffer.str();
	  m_h_TransCheck_emJetElements->GetYaxis()->SetBinLabel((i+1), name.c_str());
	  m_h_TransCheck_emJetElements->GetYaxis()->SetBinLabel((i+1+18), name.c_str());
	  m_h_TransCheck_hadJetElements->GetYaxis()->SetBinLabel((i+1), name.c_str());
	  m_h_TransCheck_hadJetElements->GetYaxis()->SetBinLabel((i+1+18), name.c_str());
	  
	}
      m_h_TransCheck_emJetElements->GetYaxis()->SetBinLabel(17, "Crate 0");
      m_h_TransCheck_emJetElements->GetYaxis()->SetBinLabel(35, "Crate 1");
      m_h_TransCheck_hadJetElements->GetYaxis()->SetBinLabel(17, "Crate 0");
      m_h_TransCheck_hadJetElements->GetYaxis()->SetBinLabel(35, "Crate 1");
      
      
      int k=1;
      for (int i=2; i<3;i++)
	{
	  for (int j=0; j<m_NoLUTSlices; j++)
	    {
	      buffer.str("");
	      buffer<<i;
	      buffer2.str("");
	      buffer2<<j;
	      name = buffer2.str() +"|" + buffer.str();
	      
	      m_h_TransCheck_emJetElements->GetXaxis()->SetBinLabel(k, name.c_str());
	      m_h_TransCheck_hadJetElements->GetXaxis()->SetBinLabel(k, name.c_str());
	      k++;
	    }
	}
      
      mLog << MSG::DEBUG << "NoLUTSlices "<< m_NoLUTSlices << endreq ;      
      
      // Comparison with Simulation
      if (m_CompareWithSimulation ==1){ 	
	
	m_h_SimBSMon_JEP=JEM_Booker.book2F("JEP_Calc_Error", "JEP Hardware Output compared to Simulation: Differences", 4,0.5,4.5,37,0.5,37.5, "", "");
	m_h_SimBSMon_JEP->SetStats(kFALSE);
	//m_h_SimBSMon_JEP->SetOption(colz);
	
	m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(1, "CalcErrors E_{x}, E_{y}");
	m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(2, "CalcErrors E_{t}");
	m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(3, "CalcErrors Hits");
	m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(4, "CalcErrors RoIs");
	
	for (int i = 0; i < 16; i++)
	  {
	    buffer.str("");
	    buffer<<i;
	    
	    name = "JEM " + buffer.str();
	    m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel((i+1), name.c_str());
	    
	    buffer.str("");
	    buffer<<i;
	    
	    name = "JEM " + buffer.str();
	    m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel((i+1+19), name.c_str());
	  }
	m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel(17, "CMM");
	m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel(18, "Crate 0: ");
	m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel(36, "CMM");
	m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel(37, "Crate 1: ");
      }
      
      
      
      
      
      //---------------------------------- Backplane transmission checks -----------------------------
      m_h_TransCheck_JEP=transmission_Booker.book2F("JEP_TransCheck", "JEP Backplane Transmission Check JEM -> CMM per Module and Crate", 2,0.5,2.5,36,0.5,36.5, "", "");
      
      m_h_TransCheck_JEP->GetXaxis()->SetBinLabel(1, "TransErrors JetHits");
      m_h_TransCheck_JEP->GetXaxis()->SetBinLabel(2, "TransErrors E_{x}, E_{y}, E_{t}");
      m_h_TransCheck_JEP->SetStats(kFALSE);
      //m_h_TransCheck_JEP->SetOption(colz);
      for (int i = 0; i < 16; i++)
	{
	  buffer.str("");
	  buffer<<i;
	  
	  name = "JEM " + buffer.str();
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel((i+1), name.c_str());
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel((i+1+19), name.c_str());
	}
      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(17, "CMM");
      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(18, "Crate 0: ");
      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(36, "Crate 1: ");


    }
  return StatusCode( StatusCode::SUCCESS );
}


/*---------------------------------------------------------*/
StatusCode JEPTransPerfMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  Helper Help;
 
  m_NoEvents++;

  // Error vector for global overview
  std::vector<int> overview(2);

  // =============================================================================================
  // ================= TriggerTowers -> JetElements ==============================================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Transmission & Functionality
  //            JetElements(PPrDAQ.TriggerTowers) = JEMDAQ.JetElements
  // ---------------------------------------------------------------------------------------------
  
  // retrieve JetElements
  const JECollection* jetElements;
  StatusCode sc = m_storeGate->retrieve(jetElements, m_BS_JetElementLocation);
  
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JetElements found in TES at " << m_BS_JetElementLocation << endreq ;
      return StatusCode::SUCCESS;
    }


  if (m_CompareWithSimulation ==1){

    mLog << MSG::DEBUG << "==== TriggerTowers -> JetElements: Transmission & Functionality ===="<< endreq ;
    
    // retrieve JetElements
    const JECollection* jetElements;
    StatusCode sc = m_storeGate->retrieve(jetElements, m_BS_JetElementLocation);
    
    if( (sc==StatusCode::FAILURE) ) 
      {
      	mLog << MSG::INFO << "No JetElements found in TES at " << m_BS_JetElementLocation << endreq ;
      	return StatusCode::SUCCESS;
      }
    
    // retrieve JetElements
    const JECollection* TT_jetElements;
    sc = m_storeGate->retrieve(TT_jetElements, m_BS_TriggerTowerLocation);
    
    if( (sc==StatusCode::FAILURE) ) 
      {
	mLog << MSG::INFO << "No JetElements found in TES at " << m_BS_TriggerTowerLocation << endreq ;
	return StatusCode::SUCCESS;
      }
    
    JECollection::const_iterator it_je ;
    
    
    // Is this loop a no-op?
    for( it_je = jetElements ->begin(); it_je < jetElements->end(); ++it_je )
      {
	LVL1::CoordToHardware ToHW;
	LVL1::Coordinate coord((*it_je)->phi(),(*it_je)->eta());
	
	/*int crate =*/ ToHW.jepCrate(coord);
	/*int module=*/ToHW.jepModule(coord);
      }
     
    int k=1;
    int i=2; //Slice number of TT_JE (normally 2)
    
    for (int j=0; j<m_NoLUTSlices; j++)  //Slice number of TT_TS (overwritten by joboptions)
      {
	TimeSliceMatch(k, j, TT_jetElements, i, jetElements, overview, &mLog);
	k++;
      }
  }

  // =============================================================================================
  // ================= JetElements -> JEM Jet Hits ===============================================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            JEMJetHits(JEMDAQ.JetElements) = JEMDAQ.JEMJetHits
  // ---------------------------------------------------------------------------------------------
  
  const JemRoiCollection* BS_JEMRoI;
  sc = m_storeGate->retrieve(BS_JEMRoI, m_BS_JEMRoILocation);
  if( (sc==StatusCode::FAILURE) ) {
    mLog << MSG::INFO << "No JEMRoI found in TES at "<< m_BS_JEMRoILocation << endreq ;
    return StatusCode::SUCCESS;
  }
  
      
  const JEMHitsCollection* BS_JEMHits;
  sc = m_storeGate->retrieve(BS_JEMHits, m_BS_JEMHitsLocation);
  if( (sc==StatusCode::FAILURE) ) {
    mLog << MSG::INFO << "No JEMHits found in TES at "<< m_BS_JEMHitsLocation << endreq ;
    return StatusCode::SUCCESS;
    }
	
  
  if (m_CompareWithSimulation ==1){

    mLog << MSG::DEBUG << "==== JetElements -> JEM Jet Hits: Functionality ===="<< endreq ;
      
    const JEMHitsCollection* BS_JEMHits;
    sc = m_storeGate->retrieve(BS_JEMHits, m_BS_JEMHitsLocation);
    if( (sc==StatusCode::FAILURE) ) {
      mLog << MSG::INFO << "No JEMHits found in TES at "<< m_BS_JEMHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
      
    JEMHitsCollection* Sim_JEMHits = 0;
    if (BS_JEMRoI) {
      Sim_JEMHits = new JEMHitsCollection;
      mLog<<MSG::DEBUG<<"Simulate JEM Hits from JEM RoIs"<<endreq;
      m_jepHitsTool->formJEMHits(BS_JEMRoI, Sim_JEMHits);
    }
     
      
    std::vector <LVL1::JEMHits*>  vBS_JEMHits;
    std::vector <LVL1::JEMHits*>  vSim_JEMHits;
    
       //use standard vector instead of datavector for transmissioncheck:
      //datavector erases not only the pointer in the vector, but also the referenced object
      //-> segmentation fault!
      
      //fill vectors
    JEMHitsCollection::const_iterator it_JEMHits ;
    int  foundModule;
    int noMatchfound;
      
    for( it_JEMHits  = BS_JEMHits->begin(); it_JEMHits != BS_JEMHits->end(); ++it_JEMHits )
    {
      vBS_JEMHits.push_back(*it_JEMHits);	       
    }   
    for( it_JEMHits  = Sim_JEMHits->begin(); it_JEMHits != Sim_JEMHits->end(); ++it_JEMHits )
    {
      vSim_JEMHits.push_back(*it_JEMHits);	       
    }   
      
      // Step over all cells and compare
    std::vector <LVL1::JEMHits*>::iterator it_BS_JEMHits;
    std::vector <LVL1::JEMHits*>::iterator it_Sim_JEMHits;
      
      
    it_BS_JEMHits=vBS_JEMHits.begin();
      
    mLog<<MSG::DEBUG<<"JEMHits Calculation difference for"<<endreq;
    while (it_BS_JEMHits != vBS_JEMHits.end()) {
      foundModule = 0;
      noMatchfound=0;
      it_Sim_JEMHits=vSim_JEMHits.begin();
	  
	  while ((foundModule==0)and(it_Sim_JEMHits != vSim_JEMHits.end()))
	    {
	      
	      if (((*it_BS_JEMHits)->crate()==(*it_Sim_JEMHits)->crate())
		  and((*it_BS_JEMHits)->module()==(*it_Sim_JEMHits)->module()))
		{
		  foundModule=1;
		  
		  if ((*it_BS_JEMHits)->JetHits()!=(*it_Sim_JEMHits)->JetHits())
		    {
		      mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_JEMHits)->crate()<<" Module "<<(*it_BS_JEMHits)->module()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS:  Hits "<<Help.Binary((*it_BS_JEMHits)->JetHits(),24)<<endreq;
		      mLog<<MSG::VERBOSE<<"Sim: Hits "<<Help.Binary((*it_Sim_JEMHits)->JetHits(),24)<<endreq;
		      noMatchfound=1;
		    }
		  m_h_SimBSMon_JEP->Fill(3,((*it_BS_JEMHits)->crate()*19+(*it_BS_JEMHits)->module()+1),noMatchfound);
		  overview[(*it_BS_JEMHits)->crate()] |= (noMatchfound << 2);
		  
		  it_BS_JEMHits=vBS_JEMHits.erase(it_BS_JEMHits);
		  it_Sim_JEMHits=vSim_JEMHits.erase(it_Sim_JEMHits);
		}
	      else ++it_Sim_JEMHits;
	    }
	  if (foundModule==0)++it_BS_JEMHits;
	}
      
      if (vBS_JEMHits.size()!=0)
	{
	  mLog<<MSG::DEBUG<<"JEMHits: additional BS data for"<<endreq;
	  
	  //fill errorcounter
	  for( it_BS_JEMHits  = vBS_JEMHits.begin(); it_BS_JEMHits != vBS_JEMHits. end(); ++it_BS_JEMHits )
	    {
	      mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_JEMHits)->crate()<<" Module "<<(*it_BS_JEMHits)->module()<<endreq;
	      mLog<<MSG::VERBOSE<<"BS: Hits "<<Help.Binary((*it_BS_JEMHits)->JetHits(),24)<<endreq;
	      
	      if ((*it_BS_JEMHits)->JetHits()!=0) {
              m_h_SimBSMon_JEP->Fill(3,((*it_BS_JEMHits)->crate()*19+(*it_BS_JEMHits)->module()+1),1);	    
              overview[(*it_BS_JEMHits)->crate()] |= (1 << 2);
	      }
	    }
	}
      
      if (vSim_JEMHits.size()!=0)
	{
	  mLog<<MSG::DEBUG<<"JEMHits: additional Sim data for"<<endreq;
	  mLog<<MSG::DEBUG<<"fill error counter but remember simulation is not zero suppressing!"<<endreq;

	  //fill error counter: but remember simulation is not zero suppressing!
	  for( it_Sim_JEMHits  = vSim_JEMHits.begin(); it_Sim_JEMHits != vSim_JEMHits. end(); ++it_Sim_JEMHits )
	    {
	      mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_JEMHits)->crate()<<" Module "<<(*it_Sim_JEMHits)->module()<<endreq;
	      mLog<<MSG::VERBOSE<<"Sim: Hits "<<Help.Binary((*it_Sim_JEMHits)->JetHits(),24)<<endreq;
	      
	      if ((*it_Sim_JEMHits)->JetHits()!=0) {
	      m_h_SimBSMon_JEP->Fill(3,((*it_Sim_JEMHits)->crate()*19+(*it_Sim_JEMHits)->module()+1),1);	    
	      overview[(*it_Sim_JEMHits)->crate()] |= (1 << 2);
	      }
	    }
	}
     delete Sim_JEMHits;	
    }


  // =============================================================================================
  // ================= JetElements -> JEM Energy Sums ============================================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            JEMEnergySums(JEMDAQ.JetElements) = JEMDAQ.JEMEnergySums
  // ---------------------------------------------------------------------------------------------
  // retrieve JEMEtSums information from simulation and bytestream... 

      const JEMEtSumsCollection* BS_JEMEtSums;
      sc = m_storeGate->retrieve(BS_JEMEtSums, m_BS_JEMEtSumsLocation);
      if( (sc==StatusCode::FAILURE) ) {
	  mLog << MSG::INFO << "No JEMEtSums found in TES at "<< m_BS_JEMEtSumsLocation << endreq ;
	  return StatusCode::SUCCESS;
      }
	
     
  
  if (m_CompareWithSimulation ==1)
    {
      
      mLog << MSG::DEBUG << "==== JetElements -> JEM Energy Sums: Functionality ===="<< endreq ;
      
      
      JEMEtSumsCollection* Sim_JEMEtSums = 0;
      if (BS_JEMEtSums) {
        Sim_JEMEtSums = new JEMEtSumsCollection;
	mLog<<MSG::DEBUG<<"Simulate JEMEtSum from JetElements"<<endreq;
	m_etSumsTool->formJEMEtSums(jetElements, Sim_JEMEtSums);
	
      }
      
      
      // compare simulation-vector with bytestream-vector, erase corresponding entries
      // -> the remaining vectors contain faulty transmissions! 
      // use standard vector instead of datavector for transmissioncheck:
      // datavector erases not only the pointer  in the vector, but also the referenced object
      // -> segmentation fault!
      std::vector <LVL1::JEMEtSums*>  vBS_JEMEtSums;
      std::vector <LVL1::JEMEtSums*>  vSim_JEMEtSums;
      
      //fill vectors with simulation and bytestream data
      JEMEtSumsCollection::const_iterator it_JEMEtSums ;
      for( it_JEMEtSums  = BS_JEMEtSums->begin(); it_JEMEtSums != BS_JEMEtSums->end(); ++it_JEMEtSums )
	{
	  vBS_JEMEtSums.push_back(*it_JEMEtSums);	       
	}   
      for( it_JEMEtSums  = Sim_JEMEtSums->begin(); it_JEMEtSums != Sim_JEMEtSums->end(); ++it_JEMEtSums )
	{
	  vSim_JEMEtSums.push_back(*it_JEMEtSums);	       
	}   
      
      // Step over all cells and compare
      //int noMatchfound used to fill: "0" if a match is found; "1" if no match is found
      std::vector <LVL1::JEMEtSums*>::iterator it_BS_JEMEtSums;
      std::vector <LVL1::JEMEtSums*>::iterator it_Sim_JEMEtSums;
      int  foundModule;
      int noMatchfound1;
      int noMatchfound2;
      

      it_BS_JEMEtSums=vBS_JEMEtSums.begin();
      
      mLog<<MSG::DEBUG<<"JEMEtSums Calculation difference for"<<endreq;  
      while (it_BS_JEMEtSums != vBS_JEMEtSums.end())
	{
	  foundModule = 0;
	  it_Sim_JEMEtSums=vSim_JEMEtSums.begin();
	  
	  while ((foundModule==0)and(it_Sim_JEMEtSums != vSim_JEMEtSums.end()))
	    {	  	  
	      if (((*it_BS_JEMEtSums)->crate()==(*it_Sim_JEMEtSums)->crate())
		  and((*it_BS_JEMEtSums)->module()==(*it_Sim_JEMEtSums)->module()))
		{
		  foundModule=1;
		  noMatchfound1 = 0;
		  noMatchfound2 = 0;
		 
		      if (((*it_BS_JEMEtSums)->Ex()!=(*it_Sim_JEMEtSums)->Ex())
			  or ((*it_BS_JEMEtSums)->Ey()!=(*it_Sim_JEMEtSums)->Ey()))
			  {
			  mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_JEMEtSums)->crate()<<" Module "<<(*it_BS_JEMEtSums)->module()<<endreq;
			  mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_JEMEtSums)->Ex()<<endreq;
			  mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_JEMEtSums)->Ey()<<endreq;
			  mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_JEMEtSums)->Et()<<endreq;
			  			  
			  mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_JEMEtSums)->Ex()<<endreq;
			  mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_JEMEtSums)->Ey()<<endreq;
			  mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_JEMEtSums)->Et()<<endreq;
			  			  
			  noMatchfound1 = 1;
			}

		      m_h_SimBSMon_JEP->Fill(1,(*it_BS_JEMEtSums)->crate()*19 + ((*it_BS_JEMEtSums)->module()+1),noMatchfound1);
		      overview[(*it_BS_JEMEtSums)->crate()] |= (noMatchfound1 << 3);

		      if ((*it_BS_JEMEtSums)->Et()!=(*it_Sim_JEMEtSums)->Et())
			{
			  mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_JEMEtSums)->crate()<<" Module "<<(*it_BS_JEMEtSums)->module()<<endreq;
			  mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_JEMEtSums)->Et()<<endreq;
			  
		      	  mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_JEMEtSums)->Et()<<endreq;
			  
			  noMatchfound2 = 1;
			}	 

   		  
		     m_h_SimBSMon_JEP->Fill(2,(*it_BS_JEMEtSums)->crate()*19 + ((*it_BS_JEMEtSums)->module()+1),noMatchfound2);
		     overview[(*it_BS_JEMEtSums)->crate()] |= (noMatchfound2 << 3);
		  	  
		  it_BS_JEMEtSums=vBS_JEMEtSums.erase(it_BS_JEMEtSums);
		  it_Sim_JEMEtSums=vSim_JEMEtSums.erase(it_Sim_JEMEtSums);
		}
	      else ++it_Sim_JEMEtSums;
	    }
	  if (foundModule==0)++it_BS_JEMEtSums;
	}
      
      if (vBS_JEMEtSums.size()!=0)
	{
	  mLog<<MSG::DEBUG<<"JEMEtSums: additional BS data for"<<endreq;
	  
	  //fill errorcounter
	  for( it_BS_JEMEtSums  = vBS_JEMEtSums.begin(); it_BS_JEMEtSums != vBS_JEMEtSums. end(); ++it_BS_JEMEtSums )
	    {
	      mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_JEMEtSums)->crate()<<" Module "<<(*it_BS_JEMEtSums)->module()<<endreq;
	      mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_JEMEtSums)->Ex()<<endreq;
	      mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_JEMEtSums)->Ey()<<endreq;
	      mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_JEMEtSums)->Et()<<endreq;
	      
	      //fill both error bins as additional data was found for all energies
	      m_h_SimBSMon_JEP->Fill(1,(*it_BS_JEMEtSums)->crate()*19 + ((*it_BS_JEMEtSums)->module()+1),1);	  
	      m_h_SimBSMon_JEP->Fill(2,(*it_BS_JEMEtSums)->crate()*19 + ((*it_BS_JEMEtSums)->module()+1),1); 
             overview[(*it_BS_JEMEtSums)->crate()] |= (1 << 3);
	    }
	}
      
      if (vSim_JEMEtSums.size()!=0)
	{
	   mLog<<MSG::DEBUG<<"JEMEtSums: additional Sim data for"<<endreq;
	   mLog<<MSG::DEBUG<<"fill error counter but remember simulation is not zero suppressing!"<<endreq;

	  //fill errorcounter
	  for( it_Sim_JEMEtSums  = vSim_JEMEtSums.begin(); it_Sim_JEMEtSums != vSim_JEMEtSums. end(); ++it_Sim_JEMEtSums )
	    {
	      mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_JEMEtSums)->crate()<<" Module "<<(*it_Sim_JEMEtSums)->module()<<endreq;
	      mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_JEMEtSums)->Ex()<<endreq;
	      mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_JEMEtSums)->Ey()<<endreq;
	      mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_JEMEtSums)->Et()<<endreq;

	      //fill both error bins as additional data was found for all energies (only if additional data !=0)
	      if(((*it_Sim_JEMEtSums))->Et()!=0 or ((*it_Sim_JEMEtSums))->Ey()!=0 or ((*it_Sim_JEMEtSums)->Ex())) 
		{
	      m_h_SimBSMon_JEP->Fill(1,(*it_Sim_JEMEtSums)->crate()*19 + ((*it_Sim_JEMEtSums)->module()+1),1);
	      m_h_SimBSMon_JEP->Fill(2,(*it_Sim_JEMEtSums)->crate()*19 + ((*it_Sim_JEMEtSums)->module()+1),1);
              overview[(*it_Sim_JEMEtSums)->crate()] |= (1 << 3);
		}
	    }
	}
      delete Sim_JEMEtSums;
    }


  // =============================================================================================
  // ================= JEM Jet Hits -> (Crate/System) CMM Jet Hits ===============================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Transmission
  //            JEMDAQ.JEMJetHits = (Crate/System)CMMJetDAQ.JEMJetHits
  // ---------------------------------------------------------------------------------------------

  mLog << MSG::DEBUG << "==== JEM Jet Hits -> (Crate/System) CMM Jet Hits: Transmission ===="<< endreq ;

  // retrieve CMM Jet Hits from Storegate
  const CMMJetHitsCollection* CMMJetHits;
  sc = m_storeGate->retrieve(CMMJetHits, m_BS_CMMJetHitsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMM JetHits found in TES at "  << m_BS_CMMJetHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }  
  CMMJetHitsCollection::const_iterator it_CMMJetHits;
          
  //retrieve JEMHits from storegate for comparison with transmitted data stored in CMMJetHits
  const JEMHitsCollection* JEMHits;
  JEMHitsCollection::const_iterator it_JEMHits ;
  sc = m_storeGate->retrieve(JEMHits, m_BS_JEMHitsLocation); 
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO<< "No JEMHits found in TES at "<< m_BS_JEMHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  // compare JEMHits-vector with CMMJEMHits-vector, erase corresponding entries
  // -> the remaining vectors contain faulty transmissions! 
  // use standard vector instead of datavector for transmissioncheck:
  // datavector erases not only the pointer  in the vector, but also the referenced object
  // -> segmentation fault!
  std::vector <LVL1::CMMJetHits*>  vCMMJEMHits;
  std::vector <LVL1::JEMHits*>  vJEMHits;
  int CrateCMMMainHits=0;
  int CrateCMMFwdHits=0;
  int SystemRemoteCMMMainHits=0;
  int SystemRemoteCMMFwdHits=0;
 
  
  // put CMM input data (JEMHits) into CMMJEMHits vector
  for( it_CMMJetHits  = CMMJetHits ->begin(); it_CMMJetHits != CMMJetHits -> end(); ++it_CMMJetHits )
    {
      // JEM information for transmission check JEMs -> CMMs 
      if ((*it_CMMJetHits)->dataID()<16)
	{
	  vCMMJEMHits.push_back(*it_CMMJetHits );
	}
      
      // CMM information for transmission check crate CMM -> system CMM
      // local main hits of crate CMM
      if (((*it_CMMJetHits)->dataID()==17)and((*it_CMMJetHits)->crate()==0)) CrateCMMMainHits=(*it_CMMJetHits)->Hits();
      // local fwd hits of crate CMM
      if (((*it_CMMJetHits)->dataID()==20)and((*it_CMMJetHits)->crate()==0)) CrateCMMFwdHits=(*it_CMMJetHits)->Hits();
      // remote main hits of system CMM
      if (((*it_CMMJetHits)->dataID()==16)and((*it_CMMJetHits)->crate()==1)) SystemRemoteCMMMainHits=(*it_CMMJetHits)->Hits();
      // remote fwd hits of crate CMM
      if (((*it_CMMJetHits)->dataID()==19)and((*it_CMMJetHits)->crate()==1)) SystemRemoteCMMFwdHits=(*it_CMMJetHits)->Hits();
    }
  for( it_JEMHits  =  JEMHits->begin(); it_JEMHits <  JEMHits-> end(); ++it_JEMHits )
    {
      vJEMHits.push_back(*it_JEMHits);
    }
  
  bool found;
  std::vector <LVL1::CMMJetHits*>::iterator it_vCMMJEMHits;
  std::vector <LVL1::JEMHits*>::iterator it_vJEMHits;
  
  int noMatchfound;
  
  //---------------------------------- backplane transmission JEMs -> CMM -----------------------------
  it_vCMMJEMHits=vCMMJEMHits.begin();
  // step through both vectors and compare...
  mLog<<MSG::DEBUG<<"JEMHits Transmission error JEM -> CMM"<<endreq;
  while (it_vCMMJEMHits != vCMMJEMHits.end())
    {
      found = 0;
      noMatchfound=0;
      it_vJEMHits=vJEMHits.begin();
      
      while ((found==0)and(it_vJEMHits != vJEMHits.end()))
	{
	  if (((*it_vCMMJEMHits)->crate()==(*it_vJEMHits)->crate())
	      
	      and((*it_vCMMJEMHits)->dataID()==(*it_vJEMHits)->module()))
	      {
	      if ((*it_vCMMJEMHits)->Hits()!=(*it_vJEMHits)->JetHits())
		{
		  mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMHits)->crate()<<" Module "<<(*it_vJEMHits)->module()<<endreq;
		  mLog<<MSG::VERBOSE<<"JEM Hit "<<Help.Binary((*it_vJEMHits)->JetHits(),24)<<endreq;
		  mLog<<MSG::VERBOSE<<"CMM Hit "<<Help.Binary((*it_vCMMJEMHits)->Hits(),24)<<endreq;
		  
		  noMatchfound=1;
		}
	      // if CMMJEMHits and JEMHits didn't match, a "1" is entered at the specific module;
	      // if a match is found, a "0" is entered (just to see that the histogram is filled
	      // since there hopefully are only zeros in there)
	      m_h_TransCheck_JEP->Fill(1,(*it_vJEMHits)->crate()*19+(*it_vJEMHits)->module()+1,noMatchfound);
	      overview[(*it_vJEMHits)->crate()] |= (noMatchfound << 4);
	      
	      it_vCMMJEMHits=vCMMJEMHits.erase(it_vCMMJEMHits);
	      it_vJEMHits=vJEMHits.erase(it_vJEMHits);
	      found=1;
	    }
	  else ++it_vJEMHits;
	}// end while (not found and not JEMHits.end)
      
      // only step further in CMMJEMHits if no corresponding entry has been found in JEMHits
      if (found==0)++it_vCMMJEMHits;
    }  // end while (not CMMJEMHits.end)
  
  // if the CMMJEMHits vector isn't empty, fill the error counter!
  if (vCMMJEMHits.size()!=0)
    {
      mLog<<MSG::DEBUG<<vCMMJEMHits.size()<<"JEMHits Transmission: additional CMM information"<<endreq;
      for( it_vCMMJEMHits  = vCMMJEMHits.begin(); it_vCMMJEMHits != vCMMJEMHits. end(); ++it_vCMMJEMHits )
	{
	 
	  mLog<<MSG::DEBUG<<"Crate "<<(*it_vCMMJEMHits)->crate()<<" Module "<<((*it_vCMMJEMHits)->dataID())<<endreq;
	  mLog<<MSG::VERBOSE<<"Hit "<<Help.Binary((*it_vCMMJEMHits)->Hits(),24)<<endreq;
	  
	  
	  if ((*it_vCMMJEMHits)->Hits()>0)
	    {
	      m_h_TransCheck_JEP->Fill(1,(*it_vCMMJEMHits)->crate()*19+(*it_vCMMJEMHits)->dataID()+1,1);
	      overview[(*it_vCMMJEMHits)->crate()] |= (1 << 4);
	    }
	}
    }
  
  // if the JEMHits vector isn't empty, fill the error counter!
  if (vJEMHits.size()!=0)
    {
      mLog<<MSG::DEBUG<<vJEMHits.size()<<" JEMHits Transmission: addtional JEM information"<<endreq;
      for( it_vJEMHits  =  vJEMHits.begin(); it_vJEMHits != vJEMHits. end(); ++it_vJEMHits )
	{
	  mLog<<MSG::DEBUG<<"Crate "<<(*it_vJEMHits)->crate()<<" Module "<<(*it_vJEMHits)->module()<<endreq;
	  mLog<<MSG::VERBOSE<<"Hit "<<Help.Binary((*it_vJEMHits)->JetHits(),24)<<endreq;
	  
	  if ((*it_vJEMHits)->JetHits()>0) 
	    {
	      m_h_TransCheck_JEP->Fill(1,(*it_vJEMHits)->crate()*19+(*it_vJEMHits)->module()+1,1);
	      overview[(*it_vJEMHits)->crate()] |= (1 << 4);
	    }
	}
    }



  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            CMMJetHits(JEMDAQ.JEMJetHits) = (Crate/System)CMMJetDAQ.CrateCMMJetHits
  // ---------------------------------------------------------------------------------------------
  const CMMJetHitsCollection* BS_CMMJetHits;
  sc = m_storeGate->retrieve(BS_CMMJetHits, m_BS_CMMJetHitsLocation);
  if( (sc==StatusCode::FAILURE) ) {
    mLog << MSG::INFO << "No CMMJetHits found in TES at "<< m_BS_CMMJetHitsLocation << endreq ;
    return StatusCode::SUCCESS;
  }
  
  if (m_CompareWithSimulation ==1)
    {
      
      mLog << MSG::DEBUG << "==== JEM Jet Hits -> (Crate/System) CMM Jet Hits: Functionality ===="<< endreq ;
      
      CMMJetHitsCollection* Sim_CMMJetHits = 0;
      if (BS_JEMHits) {
        Sim_CMMJetHits = new CMMJetHitsCollection;
	mLog<<MSG::DEBUG<<"Simulate CMMJetHits from JEMHits"<<endreq;
	m_jepHitsTool->formCMMJetHits(BS_JEMHits, Sim_CMMJetHits);
      }
	
    std::vector <LVL1::CMMJetHits*>  vBS_CMMJetHits;
    std::vector <LVL1::CMMJetHits*>  vSim_CMMJetHits;
      //use standard vector instead of datavector for transmissioncheck:
      //datavector erases not only the pointer  in the vector, but also the referenced object
      //-> segmentation fault!
      
      //fill vectors
      //CMMJetHitsCollection::const_iterator it_CMMJetHits ;
      
      for( it_CMMJetHits  = BS_CMMJetHits->begin(); it_CMMJetHits != BS_CMMJetHits->end(); ++it_CMMJetHits )
	{
	  if ((*it_CMMJetHits)->dataID()>15)
	    {
	      vBS_CMMJetHits.push_back(*it_CMMJetHits);	       
	    }
	}   
      for( it_CMMJetHits  = Sim_CMMJetHits->begin(); it_CMMJetHits != Sim_CMMJetHits->end(); ++it_CMMJetHits )
	{
	  if ((*it_CMMJetHits)->dataID()>15)
	    {
	      vSim_CMMJetHits.push_back(*it_CMMJetHits);	
	    }       
	}   
      
      // Step over all cells and compare
      std::vector <LVL1::CMMJetHits*>::iterator it_BS_CMMJetHits;
      std::vector <LVL1::CMMJetHits*>::iterator it_Sim_CMMJetHits;
      int  foundModule;
      int noMatchfound;
      
      it_BS_CMMJetHits=vBS_CMMJetHits.begin();
      
      mLog<<MSG::DEBUG<<"CMMJetHits Calculation difference for"<<endreq;
      while (it_BS_CMMJetHits != vBS_CMMJetHits.end())
	{
	  foundModule = 0;
	  noMatchfound=0;
	  it_Sim_CMMJetHits=vSim_CMMJetHits.begin();
	  
	  while ((foundModule==0)and(it_Sim_CMMJetHits != vSim_CMMJetHits.end()))
	    {
	      if (((*it_BS_CMMJetHits)->crate()==(*it_Sim_CMMJetHits)->crate())
		  and((*it_BS_CMMJetHits)->dataID()==(*it_Sim_CMMJetHits)->dataID()))
		{
		  foundModule=1;
		  
		  if ((*it_BS_CMMJetHits)->Hits()!=(*it_Sim_CMMJetHits)->Hits())
		    {
		      mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_CMMJetHits)->crate()<<" DataId "<<(*it_BS_CMMJetHits)->dataID()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Hits"<<Help.Binary((*it_BS_CMMJetHits)->Hits(),24)<<endreq;
		      mLog<<MSG::VERBOSE<<"Sim: Hits"<<Help.Binary((*it_Sim_CMMJetHits)->Hits(),24)<<endreq;
		      
		      noMatchfound=1;
		    }
		  m_h_SimBSMon_JEP->Fill(3,((*it_BS_CMMJetHits)->crate()*19+16+1),noMatchfound);
		  overview[(*it_BS_CMMJetHits)->crate()] |= (noMatchfound << 5);
		  
		  it_BS_CMMJetHits=vBS_CMMJetHits.erase(it_BS_CMMJetHits);
		  it_Sim_CMMJetHits=vSim_CMMJetHits.erase(it_Sim_CMMJetHits);
		}
	      else ++it_Sim_CMMJetHits;
	    }
	  if (foundModule==0)++it_BS_CMMJetHits;
	}
      
      if (vBS_CMMJetHits.size()!=0)
	{
	  mLog<<MSG::DEBUG<<"CMMJetHits: additional BS data for"<<endreq;
	  
	  //fill errorcounter
	  for( it_BS_CMMJetHits  = vBS_CMMJetHits.begin(); it_BS_CMMJetHits != vBS_CMMJetHits. end(); ++it_BS_CMMJetHits )
	    {
	      mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_CMMJetHits)->crate()<<" DataId "<<(*it_BS_CMMJetHits)->dataID()<<endreq;
	      mLog<<MSG::VERBOSE<<"BS: Hits"<<Help.Binary((*it_BS_CMMJetHits)->Hits(),24)<<endreq;
	      
	      m_h_SimBSMon_JEP->Fill(3,((*it_BS_CMMJetHits)->crate()*19+16+1),1);
              overview[(*it_BS_CMMJetHits)->crate()] |= (1 << 5);
	    }
	}
      
      if (vSim_CMMJetHits.size()!=0)
	{
	  mLog<<MSG::DEBUG<<"CMMJetHits: additional Sim data "<<endreq;
	  mLog<<MSG::DEBUG<<"fill error counter but remember simulation is not zero suppressing!"<<endreq;

	  //fill errorcounter: simulation isn't zero suppressing
	  for( it_Sim_CMMJetHits  = vSim_CMMJetHits.begin(); it_Sim_CMMJetHits != vSim_CMMJetHits. end(); ++it_Sim_CMMJetHits )
	    {
	      mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_CMMJetHits)->crate()<<" DataId "<<(*it_Sim_CMMJetHits)->dataID()<<endreq;
	      mLog<<MSG::VERBOSE<<"Sim: Hits"<<Help.Binary((*it_Sim_CMMJetHits)->Hits(),24)<<endreq;
	      
	      if((*it_Sim_CMMJetHits)->Hits()!=0) 
		{
		  m_h_SimBSMon_JEP->Fill(3,((*it_Sim_CMMJetHits)->crate()*19+16+1),1);
                  overview[(*it_Sim_CMMJetHits)->crate()] |= (1 << 5);
		}
	    }
	}
      delete Sim_CMMJetHits;
    }
	
    

  
  // =============================================================================================
  // ================= JEMEnergySums -> (Crate/System) CMM Energy Sums ===========================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Transmission
  //            JEMDAQ.JEMEnergySums = (Crate/System)CMMEnergySums.JEMEnergySums
  // ---------------------------------------------------------------------------------------------

  mLog << MSG::DEBUG << "==== JEMEnergySums -> (Crate/System) CMM Energy Sums: Transmission ===="<< endreq ;

  // retrieve CMM Et Sums from Storegate
  const CMMEtSumsCollection* CMMEtSums;
  sc = m_storeGate->retrieve(CMMEtSums, m_BS_CMMEtSumsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMMEtSums found in TES at " << m_BS_CMMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  CMMEtSumsCollection::const_iterator it_CMMEtSums ;
  
  const JEMEtSumsCollection* JEMEtSums;
  JEMEtSumsCollection::const_iterator it_JEMEtSums ;
  sc = m_storeGate->retrieve(JEMEtSums, m_BS_JEMEtSumsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMEtSums found in TES at "<< m_BS_JEMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
 
  
  // compare JEMHits-vector with CMMJEMHits-vector, erase corresponding entries
  // -> the remaining vectors contain faulty transmissions! 
  // use standard vector instead of datavector for transmissioncheck:
  // datavector erases not only the pointer  in the vector, but also the referenced object
  // -> segmentation fault!
  std::vector <LVL1::CMMEtSums*>  vCMMJEMEtSums;
  std::vector <LVL1::JEMEtSums*>  vJEMEtSums;
  int CrateEx=0, CrateEy=0, CrateEt=0;
  int SystemEx=0, SystemEy=0, SystemEt=0;
  noMatchfound=0;
  for( it_CMMEtSums  = CMMEtSums ->begin(); it_CMMEtSums != CMMEtSums -> end(); ++it_CMMEtSums )
    {	  
      if ((*it_CMMEtSums)->dataID()<16)
	{
	  vCMMJEMEtSums.push_back(*it_CMMEtSums );
	}
      // CMM information for transmission check crate CMM -> system CMM
      // local energies of crate CMM
      if (((*it_CMMEtSums)->dataID()==17)and((*it_CMMEtSums)->crate()==0)) 
	{
	  CrateEx=(*it_CMMEtSums)->Ex();
	  CrateEy=(*it_CMMEtSums)->Ey();
	  CrateEt=(*it_CMMEtSums)->Et();
	}
      // remote energies of system CMM
      if (((*it_CMMEtSums)->dataID()==16)and((*it_CMMEtSums)->crate()==1)) 
	{
	  SystemEx=(*it_CMMEtSums)->Ex();
	  SystemEy=(*it_CMMEtSums)->Ey();
	  SystemEt=(*it_CMMEtSums)->Et();
	}
    }
  for( it_JEMEtSums  =  JEMEtSums->begin(); it_JEMEtSums != JEMEtSums-> end(); ++it_JEMEtSums )
    {
      vJEMEtSums.push_back(*it_JEMEtSums);
    }
  
  std::vector <LVL1::CMMEtSums*>::iterator it_vCMMJEMEtSums;
  std::vector <LVL1::JEMEtSums*>::iterator it_vJEMEtSums;
  
  //---------------------------------- backplane transmission JEMs -> crate CMM -----------------------------
  it_vCMMJEMEtSums=vCMMJEMEtSums.begin();
  mLog<<MSG::DEBUG<<"CMMEtSums Transmission errors JEM -> CMM"<<endreq;
  while (it_vCMMJEMEtSums != vCMMJEMEtSums.end())
    {
      bool found = 0;
      noMatchfound = 0;
      it_vJEMEtSums=vJEMEtSums.begin();
      
      while ((found==0)and(it_vJEMEtSums != vJEMEtSums.end()))
	{
	  if (((*it_vCMMJEMEtSums)->crate()==(*it_vJEMEtSums)->crate())
	      and((*it_vCMMJEMEtSums)->dataID()==(*it_vJEMEtSums)->module()))
	    {
	      if (((*it_vCMMJEMEtSums)->Ex()!=(*it_vJEMEtSums)->Ex())
		  or ((*it_vCMMJEMEtSums)->Ey()!=(*it_vJEMEtSums)->Ey())
		  or ((*it_vCMMJEMEtSums)->Et()!=(*it_vJEMEtSums)->Et()))
		{
		  mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMEtSums)->crate()<<" Module "<<(*it_vJEMEtSums)->module()<<endreq;
		  mLog<<MSG::VERBOSE<<"JEM Ex "<<(*it_vJEMEtSums)->Ex()<<endreq;
		  mLog<<MSG::VERBOSE<<"JEM Ey "<<(*it_vJEMEtSums)->Ey()<<endreq;
		  mLog<<MSG::VERBOSE<<"JEM Et "<<(*it_vJEMEtSums)->Et()<<endreq;
		  
		  mLog<<MSG::VERBOSE<<"CMM Ex "<<(*it_vCMMJEMEtSums)->Ex()<<endreq;
		  mLog<<MSG::VERBOSE<<"CMM Ey "<<(*it_vCMMJEMEtSums)->Ey()<<endreq;
		  mLog<<MSG::VERBOSE<<"CMM Et "<<(*it_vCMMJEMEtSums)->Et()<<endreq;
		  
		  noMatchfound=1;
		}
	      // if CMMJEMEtSums and JEMEtSums didn't match, a "1" is entered at the specific module;
	      // if a match is found, a "0" is entered (just to see that the histogram is filled
	      // since there hopefully are only zeros in there)
	      
	      m_h_TransCheck_JEP->Fill(2,((*it_vJEMEtSums)->crate()*19+(*it_vJEMEtSums)->module() + 1),noMatchfound);
	      overview[(*it_vJEMEtSums)->crate()] |= (noMatchfound << 6);
	      
	      it_vCMMJEMEtSums=vCMMJEMEtSums.erase(it_vCMMJEMEtSums);
	      it_vJEMEtSums=vJEMEtSums.erase(it_vJEMEtSums);
	      found=1;
	    }
	  else ++it_vJEMEtSums;
	}
      if (found==0)++it_vCMMJEMEtSums;
    }
  
  // if the CMMJEMEtSums vector isn't empty, fill the error counter!
  if (vCMMJEMEtSums.size()!=0)
    {
      mLog<<MSG::DEBUG<<vCMMJEMEtSums.size()<<"CMMEtSums Transmission: additional CMM information"<<endreq;
      for( it_vCMMJEMEtSums  = vCMMJEMEtSums.begin(); it_vCMMJEMEtSums != vCMMJEMEtSums. end(); ++it_vCMMJEMEtSums )
	{
	  mLog<<MSG::DEBUG<<" Crate "<<(*it_vCMMJEMEtSums)->crate()<<" Module "<<(*it_vCMMJEMEtSums)->dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Ex "<<(*it_vCMMJEMEtSums)->Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Ey "<<(*it_vCMMJEMEtSums)->Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Et "<<(*it_vCMMJEMEtSums)->Et()<<endreq;
	  
	  if ((*it_vCMMJEMEtSums)->Ex()!=0 or (*it_vCMMJEMEtSums)->Ey()!=0 or (*it_vCMMJEMEtSums)->Et()!=0)
	    {
	      m_h_TransCheck_JEP->Fill(2,(*it_vCMMJEMEtSums)->crate()*19+(*it_vCMMJEMEtSums)->dataID() + 1,1);
	      overview[(*it_vCMMJEMEtSums)->crate()] |= (1 << 6);
	    }
  	}
    }
  
  // if the JEMEtSums vector isn't empty, fill the error counter!
  if (vJEMEtSums.size()!=0)
    {
      mLog<<MSG::DEBUG<<vJEMEtSums.size()<<"CMMEtSums Transmission: addtional JEM information"<<endreq;
      for( it_vJEMEtSums  =  vJEMEtSums.begin(); it_vJEMEtSums != vJEMEtSums. end(); ++it_vJEMEtSums )
	{
	  mLog<<MSG::DEBUG<<"Crate "<<(*it_vJEMEtSums)->crate()<<" Module "<<(*it_vJEMEtSums)->module()<<endreq;
	  mLog<<MSG::VERBOSE<<"JEM Ex "<<(*it_vJEMEtSums)->Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"JEM Ey "<<(*it_vJEMEtSums)->Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"JEM Et "<<(*it_vJEMEtSums)->Et()<<endreq;
	  
	  m_h_TransCheck_JEP->Fill(2,(*it_vJEMEtSums)->crate()*19+(*it_vJEMEtSums)->module() + 1,1);
	  overview[(*it_vJEMEtSums)->crate()] |= (1 << 6);
	  	}
    }
  

  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            CMMEnergySums(JEMDAQ.JEMEnergySums) = (Crate/System)CMMEnergyDAQ.CrateCMMEnergySums
  // ---------------------------------------------------------------------------------------------
    
    const CMMEtSumsCollection* BS_CMMEtSums;
      sc = m_storeGate->retrieve(BS_CMMEtSums, m_BS_CMMEtSumsLocation);
      if( (sc==StatusCode::FAILURE) ) 
	{
	  mLog << MSG::INFO << "No CMMEtSums found in TES at "<< m_BS_CMMEtSumsLocation << endreq ;
	  return StatusCode::SUCCESS;
	}
    
    if (m_CompareWithSimulation ==1)
    {

      mLog << MSG::DEBUG << "==== JEMEnergySums -> (Crate/System) CMM Energy Sums: Functionality ===="<< endreq ;
        
      
      CMMEtSumsCollection* Sim_CMMEtSums = 0;
      if (BS_JEMEtSums) {
        Sim_CMMEtSums = new CMMEtSumsCollection;
	mLog<<MSG::DEBUG<<"Simulate CMMEtSums from JEMEtSum"<<endreq;
	m_etSumsTool->formCMMEtSums(BS_JEMEtSums, Sim_CMMEtSums);
      }        
      
      std::vector <LVL1::CMMEtSums*>  vBS_CMMEtSums;
      std::vector <LVL1::CMMEtSums*>  vSim_CMMEtSums;
      //use standard vector instead of datavector for transmissioncheck:
      //datavector erases not only the pointer  in the vector, but also the referenced object
      //-> segmentation fault!
      
      //fill vectors
      //CMMEtSumsCollection::const_iterator it_CMMEtSums ;
      
      for( it_CMMEtSums  = BS_CMMEtSums->begin(); it_CMMEtSums != BS_CMMEtSums->end(); ++it_CMMEtSums )
	{
	  if ((*it_CMMEtSums)->dataID()>15)
	    {
	      vBS_CMMEtSums.push_back(*it_CMMEtSums);	       
	    }
	}   
      for( it_CMMEtSums  = Sim_CMMEtSums->begin(); it_CMMEtSums != Sim_CMMEtSums->end(); ++it_CMMEtSums )
	{
	  if ((*it_CMMEtSums)->dataID()>15)
	    {
	      vSim_CMMEtSums.push_back(*it_CMMEtSums);	
	    }       
	}   
      
      // Step over all cells and compare
      std::vector <LVL1::CMMEtSums*>::iterator it_BS_CMMEtSums;
      std::vector <LVL1::CMMEtSums*>::iterator it_Sim_CMMEtSums;
      int  foundModule;
      int noMatchfound1;
      int noMatchfound2;

      it_BS_CMMEtSums=vBS_CMMEtSums.begin();
      
      mLog<<MSG::DEBUG<<"CMMEtSums Calculation difference for"<<endreq;
      while (it_BS_CMMEtSums != vBS_CMMEtSums.end())
	{
	  foundModule = 0;
	  noMatchfound1=0;
	  noMatchfound2=0;

	  it_Sim_CMMEtSums=vSim_CMMEtSums.begin();
	  
	  while ((foundModule==0)and(it_Sim_CMMEtSums != vSim_CMMEtSums.end()))
	    {
	      if (((*it_BS_CMMEtSums)->crate()==(*it_Sim_CMMEtSums)->crate())
		  and((*it_BS_CMMEtSums)->dataID()==(*it_Sim_CMMEtSums)->dataID()))
		{
		  foundModule=1;
		  
		  if (((*it_BS_CMMEtSums)->Ex()!=(*it_Sim_CMMEtSums)->Ex())
		      or ((*it_BS_CMMEtSums)->Ey()!=(*it_Sim_CMMEtSums)->Ey()))
		    {
		      mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_CMMEtSums)->crate()<<" DataId "<<(*it_BS_CMMEtSums)->dataID()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_CMMEtSums)->Ex()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_CMMEtSums)->Ey()<<endreq;
		      		      
		      mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_CMMEtSums)->Ex()<<endreq;
		      mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_CMMEtSums)->Ey()<<endreq;
		      
		      
		      noMatchfound1=1;    
		    }
		  m_h_SimBSMon_JEP->Fill(1,(*it_BS_CMMEtSums)->crate()*19 + 16 + 1,noMatchfound1);
		  overview[(*it_BS_CMMEtSums)->crate()] |= (noMatchfound1 << 7);

		  if((*it_BS_CMMEtSums)->Et()!=(*it_Sim_CMMEtSums)->Et())
		    {
		      mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_CMMEtSums)->crate()<<" DataId "<<(*it_BS_CMMEtSums)->dataID()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_CMMEtSums)->Et()<<endreq;
		      mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_CMMEtSums)->Et()<<endreq;

		      noMatchfound2=1;
		    }    
		  m_h_SimBSMon_JEP->Fill(2,(*it_BS_CMMEtSums)->crate()*19 + 16 + 1,noMatchfound2);
		  overview[(*it_BS_CMMEtSums)->crate()] |= (noMatchfound2 << 7);
		  

		  it_BS_CMMEtSums=vBS_CMMEtSums.erase(it_BS_CMMEtSums);
		  it_Sim_CMMEtSums=vSim_CMMEtSums.erase(it_Sim_CMMEtSums);
		}
	      else ++it_Sim_CMMEtSums;
	    }
	  if (foundModule==0)++it_BS_CMMEtSums;
	}
      
      if (vBS_CMMEtSums.size()!=0)
	{
	  mLog<<MSG::DEBUG<<"CMMEtSums: additional BS data for"<<endreq;
	  
	  //fill errorcounter
	  for( it_BS_CMMEtSums  = vBS_CMMEtSums.begin(); it_BS_CMMEtSums != vBS_CMMEtSums. end(); ++it_BS_CMMEtSums )
	    {
	      mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_CMMEtSums)->crate()<<" DataId "<<(*it_BS_CMMEtSums)->dataID()<<endreq;
	      mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_CMMEtSums)->Ex()<<endreq;
	      mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_CMMEtSums)->Ey()<<endreq;
	      mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_CMMEtSums)->Et()<<endreq;
	      //fill both error bins for such a case, surpress zero data
	      if(((*it_BS_CMMEtSums)->Ex()!=0) or ((*it_BS_CMMEtSums)->Ey()!=0) or ((*it_BS_CMMEtSums)->Et()!=0))
		{
		 m_h_SimBSMon_JEP->Fill(1,(*it_BS_CMMEtSums)->crate()*19 + 16 + 1,1);
	         m_h_SimBSMon_JEP->Fill(2,(*it_BS_CMMEtSums)->crate()*19 + 16 + 1,1);
		 overview[(*it_BS_CMMEtSums)->crate()] |= (1 << 7);
	    	}
	}
      
      if (vSim_CMMEtSums.size()!=0)
	{
	  mLog<<MSG::DEBUG<<"CMMEtSums: additional Sim data for"<<endreq;
	  mLog<<MSG::DEBUG<<"fill error counter but remember simulation is not zero suppressing!"<<endreq;
	  
	  //fill errorcounter
	  for( it_Sim_CMMEtSums  = vSim_CMMEtSums.begin(); it_Sim_CMMEtSums != vSim_CMMEtSums. end(); ++it_Sim_CMMEtSums )
	    {
	      mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_CMMEtSums)->crate()<<" DataId "<<(*it_Sim_CMMEtSums)->dataID()<<endreq;
	      mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_CMMEtSums)->Ex()<<endreq;
	      mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_CMMEtSums)->Ey()<<endreq;
	      mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_CMMEtSums)->Et()<<endreq;

	     //fill both error bins for such a case 
	      if(((*it_Sim_CMMEtSums)->Ex()!=0) or ((*it_Sim_CMMEtSums)->Ey()!=0) or ((*it_Sim_CMMEtSums)->Et()!=0))
		{
	      
		  m_h_SimBSMon_JEP->Fill(1,(*it_Sim_CMMEtSums)->crate()*19 + 16 + 1,1);
		  m_h_SimBSMon_JEP->Fill(2,(*it_Sim_CMMEtSums)->crate()*19 + 16 + 1,1);
		  overview[(*it_Sim_CMMEtSums)->crate()] |= (1 << 7);
		}
	     }
         }
      }
     delete Sim_CMMEtSums; 
    }
      
  // =============================================================================================
  // ================= Crate CMM Jet Hits -> System CMM Jet Hits =================================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Transmission
  //            CrateCMMJetDAQ.CrateCMMJetHits = SystemCMMJetDAQ.RemoteCMMJetHits
  // ---------------------------------------------------------------------------------------------

  mLog << MSG::DEBUG << "==== Crate CMM Jet Hits -> System CMM Jet Hits: Transmission ===="<< endreq ;

  noMatchfound=0;
  
  if ((CrateCMMMainHits!=SystemRemoteCMMMainHits) or (CrateCMMFwdHits!=SystemRemoteCMMFwdHits))
    {
      noMatchfound=1;
      mLog<<MSG::DEBUG<<"CMMJetHits: Transmission error between crate and system CMM"<<endreq;
    }
  m_h_TransCheck_JEP->Fill(1,17,noMatchfound);
  overview[0] |= (noMatchfound << 8);


  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            TotalCMMJetHits(CrateCMMJetDAQ.CrateCMMJetHits+ SystemCMMJetDAQ.CrateCMMJetHits) = SystemCMMJetDAQ.SystemCMMJetHits
  // ---------------------------------------------------------------------------------------------
  //s. Check	Functionality
  //            CMMJetHits(JEMDAQ.JEMJetHits) = (Crate/System)CMMJetDAQ.CrateCMMJetHits




  // =============================================================================================
  // ================= Crate CMM Energy Sums -> System CMM Energy Sums ===========================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Transmission
  //            CrateCMMEnergyDAQ.CrateCMMEnergySums = SystemCMMEnergyDAQ.RemoteCMMEnergySums
  // ---------------------------------------------------------------------------------------------

  mLog << MSG::DEBUG << "==== Crate CMM Energy Sums -> System CMM Energy Sums: Transmission ===="<< endreq ;

  noMatchfound=0;
  if ((CrateEx!=SystemEx)or (CrateEy!=SystemEy)or (CrateEt!=SystemEt))
    {
      noMatchfound=1;
      mLog<<MSG::DEBUG<<"CMMEnergySums: Transmission error between crate and system CMM"<<endreq;
    }

  m_h_TransCheck_JEP->Fill(2,17,noMatchfound);
  overview[0] |= (noMatchfound << 9);


  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            TotalCMMEnergySums(CrateCMMEnergyDAQ.CrateCMMEnergySums+ SystemCMMEnergyDAQ.CrateCMMEnergySums) = SystemCMMEnergyDAQ.SystemCMMEnergySums
  // ---------------------------------------------------------------------------------------------
  // s. Check	Functionality
  //            CMMEnergySums(JEMDAQ.JEMEnergySums) = (Crate/System)CMMEnergyDAQ.CrateCMMEnergySums


  // =============================================================================================
  // ================= RoI and Level-2 ===========================================================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Transmission
  //            JEPRoI.TotalJetEtHits = SystemCMMJetDAQ.TotalJetEtHits
  //            JEPRoI.SystemCMMEnergySums = SystemCMMEnergyDAQ.SystemCMMEnergySums
  // ---------------------------------------------------------------------------------------------

  // ---------------------------------------------------------------------------------------------
  // Check	Transmission to Level-2
  // ---------------------------------------------------------------------------------------------
  
  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            JEMJetRoI(JEMDAQ.JetElements) = JEPRoI.JEMRoIJetHits
  // ---------------------------------------------------------------------------------------------
  if (m_CompareWithSimulation ==1)
    {

      mLog << MSG::DEBUG << "==== JEM RoI: Functionality ===="<< endreq ;
      
      
      JemRoiCollection* Sim_JEMRoI = 0;
      if (jetElements) {
        Sim_JEMRoI = new JemRoiCollection;
	mLog<<MSG::DEBUG<<"Simulate RoIs from JetElements"<<endreq;
	
	InternalJemRoi* intRois = new InternalJemRoi;
	
	m_jetTool->findRoIs(jetElements, intRois);
	
	InternalJemRoi::iterator roiIter  = intRois->begin();
        InternalJemRoi::iterator roiIterE = intRois->end();
        for (; roiIter != roiIterE; ++roiIter) {
          LVL1::JEMRoI* roi = new LVL1::JEMRoI((*roiIter)->RoIWord());
	  mLog << MSG::DEBUG << "roi "<<(*roiIter)->RoIWord()<<endreq;
          Sim_JEMRoI->push_back(roi);
        }
        delete intRois;	
      }	       
       
            
      std::vector <LVL1::JEMRoI*>  vBS_JEMRoI;
      std::vector <LVL1::JEMRoI*>  vSim_JEMRoI;
      //use standard vector instead of datavector for transmissioncheck:
      //datavector erases not only the pointer  in the vector, but also the referenced object
      //-> segmentation fault!
      int foundModule;
      
      //fill vectors
      JemRoiCollection::const_iterator it_JEMRoI ;
      
           
      for( it_JEMRoI  = BS_JEMRoI->begin(); it_JEMRoI != BS_JEMRoI->end(); ++it_JEMRoI )
	{
	  vBS_JEMRoI.push_back(*it_JEMRoI);	       
	}   
      for( it_JEMRoI  = Sim_JEMRoI->begin(); it_JEMRoI != Sim_JEMRoI->end(); ++it_JEMRoI )
	{
	  vSim_JEMRoI.push_back(*it_JEMRoI);	       
	}   
      
      // Step over all cells and compare
      std::vector <LVL1::JEMRoI*>::iterator it_BS_JEMRoI;
      std::vector <LVL1::JEMRoI*>::iterator it_Sim_JEMRoI;
      
      it_BS_JEMRoI=vBS_JEMRoI.begin();
      
      mLog<<MSG::DEBUG<<"JEMRoI Calculation difference for"<<endreq;
      while (it_BS_JEMRoI != vBS_JEMRoI.end())
	{
	  foundModule = 0;
	  noMatchfound=0;
	  it_Sim_JEMRoI=vSim_JEMRoI.begin();
	  
	  while ((foundModule==0)and(it_Sim_JEMRoI != vSim_JEMRoI.end()))
	    {
	      
	      if (((*it_BS_JEMRoI)->crate()==(*it_Sim_JEMRoI)->crate())
		  and((*it_BS_JEMRoI)->jem()==(*it_Sim_JEMRoI)->jem())
		  and ((*it_BS_JEMRoI)->frame()==(*it_Sim_JEMRoI)->frame())
		  and ((*it_BS_JEMRoI)->location()==(*it_Sim_JEMRoI)->location())
		  )
		{
		  foundModule=1;
		  if ((*it_BS_JEMRoI)->hits()!=(*it_Sim_JEMRoI)->hits())
		    {
		      mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_JEMRoI)->crate()<<" Module "<<(*it_BS_JEMRoI)->jem()<<" frame "<<(*it_BS_JEMRoI)->frame() <<" location "<<(*it_BS_JEMRoI)->location() <<endreq;
		      mLog<<MSG::VERBOSE<<"BS:  RoI "<<Help.Binary((*it_BS_JEMRoI)->hits(),8)<<", word "<<Help.Binary((*it_BS_JEMRoI)->roiWord(),32)<<endreq;
		      mLog<<MSG::VERBOSE<<"Sim: RoI "<<Help.Binary((*it_Sim_JEMRoI)->hits(),8)<<", word "<<Help.Binary((*it_Sim_JEMRoI)->roiWord(),32)<<endreq;
		     
		      noMatchfound=1;
		    }
		  m_h_SimBSMon_JEP->Fill(4,((*it_BS_JEMRoI)->crate()*19+(*it_BS_JEMRoI)->jem()+1),noMatchfound);
		  overview[(*it_BS_JEMRoI)->crate()] |= (noMatchfound << 10);
		  
		  it_BS_JEMRoI=vBS_JEMRoI.erase(it_BS_JEMRoI);
		  it_Sim_JEMRoI=vSim_JEMRoI.erase(it_Sim_JEMRoI);
		}
	      else ++it_Sim_JEMRoI;
	    }
	  if (foundModule==0)++it_BS_JEMRoI;
	}
      
      if (vBS_JEMRoI.size()!=0)
	{
	  mLog<<MSG::DEBUG<<"JEMRoI: additional BS data for "<<endreq;
	  
	  //fill error counter
	  for( it_BS_JEMRoI  = vBS_JEMRoI.begin(); it_BS_JEMRoI != vBS_JEMRoI. end(); ++it_BS_JEMRoI )
	    {
	      mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_JEMRoI)->crate()<<" Module "<<(*it_BS_JEMRoI)->jem()<<" frame "<<(*it_BS_JEMRoI)->frame() <<" location "<<(*it_BS_JEMRoI)->location() <<endreq;
	      mLog<<MSG::VERBOSE<<"BS: RoI "<<Help.Binary((*it_BS_JEMRoI)->hits(),24)<<endreq;
	      
	      m_h_SimBSMon_JEP->Fill(4,((*it_BS_JEMRoI)->crate()*19+(*it_BS_JEMRoI)->jem()+1),1);
              overview[(*it_BS_JEMRoI)->crate()] |= (1 << 10);
	    }
	}
      
      if (vSim_JEMRoI.size()!=0)
	{
	  mLog<<MSG::DEBUG<<"JEMRoI: additional Sim data for"<<endreq;
	  mLog<<MSG::DEBUG<<"fill error counter but remember simulation is not zero suppressing"<<endreq;
	 
	  //fill error counter if simulation data is !=0
	  for( it_Sim_JEMRoI  = vSim_JEMRoI.begin(); it_Sim_JEMRoI != vSim_JEMRoI. end(); ++it_Sim_JEMRoI )
	    {
	      mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_JEMRoI)->crate()<<" Module "<<(*it_Sim_JEMRoI)->jem()<<" frame "<<(*it_Sim_JEMRoI)->frame() <<" location "<<(*it_Sim_JEMRoI)->location() <<endreq;
	      mLog<<MSG::VERBOSE<<"Sim: RoI "<<Help.Binary((*it_Sim_JEMRoI)->hits(),24)<<endreq;
	      
	      if( (*it_Sim_JEMRoI)->hits()!=0) 
		{
		  m_h_SimBSMon_JEP->Fill(4,((*it_Sim_JEMRoI)->crate()*19+(*it_Sim_JEMRoI)->jem()+1),1);
                  overview[(*it_Sim_JEMRoI)->crate()] |= (1 << 10);
		}
	    }
	}
      delete Sim_JEMRoI;
    }
  
  // ================= CMM RoI ========================================================
  
  // retrieve RoI information from Storegate
      const LVL1::CMMRoI* BS_CR = 0;
      sc = m_storeGate->retrieve (BS_CR, m_BS_CMMRoILocation);
      if (sc==StatusCode::FAILURE) {
	  mLog <<MSG::INFO<<"No BS CMM RoI found in TES at "<< m_BS_CMMRoILocation<<endreq;
	  return StatusCode::SUCCESS;    
      }
  
    if (m_CompareWithSimulation ==1) {

      mLog << MSG::DEBUG << "==== CMM RoI: Functionality ===="<< endreq ;
      
      CMMJetHitsCollection* Sim_JetEt = 0;
      if(CMMJetHits) {
      Sim_JetEt = new CMMJetHitsCollection;
      mLog<<MSG::DEBUG<<"Simulate JetEt HitMap from CMMJetHits"<<endreq;
      m_jepHitsTool->formCMMJetHitsEtMap(CMMJetHits, Sim_JetEt);
      }
      
    bool empty = false;
      
    
     for( it_CMMJetHits  = Sim_JetEt ->begin(); it_CMMJetHits < Sim_JetEt -> end(); ++it_CMMJetHits ) {
       empty = true;
       mLog<<MSG::DEBUG<<"Sim JetHits found"<<endreq;
     }
    	
    it_CMMJetHits = Sim_JetEt -> begin();
    
    if (empty) {
      
      mLog<<MSG::DEBUG<<(*BS_CR).jetEtHits()<<endreq;
      mLog<<MSG::DEBUG<<(*it_CMMJetHits)->Hits()<<endreq;
	           
      noMatchfound=0;  
      if ( (*BS_CR).jetEtHits()!= (int)(*it_CMMJetHits)->Hits() )  //signed int compare to unsigned int
	{
	  mLog<<MSG::DEBUG<<"CMMRoI: No Match found between BS and Sim for JetEtHits"<<endreq;
	  mLog<<MSG::VERBOSE<<"BS: JetEtHits: "<<(*BS_CR).jetEtHits()<<endreq;
	      
	  mLog<<MSG::VERBOSE<<"Sim: JetEtHits: "<<(*it_CMMJetHits)->Hits()<<endreq;
	  
	  noMatchfound=1; 
	}
    }
      
    if (!empty and (*BS_CR).jetEtHits()!=0) {
      noMatchfound=1;
      mLog<<MSG::DEBUG<<"no sim but... "<<(*BS_CR).jetEtHits()<<endreq;
    }
     
      
      CMMEtSumsCollection* Sim_Et = 0;
      if(CMMEtSums) {
      Sim_Et = new CMMEtSumsCollection;
      mLog<<MSG::DEBUG<<"Simulate SumEt/MissingEt HitMap from CMMEtSums"<<endreq;
      m_etSumsTool->formCMMEtSumsEtMaps(CMMEtSums, Sim_Et);      
      }
      
      CMMEtSumsCollection::const_iterator it_CMMEtSum;	
	  	     
      
      empty = false;
      
      for( it_CMMEtSum  = Sim_Et ->begin(); it_CMMEtSum < Sim_Et -> end(); ++it_CMMEtSum )
      {  
        empty = true;
        mLog<<MSG::DEBUG<<"Sim SumEtHits found"<<endreq;
      }
      
      it_CMMEtSum = Sim_Et->begin();
      
      if (empty) {
        if ((*it_CMMEtSum)-> dataID() == 19 and (*it_CMMEtSum)->crate()==1) {
          if ( ((BS_CR)->missingEtHits()) != (int)(*it_CMMEtSum)->Et()) {
	    mLog<<MSG::DEBUG<<"CMMRoI: No Match found between BS and Sim for MissingEtHits"<<endreq;
	    mLog<<MSG::VERBOSE<<"BS: MissingEtHits: "<<(BS_CR)->missingEtHits()<<endreq;	 
    	    mLog<<MSG::VERBOSE<<"Sim: MissingEtHits: "<<(*it_CMMEtSum)->Et()<<endreq;	

	    noMatchfound=1; 
	  }
        }
	
        if ((*it_CMMEtSum)-> dataID() == 20 and (*it_CMMEtSum)->crate()==1) {
          if ( (BS_CR)->sumEtHits()!=(int)(*it_CMMEtSum)->Et()) {
	    mLog<<MSG::DEBUG<<"CMMRoI: No Match found between BS and Sim for SumEtHits"<<endreq;
	    mLog<<MSG::VERBOSE<<"BS: SumEtHits: "<<(*BS_CR).sumEtHits()<<endreq;	 
            mLog<<MSG::VERBOSE<<"Sim: SumEtHits: "<<(*it_CMMEtSum)->Et()<<endreq;	

	    noMatchfound=1; 
	  }
        }
       
        if ((*it_CMMEtSum)-> dataID() == 18 and (*it_CMMEtSum)->crate()==1) {      
          if ( ((BS_CR)->ex()!=(int)(*it_CMMEtSum)->Ex())or((BS_CR)->ey()!=(int)(*it_CMMEtSum)->Ey())or((BS_CR)->et()!=(int)(*it_CMMEtSum)->Et())) {
	    mLog<<MSG::DEBUG<<"CMMRoI: No Match found between BS and Sim for Ex, Ey or Et"<<endreq;
	    mLog<<MSG::VERBOSE<<"BS: Et: "<<(*BS_CR).et()<<endreq;
	    mLog<<MSG::VERBOSE<<"BS: Ex: "<<(*BS_CR).ex()<<endreq;
	    mLog<<MSG::VERBOSE<<"BS: Ey: "<<(*BS_CR).ey()<<endreq;	 
    
	    mLog<<MSG::VERBOSE<<"Sim: Et: "<<(*it_CMMEtSum)->Et()<<endreq;
	    mLog<<MSG::VERBOSE<<"Sim: Ex: "<<(*it_CMMEtSum)->Ex()<<endreq;
	    mLog<<MSG::VERBOSE<<"Sim: Ey: "<<(*it_CMMEtSum)->Ey()<<endreq;
	  
	    noMatchfound=1; 
	  }
        }
      }
     
      if ((!empty) and (((BS_CR)->missingEtHits()!=0) or ((BS_CR)->sumEtHits()!=0))) {
        noMatchfound=1;
        mLog<<MSG::VERBOSE<<"no sim but additional data for XE"<<(BS_CR)->missingEtHits()<<endreq;
        mLog<<MSG::VERBOSE<<"no sim but additional data for TE"<<(*BS_CR).sumEtHits()<<endreq;
      }
     
     
      m_h_SimBSMon_JEP->Fill(4,(19+16+1),noMatchfound);
      overview[1] |= (noMatchfound << 11);
    
      delete Sim_JetEt;
      delete Sim_Et;
    }

  // Write overview vector to StoreGate
  std::vector<int>* save = new std::vector<int>(overview);
  sc = m_storeGate->record(save, "L1CaloJEMMismatchVector");
  if (sc != StatusCode::SUCCESS)
    {
      mLog << MSG::ERROR << "Error recording JEM mismatch vector in TES "
           << endreq;
      return sc;
    }

  return StatusCode( StatusCode::SUCCESS );
}

/*---------------------------------------------------------*/
StatusCode JEPTransPerfMon::procHistograms( bool isEndOfEventsBlock, 
				  bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  
  mLog << MSG::DEBUG << "in procHistograms" << endreq ;
  
  
  if( isEndOfEventsBlock || isEndOfLumiBlock ) 
    {  
    }
 
  if(m_Offline==1)
    {      
      if( isEndOfRun ) 
	{
	  std::stringstream buffer;
	  buffer.str("");
	  buffer<<m_NoEvents;
	  std::string title;
	  std::stringstream runNumStr;
	  runNumStr.str("");
	  runNumStr<<m_runNum;
	  	  
	  if (m_CompareWithSimulation ==1)
	    {
	      title = m_h_SimBSMon_JEP-> GetTitle();
	      title="Run "+runNumStr.str()+": "+title+ " | #events: " + buffer.str();
	      m_h_SimBSMon_JEP->SetTitle(title.c_str());
	    }

	  title = m_h_TransCheck_JEP-> GetTitle();
	  title="Run "+runNumStr.str()+": "+title + " | #events: " + buffer.str();
	  m_h_TransCheck_JEP->SetTitle(title.c_str());
	  
	  title = m_h_TransCheck_emJetElements-> GetTitle();
	  title="Run "+runNumStr.str()+": "+title + " | #events: " + buffer.str();
	  m_h_TransCheck_emJetElements->SetTitle(title.c_str());
	  
	  title = m_h_TransCheck_hadJetElements-> GetTitle();
	  title="Run "+runNumStr.str()+": "+title + " | #events: " + buffer.str();
	  m_h_TransCheck_hadJetElements->SetTitle(title.c_str());
	}
  }
  
  return StatusCode( StatusCode::SUCCESS );
}


//---------------------------Functionality: Trigger Towers -> JetElements-----------------

/*---------------------------------------------------------*/
void JEPTransPerfMon::TimeSliceMatch(int k, int TT_TS, const JECollection* TT_jetElements, int JE_TS, const JECollection* jetElements, std::vector<int>& overview, MsgStream::MsgStream* mLog)
/*---------------------------------------------------------*/
{

  *mLog << MSG::VERBOSE << "TT="<<TT_TS<<", JE="<<JE_TS << endreq ;

  
  int JEFound;
  int em_NoJEMatchFound;
  int had_NoJEMatchFound;

  JECollection::const_iterator it_je ;
  JECollection::const_iterator it_TT_je ;
  
    
  for( it_je = jetElements ->begin(); it_je < jetElements->end(); ++it_je )
    {	  
      it_TT_je = TT_jetElements ->begin();
      JEFound = 0;
      em_NoJEMatchFound = 1;
      had_NoJEMatchFound = 1;
      

      LVL1::CoordToHardware ToHW;
      LVL1::Coordinate coord((*it_je)->phi(),(*it_je)->eta());
      
      int crate = ToHW.jepCrate(coord);
      int module=ToHW.jepModule(coord);
      
           
      while ((JEFound==0)and(it_TT_je < TT_jetElements->end()))
	{   
	      
	  if ((*it_TT_je)->key()==(*it_je)->key())
	    {
	      JEFound=1;
	      *mLog << MSG::VERBOSE << " HW: Module: " << module<<" crate: "<<crate <<" eta "<<(*it_je)->eta()<<" phi "<< (*it_je)->phi()<<endreq ;
              *mLog << MSG::VERBOSE << " Simulation: Module: " << module<<" crate: "<<crate <<" eta "<<(*it_TT_je)->eta()<<" phi "<< (*it_TT_je)->phi()<<endreq ;
	      
	       // Simulation may have energy vector values greater than saturation
       	      const int sat = 511; // layer saturation
	      if ((*it_TT_je)->emEnergyVec()[TT_TS]==(*it_je)->emEnergyVec()[JE_TS] ||((*it_TT_je)->emEnergyVec()[TT_TS]>=sat && (*it_je)->emEnergyVec()[JE_TS]>=sat)) em_NoJEMatchFound=0;        
              if ((*it_TT_je)->hadEnergyVec()[TT_TS]==(*it_je)->hadEnergyVec()[JE_TS] ||((*it_TT_je)->hadEnergyVec()[TT_TS]>=sat && (*it_je)->hadEnergyVec()[JE_TS]>=sat)) had_NoJEMatchFound=0;  
	      
	     
	    }
	  if ((JEFound==0)and(it_TT_je < TT_jetElements->end())) it_TT_je ++;
	}

      if (em_NoJEMatchFound==1 and JEFound==1)
	{
	  *mLog << MSG::VERBOSE << "no em Match Found eta: "<<(*it_TT_je)->eta()<<" phi: " <<(*it_TT_je)->phi()<<" energy HW: " << (*it_je)->emEnergyVec()[JE_TS]<<" SW: " << (*it_TT_je)->emEnergyVec()[TT_TS]<< endreq ;
	  *mLog << MSG::VERBOSE << "PPM data from Crate: "<<crate<<" Module: " <<module<< endreq ;
	}

      if (had_NoJEMatchFound==1 and JEFound==1)
	{
	  *mLog << MSG::VERBOSE << "no had Match Found eta: "<<(*it_TT_je)->eta()<<" phi: " <<(*it_TT_je)->phi()<<" energy HW: " << (*it_je)->hadEnergyVec()[JE_TS]<<" SW: " << (*it_TT_je)->hadEnergyVec()[TT_TS]<< endreq ;
	  *mLog << MSG::VERBOSE << "PPM data from Crate: "<<crate<<" Module: " <<module<< endreq ;
	}

      if (JEFound==0)
      {
   	// Not a mismatch if energy zero
   	if ((*it_je)->emEnergyVec()[JE_TS] == 0) em_NoJEMatchFound=0;
   	if ((*it_je)->hadEnergyVec()[JE_TS] == 0) had_NoJEMatchFound=0;
   	if (em_NoJEMatchFound || had_NoJEMatchFound) {
     	  *mLog <<MSG::INFO   <<"JEFound=0"<<endreq;
     	  *mLog << MSG::INFO    << " Module: " << module<<" crate: "<<crate <<" eta "<<(*it_je)->eta()<<" phi "<< (*it_je)->phi()<<endreq ;
   	  }
      }
	
      m_h_TransCheck_emJetElements->Fill(k,(crate*18+module+1),em_NoJEMatchFound);
      m_h_TransCheck_hadJetElements->Fill(k,(crate*18+module+1),had_NoJEMatchFound);
      overview[crate] |= em_NoJEMatchFound;
      overview[crate] |= (had_NoJEMatchFound << 1);

 
 *mLog <<MSG::VERBOSE<<em_NoJEMatchFound<<" and "<<had_NoJEMatchFound<<endreq;

    }
}

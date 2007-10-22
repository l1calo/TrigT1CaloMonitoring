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
// Check	Funkionalität
//              TotalCMMJetHits(CrateCMMJetDAQ.CrateCMMJetHits+ SystemCMMJetDAQ.CrateCMMJetHits) = SystemCMMJetDAQ.SystemCMMJetHits

// ================= Crate CMM Energy Sums -> System CMM Energy Sums ===========================
// Check	Transmission
//              CrateCMMEnergyDAQ.CrateCMMEnergySums = SystemCMMEnergyDAQ.RemoteCMMEnergySums
// Check	Funkionalität
//              TotalCMMEnergySums(CrateCMMEnergyDAQ.CrateCMMEnergySums+ SystemCMMEnergyDAQ.CrateCMMEnergySums) = SystemCMMEnergyDAQ.SystemCMMEnergySums

// ================= RoI and Level-2 ===========================================================
// Check	Transmission
//              JEPRoI.TotalJetEtHits = SystemCMMJetDAQ.TotalJetEtHits
//              JEPRoI.SystemCMMEnergySums = SystemCMMEnergyDAQ.SystemCMMEnergySums
// Check	Funktionality
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

#include "TrigT1CaloMonitoring/JEPTransPerfMon.h"
#include "TrigT1CaloMonitoring/JEMMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

//#include "TrigT1Calo/EnergyTrigger.h"
#include "TrigT1Calo/LVL1TriggerMenuDefs.h"
//#include "TrigT1Calo/LVL1TriggerMenu.h"
//#include "TrigT1Calo/InternalJetROI.h"
#include "TrigT1Calo/CMMRoI.h"
#include "TrigT1Calo/JEMRoI.h"
#include "TrigT1Calo/QuadLinear.h"
#include "TrigT1Calo/DataError.h"
#include "TrigT1Calo/JetElementMaker.h"
#include "TrigT1Calo/TriggerTowerCollection.h"

#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"
#include "TrigT1Calo/CoordToHardware.h"

#include "AthenaMonitoring/AthenaMonManager.h"




#include "Identifier/HWIdentifier.h"




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
  : ManagedMonitorToolBase( type, name, parent )//,
    //m_JetElementTool("LVL1::L1JetElementTools/L1JetElementTools")
/*---------------------------------------------------------*/
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file
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
  declareProperty( "EventNoInHistoTitle", m_EventNoInHisto = 1) ;


  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/TransAndPerf") ;
  /*
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
  */
}


/*---------------------------------------------------------*/
JEPTransPerfMon::~JEPTransPerfMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode JEPTransPerfMon::bookHistograms( bool isNewEventsBlock, 
				   bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in JEPTransPerfMon::bookHistograms" << endreq;
  
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
  














  ManagedMonitorToolBase::LevelOfDetail_t LevelOfDetail=shift;
  if (m_DataType=="Sim") LevelOfDetail = expert;

  MonGroup CMM_transmission ( this, (m_PathInRootFile ).c_str(), shift, eventsBlock );
  HistoBooker* transmission_Booker = new HistoBooker(&CMM_transmission, &mLog, "");
  
  MonGroup SimBSComparison_JEM ( this, (m_PathInRootFile).c_str(), shift, eventsBlock );
  HistoBooker* JEM_Booker = new HistoBooker(&SimBSComparison_JEM, &mLog, m_DataType);

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  if ( isNewEventsBlock|| isNewLumiBlock) { }
  
  if( isNewRun ) 
    {	
      m_NoEvents=0;

      std::string name;
      std::stringstream buffer,buffer2;

      m_h_usedModules=transmission_Booker->book2F("JEM_usedModules", "JEM used modules", 4,0.5,4.5,35,0.5,35.5, "", "");

      m_h_TransCheck_emJetElements=transmission_Booker->book2F("emJE_TransPerfCheck", "em JE Transmission and Performance Check for all Timeslices (TT|JE)", (5*m_NoLUTSlices),0.5,(5*m_NoLUTSlices+.5), 35,0.5,35.5, "", "");
      m_h_TransCheck_hadJetElements=transmission_Booker->book2F("hadJE_TransPerfCheck", "had JE Transmission and Performance Check for all Timeslices (TT|JE)",(5*m_NoLUTSlices),0.5,(5*m_NoLUTSlices+.5), 35,0.5,35.5, "", "");

      for (int i = 0; i < 16; i++)
	{
	  buffer.str("");
	  buffer<<i;
	  
	  name = "JEM " + buffer.str();
	  m_h_TransCheck_emJetElements->GetYaxis()->SetBinLabel((i+1), name.c_str());
	  m_h_TransCheck_emJetElements->GetYaxis()->SetBinLabel((i+1+18), name.c_str());
	  m_h_TransCheck_hadJetElements->GetYaxis()->SetBinLabel((i+1), name.c_str());
	  m_h_TransCheck_hadJetElements->GetYaxis()->SetBinLabel((i+1+18), name.c_str());
	  m_h_usedModules->GetYaxis()->SetBinLabel((i+1), name.c_str());
	  m_h_usedModules->GetYaxis()->SetBinLabel((i+1+18), name.c_str());
	}
      m_h_TransCheck_emJetElements->GetYaxis()->SetBinLabel(17, "Crate 0");
      m_h_TransCheck_emJetElements->GetYaxis()->SetBinLabel(35, "Crate 1");
      m_h_TransCheck_hadJetElements->GetYaxis()->SetBinLabel(17, "Crate 0");
      m_h_TransCheck_hadJetElements->GetYaxis()->SetBinLabel(35, "Crate 1");
      m_h_usedModules->GetYaxis()->SetBinLabel(17, "Crate 0");
      m_h_usedModules->GetYaxis()->SetBinLabel(35, "Crate 1");

      int k=1;
      for (int i=0; i<5;i++)
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

      m_h_usedModules->GetXaxis()->SetBinLabel(1, "HW em");
      m_h_usedModules->GetXaxis()->SetBinLabel(2, "HW had");
      m_h_usedModules->GetXaxis()->SetBinLabel(3, "SW em");
      m_h_usedModules->GetXaxis()->SetBinLabel(4, "SW had");

  m_h_SimBSMon_JEP=JEM_Booker->book2F("JEP_Calc_Error", "JEP Hardware Output compared to Simulation: Differences", 3,0.5,3.5,37,0.5,37.5, "", "");
  //m_h_SimBSMon_JEP-> SetOption ("text");

  m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(1, "Energy");
  m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(2, "Hits");
  m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(3, "RoI");

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

 
      //---------------------------------- Backplane transmission checks -----------------------------
      m_h_TransCheck_JEP=transmission_Booker->book2F("JEP_TransCheck", "JEP Backplane Transmission Check JEM -> CMM per Module and Crate", 2,0.5,2.5,37,0.5,37.5, "", "");
      //m_h_TransCheck_JEP-> SetOption ("text");
      m_h_TransCheck_JEP->GetXaxis()->SetBinLabel(1, "Hits");
      m_h_TransCheck_JEP->GetXaxis()->SetBinLabel(2, "Energy");
      
      for (int i = 0; i < 16; i++)
	{
	  buffer.str("");
	  buffer<<i;
	  
	  name = "JEM " + buffer.str();
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel((i+1), name.c_str());
	  	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel((i+1+20), name.c_str());
	}
      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(17, "C E CMM");
      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(18, "C J CMM");
      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(19, "Crate 0: ");
      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(37, "Crate 1: ");


   }
  return StatusCode( StatusCode::SUCCESS );
}


/*---------------------------------------------------------*/
StatusCode JEPTransPerfMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  Helper* Help = new Helper();
 
  m_NoEvents++;

  // =============================================================================================
  // ================= TriggerTowers -> JetElements ==============================================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Transmission & Functionality
  //            JetElements(PPrDAQ.TriggerTowers) = JEMDAQ.JetElements
  // ---------------------------------------------------------------------------------------------

  mLog << MSG::DEBUG << "==== TriggerTowers -> JetElements: Transmission & Functionality ===="<< endreq ;

  // retrieve JetElements
  const JECollection* jetElements;
  StatusCode sc = m_storeGate->retrieve(jetElements, m_BS_JetElementLocation);
  
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JetElements found in TES at " << m_BS_JetElementLocation << endreq ;
      return StatusCode::SUCCESS;
    }




  //--------------------- mask out for M5 ------------------------------------------------------------
  // retrieve TriggerTower
  /*
  const TriggerTowerCollection* TriggerTowers;
  sc = m_storeGate->retrieve(TriggerTowers, m_BS_TriggerTowerLocation);
  
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No TriggerTowers found in TES at " << m_BS_TriggerTowerLocation << endreq ;
      return StatusCode::SUCCESS;
    }

  // make JetElements from TriggerTowers
  sc = m_JetElementTool.retrieve();
  if (sc==StatusCode::FAILURE) mLog << MSG::ERROR << "Problem retrieving JetElementTool" << endreq;
  
  DataVector<JetElement>* TT_jetElements = new DataVector<JetElement>;
  m_JetElementTool->makeJetElements(TriggerTowers, TT_jetElements);
  */
  //--------------------------------------------------------------------------------------------------


  //----------------------- insert for M5 ------------------------------------------------------------
  // retrieve JetElements
  const JECollection* TT_jetElements;
  sc = m_storeGate->retrieve(TT_jetElements, m_BS_TriggerTowerLocation);
  
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JetElements found in TES at " << m_BS_TriggerTowerLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  //--------------------------------------------------------------------------------------------------

 
  JECollection::const_iterator it_je, it_TT_je ;
  
  for( it_je = jetElements ->begin(); it_je < jetElements->end(); ++it_je )
    {
      LVL1::CoordToHardware ToHW;
      LVL1::Coordinate coord((*it_je)->phi(),(*it_je)->eta());
      
      int crate = ToHW.jepCrate(coord);
      int module=ToHW.jepModule(coord);

      // Histograms for the Modules actually used during this run
      if ((*it_je)->emEnergyVec()[0]!=0 or (*it_je)->emEnergyVec()[1]!=0 or (*it_je)->emEnergyVec()[2]!=0 or (*it_je)->emEnergyVec()[3]!=0 or (*it_je)->emEnergyVec()[4]!=0) m_h_usedModules->Fill(1,(crate*18+module+1),1);

      if ((*it_je)->hadEnergyVec()[0]!=0 or (*it_je)->hadEnergyVec()[1]!=0 or (*it_je)->hadEnergyVec()[2]!=0 or (*it_je)->hadEnergyVec()[3]!=0 or (*it_je)->hadEnergyVec()[4]!=0) m_h_usedModules->Fill(2,(crate*18+module+1),1);
    }

  for( it_TT_je = TT_jetElements ->begin(); it_TT_je < TT_jetElements->end(); ++it_TT_je )
    {
      LVL1::CoordToHardware ToHW;
      LVL1::Coordinate coord((*it_TT_je)->phi(),(*it_TT_je)->eta());
      
      int crate = ToHW.jepCrate(coord);
      int module=ToHW.jepModule(coord);

      if ((*it_TT_je)->emEnergyVec()[0]!=0) m_h_usedModules->Fill(3,(crate*18+module+1),1);
      if ((*it_TT_je)->hadEnergyVec()[0]!=0) m_h_usedModules->Fill(4,(crate*18+module+1),1);
      
      //mLog<<MSG::DEBUG<<"Crate: "<<(*it_TT_je)->()[i]<<endreq;
	  //mLog<<MSG::DEBUG<<"TT JE em Energies: "<<(*it_TT_je)->emEnergyVec()[i]<<endreq;
	  //mLog<<MSG::DEBUG<<"TT JE had Energies: "<<(*it_TT_je)->hadEnergyVec()[i]<<endreq;
    }



  //JEPTransPerfMon::TimeSliceMatch(int k, int TT_TS, const JECollection* TT_jetElements, int JE_TS, const JECollection* jetElements, MsgStream::MsgStream* mLog)
  int k=1;
  for (int i=0; i<5;i++)
    {
      for (int j=0; j<m_NoLUTSlices; j++)
	{
	  TimeSliceMatch(k, j, TT_jetElements, i, jetElements, &mLog);
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

  mLog << MSG::DEBUG << "==== JetElements -> JEM Jet Hits: Functionality ===="<< endreq ;


  const JEMHitsCollection* BS_JEMHits;
  sc = m_storeGate->retrieve(BS_JEMHits, m_BS_JEMHitsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMHits found in TES at "<< m_BS_JEMHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  const JEMHitsCollection* Sim_JEMHits;
  sc = m_storeGate->retrieve(Sim_JEMHits, m_Sim_JEMHitsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMHits found in TES at "<< m_Sim_JEMHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  std::vector <LVL1::JEMHits>  vBS_JEMHits;
  std::vector <LVL1::JEMHits>  vSim_JEMHits;
  //use standard vector instead of datavector for transmissioncheck:
  //datavector erases not only the pointer  in the vector, but also the referenced object
  //-> segmentation fault!
  
  //fill vectors
  JEMHitsCollection::const_iterator it_JEMHits ;
  int  foundModule;
  int noMatchfound;

  for( it_JEMHits  = BS_JEMHits->begin(); it_JEMHits < BS_JEMHits->end(); ++it_JEMHits )
    {
      vBS_JEMHits.push_back(**it_JEMHits);	       
    }   
  for( it_JEMHits  = Sim_JEMHits->begin(); it_JEMHits < Sim_JEMHits->end(); ++it_JEMHits )
    {
      vSim_JEMHits.push_back(**it_JEMHits);	       
    }   
  
  // Step over all cells and compare
  std::vector <LVL1::JEMHits>::iterator it_BS_JEMHits;
  std::vector <LVL1::JEMHits>::iterator it_Sim_JEMHits;
  
  
  it_BS_JEMHits=vBS_JEMHits.begin();
  
  mLog<<MSG::DEBUG<<"JEMHits Calculation difference for"<<endreq;
  while (it_BS_JEMHits<vBS_JEMHits.end())
    {
      foundModule = 0;
      noMatchfound=0;
      it_Sim_JEMHits=vSim_JEMHits.begin();
      
      while ((foundModule==0)and(it_Sim_JEMHits<vSim_JEMHits.end()))
	{
  
	  if (((*it_BS_JEMHits).crate()==(*it_Sim_JEMHits).crate())
	      and((*it_BS_JEMHits).module()==(*it_Sim_JEMHits).module()))
	    {
	      foundModule=1;
	      
		if ((*it_BS_JEMHits).JetHits()!=(*it_Sim_JEMHits).JetHits())
		  {
		    mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_JEMHits).crate()<<" Module "<<(*it_BS_JEMHits).module()<<endreq;
		    mLog<<MSG::VERBOSE<<"BS:  Hits "<<Help->Binary((*it_BS_JEMHits).JetHits(),24)<<endreq;
		    mLog<<MSG::VERBOSE<<"Sim: Hits "<<Help->Binary((*it_Sim_JEMHits).JetHits(),24)<<endreq;
		    noMatchfound=1;
		  }
		m_h_SimBSMon_JEP->Fill(2,((*it_BS_JEMHits).crate()*19+(*it_BS_JEMHits).module()+1),noMatchfound);
		
		vBS_JEMHits.erase(it_BS_JEMHits);
		vSim_JEMHits.erase(it_Sim_JEMHits);
	    }
	  else it_Sim_JEMHits=it_Sim_JEMHits+1;
	}
      if (foundModule==0)it_BS_JEMHits=it_BS_JEMHits+1;
    }
  
  if (vBS_JEMHits.size()!=0)
    {
      mLog<<MSG::DEBUG<<"JEMHits: additional BS data for"<<endreq;
      
      //fill errorcounter
      for( it_BS_JEMHits  = vBS_JEMHits.begin(); it_BS_JEMHits <  vBS_JEMHits. end(); ++it_BS_JEMHits )
	{
	  mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_JEMHits).crate()<<" Module "<<(*it_BS_JEMHits).module()<<endreq;
	  mLog<<MSG::VERBOSE<<"BS: Hits "<<Help->Binary((*it_BS_JEMHits).JetHits(),24)<<endreq;
	  
	  if ((*it_BS_JEMHits).JetHits()>0) m_h_SimBSMon_JEP->Fill(2,((*it_BS_JEMHits).crate()*19+(*it_BS_JEMHits).module()+1),1);	    
	}
    }
  
  if (vSim_JEMHits.size()!=0)
    {
      mLog<<MSG::DEBUG<<"JEMHits: additional Sim data for"<<endreq;
      
      //fill errorcounter
      for( it_Sim_JEMHits  = vSim_JEMHits.begin(); it_Sim_JEMHits <  vSim_JEMHits. end(); ++it_Sim_JEMHits )
	{
	  mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_JEMHits).crate()<<" Module "<<(*it_Sim_JEMHits).module()<<endreq;
	  mLog<<MSG::VERBOSE<<"Sim: Hits "<<Help->Binary((*it_Sim_JEMHits).JetHits(),24)<<endreq;
	  
	  //m_h_SimBSMon_JEP->Fill(2,((*it_Sim_JEMHits).crate()*19+(*it_Sim_JEMHits).module()+1),1);	    
	}
    }
  

  // =============================================================================================
  // ================= JetElements -> JEM Energy Sums ============================================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            JEMEnergySums(JEMDAQ.JetElements) = JEMDAQ.JEMEnergySums
  // ---------------------------------------------------------------------------------------------

  mLog << MSG::DEBUG << "==== JetElements -> JEM Energy Sums: Functionality ===="<< endreq ;

  // retrieve JEMEtSums information from simulation and bytestream... 
  const JEMEtSumsCollection* BS_JEMEtSums;
  sc = m_storeGate->retrieve(BS_JEMEtSums, m_BS_JEMEtSumsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMEtSums found in TES at "<< m_BS_JEMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  const JEMEtSumsCollection* Sim_JEMEtSums;
  sc = m_storeGate->retrieve(Sim_JEMEtSums, m_Sim_JEMEtSumsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMEtSums found in TES at "<< m_Sim_JEMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  // compare simulation-vector with bytestream-vector, erase corresponding entries
  // -> the remaining vectors contain faulty transmissions! 
  // use standard vector instead of datavector for transmissioncheck:
  // datavector erases not only the pointer  in the vector, but also the referenced object
  // -> segmentation fault!
  std::vector <LVL1::JEMEtSums>  vBS_JEMEtSums;
  std::vector <LVL1::JEMEtSums>  vSim_JEMEtSums;
  
  //fill vectors with simulation and bytestream data
  JEMEtSumsCollection::const_iterator it_JEMEtSums ;
  for( it_JEMEtSums  = BS_JEMEtSums->begin(); it_JEMEtSums < BS_JEMEtSums->end(); ++it_JEMEtSums )
    {
      vBS_JEMEtSums.push_back(**it_JEMEtSums);	       
    }   
  for( it_JEMEtSums  = Sim_JEMEtSums->begin(); it_JEMEtSums < Sim_JEMEtSums->end(); ++it_JEMEtSums )
    {
      vSim_JEMEtSums.push_back(**it_JEMEtSums);	       
    }   
  
  // Step over all cells and compare
  //bool foundModule;
  //int noMatchfound = 0; // used to fill: "0" if a match is found; "1" if no match is found
  std::vector <LVL1::JEMEtSums>::iterator it_BS_JEMEtSums;
  std::vector <LVL1::JEMEtSums>::iterator it_Sim_JEMEtSums;
  
  it_BS_JEMEtSums=vBS_JEMEtSums.begin();
  
  mLog<<MSG::DEBUG<<"JEMEtSums Calculation difference for"<<endreq;  
  while (it_BS_JEMEtSums<vBS_JEMEtSums.end())
    {
      foundModule = 0;
      it_Sim_JEMEtSums=vSim_JEMEtSums.begin();
      
      while ((foundModule==0)and(it_Sim_JEMEtSums<vSim_JEMEtSums.end()))
	{	  	  
	  if (((*it_BS_JEMEtSums).crate()==(*it_Sim_JEMEtSums).crate())
	      and((*it_BS_JEMEtSums).module()==(*it_Sim_JEMEtSums).module()))
	    {
	      foundModule=1;
	      noMatchfound = 0;
	      
	      if ((*it_BS_JEMEtSums).crate()==0)
		{
		  if (((*it_BS_JEMEtSums).Ex()!=(*it_Sim_JEMEtSums).Ex())
		      or ((*it_BS_JEMEtSums).Ey()!=(*it_Sim_JEMEtSums).Ey())
		      or ((*it_BS_JEMEtSums).Et()!=(*it_Sim_JEMEtSums).Et()))
		    {
		      ;
		      mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_JEMEtSums).crate()<<" Module "<<(*it_BS_JEMEtSums).module()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_JEMEtSums).Ex()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_JEMEtSums).Ey()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_JEMEtSums).Et()<<endreq;
		      
		      mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_JEMEtSums).Ex()<<endreq;
		      mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_JEMEtSums).Ey()<<endreq;
		      mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_JEMEtSums).Et()<<endreq;
		      
		      noMatchfound = 1;
		    }
		}
	      else
		{
		  if (((*it_BS_JEMEtSums).Ex()!=(*it_Sim_JEMEtSums).Ey())
		      or ((*it_BS_JEMEtSums).Ey()!=(*it_Sim_JEMEtSums).Ex())
		      or ((*it_BS_JEMEtSums).Et()!=(*it_Sim_JEMEtSums).Et()))
		    {
		      ;
		      mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_JEMEtSums).crate()<<" Module "<<(*it_BS_JEMEtSums).module()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_JEMEtSums).Ex()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_JEMEtSums).Ey()<<endreq;
		      mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_JEMEtSums).Et()<<endreq;
		      
		      mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_JEMEtSums).Ex()<<endreq;
		      mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_JEMEtSums).Ey()<<endreq;
		      mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_JEMEtSums).Et()<<endreq;
		      
		      noMatchfound = 1;
		    }
		}
	      m_h_SimBSMon_JEP->Fill(1,(*it_BS_JEMEtSums).crate()*19 + ((*it_BS_JEMEtSums).module()+1),noMatchfound);
	      
	      vBS_JEMEtSums.erase(it_BS_JEMEtSums);
	      vSim_JEMEtSums.erase(it_Sim_JEMEtSums);
	    }
	  else it_Sim_JEMEtSums=it_Sim_JEMEtSums+1;
	}
      if (foundModule==0)it_BS_JEMEtSums=it_BS_JEMEtSums+1;
    }
  
  if (vBS_JEMEtSums.size()!=0)
    {
      mLog<<MSG::DEBUG<<"JEMEtSums: additional BS data for"<<endreq;
      
      //fill errorcounter
      for( it_BS_JEMEtSums  = vBS_JEMEtSums.begin(); it_BS_JEMEtSums <  vBS_JEMEtSums. end(); ++it_BS_JEMEtSums )
	{
	  mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_JEMEtSums).crate()<<" Module "<<(*it_BS_JEMEtSums).module()<<endreq;
	  mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_JEMEtSums).Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_JEMEtSums).Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_JEMEtSums).Et()<<endreq;
	  
	  m_h_SimBSMon_JEP->Fill(1,(*it_BS_JEMEtSums).crate()*19 + ((*it_BS_JEMEtSums).module()+1),1);	  
	}
    }
  
  if (vSim_JEMEtSums.size()!=0)
    {
      mLog<<MSG::DEBUG<<"JEMEtSums: additional Sim data for"<<endreq;

      //fill errorcounter
      for( it_Sim_JEMEtSums  = vSim_JEMEtSums.begin(); it_Sim_JEMEtSums <  vSim_JEMEtSums. end(); ++it_Sim_JEMEtSums )
	{
	  mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_JEMEtSums).crate()<<" Module "<<(*it_Sim_JEMEtSums).module()<<endreq;
	  mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_JEMEtSums).Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_JEMEtSums).Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_JEMEtSums).Et()<<endreq;
	  
	  //m_h_SimBSMon_JEP->Fill(1,(*it_Sim_JEMEtSums).crate()*19 + ((*it_Sim_JEMEtSums).module()+1),1);
	}
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
  CMMJetHitsCollection::const_iterator it_CMMJetHits ;
          
  //retrieve JEMHits from storegate for comparison with transmitted data stored in CMMJetHits
  const JEMHitsCollection* JEMHits;
  //JEMHitsCollection::const_iterator it_JEMHits ;
  sc = m_storeGate->retrieve(JEMHits, m_BS_JEMHitsLocation); 
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG:: INFO<< "No JEMHits found in TES at "<< m_BS_JEMHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  // compare JEMHits-vector with CMMJEMHits-vector, erase corresponding entries
  // -> the remaining vectors contain faulty transmissions! 
  // use standard vector instead of datavector for transmissioncheck:
  // datavector erases not only the pointer  in the vector, but also the referenced object
  // -> segmentation fault!
  std::vector <LVL1::CMMJetHits>  vCMMJEMHits;
  std::vector <LVL1::JEMHits>  vJEMHits;
  int CrateCMMMainHits=0;
  int CrateCMMFwdHits=0;
  int SystemRemoteCMMMainHits=0;
  int SystemRemoteCMMFwdHits=0;
  //int noMatchfound;
  
  // put CMM input data (JEMHits) into CMMJEMHits vector
  for( it_CMMJetHits  = CMMJetHits ->begin(); it_CMMJetHits < CMMJetHits -> end(); ++it_CMMJetHits )
    {
      // JEM information for transmission check JEMs -> CMMs 
      if ((*it_CMMJetHits)->dataID()<16)
	{
	  vCMMJEMHits.push_back(**it_CMMJetHits );
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
      vJEMHits.push_back(**it_JEMHits);
    }
  
  bool found;
  std::vector <LVL1::CMMJetHits>::iterator it_vCMMJEMHits;
  std::vector <LVL1::JEMHits>::iterator it_vJEMHits;
  //int noMatchfound;

  //---------------------------------- backplane transmission JEMs -> CMM -----------------------------
  it_vCMMJEMHits=vCMMJEMHits.begin();
  // step through both vectors and compare...
  mLog<<MSG::DEBUG<<"JEMHits Transmission error JEM -> CMM"<<endreq;
  while (it_vCMMJEMHits<vCMMJEMHits.end())
    {
      found = 0;
      noMatchfound=0;
      it_vJEMHits=vJEMHits.begin();
      
      while ((found==0)and(it_vJEMHits<vJEMHits.end()))
	{
	  if (((*it_vCMMJEMHits).crate()==(*it_vJEMHits).crate())
	      and((*it_vCMMJEMHits).dataID()==(*it_vJEMHits).module()))
	    {
	      if ((*it_vCMMJEMHits).Hits()!=(*it_vJEMHits).JetHits())
		{
		  mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMHits).crate()<<" Module "<<(*it_vJEMHits).module()<<endreq;
		  mLog<<MSG::VERBOSE<<"JEM Hit "<<Help->Binary((*it_vJEMHits).JetHits(),24)<<endreq;
		  mLog<<MSG::VERBOSE<<"CMM Hit "<<Help->Binary((*it_vCMMJEMHits).Hits(),24)<<endreq;
		  
		  noMatchfound=1;
		}
	      // if CMMJEMHits and JEMHits didn't match, a "1" is entered at the specific module;
	      // if a match is found, a "0" is entered (just to see that the histogram is filled
	      // since there hopefully are only zeros in there)
	      m_h_TransCheck_JEP->Fill(1,(*it_vJEMHits).crate()*20+(*it_vJEMHits).module()+1,noMatchfound);
	      
	      vCMMJEMHits.erase(it_vCMMJEMHits);
	      vJEMHits.erase(it_vJEMHits);
	      found=1;
	    }
	  else it_vJEMHits=it_vJEMHits+1;
	}// end while (not found and not JEMHits.end)
      
      // only step further in CMMJEMHits if no corresponding entry has been found in JEMHits
      if (found==0)it_vCMMJEMHits=it_vCMMJEMHits+1;
    }  // end while (not CMMJEMHits.end)
  
  // if the CMMJEMHits vector isn't empty, fill the error counter!
  if (vCMMJEMHits.size()!=0)
    {
      mLog<<MSG::DEBUG<<vCMMJEMHits.size()<<"JEMHits Transmission: additional CMM information"<<endreq;
      for( it_vCMMJEMHits  = vCMMJEMHits.begin(); it_vCMMJEMHits <  vCMMJEMHits. end(); ++it_vCMMJEMHits )
	{
	  mLog<<MSG::DEBUG<<"Crate "<<(*it_vCMMJEMHits).crate()<<" Module "<<(*it_vCMMJEMHits).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"Hit "<<(*it_vCMMJEMHits).Hits()<<endreq;
	  
	  if ((*it_vCMMJEMHits).Hits()>0) m_h_TransCheck_JEP->Fill(1,(*it_vCMMJEMHits).crate()*20+(*it_vCMMJEMHits).dataID()+1,1);
	}
    }
  
  // if the JEMHits vector isn't empty, fill the error counter!
  if (vJEMHits.size()!=0)
    {
      mLog<<MSG::DEBUG<<vJEMHits.size()<<" JEMHits Transmission: addtional JEM information"<<endreq;
      for( it_vJEMHits  =  vJEMHits.begin(); it_vJEMHits <  vJEMHits. end(); ++it_vJEMHits )
	{
	  mLog<<MSG::DEBUG<<"Crate "<<(*it_vJEMHits).crate()<<" Module "<<(*it_vJEMHits).module()<<endreq;
	  mLog<<MSG::VERBOSE<<"Hit "<<(*it_vJEMHits).JetHits()<<endreq;
	  
	  if ((*it_vJEMHits).JetHits()>0) m_h_TransCheck_JEP->Fill(1,(*it_vJEMHits).crate()*20+(*it_vJEMHits).module()+1,1);
	}
    }


  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            CMMJetHits(JEMDAQ.JEMJetHits) = (Crate/System)CMMJetDAQ.CrateCMMJetHits
  // ---------------------------------------------------------------------------------------------

  mLog << MSG::DEBUG << "==== JEM Jet Hits -> (Crate/System) CMM Jet Hits: Functionality ===="<< endreq ;

  const CMMJetHitsCollection* BS_CMMJetHits;
  sc = m_storeGate->retrieve(BS_CMMJetHits, m_BS_CMMJetHitsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMMJetHits found in TES at "<< m_BS_CMMJetHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  const CMMJetHitsCollection* Sim_CMMJetHits;
  sc = m_storeGate->retrieve(Sim_CMMJetHits, m_Sim_CMMJetHitsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMMJetHits found in TES at "<< m_Sim_CMMJetHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  std::vector <LVL1::CMMJetHits>  vBS_CMMJetHits;
  std::vector <LVL1::CMMJetHits>  vSim_CMMJetHits;
  //use standard vector instead of datavector for transmissioncheck:
  //datavector erases not only the pointer  in the vector, but also the referenced object
  //-> segmentation fault!
  
  //fill vectors
  //CMMJetHitsCollection::const_iterator it_CMMJetHits ;
  
  for( it_CMMJetHits  = BS_CMMJetHits->begin(); it_CMMJetHits < BS_CMMJetHits->end(); ++it_CMMJetHits )
    {
      if ((*it_CMMJetHits)->dataID()>15)
	{
	  vBS_CMMJetHits.push_back(**it_CMMJetHits);	       
	}
    }   
  for( it_CMMJetHits  = Sim_CMMJetHits->begin(); it_CMMJetHits < Sim_CMMJetHits->end(); ++it_CMMJetHits )
    {
      if ((*it_CMMJetHits)->dataID()>15)
	{
	  vSim_CMMJetHits.push_back(**it_CMMJetHits);	
	}       
    }   
  
  // Step over all cells and compare
  std::vector <LVL1::CMMJetHits>::iterator it_BS_CMMJetHits;
  std::vector <LVL1::CMMJetHits>::iterator it_Sim_CMMJetHits;
  
  it_BS_CMMJetHits=vBS_CMMJetHits.begin();
  
  mLog<<MSG::DEBUG<<"CMMJetHits Calculation difference for"<<endreq;
  while (it_BS_CMMJetHits<vBS_CMMJetHits.end())
    {
      foundModule = 0;
      noMatchfound=0;
      it_Sim_CMMJetHits=vSim_CMMJetHits.begin();
      
      while ((foundModule==0)and(it_Sim_CMMJetHits<vSim_CMMJetHits.end()))
	{
	  if (((*it_BS_CMMJetHits).crate()==(*it_Sim_CMMJetHits).crate())
	      and((*it_BS_CMMJetHits).dataID()==(*it_Sim_CMMJetHits).dataID()))
	    {
	      foundModule=1;
	      
	      if ((*it_BS_CMMJetHits).Hits()!=(*it_Sim_CMMJetHits).Hits())
		{
		  mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_CMMJetHits).crate()<<" DataId "<<(*it_BS_CMMJetHits).dataID()<<endreq;
		  mLog<<MSG::VERBOSE<<"BS: Hits"<<Help->Binary((*it_BS_CMMJetHits).Hits(),24)<<endreq;
		  mLog<<MSG::VERBOSE<<"Sim: Hits"<<Help->Binary((*it_Sim_CMMJetHits).Hits(),24)<<endreq;
		  
		  noMatchfound=1;
		}
	      m_h_SimBSMon_JEP->Fill(2,((*it_BS_CMMJetHits).crate()*19+16+1),noMatchfound);
	      
	      vBS_CMMJetHits.erase(it_BS_CMMJetHits);
	      vSim_CMMJetHits.erase(it_Sim_CMMJetHits);
	    }
	  else it_Sim_CMMJetHits=it_Sim_CMMJetHits+1;
	}
      if (foundModule==0)it_BS_CMMJetHits=it_BS_CMMJetHits+1;
    }
  
  if (vBS_CMMJetHits.size()!=0)
    {
      mLog<<MSG::DEBUG<<"CMMJetHits: additional BS data for"<<endreq;
      
      //fill errorcounter
      for( it_BS_CMMJetHits  = vBS_CMMJetHits.begin(); it_BS_CMMJetHits <  vBS_CMMJetHits. end(); ++it_BS_CMMJetHits )
	{
	  mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_CMMJetHits).crate()<<" DataId "<<(*it_BS_CMMJetHits).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"BS: Hits"<<Help->Binary((*it_BS_CMMJetHits).Hits(),24)<<endreq;
	  
	  m_h_SimBSMon_JEP->Fill(2,((*it_BS_CMMJetHits).crate()*19+16+1),1);
	}
    }
  
  if (vSim_CMMJetHits.size()!=0)
    {
      mLog<<MSG::DEBUG<<"CMMJetHits: additional Sim data "<<endreq;
      
      //fill errorcounter
      for( it_Sim_CMMJetHits  = vSim_CMMJetHits.begin(); it_Sim_CMMJetHits <  vSim_CMMJetHits. end(); ++it_Sim_CMMJetHits )
	{
	  mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_CMMJetHits).crate()<<" DataId "<<(*it_Sim_CMMJetHits).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"Sim: Hits"<<Help->Binary((*it_Sim_CMMJetHits).Hits(),24)<<endreq;
	  
	  m_h_SimBSMon_JEP->Fill(2,((*it_Sim_CMMJetHits).crate()*19+16+1),1);
	}
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
  //JEMEtSumsCollection::const_iterator it_JEMEtSums ;
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
  std::vector <LVL1::CMMEtSums>  vCMMJEMEtSums;
  std::vector <LVL1::JEMEtSums>  vJEMEtSums;
  int CrateEx=0, CrateEy=0, CrateEt=0;
  int SystemEx=0, SystemEy=0, SystemEt=0;
  noMatchfound=0;
  for( it_CMMEtSums  = CMMEtSums ->begin(); it_CMMEtSums < CMMEtSums -> end(); ++it_CMMEtSums )
    {	  
      if ((*it_CMMEtSums)->dataID()<16)
	{
	  vCMMJEMEtSums.push_back(**it_CMMEtSums );
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
  for( it_JEMEtSums  =  JEMEtSums->begin(); it_JEMEtSums <  JEMEtSums-> end(); ++it_JEMEtSums )
    {
      vJEMEtSums.push_back(**it_JEMEtSums);
    }
  
  std::vector <LVL1::CMMEtSums>::iterator it_vCMMJEMEtSums;
  std::vector <LVL1::JEMEtSums>::iterator it_vJEMEtSums;
  
  //---------------------------------- backplane transmission JEMs -> crate CMM -----------------------------
  it_vCMMJEMEtSums=vCMMJEMEtSums.begin();
  mLog<<MSG::DEBUG<<"CMMEtSums Transmission errors JEM -> CMM"<<endreq;
  while (it_vCMMJEMEtSums<vCMMJEMEtSums.end())
    {
      bool found = 0;
      noMatchfound = 0;
      it_vJEMEtSums=vJEMEtSums.begin();
      
      while ((found==0)and(it_vJEMEtSums<vJEMEtSums.end()))
	{
	  if (((*it_vCMMJEMEtSums).crate()==(*it_vJEMEtSums).crate())
	      and((*it_vCMMJEMEtSums).dataID()==(*it_vJEMEtSums).module()))
	    {
	      if (((*it_vCMMJEMEtSums).Ex()!=(*it_vJEMEtSums).Ex())
		  or ((*it_vCMMJEMEtSums).Ey()!=(*it_vJEMEtSums).Ey())
		  or ((*it_vCMMJEMEtSums).Et()!=(*it_vJEMEtSums).Et()))
		{
		  mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMEtSums).crate()<<" Module "<<(*it_vJEMEtSums).module()<<endreq;
		  mLog<<MSG::VERBOSE<<"JEM Ex "<<(*it_vJEMEtSums).Ex()<<endreq;
		  mLog<<MSG::VERBOSE<<"JEM Ey "<<(*it_vJEMEtSums).Ey()<<endreq;
		  mLog<<MSG::VERBOSE<<"JEM Et "<<(*it_vJEMEtSums).Et()<<endreq;
		  
		  mLog<<MSG::VERBOSE<<"CMM Ex "<<(*it_vCMMJEMEtSums).Ex()<<endreq;
		  mLog<<MSG::VERBOSE<<"CMM Ey "<<(*it_vCMMJEMEtSums).Ey()<<endreq;
		  mLog<<MSG::VERBOSE<<"CMM Et "<<(*it_vCMMJEMEtSums).Et()<<endreq;
		  
		  noMatchfound=1;
		}
	      // if CMMJEMEtSums and JEMEtSums didn't match, a "1" is entered at the specific module;
	      // if a match is found, a "0" is entered (just to see that the histogram is filled
	      // since there hopefully are only zeros in there)
	      
	      m_h_TransCheck_JEP->Fill(2,((*it_vJEMEtSums).crate()*20+(*it_vJEMEtSums).module() + 1),noMatchfound);
	      
	      vCMMJEMEtSums.erase(it_vCMMJEMEtSums);
	      vJEMEtSums.erase(it_vJEMEtSums);
	      found=1;
	    }
	  else it_vJEMEtSums=it_vJEMEtSums+1;
	}
      if (found==0)it_vCMMJEMEtSums=it_vCMMJEMEtSums+1;
    }
  
  // if the CMMJEMEtSums vector isn't empty, fill the error counter!
  if (vCMMJEMEtSums.size()!=0)
    {
      mLog<<MSG::DEBUG<<vCMMJEMEtSums.size()<<"CMMEtSums Transmission: additional CMM information"<<endreq;
      for( it_vCMMJEMEtSums  = vCMMJEMEtSums.begin(); it_vCMMJEMEtSums <  vCMMJEMEtSums. end(); ++it_vCMMJEMEtSums )
	{
	  mLog<<MSG::DEBUG<<" Crate "<<(*it_vCMMJEMEtSums).crate()<<" Module "<<(*it_vCMMJEMEtSums).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Ex "<<(*it_vCMMJEMEtSums).Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Ey "<<(*it_vCMMJEMEtSums).Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Et "<<(*it_vCMMJEMEtSums).Et()<<endreq;
	  
	  m_h_TransCheck_JEP->Fill(2,(*it_vCMMJEMEtSums).crate()*20+(*it_vCMMJEMEtSums).dataID() + 1,noMatchfound);
	}
    }
  
  // if the JEMEtSums vector isn't empty, fill the error counter!
  if (vJEMEtSums.size()!=0)
    {
      mLog<<MSG::DEBUG<<vJEMEtSums.size()<<"CMMEtSums Transmission: addtional JEM information"<<endreq;
      for( it_vJEMEtSums  =  vJEMEtSums.begin(); it_vJEMEtSums <  vJEMEtSums. end(); ++it_vJEMEtSums )
	{
	  mLog<<MSG::DEBUG<<"Crate "<<(*it_vJEMEtSums).crate()<<" Module "<<(*it_vJEMEtSums).module()<<endreq;
	  mLog<<MSG::VERBOSE<<"JEM Ex "<<(*it_vJEMEtSums).Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"JEM Ey "<<(*it_vJEMEtSums).Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"JEM Et "<<(*it_vJEMEtSums).Et()<<endreq;
	  
	  m_h_TransCheck_JEP->Fill(2,(*it_vJEMEtSums).crate()*20+(*it_vJEMEtSums).module() + 1,noMatchfound);
	}
    }
  

  // ---------------------------------------------------------------------------------------------
  // Check	Functionality
  //            CMMEnergySums(JEMDAQ.JEMEnergySums) = (Crate/System)CMMEnergyDAQ.CrateCMMEnergySums
  // ---------------------------------------------------------------------------------------------
  
  mLog << MSG::DEBUG << "==== JEMEnergySums -> (Crate/System) CMM Energy Sums: Functionality ===="<< endreq ;

  const CMMEtSumsCollection* BS_CMMEtSums;
  sc = m_storeGate->retrieve(BS_CMMEtSums, m_BS_CMMEtSumsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMMEtSums found in TES at "<< m_BS_CMMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  const CMMEtSumsCollection* Sim_CMMEtSums;
  sc = m_storeGate->retrieve(Sim_CMMEtSums, m_Sim_CMMEtSumsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMMEtSums found in TES at "<< m_Sim_CMMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  std::vector <LVL1::CMMEtSums>  vBS_CMMEtSums;
  std::vector <LVL1::CMMEtSums>  vSim_CMMEtSums;
  //use standard vector instead of datavector for transmissioncheck:
  //datavector erases not only the pointer  in the vector, but also the referenced object
  //-> segmentation fault!
  
  //fill vectors
  //CMMEtSumsCollection::const_iterator it_CMMEtSums ;
  
  for( it_CMMEtSums  = BS_CMMEtSums->begin(); it_CMMEtSums < BS_CMMEtSums->end(); ++it_CMMEtSums )
    {
      if ((*it_CMMEtSums)->dataID()>15)
	{
	  vBS_CMMEtSums.push_back(**it_CMMEtSums);	       
	}
    }   
  for( it_CMMEtSums  = Sim_CMMEtSums->begin(); it_CMMEtSums < Sim_CMMEtSums->end(); ++it_CMMEtSums )
    {
      if ((*it_CMMEtSums)->dataID()>15)
	{
	  vSim_CMMEtSums.push_back(**it_CMMEtSums);	
	}       
    }   
  
  // Step over all cells and compare
  std::vector <LVL1::CMMEtSums>::iterator it_BS_CMMEtSums;
  std::vector <LVL1::CMMEtSums>::iterator it_Sim_CMMEtSums;
  
  it_BS_CMMEtSums=vBS_CMMEtSums.begin();
  
  mLog<<MSG::DEBUG<<"CMMEtSums Calculation difference for"<<endreq;
  while (it_BS_CMMEtSums<vBS_CMMEtSums.end())
      {
	foundModule = 0;
	noMatchfound=0;
	it_Sim_CMMEtSums=vSim_CMMEtSums.begin();
	
	while ((foundModule==0)and(it_Sim_CMMEtSums<vSim_CMMEtSums.end()))
	  {
	    if (((*it_BS_CMMEtSums).crate()==(*it_Sim_CMMEtSums).crate())
		and((*it_BS_CMMEtSums).dataID()==(*it_Sim_CMMEtSums).dataID()))
	      {
		foundModule=1;
		
		if (((*it_BS_CMMEtSums).Ex()!=(*it_Sim_CMMEtSums).Ex())
		    or ((*it_BS_CMMEtSums).Ey()!=(*it_Sim_CMMEtSums).Ey())
		    or ((*it_BS_CMMEtSums).Et()!=(*it_Sim_CMMEtSums).Et()))
		  {
		    mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_CMMEtSums).crate()<<" DataId "<<(*it_BS_CMMEtSums).dataID()<<endreq;
		    mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_CMMEtSums).Ex()<<endreq;
		    mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_CMMEtSums).Ey()<<endreq;
		    mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_CMMEtSums).Et()<<endreq;
		    
		    mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_CMMEtSums).Ex()<<endreq;
		    mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_CMMEtSums).Ey()<<endreq;
		    mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_CMMEtSums).Et()<<endreq;
		    
		    noMatchfound=1;    
		  }
		m_h_SimBSMon_JEP->Fill(1,(*it_BS_CMMEtSums).crate()*19 + 16 + 1,noMatchfound);
		
		vBS_CMMEtSums.erase(it_BS_CMMEtSums);
		vSim_CMMEtSums.erase(it_Sim_CMMEtSums);
	      }
	    else it_Sim_CMMEtSums=it_Sim_CMMEtSums+1;
	  }
	if (foundModule==0)it_BS_CMMEtSums=it_BS_CMMEtSums+1;
      }
  
  if (vBS_CMMEtSums.size()!=0)
    {
      mLog<<MSG::DEBUG<<"CMMEtSums: additional BS data for"<<endreq;
      
      //fill errorcounter
      for( it_BS_CMMEtSums  = vBS_CMMEtSums.begin(); it_BS_CMMEtSums <  vBS_CMMEtSums. end(); ++it_BS_CMMEtSums )
	{
	  mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_CMMEtSums).crate()<<" DataId "<<(*it_BS_CMMEtSums).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_CMMEtSums).Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_CMMEtSums).Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_CMMEtSums).Et()<<endreq;
	  
	  m_h_SimBSMon_JEP->Fill(1,(*it_BS_CMMEtSums).crate()*19 + 16 + 1,1);
	  
	}
    }
  
  if (vSim_CMMEtSums.size()!=0)
    {
      mLog<<MSG::DEBUG<<"CMMEtSums: additional Sim data for"<<endreq;
      
      //fill errorcounter
      for( it_Sim_CMMEtSums  = vSim_CMMEtSums.begin(); it_Sim_CMMEtSums <  vSim_CMMEtSums. end(); ++it_Sim_CMMEtSums )
	{
	  mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_CMMEtSums).crate()<<" DataId "<<(*it_Sim_CMMEtSums).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_CMMEtSums).Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_CMMEtSums).Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_CMMEtSums).Et()<<endreq;

	  m_h_SimBSMon_JEP->Fill(1,(*it_Sim_CMMEtSums).crate()*19 + 16 + 1,1);
	}
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
  m_h_TransCheck_JEP->Fill(1,18,noMatchfound);


  // ---------------------------------------------------------------------------------------------
  // Check	Funkionalität
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


  // ---------------------------------------------------------------------------------------------
  // Check	Funkionalität
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
  // Check	Funktionality
  //            JEMJetRoI(JEMDAQ.JetElements) = JEPRoI.JEMRoIJetHits
  // ---------------------------------------------------------------------------------------------

  mLog << MSG::DEBUG << "==== JEM RoI: Funktionality ===="<< endreq ;

  const JemRoiCollection* BS_JEMRoI;
  sc = m_storeGate->retrieve(BS_JEMRoI, m_BS_JEMRoILocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMRoI found in TES at "<< m_BS_JEMRoILocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  const JemRoiCollection* Sim_JEMRoI;
  sc = m_storeGate->retrieve(Sim_JEMRoI, m_Sim_JEMRoILocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMRoI found in TES at "<< m_Sim_JEMRoILocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  std::vector <LVL1::JEMRoI>  vBS_JEMRoI;
  std::vector <LVL1::JEMRoI>  vSim_JEMRoI;
  //use standard vector instead of datavector for transmissioncheck:
  //datavector erases not only the pointer  in the vector, but also the referenced object
  //-> segmentation fault!
  
  //fill vectors
  JemRoiCollection::const_iterator it_JEMRoI ;
  
  for( it_JEMRoI  = BS_JEMRoI->begin(); it_JEMRoI < BS_JEMRoI->end(); ++it_JEMRoI )
    {
      vBS_JEMRoI.push_back(**it_JEMRoI);	       
    }   
  for( it_JEMRoI  = Sim_JEMRoI->begin(); it_JEMRoI < Sim_JEMRoI->end(); ++it_JEMRoI )
    {
      vSim_JEMRoI.push_back(**it_JEMRoI);	       
    }   
  
  // Step over all cells and compare
  std::vector <LVL1::JEMRoI>::iterator it_BS_JEMRoI;
  std::vector <LVL1::JEMRoI>::iterator it_Sim_JEMRoI;
  
  it_BS_JEMRoI=vBS_JEMRoI.begin();
  
  mLog<<MSG::DEBUG<<"JEMRoI Calculation difference for"<<endreq;
  while (it_BS_JEMRoI<vBS_JEMRoI.end())
    {
      foundModule = 0;
      noMatchfound=0;
      it_Sim_JEMRoI=vSim_JEMRoI.begin();
      
      while ((foundModule==0)and(it_Sim_JEMRoI<vSim_JEMRoI.end()))
	{
	  
	  if (((*it_BS_JEMRoI).crate()==(*it_Sim_JEMRoI).crate())
	      and((*it_BS_JEMRoI).jem()==(*it_Sim_JEMRoI).jem())
	      and ((*it_BS_JEMRoI).frame()==(*it_Sim_JEMRoI).frame())
	      and ((*it_BS_JEMRoI).location()==(*it_Sim_JEMRoI).location())
)
	    {
	      foundModule=1;
	      
	      if ((*it_BS_JEMRoI).hits()!=(*it_Sim_JEMRoI).hits())
		{
		  mLog<<MSG::DEBUG<<"Crate "<<(*it_BS_JEMRoI).crate()<<" Module "<<(*it_BS_JEMRoI).jem()<<" frame "<<(*it_BS_JEMRoI).frame() <<" location "<<(*it_BS_JEMRoI).location() <<endreq;
		  mLog<<MSG::VERBOSE<<"BS:  RoI "<<Help->Binary((*it_BS_JEMRoI).hits(),8)<<endreq;
		  mLog<<MSG::VERBOSE<<"Sim: RoI "<<Help->Binary((*it_Sim_JEMRoI).hits(),8)<<endreq;
		  
		  noMatchfound=1;
		}
	      m_h_SimBSMon_JEP->Fill(3,((*it_BS_JEMRoI).crate()*19+(*it_BS_JEMRoI).jem()+1),1);
	      
	      vBS_JEMRoI.erase(it_BS_JEMRoI);
	      vSim_JEMRoI.erase(it_Sim_JEMRoI);
	    }
	  else it_Sim_JEMRoI=it_Sim_JEMRoI+1;
	}
      if (foundModule==0)it_BS_JEMRoI=it_BS_JEMRoI+1;
    }
  
  if (vBS_JEMRoI.size()!=0)
    {
      mLog<<MSG::DEBUG<<"JEMRoI: additional BS data for "<<endreq;
      
      //fill errorcounter
      for( it_BS_JEMRoI  = vBS_JEMRoI.begin(); it_BS_JEMRoI <  vBS_JEMRoI. end(); ++it_BS_JEMRoI )
	{
	  mLog<<MSG::DEBUG<<"BS: Crate "<<(*it_BS_JEMRoI).crate()<<" Module "<<(*it_BS_JEMRoI).jem()<<" frame "<<(*it_BS_JEMRoI).frame() <<" location "<<(*it_BS_JEMRoI).location() <<endreq;
	  mLog<<MSG::VERBOSE<<"BS: RoI "<<Help->Binary((*it_BS_JEMRoI).hits(),8)<<endreq;
	  
	  m_h_SimBSMon_JEP->Fill(3,((*it_BS_JEMRoI).crate()*19+(*it_BS_JEMRoI).jem()+1),1);
	}
    }
  
  if (vSim_JEMRoI.size()!=0)
    {
      mLog<<MSG::DEBUG<<"JEMRoI: additional Sim data for"<<endreq;
      
      //fill errorcounter
      for( it_Sim_JEMRoI  = vSim_JEMRoI.begin(); it_Sim_JEMRoI <  vSim_JEMRoI. end(); ++it_Sim_JEMRoI )
	{
	  mLog<<MSG::DEBUG<<"Sim: Crate "<<(*it_Sim_JEMRoI).crate()<<" Module "<<(*it_Sim_JEMRoI).jem()<<" frame "<<(*it_Sim_JEMRoI).frame() <<" location "<<(*it_Sim_JEMRoI).location() <<endreq;
	  mLog<<MSG::VERBOSE<<"Sim: RoI "<<Help->Binary((*it_Sim_JEMRoI).hits(),8)<<endreq;
	  
	  //m_h_SimBSMon_JEP->Fill(3,((*it_Sim_JEMRoI).crate()*19+(*it_Sim_JEMRoI).jem()+1),1);
	}
    }
  
  // ================= CMM RoI ========================================================
  
  mLog << MSG::DEBUG << "==== CMM RoI: Funktionality ===="<< endreq ;

  // retrieve RoI information from Storegate
  LVL1::CMMRoI* BS_CR = new LVL1::CMMRoI ;
  sc = m_storeGate->retrieve (BS_CR, m_BS_CMMRoILocation);
  if (sc==StatusCode::FAILURE)
    {
      mLog <<MSG::INFO<<"No BS CMM RoI found in TES at "<< m_BS_CMMRoILocation<<endreq;
      return StatusCode::SUCCESS;    
    }

  LVL1::CMMRoI* Sim_CR = new LVL1::CMMRoI ;
  sc = m_storeGate->retrieve (Sim_CR, m_Sim_CMMRoILocation);
  if (sc==StatusCode::FAILURE)
    {
      mLog <<MSG::INFO<<"No Sim CMM RoI found in TES at "<< m_Sim_CMMRoILocation<<endreq;
      return StatusCode::SUCCESS;    
    }

  noMatchfound=0;  
  if ( ((*BS_CR).jetEtHits()!=(*Sim_CR).jetEtHits())or((*BS_CR).sumEtHits()!=(*Sim_CR).sumEtHits())or((*BS_CR).missingEtHits()!=(*Sim_CR).missingEtHits()))  
    {
      mLog<<MSG::DEBUG<<"CMMRoI: No Match found between BS and Sim for JetEtHits, SumEtHits or MissingEtHits"<<endreq;
      noMatchfound=1; 
    }
  m_h_SimBSMon_JEP->Fill(3,(19+16+1),noMatchfound);


  noMatchfound=0; 
  if ( ((*BS_CR).ex()!=(*Sim_CR).ex())or((*BS_CR).ey()!=(*Sim_CR).ey())or((*BS_CR).et()!=(*Sim_CR).et()))  
    {
      mLog<<MSG::DEBUG<<"CMMRoI: No Match found between BS and Sim for Ex, Ey or Et"<<endreq;
      noMatchfound=1; 
    }
  m_h_SimBSMon_JEP->Fill(3,(19+16+1),noMatchfound);
  
  


  // =============================================================================================
  // ================= JEP -> CTP ================= ==============================================
  // =============================================================================================

  // ---------------------------------------------------------------------------------------------
  // Check	Transmission
  //              SystemCMMJetDAQ.SystemCMMJetHits = JEPCTP.SystemCMMJetHits
  //              SystemCMMEnergyDAQ.TotalEtMissHits. = JEPCTP.TotalEtMissHits
  //              SystemCMMEnergyDAQ.TotalEtSumHits. = JEPCTP.TotalEtSumHits
  // ---------------------------------------------------------------------------------------------

  return StatusCode( StatusCode::SUCCESS );
}

/*---------------------------------------------------------*/
StatusCode JEPTransPerfMon::procHistograms( bool isEndOfEventsBlock, 
				  bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock ) 
    {  
    }
 
  if(m_EventNoInHisto==1)
    {      
      if( isEndOfRun ) 
	{
	  std::stringstream buffer;
	  buffer.str("");
	  buffer<<m_NoEvents;
	  std::string title;
	  
	  title = m_h_SimBSMon_JEP-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_SimBSMon_JEP->SetTitle(title.c_str());
	  
	  title = m_h_TransCheck_JEP-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TransCheck_JEP->SetTitle(title.c_str());
	  
	  title = m_h_TransCheck_emJetElements-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TransCheck_emJetElements->SetTitle(title.c_str());
	  
	  title = m_h_TransCheck_hadJetElements-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TransCheck_hadJetElements->SetTitle(title.c_str());
	}
  }
  
  return StatusCode( StatusCode::SUCCESS );
}

/*---------------------------------------------------------*/
void JEPTransPerfMon::TimeSliceMatch(int k, int TT_TS, const JECollection* TT_jetElements, int JE_TS, const JECollection* jetElements, MsgStream::MsgStream* mLog)
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

      *mLog <<MSG::VERBOSE<<"JE Crate: "<<crate<<" Module: "<<module<<endreq;

      
      while ((JEFound==0)and(it_TT_je < TT_jetElements->end()))
	{
	  if (((*it_TT_je)->phi()==(*it_je)->phi())and((*it_TT_je)->eta()==(*it_je)->eta()))
	    {
	      JEFound=1;
	      
	      if ((*it_TT_je)->emEnergyVec()[TT_TS]==(*it_je)->emEnergyVec()[JE_TS])
		{
		  em_NoJEMatchFound=0; 
		  //*mLog << MSG::DEBUG << "em Match Found eta: "<<(*it_TT_je)->eta()<<" phi: " <<(*it_TT_je)->phi()<<" energy HW: " << (*it_je)->emEnergyVec()[JE_TS]<<" SW: " << (*it_TT_je)->emEnergyVec()[TT_TS]<< endreq ;

		}
	      
	      if ((*it_TT_je)->hadEnergyVec()[TT_TS]==(*it_je)->hadEnergyVec()[JE_TS])
		{
		  had_NoJEMatchFound=0;     

		}	      
	    }
	  if ((JEFound==0)and(it_TT_je < TT_jetElements->end())) it_TT_je ++;
	}

      if (em_NoJEMatchFound==1 and JEFound==1)
	{





      Identifier EmTowerId,HadTowerId;

      int detside= m_l1CaloTTIdTools->pos_neg_z((*it_TT_je)->eta());
      int detregion = m_l1CaloTTIdTools->regionIndex((*it_TT_je)->eta());
      int eta=m_l1CaloTTIdTools->etaIndex((*it_TT_je)->eta());
      int phi=m_l1CaloTTIdTools->phiIndex((*it_TT_je)->eta(), (*it_TT_je)->phi());
      
      //---------------------------- EM Energy -----------------------------
      EmTowerId = m_lvl1Helper->tower_id(detside, 0, detregion,eta,phi );  
      HadTowerId = m_lvl1Helper->tower_id(detside, 1, detregion,eta,phi );  
      int crate, module;
      try 
	{
	  HWIdentifier ttOnlId = m_ttSvc->createTTChannelID(EmTowerId);
	  crate     = m_l1ttonlineHelper->crate(ttOnlId);
	  module    = m_l1ttonlineHelper->module(ttOnlId);
	  //mLog << MSG::DEBUG << "em PPM crate: " << crate<<"  module: "<<module << endreq ;
	} 
      catch(LArID_Exception& except) 
	{
	  *mLog << MSG::ERROR << "LArID_Exception " << (std::string) except << endreq ;
	}









	  *mLog << MSG::VERBOSE << "no em Match Found eta: "<<(*it_TT_je)->eta()<<" phi: " <<(*it_TT_je)->phi()<<" energy HW: " << (*it_je)->emEnergyVec()[JE_TS]<<" SW: " << (*it_TT_je)->emEnergyVec()[TT_TS]<< endreq ;
	  *mLog << MSG::VERBOSE << "PPM data from Crate: "<<crate<<" Module: " <<module<< endreq ;

	}
      if (had_NoJEMatchFound==1 and JEFound==1)
	{
	  *mLog << MSG::VERBOSE << "no had Match Found eta: "<<(*it_TT_je)->eta()<<" phi: " <<(*it_TT_je)->phi()<<" energy HW: " << (*it_je)->hadEnergyVec()[JE_TS]<<" SW: " << (*it_TT_je)->hadEnergyVec()[TT_TS]<< endreq ;
	}

      if (JEFound==1)
	{
	  //*mLog <<MSG::VERBOSE<<"JEFound=1"<<endreq;
	  //*mLog << MSG::VERBOSE << "em NoJEMatchFound "<<em_NoJEMatchFound<<" Module: " << module<<" crate: "<<crate <<" energy HW: " << (*it_je)->emEnergyVec()[JE_TS]<<" SW: " << (*it_TT_je)->emEnergyVec()[TT_TS]<< endreq ;
	  
	  //*mLog << MSG::VERBOSE << "had NoJEMatchFound "<<had_NoJEMatchFound<<" Module: " << module<<" crate: "<<crate <<" energy HW: " << (*it_je)->hadEnergyVec()[JE_TS]<<" SW: " << (*it_TT_je)->hadEnergyVec()[TT_TS]<< endreq ;
	  *mLog <<MSG::VERBOSE<<endreq;
	}
      else
	{
	  *mLog <<MSG::VERBOSE<<"JEFound=0"<<endreq;
	  *mLog << MSG::VERBOSE << " Module: " << module<<" crate: "<<crate <<" eta "<<(*it_je)->eta()<<" phi "<< (*it_je)->phi()<<endreq ;
	  
	}
      m_h_TransCheck_emJetElements->Fill(k,(crate*18+module+1),em_NoJEMatchFound);
      m_h_TransCheck_hadJetElements->Fill(k,(crate*18+module+1),had_NoJEMatchFound);
    }
}

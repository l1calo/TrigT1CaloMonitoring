// ********************************************************************
//
// NAME:        CMMMon.cxx
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on CMM level
//
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

#include "TrigT1CaloMonitoring/CMMMon.h"
#include "TrigT1CaloMonitoring/JEMMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

//#include "TrigT1Calo/EnergyTrigger.h"
#include "TrigT1Calo/LVL1TriggerMenuDefs.h"
#include "TrigT1Calo/LVL1TriggerMenu.h"
#include "TrigT1Calo/InternalJetROI.h"
#include "TrigT1Calo/CMMRoI.h"
#include "TrigT1Calo/QuadLinear.h"
#include "TrigT1Calo/DataError.h"

#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "AthenaMonitoring/AthenaMonManager.h"

namespace LVL1 {
  class CMMRoI;
}


// *********************************************************************
// Public Methods
// *********************************************************************

/*---------------------------------------------------------*/
CMMMon::CMMMon( const std::string & type, const std::string & name,
		const IInterface* parent )
  : ManagedMonitorToolBase( type, name, parent )
/*---------------------------------------------------------*/
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "CMMJetHitsLocation", m_CMMJetHitsLocation =  LVL1::TrigT1CaloDefs::CMMJetHitsLocation) ;
  declareProperty( "CMMEtSumsLocation", m_CMMEtSumsLocation =  LVL1::TrigT1CaloDefs::CMMEtSumsLocation) ;  
  declareProperty( "CMMRoILocation", m_CMMRoILocation =  LVL1::TrigT1CaloDefs::CMMRoILocation) ;

  declareProperty( "JEMHitsLocation", m_JEMHitsLocation =  LVL1::TrigT1CaloDefs::JEMHitsLocation) ;
  declareProperty( "JEMEtSumsLocation", m_JEMEtSumsLocation=   LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;

  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/CMM") ;
  declareProperty( "ErrorPathInRootFile", m_ErrorPathInRootFile="Stats/L1Calo/Errors") ;
  declareProperty( "DataType", m_DataType="") ;
}


/*---------------------------------------------------------*/
CMMMon::~CMMMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode CMMMon::bookHistograms( bool isNewEventsBlock, 
				   bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in CMMMon::bookHistograms" << endreq;
  
  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  
  if (sc.isFailure()) 
    {
      mLog << MSG::ERROR
	   << "Unable to retrieve pointer to StoreGateSvc"
	   << endreq;
      return sc;
    }

  ManagedMonitorToolBase::LevelOfDetail_t LevelOfDetail=shift;
  if (m_DataType=="Sim") LevelOfDetail = expert;

  MonGroup CMM_DAQ ( this, (m_PathInRootFile+"/DAQ").c_str(), expert, eventsBlock );
  HistoBooker* DAQ_Booker = new HistoBooker(&CMM_DAQ, &mLog, m_DataType);

  MonGroup CMM_input ( this, (m_PathInRootFile + "/input").c_str(), expert, eventsBlock );
  HistoBooker* input_Booker = new HistoBooker(&CMM_input, &mLog, m_DataType);
  
  MonGroup CMM_RoI ( this, (m_PathInRootFile + "/RoI").c_str(), LevelOfDetail, eventsBlock );
  HistoBooker* RoI_Booker = new HistoBooker(&CMM_RoI, &mLog, m_DataType);
  
  MonGroup CMM_transmission ( this, (m_ErrorPathInRootFile ).c_str(), shift, eventsBlock );
  HistoBooker* transmission_Booker = new HistoBooker(&CMM_transmission, &mLog, "");
  
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }
  
  if( isNewEventsBlock || isNewLumiBlock ) 
    {	

      //----------------------------------  CMM Input data from JEMs -----------------------------
      m_h_CMMJetHits_JEM_MainHits=input_Booker->book1F("MainHits_CMM_input", "Main Jet Hit Multiplicity per Threshold  --  CMM input", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMJetHits_JEM_FwdHitsRight=input_Booker->book1F("FwdHitsRight_CMM_input", "Forward Right Jet Hit Multiplicity per Threshold  --  CMM input",4 , -0.5, 3.5, "Threshold No.", "#");
      m_h_CMMJetHits_JEM_FwdHitsLeft=input_Booker->book1F("FwdHitsLeft_CMM_input", "Forward Left Jet Hit Multiplicity per Threshold  --  CMM input", 4 , -0.5, 3.5,  "Threshold No.", "#");

      m_h_CMMEtSums_JEM_Ex=input_Booker->book1F("Ex_CMM_input", "CMM E_{x}  --  CMM input", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_JEM_Ey=input_Booker->book1F("Ey_CMM_input", "CMM E_{y}  --  CMM input", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_JEM_Et=input_Booker->book1F("Et_CMM_input", "CMM E_{t}  --  CMM input", 250, 0,250, "Ex [GeV]", "#");
      

      //---------------------------------- CMM output to DAQ -----------------------------
      m_h_CMMJetHits_MainJets = DAQ_Booker->book1F("TotalMainHits_CMM_DAQ", "Main Jet Hit Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMJetHits_FwdJetsRight = DAQ_Booker->book1F("TotalFwdHitsRight_CMM_DAQ", "Forward Right Jet Hit Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5, "Threshold No.", "#");
      m_h_CMMJetHits_FwdJetsLeft = DAQ_Booker->book1F("TotalFwdHitsLeft_CMM_DAQ", "Forward Left Jet Hit Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5,  "Threshold No.", "#");
      m_h_CMMJetHits_EtMap = DAQ_Booker->book1F("JetEtHits_CMM_DAQ", "JetEt Hit Multiplicity per Threshold  --  CMM DAQ", 4 ,-0.5, 3.5, "Threshold No.", "#");
      m_h_CMMEtSums_MissingEtMap = DAQ_Booker->book1F("MissingEtHits_CMM_DAQ", "MissingEt Hit Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMEtSums_SumEtMap = DAQ_Booker->book1F("SumEtHits_CMM_DAQ", "SumEt Hit Multiplicity per Threshold  --  CMM DAQ", 4, -0.5, 3.5, "Threshold No.", "#");

      m_h_CMMEtSums_Ex = DAQ_Booker->book1F("Ex_CMM_DAQ", "CMM E_{x}  --  CMM DAQ", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_Ey = DAQ_Booker->book1F("Ey_CMM_DAQ", "CMM E_{y}  --  CMM DAQ", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_Et = DAQ_Booker->book1F("Et_CMM_DAQ", "CMM E_{t}  --  CMM DAQ", 250, 0,250, "Ex [GeV]", "#");


      //---------------------------------- CMM output to RoI -----------------------------
      m_h_CMMRoI_JetEtHits =RoI_Booker->book1F("JetEtHits_CMM_RoI","JetEt Hit Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","#");
      m_h_CMMRoI_MissingEtHits =RoI_Booker->book1F("MissingEtHits_CMM_RoI","MissingEt Hit Multiplicity per Threshold  --  CMM RoI", 8, -0.5,7.5,"Threshold No.","#");
      m_h_CMMRoI_SumEtHits =RoI_Booker->book1F("SumEtHits_CMM_RoI","SumEt Hit Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","#");

      m_h_CMMRoI_Ex = RoI_Booker->book1F("Ex_CMM_RoI", "CMM E_{x}  --  CMM RoI", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMRoI_Ey = RoI_Booker->book1F("Ey_CMM_RoI", "CMM E_{y}  --  CMM RoI", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMRoI_Et = RoI_Booker->book1F("Et_CMM_RoI", "CMM E_{t}  --  CMM RoI", 250, 0,250, "Ex [GeV]", "#");


      if (m_DataType=="BS")
	{
	  //---------------------------------- S-Link errors -----------------------------
	  m_h_CMMJet_error=transmission_Booker->book2F("CMMJet_errors", "CMM Jet S-Link Errors per per Module and Crate",10,0.5,10.5,37,0.5,37.5,"","");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(1, "Parity");

	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(3, "GLinkParity");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(4, "GLinkProtocol");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(5, "BCNMismatch");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(6, "FIFOOverflow");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(7, "ModuleError");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(8, "GLinkDown");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(9, "GLinkTimeout");
	  m_h_CMMJet_error->GetXaxis()->SetBinLabel(10, "FailingBCN");

	  m_h_CMMEnergy_error=transmission_Booker->book2F("CMMEnergy_errors", "CMM Energy S-Link Errors per per Module and Crate",10,0.5,10.5,37,0.5,37.5,"","");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(1, "Parity");

	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(3, "GLinkParity");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(4, "GLinkProtocol");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(5, "BCNMismatch");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(6, "FIFOOverflow");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(7, "ModuleError");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(8, "GLinkDown");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(9, "GLinkTimeout");
	  m_h_CMMEnergy_error->GetXaxis()->SetBinLabel(10, "FailingBCN");

 	  std::string name;
	  std::stringstream buffer;
     
	  for (int i = 0; i < 16; i++)
	    {
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_CMMJet_error->GetYaxis()->SetBinLabel((i+1), name.c_str());
	      m_h_CMMEnergy_error->GetYaxis()->SetBinLabel((i+1), name.c_str());
	      
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_CMMJet_error->GetYaxis()->SetBinLabel((i+1+19), name.c_str());
	      m_h_CMMEnergy_error->GetYaxis()->SetBinLabel((i+1+19), name.c_str());
	    }
	  m_h_CMMJet_error->GetYaxis()->SetBinLabel(17, "C J CMM ");
	  m_h_CMMJet_error->GetYaxis()->SetBinLabel(18, "Crate 0: ");
	  m_h_CMMJet_error->GetYaxis()->SetBinLabel(36, "S J CMM ");
	  m_h_CMMJet_error->GetYaxis()->SetBinLabel(37, "Crate 1: ");

	  m_h_CMMEnergy_error->GetYaxis()->SetBinLabel(17, "C E CMM ");
	  m_h_CMMEnergy_error->GetYaxis()->SetBinLabel(18, "Crate 0: ");
	  m_h_CMMEnergy_error->GetYaxis()->SetBinLabel(36, "C E CMM ");
	  m_h_CMMEnergy_error->GetYaxis()->SetBinLabel(37, "Crate 1: ");
	

	  m_h_CMMRoI_error=transmission_Booker->book1F("CMMRoI_errors", "CMM RoI S-Link Parity and Overflow",8,0.5,8.5,"");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(1, "Parity (Ex)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(2, "Parity (Ey,SumEtMap)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(3, "Parity (Et,MissingEtMap)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(4, "Parity (JetEtMap)");

	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(6, "Overflow (Ex)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(7, "Overflow (Ey)");
	  m_h_CMMRoI_error->GetXaxis()->SetBinLabel(8, "Overflow (Et)");

	  //---------------------------------- Backplane transmission checks -----------------------------
	  m_h_TransCheck_JEP=transmission_Booker->book2F("JEP_TransCheck", "JEP Backplane Transmission Check JEM -> CMM per Module and Crate", 2,0.5,2.5,37,0.5,3.5, "", "");
	  //m_h_TransCheck_JEP-> SetOption ("text");
	  m_h_TransCheck_JEP->GetXaxis()->SetBinLabel(1, "Hits");
	  m_h_TransCheck_JEP->GetXaxis()->SetBinLabel(2, "Energy");
      
	  for (int i = 0; i < 16; i++)
	    {
	      buffer.str("");
	      buffer<<i;
	      
	      name = "JEM " + buffer.str();
	      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel((i+1), name.c_str());
	      
	      buffer.str("");
	      buffer<<i;
	      
	      name = "JEM " + buffer.str();
	      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel((i+1+20), name.c_str());
	    }
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(17, "C E CMM");
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(18, "C J CMM");
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(19, "Crate 0: ");
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(37, "Crate 1: ");
	}
    }
  if( isNewRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}


/*---------------------------------------------------------*/
StatusCode CMMMon::fillHistograms()
  /*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  Helper* Help = new Helper();
  
  // =============================================================================================
  // ================= Container: CMM Jet Hits ===================================================
  // =============================================================================================
  // retrieve CMM Jet Hits from Storegate
  const CMMJetHitsCollection* CMMJetHits;
  StatusCode sc = m_storeGate->retrieve(CMMJetHits, m_CMMJetHitsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMM JetHits found in TES at "  << m_CMMJetHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }  
  mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM Jet Hits ---------------"<<endreq;  
  CMMJetHitsCollection::const_iterator it_CMMJetHits ;
  
  // Step over all cells
  for( it_CMMJetHits  = CMMJetHits ->begin(); it_CMMJetHits < CMMJetHits -> end(); ++it_CMMJetHits )
    {	  
      //put CMM Jet Hit into string (for further processing)
      std::string CMMHit = Help->Binary((*it_CMMJetHits)->Hits(),24);
      
      mLog<<MSG::DEBUG<<"CMMJetHits dataID: "<< (*it_CMMJetHits)-> dataID() <<"   Hits: "<< (*it_CMMJetHits)-> Hits()
	  << " Hits(binary): " << CMMHit <<endreq;
      
      // ------------------------------------------------------------------------------------------
      // ----------------- Histos with distribution of JEM Hit Multiplicities ---------------------
      // ------------------------------------------------------------------------------------------
      //input data from JEMs have dataID 0..15
      if ((*it_CMMJetHits)-> dataID()<16) 
	{
	  //Fwd Hits left have dataID 0 or 8
	  if (((*it_CMMJetHits)-> dataID()==0)or((*it_CMMJetHits)-> dataID()==8) ) 
	    {
	      Help->FillHitsHisto(m_h_CMMJetHits_JEM_FwdHitsLeft, CMMHit, 0, 4, 8, 2, &mLog);
	      Help->FillHitsHisto(m_h_CMMJetHits_JEM_MainHits, CMMHit, 0, 8, 0, 2, &mLog);
	    }
	  else
	    {
	      //Fwd Hits right have dataID 7 or 15
	      if (((*it_CMMJetHits)-> dataID()==7)or((*it_CMMJetHits)-> dataID()==15) ) 
		{
		  Help->FillHitsHisto(m_h_CMMJetHits_JEM_FwdHitsRight, CMMHit, 0, 4, 8, 2, &mLog);
		  Help->FillHitsHisto(m_h_CMMJetHits_JEM_MainHits, CMMHit, 0, 8, 0, 2, &mLog);
		}
	      //Main Hits for all other modules
	      else 
		{
		  Help->FillHitsHisto(m_h_CMMJetHits_JEM_MainHits, CMMHit, 0, 8, 0, 3, &mLog);
		}
	    }
	}
      // ------------------------------------------------------------------------------------------
      // ----------------- Histos with S-Link errors -----------------------------
      // ------------------------------------------------------------------------------------------
      //only for Bytestream data
      if (m_DataType=="BS")
	{
	  LVL1::DataError err((*it_CMMJetHits)-> Error());
	  //input data from JEMs have dataID 0..15   ---  fill only parity errors
	  
	  int crate = (*it_CMMJetHits)->crate();
	  int module = (*it_CMMJetHits)-> dataID();
	  
	  if (module<16)
	    {
	      // Parity
	      m_h_CMMJet_error->Fill(1,(crate*19 + 1 + module),err.get(1));
	    }
	  
	  // errors from crate CMM
	  // fill parity and L1CaloSubStatus 
	  if (module==16)
	    {
	      // Parity; set only for crate CMM -> system CMM 
	      if (crate==1) m_h_CMMJet_error->Fill(1,(1 + 16),err.get(1));
	      
	      // set L1CaloSubStatus for both Crate and System CMM
	      // GLinkParity
	      m_h_CMMJet_error->Fill(5,(crate*19 + 1 + 16),err.get(16));
	      // GLinkProtocol
	      m_h_CMMJet_error->Fill(6,(crate*19 + 1 + 16),err.get(17));
	      // BCNMismatch
	      m_h_CMMJet_error->Fill(7,(crate*19 + 1 + 16),err.get(18));
	      // FIFOOverflow
	      m_h_CMMJet_error->Fill(8,(crate*19 + 1 + 16),err.get(19));
	      // ModuleError
	      m_h_CMMJet_error->Fill(9,(crate*19 + 1 + 16),err.get(20));
	      
	      // GLinkDown
	      m_h_CMMJet_error->Fill(10,(crate*19 + 1 + 16),err.get(22));
	      // GLinkTimeout
	      m_h_CMMJet_error->Fill(11,(crate*19 + 1 + 16),err.get(23));
	      // FailingBCN
	      if (err.get(24)!=0) m_h_CMMJet_error->Fill(12,(crate*19 + 1 + 16),1);
	    }
	  // fill parity for remote fwd hits
	  if ((module==19)and (crate==1))
	    {
	      // Parity
	      m_h_CMMJet_error->Fill(1,( 1 + 16),err.get(1));
	    }
	}
      
      // ------------------------------------------------------------------------------------------
      // ----------------- Histos with distribution of CMM Hit Multiplicities ---------------------
      // ------------------------------------------------------------------------------------------
      //main total jets have dataID 18
      if ((*it_CMMJetHits)-> dataID() == 18)  
	{
	  Help->FillHitsHisto(m_h_CMMJetHits_MainJets, CMMHit, 0, 8, 0, 3, &mLog);
	}
      
      //fwd total jets 
      if ((*it_CMMJetHits)-> dataID() == 21)  
	{
	  CMMHit= Help->Binary((*it_CMMJetHits)->Hits(),16); //total fwd jets only 16 bit long!
	  mLog<<MSG::DEBUG<<"Right|Left Total Jets  Hits: " << CMMHit <<endreq;
	  
	  Help->FillHitsHisto(m_h_CMMJetHits_FwdJetsLeft, CMMHit, 0, 4, 0, 2, &mLog);
	  Help->FillHitsHisto(m_h_CMMJetHits_FwdJetsRight, CMMHit, 0, 4, 4, 2, &mLog);
	}
      
      //JetEtSum Hitmap
      if ((*it_CMMJetHits)-> dataID() == 22)  
	{
	  CMMHit= Help->Binary((*it_CMMJetHits)->Hits(),4);
	  mLog<<MSG::DEBUG<<"JetEt Hits: " << CMMHit <<endreq;
	  
	  Help->FillHitsHisto(m_h_CMMJetHits_EtMap, CMMHit, 0, 4, 0, 1, &mLog);
	}
    }  
  
  // ------------------------------------------------------------------------------------------
  // ----------------- Backplane transmission check: JEMs -> CMM ------------------------------
  // ------------------------------------------------------------------------------------------
  // only for BS data
  if (m_DataType=="BS")
    {
      mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM Jet Hits transmission ---------------"<<endreq;
      
      //retrieve JEMHits from storegate for comparison with transmitted data stored in CMMJetHits
      const JEMHitsCollection* JEMHits;
      JEMHitsCollection::const_iterator it_JEMHits ;
      sc = m_storeGate->retrieve(JEMHits, m_JEMHitsLocation); 
      if( (sc==StatusCode::FAILURE) ) 
	{
	  mLog << MSG:: INFO<< "No JEMHits found in TES at "<< m_JEMHitsLocation << endreq ;
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
      int noMatchfound;

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
      
      //---------------------------------- backplane transmission JEMs -> CMM -----------------------------
      it_vCMMJEMHits=vCMMJEMHits.begin();
      // step through both vectors and compare...
      while (it_vCMMJEMHits<vCMMJEMHits.end())
	{
	  mLog<<MSG::VERBOSE<<"CMM Crate "<<(*it_vCMMJEMHits).crate()<<" Module "<<(*it_vCMMJEMHits).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Hit "<<(*it_vCMMJEMHits).Hits()<<endreq;
	  
	  found = 0;
	  noMatchfound=0;
	  it_vJEMHits=vJEMHits.begin();
	  
	  while ((found==0)and(it_vJEMHits<vJEMHits.end()))
	    {
	      mLog<<MSG::VERBOSE<<"JEM Crate "<<(*it_vJEMHits).crate()<<" Module "<<(*it_vJEMHits).module()<<endreq;
	      mLog<<MSG::VERBOSE<<"JEM Hit "<<(*it_vJEMHits).JetHits()<<endreq;
	      
	      if (((*it_vCMMJEMHits).crate()==(*it_vJEMHits).crate())
		  and((*it_vCMMJEMHits).dataID()==(*it_vJEMHits).module()))
		{
		  if ((*it_vCMMJEMHits).Hits()!=(*it_vJEMHits).JetHits())
		    {
		      mLog<<MSG::INFO<<"Transmission errors between JEM and CMM Hits: change of value"<<endreq;
		      mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMHits).crate()<<" Module "<<(*it_vJEMHits).module()<<endreq;
		      mLog<<MSG::DEBUG<<"JEM Hit "<<Help->Binary((*it_vJEMHits).JetHits(),24)<<endreq;
		      mLog<<MSG::DEBUG<<"CMM Hit "<<Help->Binary((*it_vCMMJEMHits).Hits(),24)<<endreq;
		      
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
	  mLog<<MSG::INFO<<vCMMJEMHits.size()<<" Transmission errors between JEM and CMM Hits: additional CMM information"<<endreq;
	  for( it_vCMMJEMHits  = vCMMJEMHits.begin(); it_vCMMJEMHits <  vCMMJEMHits. end(); ++it_vCMMJEMHits )
	    {
	      mLog<<MSG::DEBUG<<"CMM Crate "<<(*it_vCMMJEMHits).crate()<<" Module "<<(*it_vCMMJEMHits).dataID()<<endreq;
	      mLog<<MSG::DEBUG<<"CMM Hit "<<(*it_vCMMJEMHits).Hits()<<endreq;
	      
	      m_h_TransCheck_JEP->Fill(1,(*it_vCMMJEMHits).crate()*20+(*it_vCMMJEMHits).dataID()+1,1);
	    }
	}
      
      // if the JEMHits vector isn't empty, fill the error counter!
      if (vJEMHits.size()!=0)
	{
	  mLog<<MSG::INFO<<vJEMHits.size()<<" Transmission errors between JEM and CMM Hits: addtional JEM information"<<endreq;
	  for( it_vJEMHits  =  vJEMHits.begin(); it_vJEMHits <  vJEMHits. end(); ++it_vJEMHits )
	    {
	      mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMHits).crate()<<" Module "<<(*it_vJEMHits).module()<<endreq;
	      mLog<<MSG::DEBUG<<"JEM Hit "<<(*it_vJEMHits).JetHits()<<endreq;
	      
	      m_h_TransCheck_JEP->Fill(1,(*it_vJEMHits).crate()*20+(*it_vJEMHits).module()+1,1);
	    }
	}
      
      //---------------------------------- backplane transmission crate CMM -> system CMM -----------------------------
      noMatchfound=0;
      
      if ((CrateCMMMainHits!=SystemRemoteCMMMainHits) or (CrateCMMFwdHits!=SystemRemoteCMMFwdHits))
	{
	  noMatchfound=1;
	  mLog<<MSG::INFO<<"Transmission error between crate and system CMM: Hits"<<endreq;
	}
      m_h_TransCheck_JEP->Fill(1,18,noMatchfound);
    }
  
  // =============================================================================================
  // ================= Container: CMM Et Sums ====================================================
  // =============================================================================================
  LVL1::QuadLinear expand;
  int ex;
  int ey;
  int et;
  
  // retrieve CMM Et Sums from Storegate
  const CMMEtSumsCollection* CMMEtSums;
  sc = m_storeGate->retrieve(CMMEtSums, m_CMMEtSumsLocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMMEtSums found in TES at " << m_CMMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM Et Sums ---------------"<<endreq;
  
  // Step over all cells 
  CMMEtSumsCollection::const_iterator it_CMMEtSums ;
  for( it_CMMEtSums  = CMMEtSums ->begin(); it_CMMEtSums < CMMEtSums -> end(); ++it_CMMEtSums )
    {	  
     
      // ------------------------------------------------------------------------------------------
      // ----------------- Histos with distribution of JEM Energies -------------------------------
      // ------------------------------------------------------------------------------------------
      // JEM energy sums, dataID < 16
      if  ((*it_CMMEtSums)-> dataID()<16) 
	{
	  // note: JEM energies are compressed -> use QuadLinear to expand!
	  ex=expand.Expand( (*it_CMMEtSums)-> Ex());
	  ey=expand.Expand( (*it_CMMEtSums)-> Ey());
	  et=expand.Expand( (*it_CMMEtSums)-> Et() );
	  
	  m_h_CMMEtSums_JEM_Ex -> Fill( ex, 1.);
	  m_h_CMMEtSums_JEM_Ey -> Fill( ey, 1.);
	  m_h_CMMEtSums_JEM_Et -> Fill( et, 1.);
	}
      
      // -----------------------------------------------------------------------------------------
      // ----------------- Histos with distribution of total Energie per system -------------------
      // ------------------------------------------------------------------------------------------
      // total energy sums
      if  ((*it_CMMEtSums)-> dataID()==18) 
	{
	  // the 0th bit of Ex and Ey is their sign
	  // -> remove it, only the value is important
	  int Ex,Ey;
	  std::string temp;
	  temp= Help->Binary((*it_CMMEtSums)-> Ex(),15);
	  temp.assign(temp,temp.length()-14,14);
	  Ex=Help->Multiplicity(temp,0,14);
	  
	  temp= Help->Binary((*it_CMMEtSums)-> Ey(),15);
	  temp.assign(temp,temp.length()-14,14);
	  Ey=Help->Multiplicity(temp,0,14);

	  m_h_CMMEtSums_Ex -> Fill( Ex, 1.);
	  m_h_CMMEtSums_Ey -> Fill( Ey, 1.);
	  m_h_CMMEtSums_Et -> Fill( (*it_CMMEtSums)-> Et(), 1.);
	}
      
      //MissingEt Hitmap
      if ((*it_CMMEtSums)-> dataID() == 19)  
	{
	  std::string CMMHit= Help->Binary((*it_CMMEtSums)->Et(),8);
	  mLog<<MSG::DEBUG<<"MissingEt Hits: " << CMMHit <<endreq;
	  
	  Help->FillHitsHisto(m_h_CMMEtSums_MissingEtMap, CMMHit, 0, 8, 0, 1, &mLog);
	}
      
      //SumEt Hitmap
      if ((*it_CMMEtSums)-> dataID() == 20)  
	{
	  std::string CMMHit= Help->Binary((*it_CMMEtSums)->Et(),4);
	  mLog<<MSG::DEBUG<<"SumEt Hits: " << CMMHit <<endreq;
	  
	  Help->FillHitsHisto(m_h_CMMEtSums_SumEtMap, CMMHit, 0, 4, 0, 1, &mLog);
	}
      
      //only for Bytestream data
      if (m_DataType=="BS")
	{
	  // ------------------------------------------------------------------------------------------
	  // ----------------- Histos with S-Link errors -----------------------------
	  // ------------------------------------------------------------------------------------------
	  LVL1::DataError exerr((*it_CMMEtSums)-> ExError());
	  LVL1::DataError eyerr((*it_CMMEtSums)-> EyError());
	  LVL1::DataError eterr((*it_CMMEtSums)-> EtError());
	  int error=0;
	  
	  int crate = (*it_CMMEtSums)->crate();
	  int module = (*it_CMMEtSums)-> dataID();

	  //input data from JEMs have dataID 0..15   ---  fill only parity errors
	  if (module<16)
	    {
	      // Parity
	      error=0;
	      if ((exerr.get(1)==1)or(eyerr.get(1)==1)or(eterr.get(1)==1)) error=1;
	      m_h_CMMEnergy_error->Fill(1,(crate*19 + 1 + module),error);
	    }
	  
	  // errors from crate CMM
	  // fill parity and L1CaloSubStatus
	  if (module==16)
	    {
	      // Parity
	      error=0;
	      if ((exerr.get(1)==1)or(eyerr.get(1)==1)or(eterr.get(1)==1)) error=1;
	      m_h_CMMEnergy_error->Fill(1,(crate*19 + 1 + 16),error);
	      
	      // GLinkParity
	      error=0;
	      if ((exerr.get(16)==1)or(eyerr.get(16)==1)or(eterr.get(16)==1)) error=1;
	      m_h_CMMEnergy_error->Fill(5,(crate*19 + 1 + 16),error);
	      // GLinkProtocol
	      error=0;
	      if ((exerr.get(17)==1)or(eyerr.get(17)==1)or(eterr.get(17)==1)) error=1;
	      m_h_CMMEnergy_error->Fill(6,(crate*19 + 1 + 16),error);
	      // BCNMismatch
	      error=0;
	      if ((exerr.get(18)==1)or(eyerr.get(18)==1)or(eterr.get(18)==1)) error=1;
	      m_h_CMMEnergy_error->Fill(7,(crate*19 + 1 + 16),error);
	      // FIFOOverflow
	      error=0;
	      if ((exerr.get(19)==1)or(eyerr.get(19)==1)or(eterr.get(19)==1)) error=1;
	      m_h_CMMEnergy_error->Fill(8,(crate*19 + 1 + 16),error);
	      // ModuleError
	      error=0;
	      if ((exerr.get(20)==1)or(eyerr.get(20)==1)or(eterr.get(20)==1)) error=1;
	      m_h_CMMEnergy_error->Fill(9,(crate*19 + 1 + 16),error);
	      
	      // GLinkDown
	      error=0;
	      if ((exerr.get(22)==1)or(eyerr.get(22)==1)or(eterr.get(22)==1)) error=1;
	      m_h_CMMEnergy_error->Fill(10,(crate*19 + 1 + 16),error);
	      // GLinkTimeout
	      error=0;
	      if ((exerr.get(23)==1)or(eyerr.get(23)==1)or(eterr.get(23)==1)) error=1;
	      m_h_CMMEnergy_error->Fill(11,(crate*19 + 1 + 16),error);
	      // FailingBCN
	      error=0;
	      if ((exerr.get(24)!=0)or(eyerr.get(24)!=0)or(eterr.get(24)!=0)) error=1;
	      m_h_CMMEnergy_error->Fill(12,(crate*19 + 1 + 16),error);
	    }
	}
    }
  // ------------------------------------------------------------------------------------------
  // ----------------- Backplane transmission check: JEMs -> CMM ------------------------------
  // ------------------------------------------------------------------------------------------
  //only for Bytestream data
  if (m_DataType=="BS")
    {
      mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM EtSums transmission ---------------"<<endreq;
      const JEMEtSumsCollection* JEMEtSums;
      JEMEtSumsCollection::const_iterator it_JEMEtSums ;
      sc = m_storeGate->retrieve(JEMEtSums, m_JEMEtSumsLocation);
      if( (sc==StatusCode::FAILURE) ) 
	{
	  mLog << MSG::INFO << "No JEMEtSums found in TES at "<< m_JEMEtSumsLocation << endreq ;
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
      int noMatchfound=0;
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
      while (it_vCMMJEMEtSums<vCMMJEMEtSums.end())
	{
	  mLog<<MSG::VERBOSE<<"CMM Crate "<<(*it_vCMMJEMEtSums).crate()<<" Module "<<(*it_vCMMJEMEtSums).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Ex "<<(*it_vCMMJEMEtSums).Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Ey "<<(*it_vCMMJEMEtSums).Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Et "<<(*it_vCMMJEMEtSums).Et()<<endreq;
	  
	  bool found = 0;
	  noMatchfound = 0;
	  it_vJEMEtSums=vJEMEtSums.begin();
	  
	  while ((found==0)and(it_vJEMEtSums<vJEMEtSums.end()))
	    {
	      mLog<<MSG::VERBOSE<<"JEM Crate "<<(*it_vJEMEtSums).crate()<<" Module "<<(*it_vJEMEtSums).module()<<endreq;
	      mLog<<MSG::VERBOSE<<"JEM Ex "<<(*it_vJEMEtSums).Ex()<<endreq;
	      mLog<<MSG::VERBOSE<<"JEM Ey "<<(*it_vJEMEtSums).Ey()<<endreq;
	      mLog<<MSG::VERBOSE<<"JEM Et "<<(*it_vJEMEtSums).Et()<<endreq;
	      
	      if (((*it_vCMMJEMEtSums).crate()==(*it_vJEMEtSums).crate())
		  and((*it_vCMMJEMEtSums).dataID()==(*it_vJEMEtSums).module()))
		{
		  if (((*it_vCMMJEMEtSums).Ex()!=(*it_vJEMEtSums).Ex())
		      or ((*it_vCMMJEMEtSums).Ey()!=(*it_vJEMEtSums).Ey())
		      or ((*it_vCMMJEMEtSums).Et()!=(*it_vJEMEtSums).Et()))
		    {
		      mLog<<MSG::INFO<<"Transmission errors between JEM and CMM EtSums: change of value"<<endreq;
		      mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMEtSums).crate()<<" Module "<<(*it_vJEMEtSums).module()<<endreq;
		      mLog<<MSG::DEBUG<<"JEM Ex "<<(*it_vJEMEtSums).Ex()<<endreq;
		      mLog<<MSG::DEBUG<<"JEM Ey "<<(*it_vJEMEtSums).Ey()<<endreq;
		      mLog<<MSG::DEBUG<<"JEM Et "<<(*it_vJEMEtSums).Et()<<endreq;
		      
		      mLog<<MSG::DEBUG<<"CMM Ex "<<(*it_vCMMJEMEtSums).Ex()<<endreq;
		      mLog<<MSG::DEBUG<<"CMM Ey "<<(*it_vCMMJEMEtSums).Ey()<<endreq;
		      mLog<<MSG::DEBUG<<"CMM Et "<<(*it_vCMMJEMEtSums).Et()<<endreq;
		      
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
	  mLog<<MSG::INFO<<vCMMJEMEtSums.size()<<" Transmission errors between JEM and CMM EtSums: additional CMM information"<<endreq;
	  for( it_vCMMJEMEtSums  = vCMMJEMEtSums.begin(); it_vCMMJEMEtSums <  vCMMJEMEtSums. end(); ++it_vCMMJEMEtSums )
	    {
	      mLog<<MSG::DEBUG<<"CMM Crate "<<(*it_vCMMJEMEtSums).crate()<<" Module "<<(*it_vCMMJEMEtSums).dataID()<<endreq;
	      mLog<<MSG::DEBUG<<"CMM Ex "<<(*it_vCMMJEMEtSums).Ex()<<endreq;
	      mLog<<MSG::DEBUG<<"CMM Ey "<<(*it_vCMMJEMEtSums).Ey()<<endreq;
	      mLog<<MSG::DEBUG<<"CMM Et "<<(*it_vCMMJEMEtSums).Et()<<endreq;
	      
	      m_h_TransCheck_JEP->Fill(2,(*it_vCMMJEMEtSums).crate()*20+(*it_vCMMJEMEtSums).dataID() + 1,noMatchfound);
	    }
	}
      
      // if the JEMEtSums vector isn't empty, fill the error counter!
      if (vJEMEtSums.size()!=0)
	{
	  mLog<<MSG::INFO<<vJEMEtSums.size()<<" Transmission errors between JEM and CMM EtSums: addtional JEM information"<<endreq;
	  for( it_vJEMEtSums  =  vJEMEtSums.begin(); it_vJEMEtSums <  vJEMEtSums. end(); ++it_vJEMEtSums )
	    {
	      mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMEtSums).crate()<<" Module "<<(*it_vJEMEtSums).module()<<endreq;
	      mLog<<MSG::DEBUG<<"JEM Ex "<<(*it_vJEMEtSums).Ex()<<endreq;
	      mLog<<MSG::DEBUG<<"JEM Ey "<<(*it_vJEMEtSums).Ey()<<endreq;
	      mLog<<MSG::DEBUG<<"JEM Et "<<(*it_vJEMEtSums).Et()<<endreq;
	      
	      m_h_TransCheck_JEP->Fill(2,(*it_vJEMEtSums).crate()*20+(*it_vJEMEtSums).module() + 1,noMatchfound);
	    }
	}
      //---------------------------------- backplane transmission crate CMM -> system CMM -----------------------------
      noMatchfound=0;
      if ((CrateEx!=SystemEx)or (CrateEy!=SystemEy)or (CrateEt!=SystemEt))
	{
	  noMatchfound=1;
	  mLog<<MSG::INFO<<"Transmission error between crate and system CMM: Energy"<<endreq;
	}
      m_h_TransCheck_JEP->Fill(2,17,noMatchfound);
    }

  // =============================================================================================
  // ================= Container: CMM RoI ========================================================
  // =============================================================================================
  
  // retrieve RoI information from Storegate
  LVL1::CMMRoI* CR = new LVL1::CMMRoI ;
  sc = m_storeGate->retrieve (CR, m_CMMRoILocation);
  if (sc==StatusCode::FAILURE)
    {
      mLog <<MSG::INFO<<"No CMM RoI found in TES at "<< m_CMMRoILocation<<endreq;
      return StatusCode::SUCCESS;    
    }

  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" CMM RoI ---------------"<<endreq;

  // ------------------------------------------------------------------------------------------
  // ----------------- Histos filled with CMM RoI information ---------------------------------
  // ------------------------------------------------------------------------------------------

  mLog<<MSG::DEBUG<<"JetEtHits: "<<Help->Binary((CR)->jetEtHits(),4)<<"; SumEtHits: "<<Help->Binary((CR)->sumEtHits(),4)<<"; MissingEtHits: "<<Help->Binary((CR)->missingEtHits(),8)<<endreq;

  // Jet Et Hits
  std::string CMMRoIHit = Help->Binary((CR)-> jetEtHits(),4);
  Help->FillHitsHisto(m_h_CMMRoI_JetEtHits, CMMRoIHit, 0, 4, 0, 1, &mLog);

  // Sum Et Hits
  CMMRoIHit = Help->Binary((CR)-> sumEtHits(),4);
  Help->FillHitsHisto(m_h_CMMRoI_SumEtHits, CMMRoIHit, 0, 4, 0, 1, &mLog);

  // Missing Et Hits
  CMMRoIHit = Help->Binary((CR)-> missingEtHits(),8);
  Help->FillHitsHisto(m_h_CMMRoI_MissingEtHits , CMMRoIHit, 0, 8, 0, 1, &mLog);
 
  // total Ex, Ey, Et
  // the 0th bit of Ex and Ey is their sign
  // -> remove it, only the value is important
  int Ex,Ey;
  Ex=abs((CR)-> ex());
  Ey=abs((CR)-> ey());
  
  mLog<<MSG::DEBUG<<" Ex: "<<Ex<<"; Ey: "<<Ey<<"; Et "<<(CR)->et()<<endreq;
  
  m_h_CMMRoI_Ex -> Fill( Ex,1);
  m_h_CMMRoI_Ey -> Fill( Ey,1);
  m_h_CMMRoI_Et -> Fill( (CR)->et(),1);
  
  // errors
  if (m_DataType=="BS")
    {
      LVL1::DataError exerr((CR)-> exError());
      LVL1::DataError eyerr((CR)-> eyError());
      LVL1::DataError eterr((CR)-> etError());
      LVL1::DataError jetEterr((CR)-> jetEtError());

      // Parity (Ex)
      m_h_CMMRoI_error->Fill(1,exerr.get(1));
      // Parity (Ey,SumEtMap)
      m_h_CMMRoI_error->Fill(2,eyerr.get(1));
      // Parity (Et,MissingEtMap)
      m_h_CMMRoI_error->Fill(3,eterr.get(1));
      // Parity (JetEtMap)
      m_h_CMMRoI_error->Fill(4,jetEterr.get(1));
      
      // Overflow (Ex)
      m_h_CMMRoI_error->Fill(6,exerr.get(0));
      // Overflow (Ey)
       m_h_CMMRoI_error->Fill(7,eyerr.get(0));
     // Overflow (Et)
      m_h_CMMRoI_error->Fill(8,eterr.get(0)); 
    }

  return StatusCode( StatusCode::SUCCESS );
}

/*---------------------------------------------------------*/
StatusCode CMMMon::procHistograms( bool isEndOfEventsBlock, 
				  bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  if( isEndOfEventsBlock || isEndOfLumiBlock ) 
    {  
    }
	
  if( isEndOfRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}

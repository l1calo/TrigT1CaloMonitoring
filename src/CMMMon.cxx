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
#include "TrigT1CaloMonitoring/MonHelper.h"

//#include "TrigT1Calo/EnergyTrigger.h"
#include "TrigT1Calo/LVL1TriggerMenuDefs.h"
//#include "TrigT1Calo/LVL1TriggerMenu.h"
//#include "TrigT1Calo/InternalJetROI.h"
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
  declareProperty( "EventNoInHistoTitle", m_EventNoInHisto = 1) ;

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

  MonGroup CMM_DAQ ( this, (m_PathInRootFile+"_DAQ").c_str(), expert, eventsBlock );
  HistoBooker* DAQ_Booker = new HistoBooker(&CMM_DAQ, &mLog, m_DataType);

  MonGroup CMM_input ( this, (m_PathInRootFile + "_input").c_str(), expert, eventsBlock );
  HistoBooker* input_Booker = new HistoBooker(&CMM_input, &mLog, m_DataType);
  
  MonGroup CMM_RoI ( this, (m_PathInRootFile + "_RoI").c_str(), LevelOfDetail, eventsBlock );
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
      m_NoEvents=0;
      //----------------------------------  CMM Input data from JEMs -----------------------------
      m_h_CMMJetHits_JEM_MainHits=input_Booker->book1F("MainHits_CMM_input", "Main Jet Hit Multiplicity per Threshold  --  CMM input", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMJetHits_JEM_FwdHitsRight=input_Booker->book1F("FwdHitsRight_CMM_input", "Forward Right Jet Hit Multiplicity per Threshold  --  CMM input",4 , -0.5, 3.5, "Threshold No.", "#");
      m_h_CMMJetHits_JEM_FwdHitsLeft=input_Booker->book1F("FwdHitsLeft_CMM_input", "Forward Left Jet Hit Multiplicity per Threshold  --  CMM input", 4 , -0.5, 3.5,  "Threshold No.", "#");

      m_h_CMMEtSums_JEM_Ex=input_Booker->book1F("Ex_CMM_input", "CMM E_{x}^{JEM}  --  CMM input", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_JEM_Ey=input_Booker->book1F("Ey_CMM_input", "CMM E_{y}^{JEM}  --  CMM input", 250, 0,250, "Ey [GeV]", "#");
      m_h_CMMEtSums_JEM_Et=input_Booker->book1F("Et_CMM_input", "CMM E_{t}^{JEM}  --  CMM input", 250, 0,250, "Et [GeV]", "#");
      

      //---------------------------------- CMM output to DAQ -----------------------------
      m_h_CMMJetHits_MainJets = DAQ_Booker->book1F("TotalMainHits_CMM_DAQ", "Main Jet Hit Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMJetHits_FwdJetsRight = DAQ_Booker->book1F("TotalFwdHitsRight_CMM_DAQ", "Forward Right Jet Hit Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5, "Threshold No.", "#");
      m_h_CMMJetHits_FwdJetsLeft = DAQ_Booker->book1F("TotalFwdHitsLeft_CMM_DAQ", "Forward Left Jet Hit Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5,  "Threshold No.", "#");
      m_h_CMMJetHits_EtMap = DAQ_Booker->book1F("JetEtHits_CMM_DAQ", "JetEt Hit Multiplicity per Threshold  --  CMM DAQ", 4 ,-0.5, 3.5, "Threshold No.", "#");
      m_h_CMMEtSums_MissingEtMap = DAQ_Booker->book1F("MissingEtHits_CMM_DAQ", "MissingEt Hit Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMEtSums_SumEtMap = DAQ_Booker->book1F("SumEtHits_CMM_DAQ", "SumEt Hit Multiplicity per Threshold  --  CMM DAQ", 4, -0.5, 3.5, "Threshold No.", "#");

      m_h_CMMEtSums_Ex = DAQ_Booker->book1F("Ex_CMM_DAQ", "CMM E_{x}^{CMM}  --  CMM DAQ", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_Ey = DAQ_Booker->book1F("Ey_CMM_DAQ", "CMM E_{y}^{CMM}  --  CMM DAQ", 250, 0,250, "Ey [GeV]", "#");
      m_h_CMMEtSums_Et = DAQ_Booker->book1F("Et_CMM_DAQ", "CMM E_{t}^{CMM}  --  CMM DAQ", 250, 0,250, "Et [GeV]", "#");


      //---------------------------------- CMM output to RoI -----------------------------
      m_h_CMMRoI_JetEtHits =RoI_Booker->book1F("JetEtHits_CMM_RoI","JetEt Hit Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","#");
      m_h_CMMRoI_MissingEtHits =RoI_Booker->book1F("MissingEtHits_CMM_RoI","MissingEt Hit Multiplicity per Threshold  --  CMM RoI", 8, -0.5,7.5,"Threshold No.","#");
      m_h_CMMRoI_SumEtHits =RoI_Booker->book1F("SumEtHits_CMM_RoI","SumEt Hit Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","#");

      m_h_CMMRoI_Ex = RoI_Booker->book1F("Ex_CMM_RoI", "CMM E_{x}^{CMM}  --  CMM RoI", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMRoI_Ey = RoI_Booker->book1F("Ey_CMM_RoI", "CMM E_{y}^{CMM}  --  CMM RoI", 250, 0,250, "Ey [GeV]", "#");
      m_h_CMMRoI_Et = RoI_Booker->book1F("Et_CMM_RoI", "CMM E_{t}^{CMM}  --  CMM RoI", 250, 0,250, "Et [GeV]", "#");


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
  m_NoEvents++;

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
	  mLog<<MSG::DEBUG<<"    Ex: "<<Ex<<"; Ey: "<<Ey<<"; Et "<<(*it_CMMEtSums)->Et()<<endreq;
	  mLog<<MSG::DEBUG<<"Roh Ex: "<<(*it_CMMEtSums)->Ex()<<"; Ey: "<<(*it_CMMEtSums)->Ey()<<"; Et "<<(*it_CMMEtSums)->Et()<<endreq;

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
 

  // the 0th bit of Ex and Ey is their sign
  // -> remove it, only the value is important
  int Ex,Ey;
  std::string temp;
  temp= Help->Binary((CR)-> ex(),15);
  temp.assign(temp,temp.length()-14,14);
  Ex=Help->Multiplicity(temp,0,14);
  
  temp= Help->Binary((CR)-> ey(),15);
  temp.assign(temp,temp.length()-14,14);
  Ey=Help->Multiplicity(temp,0,14);
  
  mLog<<MSG::DEBUG<<"    Ex: "<<Ex<<"; Ey: "<<Ey<<"; Et "<<(CR)->et()<<endreq;
  mLog<<MSG::DEBUG<<"Roh Ex: "<<(CR)-> ex()<<"; Ey: "<<(CR)-> ey()<<"; Et "<<(CR)-> et()<<endreq;

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
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock ) 
    {  
    }
	
  if(m_EventNoInHisto==1)
    {
      if (m_DataType=="BS")
	{
	  if( isEndOfRun ) { 
	    std::stringstream buffer;
	    buffer.str("");
	    buffer<<m_NoEvents;
	    std::string title;
	    
	    title = m_h_CMMJet_error-> GetTitle();
	    title=title + " | #events: " + buffer.str();
	    m_h_CMMJet_error->SetTitle(title.c_str());
	    
	    title = m_h_CMMEnergy_error-> GetTitle();
	    title=title + " | #events: " + buffer.str();
	    m_h_CMMEnergy_error->SetTitle(title.c_str());
	    
	    title = m_h_CMMRoI_error-> GetTitle();
	    title=title + " | #events: " + buffer.str();
	    m_h_CMMRoI_error->SetTitle(title.c_str());
	  }
	}
    }
  return StatusCode( StatusCode::SUCCESS );
}

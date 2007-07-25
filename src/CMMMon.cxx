
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


#include "TrigT1CaloMonitoring/CMMMon.h"
#include "TrigT1CaloMonitoring/JEMMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

//#include "TrigT1Calo/EnergyTrigger.h"
#include "TrigT1Calo/LVL1TriggerMenuDefs.h"
#include "TrigT1Calo/LVL1TriggerMenu.h"
#include "TrigT1Calo/InternalJetROI.h"
#include "TrigT1Calo/CMMRoI.h"
#include "TrigT1Calo/QuadLinear.h"


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
  
  MonGroup CMM_transmission ( this, (m_PathInRootFile + "/TransmissionCheck").c_str(), shift, eventsBlock );
  HistoBooker* transmission_Booker = new HistoBooker(&CMM_transmission, &mLog, m_DataType);
  
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }
  
  if( isNewEventsBlock || isNewLumiBlock ) 
    {	

      // CMM Input data
      m_h_CMMJetHits_JEM_MainHits=input_Booker->book1F("MainHits_CMM_input", "Main Jet Hit Multiplicity per Threshold  --  CMM input", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMJetHits_JEM_FwdHitsRight=input_Booker->book1F("FwdHitsRight_CMM_input", "Forward Right Jet Hit Multiplicity per Threshold  --  CMM input",4 , -0.5, 3.5, "Threshold No.", "#");
      m_h_CMMJetHits_JEM_FwdHitsLeft=input_Booker->book1F("FwdHitsLeft_CMM_input", "Forward Left Jet Hit Multiplicity per Threshold  --  CMM input", 4 , -0.5, 3.5,  "Threshold No.", "#");

      m_h_CMMEtSums_JEM_Ex=input_Booker->book1F("Ex_CMM_input", "CMM E_{x}  --  CMM input", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_JEM_Ey=input_Booker->book1F("Ey_CMM_input", "CMM E_{y}  --  CMM input", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_JEM_Et=input_Booker->book1F("Et_CMM_input", "CMM E_{t}  --  CMM input", 250, 0,250, "Ex [GeV]", "#");
      

      // CMM DAQ data
      m_h_CMMJetHits_MainJets = DAQ_Booker->book1F("TotalMainHits_CMM_DAQ", "Main Jet Hit Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMJetHits_FwdJetsRight = DAQ_Booker->book1F("TotalFwdHitsRight_CMM_DAQ", "Forward Right Jet Hit Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5, "Threshold No.", "#");
      m_h_CMMJetHits_FwdJetsLeft = DAQ_Booker->book1F("TotalFwdHitsLeft_CMM_DAQ", "Forward Left Jet Hit Multiplicity per Threshold  --  CMM DAQ", 4 , -0.5, 3.5,  "Threshold No.", "#");
      m_h_CMMJetHits_EtMap = DAQ_Booker->book1F("JetEtHits_CMM_DAQ", "JetEt Hit Multiplicity per Threshold  --  CMM DAQ", 4 ,-0.5, 3.5, "Threshold No.", "#");
      m_h_CMMEtSums_MissingEtMap = DAQ_Booker->book1F("MissingEtHits_CMM_DAQ", "MissingEt Hit Multiplicity per Threshold  --  CMM DAQ", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMEtSums_SumEtMap = DAQ_Booker->book1F("SumEtHits_CMM_DAQ", "SumEt Hit Multiplicity per Threshold  --  CMM DAQ", 4, -0.5, 3.5, "Threshold No.", "#");

      m_h_CMMEtSums_Ex = DAQ_Booker->book1F("Ex_CMM_DAQ", "CMM E_{x}  --  CMM DAQ", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_Ey = DAQ_Booker->book1F("Ey_CMM_DAQ", "CMM E_{y}  --  CMM DAQ", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_Et = DAQ_Booker->book1F("Et_CMM_DAQ", "CMM E_{t}  --  CMM DAQ", 250, 0,250, "Ex [GeV]", "#");


      //CMM RoI data
      m_h_CMMRoI_JetEtHits =RoI_Booker->book1F("JetEtHits_CMM_RoI","JetEt Hit Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","#");
      m_h_CMMRoI_MissingEtHits =RoI_Booker->book1F("MissingEtHits_CMM_RoI","MissingEt Hit Multiplicity per Threshold  --  CMM RoI", 8, -0.5,7.5,"Threshold No.","#");
      m_h_CMMRoI_SumEtHits =RoI_Booker->book1F("SumEtHits_CMM_RoI","SumEt Hit Multiplicity per Threshold  --  CMM RoI", 4, -0.5,3.5,"Threshold No.","#");

      m_h_CMMRoI_Ex = RoI_Booker->book1F("Ex_CMM_RoI", "CMM E_{x}  --  CMM RoI", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMRoI_Ey = RoI_Booker->book1F("Ey_CMM_RoI", "CMM E_{y}  --  CMM RoI", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMRoI_Et = RoI_Booker->book1F("Et_CMM_RoI", "CMM E_{t}  --  CMM RoI", 250, 0,250, "Ex [GeV]", "#");


      if (m_DataType=="BS")
	{
	  //Parity Errors
	  m_h_CMMJetHits_JEM_Crate0ParityError = transmission_Booker->book1F("Crate0JEMHitsParityError_CMM_transmission", "Crate0 JEM Hits Parity Error  --  CMM Transmission",
									     16, -0.5, 15.5, "Module No.", "#");
	  m_h_CMMJetHits_JEM_Crate1ParityError = transmission_Booker->book1F("Crate1JEMHitsParityError_CMM_transmission", "Crate1 JEM Hits Parity Error  --  CMM Transmission",
									     16, -0.5, 15.5, "Module No.", "#");
	  
	  //transmission checks
	  m_h_TransCheck_JEP=transmission_Booker->book2F("JEP_TransCheck", "JEP Backplane Transmission Check JEM -> CMM per Module and Crate", 2,0.5,2.5,40,-0.5,39.5, "", "");
	  //m_h_TransCheck_JEP-> SetOption ("text");
	  m_h_TransCheck_JEP->GetXaxis()->SetBinLabel(1, "Energy");
	  m_h_TransCheck_JEP->GetXaxis()->SetBinLabel(2, "Hits");
      
	  std::string name;
	  std::stringstream buffer;
      
	  for (int i = 0; i < 16; i++)
	    {
	      buffer.str("");
	      buffer<<i;
	      
	      name = "JEM " + buffer.str();
	      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel((i+1), name.c_str());
	      
	      buffer.str("");
	      buffer<<i;
	      
	      name = "JEM " + buffer.str();
	      m_h_TransCheck_JEP->GetYaxis()->SetBinLabel((i+1+19), name.c_str());
	    }
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(17, "Crate Energy CMM");
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(18, "Crate Jet CMM");
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(19, "Crate 1: ");
	  
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(36, "Crate Energy CMM");
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(37, "System Energy CMM");
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(38, "Crate Jet CMM");
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(39, "System Jet CMM");
	  m_h_TransCheck_JEP->GetYaxis()->SetBinLabel(40, "Crate 0: ");
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
  
  
  //******************************* 
  //********* CMM Jet Hits ******** 
  //******************************* 

  //fill  Histos
  const CMMJetHitsCollection* CMMJetHits;

  StatusCode sc = m_storeGate->retrieve(CMMJetHits, m_CMMJetHitsLocation);
  
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No CMM JetHits found in TES at "  << m_CMMJetHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }
  
  mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM Jet Hits ---------------"<<endreq;
  
  // Step over all cells and put into hist
  CMMJetHitsCollection::const_iterator it_CMMJetHits ;
  
  for( it_CMMJetHits  = CMMJetHits ->begin(); it_CMMJetHits < CMMJetHits -> end(); ++it_CMMJetHits )
    {	  
      std::string CMMHit = Help->Binary((*it_CMMJetHits)->Hits(),24);

      mLog<<MSG::DEBUG<<"CMMJetHits dataID: "<< (*it_CMMJetHits)-> dataID() <<"   Hits: "<< (*it_CMMJetHits)-> Hits()
	  << " Hits(binary): " << CMMHit <<endreq;

      if ((*it_CMMJetHits)-> dataID()<16) //input data from JEMs
	{
	  if (((*it_CMMJetHits)-> dataID()==0)or((*it_CMMJetHits)-> dataID()==8) ) //Fwd Hits left
	    {
	      Help->FillHitsHisto(m_h_CMMJetHits_JEM_FwdHitsLeft, CMMHit, 0, 4, 8, 2, &mLog);
	      Help->FillHitsHisto(m_h_CMMJetHits_JEM_MainHits, CMMHit, 0, 8, 0, 2, &mLog);
	    }
	  else
	    {
	      if (((*it_CMMJetHits)-> dataID()==7)or((*it_CMMJetHits)-> dataID()==15) ) //Fwd Hits right
		{
		  Help->FillHitsHisto(m_h_CMMJetHits_JEM_FwdHitsRight, CMMHit, 0, 4, 8, 2, &mLog);
		   Help->FillHitsHisto(m_h_CMMJetHits_JEM_MainHits, CMMHit, 0, 8, 0, 2, &mLog);
		}
	      else //Main Hits
		{
		  Help->FillHitsHisto(m_h_CMMJetHits_JEM_MainHits, CMMHit, 0, 8, 0, 3, &mLog);
		}
	    }

	  if (m_DataType=="BS")
	    {
	      if ((*it_CMMJetHits)-> crate()==0) //Fill ErrorCounter
		{
		  m_h_CMMJetHits_JEM_Crate0ParityError -> Fill(((*it_CMMJetHits)-> dataID()),((*it_CMMJetHits)-> Error()));
		}
	      if ((*it_CMMJetHits)-> crate()==1) //Fill ErrorCounter
		{
		  m_h_CMMJetHits_JEM_Crate1ParityError -> Fill(((*it_CMMJetHits)-> dataID()),((*it_CMMJetHits)-> Error()));
		}
	    }
	}
      if ((*it_CMMJetHits)-> dataID() == 18)  //main total jets
	{
	  Help->FillHitsHisto(m_h_CMMJetHits_MainJets, CMMHit, 0, 8, 0, 3, &mLog);
	}

      if ((*it_CMMJetHits)-> dataID() == 21)  //fwd total jets
	{
	  CMMHit= Help->Binary((*it_CMMJetHits)->Hits(),16); //total fwd jets only 16 bit long!
	  mLog<<MSG::DEBUG<<"Right|Left Total Jets  Hits: " << CMMHit <<endreq;
	  
	  Help->FillHitsHisto(m_h_CMMJetHits_FwdJetsLeft, CMMHit, 0, 4, 0, 2, &mLog);
	  Help->FillHitsHisto(m_h_CMMJetHits_FwdJetsRight, CMMHit, 0, 4, 4, 2, &mLog);
	}

      if ((*it_CMMJetHits)-> dataID() == 22)  //JetEtSum Hitmap
	{
	  CMMHit= Help->Binary((*it_CMMJetHits)->Hits(),4);
	  mLog<<MSG::DEBUG<<"JetEt Hits: " << CMMHit <<endreq;

	  Help->FillHitsHisto(m_h_CMMJetHits_EtMap, CMMHit, 0, 4, 0, 1, &mLog);
	}
       
    }    

  if (m_DataType=="BS")
    {

      mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM Jet Hits transmission ---------------"<<endreq;
      //Transmission check
      const JEMHitsCollection* JEMHits;
      JEMHitsCollection::const_iterator it_JEMHits ;
      sc = m_storeGate->retrieve(JEMHits, m_JEMHitsLocation);
      
      if( (sc==StatusCode::FAILURE) ) 
	{
	  mLog << MSG:: INFO<< "No JEMHits found in TES at "<< m_JEMHitsLocation << endreq ;
	  return StatusCode::SUCCESS;
	}
      
      std::vector <LVL1::CMMJetHits>  vCMMJEMHits;
      std::vector <LVL1::JEMHits>  vJEMHits;
      //use standard vector instead of datavector for transmissioncheck:
      //datavector erases not only the pointer  in the vector, but also the referenced object
      //-> segmentation fault!
      
      for( it_CMMJetHits  = CMMJetHits ->begin(); it_CMMJetHits < CMMJetHits -> end(); ++it_CMMJetHits )
	{	  
	  if ((*it_CMMJetHits)->dataID()<16)
	    {
	      vCMMJEMHits.push_back(**it_CMMJetHits );
	    }
	}
      for( it_JEMHits  =  JEMHits->begin(); it_JEMHits <  JEMHits-> end(); ++it_JEMHits )
	{
	  vJEMHits.push_back(**it_JEMHits);
	}
      
      bool found;
      std::vector <LVL1::CMMJetHits>::iterator it_vCMMJEMHits;
      std::vector <LVL1::JEMHits>::iterator it_vJEMHits;
      
      it_vCMMJEMHits=vCMMJEMHits.begin();
      while (it_vCMMJEMHits<vCMMJEMHits.end())
	{
	  mLog<<MSG::VERBOSE<<"CMM Crate "<<(*it_vCMMJEMHits).crate()<<" Module "<<(*it_vCMMJEMHits).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Hit "<<(*it_vCMMJEMHits).Hits()<<endreq;
	  
	  found = 0;
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
		      
		      if ((*it_vJEMHits).crate()==0) 
			{
			  m_h_TransCheck_JEP->Fill(1,(*it_vJEMHits).module(),1);
			}
		      else 
			{
			  m_h_TransCheck_JEP->Fill(1,((*it_vJEMHits).module()+19),1);
			}
		    }
		  vCMMJEMHits.erase(it_vCMMJEMHits);
		  vJEMHits.erase(it_vJEMHits);
		  found=1;
		}
	      else it_vJEMHits=it_vJEMHits+1;
	    }
	  if (found==0)it_vCMMJEMHits=it_vCMMJEMHits+1;
	}
      
      if (vCMMJEMHits.size()!=0)
	{
	  mLog<<MSG::INFO<<vCMMJEMHits.size()<<" Transmission errors between JEM and CMM Hits: additional CMM information"<<endreq;
	  
	  //fill errorcounter
	  for( it_vCMMJEMHits  = vCMMJEMHits.begin(); it_vCMMJEMHits <  vCMMJEMHits. end(); ++it_vCMMJEMHits )
	    {
	      mLog<<MSG::DEBUG<<"CMM Crate "<<(*it_vCMMJEMHits).crate()<<" Module "<<(*it_vCMMJEMHits).dataID()<<endreq;
	      mLog<<MSG::DEBUG<<"CMM Hit "<<(*it_vCMMJEMHits).Hits()<<endreq;
	      
	      if ((*it_vCMMJEMHits).crate()==0) 
		{
		  m_h_TransCheck_JEP->Fill(1,(*it_vCMMJEMHits).dataID(),1);
		}
	      else 
		{
		  m_h_TransCheck_JEP->Fill(1,((*it_vCMMJEMHits).dataID()+19),1);
		}
	      
	    }
	}
      
      if (vJEMHits.size()!=0)
	{
	  mLog<<MSG::INFO<<vJEMHits.size()<<" Transmission errors between JEM and CMM Hits: addtional JEM information"<<endreq;
	  
	  //fill errorcounter
	  for( it_vJEMHits  =  vJEMHits.begin(); it_vJEMHits <  vJEMHits. end(); ++it_vJEMHits )
	    {
	      mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMHits).crate()<<" Module "<<(*it_vJEMHits).module()<<endreq;
	      mLog<<MSG::DEBUG<<"JEM Hit "<<(*it_vJEMHits).JetHits()<<endreq;
	      
	      if ((*it_vJEMHits).crate()==0) 
		{
		  m_h_TransCheck_JEP->Fill(1,(*it_vJEMHits).module(),1);
		}
	      else 
		{
		  m_h_TransCheck_JEP->Fill(1,((*it_vJEMHits).module()+19),1);
		}
	      
	    }
	}
    }
 	
//******************************* 
//********* CMM Et Sums ********* 
//******************************* 
  LVL1::QuadLinear expand;
  int ex;
  int ey;
  int et;

  //fill  Histos
  const CMMEtSumsCollection* CMMEtSums;
  sc = m_storeGate->retrieve(CMMEtSums, m_CMMEtSumsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO
	   << "No CMMEtSums found in TES at "
	   << m_CMMEtSumsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM Et Sums ---------------"<<endreq;

  // Step over all cells and put into hist
  CMMEtSumsCollection::const_iterator it_CMMEtSums ;

  for( it_CMMEtSums  = CMMEtSums ->begin(); it_CMMEtSums < CMMEtSums -> end(); ++it_CMMEtSums )
    {	  
     // fill histograms
      if  ((*it_CMMEtSums)-> dataID()<16) // JEM energy sums  important: energies are compressed
	{
	  ex=expand.Expand( (*it_CMMEtSums)-> Ex());
	  ey=expand.Expand( (*it_CMMEtSums)-> Ey());
	  et=expand.Expand( (*it_CMMEtSums)-> Et() );


	  m_h_CMMEtSums_JEM_Ex -> Fill( ex, 1.);
	  m_h_CMMEtSums_JEM_Ey -> Fill( ey, 1.);
	  m_h_CMMEtSums_JEM_Et -> Fill( et, 1.);
	}

      if  ((*it_CMMEtSums)-> dataID()==18) // total energy sums
	{
	  m_h_CMMEtSums_Ex -> Fill( (*it_CMMEtSums)-> Ex(), 1.);
	  m_h_CMMEtSums_Ey -> Fill( (*it_CMMEtSums)-> Ey(), 1.);
	  m_h_CMMEtSums_Et -> Fill( (*it_CMMEtSums)-> Et(), 1.);
	}
      
      if ((*it_CMMEtSums)-> dataID() == 19)  //MissingEt Hitmap
	{
	  std::string CMMHit= Help->Binary((*it_CMMEtSums)->Et(),8);
	  mLog<<MSG::DEBUG<<"MissingEt Hits: " << CMMHit <<endreq;

	  Help->FillHitsHisto(m_h_CMMEtSums_MissingEtMap, CMMHit, 0, 8, 0, 1, &mLog);
	}

      if ((*it_CMMEtSums)-> dataID() == 20)  //SumEt Hitmap
	{
	  std::string CMMHit= Help->Binary((*it_CMMEtSums)->Et(),4);
	  mLog<<MSG::DEBUG<<"SumEt Hits: " << CMMHit <<endreq;

	  Help->FillHitsHisto(m_h_CMMEtSums_SumEtMap, CMMHit, 0, 4, 0, 1, &mLog);
	}

      mLog<<MSG::VERBOSE<<"CMMEtSums crate: "<< (*it_CMMEtSums)-> crate()<<" dataID: "<< (*it_CMMEtSums)-> dataID()
	  <<" Ex: "<<  (*it_CMMEtSums)-> Ex()
	  <<" Ey: "<<  (*it_CMMEtSums)-> Ey()
	  <<" Et: "<<  (*it_CMMEtSums)-> Et()<<endreq;
    }     	

  if (m_DataType=="BS")
    {
      
      mLog<<MSG::DEBUG<<"--------------  "<< m_DataType<<" CMM EtSums transmission ---------------"<<endreq;
      //Transmission check
      const JEMEtSumsCollection* JEMEtSums;
      JEMEtSumsCollection::const_iterator it_JEMEtSums ;
      sc = m_storeGate->retrieve(JEMEtSums, m_JEMEtSumsLocation);
      if( (sc==StatusCode::FAILURE) ) 
	{
	  mLog << MSG::INFO << "No JEMEtSums found in TES at "<< m_JEMEtSumsLocation << endreq ;
	  return StatusCode::SUCCESS;
	}
      
      std::vector <LVL1::CMMEtSums>  vCMMJEMEtSums;
      std::vector <LVL1::JEMEtSums>  vJEMEtSums;
      //use standard vector instead of datavector for transmissioncheck:
      //datavector erases not only the pointer  in the vector, but also the referenced object
      //-> segmentation fault!
      
      for( it_CMMEtSums  = CMMEtSums ->begin(); it_CMMEtSums < CMMEtSums -> end(); ++it_CMMEtSums )
	{	  
	  if ((*it_CMMEtSums)->dataID()<16)
	    {
	      vCMMJEMEtSums.push_back(**it_CMMEtSums );
	    }
	}
      for( it_JEMEtSums  =  JEMEtSums->begin(); it_JEMEtSums <  JEMEtSums-> end(); ++it_JEMEtSums )
	{
	  vJEMEtSums.push_back(**it_JEMEtSums);
	}
      
      std::vector <LVL1::CMMEtSums>::iterator it_vCMMJEMEtSums;
      std::vector <LVL1::JEMEtSums>::iterator it_vJEMEtSums;
      
      it_vCMMJEMEtSums=vCMMJEMEtSums.begin();
      while (it_vCMMJEMEtSums<vCMMJEMEtSums.end())
	{
	  mLog<<MSG::VERBOSE<<"CMM Crate "<<(*it_vCMMJEMEtSums).crate()<<" Module "<<(*it_vCMMJEMEtSums).dataID()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Ex "<<(*it_vCMMJEMEtSums).Ex()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Ey "<<(*it_vCMMJEMEtSums).Ey()<<endreq;
	  mLog<<MSG::VERBOSE<<"CMM Et "<<(*it_vCMMJEMEtSums).Et()<<endreq;
	  
	  bool found = 0;
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
		      
		      if ((*it_vJEMEtSums).crate()==0)
			{
			  m_h_TransCheck_JEP->Fill(2,(*it_vJEMEtSums).module(),1);
			}
		      else 
			{
			  m_h_TransCheck_JEP->Fill(2,((*it_vJEMEtSums).module()+19),1);
			}
		      
		      
		    }
		  vCMMJEMEtSums.erase(it_vCMMJEMEtSums);
		  vJEMEtSums.erase(it_vJEMEtSums);
		  found=1;
		}
	      else it_vJEMEtSums=it_vJEMEtSums+1;
	    }
	  if (found==0)it_vCMMJEMEtSums=it_vCMMJEMEtSums+1;
	}
      
      if (vCMMJEMEtSums.size()!=0)
	{
	  mLog<<MSG::INFO<<vCMMJEMEtSums.size()<<" Transmission errors between JEM and CMM EtSums: additional CMM information"<<endreq;
	  
	  //fill errorcounter
	  for( it_vCMMJEMEtSums  = vCMMJEMEtSums.begin(); it_vCMMJEMEtSums <  vCMMJEMEtSums. end(); ++it_vCMMJEMEtSums )
	    {
	      mLog<<MSG::DEBUG<<"CMM Crate "<<(*it_vCMMJEMEtSums).crate()<<" Module "<<(*it_vCMMJEMEtSums).dataID()<<endreq;
	      mLog<<MSG::DEBUG<<"CMM Ex "<<(*it_vCMMJEMEtSums).Ex()<<endreq;
	      mLog<<MSG::DEBUG<<"CMM Ey "<<(*it_vCMMJEMEtSums).Ey()<<endreq;
	      mLog<<MSG::DEBUG<<"CMM Et "<<(*it_vCMMJEMEtSums).Et()<<endreq;
	      
	      if ((*it_vCMMJEMEtSums).crate()==0)
		{
		  m_h_TransCheck_JEP->Fill(2,(*it_vCMMJEMEtSums).dataID(),1);
		}
	      else 
		{
		  m_h_TransCheck_JEP->Fill(2,((*it_vCMMJEMEtSums).dataID()+19),1);
		}
	      
	    }
	}
      
      if (vJEMEtSums.size()!=0)
	{
	  mLog<<MSG::INFO<<vJEMEtSums.size()<<" Transmission errors between JEM and CMM EtSums: addtional JEM information"<<endreq;
	  
	  //fill errorcounter
	  for( it_vJEMEtSums  =  vJEMEtSums.begin(); it_vJEMEtSums <  vJEMEtSums. end(); ++it_vJEMEtSums )
	    {
	      mLog<<MSG::DEBUG<<"JEM Crate "<<(*it_vJEMEtSums).crate()<<" Module "<<(*it_vJEMEtSums).module()<<endreq;
	      mLog<<MSG::DEBUG<<"JEM Ex "<<(*it_vJEMEtSums).Ex()<<endreq;
	      mLog<<MSG::DEBUG<<"JEM Ey "<<(*it_vJEMEtSums).Ey()<<endreq;
	      mLog<<MSG::DEBUG<<"JEM Et "<<(*it_vJEMEtSums).Et()<<endreq;
	      
	      if ((*it_vJEMEtSums).crate()==0)
		{
		  m_h_TransCheck_JEP->Fill(2,(*it_vJEMEtSums).module(),1);
		}
	      else 
		{
		  m_h_TransCheck_JEP->Fill(2,((*it_vJEMEtSums).module()+19),1);
		}
	      
	    }
	}
    }


//******************************* 
//********* CMM RoI ************* 
//******************************* 

  LVL1::CMMRoI* CR = new LVL1::CMMRoI ;

  sc = m_storeGate->retrieve (CR, m_CMMRoILocation);
    	 
  if (sc==StatusCode::FAILURE)
    {
      mLog <<MSG::INFO<<"No CMM RoI found in TES at "<< m_CMMRoILocation<<endreq;
      return StatusCode::SUCCESS;    
    }

  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" CMM RoI ---------------"<<endreq;


  mLog<<"JetEtHits: "<<Help->Binary((CR)->jetEtHits(),4)<<"; SumEtHits: "<<Help->Binary((CR)->sumEtHits(),4)<<"; MissingEtHits: "<<Help->Binary((CR)->missingEtHits(),8)<<endreq;
  mLog<<"Ex: "<<(CR)->ex()<<"; Ey: "<<(CR)->ey()<<"; Et: "<<(CR)->et()<<endreq;

  std::string CMMRoIHit = Help->Binary((CR)-> jetEtHits(),4);
  Help->FillHitsHisto(m_h_CMMRoI_JetEtHits, CMMRoIHit, 0, 4, 0, 1, &mLog);

  CMMRoIHit = Help->Binary((CR)-> sumEtHits(),4);
  Help->FillHitsHisto(m_h_CMMRoI_SumEtHits, CMMRoIHit, 0, 4, 0, 1, &mLog);

  CMMRoIHit = Help->Binary((CR)-> missingEtHits(),8);
  Help->FillHitsHisto(m_h_CMMRoI_MissingEtHits , CMMRoIHit, 0, 8, 0, 1, &mLog);
 
  m_h_CMMRoI_Ex -> Fill( (CR)->ex(),1);
  m_h_CMMRoI_Ey -> Fill( (CR)->ey(),1);
  m_h_CMMRoI_Et -> Fill( (CR)->et(),1);
 
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

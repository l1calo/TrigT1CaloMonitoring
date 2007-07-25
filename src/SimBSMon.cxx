
// ********************************************************************
//
// NAME:        SimBSMon.cxx
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


#include "TrigT1CaloMonitoring/SimBSMon.h"
#include "TrigT1CaloMonitoring/JEMMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

//#include "TrigT1Calo/EnergyTrigger.h"
#include "TrigT1Calo/LVL1TriggerMenuDefs.h"
#include "TrigT1Calo/LVL1TriggerMenu.h"
#include "TrigT1Calo/InternalJetROI.h"
#include "TrigT1Calo/JEMRoI.h"
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
SimBSMon::SimBSMon( const std::string & type, const std::string & name,
		const IInterface* parent )
  : ManagedMonitorToolBase( type, name, parent )
/*---------------------------------------------------------*/
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "BS_CMMJetHitsLocation", m_BS_CMMJetHitsLocation =  LVL1::TrigT1CaloDefs::CMMJetHitsLocation) ;
  declareProperty( "BS_CMMEtSumsLocation", m_BS_CMMEtSumsLocation =  LVL1::TrigT1CaloDefs::CMMEtSumsLocation) ;  
  declareProperty( "BS_CMMRoILocation", m_BS_CMMRoILocation =  LVL1::TrigT1CaloDefs::CMMRoILocation) ;

  declareProperty( "BS_JEMHitsLocation", m_BS_JEMHitsLocation =  LVL1::TrigT1CaloDefs::JEMHitsLocation) ;
  declareProperty( "BS_JEMEtSumsLocation", m_BS_JEMEtSumsLocation=   LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;
  declareProperty( "BS_JEMRoILocation", m_BS_JEMRoILocation=   LVL1::TrigT1CaloDefs::JEMRoILocation) ;

  declareProperty( "Sim_CMMJetHitsLocation", m_Sim_CMMJetHitsLocation =  LVL1::TrigT1CaloDefs::CMMJetHitsLocation) ;
  declareProperty( "Sim_CMMEtSumsLocation", m_Sim_CMMEtSumsLocation =  LVL1::TrigT1CaloDefs::CMMEtSumsLocation) ;  
  declareProperty( "Sim_CMMRoILocation", m_Sim_CMMRoILocation =  LVL1::TrigT1CaloDefs::CMMRoILocation) ;

  declareProperty( "Sim_JEMHitsLocation", m_Sim_JEMHitsLocation =  LVL1::TrigT1CaloDefs::JEMHitsLocation) ;
  declareProperty( "Sim_JEMEtSumsLocation", m_Sim_JEMEtSumsLocation=   LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;
  declareProperty( "Sim_JEMRoILocation", m_Sim_JEMRoILocation=   LVL1::TrigT1CaloDefs::JEMRoILocation) ;

  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/CMM") ;
  declareProperty( "DataType", m_DataType="") ;
}


/*---------------------------------------------------------*/
SimBSMon::~SimBSMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode SimBSMon::bookHistograms( bool isNewEventsBlock, 
				   bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in SimBSMon::bookHistograms" << endreq;
  
  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  
  if (sc.isFailure()) 
    {
      mLog << MSG::ERROR << "Unable to retrieve pointer to StoreGateSvc" << endreq;
      return sc;
    }
  MonGroup SimBSComparison_JEM ( this, (m_PathInRootFile).c_str(), shift, eventsBlock );
  HistoBooker* JEM_Booker = new HistoBooker(&SimBSComparison_JEM, &mLog, m_DataType);

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }
  
  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
  
  m_h_SimBSMon_JEP=JEM_Booker->book2F("JEP_Calc_Error", "JEP Hardware Output compared to Simulation: Differences", 3,0.5,3.5,36,-0.5,35.5, "Calculated Magnitude", "");
  //m_h_SimBSMon_JEP-> SetOption ("text");

  m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(1, "Energy");
  m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(2, "Hits");
  m_h_SimBSMon_JEP->GetXaxis()->SetBinLabel(3, "RoI");

  std::string name;
  std::stringstream buffer;
    
  for (int i = 0; i < 16; i++)
    {
      buffer.str("");
      buffer<<i;
      
      name = "JEM " + buffer.str();
      m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel((i+1), name.c_str());
      
      buffer.str("");
      buffer<<i;
      
      name = "JEM " + buffer.str();
      m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel((i+1+18), name.c_str());
    }
      m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel(17, "CMM");
      m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel(18, "Crate 1: ");
      m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel(35, "CMM");
      m_h_SimBSMon_JEP->GetYaxis()->SetBinLabel(36, "Crate 0: ");

    }

  if( isNewRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}


/*---------------------------------------------------------*/
StatusCode SimBSMon::fillHistograms()
  /*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  Helper* Help = new Helper();
  
  //*******************************************************
  //************** JEMEtSums ******************************
  //*******************************************************

  const JEMEtSumsCollection* BS_JEMEtSums;
  StatusCode sc = m_storeGate->retrieve(BS_JEMEtSums, m_BS_JEMEtSumsLocation);
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

  std::vector <LVL1::JEMEtSums>  vBS_JEMEtSums;
  std::vector <LVL1::JEMEtSums>  vSim_JEMEtSums;
  //use standard vector instead of datavector for transmissioncheck:
  //datavector erases not only the pointer  in the vector, but also the referenced object
  //-> segmentation fault!

  //fill vectors
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
    bool foundModule;
    std::vector <LVL1::JEMEtSums>::iterator it_BS_JEMEtSums;
    std::vector <LVL1::JEMEtSums>::iterator it_Sim_JEMEtSums;

    it_BS_JEMEtSums=vBS_JEMEtSums.begin();

    while (it_BS_JEMEtSums<vBS_JEMEtSums.end())
      {
	mLog<<MSG::VERBOSE<<"BS: Crate "<<(*it_BS_JEMEtSums).crate()<<" Module "<<(*it_BS_JEMEtSums).module()<<endreq;
	mLog<<MSG::VERBOSE<<"BS: Et "<<(*it_BS_JEMEtSums).Et()<<endreq;
      
	foundModule = 0;
	it_Sim_JEMEtSums=vSim_JEMEtSums.begin();

	while ((foundModule==0)and(it_Sim_JEMEtSums<vSim_JEMEtSums.end()))
	  {
	    mLog<<MSG::VERBOSE<<"Sim: Crate "<<(*it_Sim_JEMEtSums).crate()<<" Modul "<<(*it_Sim_JEMEtSums).module()<<endreq;
	    mLog<<MSG::VERBOSE<<"Sim: Et "<<(*it_Sim_JEMEtSums).Et()<<endreq;


	    if (((*it_BS_JEMEtSums).crate()==(*it_Sim_JEMEtSums).crate())
		and((*it_BS_JEMEtSums).module()==(*it_Sim_JEMEtSums).module()))
	      {
		foundModule=1;

		if (((*it_BS_JEMEtSums).Ex()!=(*it_Sim_JEMEtSums).Ex())
		    or ((*it_BS_JEMEtSums).Ey()!=(*it_Sim_JEMEtSums).Ey())
		    or ((*it_BS_JEMEtSums).Et()!=(*it_Sim_JEMEtSums).Et()))
		  {
		    mLog<<MSG::INFO<<"JEMEtSums Calculation difference between ByteStream and Simulation JEM Energy"<<endreq;
		    mLog<<MSG::INFO<<"Crate "<<(*it_BS_JEMEtSums).crate()<<" Module "<<(*it_BS_JEMEtSums).module()<<endreq;
		    mLog<<MSG::DEBUG<<"BS: Ex (compressed)"<<(*it_BS_JEMEtSums).Ex()<<endreq;
		    mLog<<MSG::DEBUG<<"BS: Ey (compressed)"<<(*it_BS_JEMEtSums).Ey()<<endreq;
		    mLog<<MSG::DEBUG<<"BS: Et (compressed)"<<(*it_BS_JEMEtSums).Et()<<endreq;

		    mLog<<MSG::DEBUG<<"Sim: Ex (compressed)"<<(*it_Sim_JEMEtSums).Ex()<<endreq;
		    mLog<<MSG::DEBUG<<"Sim: Ey (compressed)"<<(*it_Sim_JEMEtSums).Ey()<<endreq;
		    mLog<<MSG::DEBUG<<"Sim: Et (compressed)"<<(*it_Sim_JEMEtSums).Et()<<endreq;

		    if ((*it_BS_JEMEtSums).crate()==0) 
		    {
		      m_h_SimBSMon_JEP->Fill(1,((*it_BS_JEMEtSums).module()+18),1);
		    }
		    else 
		    {
		      m_h_SimBSMon_JEP->Fill(1,((*it_BS_JEMEtSums).module()+18),1);
		    }
		  }

		vBS_JEMEtSums.erase(it_BS_JEMEtSums);
		vSim_JEMEtSums.erase(it_Sim_JEMEtSums);
	      }
	    else it_Sim_JEMEtSums=it_Sim_JEMEtSums+1;
	  }
	if (foundModule==0)it_BS_JEMEtSums=it_BS_JEMEtSums+1;
      }

    if (vBS_JEMEtSums.size()!=0)
      {
	mLog<<MSG::INFO<<"JEMEtSums No Match found between BS and Sim for"<<endreq;

	//fill errorcounter
	for( it_BS_JEMEtSums  = vBS_JEMEtSums.begin(); it_BS_JEMEtSums <  vBS_JEMEtSums. end(); ++it_BS_JEMEtSums )
	  {
	    mLog<<MSG::INFO<<"BS: Crate "<<(*it_BS_JEMEtSums).crate()<<" Module "<<(*it_BS_JEMEtSums).module()<<endreq;
	    mLog<<MSG::DEBUG<<"BS: Ex (compressed)"<<(*it_BS_JEMEtSums).Ex()<<endreq;
	    mLog<<MSG::DEBUG<<"BS: Ey (compressed)"<<(*it_BS_JEMEtSums).Ey()<<endreq;
	    mLog<<MSG::DEBUG<<"BS: Et (compressed)"<<(*it_BS_JEMEtSums).Et()<<endreq;
	    
	    if ((*it_BS_JEMEtSums).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(1,((*it_BS_JEMEtSums).module()+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(1,((*it_BS_JEMEtSums).module()),1);
	      }
	  }
      }

    if (vSim_JEMEtSums.size()!=0)
      {
	mLog<<MSG::INFO<<"JEMEtSums No Match found between BS and Sim for"<<endreq;

	//fill errorcounter
	for( it_Sim_JEMEtSums  = vSim_JEMEtSums.begin(); it_Sim_JEMEtSums <  vSim_JEMEtSums. end(); ++it_Sim_JEMEtSums )
	  {
	    mLog<<MSG::INFO<<"Sim: Crate "<<(*it_Sim_JEMEtSums).crate()<<" Module "<<(*it_Sim_JEMEtSums).module()<<endreq;
	    mLog<<MSG::DEBUG<<"Sim: Ex (compressed)"<<(*it_Sim_JEMEtSums).Ex()<<endreq;
	    mLog<<MSG::DEBUG<<"Sim: Ey (compressed)"<<(*it_Sim_JEMEtSums).Ey()<<endreq;
	    mLog<<MSG::DEBUG<<"Sim: Et (compressed)"<<(*it_Sim_JEMEtSums).Et()<<endreq;
	    
	    if ((*it_BS_JEMEtSums).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(1,((*it_BS_JEMEtSums).module()+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(1,((*it_BS_JEMEtSums).module()),1);
	      }
	  }
      }
  //*******************************************************
  //************** CMMEtSums ******************************
  //*******************************************************

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
  CMMEtSumsCollection::const_iterator it_CMMEtSums ;

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

    while (it_BS_CMMEtSums<vBS_CMMEtSums.end())
      {
	mLog<<MSG::VERBOSE<<"BS: Crate "<<(*it_BS_CMMEtSums).crate()<<" dataID "<<(*it_BS_CMMEtSums).dataID()<<endreq;
	mLog<<MSG::VERBOSE<<"BS: Ex (compressed)"<<(*it_BS_CMMEtSums).Ex()<<endreq;
	mLog<<MSG::VERBOSE<<"BS: Ey (compressed)"<<(*it_BS_CMMEtSums).Ey()<<endreq;
	mLog<<MSG::VERBOSE<<"BS: Et (compressed)"<<(*it_BS_CMMEtSums).Et()<<endreq;
     
	foundModule = 0;
	it_Sim_CMMEtSums=vSim_CMMEtSums.begin();

	while ((foundModule==0)and(it_Sim_CMMEtSums<vSim_CMMEtSums.end()))
	  {
	    mLog<<MSG::VERBOSE<<"Sim: Crate "<<(*it_Sim_CMMEtSums).crate()<<" DataID "<<(*it_Sim_CMMEtSums).dataID()<<endreq;
	    mLog<<MSG::VERBOSE<<"Sim: Ex (compressed)"<<(*it_Sim_CMMEtSums).Ex()<<endreq;
	    mLog<<MSG::VERBOSE<<"Sim: Ey (compressed)"<<(*it_Sim_CMMEtSums).Ey()<<endreq;
	    mLog<<MSG::VERBOSE<<"Sim: Et (compressed)"<<(*it_Sim_CMMEtSums).Et()<<endreq;

	    if (((*it_BS_CMMEtSums).crate()==(*it_Sim_CMMEtSums).crate())
		and((*it_BS_CMMEtSums).dataID()==(*it_Sim_CMMEtSums).dataID()))
	      {
		foundModule=1;

		if (((*it_BS_CMMEtSums).Ex()!=(*it_Sim_CMMEtSums).Ex())
		    or ((*it_BS_CMMEtSums).Ey()!=(*it_Sim_CMMEtSums).Ey())
		    or ((*it_BS_CMMEtSums).Et()!=(*it_Sim_CMMEtSums).Et()))
		  {
		    mLog<<MSG::INFO<<"CMMEtSums Calculation difference between ByteStream and Simulation CMM Energy"<<endreq;
		    mLog<<MSG::INFO<<"Crate "<<(*it_BS_CMMEtSums).crate()<<" DataId "<<(*it_BS_CMMEtSums).dataID()<<endreq;
		    mLog<<MSG::DEBUG<<"BS: Ex (compressed)"<<(*it_BS_CMMEtSums).Ex()<<endreq;
		    mLog<<MSG::DEBUG<<"BS: Ey (compressed)"<<(*it_BS_CMMEtSums).Ey()<<endreq;
		    mLog<<MSG::DEBUG<<"BS: Et (compressed)"<<(*it_BS_CMMEtSums).Et()<<endreq;

		    mLog<<MSG::DEBUG<<"Sim: Ex (compressed)"<<(*it_Sim_CMMEtSums).Ex()<<endreq;
		    mLog<<MSG::DEBUG<<"Sim: Ey (compressed)"<<(*it_Sim_CMMEtSums).Ey()<<endreq;
		    mLog<<MSG::DEBUG<<"Sim: Et (compressed)"<<(*it_Sim_CMMEtSums).Et()<<endreq;

		    if ((*it_BS_JEMEtSums).crate()==0) 
		      {
			m_h_SimBSMon_JEP->Fill(1,(16+18),1);
		      }
		    else 
		      {
			m_h_SimBSMon_JEP->Fill(1,16,1);
		      }
		    
		  }

		vBS_CMMEtSums.erase(it_BS_CMMEtSums);
		vSim_CMMEtSums.erase(it_Sim_CMMEtSums);
	      }
	    else it_Sim_CMMEtSums=it_Sim_CMMEtSums+1;
	  }
	if (foundModule==0)it_BS_CMMEtSums=it_BS_CMMEtSums+1;
      }

    if (vBS_CMMEtSums.size()!=0)
      {
	mLog<<MSG::INFO<<"CMMEtSums No Match found between BS and Sim for"<<endreq;

	//fill errorcounter
	for( it_BS_CMMEtSums  = vBS_CMMEtSums.begin(); it_BS_CMMEtSums <  vBS_CMMEtSums. end(); ++it_BS_CMMEtSums )
	  {
	    mLog<<MSG::INFO<<"BS: Crate "<<(*it_BS_CMMEtSums).crate()<<" DataId "<<(*it_BS_CMMEtSums).dataID()<<endreq;
	    mLog<<MSG::DEBUG<<"BS: Ex (compressed)"<<(*it_BS_CMMEtSums).Ex()<<endreq;
	    mLog<<MSG::DEBUG<<"BS: Ey (compressed)"<<(*it_BS_CMMEtSums).Ey()<<endreq;
	    mLog<<MSG::DEBUG<<"BS: Et (compressed)"<<(*it_BS_CMMEtSums).Et()<<endreq;
	    
	    if ((*it_BS_JEMEtSums).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(1,(16+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(1,16,1);
	      }

	  }
      }

    if (vSim_CMMEtSums.size()!=0)
      {
	mLog<<MSG::INFO<<"CMMEtSums No Match found between BS and Sim for"<<endreq;

	//fill errorcounter
	for( it_Sim_CMMEtSums  = vSim_CMMEtSums.begin(); it_Sim_CMMEtSums <  vSim_CMMEtSums. end(); ++it_Sim_CMMEtSums )
	  {
	    mLog<<MSG::INFO<<"Sim: Crate "<<(*it_Sim_CMMEtSums).crate()<<" DataId "<<(*it_Sim_CMMEtSums).dataID()<<endreq;
	    mLog<<MSG::DEBUG<<"Sim: Ex (compressed)"<<(*it_Sim_CMMEtSums).Ex()<<endreq;
	    mLog<<MSG::DEBUG<<"Sim: Ey (compressed)"<<(*it_Sim_CMMEtSums).Ey()<<endreq;
	    mLog<<MSG::DEBUG<<"Sim: Et (compressed)"<<(*it_Sim_CMMEtSums).Et()<<endreq;

	    if ((*it_BS_JEMEtSums).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(1,(16+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(1,16,1);
	      }

	  }
      }

  //*******************************************************
  //************** JEMHits ******************************
  //*******************************************************
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

    while (it_BS_JEMHits<vBS_JEMHits.end())
      {
	mLog<<MSG::VERBOSE<<"BS: Crate "<<(*it_BS_JEMHits).crate()<<" Module "<<(*it_BS_JEMHits).module()<<endreq;
	mLog<<MSG::VERBOSE<<"BS: Hits "<<Help->Binary((*it_BS_JEMHits).JetHits(),24)<<endreq;
      
	foundModule = 0;
	it_Sim_JEMHits=vSim_JEMHits.begin();

	while ((foundModule==0)and(it_Sim_JEMHits<vSim_JEMHits.end()))
	  {
	    mLog<<MSG::VERBOSE<<"Sim: Crate "<<(*it_Sim_JEMHits).crate()<<" Modul "<<(*it_Sim_JEMHits).module()<<endreq;
	    mLog<<MSG::VERBOSE<<"Sim: Hits "<<Help->Binary((*it_Sim_JEMHits).JetHits(),24)<<endreq;


	    if (((*it_BS_JEMHits).crate()==(*it_Sim_JEMHits).crate())
		and((*it_BS_JEMHits).module()==(*it_Sim_JEMHits).module()))
	      {
		foundModule=1;

		if ((*it_BS_JEMHits).JetHits()!=(*it_Sim_JEMHits).JetHits())
		  {
		    mLog<<MSG::INFO<<"Hits Calculation difference between ByteStream and Simulation JEM Hits"<<endreq;
		    mLog<<MSG::INFO<<"Crate "<<(*it_BS_JEMHits).crate()<<" Module "<<(*it_BS_JEMHits).module()<<endreq;
		    mLog<<MSG::DEBUG<<"BS: Hits "<<Help->Binary((*it_BS_JEMHits).JetHits(),24)<<endreq;
		    mLog<<MSG::DEBUG<<"Sim: Hits "<<Help->Binary((*it_Sim_JEMHits).JetHits(),24)<<endreq;
	
		    if ((*it_BS_JEMHits).crate()==0) 
		      {
			m_h_SimBSMon_JEP->Fill(2,((*it_BS_JEMHits).module()+18),1);
		      }
		    else 
		      {
			m_h_SimBSMon_JEP->Fill(2,(*it_BS_JEMHits).module(),1);
		      }

		  }

		vBS_JEMHits.erase(it_BS_JEMHits);
		vSim_JEMHits.erase(it_Sim_JEMHits);
	      }
	    else it_Sim_JEMHits=it_Sim_JEMHits+1;
	  }
	if (foundModule==0)it_BS_JEMHits=it_BS_JEMHits+1;
      }

    if (vBS_JEMHits.size()!=0)
      {
	mLog<<MSG::INFO<<"JEMHits No Match found in JEMHits between BS and Sim for "<<endreq;

	//fill errorcounter
	for( it_BS_JEMHits  = vBS_JEMHits.begin(); it_BS_JEMHits <  vBS_JEMHits. end(); ++it_BS_JEMHits )
	  {
	    mLog<<MSG::INFO<<"BS: Crate "<<(*it_BS_JEMHits).crate()<<" Module "<<(*it_BS_JEMHits).module()<<endreq;
	    mLog<<MSG::DEBUG<<"BS: Hits "<<Help->Binary((*it_BS_JEMHits).JetHits(),24)<<endreq;
	    
	    if ((*it_BS_JEMHits).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(2,((*it_BS_JEMHits).module()+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(2,(*it_BS_JEMHits).module(),1);
	      }
	    
	  }
      }

    if (vSim_JEMHits.size()!=0)
      {
	mLog<<MSG::INFO<<"JEMHits No Match found between BS and Sim for"<<endreq;

	//fill errorcounter
	for( it_Sim_JEMHits  = vSim_JEMHits.begin(); it_Sim_JEMHits <  vSim_JEMHits. end(); ++it_Sim_JEMHits )
	  {
	    mLog<<MSG::INFO<<"Sim: Crate "<<(*it_Sim_JEMHits).crate()<<" Module "<<(*it_Sim_JEMHits).module()<<endreq;
	    mLog<<MSG::DEBUG<<"Sim: Hits "<<Help->Binary((*it_Sim_JEMHits).JetHits(),24)<<endreq;
	    
	    if ((*it_BS_JEMHits).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(2,((*it_BS_JEMHits).module()+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(2,(*it_BS_JEMHits).module(),1);
	      }
	    
	  }
      }

  //*******************************************************
  //************** CMMJetHits ******************************
  //*******************************************************

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
  CMMJetHitsCollection::const_iterator it_CMMJetHits ;

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

    while (it_BS_CMMJetHits<vBS_CMMJetHits.end())
      {
	mLog<<MSG::VERBOSE<<"BS: Crate "<<(*it_BS_CMMJetHits).crate()<<" dataID "<<(*it_BS_CMMJetHits).dataID()<<endreq;
	mLog<<MSG::VERBOSE<<"BS: Hits"<<Help->Binary((*it_BS_CMMJetHits).Hits(),24)<<endreq;
     
	foundModule = 0;
	it_Sim_CMMJetHits=vSim_CMMJetHits.begin();

	while ((foundModule==0)and(it_Sim_CMMJetHits<vSim_CMMJetHits.end()))
	  {
	    mLog<<MSG::VERBOSE<<"Sim: Crate "<<(*it_Sim_CMMJetHits).crate()<<" DataID "<<(*it_Sim_CMMJetHits).dataID()<<endreq;
	    mLog<<MSG::VERBOSE<<"Sim: Hits"<<Help->Binary((*it_Sim_CMMJetHits).Hits(),24)<<endreq;

	    if (((*it_BS_CMMJetHits).crate()==(*it_Sim_CMMJetHits).crate())
		and((*it_BS_CMMJetHits).dataID()==(*it_Sim_CMMJetHits).dataID()))
	      {
		foundModule=1;

		if ((*it_BS_CMMJetHits).Hits()!=(*it_Sim_CMMJetHits).Hits())
		  {
		    mLog<<MSG::INFO<<"CMMJetHits Calculation difference between ByteStream and Simulation CMM Hits"<<endreq;
		    mLog<<MSG::INFO<<"Crate "<<(*it_BS_CMMJetHits).crate()<<" DataId "<<(*it_BS_CMMJetHits).dataID()<<endreq;
		    mLog<<MSG::DEBUG<<"BS: Hits"<<Help->Binary((*it_BS_CMMJetHits).Hits(),24)<<endreq;

		    mLog<<MSG::DEBUG<<"Sim: Hits"<<Help->Binary((*it_Sim_CMMJetHits).Hits(),24)<<endreq;

		    if ((*it_BS_JEMHits).crate()==0) 
		      {
			m_h_SimBSMon_JEP->Fill(2,(16+18),1);
		      }
		    else 
		      {
			m_h_SimBSMon_JEP->Fill(2,16,1);
		      }

		  }

		vBS_CMMJetHits.erase(it_BS_CMMJetHits);
		vSim_CMMJetHits.erase(it_Sim_CMMJetHits);
	      }
	    else it_Sim_CMMJetHits=it_Sim_CMMJetHits+1;
	  }
	if (foundModule==0)it_BS_CMMJetHits=it_BS_CMMJetHits+1;
      }

    if (vBS_CMMJetHits.size()!=0)
      {
	mLog<<MSG::INFO<<"CMMJetHits No Match found between BS and Sim for"<<endreq;

	//fill errorcounter
	for( it_BS_CMMJetHits  = vBS_CMMJetHits.begin(); it_BS_CMMJetHits <  vBS_CMMJetHits. end(); ++it_BS_CMMJetHits )
	  {
	    mLog<<MSG::INFO<<"BS: Crate "<<(*it_BS_CMMJetHits).crate()<<" DataId "<<(*it_BS_CMMJetHits).dataID()<<endreq;
	    mLog<<MSG::DEBUG<<"BS: Hits"<<Help->Binary((*it_BS_CMMJetHits).Hits(),24)<<endreq;
	    
	    if ((*it_BS_JEMHits).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(2,(16+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(2,16,1);
	      }
	    
	  }
      }

    if (vSim_CMMJetHits.size()!=0)
      {
	mLog<<MSG::INFO<<"CMMJetHits No Match found between BS and Sim for"<<endreq;

	//fill errorcounter
	for( it_Sim_CMMJetHits  = vSim_CMMJetHits.begin(); it_Sim_CMMJetHits <  vSim_CMMJetHits. end(); ++it_Sim_CMMJetHits )
	  {
	    mLog<<MSG::INFO<<"Sim: Crate "<<(*it_Sim_CMMJetHits).crate()<<" DataId "<<(*it_Sim_CMMJetHits).dataID()<<endreq;
	    mLog<<MSG::DEBUG<<"Sim: Hits"<<Help->Binary((*it_Sim_CMMJetHits).Hits(),24)<<endreq;
	    
	    if ((*it_BS_JEMHits).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(2,(16+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(2,16,1);
	      }

	  }
      }

  //*******************************************************
  //************** JEM RoI ********************************
  //*******************************************************
  const JEMRoICollection* BS_JEMRoI;
  sc = m_storeGate->retrieve(BS_JEMRoI, m_BS_JEMRoILocation);
  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMRoI found in TES at "<< m_BS_JEMRoILocation << endreq ;
      return StatusCode::SUCCESS;
    }

  const JEMRoICollection* Sim_JEMRoI;
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
  JEMRoICollection::const_iterator it_JEMRoI ;

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

    while (it_BS_JEMRoI<vBS_JEMRoI.end())
      {
	mLog<<MSG::VERBOSE<<"BS: Crate "<<(*it_BS_JEMRoI).crate()<<" Module "<<(*it_BS_JEMRoI).jem()<<endreq;
	mLog<<MSG::VERBOSE<<"BS: RoI "<<Help->Binary((*it_BS_JEMRoI).hits(),8)<<endreq;
      
	foundModule = 0;
	it_Sim_JEMRoI=vSim_JEMRoI.begin();

	while ((foundModule==0)and(it_Sim_JEMRoI<vSim_JEMRoI.end()))
	  {
	    mLog<<MSG::VERBOSE<<"Sim: Crate "<<(*it_Sim_JEMRoI).crate()<<" Modul "<<(*it_Sim_JEMRoI).jem()<<endreq;
	    mLog<<MSG::VERBOSE<<"Sim: RoI "<<Help->Binary((*it_Sim_JEMRoI).hits(),8)<<endreq;


	    if (((*it_BS_JEMRoI).crate()==(*it_Sim_JEMRoI).crate())
		and((*it_BS_JEMRoI).jem()==(*it_Sim_JEMRoI).jem()))
	      {
		foundModule=1;

		if ((*it_BS_JEMRoI).hits()!=(*it_Sim_JEMRoI).hits())
		  {
		    mLog<<MSG::INFO<<"JEMRoI Calculation difference between ByteStream and Simulation JEM RoI"<<endreq;
		    mLog<<MSG::INFO<<"Crate "<<(*it_BS_JEMRoI).crate()<<" Module "<<(*it_BS_JEMRoI).jem()<<endreq;
		    mLog<<MSG::DEBUG<<"BS: RoI "<<Help->Binary((*it_BS_JEMRoI).hits(),8)<<endreq;
		    mLog<<MSG::DEBUG<<"Sim: RoI "<<Help->Binary((*it_Sim_JEMRoI).hits(),8)<<endreq;
	
		    if ((*it_BS_JEMRoI).crate()==0) 
		      {
			m_h_SimBSMon_JEP->Fill(3,((*it_BS_JEMRoI).jem()+18),1);
		      }
		    else 
		      {
			m_h_SimBSMon_JEP->Fill(3,(*it_BS_JEMRoI).jem(),1);
		      }
		    
		  }

		vBS_JEMRoI.erase(it_BS_JEMRoI);
		vSim_JEMRoI.erase(it_Sim_JEMRoI);
	      }
	    else it_Sim_JEMRoI=it_Sim_JEMRoI+1;
	  }
	if (foundModule==0)it_BS_JEMRoI=it_BS_JEMRoI+1;
      }

    if (vBS_JEMRoI.size()!=0)
      {
	mLog<<MSG::INFO<<"JEMRoI No Match found in JEMRoI between BS and Sim for "<<endreq;

	//fill errorcounter
	for( it_BS_JEMRoI  = vBS_JEMRoI.begin(); it_BS_JEMRoI <  vBS_JEMRoI. end(); ++it_BS_JEMRoI )
	  {
	    mLog<<MSG::INFO<<"BS: Crate "<<(*it_BS_JEMRoI).crate()<<" Module "<<(*it_BS_JEMRoI).jem()<<endreq;
	    mLog<<MSG::DEBUG<<"BS: RoI "<<Help->Binary((*it_BS_JEMRoI).hits(),8)<<endreq;
	    
	    if ((*it_BS_JEMRoI).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(3,((*it_BS_JEMRoI).jem()+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(3,(*it_BS_JEMRoI).jem(),1);
	      }
	    
	  }
      }

    if (vSim_JEMRoI.size()!=0)
      {
	mLog<<MSG::INFO<<"JEMRoI No Match found between BS and Sim for"<<endreq;

	//fill errorcounter
	for( it_Sim_JEMRoI  = vSim_JEMRoI.begin(); it_Sim_JEMRoI <  vSim_JEMRoI. end(); ++it_Sim_JEMRoI )
	  {
	    mLog<<MSG::INFO<<"Sim: Crate "<<(*it_Sim_JEMRoI).crate()<<" Module "<<(*it_Sim_JEMRoI).jem()<<endreq;
	    mLog<<MSG::DEBUG<<"Sim: RoI "<<Help->Binary((*it_Sim_JEMRoI).hits(),8)<<endreq;
	    
	    if ((*it_BS_JEMRoI).crate()==0) 
	      {
		m_h_SimBSMon_JEP->Fill(3,((*it_BS_JEMRoI).jem()+18),1);
	      }
	    else 
	      {
		m_h_SimBSMon_JEP->Fill(3,(*it_BS_JEMRoI).jem(),1);
	      }

	  }
      }

  //*******************************************************
  //************** CMM RoI ********************************
  //*******************************************************





  return StatusCode( StatusCode::SUCCESS );
}

/*---------------------------------------------------------*/
StatusCode SimBSMon::procHistograms( bool isEndOfEventsBlock, 
				  bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
        if( isEndOfEventsBlock || isEndOfLumiBlock ) 
	  {

	}
	
	if( isEndOfRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}

// ********************************************************************
//
// NAME:        JEMMon.cxx
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on JEM level
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

#include "TrigT1CaloMonitoring/JEMMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

//#include "TrigT1Calo/JEMEtSums.h"
#include "TrigT1Calo/EnergyTrigger.h"
#include "TrigT1Calo/JEMRoI.h"

#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "AthenaMonitoring/AthenaMonManager.h"


// *********************************************************************
// Public Methods
// *********************************************************************

/*---------------------------------------------------------*/
JEMMon::JEMMon( const std::string & type, const std::string & name,
		const IInterface* parent )
  : ManagedMonitorToolBase( type, name, parent )
/*---------------------------------------------------------*/
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "JEMHitsLocation", m_JEMHitsLocation =  LVL1::TrigT1CaloDefs::JEMHitsLocation) ;
  declareProperty( "JEMEtSumsLocation", m_JEMEtSumsLocation=   LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;
  declareProperty( "JEMRoILocation", m_JEMRoILocation =  LVL1::TrigT1CaloDefs::JEMRoILocation) ;

  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/JEM") ;
  declareProperty( "DataType", m_DataType="") ;

}


/*---------------------------------------------------------*/
JEMMon::~JEMMon()
/*---------------------------------------------------------*/
{
}


/*---------------------------------------------------------*/
StatusCode JEMMon::bookHistograms( bool isNewEventsBlock, 
				   bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in JEMMon::bookHistograms" << endreq;

  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) 
    {
      mLog << MSG::ERROR
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
	
  MonGroup JEM ( this, m_PathInRootFile, expert, eventsBlock );
  HistoBooker* expert_Booker = new HistoBooker(&JEM, &mLog, m_DataType);

  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
	  
      m_h_JEMHits_MainHits = expert_Booker->book1F("MainHits", "JEM  --  Main Hits)", 8, -0.5, 7.5,  "Threshold No.", "#");
      m_h_JEMHits_FwdHitsRight = expert_Booker->book1F("FwdHits_Right", "JEM  --  Forward Hits Right", 4, -0.5, 3.5, "Threshold No.", "#");
      m_h_JEMHits_FwdHitsLeft = expert_Booker->book1F("FwdHits_Left", "JEM  --  Forward HitsLeft", 4, -0.5, 3.5, "Threshold No.", "#");

      m_h_JEMEtSums_Ex = expert_Booker->book1F("Ex", "JEM  --  Ex", 250, 0,250, "Ex [GeV]", "#");
      m_h_JEMEtSums_Ey = expert_Booker->book1F("Ey", "JEM  --  Ey", 250, 0,250, "Ex [GeV]", "#");
      m_h_JEMEtSums_Et = expert_Booker->book1F("Et", "JEM  --  Et", 250, 0,250, "Ex [GeV]", "#");

      // JEM RoI
      m_h_JEMRoI_Hits = expert_Booker->book1F("Hits","JEM  --  RoI Hits",8, -0.5,7.5, "Threshold No." , "Number");

    }
  
  if( isNewRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}



/*---------------------------------------------------------*/
StatusCode JEMMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  Helper* Help = new Helper();

  //******************************* 
  //********* JEM Hits  *********** 
  //******************************* 

  // fill  Histos
  const JEMHitsCollection* JEMHits;
  StatusCode sc = m_storeGate->retrieve(JEMHits, m_JEMHitsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No JEMHits found in TES at "
	   << m_JEMHitsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<endreq;
  mLog<<MSG::DEBUG<<"-------------- "<<m_DataType<<" JEM Hits ---------------"<<endreq;

  // Step over all cells and put into hist
  JEMHitsCollection::const_iterator it_JEMHits ;

  for( it_JEMHits  = JEMHits ->begin(); it_JEMHits < JEMHits -> end(); ++it_JEMHits )
    {	  
      int j=0;
      int JetMult=0;
      std::string JEMHit = Help->Binary((*it_JEMHits)-> JetHits(),24);
      // the binary hit information is represented by an integer number, that is converted to a string in order
      // to get the real binary information
      // later the multiplicities of the several thresholds are retrieved from this string

      mLog<<MSG::DEBUG<<" JEMHits Crate: "<< (*it_JEMHits)->crate()<<"  Module: "<<(*it_JEMHits)->module()
	  << "  JetHits: "<<(*it_JEMHits)-> JetHits()
	  << "   Hits(binary): "  <<JEMHit<<   endreq;
      
      if  ((*it_JEMHits)->forward()==0) //Main Jets
	{
	  for(j=0;j<8;j++)
	    {
	      JetMult=Help->Multiplicity(JEMHit,j,3);

	      mLog<<MSG::DEBUG<<" MainThresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_JEMHits_MainHits -> Fill( j, JetMult);
	    }
	}

      if  ((*it_JEMHits)->forward()==1) //fwd jets a bit complicated!
	// fwd and main hits are contained in the same hitword
	{
	  for(j=0;j<4;j++)  //fwd hits
	    {
	      JetMult=Help->Multiplicity(JEMHit,j,2); // only 2 bits per thresh

	      if (((*it_JEMHits)-> module()==0) or((*it_JEMHits)-> module()==8) )//left fwd hits
		// JEMs No 0 and 8 are processing forward left hits,
		// JEMs No 7 and 15 forward right hits
		{
		  mLog<<MSG::DEBUG<<" LeftFwd Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
		  m_h_JEMHits_FwdHitsLeft -> Fill( j, JetMult);
		}
	      if (((*it_JEMHits)-> module()==7) or((*it_JEMHits)-> module()==15) )//right fwd hits
		{
		  mLog<<MSG::DEBUG<<" RightFwd Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
		  m_h_JEMHits_FwdHitsRight -> Fill( j, JetMult);
		}


	    }
	  for(j=4;j<12;j++) //main hits
	    {
	      JetMult=Help->Multiplicity(JEMHit,j,2);

	      mLog<<MSG::DEBUG<<" Main Thresh.No: "<<(j-4)<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_JEMHits_MainHits -> Fill( (j-4), JetMult);
	    }

	}


    }     	
 
  //******************************* 
  //********* JEM Et Sums ********* 
  //******************************* 

  //fill  Histos
  const JEMEtSumsCollection* JEMEtSums;
  sc = m_storeGate->retrieve(JEMEtSums, m_JEMEtSumsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No JEMEtSums found in TES at "
	   << m_JEMEtSumsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<"-------------- JEM Et Sums ---------------"<<endreq;

  // Step over all cells and put into hist
  JEMEtSumsCollection::const_iterator it_JEMEtSums ;

  for( it_JEMEtSums  = JEMEtSums ->begin(); it_JEMEtSums < JEMEtSums -> end(); ++it_JEMEtSums )
    {	  
     
     // fill histograms
      m_h_JEMEtSums_Ex -> Fill( (*it_JEMEtSums)-> Ex(), 1.); 
      m_h_JEMEtSums_Ey -> Fill( (*it_JEMEtSums)-> Ey(), 1.); 
      m_h_JEMEtSums_Et -> Fill( (*it_JEMEtSums)-> Et() , 1.); 
      mLog <<MSG::DEBUG<< " JEMEtSums Crate: "<<(*it_JEMEtSums)->crate()<<"  Module: "<<(*it_JEMEtSums)->module()
	   <<"   Ex: "<<  (*it_JEMEtSums)-> Ex() 
	   <<"   Ey: "<<  (*it_JEMEtSums)-> Ey() 
	   <<"   Et: "<<  (*it_JEMEtSums)-> Et() <<endreq;
    }   


  //******************************* 
  //********* JEM RoI ************* 
  //******************************* 
  mLog<<MSG::DEBUG<<"in fill Histograms()"<<endreq;

  const JemRoiCollection* JEMRoIs = 0;
  sc = m_storeGate->retrieve (JEMRoIs, m_JEMRoILocation);
    	 
  if (sc==StatusCode::FAILURE)
    {
      mLog <<MSG::ERROR<<"No JEM RoIs found in TES at"<< m_JEMRoILocation<<endreq;
      return StatusCode::SUCCESS;    
    }

  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" JEM RoIs ---------------"<<endreq;


  // Step over all cells and put into hist
  JemRoiCollection::const_iterator it_JEMRoIs ;

  for( it_JEMRoIs  = JEMRoIs ->begin(); it_JEMRoIs < JEMRoIs -> end(); ++it_JEMRoIs )
    {	  
      mLog<<"JEMRoI Word: "<<(*it_JEMRoIs)->roiWord()<<endreq;
      mLog<<"Crate: "<<(*it_JEMRoIs)->crate()<<"; JEM: "<<(*it_JEMRoIs)->jem()<<"; Frame: "<<(*it_JEMRoIs)->frame()<<endreq;
      mLog<<"Location: "<<(*it_JEMRoIs)->location()<<"; forward: "<<(*it_JEMRoIs)->forward()<<"; Hits: "<<Help->Binary((*it_JEMRoIs)->hits(),8)<<endreq;

      int j=0;
      int JetMult=0;
      std::string JEMRoI = Help->Binary((*it_JEMRoIs)-> hits(),8);

      for(j=0;j<8;j++)
	{
	  JetMult=Help->Multiplicity(JEMRoI,j,1);
	  
	  mLog<<MSG::VERBOSE<<"Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	  m_h_JEMRoI_Hits -> Fill( j, JetMult);
	}
    }

      //const LVL1::JEMRoI* const roi = *pos;
      //const uint32_t key = roi->roiWord();




  	

   return StatusCode( StatusCode::SUCCESS );
}

/*---------------------------------------------------------*/
StatusCode JEMMon::procHistograms( bool isEndOfEventsBlock, 
				   bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
        if( isEndOfEventsBlock || isEndOfLumiBlock ) 
	  {

	}
	
	if( isEndOfRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}

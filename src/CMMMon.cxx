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
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "TrigT1Calo/EnergyTrigger.h"
#include "TrigT1Calo/LVL1TriggerMenuDefs.h"
#include "TrigT1Calo/LVL1TriggerMenu.h"
#include "TrigT1Calo/InternalJetROI.h"
#include "TrigT1Calo/CMMRoI.h"


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
  MonGroup CMM ( this, m_PathInRootFile, shift, eventsBlock );
  HistoBooker* shift_Booker = new HistoBooker(&CMM, &mLog, m_DataType);
  
  
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }
  
  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
      // CMM Jet Hits
      m_h_CMMJetHits_MainJets = shift_Booker->book1F("TotalMainJets", "CMM  --  Total Main Jets", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMJetHits_FwdJetsRight = shift_Booker->book1F("TotalFwdJets_Right", "CMM  --  Total Forward Jets Right",4 , -0.5, 3.5, "Threshold No.", "#");
      m_h_CMMJetHits_FwdJetsLeft = shift_Booker->book1F("TotalFwdJets_Left", "CMM  --  Total Forward Jets Left",4 , -0.5, 3.5,  "Threshold No.", "#");
      
      m_h_CMMJetHits_EtMap = shift_Booker->book1F("EtMap", "CMM  --  Et Map", 4 ,-0.5, 3.5, "Threshold No.", "#");
            
      // CMM Et Sums
      m_h_CMMEtSums_Ex = shift_Booker->book1F("Ex", "CMM  --  Ex", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_Ey = shift_Booker->book1F("Ey", "CMM  --  Ey", 250, 0,250, "Ex [GeV]", "#");
      m_h_CMMEtSums_Et = shift_Booker->book1F("Et", "CMM  --  Et", 250, 0,250, "Ex [GeV]", "#");
      
      m_h_CMMEtSums_MissingEtMap = shift_Booker->book1F("MissingEtMap", "CMM  --  Missing Et Map", 8, -0.5, 7.5, "Threshold No.", "#");
      m_h_CMMEtSums_SumEtMap = shift_Booker->book1F("SumEtMap", " CMM  --  Sum Et Map", 4, -0.5, 3.5, "Threshold No.", "#");
      
      // CMM RoI
      m_h_CMMRoI_JetEtHits =shift_Booker->book1F("JetEtHits"," CMM  --  RoI JetEtHits",4, -0.5,3.5,"Threshold No.","#");
      
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
      mLog << MSG::DEBUG
	   << "No CMM JetHits found in TES at "
	   << m_CMMJetHitsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }
  
  mLog<<MSG::DEBUG<<endreq;
  mLog<<MSG::DEBUG<<"-------------- CMM Jet Hits ---------------"<<endreq;
  
  // Step over all cells and put into hist
  CMMJetHitsCollection::const_iterator it_CMMJetHits ;
  
  for( it_CMMJetHits  = CMMJetHits ->begin(); it_CMMJetHits < CMMJetHits -> end(); ++it_CMMJetHits )
    {	  
      std::string CMMHit = Help->Binary((*it_CMMJetHits)->Hits(),24);
      int j=0;
      int JetMult=0;
      // the binary hit information is represented by an integer number, that is converted to a string in order
      // to get the real binary information
      // later the multiplicities of the several thresholds are retrieved from this string
      

      mLog<<MSG::DEBUG<<"CMMJetHits dataID: "<< (*it_CMMJetHits)-> dataID() <<"   Hits: "<< (*it_CMMJetHits)-> Hits()
	  << " Hits(binary): " << CMMHit <<endreq;

      if ((*it_CMMJetHits)-> dataID() == 18)  //main total jets
	{
	  for(int j=0;j<8;j++)
	    {
	      JetMult=Help->Multiplicity(CMMHit,j,3);

	      mLog<<MSG::DEBUG<<"Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_CMMJetHits_MainJets -> Fill( j ,JetMult);
	    }
	}

      if ((*it_CMMJetHits)-> dataID() == 21)  //fwd total jets
	{
	  CMMHit= Help->Binary((*it_CMMJetHits)->Hits(),16); //total fwd jets only 16 bit long!
	  mLog<<MSG::DEBUG<<"FwdCMMHit: "<<CMMHit<<endreq;
	  
	  for(j=0;j<4;j++) //fwd Right jets
	    {
	      JetMult=Help->Multiplicity(CMMHit,j,2);

	      m_h_CMMJetHits_FwdJetsRight -> Fill( j ,JetMult);
	      mLog<<MSG::DEBUG<<"Right FwdThresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	    }
	  for(j=4;j<8;j++) //fwd left jets
	    {
	      JetMult=Help->Multiplicity(CMMHit,j,2);

	      m_h_CMMJetHits_FwdJetsLeft-> Fill( (j-4) ,JetMult);
	      mLog<<MSG::DEBUG<<"Left fwdThresh.No: "<<(j-4)<<" Multiplicity: "<< JetMult<<endreq;
	    }
	}

      if ((*it_CMMJetHits)-> dataID() == 22)  //JetEtSum Hitmap
	{
	  CMMHit= Help->Binary((*it_CMMJetHits)->Hits(),4);

	  for(int j=0;j<4;j++)
	    {
	      JetMult=Help->Multiplicity(CMMHit,j,1);

	      mLog<<MSG::DEBUG<<"JetEtSum Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_CMMJetHits_EtMap -> Fill( j ,JetMult);
	    }
	}
       
    }     	
//******************************* 
//********* CMM Et Sums ********* 
//******************************* 

  //fill  Histos
  const CMMEtSumsCollection* CMMEtSums;
  sc = m_storeGate->retrieve(CMMEtSums, m_CMMEtSumsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No CMMEtSums found in TES at "
	   << m_CMMEtSumsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<"-------------- CMM Et Sums ---------------"<<endreq;

  // Step over all cells and put into hist
  CMMEtSumsCollection::const_iterator it_CMMEtSums ;

  for( it_CMMEtSums  = CMMEtSums ->begin(); it_CMMEtSums < CMMEtSums -> end(); ++it_CMMEtSums )
    {	  
     // fill histograms
      if  ((*it_CMMEtSums)-> dataID()==18) // total energy sums
	{
	  m_h_CMMEtSums_Ex -> Fill( (*it_CMMEtSums)-> Ex(), 1.);
	  m_h_CMMEtSums_Ey -> Fill( (*it_CMMEtSums)-> Ey(), 1.);
	  m_h_CMMEtSums_Et -> Fill( (*it_CMMEtSums)-> Et(), 1.);
	}
      
      if ((*it_CMMEtSums)-> dataID() == 19)  //SumEt Hitmap
	{
	  std::string CMMHit= Help->Binary((*it_CMMEtSums)->Et(),8);

	  for(int j=0;j<8;j++)
	    {
	      int JetMult=Help->Multiplicity(CMMHit,j,1);

	      mLog<<MSG::DEBUG<<"JetEtSum Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_CMMEtSums_SumEtMap -> Fill( j ,JetMult);
	    }
	}

      if ((*it_CMMEtSums)-> dataID() == 20)  //SumEt Hitmap
	{
	  std::string CMMHit= Help->Binary((*it_CMMEtSums)->Et(),4);

	  for(int j=0;j<4;j++)
	    {
	      int JetMult=Help->Multiplicity(CMMHit,j,1);

	      mLog<<MSG::DEBUG<<"JetEtSum Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_CMMEtSums_SumEtMap -> Fill( j ,JetMult);
	    }
	}

      mLog<<MSG::DEBUG<<" CMMEtSums crate: "<< (*it_CMMEtSums)-> crate()<<" dataID: "<< (*it_CMMEtSums)-> dataID()
	  <<" Ex: "<<  (*it_CMMEtSums)-> Ex()
	  <<" Ey: "<<  (*it_CMMEtSums)-> Ey()
	  <<" Et: "<<  (*it_CMMEtSums)-> Et()<<endreq;
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

  int j=0;
  int JetMult=0;
  std::string CMMRoIHit = Help->Binary((CR)-> jetEtHits(),4);

  for(j=0;j<4;j++)
    {
      JetMult=Help->Multiplicity(CMMRoIHit,j,1);
      
      mLog<<MSG::VERBOSE<<"Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
      m_h_CMMRoI_JetEtHits -> Fill( j, JetMult);
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

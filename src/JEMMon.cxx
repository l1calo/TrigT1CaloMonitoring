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
//#include "TrigT1Calo/JEMEtSums.h"
#include "TrigT1Calo/EnergyTrigger.h"
#include "TrigT1CaloMonitoring/MonHelpers.h"


#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "AthenaMonitoring/AthenaMonManager.h"


// *********************************************************************
// Public Methods
// *********************************************************************

JEMMon::
JEMMon( const std::string & type, const std::string & name,
                        const IInterface* parent )
	: ManagedMonitorToolBase( type, name, parent )
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "BS_JEMHitsLocation", m_BS_JEMHitsLocation =  LVL1::TrigT1CaloDefs::JEMHitsLocation) ;
  declareProperty( "BS_JEMEtSumsLocation", m_BS_JEMEtSumsLocation=   LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;

  declareProperty( "Sim_JEMHitsLocation", m_Sim_JEMHitsLocation =  LVL1::TrigT1CaloDefs::JEMHitsLocation) ;
  declareProperty( "Sim_JEMEtSumsLocation", m_Sim_JEMEtSumsLocation=   LVL1::TrigT1CaloDefs::JEMEtSumsLocation) ;

}


JEMMon::
~JEMMon()
{
}


//_______________________________ book Histograms ___________________________________________
StatusCode
JEMMon::
bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
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
	
  MonGroup BS_JEM ( this, "Stats/BS_JEM", expert, eventsBlock );

  MonGroup Sim_JEM ( this, "Stats/Sim_JEM", expert, eventsBlock );

  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
	  
      // JEM Hits
      //BS Histos
      m_h_BS_JEMHits_MainHits = new TH1F("MainHits", "BS JEMHits  --  Main Hits)", 8, -0.5,7.5);
      m_h_BS_JEMHits_MainHits -> GetXaxis() -> SetTitle("Threshold No.");
      BS_JEM.regHist (m_h_BS_JEMHits_MainHits);

      m_h_BS_JEMHits_FwdHitsRight = new TH1F("FwdHitsRight", "BS JEMHits  --  Forward Hits Right", 4, -0.5,3.5);
      m_h_BS_JEMHits_FwdHitsRight -> GetXaxis() -> SetTitle("Threshold No.");
      BS_JEM.regHist (m_h_BS_JEMHits_FwdHitsRight);

      m_h_BS_JEMHits_FwdHitsLeft = new TH1F("FwdHitsLeft", "BS JEMHits  --  Forward HitsLeft", 4, -0.5,3.5);
      m_h_BS_JEMHits_FwdHitsLeft -> GetXaxis() -> SetTitle("Threshold No.");
      BS_JEM.regHist (m_h_BS_JEMHits_FwdHitsLeft);

      //Sim Histos
      m_h_Sim_JEMHits_MainHits = new TH1F("MainHits", "Sim JEMHits  --  Main Hits", 8, -0.5,7.5);
      m_h_Sim_JEMHits_MainHits -> GetXaxis() -> SetTitle("Threshold No.");
      Sim_JEM.regHist (m_h_Sim_JEMHits_MainHits);

      m_h_Sim_JEMHits_FwdHitsRight = new TH1F("FwdHitsRight", "Sim JEMHits  --  Forward Hits Right", 4, -0.5,3.5);
      m_h_Sim_JEMHits_FwdHitsRight -> GetXaxis() -> SetTitle("Threshold No.");
      Sim_JEM.regHist (m_h_Sim_JEMHits_FwdHitsRight);

      m_h_Sim_JEMHits_FwdHitsLeft = new TH1F("FwdHitsLeft", "Sim JEMHits  --  Forward Hits Left", 4, -0.5,3.5);
      m_h_Sim_JEMHits_FwdHitsLeft -> GetXaxis() -> SetTitle("Threshold No.");
      Sim_JEM.regHist (m_h_Sim_JEMHits_FwdHitsLeft);

      // JEM Et Sums
      //BS Histos
      m_h_BS_JEMEtSums_Ex = new TH1F("JEM Ex", "BS JEMEtSums  --  Ex", 250, 0,250);
      m_h_BS_JEMEtSums_Ex -> GetXaxis() -> SetTitle("Ex [GeV]");
      BS_JEM.regHist (m_h_BS_JEMEtSums_Ex);

      m_h_BS_JEMEtSums_Ey = new TH1F("JEM Ey", "BS JEMEtSums  --  Ey", 250, 0,250);
      m_h_BS_JEMEtSums_Ey -> GetXaxis() -> SetTitle("Ey [GeV]");
      BS_JEM.regHist (m_h_BS_JEMEtSums_Ey);

      m_h_BS_JEMEtSums_Et = new TH1F("JEM Et", "BS JEMEtSums  --  Et", 250, 0,250);
      m_h_BS_JEMEtSums_Et -> GetXaxis() -> SetTitle("Et [GeV]");
      BS_JEM.regHist (m_h_BS_JEMEtSums_Et);

      //Sim Histos
      m_h_Sim_JEMEtSums_Ex = new TH1F("JEM Ex", "Sim JEMEtSums  --  Ex", 250, 0,250);
      m_h_Sim_JEMEtSums_Ex -> GetXaxis() -> SetTitle("Ex [GeV]");
      Sim_JEM.regHist (m_h_Sim_JEMEtSums_Ex);

      m_h_Sim_JEMEtSums_Ey = new TH1F("JEM Ey", "Sim JEMEtSums  --  Ey", 250, 0,250);
      m_h_Sim_JEMEtSums_Ey -> GetXaxis() -> SetTitle("Ey [GeV]");
      Sim_JEM.regHist (m_h_Sim_JEMEtSums_Ey);

      m_h_Sim_JEMEtSums_Et = new TH1F("JEM Et", "Sim JEMEtSums  --  Et", 250, 0,250);
      m_h_Sim_JEMEtSums_Et -> GetXaxis() -> SetTitle("Et [GeV]");
      Sim_JEM.regHist (m_h_Sim_JEMEtSums_Et);

    }
  
  if( isNewRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}



//_______________________________ fill Histograms ___________________________________________
StatusCode
JEMMon::
fillHistograms()
{
  MsgStream mLog( msgSvc(), name() );

  //******************************* 
  //********* JEM Hits  *********** 
  //******************************* 

  // fill BS Histos
  const JEMHitsCollection* JEMHits;
  StatusCode sc = m_storeGate->retrieve(JEMHits, m_BS_JEMHitsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No JEMHits found in TES at "
	   << m_BS_JEMHitsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<endreq;
  mLog<<MSG::DEBUG<<"-------------- JEM Hits ---------------"<<endreq;

  // Step over all cells and put into hist
  JEMHitsCollection::const_iterator it_JEMHits ;

  for( it_JEMHits  = JEMHits ->begin(); it_JEMHits < JEMHits -> end(); ++it_JEMHits )
    {	  
      int j=0;
      int JetMult=0;
      std::string JEMHit = Binary((*it_JEMHits)-> JetHits(),24);
      // the binary hit information is represented by an integer number, that is converted to a string in order
      // to get the real binary information
      // later the multiplicities of the several thresholds are retrieved from this string

      mLog<<MSG::DEBUG<<"BS JEMHits Crate: "<< (*it_JEMHits)->crate()<<"  Module: "<<(*it_JEMHits)->module()
	  << "  JetHits: "<<(*it_JEMHits)-> JetHits()
	  << "   Hits(binary): "  <<JEMHit<<   endreq;
      
      if  ((*it_JEMHits)->forward()==0) //Main Jets
	{
	  for(j=0;j<8;j++)
	    {
	      JetMult=Multiplicity(JEMHit,j,3);

	      mLog<<MSG::DEBUG<<"BS MainThresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_BS_JEMHits_MainHits -> Fill( j, JetMult);
	    }
	}

      if  ((*it_JEMHits)->forward()==1) //fwd jets a bit complicated!
	// fwd and main hits are contained in the same hitword
	{
	  for(j=0;j<4;j++)  //fwd hits
	    {
	      JetMult=Multiplicity(JEMHit,j,2); // only 2 bits per thresh

	      if (((*it_JEMHits)-> module()==0) or((*it_JEMHits)-> module()==8) )//left fwd hits
		// JEMs No 0 and 8 are processing forward left hits,
		// JEMs No 7 and 15 forward right hits
		{
		  mLog<<MSG::DEBUG<<"BS LeftFwd Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
		  m_h_BS_JEMHits_FwdHitsLeft -> Fill( j, JetMult);
		}
	      if (((*it_JEMHits)-> module()==7) or((*it_JEMHits)-> module()==15) )//right fwd hits
		{
		  mLog<<MSG::DEBUG<<"BS RightFwd Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
		  m_h_BS_JEMHits_FwdHitsRight -> Fill( j, JetMult);
		}


	    }
	  for(j=4;j<12;j++) //main hits
	    {
	      JetMult=Multiplicity(JEMHit,j,2);

	      mLog<<MSG::DEBUG<<"BS Main Thresh.No: "<<(j-4)<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_BS_JEMHits_MainHits -> Fill( (j-4), JetMult);
	    }

	}


    }     	
 

  mLog<<MSG::DEBUG<<"Sim"<<endreq;


  //fill Sim Histos
  sc = m_storeGate->retrieve(JEMHits, m_Sim_JEMHitsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No JEMHits found in TES at "
	   << m_Sim_JEMHitsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  // Step over all cells and put into hist
  for( it_JEMHits  = JEMHits ->begin(); it_JEMHits < JEMHits -> end(); ++it_JEMHits )
    {	  
     // fill histograms
      int j=0;
      int JetMult=0;
      std::string JEMHit = Binary((*it_JEMHits)-> JetHits(),24);

      mLog<<MSG::DEBUG<<"Sim JEMHits Crate: "<< (*it_JEMHits)->crate()<<"  Module: "<<(*it_JEMHits)->module()
	  << "  JetHits: "<<(*it_JEMHits)-> JetHits()
	  << " Hits(binary): " <<JEMHit <<   endreq;

      if  ((*it_JEMHits)->forward()==0)   //main jets
	{
	  for(j=0;j<8;j++)
	    {
	      JetMult=Multiplicity(JEMHit,j,3);

	      mLog<<MSG::DEBUG<<"Sim Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_Sim_JEMHits_MainHits -> Fill( j, JetMult);
	    }
	}

      if  ((*it_JEMHits)->forward()==1) //fwd jets
	{
	  for(j=0;j<4;j++)
	    {

	      JetMult=Multiplicity(JEMHit,j,2);

	      if (((*it_JEMHits)-> module()==0) or((*it_JEMHits)-> module()==8) )//left fwd hits
		{
		  mLog<<MSG::DEBUG<<"Sim LeftFwd Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
		  m_h_Sim_JEMHits_FwdHitsLeft -> Fill( j, JetMult);
		}
	      if (((*it_JEMHits)-> module()==7) or((*it_JEMHits)-> module()==15) )//right fwd hits
		{
		  mLog<<MSG::DEBUG<<"Sim RightFwd Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
		  m_h_Sim_JEMHits_FwdHitsRight -> Fill( j, JetMult);
		}

	    }
	  for(j=4;j<12;j++)
	    {
	      JetMult=Multiplicity(JEMHit,j,2);

	      mLog<<MSG::DEBUG<<"Sim Main Thresh.No: "<<(j-4)<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_Sim_JEMHits_MainHits -> Fill( (j-4), JetMult);
	    }

	}


    } 

  //******************************* 
  //********* JEM Et Sums ********* 
  //******************************* 

  //fill BS Histos
  const JEMEtSumsCollection* JEMEtSums;
  sc = m_storeGate->retrieve(JEMEtSums, m_BS_JEMEtSumsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No JEMEtSums found in TES at "
	   << m_BS_JEMEtSumsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<"-------------- JEM Et Sums ---------------"<<endreq;

  // Step over all cells and put into hist
  JEMEtSumsCollection::const_iterator it_JEMEtSums ;

  for( it_JEMEtSums  = JEMEtSums ->begin(); it_JEMEtSums < JEMEtSums -> end(); ++it_JEMEtSums )
    {	  
     
     // fill histograms
      m_h_BS_JEMEtSums_Ex -> Fill( (*it_JEMEtSums)-> Ex(), 1.); 
      m_h_BS_JEMEtSums_Ey -> Fill( (*it_JEMEtSums)-> Ey(), 1.); 
      m_h_BS_JEMEtSums_Et -> Fill( (*it_JEMEtSums)-> Et() , 1.); 
      mLog <<MSG::DEBUG<< "BS JEMEtSums Crate: "<<(*it_JEMEtSums)->crate()<<"  Module: "<<(*it_JEMEtSums)->module()
	   <<"   Ex: "<<  (*it_JEMEtSums)-> Ex() 
	   <<"   Ey: "<<  (*it_JEMEtSums)-> Ey() 
	   <<"   Et: "<<  (*it_JEMEtSums)-> Et() <<endreq;
    }     	


  //fill Sim Histos
  sc = m_storeGate->retrieve(JEMEtSums, m_Sim_JEMEtSumsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No JEMEtSums found in TES at "
	   << m_Sim_JEMEtSumsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  // Step over all cells and put into hist

  for( it_JEMEtSums  = JEMEtSums ->begin(); it_JEMEtSums < JEMEtSums -> end(); ++it_JEMEtSums )
    {	  
     // fill histograms
      m_h_Sim_JEMEtSums_Ex -> Fill( (*it_JEMEtSums)-> Ex(), 1.); 
      m_h_Sim_JEMEtSums_Ey -> Fill( (*it_JEMEtSums)-> Ey(), 1.); 
      m_h_Sim_JEMEtSums_Et -> Fill( (*it_JEMEtSums)-> Et(), 1.); 
      mLog <<MSG::DEBUG<< "Sim JEMEtSums Crate: "<<(*it_JEMEtSums)->crate()<<"  Module: "<<(*it_JEMEtSums)->module()
	   <<"   Ex: "<<  (*it_JEMEtSums)-> Ex() 
	   <<"   Ey: "<<  (*it_JEMEtSums)-> Ey() 
	   <<"   Et: "<<  (*it_JEMEtSums)-> Et() <<endreq;
    }     	 

   return StatusCode( StatusCode::SUCCESS );
}

//_______________________________ proc  Histograms ___________________________________________
StatusCode
JEMMon::
procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
{
        if( isEndOfEventsBlock || isEndOfLumiBlock ) 
	  {

	}
	
	if( isEndOfRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}

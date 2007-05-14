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
#include "TrigT1Calo/EnergyTrigger.h"
#include "TrigT1Calo/LVL1TriggerMenuDefs.h"
#include "TrigT1Calo/LVL1TriggerMenu.h"
#include "TrigT1Calo/InternalJetROI.h"
#include "TrigT1CaloMonitoring/MonHelpers.h"


#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "AthenaMonitoring/AthenaMonManager.h"


// *********************************************************************
// Public Methods
// *********************************************************************

CMMMon::
CMMMon( const std::string & type, const std::string & name,
                        const IInterface* parent )
	: ManagedMonitorToolBase( type, name, parent )
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "BS_CMMJetHitsLocation", m_BS_CMMJetHitsLocation =  LVL1::TrigT1CaloDefs::CMMJetHitsLocation) ;
  declareProperty( "BS_CMMEtSumsLocation", m_BS_CMMEtSumsLocation =  LVL1::TrigT1CaloDefs::CMMEtSumsLocation) ;

  declareProperty( "Sim_CMMJetHitsLocation", m_Sim_CMMJetHitsLocation =  LVL1::TrigT1CaloDefs::CMMJetHitsLocation) ;
  declareProperty( "Sim_CMMEtSumsLocation", m_Sim_CMMEtSumsLocation =  LVL1::TrigT1CaloDefs::CMMEtSumsLocation) ;
  

}


CMMMon::
~CMMMon()
{
}


//_______________________________ book Histograms ___________________________________________
StatusCode
CMMMon::
bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
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
 
 if( m_environment == AthenaMonManager::online ) {
   // book histograms that are only made in the online environment...
 }
	
 if( m_dataType == AthenaMonManager::cosmics ) {
   // book histograms that are only relevant for cosmics data...
 }
	
 MonGroup BS_CMM ( this, "Stats/BS_CMM", shift, eventsBlock );

 MonGroup Sim_CMM ( this, "Stats/Sim_CMM", shift, eventsBlock );

 if( isNewEventsBlock || isNewLumiBlock ) 
   {	
     
     // CMM Jet Hits
     //BS Histos
     m_h_BS_CMMJetHits_MainJets = new TH1F("TotalMainJets", "BS CMMJetHits  --  Total Main Jets", 8, -0.5,7.5);
     m_h_BS_CMMJetHits_MainJets -> GetXaxis() -> SetTitle("Threshold No.");
     BS_CMM.regHist (m_h_BS_CMMJetHits_MainJets);

     m_h_BS_CMMJetHits_FwdJetsRight = new TH1F("TotalFwdJets Right", "BS CMMJetHits  --  Total Forward Jets Right",4 , -0.5,3.5);
     m_h_BS_CMMJetHits_FwdJetsRight -> GetXaxis() -> SetTitle("Threshold No.");
     BS_CMM.regHist (m_h_BS_CMMJetHits_FwdJetsRight);

     m_h_BS_CMMJetHits_FwdJetsLeft = new TH1F("TotalFwdJets Left", "BS CMMJetHits  --  Total Forward Jets Left",4 , -0.5,3.5);
     m_h_BS_CMMJetHits_FwdJetsLeft -> GetXaxis() -> SetTitle("Threshold No.");
     BS_CMM.regHist (m_h_BS_CMMJetHits_FwdJetsLeft);

     m_h_BS_CMMJetHits_EtMap = new TH1F("EtMap", "BS CMMJetHits  --  Et Map",4 ,-0.5,3.5);
     m_h_BS_CMMJetHits_EtMap -> GetXaxis() -> SetTitle("Threshold No.");
     BS_CMM.regHist (m_h_BS_CMMJetHits_EtMap);

     //Sim Histos
     m_h_Sim_CMMJetHits_MainJets = new TH1F("TotalMainJets", "Sim CMMJetHits  --  Total MainJets", 8, -0.5,7.5);
     m_h_Sim_CMMJetHits_MainJets -> GetXaxis() -> SetTitle("Threshold No.");
     Sim_CMM.regHist (m_h_Sim_CMMJetHits_MainJets);

     m_h_Sim_CMMJetHits_FwdJetsRight = new TH1F("TotalFwdJets Right", "Sim CMMJetHits  --  Total Forward Jets Right",4 , -0.5,3.5);
     m_h_Sim_CMMJetHits_FwdJetsRight -> GetXaxis() -> SetTitle("Threshold No.");
     Sim_CMM.regHist (m_h_Sim_CMMJetHits_FwdJetsRight);

     m_h_Sim_CMMJetHits_FwdJetsLeft = new TH1F("TotalFwdJets Left", "Sim CMMJetHits  --  Total Forward Jets Left",4 , -0.5,3.5);
     m_h_Sim_CMMJetHits_FwdJetsLeft -> GetXaxis() -> SetTitle("Threshold No.");
     Sim_CMM.regHist (m_h_Sim_CMMJetHits_FwdJetsLeft);

     m_h_Sim_CMMJetHits_EtMap = new TH1F("EtMap", "Sim CMMJetHits  --  Et Map",4 , -0.5,3.5);
     m_h_Sim_CMMJetHits_EtMap -> GetXaxis() -> SetTitle("Threshold No.");
     Sim_CMM.regHist (m_h_Sim_CMMJetHits_EtMap);

     // CMM Et Sums
     //BS Histos
     m_h_BS_CMMEtSums_Ex = new TH1F("CMM Ex", "BS CMMEtSums  --  Ex", 250, 0,250);
     m_h_BS_CMMEtSums_Ex -> GetXaxis() -> SetTitle("Ex [GeV]");
     BS_CMM.regHist (m_h_BS_CMMEtSums_Ex);

     m_h_BS_CMMEtSums_Ey = new TH1F("CMM Ey", "BS CMMEtSums  --  Ey", 250, 0,250);
     m_h_BS_CMMEtSums_Ey -> GetXaxis() -> SetTitle("Ey [GeV]");
     BS_CMM.regHist (m_h_BS_CMMEtSums_Ey);

     m_h_BS_CMMEtSums_Et = new TH1F("CMM Et", "BS CMMEtSums  --  Et", 250, 0,250);
     m_h_BS_CMMEtSums_Et -> GetXaxis() -> SetTitle("Et [GeV]");
     BS_CMM.regHist (m_h_BS_CMMEtSums_Et);

     m_h_BS_CMMEtSums_MissingEtMap = new TH1F("CMM MissingEt Map", "BS CMMEtSums  --  Missing Et Map", 8, -0.5,7.5);
     m_h_BS_CMMEtSums_MissingEtMap -> GetXaxis() -> SetTitle("Threshold No.");
     BS_CMM.regHist (m_h_BS_CMMEtSums_MissingEtMap);

     m_h_BS_CMMEtSums_SumEtMap = new TH1F("CMM SumEt Map", "BS CMMEtSums  --  Sum Et Map", 4, -0.5,3.5);
     m_h_BS_CMMEtSums_SumEtMap -> GetXaxis() -> SetTitle("Threshold No.");
     BS_CMM.regHist (m_h_BS_CMMEtSums_SumEtMap);

     //Sim Histos
     m_h_Sim_CMMEtSums_Ex = new TH1F("CMM Ex", "Sim CMMEtSums  --  Ex", 250, 0,250);
     m_h_Sim_CMMEtSums_Ex -> GetXaxis() -> SetTitle("Ex [GeV]");
     Sim_CMM.regHist (m_h_Sim_CMMEtSums_Ex);

     m_h_Sim_CMMEtSums_Ey = new TH1F("CMM Ey", "Sim CMMEtSums  --  Ey", 250, 0,250);
     m_h_Sim_CMMEtSums_Ey -> GetXaxis() -> SetTitle("Ey [GeV]");
     Sim_CMM.regHist (m_h_Sim_CMMEtSums_Ey);

     m_h_Sim_CMMEtSums_Et = new TH1F("CMM Et", "Sim CMMEtSums  --  Et", 250, 0,250);
     m_h_Sim_CMMEtSums_Et -> GetXaxis() -> SetTitle("Et [GeV]");
     Sim_CMM.regHist (m_h_Sim_CMMEtSums_Et);

     m_h_Sim_CMMEtSums_MissingEtMap = new TH1F("CMM MissingEt Map", "Sim CMMEtSums  --  Missing Et Map", 8, -0.5,7.5);
     m_h_Sim_CMMEtSums_MissingEtMap -> GetXaxis() -> SetTitle("Threshold No.");
     Sim_CMM.regHist (m_h_Sim_CMMEtSums_MissingEtMap);

     m_h_Sim_CMMEtSums_SumEtMap = new TH1F("CMM SumEt Map", "Sim CMMEtSums  --  Sum Et Map", 4, -0.5,3.5);
     m_h_Sim_CMMEtSums_SumEtMap -> GetXaxis() -> SetTitle("Threshold No.");
     Sim_CMM.regHist (m_h_Sim_CMMEtSums_SumEtMap);

   }
	
 if( isNewRun ) { }
  
 return StatusCode( StatusCode::SUCCESS );
}


//_______________________________ fill Histograms ___________________________________________
StatusCode
CMMMon::
fillHistograms()
{
  MsgStream mLog( msgSvc(), name() );


//******************************* 
//********* CMM Jet Hits ******** 
//******************************* 

  //fill BS Histos
  const CMMJetHitsCollection* CMMJetHits;
  StatusCode sc = m_storeGate->retrieve(CMMJetHits, m_BS_CMMJetHitsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No CMM JetHits found in TES at "
	   << m_BS_CMMJetHitsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<endreq;
  mLog<<MSG::DEBUG<<"-------------- CMM Jet Hits ---------------"<<endreq;

  // Step over all cells and put into hist
  CMMJetHitsCollection::const_iterator it_CMMJetHits ;

  for( it_CMMJetHits  = CMMJetHits ->begin(); it_CMMJetHits < CMMJetHits -> end(); ++it_CMMJetHits )
    {	  
      std::string CMMHit = Binary((*it_CMMJetHits)->Hits(),24);
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
	      JetMult=Multiplicity(CMMHit,j,3);

	      mLog<<MSG::DEBUG<<"Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_BS_CMMJetHits_MainJets -> Fill( j ,JetMult);
	    }
	}

      if ((*it_CMMJetHits)-> dataID() == 21)  //fwd total jets
	{
	  CMMHit= Binary((*it_CMMJetHits)->Hits(),16); //total fwd jets only 16 bit long!
	  mLog<<MSG::DEBUG<<"FwdCMMHit: "<<CMMHit<<endreq;
	  
	  for(j=0;j<4;j++) //fwd Right jets
	    {
	      JetMult=Multiplicity(CMMHit,j,2);

	      m_h_BS_CMMJetHits_FwdJetsRight -> Fill( j ,JetMult);
	      mLog<<MSG::DEBUG<<"Right FwdThresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	    }
	  for(j=4;j<8;j++) //fwd left jets
	    {
	      JetMult=Multiplicity(CMMHit,j,2);

	      m_h_BS_CMMJetHits_FwdJetsLeft-> Fill( (j-4) ,JetMult);
	      mLog<<MSG::DEBUG<<"Left fwdThresh.No: "<<(j-4)<<" Multiplicity: "<< JetMult<<endreq;
	    }
	}

      if ((*it_CMMJetHits)-> dataID() == 22)  //JetEtSum Hitmap
	{
	  CMMHit= Binary((*it_CMMJetHits)->Hits(),4);

	  for(int j=0;j<4;j++)
	    {
	      JetMult=Multiplicity(CMMHit,j,1);

	      mLog<<MSG::DEBUG<<"JetEtSum Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_BS_CMMJetHits_EtMap -> Fill( j ,JetMult);
	    }
	}
       
    }     	
  
  mLog<<MSG::DEBUG<<"Sim"<<endreq;

  //fill Sim Histos
  sc = m_storeGate->retrieve(CMMJetHits, m_BS_CMMJetHitsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No CMM JetHits found in TES at "
	   << m_Sim_CMMJetHitsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  // Step over all cells and put into hist
  for( it_CMMJetHits  = CMMJetHits ->begin(); it_CMMJetHits < CMMJetHits -> end(); ++it_CMMJetHits )
    {	  
     // fill histograms
      std::string CMMHit = Binary((*it_CMMJetHits)->Hits(),24);
      int j=0;
      int JetMult=0;

      mLog<<MSG::DEBUG<<"Sim CMMJetHits dataID: "<< (*it_CMMJetHits)-> dataID() <<"   Hits: "<< (*it_CMMJetHits)-> Hits()
	  << " Hits(binary): " << CMMHit <<endreq;

     if ((*it_CMMJetHits)-> dataID() == 18)  //main total jets
       {
	 for(int j=0;j<8;j++)
	   {
	     JetMult=Multiplicity(CMMHit,j,3);

	     mLog<<MSG::DEBUG<<"Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	     m_h_Sim_CMMJetHits_MainJets -> Fill( j ,JetMult);
	   }
       }

      if ((*it_CMMJetHits)-> dataID() == 21)  //fwd total jets
	{
	  CMMHit= Binary((*it_CMMJetHits)->Hits(),16); //total fwd jets only 16 bit long!
	  mLog<<MSG::DEBUG<<"FwdCMMHit: "<<CMMHit<<endreq;
	  
	  for(j=0;j<4;j++) //fwd Right jets
	    {
	      JetMult=Multiplicity(CMMHit,j,2);

	      m_h_Sim_CMMJetHits_FwdJetsRight -> Fill( j ,JetMult);
	      mLog<<MSG::DEBUG<<"Right FwdThresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	    }
	  for(j=4;j<8;j++) //fwd left jets
	    {
	      JetMult=Multiplicity(CMMHit,j,2);

	      m_h_Sim_CMMJetHits_FwdJetsLeft-> Fill( (j-4) ,JetMult);
	      mLog<<MSG::DEBUG<<"Left fwdThresh.No: "<<(j-4)<<" Multiplicity: "<< JetMult<<endreq;
	    }
	}

      if ((*it_CMMJetHits)-> dataID() == 22)  //JetEtSum Hitmap
	{
	  CMMHit= Binary((*it_CMMJetHits)->Hits(),4);

	  for(int j=0;j<4;j++)
	    {
	      JetMult=Multiplicity(CMMHit,j,1);

	      mLog<<MSG::DEBUG<<"JetEtSum Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_Sim_CMMJetHits_EtMap -> Fill( j ,JetMult);
	    }
	}

    }
//******************************* 
//********* CMM Et Sums ********* 
//******************************* 

  //fill BS Histos
  const CMMEtSumsCollection* CMMEtSums;
  sc = m_storeGate->retrieve(CMMEtSums, m_BS_CMMEtSumsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No CMMEtSums found in TES at "
	   << m_BS_CMMEtSumsLocation
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
	  m_h_BS_CMMEtSums_Ex -> Fill( (*it_CMMEtSums)-> Ex(), 1.);
	  m_h_BS_CMMEtSums_Ey -> Fill( (*it_CMMEtSums)-> Ey(), 1.);
	  m_h_BS_CMMEtSums_Et -> Fill( (*it_CMMEtSums)-> Et(), 1.);
	}
      
      if ((*it_CMMEtSums)-> dataID() == 19)  //SumEt Hitmap
	{
	  std::string CMMHit= Binary((*it_CMMEtSums)->Et(),8);

	  for(int j=0;j<8;j++)
	    {
	      int JetMult=Multiplicity(CMMHit,j,1);

	      mLog<<MSG::DEBUG<<"JetEtSum Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_BS_CMMEtSums_SumEtMap -> Fill( j ,JetMult);
	    }
	}

      if ((*it_CMMEtSums)-> dataID() == 20)  //SumEt Hitmap
	{
	  std::string CMMHit= Binary((*it_CMMEtSums)->Et(),4);

	  for(int j=0;j<4;j++)
	    {
	      int JetMult=Multiplicity(CMMHit,j,1);

	      mLog<<MSG::DEBUG<<"JetEtSum Thresh.No: "<<j<<" Multiplicity: "<< JetMult<<endreq;
	      m_h_BS_CMMEtSums_SumEtMap -> Fill( j ,JetMult);
	    }
	}

      mLog<<MSG::DEBUG<<"BS CMMEtSums crate: "<< (*it_CMMEtSums)-> crate()<<" dataID: "<< (*it_CMMEtSums)-> dataID()
	  <<" Ex: "<<  (*it_CMMEtSums)-> Ex()
	  <<" Ey: "<<  (*it_CMMEtSums)-> Ey()
	  <<" Et: "<<  (*it_CMMEtSums)-> Et()<<endreq;
    }     	


  //fill Sim Histos
  sc = m_storeGate->retrieve(CMMEtSums, m_Sim_CMMEtSumsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No CMMEtSums found in TES at "
	   << m_Sim_CMMEtSumsLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  // Step over all cells and put into hist
  for( it_CMMEtSums  = CMMEtSums ->begin(); it_CMMEtSums < CMMEtSums -> end(); ++it_CMMEtSums )
    {	  
     // fill histograms
      if  ((*it_CMMEtSums)-> dataID()==18) //only total energy sums
	{
	  m_h_Sim_CMMEtSums_Ex -> Fill( (*it_CMMEtSums)-> Ex(), 1.);
	  m_h_Sim_CMMEtSums_Ey -> Fill( (*it_CMMEtSums)-> Ey(), 1.);
	  m_h_Sim_CMMEtSums_Et -> Fill( (*it_CMMEtSums)-> Et(), 1.);
	}
      
      if  ((*it_CMMEtSums)-> dataID()==19) m_h_Sim_CMMEtSums_MissingEtMap -> Fill( (*it_CMMEtSums)-> Et(), 1.);
      if  ((*it_CMMEtSums)-> dataID()==20) m_h_Sim_CMMEtSums_SumEtMap -> Fill( (*it_CMMEtSums)-> Et(), 1.);

      mLog<<MSG::DEBUG<<"Sim CMMEtSums crate: "<< (*it_CMMEtSums)-> crate()<<" dataID: "<< (*it_CMMEtSums)-> dataID()
	  <<" Ex: "<<  (*it_CMMEtSums)-> Ex()
	  <<" Ey: "<<  (*it_CMMEtSums)-> Ey()
	  <<" Et: "<<  (*it_CMMEtSums)-> Et()<<endreq;

    } 

  return StatusCode( StatusCode::SUCCESS );
}

//_______________________________ proc  Histograms ___________________________________________
StatusCode
CMMMon::
procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
{
        if( isEndOfEventsBlock || isEndOfLumiBlock ) 
	  {

	}
	
	if( isEndOfRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}

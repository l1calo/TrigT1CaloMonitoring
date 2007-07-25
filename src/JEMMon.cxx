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
#include "TrigT1Calo/QuadLinear.h"

#include "TrigT1Interfaces/JEPRoIDecoder.h"
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
	
  ManagedMonitorToolBase::LevelOfDetail_t LevelOfDetail=shift;
  if (m_DataType=="Sim") LevelOfDetail = expert;

  MonGroup JEM_DAQ ( this, (m_PathInRootFile+"/DAQ").c_str(), expert, eventsBlock );
  HistoBooker* DAQ_Booker = new HistoBooker(&JEM_DAQ, &mLog, m_DataType);

  MonGroup JEM_RoI ( this, (m_PathInRootFile+"/RoI").c_str(), LevelOfDetail, eventsBlock );
  HistoBooker* RoI_Booker = new HistoBooker(&JEM_RoI, &mLog, m_DataType);

  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
      Helper* Help = new Helper();
	  
      m_h_JEMHits_MainHits = DAQ_Booker->book1F("MainHits_JEM_DAQ", "Main Jet Hit Multiplicity per Threshold  --  JEM DAQ", 8, -0.5, 7.5,  "Threshold No.", "#");
      m_h_JEMHits_FwdHitsRight = DAQ_Booker->book1F("FwdHitsRight_JEM_DAQ", "Forward Right Jet Hit Multiplicity per Threshold  --  JEM DAQ", 4, -0.5, 3.5, "Threshold No.", "#");
      m_h_JEMHits_FwdHitsLeft = DAQ_Booker->book1F("FwdHitsLeft_JEM_DAQ", "Forward Left Jet Hit Multiplicity per Threshold  --  JEM DAQ", 4, -0.5, 3.5, "Threshold No.", "#");

      m_h_JEMEtSums_Ex = DAQ_Booker->book1F("Ex_JEM_DAQ", "JEM E_{x}  --  JEM DAQ", 250, 0,250, "Ex [GeV]", "#");
      m_h_JEMEtSums_Ey = DAQ_Booker->book1F("Ey_JEM_DAQ", "JEM E_{y}  --  JEM DAQ", 250, 0,250, "Ex [GeV]", "#");
      m_h_JEMEtSums_Et = DAQ_Booker->book1F("Et_JEM_DAQ", "JEM E_{t}  --  JEM DAQ", 250, 0,250, "Ex [GeV]", "#");

      // JEM RoI
      m_h_JEMRoI_MainHits = RoI_Booker->book1F("MainHits_JEM_RoI", "Main Jet Hit Multiplicity per Threshold  --  JEM RoI", 8, -0.5, 7.5,  "Threshold No.", "#");      
      m_h_JEMRoI_FwdHitsRight = RoI_Booker->book1F("FwdHitsRight_JEM_RoI", "Forward Right Jet Hit Multiplicity per Threshold  --  JEM RoI", 4, -0.5, 3.5, "Threshold No.", "#");
      m_h_JEMRoI_FwdHitsLeft = RoI_Booker->book1F("FwdHitsLeft_JEM_RoI", "Forward Left Jet Hit Multiplicity per Threshold  --  JEM RoI", 4, -0.5, 3.5, "Threshold No.", "#");

      m_h_JEMRoI_Thresh1_etaphi = RoI_Booker->book2F("Thresh1_eta-phi_JEM_RoI", "#eta - #phi Map of Hits passing Threshold <ThreshNo>  --  JEM RoI", 50, -5, 5, 32,0,6.4, "#eta", "#phi");
      m_h_JEMRoI_Thresh1_etaphi->SetBins(32,Help->JEEtaBinning(),32,Help->JEPhiBinning());

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
      mLog << MSG::INFO << "No JEMHits found in TES at "
	   << m_JEMHitsLocation << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" JEM Hits ---------------"<<endreq;

  // Step over all cells and put into hist
  JEMHitsCollection::const_iterator it_JEMHits ;

  for( it_JEMHits  = JEMHits ->begin(); it_JEMHits < JEMHits -> end(); ++it_JEMHits )
    {	  
      std::string JEMHit = Help->Binary((*it_JEMHits)-> JetHits(),24);
      // the binary hit information is represented by an integer number, 
      //that is converted to a string in order
      // to get the real binary information
      // later the multiplicities of the several thresholds are retrieved from this string

      mLog<<MSG::DEBUG<<"Crate: "<< (*it_JEMHits)->crate()<<"  Module: "<<(*it_JEMHits)->module()
	  << "  JetHits: " <<JEMHit<<   endreq;
      
      if  ((*it_JEMHits)->forward()==0) //Main Jets
	{
	  Help->FillHitsHisto(m_h_JEMHits_MainHits,JEMHit, 0, 8, 0, 3, &mLog);
	}

      if  ((*it_JEMHits)->forward()==1) //fwd jets a bit complicated!
	// fwd and main hits are contained in the same hitword
	{
	  if (((*it_JEMHits)-> module()==0) or((*it_JEMHits)-> module()==8) )//left fwd hits
	    // JEMs No 0 and 8 are processing forward left hits,
	    // JEMs No 7 and 15 forward right hits
	    {
	      Help->FillHitsHisto(m_h_JEMHits_FwdHitsLeft, JEMHit, 0, 4, 8, 2, &mLog);
	    }
	  if (((*it_JEMHits)-> module()==7) or((*it_JEMHits)-> module()==15) )//right fwd hits
	    {
	      Help->FillHitsHisto(m_h_JEMHits_FwdHitsRight, JEMHit, 0, 4, 8, 2, &mLog);
	    }
	  Help->FillHitsHisto(m_h_JEMHits_MainHits , JEMHit, 0, 8, 0, 2, &mLog);
	}
    }

  //******************************* 
  //********* JEM Et Sums ********* 
  //******************************* 

  const JEMEtSumsCollection* JEMEtSums;
  sc = m_storeGate->retrieve(JEMEtSums, m_JEMEtSumsLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JEMEtSums found in TES at "<< m_JEMEtSumsLocation << endreq ;
      return StatusCode::SUCCESS;
    }

  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" JEM Et Sums ---------------"<<endreq;

  // Step over all cells and put into hist
  JEMEtSumsCollection::const_iterator it_JEMEtSums ;
  LVL1::QuadLinear expand;
  //  expand = new QuadLinear();
  int ex;
  int ey;
  int et;

  for( it_JEMEtSums  = JEMEtSums ->begin(); it_JEMEtSums < JEMEtSums -> end(); ++it_JEMEtSums )
    {	       
      ex=expand.Expand( (*it_JEMEtSums)-> Ex());
      ey=expand.Expand( (*it_JEMEtSums)-> Ey());
      et=expand.Expand( (*it_JEMEtSums)-> Et() );

      m_h_JEMEtSums_Ex -> Fill(ex, 1.); 
      m_h_JEMEtSums_Ey -> Fill(ey, 1.); 
      m_h_JEMEtSums_Et -> Fill(et, 1.); 
      mLog <<MSG::DEBUG<< " JEMEtSums Crate: "<<(*it_JEMEtSums)->crate()<<"  Module: "<<(*it_JEMEtSums)->module()
	   <<"   Ex: "<<  ex
	   <<"   Ey: "<<  ey 
	   <<"   Et: "<<  et  << "    Et compressed: "<< (*it_JEMEtSums)-> Et() <<endreq;
    }   


  //******************************* 
  //********* JEM RoI ************* 
  //******************************* 
  const JemRoiCollection* JEMRoIs = 0;
  sc = m_storeGate->retrieve (JEMRoIs, m_JEMRoILocation);
    	 
  if (sc==StatusCode::FAILURE)
    {
      mLog <<MSG::INFO<<"No JEM RoIs found in TES at"<< m_JEMRoILocation<<endreq;
      return StatusCode::SUCCESS;    
    }

  mLog<<MSG::DEBUG<<"-------------- "<< m_DataType<<" JEM RoIs ---------------"<<endreq;

  // Step over all cells and put into hist
  JemRoiCollection::const_iterator it_JEMRoIs ;

  for( it_JEMRoIs  = JEMRoIs ->begin(); it_JEMRoIs < JEMRoIs -> end(); ++it_JEMRoIs )
    {	  
      std::string JEMRoIHits ;
      if  ((*it_JEMRoIs)->forward()==0) //Main Jets
	{
	  JEMRoIHits = Help->Binary((*it_JEMRoIs)-> hits(),8);
	  Help->FillHitsHisto(m_h_JEMRoI_MainHits , JEMRoIHits, 0, 8, 0, 1, &mLog);
	}      
      else 
	{
	  JEMRoIHits = Help->Binary((*it_JEMRoIs)-> hits(),4);

	  if (((*it_JEMRoIs)-> jem()==0) or((*it_JEMRoIs)-> jem()==8) )//left fwd hits
	    // JEMs No 0 and 8 are processing forward left hits,
	    // JEMs No 7 and 15 forward right hits
	    {
	      Help->FillHitsHisto(m_h_JEMRoI_FwdHitsLeft, JEMRoIHits, 0, 4, 0, 1, &mLog);
	    }
	  if (((*it_JEMRoIs)-> jem()==7) or((*it_JEMRoIs)-> jem()==15) )//right fwd hits
	    // JEMs No 0 and 8 are processing forward left hits,
	    // JEMs No 7 and 15 forward right hits
	    {
	      Help->FillHitsHisto(m_h_JEMRoI_FwdHitsRight, JEMRoIHits, 0, 4, 0, 1, &mLog);
	    }
	}
      
      mLog<<MSG::DEBUG<<"JEMRoI Word: "<<Help->Binary((*it_JEMRoIs)->roiWord(),32)<<endreq;
      mLog<<MSG::DEBUG<<"Crate: "<<(*it_JEMRoIs)->crate()<<"; JEM: "<<(*it_JEMRoIs)->jem()
	  <<"; forward: "<<(*it_JEMRoIs)->forward() <<"; Hits: "<<JEMRoIHits<<endreq;
      //mLog<<"Frame: "<<(*it_JEMRoIs)->frame()Location: "<<(*it_JEMRoIs)->location()<<"<<endreq;

      //RoI-Hits per Threshold and eta-phi
      LVL1::JEPRoIDecoder decoder;
      LVL1::CoordinateRange coordRange = decoder.coordinate((*it_JEMRoIs)->roiWord());
      double eta = coordRange.eta();
      double phi = coordRange.phi();

      JEMRoIHits = Help->Binary((*it_JEMRoIs)-> hits(),8);

      if ((Help->Multiplicity(JEMRoIHits,1,1))!=0)
	{
	  m_h_JEMRoI_Thresh1_etaphi->Fill(eta,phi,1);
	}
      
    }

  mLog<<MSG::DEBUG<<"--------------------------------------"<<endreq;

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

// ********************************************************************
//
// NAME:        JetElementMon.cxx
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the inputdata (JetElements) of the JEP
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


#include "TrigT1Calo/JetElementMaker.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "AthenaMonitoring/AthenaMonManager.h"
#include "TrigT1CaloMonitoring/JetElementMon.h"


// *********************************************************************
// Public Methods
// *********************************************************************

/*---------------------------------------------------------*/
JetElementMon::JetElementMon( const std::string & type, const std::string & name,
			      const IInterface* parent )
  : ManagedMonitorToolBase( type, name, parent )
/*---------------------------------------------------------*/
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "JetElementLocation", m_JetElementLocation = LVL1::TrigT1CaloDefs::JetElementLocation); 
  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/JetElements") ;
  declareProperty( "DataType", m_DataType="") ;
}

/*---------------------------------------------------------*/
JetElementMon::~JetElementMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode JetElementMon::bookHistograms( bool isNewEventsBlock, 
					  bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in JetElementMon::bookHistograms" << endreq;
  
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

  MonGroup JetElements_expert (this,m_PathInRootFile ,expert, eventsBlock);
  HistoBooker* expert_Booker = new HistoBooker(&JetElements_expert, &mLog, m_DataType);

  MonGroup JetElements_shift (this,m_PathInRootFile ,shift, eventsBlock);
  HistoBooker* shift_Booker = new HistoBooker(&JetElements_shift, &mLog, m_DataType);
  
  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
      Helper* Help = new Helper();

      // register Histograms for  data
      m_h_je_eta = expert_Booker->book1F("eta_JEM_input", "JE eta  --  JEM input" , 50, -5, 5, "#eta" , "#");
      m_h_je_eta ->SetBins(32,Help->JEEtaBinning());

      m_h_je_phi = expert_Booker->book1F("phi_JEM_input", "JE phi  --  JEM input", 32, 0, 6.4, "#phi" , "#");
      m_h_je_phi->SetBins(32,Help->JEPhiBinning());

      m_h_je_emenergy = expert_Booker->book1F("EmEnergy_JEM_input", "JE em energy  --  JEM input", 
					      100, 0, 100, "em energy [GeV]" , "#");
      m_h_je_hadenergy  = expert_Booker->book1F("HadEnergy_JEM_input", "JE had energy  --  JEM input", 
						100, 0, 100, "had energy [GeV]" , "#");
      m_h_je_energy = expert_Booker->book1F("Et_JEM_input", "JE Et (em + had)  --  JEM input", 
					    100, 0, 100, "(em+had) energy [GeV]" , "#");


      m_h_je_etaphi = shift_Booker->book2F("eta-phi_JEM_input", "JE per eta - phi  --  JEM input", 
					   50, -5, 5, 32, 0, 6.4 , "#eta", "#phi");
      m_h_je_etaphi->SetBins(32,Help->JEEtaBinning(),32,Help->JEPhiBinning());

      m_h_je_energy_etaphi = shift_Booker->book2F("Et_eta-phi_JEM_input", "JE Et per eta - phi  --  JEM input",  
						  50, -5, 5, 32, 0, 6.4, "#eta", "#phi");	  
      m_h_je_etaphi->SetBins(32,Help->JEEtaBinning(),32,Help->JEPhiBinning());
    }
	
  if( isNewRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}


/*---------------------------------------------------------*/
StatusCode JetElementMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );

  //fill Histos
  const JECollection* jetElements;
  StatusCode sc = m_storeGate->retrieve(jetElements, m_JetElementLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO
	   << "No JetElements found in TES at "
	   << m_JetElementLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  // Step over all cells and put into hist
  JECollection::const_iterator it_je ;
  for( it_je = jetElements ->begin(); it_je < jetElements->end(); ++it_je )
    {	  
      mLog << MSG::VERBOSE<<m_DataType <<" JE has coords ("<<(*it_je)->phi()<<", "<<(*it_je)->eta()
	   << " and energies : "<<(*it_je)->emEnergy()<<", "<<(*it_je)->hadEnergy()<<" (Em,Had)"<<endreq;

      // fill histograms
      m_h_je_eta -> Fill( (*it_je)-> eta(), 1.);
      m_h_je_phi->Fill( (*it_je)->phi() , 1.);
      m_h_je_emenergy->Fill( (*it_je)->emEnergy() , 1.);
      m_h_je_hadenergy->Fill( (*it_je)->hadEnergy() , 1.);
      m_h_je_energy->Fill( (*it_je)->energy() , 1.);
      
      m_h_je_etaphi->Fill( (*it_je)->eta(),(*it_je)->phi() , 1.);
      m_h_je_energy_etaphi->Fill( (*it_je)->eta(),(*it_je)->phi() ,(*it_je)->energy() );  
  
    }     	
  return StatusCode( StatusCode::SUCCESS );
}

/*---------------------------------------------------------*/
StatusCode JetElementMon::procHistograms( bool isEndOfEventsBlock, 
					  bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  if( isEndOfEventsBlock || isEndOfLumiBlock ) 
    {
      
    }
	
  if( isEndOfRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}

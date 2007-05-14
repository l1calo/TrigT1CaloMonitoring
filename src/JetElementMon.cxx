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


//#include "TrigT1Calo/JEMEtSums.h"
#include "TrigT1Calo/JetElementMaker.h"

#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "AthenaMonitoring/AthenaMonManager.h"
#include "TrigT1CaloMonitoring/JetElementMon.h"


// *********************************************************************
// Public Methods
// *********************************************************************

JetElementMon::JetElementMon( const std::string & type, const std::string & name,
                        const IInterface* parent )
	: ManagedMonitorToolBase( type, name, parent )
{
  // This is how you declare the parameters to Gaudi so that
  // they can be over-written via the job options file

  declareProperty( "BS_JetElementLocation", m_BS_JetElementLocation = LVL1::TrigT1CaloDefs::JetElementLocation); 
  declareProperty( "Sim_JetElementLocation", m_Sim_JetElementLocation  =  LVL1::TrigT1CaloDefs::JetElementLocation ) ;

}


JetElementMon::
~JetElementMon()
{
}


//_______________________________ book Histograms ___________________________________________
StatusCode
JetElementMon::
bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
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

  MonGroup BS_JetElements_expert (this, "Stats/BS_JetElements",expert, eventsBlock);
  MonGroup BS_JetElements_shift (this, "Stats/BS_JetElements",shift, eventsBlock);
  
  MonGroup Sim_JetElements_expert (this, "Stats/Sim_JetElements",expert, eventsBlock);
  MonGroup Sim_JetElements_shift (this, "Stats/Sim_JetElements",shift, eventsBlock);

  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
      // register Histograms for BS data
      m_h_BS_je_eta     = new TH1F("eta","BS JE eta", 50,-5,5);
      m_h_BS_je_eta -> GetXaxis() -> SetTitle("#eta");
      BS_JetElements_expert.regHist (m_h_BS_je_eta);
      
      m_h_BS_je_phi     = new TH1F("phi","BS JE phi",32,0,6.4);
      m_h_BS_je_phi -> GetXaxis() -> SetTitle("#phi");
      BS_JetElements_expert.regHist (m_h_BS_je_phi );

      m_h_BS_je_emenergy     = new TH1F("em energy","BS JE em energy",100,0,100);
      m_h_BS_je_emenergy -> GetXaxis() -> SetTitle("em energy [GeV]");
      BS_JetElements_expert.regHist (m_h_BS_je_emenergy);

      m_h_BS_je_hadenergy     = new TH1F("had energy","BS JE had energy",100,0,100);
      m_h_BS_je_hadenergy -> GetXaxis() -> SetTitle("had energy [GeV]");
      BS_JetElements_expert.regHist (m_h_BS_je_hadenergy);

      m_h_BS_je_energy     = new TH1F("Et","BS JE Et (em + had)",100,0,100);
      m_h_BS_je_energy -> GetXaxis() -> SetTitle("(em+had) energy [GeV]");
      BS_JetElements_expert.regHist (m_h_BS_je_energy );

      m_h_BS_je_etaphi      = new TH2F("eta phi", "BS JE per eta - phi", 50,-5,5,32,0,6.4);
      m_h_BS_je_etaphi -> SetOption ("colz");
      m_h_BS_je_etaphi -> GetXaxis() -> SetTitle("#eta");
      m_h_BS_je_etaphi -> GetYaxis() -> SetTitle("#phi");
      BS_JetElements_shift.regHist (m_h_BS_je_etaphi );

      m_h_BS_je_energy_etaphi= new TH2F("Et per eta-phi", "BS JE Et per eta - phi", 50,-5,5,32,0,6.4);
      m_h_BS_je_energy_etaphi -> SetOption ("colz");
      m_h_BS_je_energy_etaphi -> GetXaxis() -> SetTitle("#eta");
      m_h_BS_je_energy_etaphi -> GetYaxis() -> SetTitle("#phi");
      BS_JetElements_shift.regHist (m_h_BS_je_energy_etaphi);
	  
      // register Histograms for Simulation data
      m_h_Sim_je_eta     = new TH1F("eta","Sim JE eta",50,-5,5);
      m_h_Sim_je_eta -> GetXaxis() -> SetTitle("#eta");
      Sim_JetElements_expert.regHist (m_h_Sim_je_eta);
      
      m_h_Sim_je_phi     = new TH1F("phi","Sim JE phi",32,0,6.4);
      m_h_Sim_je_phi -> GetXaxis() -> SetTitle("#phi");
      Sim_JetElements_expert.regHist (m_h_Sim_je_phi );
      
      m_h_Sim_je_emenergy     = new TH1F("em energy","Sim JE em energy",100,0,100);
      m_h_Sim_je_emenergy -> GetXaxis() -> SetTitle("em energy [GeV]");
      Sim_JetElements_expert.regHist (m_h_Sim_je_emenergy);

      m_h_Sim_je_hadenergy     = new TH1F("had energy","Sim JE had energy",100,0,100);
      m_h_Sim_je_hadenergy -> GetXaxis() -> SetTitle("had energy [GeV]");
      Sim_JetElements_expert.regHist (m_h_Sim_je_hadenergy);

      m_h_Sim_je_energy     = new TH1F("Et","Sim JE Et (em + had)",100,0,100);
      m_h_Sim_je_energy -> GetXaxis() -> SetTitle("(em+had) energy [GeV]");
      Sim_JetElements_expert.regHist (m_h_Sim_je_energy );

      m_h_Sim_je_etaphi      = new TH2F("eta-phi", "Sim JE per eta - phi", 50,-5,5,32,0,6.4);
      m_h_Sim_je_etaphi -> SetOption ("colz");
      m_h_Sim_je_etaphi -> GetXaxis() -> SetTitle("#eta");
      m_h_Sim_je_etaphi -> GetYaxis() -> SetTitle("#phi");
      Sim_JetElements_shift.regHist (m_h_Sim_je_etaphi );

      m_h_Sim_je_energy_etaphi= new TH2F("Et per eta-phi", "Sim JE Et per eta - phi", 50,-5,5,32,0,6.4);
      m_h_Sim_je_energy_etaphi -> SetOption ("colz");
      m_h_Sim_je_energy_etaphi -> GetXaxis() -> SetTitle("#eta");
      m_h_Sim_je_energy_etaphi -> GetYaxis() -> SetTitle("#phi");
      Sim_JetElements_shift.regHist (m_h_Sim_je_energy_etaphi);

    }
	
  if( isNewRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}




//_______________________________ fill Histograms ___________________________________________
StatusCode
JetElementMon::
fillHistograms()
{
  MsgStream mLog( msgSvc(), name() );

  //fill BS Histos
  const JECollection* jetElements;
  StatusCode sc = m_storeGate->retrieve(jetElements, m_BS_JetElementLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::DEBUG
	   << "No jetElemtns found in TES at "
	   << m_BS_JetElementLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  // Step over all cells and put into hist
  JECollection::const_iterator it_je ;
  for( it_je = jetElements ->begin(); it_je < jetElements->end(); ++it_je )
    {	  
      mLog << MSG::VERBOSE<<"JE has coords ("<<(*it_je)->phi()<<", "<<(*it_je)->eta()
	   << " and energies : "<<(*it_je)->emEnergy()<<", "<<(*it_je)->hadEnergy()<<" (Em,Had)"<<endreq;

      // fill histograms
      m_h_BS_je_eta -> Fill( (*it_je)-> eta(), 1.);
      m_h_BS_je_phi->Fill( (*it_je)->phi() , 1.);
      m_h_BS_je_emenergy->Fill( (*it_je)->emEnergy() , 1.);
      m_h_BS_je_hadenergy->Fill( (*it_je)->hadEnergy() , 1.);
      m_h_BS_je_energy->Fill( (*it_je)->energy() , 1.);
      
      m_h_BS_je_etaphi->Fill( (*it_je)->eta(),(*it_je)->phi() , 1.);
      m_h_BS_je_energy_etaphi->Fill( (*it_je)->eta(),(*it_je)->phi() ,(*it_je)->energy() );  
  
    }     	
  
  //fill Simulation data
  sc = m_storeGate->retrieve(jetElements, m_Sim_JetElementLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::VERBOSE
	   << "No jetElemtns found in TES at "
	   << m_Sim_JetElementLocation
	   << endreq ;
      return StatusCode::SUCCESS;
    }

  // Step over all cells and put into hist
  for( it_je  =jetElements ->begin(); it_je < jetElements->end(); ++it_je )
    {	  
      mLog << MSG::VERBOSE<<"JE has coords ("<<(*it_je)->phi()<<", "<<(*it_je)->eta()
	   << " and energies : "<<(*it_je)->emEnergy()<<", "<<(*it_je)->hadEnergy()<<" (Em,Had)"<<endreq;

      // fill histograms
      m_h_Sim_je_eta -> Fill( (*it_je)-> eta(), 1.);
      m_h_Sim_je_phi->Fill( (*it_je)->phi() , 1.);
      m_h_Sim_je_emenergy->Fill( (*it_je)->emEnergy() , 1.);
      m_h_Sim_je_hadenergy->Fill( (*it_je)->hadEnergy() , 1.);
      m_h_Sim_je_energy->Fill( (*it_je)->energy() , 1.);
      
      m_h_Sim_je_etaphi->Fill( (*it_je)->eta(),(*it_je)->phi() , 1.);
      m_h_Sim_je_energy_etaphi->Fill( (*it_je)->eta(),(*it_je)->phi() ,(*it_je)->energy() );  
  
    }
  
  return StatusCode( StatusCode::SUCCESS );
}

//_______________________________ proc  Histograms ___________________________________________
StatusCode
JetElementMon::
procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
{
        if( isEndOfEventsBlock || isEndOfLumiBlock ) 
	  {

	}
	
	if( isEndOfRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}

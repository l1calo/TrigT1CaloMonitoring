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
#include "TrigT1Calo/DataError.h"
#include "TrigT1Calo/CoordToHardware.h"

#include "TrigT1Interfaces/TrigT1CaloDefs.h"
#include "TrigT1Interfaces/Coordinate.h"

#include "AthenaMonitoring/AthenaMonManager.h"
#include "TrigT1CaloMonitoring/JetElementMon.h"


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
  declareProperty( "ErrorPathInRootFile", m_ErrorPathInRootFile="Stats/L1Calo/Errors") ;
  declareProperty( "NumberOfSlices", m_SliceNo = 5);
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
      mLog << MSG::ERROR << "Unable to retrieve pointer to StoreGateSvc"<< endreq;
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

  MonGroup JetElements_expert (this,m_PathInRootFile ,expert, eventsBlock);
  HistoBooker* expert_Booker = new HistoBooker(&JetElements_expert, &mLog, m_DataType);

  MonGroup JetElements_shift (this,m_PathInRootFile ,LevelOfDetail, eventsBlock);
  HistoBooker* shift_Booker = new HistoBooker(&JetElements_shift, &mLog, m_DataType);
  
  MonGroup JE_TriggeredSlice( this, (m_PathInRootFile).c_str(), shift, eventsBlock );
  HistoBooker* TriggeredSlice_Booker = new HistoBooker(&JE_TriggeredSlice, &mLog, m_DataType);

  MonGroup JEM_Error( this, (m_ErrorPathInRootFile).c_str(), shift, eventsBlock );
  HistoBooker* Error_Booker = new HistoBooker(&JEM_Error, &mLog, "");

  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
      Helper* Help = new Helper();

      // register Histograms for  data
      m_h_je_eta = expert_Booker->book1F("eta_JEM_input", "JE distribution per #eta  --  JEM input" , 50, -5, 5, "#eta" , "#");
      m_h_je_eta ->SetBins(32,Help->JEEtaBinning());
      m_h_je_phi = expert_Booker->book1F("phi_JEM_input", "JE distribution per #phi  --  JEM input", 32, 0, 6.4, "#phi" , "#");
      m_h_je_phi->SetBins(32,Help->JEPhiBinning());


      m_h_je_emenergy = expert_Booker->book1F("EmEnergy_JEM_input", "JE EM energy distribution  --  JEM input", 100, 0, 100, "em energy [GeV]" , "#");
      m_h_je_hadenergy  = expert_Booker->book1F("HadEnergy_JEM_input", "JE HAD energy distribution  --  JEM input", 100, 0, 100, "had energy [GeV]" , "#");


      m_h_je_energy_emHitMap = shift_Booker->book2F("JE_EM_HitMap_energy_JEM_input", "#eta - #phi map of EM JE weighted with energy  --  JEM input", 50, -5, 5, 32, 0, 6.4 , "#eta", "#phi");
      m_h_je_energy_emHitMap->SetBins(32,Help->JEEtaBinning(),32,Help->JEPhiBinning());

      m_h_je_energy_hadHitMap = shift_Booker->book2F("JE_HAD_HitMap_energy_JEM_input", "#eta - #phi map of HAD JE weighted with energy  --  JEM input", 50, -5, 5, 32, 0, 6.4, "#eta", "#phi");	  
      m_h_je_energy_hadHitMap->SetBins(32,Help->JEEtaBinning(),32,Help->JEPhiBinning());

      if (m_DataType=="BS")
	{
	  std::string name,title;
	  std::stringstream buffer;
	  
	  for (int i = 0; i < m_SliceNo; i++)
	    {
	      buffer.str("");
	      buffer<<i;
	      
	      name = "JE_EM_HitMap_" + buffer.str() + "_JEM_input";
	      title = "#eta - #phi map of EM JE for Timeslice " + buffer.str() +  "  --  JEM input";
	      m_h_je_emHitMap[i]=shift_Booker->book2F(name,title,50, -5, 5, 32, 0, 6.4, "#eta", "#phi");	  
	      m_h_je_emHitMap[i]->SetBins(32,Help->JEEtaBinning(),32,Help->JEPhiBinning());
	      
	      buffer.str("");
	      buffer<<i;
	      
	      name = "JE_HAD_HitMap_" + buffer.str() + "_JEM_input";
	      title = "#eta - #phi map of HAD JE for Timeslice " + buffer.str() +  "  --  JEM input";
	      m_h_je_hadHitMap[i]=shift_Booker->book2F(name,title,50, -5, 5, 32, 0, 6.4, "#eta", "#phi");	  
	      m_h_je_hadHitMap[i]->SetBins(32,Help->JEEtaBinning(),32,Help->JEPhiBinning());
	    }

	  // ----------------------------------- Error Histos ------------------------------------------------------
	  m_h_je_error = Error_Booker->book2F("JEM_Error","JEM S-Link Error per Module and Crate",12,0.5,12.5,35,0.5,35.5,"","");
	  //m_h_je_error -> SetOption ("text");
	  m_h_je_error->GetXaxis()->SetBinLabel(1, "EM Parity");
	  m_h_je_error->GetXaxis()->SetBinLabel(2, "HAD Parity");
	  m_h_je_error->GetXaxis()->SetBinLabel(3, "PPM Link down");

	  m_h_je_error->GetXaxis()->SetBinLabel(5, "GLinkParity");
	  m_h_je_error->GetXaxis()->SetBinLabel(6, "GLinkProtocol");
	  m_h_je_error->GetXaxis()->SetBinLabel(7, "BCNMismatch");
	  m_h_je_error->GetXaxis()->SetBinLabel(8, "FIFOOverflow");
	  m_h_je_error->GetXaxis()->SetBinLabel(9, "ModuleError");
	  m_h_je_error->GetXaxis()->SetBinLabel(10, "GLinkDown");
	  m_h_je_error->GetXaxis()->SetBinLabel(11, "GLinkTimeout");
	  m_h_je_error->GetXaxis()->SetBinLabel(12, "FailingBCN");
      
	  for (int i = 0; i < 16; i++)
	    {
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_je_error->GetYaxis()->SetBinLabel((i+1), name.c_str());
	      
	      buffer.str("");
	      buffer<<i;
	      name = "JEM " + buffer.str();
	      m_h_je_error->GetYaxis()->SetBinLabel((i+19), name.c_str());
	    }
	  m_h_je_error->GetYaxis()->SetBinLabel(17, "Crate 0: ");
	  m_h_je_error->GetYaxis()->SetBinLabel(35, "Crate 1: ");
	}



      // number of triggered slice
      m_h_je_triggeredSlice=TriggeredSlice_Booker->book1F("JE_TriggeredSlice","Number of the Triggered Slice for JE",7,-0.5,6.5,"#Slice");
}
	
  if( isNewRun ) { }
  
  return StatusCode( StatusCode::SUCCESS );
}


/*---------------------------------------------------------*/
StatusCode JetElementMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );

  // =============================================================================================
  // ================= Container: JetElements ====================================================
  // =============================================================================================
  // retrieve JetElements
  const JECollection* jetElements;
  StatusCode sc = m_storeGate->retrieve(jetElements, m_JetElementLocation);

  if( (sc==StatusCode::FAILURE) ) 
    {
      mLog << MSG::INFO << "No JetElements found in TES at " << m_JetElementLocation << endreq ;
      return StatusCode::SUCCESS;
    }

  // Step over all cells 
  JECollection::const_iterator it_je ;
  for( it_je = jetElements ->begin(); it_je < jetElements->end(); ++it_je )
    {	  
      mLog << MSG::VERBOSE<<m_DataType <<" JE has coords ("<<(*it_je)->phi()<<", "<<(*it_je)->eta()
	   << " and energies : "<<(*it_je)->emEnergy()<<", "<<(*it_je)->hadEnergy()<<" (Em,Had)"<<endreq;

      m_h_je_eta -> Fill( (*it_je)-> eta(), 1.);
      m_h_je_phi->Fill( (*it_je)->phi() , 1.);

      m_h_je_emenergy->Fill( (*it_je)->emEnergy() , 1.);
      m_h_je_hadenergy->Fill( (*it_je)->hadEnergy() , 1.);
      
      m_h_je_energy_emHitMap->Fill( (*it_je)->eta(),(*it_je)->phi() , (*it_je)->emEnergy());
      m_h_je_energy_hadHitMap->Fill( (*it_je)->eta(),(*it_je)->phi() ,(*it_je)->hadEnergy() ); 

      // number of triggered slice
      m_h_je_triggeredSlice->Fill((*it_je)->peak(),1);

      // ------------------------------------------------------------------------------------------
      // ----------------- Histos filled only for BS data -----------------------------------------
      // ------------------------------------------------------------------------------------------
      if (m_DataType=="BS")
	{
	  // ----------------- HitMaps per time slice -----------------------------------------
	  for (int i = 0; i < m_SliceNo; i++)
	    {
	      if (i < ((*it_je)->emEnergyVec()).size())
		{
		  if ((*it_je)->emEnergyVec()[i] != 0) m_h_je_emHitMap[i]-> Fill( (*it_je)->eta(),(*it_je)->phi() ,1);
		  if ((*it_je)->hadEnergyVec()[i] != 0) m_h_je_hadHitMap[i]-> Fill( (*it_je)->eta(),(*it_je)->phi() ,1);
		} 
	    }

	  // ----------------- Error Histos -----------------------------------------
	  LVL1::DataError err((*it_je)->emError());
	  LVL1::DataError haderr((*it_je)->hadError());
	  LVL1::CoordToHardware ToHW;
	  LVL1::Coordinate coord((*it_je)->phi(),(*it_je)->eta());

	  int crate = ToHW.jepCrate(coord);
	  int module=ToHW.jepModule(coord);
	  // EM Parity
	  m_h_je_error->Fill(1,(crate*18 + module+1),err.get(1));
	  // HAD Parity
	  m_h_je_error->Fill(1,(crate*18 + module+1),haderr.get(1));
	  // PPM Link down
	  m_h_je_error->Fill(3,(crate*18 + module+1),err.get(2));

	  // GLinkParity
	  m_h_je_error->Fill(5,(crate*18 + module+1),err.get(16));
          // GLinkProtocol
	  m_h_je_error->Fill(6,(crate*18 + module+1),err.get(17));
	  // BCNMismatch
	  m_h_je_error->Fill(7,(crate*18 + module+1),err.get(18));
	  // FIFOOverflow
	  m_h_je_error->Fill(8,(crate*18 + module+1),err.get(19));
	  // Module Error
	  m_h_je_error->Fill(9,(crate*18 + module+1),err.get(20));

	  // GLinkDown
	  m_h_je_error->Fill(10,(crate*18 + module+1),err.get(22));
	  // GLinkTimeout
	  m_h_je_error->Fill(11,(crate*18 + module+1),err.get(23));
	  // FailingBCN
	  if (err.get(24)!=0) m_h_je_error->Fill(12,(crate*18 + module+1),1);
	}   
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

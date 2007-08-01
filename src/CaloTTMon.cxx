// ********************************************************************
//
// NAME:     CaloTTMon.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
//
// ********************************************************************

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ITHistSvc.h"

#include <TH1D.h>
#include <TH2D.h>

#include "TString.h"

#include "StoreGate/StoreGateSvc.h"

#include "CaloEvent/CaloTower.h"
#include "CaloEvent/CaloTowerContainer.h"

#include "CaloEvent/CaloCell.h"
#include "CaloEvent/CaloCellContainer.h"

#include "TrigT1CaloMonitoring/CaloTTMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "TrigT1Calo/TriggerTowerCollection.h"
#include "TrigT1Calo/TrigT1CaloDict.h"
#include "TrigT1Calo/TriggerTower_ClassDEF.h"

#include "GaudiKernel/IService.h"
#include "GaudiKernel/IToolSvc.h"


/*---------------------------------------------------------*/
CaloTTMon::CaloTTMon(const std::string & type, const std::string & name,
					 const IInterface* parent)
  : ManagedMonitorToolBase ( type, name, parent )
/*---------------------------------------------------------*/
{
  declareProperty("CaloTTContainer",  m_CaloTTContainerName = "AllCalos");
  declareProperty("Calo_HitMap_Thresh0",  m_CaloTT_HitMap_Thresh0 = 1);
  declareProperty("Calo_HitMap_Thresh1",  m_CaloTT_HitMap_Thresh1 = 3);
  declareProperty("Calo_HitMap_Thresh2",  m_CaloTT_HitMap_Thresh2 = 7);
 
  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/CMM") ;
  declareProperty( "DataType", m_DataType="") ;
}

/*---------------------------------------------------------*/
CaloTTMon::~CaloTTMon()
/*---------------------------------------------------------*/
{
}


/*---------------------------------------------------------*/
StatusCode CaloTTMon::bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "in CaloTTMon::bookHistograms" << endreq;

  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR << "Unable to retrieve pointer to StoreGateSvc" << endreq;
      return sc;
    }
  
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  MonGroup Calo_Tower( this, (m_PathInRootFile+"/CaloReadout").c_str(), shift, eventsBlock );
  HistoBooker* CaloTower_Booker = new HistoBooker(&Calo_Tower, &log, m_DataType);

  MonGroup Calo_HitMaps( this, (m_PathInRootFile+"/CaloTTHitMaps").c_str(), shift, eventsBlock );
  HistoBooker* HitMaps_Booker = new HistoBooker(&Calo_HitMaps, &log, m_DataType);

  MonGroup Calo_Dist( this, (m_PathInRootFile+"/CaloTTEnergyDistribution").c_str(), shift, eventsBlock );
  HistoBooker* EnergyDistribution_Booker = new HistoBooker(&Calo_Dist, &log, m_DataType);

  if( isNewEventsBlock || isNewLumiBlock ) 
    {	

      Helper* Help = new Helper();
      std::string name,title;
      std::stringstream buffer, etabuffer,phibuffer;


      m_h_CaloTower_et=CaloTower_Booker->book1F("CaloTower_et","Distribution of CaloTower E_{T}",256,-0.5,255.5,"E_{T} in GeV");
      m_h_CaloTower_eta=CaloTower_Booker->book1F("CaloTower_eta","#eta - Distribution of CaloTowers",100,-4.9,4.9,"#eta");
      m_h_CaloTower_phi=CaloTower_Booker->book1F("CaloTower_phi","#phi - Distribution of CaloTowers",64,-M_PI,M_PI,"#phi");      m_h_CaloTower_HitMap=CaloTower_Booker->book2F("CaloTower_HitMap","#eta - #phi Map of CaloTowers",100,-4.9,4.9, 64,-M_PI,M_PI,"#eta","#phi");

      
      m_h_CaloCell_et=CaloTower_Booker->book1F("CaloCell_et","Distribution of CaloCell E_{T}",256,-0.5,255.5,"E_{T} in GeV");
      m_h_CaloCell_eta=CaloTower_Booker->book1F("CaloCell_eta","#eta - Distribution of CaloCells",100,-4.9,4.9,"#eta");
      m_h_CaloCell_phi=CaloTower_Booker->book1F("CaloCell_phi","#phi - Distribution of CaloCells",64,-M_PI,M_PI,"#phi");
      m_h_CaloCell_HitMap=CaloTower_Booker->book2F("CaloCell_HitMap","#eta - #phi Map of CaloCells",100,-4.9,4.9, 64,-M_PI,M_PI,"#eta","#phi");


      //CaloTT Hitmaps per threshold
     buffer.str("");
     buffer<<m_CaloTT_HitMap_Thresh0;
   m_h_CaloTT_HitMap_emLUT_Thresh0=HitMaps_Booker->book2F("CaloTTEmLUT"+buffer.str(),"#eta - #phi Map of CaloTT EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_CaloTT_HitMap_emLUT_Thresh0->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
  
      buffer.str("");
      buffer<<m_CaloTT_HitMap_Thresh1;
   m_h_CaloTT_HitMap_emLUT_Thresh1=HitMaps_Booker->book2F("CaloTTEmLUT"+buffer.str(),"#eta - #phi Map of CaloTT EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_CaloTT_HitMap_emLUT_Thresh1->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());  
 
      buffer.str("");
      buffer<<m_CaloTT_HitMap_Thresh2;
   m_h_CaloTT_HitMap_emLUT_Thresh2=HitMaps_Booker->book2F("CaloTTEmLUT"+buffer.str(),"#eta - #phi Map of CaloTT EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_CaloTT_HitMap_emLUT_Thresh2->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
  

      buffer.str("");
      buffer<<m_CaloTT_HitMap_Thresh0;
   m_h_CaloTT_HitMap_hadLUT_Thresh0=HitMaps_Booker->book2F("CaloTTHadLUT"+buffer.str(),"#eta - #phi Map of CaloTT Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_CaloTT_HitMap_hadLUT_Thresh0->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());   

      buffer.str("");
      buffer<<m_CaloTT_HitMap_Thresh1;
   m_h_CaloTT_HitMap_hadLUT_Thresh1=HitMaps_Booker->book2F("CaloTTHadLUT"+buffer.str(),"#eta - #phi Map of CaloTT Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_CaloTT_HitMap_hadLUT_Thresh1->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());   

      buffer.str("");
      buffer<<m_CaloTT_HitMap_Thresh2;
   m_h_CaloTT_HitMap_hadLUT_Thresh2=HitMaps_Booker->book2F("CaloTTHadLUT"+buffer.str(),"#eta - #phi Map of CaloTT Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
   m_h_CaloTT_HitMap_hadLUT_Thresh2->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());  

   //distribution of LUT peak per detector region
   m_h_CaloTT_emLUT=EnergyDistribution_Booker->book1F("CaloTTEMLUT_Dist","CaloTT EM LUT Distribution of Peak",256,-0.5,255.5,"em LUT Peak [GeV]");
   m_h_CaloTT_emLUT_eta=EnergyDistribution_Booker->book1F("CaloTTEMLUT_eta","CaloTT EM LUT Distribution of Peak per #eta",21,-0.5,255.5,"#eta");
   m_h_CaloTT_emLUT_eta->SetBins(66,Help->TTEtaBinning());

   m_h_CaloTT_emLUT_phi=EnergyDistribution_Booker->book1F("CaloTTEMLUT_phi","CaloTT EM LUT Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
   m_h_CaloTT_emLUT_phi->SetBins(64,Help->TTPhiBinning());  


   m_h_CaloTT_hadLUT=EnergyDistribution_Booker->book1F("CaloTTHADLUT_Dist","CaloTT HAD LUT Distribution of Peak",256,-0.5,255.5,"had LUT Peak [GeV]"); 
   m_h_CaloTT_hadLUT_eta=EnergyDistribution_Booker->book1F("CaloTTHADLUT_eta","CaloTT HAD LUT Distribution of Peak per #eta",256,-0.5,255.5,"#eta");
   m_h_CaloTT_hadLUT_eta->SetBins(66,Help->TTEtaBinning());
   m_h_CaloTT_hadLUT_phi=EnergyDistribution_Booker->book1F("CaloTTHADLUT_phi","CaloTT HAD LUT Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
   m_h_CaloTT_hadLUT_phi->SetBins(64,Help->TTPhiBinning());  
    }
  if( isNewRun ) { }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CaloTTMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  
  log << MSG::DEBUG << "in fillHistograms()" << endreq;

  //Retrieve Calo Tower collection from SG
  const CaloTowerContainer* CaloTowerTES = 0; 
  StatusCode sc=m_storeGate->retrieve(CaloTowerTES, "CombinedTower"); 
  if( sc.isFailure()  ||  !CaloTowerTES ) {
    log << MSG::ERROR<< "No CaloTowerContainer found at CombinedTower" << endreq; 
  }

  //Retreive Calo Cell collection from SG
  const CaloCellContainer* CaloCellTES = 0;
  sc=m_storeGate->retrieve(CaloCellTES, "AllCalo"); 
  if( sc.isFailure()  ||  !CaloCellTES ) {
    log << MSG::ERROR<< "No CaloCellContainer found at AllCalo"<< endreq; 
  }
  

 // =============================================================================================
  // ================= CaloTowers (combined LAr and Tile data) From reconstruction ==============
  // ============================================================================================

  // CaloTowers exist irrespective of how low the Energy is, 
  // so need an Et cut to select only those towers with a deposited energy

  CaloTowerContainer::const_iterator CaloTowerIterator    = CaloTowerTES->begin();
  CaloTowerContainer::const_iterator CaloTowerIteratorEnd = CaloTowerTES->end();

  //int TTkey; //change this to key and use previous int
  //double TTtoCaloEnergyRatio;

 for (; CaloTowerIterator != CaloTowerIteratorEnd; ++CaloTowerIterator)
  {

    //select only Towers with an energy deposit
    if((*CaloTowerIterator)->et()>0.*GeV)
      {
	m_h_CaloTower_phi->Fill( (*CaloTowerIterator)->phi(), 1.);
	m_h_CaloTower_eta->Fill( (*CaloTowerIterator)->eta(), 1.);
	m_h_CaloTower_HitMap->Fill( (*CaloTowerIterator)->eta(), (*CaloTowerIterator)->phi(), 1.);
	m_h_CaloTower_et->Fill( (*CaloTowerIterator)->et(), 1.);
     }//eo calotower et>0gev
  }

  // =============================================================================================
  // ================= CaloCells  ==============
  // ============================================================================================

  // CaloCells exist irrespective of how low the Energy is, 
  // so need an Et cut to select only those towers with a deposited energy

  CaloCellContainer::const_iterator CaloCellIterator    = CaloCellTES->begin();
  CaloCellContainer::const_iterator CaloCellIteratorEnd = CaloCellTES->end();

  //int TTkey; //change this to key and use previous int
  //double TTtoCaloEnergyRatio;

 for (; CaloCellIterator != CaloCellIteratorEnd; ++CaloCellIterator)
  {

    //select only Cells with an energy deposit
    if((*CaloCellIterator)->et()>0.*GeV )
      {
	m_h_CaloCell_phi->Fill( (*CaloCellIterator)->phi(), 1.);
	m_h_CaloCell_eta->Fill( (*CaloCellIterator)->eta(), 1.);
	m_h_CaloCell_HitMap->Fill( (*CaloCellIterator)->eta(), (*CaloCellIterator)->phi(), 1.);
	m_h_CaloCell_et->Fill( (*CaloCellIterator)->et(), 1.);
     }//eo calotower et>0gev
  }









  //Retrieve Calo TriggerTowers from SG
  const TriggerTowerCollection* TriggerTowerTES = 0; 
  sc = m_storeGate->retrieve(TriggerTowerTES, m_CaloTTContainerName); 
  if( (sc==StatusCode::FAILURE) ) 
    {
      log << MSG::INFO << "No TriggerTower found in TES at "<< m_CaloTTContainerName<< endreq ;
      return StatusCode::SUCCESS;
    }


  // =========================================================================================
  // ================================== TriggerTower Plots ===================================
  // =========================================================================================

  // Calos are only filled/stored if they contain energy, 
  // so looping over all Calos will only loop over those with energy

  TriggerTowerCollection::const_iterator TriggerTowerIterator    = TriggerTowerTES->begin(); 
  TriggerTowerCollection::const_iterator TriggerTowerIteratorEnd = TriggerTowerTES->end(); 

  for (; TriggerTowerIterator != TriggerTowerIteratorEnd; ++TriggerTowerIterator) 
    {
      int EmEnergy = (*TriggerTowerIterator)->emLUT()[(*TriggerTowerIterator)->emPeak()];

      if (EmEnergy>0) 
	{
	  m_h_CaloTT_emLUT_eta-> Fill((*TriggerTowerIterator)->eta(),1);
	  m_h_CaloTT_emLUT_phi-> Fill((*TriggerTowerIterator)->phi(),1);
	  m_h_CaloTT_emLUT->Fill(EmEnergy,1);
	}
	   
      if (EmEnergy>m_CaloTT_HitMap_Thresh0)
	{
	  m_h_CaloTT_HitMap_emLUT_Thresh0->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (EmEnergy>m_CaloTT_HitMap_Thresh1)
	{
	  m_h_CaloTT_HitMap_emLUT_Thresh1->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (EmEnergy>m_CaloTT_HitMap_Thresh2)
	{
	  m_h_CaloTT_HitMap_emLUT_Thresh2->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      
      int HadEnergy = (*TriggerTowerIterator)->hadLUT()[(*TriggerTowerIterator)->hadPeak()];

      if (HadEnergy>0) 
	{
	  m_h_CaloTT_hadLUT_eta-> Fill((*TriggerTowerIterator)->eta(),1);
	  m_h_CaloTT_hadLUT_phi-> Fill((*TriggerTowerIterator)->phi(),1);
	  m_h_CaloTT_hadLUT->Fill(HadEnergy,1);
	}

	     
     if (HadEnergy>m_CaloTT_HitMap_Thresh0)
	{
	  m_h_CaloTT_HitMap_hadLUT_Thresh0->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (HadEnergy>m_CaloTT_HitMap_Thresh1)
	{
	  m_h_CaloTT_HitMap_hadLUT_Thresh1->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (HadEnergy>m_CaloTT_HitMap_Thresh2)
	{
	  m_h_CaloTT_HitMap_hadLUT_Thresh2->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}

    }
  return StatusCode( StatusCode::SUCCESS );

}
/*---------------------------------------------------------*/
StatusCode CaloTTMon::procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  if( isEndOfEventsBlock || isEndOfLumiBlock ) 
    {  
    }
	
  if( isEndOfRun ) 
    {   
    }
  
  return StatusCode( StatusCode::SUCCESS );
}


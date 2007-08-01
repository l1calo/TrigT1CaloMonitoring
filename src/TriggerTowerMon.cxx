// ********************************************************************
//
// NAME:     TrigT1CaloMonTool.cxx
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

#include "TrigT1CaloMonitoring/TriggerTowerMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "TrigT1Calo/TriggerTowerCollection.h"
#include "TrigT1Calo/TrigT1CaloDict.h"
#include "TrigT1Calo/TriggerTower_ClassDEF.h"

#include "GaudiKernel/IService.h"
#include "GaudiKernel/IToolSvc.h"


/*---------------------------------------------------------*/
TriggerTowerMon::TriggerTowerMon(const std::string & type, const std::string & name,
					 const IInterface* parent)
  : ManagedMonitorToolBase ( type, name, parent )
/*---------------------------------------------------------*/
{
  declareProperty("BS_TriggerTowerContainer",  m_TriggerTowerContainerName = "LVL1TriggerTowers");
  declareProperty("LUTHitMap_Thresh0",  m_TT_HitMap_Thresh0 = 1);
  declareProperty("LUTHitMap_Thresh1",  m_TT_HitMap_Thresh1 = 3);
  declareProperty("LUTHitMap_Thresh2",  m_TT_HitMap_Thresh2 = 7);
  declareProperty("ADCHitMap_Thresh",  m_TT_ADC_HitMap_Thresh = 15);
  declareProperty("DistPerChannel", m_TT_DistPerChannel=1);
  declareProperty("DistPerChannelAndTimeSlice", m_TT_DistPerChannelAndTimeSlice=0);
 
  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/CMM") ;
  declareProperty( "DataType", m_DataType="") ;

  m_SliceNo=5;
}

/*---------------------------------------------------------*/
TriggerTowerMon::~TriggerTowerMon()
/*---------------------------------------------------------*/
{
}


/*---------------------------------------------------------*/
StatusCode TriggerTowerMon::bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "in TriggerTowerMon::bookHistograms" << endreq;

  /** get a handle of StoreGate for access to the Event Store */
  StatusCode sc = service("StoreGateSvc", m_storeGate);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR
	   << "Unable to retrieve pointer to StoreGateSvc"
	   << endreq;
      return sc;
    }

  // Get a pointer to DetectorStore services
  sc = service("DetectorStore", m_detStore);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR << "Cannot access DetectorStore" << endreq;
      return StatusCode::FAILURE;;
    }

  // Retrieve the CaloIdManager from the detector store
  sc = m_detStore->retrieve(m_caloMgr);
  if (sc.isFailure()) 
    {
      log << MSG::ERROR << "Unable to retrieve CaloIdManager from DetectorStore" << endreq;
      return StatusCode::FAILURE;
    }
  
  // Use the CaloIdManager to get a pointer to an instance of the CaloLVL1_ID helper
  m_lvl1Helper = m_caloMgr->getLVL1_ID();
  if (!m_lvl1Helper) 
    {
      log << MSG::ERROR << "Could not access CaloLVL1_ID helper" << endreq;
      return StatusCode::FAILURE;
    }

  IToolSvc* toolSvc;
  StatusCode status = service( "ToolSvc",toolSvc  );

  if(status.isSuccess()) 
    {
      IAlgTool *algtool;
      sc = toolSvc->retrieveTool("L1CaloTTIdTools", algtool);
      log<<MSG::DEBUG<<"L1CaloTTIdTools retrieved"<<endreq;
      if (sc!=StatusCode::SUCCESS) 
	{
	  log << MSG::ERROR << " Cannot get L1CaloTTIdTools !" << endreq;
	  return sc;
	}
      m_l1CaloTTIdTools = dynamic_cast<L1CaloTTIdTools*> (algtool);
    } 
  else 
    {
      return StatusCode::FAILURE;
    }
  
  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  MonGroup TT_EmADCPeak( this, (m_PathInRootFile+"/EmADCPeak").c_str(), expert, eventsBlock );
  HistoBooker* EmADCPeak_Booker = new HistoBooker(&TT_EmADCPeak, &log, m_DataType);

  MonGroup TT_HadADCPeak( this, (m_PathInRootFile+"/HadADCPeak").c_str(), expert, eventsBlock );
  HistoBooker* HadADCPeak_Booker = new HistoBooker(&TT_HadADCPeak, &log, m_DataType);

  MonGroup TT_EmLUTPeak( this, (m_PathInRootFile+"/EmLUTPeak").c_str(), expert, eventsBlock );
  HistoBooker* EmLUTPeak_Booker = new HistoBooker(&TT_EmLUTPeak, &log, m_DataType);

  MonGroup TT_HadLUTPeak( this, (m_PathInRootFile+"/HadLUTPeak").c_str(), expert, eventsBlock );
  HistoBooker* HadLUTPeak_Booker = new HistoBooker(&TT_HadLUTPeak, &log, m_DataType);


  MonGroup TT_HitMaps( this, (m_PathInRootFile+"/LUTHitMaps").c_str(), shift, eventsBlock );
  HistoBooker* HitMaps_Booker = new HistoBooker(&TT_HitMaps, &log, m_DataType);

  MonGroup TT_ADC( this, (m_PathInRootFile+"/ADCTimeSlices").c_str(), shift, eventsBlock );
  HistoBooker* ADCTimeSlice_Booker = new HistoBooker(&TT_ADC, &log, m_DataType);

  MonGroup TT_LUTPeakDist( this, (m_PathInRootFile+"/LUTPeakDistribution").c_str(), shift, eventsBlock );
  HistoBooker* LUTPeakDistribution_Booker = new HistoBooker(&TT_LUTPeakDist, &log, m_DataType);

  MonGroup TT_TriggeredSlice( this, (m_PathInRootFile).c_str(), shift, eventsBlock );
  HistoBooker* TriggeredSlice_Booker = new HistoBooker(&TT_TriggeredSlice, &log, m_DataType);

  if( isNewEventsBlock || isNewLumiBlock ) 
    {	
      Helper* Help = new Helper();

      //---------------------------- Energy distributions (ADC and LUT) per Channel -----------------------------
      std::vector<Identifier>::const_iterator tower_it = m_lvl1Helper->tower_begin();
      std::string name,title;
      std::stringstream buffer, etabuffer,phibuffer;
      
      if (m_TT_DistPerChannel==1)
	{
	  for(;tower_it!=m_lvl1Helper->tower_end();++tower_it) 
	    {
	      Identifier towerId = (*tower_it);
	      buffer.str("");
	      etabuffer.str("");
	      phibuffer.str("");
	      buffer<<towerId;
	      etabuffer<<m_l1CaloTTIdTools->IDeta(towerId);
	      phibuffer<<m_l1CaloTTIdTools->IDphi(towerId);
	      
	      if (m_lvl1Helper->sampling(towerId)==0) //EM TT
		{
		  name = "emADCTT_" + buffer.str();
		  title = "TT EM ADC Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_EmADCPeak[towerId]=EmADCPeak_Booker->book1F(name,title,256,-0.5,255.5,"em ADC values");
		  
		  
		  name = "emLUTTT_" + buffer.str();
		  title = "TT EM LUT Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_EmLUTPeak[towerId]=EmLUTPeak_Booker->book1F(name,title,256,-0.5,255.5,"em LUT [GeV]");
		}
	      
	      if (m_lvl1Helper->sampling(towerId)==1) //Had TT
		{
		  name = "hadADCTT_" + buffer.str();
		  title = "TT HAD ADC Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HadADCPeak[towerId]=HadADCPeak_Booker->book1F(name,title,256,-0.5,255.5,"had ADC values");
		  
		  name = "hadLUTTT_" + buffer.str();
		  title = "TT EM LUT Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HadLUTPeak[towerId]=HadLUTPeak_Booker->book1F(name,title,256,-0.5,255.5,"had LUT [GeV]");
		}
	    }
	}

      //---------------------------- ADC Hitmaps per TimeSlice -----------------------------
      for(int i=0;i<m_SliceNo;i++) 
	{
	  buffer.str("");
	  buffer<<i;
	  etabuffer.str("");
	  etabuffer<<m_TT_ADC_HitMap_Thresh;
	  
	  name = "emADCHitMap_" + buffer.str();
	  title = "#eta - #phi Map of EM ADC > " + etabuffer.str() +" for timeslice " + buffer.str();
	  m_h_TT_HitMap_emADC[i]=ADCTimeSlice_Booker->book2F(name,title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	  m_h_TT_HitMap_emADC[i]->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
	  
	  name = "hadADCHitMap_" + buffer.str();
	  title = "#eta - #phi Map of HAD ADC > " + etabuffer.str() +" for timeslice " + buffer.str();
	  m_h_TT_HitMap_hadADC[i]=ADCTimeSlice_Booker->book2F(name,title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	  m_h_TT_HitMap_hadADC[i]->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
	}
      
      //---------------------------- LUT Hitmaps per threshold -----------------------------
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh0;
      m_h_TT_HitMap_emLUT_Thresh0=HitMaps_Booker->book2F("TTEmLUT"+buffer.str(),"#eta - #phi Map of EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_emLUT_Thresh0->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh1;
      m_h_TT_HitMap_emLUT_Thresh1=HitMaps_Booker->book2F("TTEmLUT"+buffer.str(),"#eta - #phi Map of EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_emLUT_Thresh1->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());  
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh2;
      m_h_TT_HitMap_emLUT_Thresh2=HitMaps_Booker->book2F("TTEmLUT"+buffer.str(),"#eta - #phi Map of EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_emLUT_Thresh2->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning()); 
      
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh0;
      m_h_TT_HitMap_hadLUT_Thresh0=HitMaps_Booker->book2F("TTHadLUT"+buffer.str(),"#eta - #phi Map of Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_hadLUT_Thresh0->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());   
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh1;
      m_h_TT_HitMap_hadLUT_Thresh1=HitMaps_Booker->book2F("TTHadLUT"+buffer.str(),"#eta - #phi Map of Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_hadLUT_Thresh1->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());   
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh2;
      m_h_TT_HitMap_hadLUT_Thresh2=HitMaps_Booker->book2F("TTHadLUT"+buffer.str(),"#eta - #phi Map of Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_hadLUT_Thresh2->SetBins(66,Help->TTEtaBinning(),64,Help->TTPhiBinning());  
      
      //---------------------------- distribution of LUT peak per detector region -----------------------------
      m_h_TT_emLUT=LUTPeakDistribution_Booker->book1F("TTEMLUT_Dist","TT EM LUT Distribution of Peak",256,-0.5,255.5,"em LUT Peak [GeV]");
      m_h_TT_emLUT_eta=LUTPeakDistribution_Booker->book1F("TTEMLUT_eta","TT EM LUT Distribution of Peak per #eta",21,-0.5,255.5,"#eta");
      m_h_TT_emLUT_eta->SetBins(66,Help->TTEtaBinning());
      
      m_h_TT_emLUT_phi=LUTPeakDistribution_Booker->book1F("TTEMLUT_phi","TT EM LUT Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
      m_h_TT_emLUT_phi->SetBins(64,Help->TTPhiBinning());  
      
      
      m_h_TT_hadLUT=LUTPeakDistribution_Booker->book1F("TTHADLUT_Dist","TT HAD LUT Distribution of Peak",256,-0.5,255.5,"had LUT Peak [GeV]"); 
      m_h_TT_hadLUT_eta=LUTPeakDistribution_Booker->book1F("TTHADLUT_eta","TT HAD LUT Distribution of Peak per #eta",256,-0.5,255.5,"#eta");
      m_h_TT_hadLUT_eta->SetBins(66,Help->TTEtaBinning());
      m_h_TT_hadLUT_phi=LUTPeakDistribution_Booker->book1F("TTHADLUT_phi","TT HAD LUT Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
      m_h_TT_hadLUT_phi->SetBins(64,Help->TTPhiBinning());  
      
      //---------------------------- number of triggered slice -----------------------------
      m_h_TT_triggeredSlice_em=TriggeredSlice_Booker->book1F("TT_EMTriggeredSlice","Number of the EM Triggered Slice",7,-0.5,6.5,"#Slice");
      m_h_TT_triggeredSlice_had=TriggeredSlice_Booker->book1F("TT_HADTriggeredSlice","Number of the HAD Triggered Slice",7,-0.5,6.5,"#Slice");
      


    }

  if( isNewRun ) { }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode TriggerTowerMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  
  log << MSG::DEBUG << "in fillHistograms()" << endreq;
  
  //Retrieve TriggerTowers from SG
  const TriggerTowerCollection* TriggerTowerTES = 0; 
  StatusCode sc = m_storeGate->retrieve(TriggerTowerTES, m_TriggerTowerContainerName); 
  if( (sc==StatusCode::FAILURE) ) 
    {
      log << MSG::INFO << "No TriggerTower found in TES at "<< m_TriggerTowerContainerName<< endreq ;
      return StatusCode::SUCCESS;
    }
  // =============================================================================================
  // ================= Container: TriggerTower ===================================================
  // =============================================================================================
  TriggerTowerCollection::const_iterator TriggerTowerIterator    = TriggerTowerTES->begin(); 
  TriggerTowerCollection::const_iterator TriggerTowerIteratorEnd = TriggerTowerTES->end(); 
  
  for (; TriggerTowerIterator != TriggerTowerIteratorEnd; ++TriggerTowerIterator) 
    {
      Identifier EmTowerId,HadTowerId;
      
      int detside= m_l1CaloTTIdTools->pos_neg_z((*TriggerTowerIterator)->eta());
      int detregion = m_l1CaloTTIdTools->regionIndex((*TriggerTowerIterator)->eta());
      int eta=m_l1CaloTTIdTools->etaIndex((*TriggerTowerIterator)->eta());
      int phi=m_l1CaloTTIdTools->phiIndex((*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi());
      
      //---------------------------- EM Energy -----------------------------
      EmTowerId = m_lvl1Helper->tower_id(detside, 0, detregion,eta,phi );  
      // em ADC Peak per channel
      int EmEnergy = (*TriggerTowerIterator)->emADC()[(*TriggerTowerIterator)->emPeak()];
      if (m_TT_DistPerChannel==1) m_h_TT_EmADCPeak[EmTowerId]->Fill(EmEnergy,1);
      // em LUT Peak per channel
      EmEnergy = (*TriggerTowerIterator)->emLUT()[(*TriggerTowerIterator)->emPeak()];
      if (m_TT_DistPerChannel==1) m_h_TT_EmLUTPeak[EmTowerId]->Fill(EmEnergy,1);

      // em energy distributions per detector region
      if (EmEnergy>0) 
	{
	  m_h_TT_emLUT_eta-> Fill((*TriggerTowerIterator)->eta(),1);
	  m_h_TT_emLUT_phi-> Fill((*TriggerTowerIterator)->phi(),1);
	  m_h_TT_emLUT->Fill(EmEnergy,1);
	}
	
      // EM LUT HitMaps   
      if (EmEnergy>m_TT_HitMap_Thresh0)
	{
	  m_h_TT_HitMap_emLUT_Thresh0->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (EmEnergy>m_TT_HitMap_Thresh1)
	{
	  m_h_TT_HitMap_emLUT_Thresh1->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (EmEnergy>m_TT_HitMap_Thresh2)
	{
	  m_h_TT_HitMap_emLUT_Thresh2->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      
       //---------------------------- EM Energy -----------------------------
      HadTowerId = m_lvl1Helper->tower_id(detside, 1, detregion,eta,phi );  
      // HAD ADC Peak per channel      
      int HadEnergy = (*TriggerTowerIterator)->hadADC()[(*TriggerTowerIterator)->hadPeak()];
      if (m_TT_DistPerChannel==1) m_h_TT_HadADCPeak[HadTowerId]->Fill(HadEnergy,1);
      // had LUT peak per channel
      HadEnergy = (*TriggerTowerIterator)->hadLUT()[(*TriggerTowerIterator)->hadPeak()];
      if (m_TT_DistPerChannel==1) m_h_TT_HadLUTPeak[HadTowerId]->Fill(HadEnergy,1);

      // had energy distribution per detector region
      if (HadEnergy>0) 
	{
	  m_h_TT_hadLUT_eta-> Fill((*TriggerTowerIterator)->eta(),1);
	  m_h_TT_hadLUT_phi-> Fill((*TriggerTowerIterator)->phi(),1);
	  m_h_TT_hadLUT->Fill(HadEnergy,1);
	}
     
      // had LUT HitMaps   
      if (HadEnergy>m_TT_HitMap_Thresh0)
	{
	  m_h_TT_HitMap_hadLUT_Thresh0->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (HadEnergy>m_TT_HitMap_Thresh1)
	{
	  m_h_TT_HitMap_hadLUT_Thresh1->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}
      if (HadEnergy>m_TT_HitMap_Thresh2)
	{
	  m_h_TT_HitMap_hadLUT_Thresh2->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	}

      // ADC HitMaps per timeslice
      for (int i=0; i<m_SliceNo;i++)
	{
	  if (i<( (*TriggerTowerIterator)->emADC()).size())
	    {
	      if (( (*TriggerTowerIterator)->emADC())[i] > m_TT_ADC_HitMap_Thresh)
		{
		  m_h_TT_HitMap_emADC[i] ->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
		}
	    }
	  if (i<( (*TriggerTowerIterator)->hadADC()).size())
	    {
	      if (( (*TriggerTowerIterator)->hadADC())[i] > m_TT_ADC_HitMap_Thresh)
		{
		  m_h_TT_HitMap_hadADC[i] ->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
		}
	    }
	}    
      
      // number of triggered slice
      m_h_TT_triggeredSlice_em->Fill((*TriggerTowerIterator)->emPeak(),1);
      m_h_TT_triggeredSlice_had->Fill((*TriggerTowerIterator)->hadPeak(),1);
    }

  return StatusCode( StatusCode::SUCCESS );

}
/*---------------------------------------------------------*/
StatusCode TriggerTowerMon::procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
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


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
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IToolSvc.h"

#include <TH1D.h>
#include <TH2D.h>
#include "TString.h"

#include "StoreGate/StoreGateSvc.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TrigT1CaloMonitoring/PPrMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"

#include "TrigT1Calo/TriggerTowerCollection.h"
#include "TrigT1Calo/TrigT1CaloDict.h"
#include "TrigT1Calo/TriggerTower_ClassDEF.h"
#include "TrigT1Calo/DataError.h"
#include "Identifier/HWIdentifier.h"


/*---------------------------------------------------------*/
PPrMon::PPrMon(const std::string & type, const std::string & name,
					 const IInterface* parent)
  : ManagedMonitorToolBase ( type, name, parent )
/*---------------------------------------------------------*/
{
  declareProperty("BS_TriggerTowerContainer",  m_TriggerTowerContainerName = "LVL1TriggerTowers");
  declareProperty("LUTHitMap_Thresh0",  m_TT_HitMap_Thresh0 = 1);
  declareProperty("LUTHitMap_Thresh1",  m_TT_HitMap_Thresh1 = 3);
  declareProperty("LUTHitMap_Thresh2",  m_TT_HitMap_Thresh2 = 7);
  declareProperty("ADCPedestal",  m_TT_ADC_Pedestal = 40);
  declareProperty("ADCHitMap_Thresh",  m_TT_ADC_HitMap_Thresh = 15);
  declareProperty("DistPerChannel", m_TT_DistPerChannel=1);
  declareProperty("DistPerChannelAndTimeSlice", m_TT_DistPerChannelAndTimeSlice=0);
  declareProperty( "MaxEnergyRange", m_MaxEnergyRange = 50) ;
  declareProperty( "Offline", m_Offline = 1) ;
  declareProperty( "ADCTimingPerChannel", m_TT_ADCTimingPerChannel=1);
  declareProperty( "HADFADCCut",  m_HADFADCCut=80);
  declareProperty( "EMFADCCut",  m_EMFADCCut=80);

  declareProperty( "PathInRootFile", m_PathInRootFile="Stats/L1Calo/PPr") ;
  declareProperty( "ErrorPathInRootFile", m_ErrorPathInRootFile="Stats/L1Calo/Errors") ;
  declareProperty( "EventPathInRootFile", m_EventPathInRootFile="Stats/L1Calo") ;
  declareProperty( "TypeOfData", m_DataType="") ;

  m_SliceNo=5;
}

/*---------------------------------------------------------*/
PPrMon::~PPrMon()
/*---------------------------------------------------------*/
{
}


/*---------------------------------------------------------*/
StatusCode PPrMon::bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  MsgStream log( msgSvc(), name() );
  log << MSG::DEBUG << "in PPrMon::bookHistograms" << endreq;

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

  
  ISvcLocator* svcLoc = Gaudi::svcLocator( );
  toolSvc = 0; // Pointer to Tool Service
  sc = svcLoc->service( "ToolSvc",toolSvc  );
  if(sc.isSuccess()) 
    {
      sc = toolSvc->retrieveTool("CaloTriggerTowerService",m_ttSvc);
      if(sc.isFailure())
	{
	  log << MSG::ERROR << "Could not retrieve CaloTriggerTowerService Tool" << endreq;
	  return StatusCode::FAILURE;
	}
    } 
  else 
    {
      log << MSG::ERROR << "Could not retrieve ToolSvc" << endreq;
      return StatusCode::FAILURE;
    }


  // Use the CaloIdManager to get a pointer to an instance of the TTOnlineID helper
  m_l1ttonlineHelper = m_caloMgr->getTTOnlineID();
  if (!m_l1ttonlineHelper ) 
    {
      log << MSG::ERROR << "Could not access TTOnlineID helper" << endreq;
      return StatusCode::FAILURE;
    }
  

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  MonGroup TT_EmADCPeak( this, (m_PathInRootFile+"_EmFADCPeak").c_str(), expert, run );
  HistoBooker EmADCPeak_Booker(&TT_EmADCPeak, &log, m_DataType);

  MonGroup TT_HadADCPeak( this, (m_PathInRootFile+"_HadFADCPeak").c_str(), expert, run );
  HistoBooker HadADCPeak_Booker(&TT_HadADCPeak, &log, m_DataType);

  MonGroup TT_HadADCTiming( this, (m_PathInRootFile+"_HadFADCTiming").c_str(), expert, run );
  HistoBooker HadADCTiming_Booker(&TT_HadADCTiming, &log, m_DataType);

  MonGroup TT_EmADCTiming( this, (m_PathInRootFile+"_EmFADCTiming").c_str(), expert, run );
  HistoBooker EmADCTiming_Booker(&TT_EmADCTiming, &log, m_DataType);

  MonGroup TT_EmLUTPeak( this, (m_PathInRootFile+"_EmLUTPeak").c_str(), expert, run );
  HistoBooker EmLUTPeak_Booker(&TT_EmLUTPeak, &log, m_DataType);

  MonGroup TT_HadLUTPeak( this, (m_PathInRootFile+"_HadLUTPeak").c_str(), expert, run );
  HistoBooker HadLUTPeak_Booker(&TT_HadLUTPeak, &log, m_DataType);

  MonGroup TT_HitMaps( this, (m_PathInRootFile+"_LUTHitMaps").c_str(), shift, run );
  HistoBooker HitMaps_Booker(&TT_HitMaps, &log, m_DataType);

  MonGroup TT_ADC( this, (m_PathInRootFile+"_FADCHitMaps").c_str(), shift, run );
  HistoBooker ADCTimeSlice_Booker(&TT_ADC, &log, m_DataType);

  MonGroup TT_LUTPeakDist( this, (m_PathInRootFile+"_LUTPeakDistribution").c_str(), shift, run );
  HistoBooker LUTPeakDistribution_Booker(&TT_LUTPeakDist, &log, m_DataType);

  MonGroup TT_Error( this, (m_ErrorPathInRootFile).c_str(), shift, run );
  HistoBooker Error_Booker(&TT_Error, &log, "");

  MonGroup NoEvents( this, (m_EventPathInRootFile).c_str(), expert, run );
  HistoBooker NoEvent_Booker(&NoEvents, &log, "");

  if ( isNewEventsBlock|| isNewLumiBlock) { }

  if( isNewRun )

    //if( isNewEventsBlock || isNewLumiBlock ) 
    {	
      Helper Help;
      m_NoEvents=0;
  

      //---------------------------- Energy distributions (ADC and LUT) per Channel -----------------------------
      std::vector<Identifier>::const_iterator tower_it = m_lvl1Helper->tower_begin();
      std::string name,title;
      std::stringstream buffer, etabuffer,phibuffer;
      
      log << MSG::DEBUG << "before book::energy dists" << endreq;
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
		  name = "emFADCTT_" + buffer.str();
		  title = "TT EM FADC Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_EmADCPeak[towerId]=EmADCPeak_Booker.book1F(name,title,256,-0.5,255.5,"em FADC values");
		  
		  
		  name = "emLUTTT_" + buffer.str();
		  title = "TT EM LUT Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_EmLUTPeak[towerId]=EmLUTPeak_Booker.book1F(name,title,256,-0.5,255.5,"em LUT [GeV]");
		}
	      
	      if (m_lvl1Helper->sampling(towerId)==1) //Had TT
		{
		  name = "hadFADCTT_" + buffer.str();
		  title = "TT HAD FADC Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HadADCPeak[towerId]=HadADCPeak_Booker.book1F(name,title,256,-0.5,255.5,"had FADC values");
		  
		  name = "hadLUTTT_" + buffer.str();
		  title = "TT EM LUT Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HadLUTPeak[towerId]=HadLUTPeak_Booker.book1F(name,title,256,-0.5,255.5,"had LUT [GeV]");
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
	  
	  name = "emFADCHitMap_" + buffer.str();
	  title = "#eta - #phi Map of EM FADC > " + etabuffer.str() +" for timeslice " + buffer.str();
	  m_h_TT_HitMap_emADC[i]=ADCTimeSlice_Booker.book2F(name,title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	  m_h_TT_HitMap_emADC[i]->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 
	  
	  name = "hadFADCHitMap_" + buffer.str();
	  title = "#eta - #phi Map of HAD FADC > " + etabuffer.str() +" for timeslice " + buffer.str();
	  m_h_TT_HitMap_hadADC[i]=ADCTimeSlice_Booker.book2F(name,title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	  m_h_TT_HitMap_hadADC[i]->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 
	}




      //---------------------------- Timing of FADC Signal -----------------------------
     if (m_TT_ADCTimingPerChannel==1)
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
	      
	      if (m_lvl1Helper->sampling(towerId)==1 )
		{
		  name = "hadFADCTiming_" + buffer.str();
		  title = "TT had FADC Timing for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HitMap_hadADCChannel_timing[towerId]= HadADCTiming_Booker.bookProfile(name, title,5,-0.5,4.5, "TimeSlice No", "had FADC counts");
		}

	      if (m_lvl1Helper->sampling(towerId)==0) 
		{
		  name = "emFADCTiming_" + buffer.str();
		  title = "TT em FADC Timing for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HitMap_emADCChannel_timing[towerId]= EmADCTiming_Booker.bookProfile(name, title,5,-0.5,4.5, "TimeSlice No", "had FADC counts");
		}
	    }
	}
      
     m_h_TT_ADC_hadTiming_signal= ADCTimeSlice_Booker.bookProfile2D("ADC_hadTiming_signal","Average Maximum TimeSlice for had Signal (TS:1-5)",100,-4.9,4.9, 64,0,2*M_PI,"#eta", "#phi");
      m_h_TT_ADC_hadTiming_signal->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 

      m_h_TT_ADC_emTiming_signal=ADCTimeSlice_Booker.bookProfile2D("ADC_emTiming_signal","em Timing",100,-4.9,4.9, 64,0,2*M_PI, "#eta", "#phi");
      m_h_TT_ADC_emTiming_signal->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 


      //---------------------------- LUT Hitmaps per threshold -----------------------------
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh0;
      m_h_TT_HitMap_emLUT_Thresh0=HitMaps_Booker.book2F("TTEmLUT"+buffer.str(),"#eta - #phi Map of EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_emLUT_Thresh0->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh1;
      m_h_TT_HitMap_emLUT_Thresh1=HitMaps_Booker.book2F("TTEmLUT"+buffer.str(),"#eta - #phi Map of EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_emLUT_Thresh1->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());  
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh2;
      m_h_TT_HitMap_emLUT_Thresh2=HitMaps_Booker.book2F("TTEmLUT"+buffer.str(),"#eta - #phi Map of EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_emLUT_Thresh2->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 
      
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh0;
      m_h_TT_HitMap_hadLUT_Thresh0=HitMaps_Booker.book2F("TTHadLUT"+buffer.str(),"#eta - #phi Map of Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_hadLUT_Thresh0->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());   
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh1;
      m_h_TT_HitMap_hadLUT_Thresh1=HitMaps_Booker.book2F("TTHadLUT"+buffer.str(),"#eta - #phi Map of Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_hadLUT_Thresh1->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());   
      
      buffer.str("");
      buffer<<m_TT_HitMap_Thresh2;
      m_h_TT_HitMap_hadLUT_Thresh2=HitMaps_Booker.book2F("TTHadLUT"+buffer.str(),"#eta - #phi Map of Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_HitMap_hadLUT_Thresh2->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());  
      
      //---------------------------- distribution of LUT peak per detector region -----------------------------
      m_h_TT_emLUT=LUTPeakDistribution_Booker.book1F("emLUT_peak","EM LUT: Distribution of Peak",m_MaxEnergyRange,0,m_MaxEnergyRange,"em LUT Peak [GeV]");
      m_h_TT_emLUT_eta=LUTPeakDistribution_Booker.book1F("emLUT_eta","EM LUT: Distribution of Peak per #eta",21,-0.5,255.5,"#eta");
      m_h_TT_emLUT_eta->SetBins(66,Help.TTEtaBinning());
      
      m_h_TT_emLUT_phi=LUTPeakDistribution_Booker.book1F("emLUT_phi","EM LUT: Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
      m_h_TT_emLUT_phi->SetBins(64,Help.TTPhiBinning());  
      
      
      m_h_TT_hadLUT=LUTPeakDistribution_Booker.book1F("hadLUT_Dist","HAD LUT: Distribution of Peak",m_MaxEnergyRange,0,m_MaxEnergyRange,"had LUT Peak [GeV]"); 
      m_h_TT_hadLUT_eta=LUTPeakDistribution_Booker.book1F("hadLUT_eta","HAD LUT: Distribution of Peak per #eta",256,-0.5,255.5,"#eta");
      m_h_TT_hadLUT_eta->SetBins(66,Help.TTEtaBinning());
      m_h_TT_hadLUT_phi=LUTPeakDistribution_Booker.book1F("hadLUT_phi","HAD LUT: Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
      m_h_TT_hadLUT_phi->SetBins(64,Help.TTPhiBinning());  
      
      //---------------------------- SubStatus Word errors -----------------------------
      m_h_TT_emerror=Error_Booker.book1F("TT_emerror","EM TT SubStatus Word errors",17,0.5,17.5,"");
      
      m_h_TT_emerror->GetXaxis()->SetBinLabel(1, "ChannelDisabled");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(2, "MCMAbsent");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(3, "Timeout");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(4, "ASICFull");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(5, "EventMismatch");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(6, "BunchMismatch");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(7, "FIFOCorrupt");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(8, "PinParity");
      
      m_h_TT_emerror->GetXaxis()->SetBinLabel(10, "GLinkParity");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(11, "GLinkProtocol");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(12, "BCNMismatch");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(13, "FIFOOverflow");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(14, "ModuleError");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(15, "GLinkDown");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(16, "GLinkTimeout");
      m_h_TT_emerror->GetXaxis()->SetBinLabel(17, "FailingBCN");

      m_h_TT_haderror=Error_Booker.book1F("TT_haderror","HAD TT SubStatus Word errors",17,0.5,17.5,"");
      
      m_h_TT_haderror->GetXaxis()->SetBinLabel(1, "ChannelDisabled");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(2, "MCMAbsent");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(3, "Timeout");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(4, "ASICFull");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(5, "EventMismatch");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(6, "BunchMismatch");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(7, "FIFOCorrupt");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(8, "PinParity");
      
      m_h_TT_haderror->GetXaxis()->SetBinLabel(10, "GLinkParity");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(11, "GLinkProtocol");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(12, "BCNMismatch");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(13, "FIFOOverflow");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(14, "ModuleError");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(15, "GLinkDown");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(16, "GLinkTimeout");
      m_h_TT_haderror->GetXaxis()->SetBinLabel(17, "FailingBCN");
      
      m_h_TT_error_Crate_03=Error_Booker.book2F("TT_error_Crate_0-3","TT SubStatus Word errors in crates 0-3",17,0.5,17.5,71,0.5,71.5,"","");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(1, "ChannelDisabled");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(2, "MCMAbsent");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(3, "Timeout");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(4, "ASICFull");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(5, "EventMismatch");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(6, "BunchMismatch");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(7, "FIFOCorrupt");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(8, "PinParity");
      
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(10, "GLinkParity");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(11, "GLinkProtocol");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(12, "BCNMismatch");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(13, "FIFOOverflow");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(14, "ModuleError");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(15, "GLinkDown");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(16, "GLinkTimeout");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(17, "FailingBCN");

      m_h_TT_error_Crate_47=Error_Booker.book2F("TT_error_Crate_4-7","TT SubStatus Word errors in crates 4-7",17,0.5,17.5,71,0.5,71.5,"","");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(1, "ChannelDisabled");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(2, "MCMAbsent");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(3, "Timeout");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(4, "ASICFull");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(5, "EventMismatch");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(6, "BunchMismatch");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(7, "FIFOCorrupt");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(8, "PinParity");
      
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(10, "GLinkParity");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(11, "GLinkProtocol");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(12, "BCNMismatch");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(13, "FIFOOverflow");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(14, "ModuleError");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(15, "GLinkDown");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(16, "GLinkTimeout");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(17, "FailingBCN");

    
      for (int i=5; i<21; i++)
	{
	  buffer.str("");
	  buffer<<i;
	  
	  name = "PPM " + buffer.str();
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel((i-4)+0*18, name.c_str());
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel((i-4)+1*18, name.c_str());
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel((i-4)+2*18, name.c_str());
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel((i-4)+3*18, name.c_str());


	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel((i-4)+(4-4)*18, name.c_str());
	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel((i-4)+(5-4)*18, name.c_str());
	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel((i-4)+(6-4)*18, name.c_str());
	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel((i-4)+(7-4)*18, name.c_str());
	  i++;
	}
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel(17+0*18, "Crate 0");
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel(17+1*18, "Crate 1");
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel(17+2*18, "Crate 2");
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel(17+3*18, "Crate 3");

	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel(17+(4-4)*18, "Crate 4");
	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel(17+(5-4)*18, "Crate 5");
	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel(17+(6-4)*18, "Crate 6");
	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel(17+(7-4)*18, "Crate 7");



      m_h_TT_em_GLinkDown=Error_Booker.book2F("TT_em_GLinkDown","EM TT GLink Down per #eta-#phi",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_em_GLinkDown->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 
      
      m_h_TT_em_GLinkTimeout=Error_Booker.book2F("TT_em_GLinkTimeout","EM TT GLink Timeout per #eta-#phi",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_em_GLinkTimeout->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 


      m_h_TT_had_GLinkDown=Error_Booker.book2F("TT_had_GLinkDown","HAD TT GLink Down per #eta-#phi",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_had_GLinkDown->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 

      m_h_TT_had_GLinkTimeout=Error_Booker.book2F("TT_had_GLinkTimeout","HAD TT GLink Timeout per #eta-#phi",100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
      m_h_TT_had_GLinkTimeout->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 

      //---------------------------- number of triggered slice -----------------------------
      m_h_TT_triggeredSlice_em=ADCTimeSlice_Booker.book1F("TT_EMTriggeredSlice","Number of the EM Triggered Slice",7,-0.5,6.5,"#Slice");
      m_h_TT_triggeredSlice_had=ADCTimeSlice_Booker.book1F("TT_HADTriggeredSlice","Number of the HAD Triggered Slice",7,-0.5,6.5,"#Slice");
      
   
      //----------------------------- number of events ----------------------------------
      m_h_NumberEvents= NoEvent_Booker.book1F("NumberEvents","Number of processed events",1,0.5,1.5,"");
      m_h_NumberEvents->GetXaxis()->SetBinLabel(1,"Number of Events");

    }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPrMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  m_NoEvents++;

  log << MSG::DEBUG << "in fillHistograms()" << endreq;
  
  //Retrieve TriggerTowers from SG
  const TriggerTowerCollection* TriggerTowerTES = 0; 
  StatusCode sc = m_storeGate->retrieve(TriggerTowerTES, m_TriggerTowerContainerName); 
  if( (sc==StatusCode::FAILURE) ) 
    {
      log << MSG::INFO << "No TriggerTower found in TES at "<< m_TriggerTowerContainerName<< endreq ;
      return StatusCode::SUCCESS;
    }

  m_h_NumberEvents->Fill(1,1);
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
      int EmEnergy = (*TriggerTowerIterator)->emADC()[(*TriggerTowerIterator)->emADCPeak()];
      //if (m_TT_DistPerChannel==1) m_h_TT_EmADCPeak[EmTowerId]->Fill(EmEnergy,1);
      // em LUT Peak per channel
      EmEnergy = (*TriggerTowerIterator)->emLUT()[(*TriggerTowerIterator)->emPeak()];
      //if (m_TT_DistPerChannel==1) m_h_TT_EmLUTPeak[EmTowerId]->Fill(EmEnergy,1);

      // em energy distributions per detector region
      if (EmEnergy>0) 
	{
	  m_h_TT_emLUT_eta-> Fill((*TriggerTowerIterator)->eta(),1);
	  m_h_TT_emLUT_phi-> Fill((*TriggerTowerIterator)->phi(),1);
	  m_h_TT_emLUT->Fill(EmEnergy,1);
	}
	
       //---------------------------- EM LUT HitMaps -----------------------------
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
      
       //---------------------------- HAD Energy -----------------------------
      HadTowerId = m_lvl1Helper->tower_id(detside, 1, detregion,eta,phi );  
      // HAD ADC Peak per channel      
      int HadEnergy = (*TriggerTowerIterator)->hadADC()[(*TriggerTowerIterator)->hadADCPeak()];
      //if (m_TT_DistPerChannel==1) m_h_TT_HadADCPeak[HadTowerId]->Fill(HadEnergy,1);
      // had LUT peak per channel
      HadEnergy = (*TriggerTowerIterator)->hadLUT()[(*TriggerTowerIterator)->hadPeak()];
      //if (m_TT_DistPerChannel==1) m_h_TT_HadLUTPeak[HadTowerId]->Fill(HadEnergy,1);


      // had energy distribution per detector region
      if (HadEnergy>0) 
	{
	  m_h_TT_hadLUT_eta-> Fill((*TriggerTowerIterator)->eta(),1);
	  m_h_TT_hadLUT_phi-> Fill((*TriggerTowerIterator)->phi(),1);
	  m_h_TT_hadLUT->Fill(HadEnergy,1);
	}
     
       //---------------------------- had LUT HitMaps -----------------------------
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

       //---------------------------- ADC HitMaps per timeslice -----------------------------
      for (int i=0; i<m_SliceNo;i++)
	{
	  if (i<static_cast<int>(( (*TriggerTowerIterator)->emADC()).size()))
	    {
	      if (( (*TriggerTowerIterator)->emADC())[i] > m_TT_ADC_HitMap_Thresh)
		{
		  m_h_TT_HitMap_emADC[i] ->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
		}
	    }

	  if (i<static_cast<int>(( (*TriggerTowerIterator)->hadADC()).size()))
	    {
	      if (( (*TriggerTowerIterator)->hadADC())[i] > m_TT_ADC_HitMap_Thresh)
		{
		  m_h_TT_HitMap_hadADC[i] ->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);	  
		}
	    }
	}    

        //---------------------------- Timing of FADC Signal -----------------------------
      int hadFADCSum=0;     
      int emFADCSum=0;     
      double max;

      emFADCSum=FADCSum((*TriggerTowerIterator)->emADC());
      hadFADCSum=FADCSum((*TriggerTowerIterator)->hadADC());

      if (m_TT_ADCTimingPerChannel==1)
	{
	  for (int i=0; i<m_SliceNo;i++)
	    {
	      if (i<static_cast<int>(( (*TriggerTowerIterator)->emADC()).size()))
		{ 
		  if (emFADCSum>m_EMFADCCut) 
		    {
		      m_h_TT_HitMap_emADCChannel_timing[EmTowerId]->Fill(i,(*TriggerTowerIterator)->emADC()[i]);
		    }
		}
	      
	      if (i<static_cast<int>(( (*TriggerTowerIterator)->hadADC()).size()))
		{
		  if (hadFADCSum>m_HADFADCCut) 
		    {
		      m_h_TT_HitMap_hadADCChannel_timing[HadTowerId]->Fill(i,(*TriggerTowerIterator)->hadADC()[i]);
		    }
		}
	    }
	}

      if (emFADCSum>m_EMFADCCut)
	{
	  max = recTime((*TriggerTowerIterator)->emADC())+1;
	  //log << MSG::INFO << "TimeSlice of Maximum "<< max<< endreq ;
	  m_h_TT_ADC_emTiming_signal->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),max);
	}

      if (hadFADCSum>m_HADFADCCut)
	{
	  max = recTime((*TriggerTowerIterator)->hadADC())+1;
	  //log << MSG::INFO << "TimeSlice of Maximum "<< max<< endreq ;
	  m_h_TT_ADC_hadTiming_signal->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),max);
	}



      

      //---------------------------- SubStatus Word errors -----------------------------
      //----------------------------- em ---------------------------------------

      LVL1::DataError emerr((*TriggerTowerIterator)-> emError());

      int crate=0;
      int module=0;
      try 
	{
	  HWIdentifier ttOnlId = m_ttSvc->createTTChannelID(EmTowerId);
	  crate     = m_l1ttonlineHelper->crate(ttOnlId);
	  module    = m_l1ttonlineHelper->module(ttOnlId);
	  //log << MSG::DEBUG << "em PPM crate: " << crate<<"  module: "<<module << endreq ;
	} 
      catch(LArID_Exception& except) 
	{
	  log << MSG::ERROR << "LArID_Exception " << (std::string) except << endreq ;
	}
      if (crate>3) log << MSG::INFO << "Wrong em Crate " << crate << endreq ;


      // ChannelDisabled
      m_h_TT_emerror->Fill(1,emerr.get(4));
      // MCMAbsent
      m_h_TT_emerror->Fill(2,emerr.get(5));
      // Timeout
      m_h_TT_emerror->Fill(3,emerr.get(6));
      // ASICFull
      m_h_TT_emerror->Fill(4,emerr.get(7));
      // EventMismatch
      m_h_TT_emerror->Fill(5,emerr.get(8));
      // BunchMismatch
      m_h_TT_emerror->Fill(6,emerr.get(9));
      // FIFOCorrupt
      m_h_TT_emerror->Fill(7,emerr.get(10));
      // PinParity
      m_h_TT_emerror->Fill(8,emerr.get(11));
      		  
      // GLinkParity
      m_h_TT_emerror->Fill(10,emerr.get(16));
      // GLinkProtocol
      m_h_TT_emerror->Fill(11,emerr.get(17));
      // BCNMismatch
      m_h_TT_emerror->Fill(12,emerr.get(18));
      // FIFOOverflow
      m_h_TT_emerror->Fill(13,emerr.get(19));
      // ModuleError
      m_h_TT_emerror->Fill(14,emerr.get(20));
      
      // GLinkDown
      m_h_TT_emerror->Fill(15,emerr.get(22));
      // GLinkTimeout
      m_h_TT_emerror->Fill(16,emerr.get(23));
      // FailingBCN
      m_h_TT_emerror->Fill(17,emerr.get(24));

      //---------------- per crate and module -------------------------  
      // ChannelDisabled
      m_h_TT_error_Crate_03->Fill(1,(module-4)+(crate*18),emerr.get(4));
      // MCMAbsent
      m_h_TT_error_Crate_03->Fill(2,(module-4)+(crate*18),emerr.get(5));
      // Timeout
      m_h_TT_error_Crate_03->Fill(3,(module-4)+(crate*18),emerr.get(6));
      // ASICFull
      m_h_TT_error_Crate_03->Fill(4,(module-4)+(crate*18),emerr.get(7));
      // EventMismatch
      m_h_TT_error_Crate_03->Fill(5,(module-4)+(crate*18),emerr.get(8));
      // BunchMismatch
      m_h_TT_error_Crate_03->Fill(6,(module-4)+(crate*18),emerr.get(9));
      // FIFOCorrupt
      m_h_TT_error_Crate_03->Fill(7,(module-4)+(crate*18),emerr.get(10));
      // PinParity
      m_h_TT_error_Crate_03->Fill(8,(module-4)+(crate*18),emerr.get(11));
      		  
      // GLinkParity
      m_h_TT_error_Crate_03->Fill(10,(module-4)+(crate*18),emerr.get(16));
      // GLinkProtocol
      m_h_TT_error_Crate_03->Fill(11,(module-4)+(crate*18),emerr.get(17));
      // BCNMismatch
      m_h_TT_error_Crate_03->Fill(12,(module-4)+(crate*18),emerr.get(18));
      // FIFOOverflow
      m_h_TT_error_Crate_03->Fill(13,(module-4)+(crate*18),emerr.get(19));
      // ModuleError
      m_h_TT_error_Crate_03->Fill(14,(module-4)+(crate*18),emerr.get(20));
      
      // GLinkDown
      m_h_TT_error_Crate_03->Fill(15,(module-4)+(crate*18),emerr.get(22));
      // GLinkTimeout
      m_h_TT_error_Crate_03->Fill(16,(module-4)+(crate*18),emerr.get(23));
      // FailingBCN
      m_h_TT_error_Crate_03->Fill(17,(module-4)+(crate*18),emerr.get(24));


    
     LVL1::DataError haderr((*TriggerTowerIterator)-> hadError());

      // ChannelDisabled
      m_h_TT_haderror->Fill(1,haderr.get(4));
      // MCMAbsent
      m_h_TT_haderror->Fill(2,haderr.get(5));
      // Timeout
      m_h_TT_haderror->Fill(3,haderr.get(6));
      // ASICFull
      m_h_TT_haderror->Fill(4,haderr.get(7));
      // EventMismatch
      m_h_TT_haderror->Fill(5,haderr.get(8));
      // BunchMismatch
      m_h_TT_haderror->Fill(6,haderr.get(9));
      // FIFOCorrupt
      m_h_TT_haderror->Fill(7,haderr.get(10));
      // PinParity
      m_h_TT_haderror->Fill(8,haderr.get(11));
      		  
      // GLinkParity
      m_h_TT_haderror->Fill(10,haderr.get(16));
      // GLinkProtocol
      m_h_TT_haderror->Fill(11,haderr.get(17));
      // BCNMismatch
      m_h_TT_haderror->Fill(12,haderr.get(18));
      // FIFOOverflow
      m_h_TT_haderror->Fill(13,haderr.get(19));
      // ModuleError
      m_h_TT_haderror->Fill(14,haderr.get(20));
      
      // GLinkDown
      m_h_TT_haderror->Fill(15,haderr.get(22));
      // GLinkTimeout
      m_h_TT_haderror->Fill(16,haderr.get(23));
      // FailingBCN
      if (haderr.get(24)!=0) m_h_TT_haderror->Fill(17,1);

      try 
	{
	  HWIdentifier ttOnlId = m_ttSvc->createTTChannelID(HadTowerId);
	  crate     = m_l1ttonlineHelper->crate(ttOnlId);
	  module    = m_l1ttonlineHelper->module(ttOnlId);
	  //log << MSG::DEBUG << "em PPM crate: " << crate<<"  module: "<<module << endreq ;
	} 
      catch(LArID_Exception& except) 
	{
	  log << MSG::ERROR << "LArID_Exception " << (std::string) except << endreq ;
	}
      if (crate<4) log << MSG::INFO << "Wrong had Crate " << crate << endreq ;

      //---------------- per crate and module -------------------------  m_h_TT_error_Crate_03
      // ChannelDisabled
      m_h_TT_error_Crate_47->Fill(1,(module-4)+((crate-4)*18),haderr.get(4));
      // MCMAbsent
      m_h_TT_error_Crate_47->Fill(2,(module-4)+((crate-4)*18),haderr.get(5));
      // Timeout
      m_h_TT_error_Crate_47->Fill(3,(module-4)+((crate-4)*18),haderr.get(6));
      // ASICFull
      m_h_TT_error_Crate_47->Fill(4,(module-4)+((crate-4)*18),haderr.get(7));
      // EventMismatch
      m_h_TT_error_Crate_47->Fill(5,(module-4)+((crate-4)*18),haderr.get(8));
      // BunchMismatch
      m_h_TT_error_Crate_47->Fill(6,(module-4)+((crate-4)*18),haderr.get(9));
      // FIFOCorrupt
      m_h_TT_error_Crate_47->Fill(7,(module-4)+((crate-4)*18),haderr.get(10));
      // PinParity
      m_h_TT_error_Crate_47->Fill(8,(module-4)+((crate-4)*18),haderr.get(11));
      		  
      // GLinkParity
      m_h_TT_error_Crate_47->Fill(10,(module-4)+((crate-4)*18),haderr.get(16));
      // GLinkProtocol
      m_h_TT_error_Crate_47->Fill(11,(module-4)+((crate-4)*18),haderr.get(17));
      // BCNMismatch
      m_h_TT_error_Crate_47->Fill(12,(module-4)+((crate-4)*18),haderr.get(18));
      // FIFOOverflow
      m_h_TT_error_Crate_47->Fill(13,(module-4)+((crate-4)*18),haderr.get(19));
      // ModuleError
      m_h_TT_error_Crate_47->Fill(14,(module-4)+((crate-4)*18),haderr.get(20));
      
      // GLinkDown
      m_h_TT_error_Crate_47->Fill(15,(module-4)+((crate-4)*18),haderr.get(22));
      // GLinkTimeout
      m_h_TT_error_Crate_47->Fill(16,(module-4)+((crate-4)*18),haderr.get(23));
      // FailingBCN
      m_h_TT_error_Crate_47->Fill(17,(module-4)+((crate-4)*18),haderr.get(24));


      m_h_TT_em_GLinkDown->Fill((*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi(), emerr.get(22)); 
      m_h_TT_em_GLinkTimeout->Fill((*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi(), emerr.get(23)); 
      m_h_TT_had_GLinkDown->Fill((*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi(), haderr.get(22)); 
      m_h_TT_had_GLinkTimeout->Fill((*TriggerTowerIterator)->eta(), (*TriggerTowerIterator)->phi(), haderr.get(23)); 



      // number of triggered slice
      m_h_TT_triggeredSlice_em->Fill((*TriggerTowerIterator)->emADCPeak(),1);
      m_h_TT_triggeredSlice_had->Fill((*TriggerTowerIterator)->hadADCPeak(),1);


    }

  return StatusCode( StatusCode::SUCCESS );

}
/*---------------------------------------------------------*/
StatusCode PPrMon::procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock ) 
    {  
    }
	
  if(m_Offline==1)
    {
      if( isEndOfRun ) 
	{   
	  
	  std::stringstream buffer;
	  buffer.str("");
	  buffer<<m_NoEvents;
	  std::string title;
	  
	  title = m_h_TT_emerror-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_emerror->SetTitle(title.c_str());
	  
	  title = m_h_TT_error_Crate_03-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_error_Crate_03->SetTitle(title.c_str());
	  
	  title = m_h_TT_error_Crate_47-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_error_Crate_47->SetTitle(title.c_str());
	  
	  title = m_h_TT_haderror-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_haderror->SetTitle(title.c_str());
	  
	  title = m_h_TT_em_GLinkDown-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_em_GLinkDown->SetTitle(title.c_str());
	  
	  title = m_h_TT_em_GLinkTimeout-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_em_GLinkTimeout->SetTitle(title.c_str());
	  
	  title = m_h_TT_had_GLinkDown-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_had_GLinkDown->SetTitle(title.c_str());
	  
	  title = m_h_TT_had_GLinkTimeout-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_had_GLinkTimeout->SetTitle(title.c_str());
	}
    }
  return StatusCode( StatusCode::SUCCESS );
}


/*---------------------------------------------------------*/
double PPrMon::recTime(const std::vector<int>& vFAdc) {
/*---------------------------------------------------------*/

  double x[3];
  double y[3];
  double binshift = 0.;
  
  int indmax=0;
  int index=0;
  std::vector<int>::const_iterator it = vFAdc.begin();
  for(;it!=vFAdc.end();++it,++index) {
    if(vFAdc[indmax]<*it) indmax=index;
  }
  
  //cout<<"indmax: "<<indmax<<endl;
  
  double max = 0.;
  
  if(indmax==0) {
    max=0.+binshift;
    //x[0] = 0 + binshift;  y[0] = vFAdc[0];
    //x[1] = 1     + binshift; y[1] = vFAdc[1];
    //x[2] = 2 + 1 +binshift; y[2] = vFAdc[2];
  } else if(indmax==4) {
    max=4.+binshift;
    //x[0] = 2 + binshift;  y[0] = vFAdc[2];
    //x[1] = 3     + binshift; y[1] = vFAdc[3];
    //x[2] = 4 + 1 +binshift; y[2] = vFAdc[4];
  } else {
    //if(indmax!=0 && indmax!=4) {
    
    x[0] = indmax - 1 + binshift;  y[0] = vFAdc[indmax-1];
    x[1] = indmax     + binshift; y[1] = vFAdc[indmax];
    x[2] = indmax + 1 +binshift; y[2] = vFAdc[indmax+1];
    
    double a = ( (x[0]-x[2])*(y[1]-y[2]) - (x[1]-x[2])*(y[0]-y[2]) ) / (
									(x[0]-x[2])*(x[1]*x[1]-x[2]*x[2]) - (x[1]-x[2])*(x[0]*x[0]-x[2]*x[2]) );
    double b = ( (y[1]-y[2])*(x[0]*x[0]-x[2]*x[2]) -
		 (y[0]-y[2])*(x[1]*x[1]-x[2]*x[2]) ) / (
							(x[1]-x[2])*(x[0]*x[0]-x[2]*x[2]) - (x[0]-x[2])*(x[1]*x[1]-x[2]*x[2]) );
    //double c = y[0] - b*x[0] - a*x[0]*x[0];
    max = -b/(2*a);
    
  }
  
  return max;
}



/*---------------------------------------------------------*/
int PPrMon::FADCSum(const std::vector<int>& vFAdc) {
  /*---------------------------------------------------------*/
  
  std::vector<int>::const_iterator it = vFAdc.begin();
  int FADCSum =0;
  int FADCMin=999999;
  for(;it!=vFAdc.end();++it) {
    FADCSum+=*it;
    if(*it<FADCMin)FADCMin=*it;
  }

  FADCSum-=vFAdc.size()*FADCMin; // pseudo pedestal subtraction...
  
  return FADCSum;
}


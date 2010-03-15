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
#include "SGTools/StlVectorClids.h"

#include "TString.h"

#include "StoreGate/StoreGateSvc.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

#include "TrigT1CaloMonitoring/PPrMon.h"
#include "TrigT1CaloMonitoring/MonHelper.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"

#include "TrigT1CaloEvent/TriggerTowerCollection.h"
#include "TrigT1CaloEvent/TriggerTower_ClassDEF.h"
#include "TrigT1CaloUtils/DataError.h"
#include "Identifier/HWIdentifier.h"


/*---------------------------------------------------------*/
PPrMon::PPrMon(const std::string & type, const std::string & name,
					 const IInterface* parent)
  : ManagedMonitorToolBase ( type, name, parent ),
    m_errorTool("TrigT1CaloMonErrorTool")
/*---------------------------------------------------------*/
{
  declareProperty("BS_TriggerTowerContainer",  m_TriggerTowerContainerName = "LVL1TriggerTowers");
  declareProperty("LUTHitMap_Thresh0",  m_TT_HitMap_Thresh0 = 1);
  declareProperty("LUTHitMap_Thresh1",  m_TT_HitMap_Thresh1 = 3);
  declareProperty("LUTHitMap_Thresh2",  m_TT_HitMap_Thresh2 = 7);
  declareProperty("LUTHitMap_ThreshMax",  m_TT_HitMap_ThreshMax = 10);
  declareProperty("LUTHitMap_LumiBlocks",  m_TT_HitMap_LumiBlocks = 10);
  declareProperty("ADCPedestal",  m_TT_ADC_Pedestal = 35);
  declareProperty("ADCHitMap_Thresh",  m_TT_ADC_HitMap_Thresh = 15);
  declareProperty("DistPerChannel", m_TT_DistPerChannel=0);
  declareProperty("DistPerChannelAndTimeSlice", m_TT_DistPerChannelAndTimeSlice=0);
  declareProperty("MaxEnergyRange", m_MaxEnergyRange = 256) ;
  declareProperty("Offline", m_Offline = 1) ;
  declareProperty("ADCTimingPerChannel", m_TT_ADCTimingPerChannel=0);
  // Next two cuts now average per timeslice
  declareProperty("HADFADCCut",  m_HADFADCCut=8);
  declareProperty("EMFADCCut",  m_EMFADCCut=8);

  declareProperty("PathInRootFile", m_PathInRootFile="L1Calo/PPM") ;
  declareProperty("ErrorPathInRootFile", m_ErrorPathInRootFile="L1Calo/PPM/Errors") ;
  declareProperty("EventPathInRootFile", m_EventPathInRootFile="L1Calo/Overview") ;
  declareProperty("TypeOfData", m_DataType="") ;
  declareProperty("OnlineTest", m_onlineTest = false,
                  "Test online code when running offline");

  // Maximum possible number of ADC slices
  m_SliceNo=15;
}

/*---------------------------------------------------------*/
PPrMon::~PPrMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode PPrMon::initialize()
/*---------------------------------------------------------*/
{
  MsgStream log( msgSvc(), name() );

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_errorTool.retrieve();
  if( sc.isFailure() ) {
    log << MSG::ERROR << "Unable to locate Tool TrigT1CaloMonErrorTool"
                      << endreq;
    return sc;
  }
  return StatusCode::SUCCESS;
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

  /*
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
  */
  MonGroup TT_EmADCPeak( this, m_PathInRootFile+"/ADC/Channels/Distributions_Em", expert, run );
  HistoBooker EmADCPeak_Booker(&TT_EmADCPeak, &log, "");

  MonGroup TT_HadADCPeak( this, m_PathInRootFile+"/ADC/Channels/Distributions_Had", expert, run );
  HistoBooker HadADCPeak_Booker(&TT_HadADCPeak, &log, "");

  MonGroup TT_HadADCTiming( this, m_PathInRootFile+"/ADC/Channels/Timing_Had", expert, run );
  HistoBooker HadADCTiming_Booker(&TT_HadADCTiming, &log, "");

  MonGroup TT_EmADCTiming( this, m_PathInRootFile+"/ADC/Channels/Timing_Em", expert, run );
  HistoBooker EmADCTiming_Booker(&TT_EmADCTiming, &log, "");

  MonGroup TT_EmLUTPeak( this, m_PathInRootFile+"/LUT/Channels/Distributions_Em", expert, run );
  HistoBooker EmLUTPeak_Booker(&TT_EmLUTPeak, &log, "");

  MonGroup TT_HadLUTPeak( this, m_PathInRootFile+"/LUT/Channels/Distributions_Had", expert, run );
  HistoBooker HadLUTPeak_Booker(&TT_HadLUTPeak, &log, "");

  MonGroup TT_HitMaps( this, m_PathInRootFile+"/LUT/EtaPhiMaps", shift, run );
  HistoBooker HitMaps_Booker(&TT_HitMaps, &log, "");

  MonGroup TT_ADC( this, m_PathInRootFile+"/ADC/EtaPhiMaps", shift, run );
  HistoBooker ADCHitMaps_Booker(&TT_ADC, &log, "");

  MonGroup TT_ADCSlices( this, m_PathInRootFile+"/ADC/Timeslices", shift, run );
  HistoBooker ADCTimeSlice_Booker(&TT_ADCSlices, &log, "");

  MonGroup TT_LUTPeakDist( this, m_PathInRootFile+"/LUT/Distributions", shift, run );
  HistoBooker LUTPeakDistribution_Booker(&TT_LUTPeakDist, &log, "");

  MonGroup TT_Error( this, m_ErrorPathInRootFile, shift, run );
  HistoBooker Error_Booker(&TT_Error, &log, "");

  MonGroup TT_ErrorDetail( this, m_ErrorPathInRootFile+"/Detail", expert, run );
  HistoBooker ErrorDetail_Booker(&TT_ErrorDetail, &log, "");

  MonGroup NoEvents( this, m_EventPathInRootFile, expert, run );
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
      std::stringstream buffer, etabuffer,phibuffer, thresbuffer;
	
      log << MSG::DEBUG << "before book::energy dists" << endreq;

      if (m_TT_DistPerChannel==1)
	{
	  for (;tower_it!=m_lvl1Helper->tower_end();++tower_it) 
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
		  //name = "emFADCTT_" + buffer.str();
		  name = "ppm_em_1d_tt_adc_Channel" + buffer.str();
		  title = "TT EM FADC Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_EmADCPeak[towerId]=EmADCPeak_Booker.book1F(name,title,256,-0.5,255.5,"em FADC values");
		  
		  
		  //name = "emLUTTT_" + buffer.str();
		  name = "ppm_em_1d_tt_lut_Channel" + buffer.str();
		  title = "TT EM LUT Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_EmLUTPeak[towerId]=EmLUTPeak_Booker.book1F(name,title,256,-0.5,255.5,"em LUT [GeV]");
		}
	    
		
	      
	      if (m_lvl1Helper->sampling(towerId)==1) //Had TT
		{
		  //name = "hadFADCTT_" + buffer.str();
		  name = "ppm_had_1d_tt_adc_Channel" + buffer.str();
		  title = "TT HAD FADC Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HadADCPeak[towerId]=HadADCPeak_Booker.book1F(name,title,256,-0.5,255.5,"had FADC values");
		  
		  //name = "hadLUTTT_" + buffer.str();
		  name = "ppm_had_1d_tt_lut_Channel" + buffer.str();
		  title = "TT EM LUT Distribution of Peak for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HadLUTPeak[towerId]=HadLUTPeak_Booker.book1F(name,title,256,-0.5,255.5,"had LUT [GeV]");
		}
	    }
	  }
	  

      //---------------------------- ADC Hitmaps for Triggered Timeslice -----------------------------
    

          thresbuffer.str("");
	  thresbuffer<<m_TT_ADC_HitMap_Thresh;

	  //title="#eta - #phi Map of EM FADC > "+ thresbuffer.str() + " for timeslice 2";
	  //m_h_TT_HitMap_emADC_00100=ADCTimeSlice_Booker.book2F("emFADCHitMap2",title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	  title="#eta - #phi Map of EM FADC > "+ thresbuffer.str() + " for triggered timeslice";
	  m_h_TT_HitMap_emADC_00100=ADCHitMaps_Booker.book2F("ppm_em_2d_etaPhi_tt_adc_HitMap",title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	  m_h_TT_HitMap_emADC_00100->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());
	  m_h_TT_HitMap_emADC_00100->SetStats(kFALSE);
	  //m_h_TT_HitMap_emADC_00100->SetOption(colz);
	  //title="#eta - #phi Profile Map of EM FADC > "+ thresbuffer.str() + " for timeslice 2";
	  //m_p_TT_HitMap_emADC_00100=ADCTimeSlice_Booker.bookProfile2Dbin("emFADCHitMap2_Profile",title,66,Help.TTEtaBinning(), 64,Help.TTPhiBinning(),"#eta","#phi");
	  title="#eta - #phi Profile Map of EM FADC > "+ thresbuffer.str() + " for triggered timeslice";
	  m_p_TT_HitMap_emADC_00100=ADCHitMaps_Booker.bookProfile2Dbin("ppm_em_2d_etaPhi_tt_adc_ProfileHitMap",title,66,Help.TTEtaBinning(), 64,Help.TTPhiBinning(),"#eta","#phi");
	  //m_p_TT_HitMap_emADC_00100->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 
	  m_p_TT_HitMap_emADC_00100->SetStats(kFALSE);
	  //m_p_TT_HitMap_emADC_00100->SetOption(colz);


          //title="#eta - #phi Map of HAD FADC > "+ thresbuffer.str() + " for timeslice 2";
	  //m_h_TT_HitMap_hadADC_00100=ADCTimeSlice_Booker.book2F("hadFADCHitMap2",title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
          title="#eta - #phi Map of HAD FADC > "+ thresbuffer.str() + " for triggered timeslice";
	  m_h_TT_HitMap_hadADC_00100=ADCHitMaps_Booker.book2F("ppm_had_2d_etaPhi_tt_adc_HitMap",title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	  m_h_TT_HitMap_hadADC_00100->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 
	  m_h_TT_HitMap_hadADC_00100->SetStats(kFALSE);
	  //m_h_TT_HitMap_hadADC_00100->SetOption(colz);
	  //title="#eta - #phi Profile Map of HAD FADC > "+ thresbuffer.str() + " for timeslice 2";
	  //m_p_TT_HitMap_hadADC_00100=ADCTimeSlice_Booker.bookProfile2Dbin("hadFADCHitMap2_Profile",title,66,Help.TTEtaBinning(), 64,Help.TTPhiBinning(),"#eta","#phi");
	  title="#eta - #phi Profile Map of HAD FADC > "+ thresbuffer.str() + " for triggered timeslice";
	  m_p_TT_HitMap_hadADC_00100=ADCHitMaps_Booker.bookProfile2Dbin("ppm_had_2d_etaPhi_tt_adc_ProfileHitMap",title,66,Help.TTEtaBinning(), 64,Help.TTPhiBinning(),"#eta","#phi");
	  //m_p_TT_HitMap_hadADC_00100->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning()); 
	  m_p_TT_HitMap_hadADC_00100->SetStats(kFALSE);
	  //m_p_TT_HitMap_hadADC_00100->SetOption(colz);
	  //m_h_dist_had_max=ADCTimeSlice_Booker.book1F("haddist_maximum"," had. Distribution of Average Maximum Timeslice",5,0.5,5.5,"time slice (1-5)");
	  //m_h_dist_em_max=ADCTimeSlice_Booker.book1F("emdist_maximum"," em. Distribution of Average Maximum Timeslice",5,0.5,5.5,"time slice (1-5)");
	  m_h_dist_had_max=ADCTimeSlice_Booker.book1F("ppm_had_1d_tt_adc_MaxTimeslice"," had. Distribution of Average Maximum Timeslice",m_SliceNo,-0.5,m_SliceNo-0.5,"time slice");
	  m_h_dist_em_max=ADCTimeSlice_Booker.book1F("ppm_em_1d_tt_adc_MaxTimeslice"," em. Distribution of Average Maximum Timeslice",m_SliceNo,-0.5,m_SliceNo-0.5,"time slice");
	  
  
      //---------------------------- Timing of FADC Signal -----------------------------
     /* Comment out for now to avoid having to branch
     if (m_TT_ADCTimingPerChannel==1)
	{
	  for(;tower_it!=m_lvl1Helper->tower_end();++tower_it) 
	    {
	      Identifier towerId = (*tower_it);

	      buffer.str("");
	      etabuffer.str("");
	      phibuffer.str("");
	      buffer << towerId;
	      etabuffer<<m_l1CaloTTIdTools->IDeta(towerId);
	      phibuffer<<m_l1CaloTTIdTools->IDphi(towerId);
	      
	      if (m_lvl1Helper->sampling(towerId)==1 )
		{
		  //name = "hadFADCTiming_" + buffer.str();
		  name = "ppm_had_1d_tt_adc_Timing" + buffer.str();
		  title = "TT had FADC Timing for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HitMap_hadADCChannel_timing[towerId.get_identifier32().get_compact()] = 
                      HadADCTiming_Booker.bookProfile(name, title,m_SliceNo,-0.5,m_SliceNo-0.5, 
                                                      "TimeSlice No", "had FADC counts");
		}

	      if (m_lvl1Helper->sampling(towerId)==0) 
		{
		  //name = "emFADCTiming_" + buffer.str();
		  name = "ppm_em_1d_tt_adc_Timing" + buffer.str();
		  title = "TT em FADC Timing for #eta = " + etabuffer.str() + " | #phi = " + phibuffer.str();
		  m_h_TT_HitMap_emADCChannel_timing[towerId.get_identifier32().get_compact()] = 
                      EmADCTiming_Booker.bookProfile(name, title,m_SliceNo,-0.5,m_SliceNo-0.5, 
                                                     "TimeSlice No", "had FADC counts");
		}
	    }
	}
     */


     //-----------------------------Average Maximum Timeslice-------------------------------------------------
      
      //m_h_TT_ADC_hadTiming_signal= ADCTimeSlice_Booker.bookProfile2Dbin("ADC_hadTiming_signal","Average Maximum TimeSlice for had Signal (TS:1-5)",66,Help.TTEtaBinning(), 64,Help.TTPhiBinning(),"#eta", "#phi");
      m_h_TT_ADC_hadTiming_signal= ADCTimeSlice_Booker.bookProfile2Dbin("ppm_had_2d_etaPhi_tt_adc_MaxTimeslice","Average Maximum TimeSlice for had Signal (TS:1-15)",66,Help.TTEtaBinning(), 64,Help.TTPhiBinning(),"#eta", "#phi");

      //m_h_TT_ADC_emTiming_signal=ADCTimeSlice_Booker.bookProfile2Dbin("ADC_emTiming_signal","Average Maximum TimeSlice for em Signal (TS:1-5)",66,Help.TTEtaBinning(), 64,Help.TTPhiBinning(), "#eta", "#phi");
      m_h_TT_ADC_emTiming_signal=ADCTimeSlice_Booker.bookProfile2Dbin("ppm_em_2d_etaPhi_tt_adc_MaxTimeslice","Average Maximum TimeSlice for em Signal (TS:1-15)",66,Help.TTEtaBinning(), 64,Help.TTPhiBinning(), "#eta", "#phi");

      //---------------------------- Signal shape ------------------------------------------

      m_h_TT_SignalProfile.clear();
      int emPart = MaxPartitions/2;
      for (int p = 0; p < MaxPartitions; ++p)
        {
          if (p < emPart) name = "ppm_em_1d_tt_adc_SignalProfile" + partitionName(p);
	  else            name = "ppm_had_1d_tt_adc_SignalProfile" + partitionName(p);
	  title = "Signal Shape Profile for " + partitionName(p);
	  m_h_TT_SignalProfile.push_back(ADCTimeSlice_Booker.bookProfile(name, title, m_SliceNo, 0, m_SliceNo, "Timeslice", ""));
        }

      //---------------------------- LUT Hitmaps per threshold -----------------------------

      // Per run and last N lumiblocks - online only
      if (m_environment == AthenaMonManager::online || m_onlineTest)
        {
          m_h_TT_HitMap_emLUT_Thresh.clear();
          m_h_TT_HitMap_hadLUT_Thresh.clear();
          for (int thresh = 0; thresh < m_TT_HitMap_ThreshMax; ++thresh)
            {
	      buffer.str("");
	      buffer<<thresh;
	      TH2F* hist = HitMaps_Booker.book2F("ppm_em_2d_etaPhi_tt_lut_Threshold"+buffer.str(),"#eta - #phi Map of EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	      hist->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());
	      m_h_TT_HitMap_emLUT_Thresh.push_back(hist);
	      hist = HitMaps_Booker.book2F("ppm_had_2d_etaPhi_tt_lut_Threshold"+buffer.str(),"#eta - #phi Map of Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	      hist->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());
	      m_h_TT_HitMap_hadLUT_Thresh.push_back(hist);
            }
	  for (int block = 0; block <= m_TT_HitMap_LumiBlocks; ++block)
	    {
              buffer.str("");
	      buffer<<block;
              MonGroup lumiGroup( this, m_PathInRootFile+"/LUT/EtaPhiMaps/lumi_"+buffer.str(), expert, run );
              HistoBooker lumiBooker(&lumiGroup, &log, "");
	      for (int thresh = 0; thresh < m_TT_HitMap_ThreshMax; ++thresh)
	        {
		  buffer.str("");
		  buffer<<thresh<<"Lumi"<<block;
		  std::string name = "ppm_em_2d_etaPhi_tt_lut_Thresh"+buffer.str();
		  buffer.str("");
		  if (block == 0) buffer<<thresh<<", Current Lumi-block";
		  else            buffer<<thresh<<", Lumi-block -"<<block;
		  std::string title = "#eta - #phi Map of EM LUT > "+buffer.str();
		  TH2F* hist = lumiBooker.book2F(name,title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
		  hist->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());
		  m_h_TT_HitMap_emLUT_Thresh.push_back(hist);
		  title = "#eta - #phi Map of Had LUT > "+buffer.str();
		  buffer.str("");
		  buffer<<thresh<<"Lumi"<<block;
		  name = "ppm_had_2d_etaPhi_tt_lut_Thresh"+buffer.str();
		  hist = lumiBooker.book2F(name,title,100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
		  hist->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());
		  m_h_TT_HitMap_hadLUT_Thresh.push_back(hist);
	        }
            }
        }
      
      //---------------------------- distribution of LUT peak per detector region -----------------------------
      //m_h_TT_emLUT=LUTPeakDistribution_Booker.book1F("emLUT_peak","EM LUT: Distribution of Peak",m_MaxEnergyRange,0,m_MaxEnergyRange,"em LUT Peak [GeV]");
      //m_h_TT_emLUT_eta=LUTPeakDistribution_Booker.book1F("emLUT_eta","EM LUT: Distribution of Peak per #eta",21,-0.5,255.5,"#eta");
      m_h_TT_emLUT=LUTPeakDistribution_Booker.book1F("ppm_em_1d_tt_lut_Et","EM LUT: Distribution of Peak",m_MaxEnergyRange-1,1,m_MaxEnergyRange,"em LUT Peak [GeV]");
      m_h_TT_emLUT_eta=LUTPeakDistribution_Booker.book1F("ppm_em_1d_tt_lut_Eta","EM LUT: Distribution of Peak per #eta",21,-0.5,255.5,"#eta");
      m_h_TT_emLUT_eta->SetBins(66,Help.TTEtaBinning());
      
      //m_h_TT_emLUT_phi=LUTPeakDistribution_Booker.book1F("emLUT_phi","EM LUT: Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
      m_h_TT_emLUT_phi=LUTPeakDistribution_Booker.book1F("ppm_em_1d_tt_lut_Phi","EM LUT: Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
      m_h_TT_emLUT_phi->SetBins(64,Help.TTPhiBinning());  
      
      
      //m_h_TT_hadLUT=LUTPeakDistribution_Booker.book1F("hadLUT_peak","HAD LUT: Distribution of Peak",m_MaxEnergyRange,0,m_MaxEnergyRange,"had LUT Peak [GeV]"); 
      //m_h_TT_hadLUT_eta=LUTPeakDistribution_Booker.book1F("hadLUT_eta","HAD LUT: Distribution of Peak per #eta",256,-0.5,255.5,"#eta");
      m_h_TT_hadLUT=LUTPeakDistribution_Booker.book1F("ppm_had_1d_tt_lut_Et","HAD LUT: Distribution of Peak",m_MaxEnergyRange-1,1,m_MaxEnergyRange,"had LUT Peak [GeV]"); 
      m_h_TT_hadLUT_eta=LUTPeakDistribution_Booker.book1F("ppm_had_1d_tt_lut_Eta","HAD LUT: Distribution of Peak per #eta",256,-0.5,255.5,"#eta");
      m_h_TT_hadLUT_eta->SetBins(66,Help.TTEtaBinning());
      //m_h_TT_hadLUT_phi=LUTPeakDistribution_Booker.book1F("hadLUT_phi","HAD LUT: Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
      m_h_TT_hadLUT_phi=LUTPeakDistribution_Booker.book1F("ppm_had_1d_tt_lut_Phi","HAD LUT: Distribution of Peak per #phi",256,-0.5,255.5,"#phi");
      m_h_TT_hadLUT_phi->SetBins(64,Help.TTPhiBinning());  

      m_h_TT_BCLUT=LUTPeakDistribution_Booker.book1F("ppm_1d_tt_lut_LutPerBCN","Num of LUT > 5 per BC",0xdec,0,0xdec,"Bunch Crossing","Num. of LUT above limit");
      m_h_TT_BCID=LUTPeakDistribution_Booker.book2F("ppm_2d_tt_lut_BcidBits","PPM: Bits of BCID Logic Word Vs. LUT",8,0.,8.,256,0.,256.,"","");
      m_h_TT_BCID->GetXaxis()->SetBinLabel(1, "none");
      m_h_TT_BCID->GetXaxis()->SetBinLabel(2, "extBC only");
      m_h_TT_BCID->GetXaxis()->SetBinLabel(3, "satBC only");
      m_h_TT_BCID->GetXaxis()->SetBinLabel(4, "extBC & satBC");
      m_h_TT_BCID->GetXaxis()->SetBinLabel(5, "peakF only");
      m_h_TT_BCID->GetXaxis()->SetBinLabel(6, "extBC & peakF");
      m_h_TT_BCID->GetXaxis()->SetBinLabel(7, "satBC & peakF");
      m_h_TT_BCID->GetXaxis()->SetBinLabel(8, "all");
      m_h_TT_BCID->SetStats(kFALSE);
      

      //-------------------------Summary of Errors-----------------------------------------------

      //m_h_TT_Error=Error_Booker.book1F("TT_errorSummary","Summary of Errors",7,0.5,7.5,""); //without Ppm fw errors
      m_h_TT_Error=Error_Booker.book1F("ppm_1d_ErrorSummary","Summary of Errors",7,0.5,7.5,""); //without Ppm fw errors

      m_h_TT_Error->GetXaxis()->SetBinLabel(1, "GLinkParity");
      m_h_TT_Error->GetXaxis()->SetBinLabel(2, "GLinkProtocol");
      m_h_TT_Error->GetXaxis()->SetBinLabel(3, "FIFOOverflow");
      m_h_TT_Error->GetXaxis()->SetBinLabel(4, "ModuleError");
      m_h_TT_Error->GetXaxis()->SetBinLabel(5, "GLinkDown");
      m_h_TT_Error->GetXaxis()->SetBinLabel(6, "GLinkTimeout");
      m_h_TT_Error->GetXaxis()->SetBinLabel(7, "BCNMismatch");

      //---------------------------- SubStatus Word errors -----------------------------
      // divided in: crate, ROD status and PPm fw errors, BCN mismatch

      //BCN mismatch
      //m_h_BCNmis_Crate_03=Error_Booker.book2F("TT_BCNmis_Crate_0-3","BCN Mismatch (crates 0-3)",1,0.5,1.5,71,0.5,71.5,"","");
      m_h_BCNmis_Crate_03=Error_Booker.book2F("ppm_2d_BcnMismatch03","BCN Mismatch (crates 0-3)",1,0.5,1.5,71,0.5,71.5,"","");
       m_h_BCNmis_Crate_03->GetXaxis()->SetBinLabel(1, "BCNMismatch");

      //m_h_BCNmis_Crate_47=Error_Booker.book2F("TT_BCNmis_Crate_4-7","BCN Mismatch (crates 4-7)",1,0.5,1.5,71,0.5,71.5,"","");
      m_h_BCNmis_Crate_47=Error_Booker.book2F("ppm_2d_BcnMismatch47","BCN Mismatch (crates 4-7)",1,0.5,1.5,71,0.5,71.5,"","");
       m_h_BCNmis_Crate_47->GetXaxis()->SetBinLabel(1, "BCNMismatch");
      
       //L1Calo Substatus word
      //m_h_TT_error_Crate_03=Error_Booker.book2F("TT_error_Crate_0-3","Errors from TT SubStatus Word (crates 0-3)",6,0.5,6.5,71,0.5,71.5,"","");
      m_h_TT_error_Crate_03=Error_Booker.book2F("ppm_2d_Status03","Errors from TT SubStatus Word (crates 0-3)",6,0.5,6.5,71,0.5,71.5,"","");
            
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(1, "GLinkParity");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(2, "GLinkProtocol");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(3, "FIFOOverflow");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(4, "ModuleError");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(5, "GLinkDown");
      m_h_TT_error_Crate_03->GetXaxis()->SetBinLabel(6, "GLinkTimeout");
      

      //m_h_TT_error_Crate_47=Error_Booker.book2F("TT_error_Crate_4-7","Errors from TT SubStatus Word (crates 4-7)",6,0.5,6.5,71,0.5,71.5,"","");
      m_h_TT_error_Crate_47=Error_Booker.book2F("ppm_2d_Status47","Errors from TT SubStatus Word (crates 4-7)",6,0.5,6.5,71,0.5,71.5,"","");
            
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(1, "GLinkParity");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(2, "GLinkProtocol");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(3, "FIFOOverflow");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(4, "ModuleError");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(5, "GLinkDown");
      m_h_TT_error_Crate_47->GetXaxis()->SetBinLabel(6, "GLinkTimeout");
      

      //error bit field from ASIC data
     //m_h_fwPpmError_Crate_03=Error_Booker.book2F("TT_fwError_Crate_0-3","Errors from ASIC error field (crates 0-3)",8,0.5,8.5,71,0.5,71.5,"","");
     m_h_fwPpmError_Crate_03=Error_Booker.book2F("ppm_2d_ErrorField03","Errors from ASIC error field (crates 0-3)",8,0.5,8.5,71,0.5,71.5,"","");
     
      m_h_fwPpmError_Crate_03->GetXaxis()->SetBinLabel(1, "ChannelDisabled");
      m_h_fwPpmError_Crate_03->GetXaxis()->SetBinLabel(2, "MCMAbsent");
      m_h_fwPpmError_Crate_03->GetXaxis()->SetBinLabel(3, "Timeout");
      m_h_fwPpmError_Crate_03->GetXaxis()->SetBinLabel(4, "ASICFull");
      m_h_fwPpmError_Crate_03->GetXaxis()->SetBinLabel(5, "EventMismatch");
      m_h_fwPpmError_Crate_03->GetXaxis()->SetBinLabel(6, "BunchMismatch");
      m_h_fwPpmError_Crate_03->GetXaxis()->SetBinLabel(7, "FIFOCorrupt");
      m_h_fwPpmError_Crate_03->GetXaxis()->SetBinLabel(8, "PinParity");
 
     //m_h_fwPpmError_Crate_47=Error_Booker.book2F("TT_fwError_Crate_4-7","Errors from ASIC error field (crates 4-7)",8,0.5,8.5,71,0.5,71.5,"","");
     m_h_fwPpmError_Crate_47=Error_Booker.book2F("ppm_2d_ErrorField47","Errors from ASIC error field (crates 4-7)",8,0.5,8.5,71,0.5,71.5,"","");

      m_h_fwPpmError_Crate_47->GetXaxis()->SetBinLabel(1, "ChannelDisabled");
      m_h_fwPpmError_Crate_47->GetXaxis()->SetBinLabel(2, "MCMAbsent");
      m_h_fwPpmError_Crate_47->GetXaxis()->SetBinLabel(3, "Timeout");
      m_h_fwPpmError_Crate_47->GetXaxis()->SetBinLabel(4, "ASICFull");
      m_h_fwPpmError_Crate_47->GetXaxis()->SetBinLabel(5, "EventMismatch");
      m_h_fwPpmError_Crate_47->GetXaxis()->SetBinLabel(6, "BunchMismatch");
      m_h_fwPpmError_Crate_47->GetXaxis()->SetBinLabel(7, "FIFOCorrupt");
      m_h_fwPpmError_Crate_47->GetXaxis()->SetBinLabel(8, "PinParity");

    
      for (int i=5; i<21; i++)
	{
	  buffer.str("");
	  buffer<<i-5;
	  
	  name = "PPM " + buffer.str();
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel((i-4)+0*18, name.c_str());
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel((i-4)+1*18, name.c_str());
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel((i-4)+2*18, name.c_str());
	  m_h_TT_error_Crate_03->GetYaxis()->SetBinLabel((i-4)+3*18, name.c_str());
          
	  m_h_fwPpmError_Crate_03->GetYaxis()->SetBinLabel((i-4)+0*18, name.c_str());
	  m_h_fwPpmError_Crate_03->GetYaxis()->SetBinLabel((i-4)+1*18, name.c_str());
	  m_h_fwPpmError_Crate_03->GetYaxis()->SetBinLabel((i-4)+2*18, name.c_str());
	  m_h_fwPpmError_Crate_03->GetYaxis()->SetBinLabel((i-4)+3*18, name.c_str());

	  
	  m_h_BCNmis_Crate_03->GetYaxis()->SetBinLabel((i-4)+0*18, name.c_str());
	  m_h_BCNmis_Crate_03->GetYaxis()->SetBinLabel((i-4)+1*18, name.c_str());
	  m_h_BCNmis_Crate_03->GetYaxis()->SetBinLabel((i-4)+2*18, name.c_str());
	  m_h_BCNmis_Crate_03->GetYaxis()->SetBinLabel((i-4)+3*18, name.c_str());

	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel((i-4)+(4-4)*18, name.c_str());
	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel((i-4)+(5-4)*18, name.c_str());
	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel((i-4)+(6-4)*18, name.c_str());
	  m_h_TT_error_Crate_47->GetYaxis()->SetBinLabel((i-4)+(7-4)*18, name.c_str());

	  m_h_fwPpmError_Crate_47->GetYaxis()->SetBinLabel((i-4)+(4-4)*18, name.c_str());
	  m_h_fwPpmError_Crate_47->GetYaxis()->SetBinLabel((i-4)+(5-4)*18, name.c_str());
	  m_h_fwPpmError_Crate_47->GetYaxis()->SetBinLabel((i-4)+(6-4)*18, name.c_str());
	  m_h_fwPpmError_Crate_47->GetYaxis()->SetBinLabel((i-4)+(7-4)*18, name.c_str());

	  m_h_BCNmis_Crate_47->GetYaxis()->SetBinLabel((i-4)+(4-4)*18, name.c_str());
	  m_h_BCNmis_Crate_47->GetYaxis()->SetBinLabel((i-4)+(5-4)*18, name.c_str());
	  m_h_BCNmis_Crate_47->GetYaxis()->SetBinLabel((i-4)+(6-4)*18, name.c_str());
	  m_h_BCNmis_Crate_47->GetYaxis()->SetBinLabel((i-4)+(7-4)*18, name.c_str());


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


	  m_h_fwPpmError_Crate_03->GetYaxis()->SetBinLabel(17+0*18, "Crate 0");
	  m_h_fwPpmError_Crate_03->GetYaxis()->SetBinLabel(17+1*18, "Crate 1");
	  m_h_fwPpmError_Crate_03->GetYaxis()->SetBinLabel(17+2*18, "Crate 2");
	  m_h_fwPpmError_Crate_03->GetYaxis()->SetBinLabel(17+3*18, "Crate 3");

	  m_h_fwPpmError_Crate_47->GetYaxis()->SetBinLabel(17+(4-4)*18, "Crate 4");
	  m_h_fwPpmError_Crate_47->GetYaxis()->SetBinLabel(17+(5-4)*18, "Crate 5");
	  m_h_fwPpmError_Crate_47->GetYaxis()->SetBinLabel(17+(6-4)*18, "Crate 6");
	  m_h_fwPpmError_Crate_47->GetYaxis()->SetBinLabel(17+(7-4)*18, "Crate 7");

	  m_h_BCNmis_Crate_03->GetYaxis()->SetBinLabel(17+0*18, "Crate 0");
	  m_h_BCNmis_Crate_03->GetYaxis()->SetBinLabel(17+1*18, "Crate 1");
	  m_h_BCNmis_Crate_03->GetYaxis()->SetBinLabel(17+2*18, "Crate 2");
	  m_h_BCNmis_Crate_03->GetYaxis()->SetBinLabel(17+3*18, "Crate 3");

	  m_h_BCNmis_Crate_47->GetYaxis()->SetBinLabel(17+(4-4)*18, "Crate 4");
	  m_h_BCNmis_Crate_47->GetYaxis()->SetBinLabel(17+(5-4)*18, "Crate 5");
	  m_h_BCNmis_Crate_47->GetYaxis()->SetBinLabel(17+(6-4)*18, "Crate 6");
	  m_h_BCNmis_Crate_47->GetYaxis()->SetBinLabel(17+(7-4)*18, "Crate 7");

      m_h_ErrorDetails.clear();
      std::vector<std::string> errNames;
      errNames.push_back("Channel0Disabled");
      errNames.push_back("Channel1Disabled");
      errNames.push_back("Channel2Disabled");
      errNames.push_back("Channel3Disabled");
      errNames.push_back("MCMAbsent");
      errNames.push_back("");
      errNames.push_back("Timeout");
      errNames.push_back("ASICFull");
      errNames.push_back("EventMismatch");
      errNames.push_back("BunchMismatch");
      errNames.push_back("FIFOCorrupt");
      errNames.push_back("PinParity");
      for (int error = 0; error < 12; error+=2) 
        {
          for (int crate = 0; crate < 8; crate+=2)
            {
	      buffer.str("");
	      buffer<<crate;
	      std::string name = "ppm_2d_"+errNames[error]+errNames[error+1]+"Crate"+buffer.str();
	      std::string title = "ASIC Errors "+errNames[error]+" "+errNames[error+1]+" for Crates "+buffer.str();
	      buffer.str("");
	      buffer<<(crate+1);
	      name += buffer.str();
	      title += "-"+buffer.str();
	      TH2F* hist = 0;
	      if (error != 4) hist = ErrorDetail_Booker.book2F(name,title,32,0,32,32,0,32,"MCM","Crate/Module");
	      else            hist = ErrorDetail_Booker.book2F(name,title,16,0,16,32,0,32,"MCM","Crate/Module");
	      m_h_ErrorDetails.push_back(hist);
	      for (int mcm = 0; mcm < 16; mcm+=2)
	        {
		  if (mcm == 0)
		    {
		      hist->GetXaxis()->SetBinLabel(1, errNames[error].c_str());
		      if (error != 4) hist->GetXaxis()->SetBinLabel(17, errNames[error+1].c_str());
                    }
                  else
		    {
		      buffer.str("");
		      buffer<<mcm;
		      hist->GetXaxis()->SetBinLabel(1+mcm, buffer.str().c_str());
		      if (error != 4) hist->GetXaxis()->SetBinLabel(17+mcm, buffer.str().c_str());
                    }
                }
              for (int cr = 0; cr < 2; ++cr)
	        {
		  for (int module = 0; module < 16; module+=2)
		    {
		      buffer.str("");
		      buffer<<(cr+crate)<<"/"<<module;
		      hist->GetYaxis()->SetBinLabel(1+cr*16+module, buffer.str().c_str());
                    }
                }
            }
        }

	  
      //---------------------------- number of triggered slice -----------------------------
      //m_h_TT_triggeredSlice_em=ADCTimeSlice_Booker.book1F("TT_EMTriggeredSlice","Number of the EM Triggered Slice",7,-0.5,6.5,"#Slice");
      //m_h_TT_triggeredSlice_had=ADCTimeSlice_Booker.book1F("TT_HADTriggeredSlice","Number of the HAD Triggered Slice",7,-0.5,6.5,"#Slice");
      m_h_TT_triggeredSlice_em=ADCTimeSlice_Booker.book1F("ppm_em_1d_tt_adc_TriggeredSlice","Number of the EM Triggered Slice",m_SliceNo,-0.5,m_SliceNo-0.5,"#Slice");
      m_h_TT_triggeredSlice_had=ADCTimeSlice_Booker.book1F("ppm_had_1d_tt_adc_TriggeredSlice","Number of the HAD Triggered Slice",m_SliceNo,-0.5,m_SliceNo-0.5,"#Slice");
      
   
      //----------------------------- number of events ----------------------------------
      //m_h_NumberEvents= NoEvent_Booker.book1F("NumberEvents","Number of processed events",1,0.5,1.5,"");
      m_h_NumberEvents= NoEvent_Booker.book1F("l1calo_1d_NumberOfEvents","Number of processed events",2,0.5,2.5,"");
      m_h_NumberEvents->GetXaxis()->SetBinLabel(1,"Processed Events");
      m_h_NumberEvents->GetXaxis()->SetBinLabel(2,"Corrupt Events Skipped");
	     
	}	

  if ( isNewLumiBlock )
    {
      //---------------------------- LUT Hitmaps per threshold -----------------------------
      if (m_environment == AthenaMonManager::online || m_onlineTest)
        {
	  // Current lumi copied to lumi-1 and so on
	  for (int block = m_TT_HitMap_LumiBlocks-1; block >= 0; --block)
	    {
	      for (int thresh = 0; thresh < m_TT_HitMap_ThreshMax; ++thresh)
	        {
		  TH2F* hist1 = m_h_TT_HitMap_emLUT_Thresh[(block+1)*m_TT_HitMap_ThreshMax + thresh];
		  TH2F* hist2 = m_h_TT_HitMap_emLUT_Thresh[(block+2)*m_TT_HitMap_ThreshMax + thresh];
		  TH2F* hist3 = m_h_TT_HitMap_hadLUT_Thresh[(block+1)*m_TT_HitMap_ThreshMax + thresh];
		  TH2F* hist4 = m_h_TT_HitMap_hadLUT_Thresh[(block+2)*m_TT_HitMap_ThreshMax + thresh];
		  hist2->Reset();
		  hist4->Reset();
		  for (int binx = 1; binx <= hist1->GetNbinsX(); ++binx)
		    {
		      for (int biny = 1; biny <= hist1->GetNbinsY(); biny++)
		        {
			  double val = hist1->GetBinContent(binx, biny);
			  if (val) hist2->SetBinContent(binx, biny, val);
			  val = hist3->GetBinContent(binx, biny);
			  if (val) hist4->SetBinContent(binx, biny, val);
                        }
                    }
		  if (block == 0)
		    {
		      hist1->Reset();
		      hist3->Reset();
		    }
                }
            }
        }
      else
        {
	  // Offline - per lumiblock - merge will give per run
          Helper Help;
          m_h_TT_HitMap_emLUT_Thresh.clear();
          m_h_TT_HitMap_hadLUT_Thresh.clear();
          MonGroup TT_LumiHitMaps( this, m_PathInRootFile+"/LUT/EtaPhiMaps", expert, lumiBlock );
          HistoBooker LumiHitMaps_Booker(&TT_LumiHitMaps, &log, "");
          std::stringstream buffer;
          for (int thresh = 0; thresh < m_TT_HitMap_ThreshMax; ++thresh)
            {
	      buffer.str("");
	      buffer<<thresh;
	      TH2F* hist = LumiHitMaps_Booker.book2F("ppm_em_2d_etaPhi_tt_lut_Threshold"+buffer.str(),"#eta - #phi Map of EM LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	      hist->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());
	      m_h_TT_HitMap_emLUT_Thresh.push_back(hist);
	      hist = LumiHitMaps_Booker.book2F("ppm_had_2d_etaPhi_tt_lut_Threshold"+buffer.str(),"#eta - #phi Map of Had LUT > "+buffer.str(),100,-4.9,4.9, 64,0,2*M_PI,"#eta","#phi");
	      hist->SetBins(66,Help.TTEtaBinning(),64,Help.TTPhiBinning());
	      m_h_TT_HitMap_hadLUT_Thresh.push_back(hist);
            }
	  
        }
    }
    
  return StatusCode::SUCCESS;
	}

/*---------------------------------------------------------*/
StatusCode PPrMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  log << MSG::DEBUG << "in fillHistograms()" << endreq;

  // Skip events believed to be corrupt

  if (m_errorTool->corrupt()) {
    m_h_NumberEvents->Fill(2,1);
    log << MSG::DEBUG << "Skipping corrupt event" << endreq;
    return StatusCode::SUCCESS;
  }
  m_h_NumberEvents->Fill(1,1);  
  m_NoEvents++;

  // Error vector for global overview
  std::vector<int> overview(8);
  
  //Retrieve TriggerTowers from SG
  const TriggerTowerCollection* TriggerTowerTES = 0; 
  StatusCode sc = m_storeGate->retrieve(TriggerTowerTES, m_TriggerTowerContainerName); 
  if( (sc==StatusCode::FAILURE) ) 
    {
      log << MSG::INFO << "No TriggerTower found in TES at "<< m_TriggerTowerContainerName<< endreq ;
      return StatusCode::SUCCESS;
    }

  // Get Bunch crossing number from EventInfo
  uint32_t bunchCrossing = 0;
  const EventInfo* evInfo = 0;
  sc = m_storeGate->retrieve(evInfo);
  if (sc.isFailure()) {
    log << MSG::DEBUG << "No EventInfo found" << endreq;
  } else {
    const EventID* evID = evInfo->event_ID();
    if (evID) bunchCrossing = evID->bunch_crossing_id();
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
          // Bunch crossing and BCID bits
          if (EmEnergy>5) m_h_TT_BCLUT->Fill(bunchCrossing);
          m_h_TT_BCID->Fill((*TriggerTowerIterator)->emBCID(), EmEnergy);
	}
	 
       //---------------------------- EM LUT HitMaps -----------------------------
      for (int thresh = 0; thresh < m_TT_HitMap_ThreshMax; ++thresh)
        {
	  if (EmEnergy > thresh)
	    {
	      m_h_TT_HitMap_emLUT_Thresh[thresh]->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	      if (m_environment == AthenaMonManager::online || m_onlineTest)
	        {
	          m_h_TT_HitMap_emLUT_Thresh[thresh+m_TT_HitMap_ThreshMax]->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	        }
            }
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
          // Bunch crossing and BCID bits
          if (HadEnergy>5) m_h_TT_BCLUT->Fill(bunchCrossing);
          m_h_TT_BCID->Fill((*TriggerTowerIterator)->hadBCID(), HadEnergy);
	}
    
       //---------------------------- had LUT HitMaps -----------------------------
      for (int thresh = 0; thresh < m_TT_HitMap_ThreshMax; ++thresh)
        {
	  if (HadEnergy > thresh)
	    {
	      m_h_TT_HitMap_hadLUT_Thresh[thresh]->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	      if (m_environment == AthenaMonManager::online || m_onlineTest)
	        {
	          m_h_TT_HitMap_hadLUT_Thresh[thresh+m_TT_HitMap_ThreshMax]->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
	        }
            }
        }
    


     //---------------------------- ADC HitMaps per timeslice -----------------------------


      int tslice=(*TriggerTowerIterator)->emADCPeak();
    
if (tslice<static_cast<int>(( (*TriggerTowerIterator)->emADC()).size()))
	    {
	      if (( (*TriggerTowerIterator)->emADC())[tslice] > m_TT_ADC_HitMap_Thresh)
		{
		  m_h_TT_HitMap_emADC_00100 ->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
		  m_p_TT_HitMap_emADC_00100 ->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),( (*TriggerTowerIterator)->emADC())[tslice]);
		}
	    }


      tslice=(*TriggerTowerIterator)->hadADCPeak();

 if (tslice<static_cast<int>(( (*TriggerTowerIterator)->hadADC()).size()))
	    {
	      if (( (*TriggerTowerIterator)->hadADC())[tslice] > m_TT_ADC_HitMap_Thresh)
		{
		  m_h_TT_HitMap_hadADC_00100 ->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),1);
		  m_p_TT_HitMap_hadADC_00100 ->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),( (*TriggerTowerIterator)->hadADC())[tslice]); 
		}

	    }

	             
        //---------------------------- Timing of FADC Signal -----------------------------
      int hadFADCSum=0;     
      int emFADCSum=0;     
      double max;

      emFADCSum=FADCSum((*TriggerTowerIterator)->emADC());
      hadFADCSum=FADCSum((*TriggerTowerIterator)->hadADC());

      /* Comment out for now to avoid having to branch
      if (m_TT_ADCTimingPerChannel==1)
	{
	  for (int i=0; i<m_SliceNo;i++)
	    {
	      if (i<static_cast<int>(( (*TriggerTowerIterator)->emADC()).size()))
		{ 
		  if (emFADCSum>m_EMFADCCut) 
		    {
		      m_h_TT_HitMap_emADCChannel_timing[EmTowerId.get_identifier32().get_compact()]->Fill(i,(*TriggerTowerIterator)->emADC()[i]);
		    }
		}
	      
	      if (i<static_cast<int>(( (*TriggerTowerIterator)->hadADC()).size()))
		{
		  if (hadFADCSum>m_HADFADCCut) 
		    {
		      m_h_TT_HitMap_hadADCChannel_timing[HadTowerId.get_identifier32().get_compact()]->Fill(i,(*TriggerTowerIterator)->hadADC()[i]);
		    }
		}
	    }
	}
      */

      if (emFADCSum>m_EMFADCCut)
	{
	  max = recTime((*TriggerTowerIterator)->emADC());
	  //log << MSG::INFO << "TimeSlice of Maximum "<< max<< endreq ;
	  if (max >= 0.) {
	    m_h_TT_ADC_emTiming_signal->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),max+1);
	    m_h_dist_em_max->Fill(max);
	  }
	}

      if (hadFADCSum>m_HADFADCCut)
	{
	  max = recTime((*TriggerTowerIterator)->hadADC());
	  //log << MSG::INFO << "TimeSlice of Maximum "<< max<< endreq ;
	  if (max >= 0.) {
	    m_h_TT_ADC_hadTiming_signal->Fill((*TriggerTowerIterator)->eta(),(*TriggerTowerIterator)->phi(),max+1);
	    m_h_dist_had_max->Fill(max);
	  }
        }

      //------------------------ Signal shape profile ---------------------

      if ((*TriggerTowerIterator)->emEnergy() > 0)
        {
          const int emPart  = partition(0, (*TriggerTowerIterator)->eta());
	  const std::vector<int>& emADC((*TriggerTowerIterator)->emADC());
	  std::vector<int>::const_iterator it  = emADC.begin();
	  std::vector<int>::const_iterator itE = emADC.end();
	  for (int slice = 0; it != itE && slice < m_SliceNo; ++it, ++slice) m_h_TT_SignalProfile[emPart]->Fill(slice, *it);
        }
      if ((*TriggerTowerIterator)->hadEnergy() > 0)
        {
          const int hadPart = partition(1, (*TriggerTowerIterator)->eta());
	  const std::vector<int>& hadADC((*TriggerTowerIterator)->hadADC());
	  std::vector<int>::const_iterator it  = hadADC.begin();
	  std::vector<int>::const_iterator itE = hadADC.end();
	  for (int slice = 0; it != itE && slice < m_SliceNo; ++it, ++slice) m_h_TT_SignalProfile[hadPart]->Fill(slice, *it);
        }
    

      //---------------------------- SubStatus Word errors -----------------------------
      //----------------------------- em ---------------------------------------

      using LVL1::DataError;
      DataError emerr((*TriggerTowerIterator)-> emError());

      int crate=0;
      int module=0;
      int submodule=0;
      int channel=0;
      try 
	{
	  HWIdentifier ttOnlId = m_ttSvc->createTTChannelID(EmTowerId);
	  crate     = m_l1ttonlineHelper->crate(ttOnlId);
	  module    = m_l1ttonlineHelper->module(ttOnlId);
	  submodule = m_l1ttonlineHelper->submodule(ttOnlId);
	  channel   = m_l1ttonlineHelper->channel(ttOnlId);

	  //log << MSG::DEBUG << "em PPM crate: " << crate<<"  module: "<<module << " submodule "<<submodule<<" channel "<< channel<< endreq ;

	 }
 
      catch(CaloID_Exception& except) 
	{
	  log << MSG::ERROR << "CaloID_Exception " << (std::string) except << endreq ;
	}
   
      //Summary

      if (emerr.get(DataError::GLinkParity))   m_h_TT_Error->Fill(1);
      if (emerr.get(DataError::GLinkProtocol)) m_h_TT_Error->Fill(2);
      if (emerr.get(DataError::FIFOOverflow))  m_h_TT_Error->Fill(3);
      if (emerr.get(DataError::ModuleError))   m_h_TT_Error->Fill(4);
      if (emerr.get(DataError::GLinkDown))     m_h_TT_Error->Fill(5);
      if (emerr.get(DataError::GLinkTimeout))  m_h_TT_Error->Fill(6);
      if (emerr.get(DataError::BCNMismatch))   m_h_TT_Error->Fill(7);
      

      // em signals Crate 0-3
      //em+had FCAL signals get processed in one crate (Crates 4-7)

      if (crate>3) {
	//---------------- per crate and module --------------------  m_h_TT_error_Crate_47
      const int ypos = (module-4)+((crate-4)*18);
      if (emerr.get(DataError::ChannelDisabled)) m_h_fwPpmError_Crate_47->Fill(1,ypos);
      if (emerr.get(DataError::MCMAbsent))       m_h_fwPpmError_Crate_47->Fill(2,ypos);
      if (emerr.get(DataError::Timeout))         m_h_fwPpmError_Crate_47->Fill(3,ypos);
      if (emerr.get(DataError::ASICFull))        m_h_fwPpmError_Crate_47->Fill(4,ypos);
      if (emerr.get(DataError::EventMismatch))   m_h_fwPpmError_Crate_47->Fill(5,ypos);
      if (emerr.get(DataError::BunchMismatch))   m_h_fwPpmError_Crate_47->Fill(6,ypos);
      if (emerr.get(DataError::FIFOCorrupt))     m_h_fwPpmError_Crate_47->Fill(7,ypos);
      if (emerr.get(DataError::PinParity))       m_h_fwPpmError_Crate_47->Fill(8,ypos);
      		  
      if (emerr.get(DataError::GLinkParity))     m_h_TT_error_Crate_47->Fill(1,ypos);
      if (emerr.get(DataError::GLinkProtocol))   m_h_TT_error_Crate_47->Fill(2,ypos);
      if (emerr.get(DataError::FIFOOverflow))    m_h_TT_error_Crate_47->Fill(3,ypos);
      if (emerr.get(DataError::ModuleError))     m_h_TT_error_Crate_47->Fill(4,ypos);
      if (emerr.get(DataError::GLinkDown))       m_h_TT_error_Crate_47->Fill(5,ypos);
      if (emerr.get(DataError::GLinkTimeout))    m_h_TT_error_Crate_47->Fill(6,ypos);
      
      if (emerr.get(DataError::BCNMismatch))     m_h_BCNmis_Crate_47->Fill(1,ypos);

	}

      else {
	//---------------- per crate and module -------------------------  
      const int ypos = (module-4)+(crate*18);
      if (emerr.get(DataError::ChannelDisabled)) m_h_fwPpmError_Crate_03->Fill(1,ypos);
      if (emerr.get(DataError::MCMAbsent))       m_h_fwPpmError_Crate_03->Fill(2,ypos);
      if (emerr.get(DataError::Timeout))         m_h_fwPpmError_Crate_03->Fill(3,ypos);
      if (emerr.get(DataError::ASICFull))        m_h_fwPpmError_Crate_03->Fill(4,ypos);
      if (emerr.get(DataError::EventMismatch))   m_h_fwPpmError_Crate_03->Fill(5,ypos);
      if (emerr.get(DataError::BunchMismatch))   m_h_fwPpmError_Crate_03->Fill(6,ypos);
      if (emerr.get(DataError::FIFOCorrupt))     m_h_fwPpmError_Crate_03->Fill(7,ypos);
      if (emerr.get(DataError::PinParity))       m_h_fwPpmError_Crate_03->Fill(8,ypos);
      		  
      if (emerr.get(DataError::GLinkParity))     m_h_TT_error_Crate_03->Fill(1,ypos);
      if (emerr.get(DataError::GLinkProtocol))   m_h_TT_error_Crate_03->Fill(2,ypos);
      if (emerr.get(DataError::FIFOOverflow))    m_h_TT_error_Crate_03->Fill(3,ypos);
      if (emerr.get(DataError::ModuleError))     m_h_TT_error_Crate_03->Fill(4,ypos);
      if (emerr.get(DataError::GLinkDown))       m_h_TT_error_Crate_03->Fill(5,ypos);
      if (emerr.get(DataError::GLinkTimeout))    m_h_TT_error_Crate_03->Fill(6,ypos);
      
      if (emerr.get(DataError::BCNMismatch))     m_h_BCNmis_Crate_03->Fill(1,ypos);

	}

      if (emerr.get(DataError::ChannelDisabled) || emerr.get(DataError::MCMAbsent)) overview[crate] |= 1;

      if (emerr.get(DataError::Timeout)       || emerr.get(DataError::ASICFull)      ||
          emerr.get(DataError::EventMismatch) || emerr.get(DataError::BunchMismatch) ||
	  emerr.get(DataError::FIFOCorrupt)   || emerr.get(DataError::PinParity)) overview[crate] |= (1 << 1);

      if (emerr.get(DataError::GLinkParity)  || emerr.get(DataError::GLinkProtocol) ||
          emerr.get(DataError::FIFOOverflow) || emerr.get(DataError::ModuleError)   ||
	  emerr.get(DataError::GLinkDown)    || emerr.get(DataError::GLinkTimeout)  ||
	  emerr.get(DataError::BCNMismatch)) overview[crate] |= (1 << 2);

      // Detailed plots by MCM
      int ypos = (crate%2)*16+module-5;
      if (emerr.get(DataError::ChannelDisabled)) m_h_ErrorDetails[(channel/2)*4+crate/2]->Fill((channel%2)*16+submodule, ypos);
      if (emerr.get(DataError::MCMAbsent))       m_h_ErrorDetails[8+crate/2]->Fill(submodule, ypos);
      if (emerr.get(DataError::Timeout))         m_h_ErrorDetails[12+crate/2]->Fill(submodule, ypos);
      if (emerr.get(DataError::ASICFull))        m_h_ErrorDetails[12+crate/2]->Fill(16+submodule, ypos);
      if (emerr.get(DataError::EventMismatch))   m_h_ErrorDetails[16+crate/2]->Fill(submodule, ypos);
      if (emerr.get(DataError::BunchMismatch))   m_h_ErrorDetails[16+crate/2]->Fill(16+submodule, ypos);
      if (emerr.get(DataError::FIFOCorrupt))     m_h_ErrorDetails[20+crate/2]->Fill(submodule, ypos);
      if (emerr.get(DataError::PinParity))       m_h_ErrorDetails[20+crate/2]->Fill(16+submodule, ypos);
    
     DataError haderr((*TriggerTowerIterator)-> hadError());

     //had signals in Crates 4-7
      		  
      if (haderr.get(DataError::GLinkParity))   m_h_TT_Error->Fill(1);
      if (haderr.get(DataError::GLinkProtocol)) m_h_TT_Error->Fill(2);
      if (haderr.get(DataError::FIFOOverflow))  m_h_TT_Error->Fill(3);
      if (haderr.get(DataError::ModuleError))   m_h_TT_Error->Fill(4);
      if (haderr.get(DataError::GLinkDown))     m_h_TT_Error->Fill(5);
      if (haderr.get(DataError::GLinkTimeout))  m_h_TT_Error->Fill(6);
      if (haderr.get(DataError::BCNMismatch))   m_h_TT_Error->Fill(7);

      try 
	{
	  HWIdentifier ttOnlId = m_ttSvc->createTTChannelID(HadTowerId);
	  crate     = m_l1ttonlineHelper->crate(ttOnlId);
	  module    = m_l1ttonlineHelper->module(ttOnlId);
	  submodule = m_l1ttonlineHelper->submodule(ttOnlId);
	  channel   = m_l1ttonlineHelper->channel(ttOnlId);

	  //log << MSG::DEBUG << "em PPM crate: " << crate<<"  module: "<<module << " submodule "<<submodule<<" channel "<< channel<< endreq ;

	} 

      catch(CaloID_Exception& except) 
	{
	  log << MSG::ERROR << "CaloID_Exception " << (std::string) except << endreq ;
	}
    

      //---------------- per crate and module --------------------  m_h_TT_error_Crate_03
      ypos = (module-4)+((crate-4)*18);
      if (haderr.get(DataError::ChannelDisabled)) m_h_fwPpmError_Crate_47->Fill(1,ypos);
      if (haderr.get(DataError::MCMAbsent))       m_h_fwPpmError_Crate_47->Fill(2,ypos);
      if (haderr.get(DataError::Timeout))         m_h_fwPpmError_Crate_47->Fill(3,ypos);
      if (haderr.get(DataError::ASICFull))        m_h_fwPpmError_Crate_47->Fill(4,ypos);
      if (haderr.get(DataError::EventMismatch))   m_h_fwPpmError_Crate_47->Fill(5,ypos);
      if (haderr.get(DataError::BunchMismatch))   m_h_fwPpmError_Crate_47->Fill(6,ypos);
      if (haderr.get(DataError::FIFOCorrupt))     m_h_fwPpmError_Crate_47->Fill(7,ypos);
      if (haderr.get(DataError::PinParity))       m_h_fwPpmError_Crate_47->Fill(8,ypos);
      		  
      if (haderr.get(DataError::GLinkParity))     m_h_TT_error_Crate_47->Fill(1,ypos);
      if (haderr.get(DataError::GLinkProtocol))   m_h_TT_error_Crate_47->Fill(2,ypos);
      if (haderr.get(DataError::FIFOOverflow))    m_h_TT_error_Crate_47->Fill(3,ypos);
      if (haderr.get(DataError::ModuleError))     m_h_TT_error_Crate_47->Fill(4,ypos);
      if (haderr.get(DataError::GLinkDown))       m_h_TT_error_Crate_47->Fill(5,ypos);
      if (haderr.get(DataError::GLinkTimeout))    m_h_TT_error_Crate_47->Fill(6,ypos);
      
      if (haderr.get(DataError::BCNMismatch))     m_h_BCNmis_Crate_47->Fill(1,ypos);

      if (haderr.get(DataError::ChannelDisabled) || haderr.get(DataError::MCMAbsent)) overview[crate] |= 1;

      if (haderr.get(DataError::Timeout)       || haderr.get(DataError::ASICFull)      ||
          haderr.get(DataError::EventMismatch) || haderr.get(DataError::BunchMismatch) ||
	  haderr.get(DataError::FIFOCorrupt)   || haderr.get(DataError::PinParity)) overview[crate] |= (1 << 1);

      if (haderr.get(DataError::GLinkParity)  || haderr.get(DataError::GLinkProtocol) ||
          haderr.get(DataError::FIFOOverflow) || haderr.get(DataError::ModuleError)   ||
	  haderr.get(DataError::GLinkDown)    || haderr.get(DataError::GLinkTimeout)  ||
	  haderr.get(DataError::BCNMismatch)) overview[crate] |= (1 << 2);

      // Detailed plots by MCM
      ypos = (crate%2)*16+module-5;
      if (haderr.get(DataError::ChannelDisabled)) m_h_ErrorDetails[(channel/2)*4+crate/2]->Fill((channel%2)*16+submodule, ypos);
      if (haderr.get(DataError::MCMAbsent))       m_h_ErrorDetails[8+crate/2]->Fill(submodule, ypos);
      if (haderr.get(DataError::Timeout))         m_h_ErrorDetails[12+crate/2]->Fill(submodule, ypos);
      if (haderr.get(DataError::ASICFull))        m_h_ErrorDetails[12+crate/2]->Fill(16+submodule, ypos);
      if (haderr.get(DataError::EventMismatch))   m_h_ErrorDetails[16+crate/2]->Fill(submodule, ypos);
      if (haderr.get(DataError::BunchMismatch))   m_h_ErrorDetails[16+crate/2]->Fill(16+submodule, ypos);
      if (haderr.get(DataError::FIFOCorrupt))     m_h_ErrorDetails[20+crate/2]->Fill(submodule, ypos);
      if (haderr.get(DataError::PinParity))       m_h_ErrorDetails[20+crate/2]->Fill(16+submodule, ypos);
      
     
      // number of triggered slice
      m_h_TT_triggeredSlice_em->Fill((*TriggerTowerIterator)->emADCPeak(),1);
      m_h_TT_triggeredSlice_had->Fill((*TriggerTowerIterator)->hadADCPeak(),1);

	}	     
	     
  // Write overview vector to StoreGate
  std::vector<int>* save = new std::vector<int>(overview);
  sc = m_storeGate->record(save, "L1CaloPPMErrorVector");
  if (sc != StatusCode::SUCCESS)
    {
      log << MSG::ERROR << "Error recording PPM error vector in TES "
          << endreq;
      return sc;
    }

  
  return StatusCode::SUCCESS;
}


   
/*---------------------------------------------------------*/
StatusCode PPrMon::procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun )
/*---------------------------------------------------------*/
{
  MsgStream mLog( msgSvc(), name() );
  mLog << MSG::DEBUG << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun ) 
    {  
    }
	
  /*
  if(m_Offline==1)
    {
      if( isEndOfRun ) 
	{   
	  
	  std::stringstream buffer;
	  buffer.str("");
	  buffer<<m_NoEvents;
	  std::string title;
	  
	  title = m_h_TT_Error-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_Error->SetTitle(title.c_str());
	  
	  title = m_h_TT_error_Crate_03-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_error_Crate_03->SetTitle(title.c_str());
	  
	  title = m_h_TT_error_Crate_47-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_TT_error_Crate_47->SetTitle(title.c_str());
	  
	  title = m_h_fwPpmError_Crate_03-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_fwPpmError_Crate_03->SetTitle(title.c_str());

	  title = m_h_fwPpmError_Crate_47-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_fwPpmError_Crate_47->SetTitle(title.c_str());

	  title = m_h_BCNmis_Crate_47-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_BCNmis_Crate_47->SetTitle(title.c_str());

	  title = m_h_BCNmis_Crate_03-> GetTitle();
	  title=title + " | #events: " + buffer.str();
	  m_h_BCNmis_Crate_03->SetTitle(title.c_str());
	  
	}
    }
  */
  return StatusCode::SUCCESS;
}


/*---------------------------------------------------------*/
double PPrMon::recTime(const std::vector<int>& vFAdc) {
/*---------------------------------------------------------*/

  /*
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
  
  double max = indmax+binshift;
  
  if (indmax!=0 && indmax!=(int)(vFAdc.size()-1)) {
       
    x[0] = indmax - 1 + binshift;  y[0] = vFAdc[indmax-1];
    x[1] = indmax     + binshift; y[1] = vFAdc[indmax];
    x[2] = indmax + 1 +binshift; y[2] = vFAdc[indmax+1];
    

    //This is a parabola fit function to find the maximum in ADC counts (asymmetric distribution) 
    //Simplified: max = indmax+(y[0]-y[2])/(2*(y[0]+y[2]-2*y[1]))
 
    double a = ( (x[0]-x[2])*(y[1]-y[2]) - (x[1]-x[2])*(y[0]-y[2]) ) / (
									(x[0]-x[2])*(x[1]*x[1]-x[2]*x[2]) - (x[1]-x[2])*(x[0]*x[0]-x[2]*x[2]) );
    double b = ( (y[1]-y[2])*(x[0]*x[0]-x[2]*x[2]) -
		 (y[0]-y[2])*(x[1]*x[1]-x[2]*x[2]) ) / (
							(x[1]-x[2])*(x[0]*x[0]-x[2]*x[2]) - (x[0]-x[2])*(x[1]*x[1]-x[2]*x[2]) );
    //double c = y[0] - b*x[0] - a*x[0]*x[0];
    if (a != 0.) max = -b/(2*a);
    
  }
  */

  double max = -1.;
  int slices = vFAdc.size();
  if (slices > 0) {
    max = 0.;
    int maxAdc = vFAdc[0];
    for (int sl = 1; sl < slices; ++sl) {
      if (vFAdc[sl] > maxAdc) {
        maxAdc = vFAdc[sl];
        max = sl;
      } else if (vFAdc[sl] == maxAdc) max = -1.;
    }
    if (maxAdc == 0) max = -1.;
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
  FADCSum/=vFAdc.size();         // average per slice
  
  return FADCSum;
}

/*---------------------------------------------------------*/
int PPrMon::partition(int layer, double eta) {
/*---------------------------------------------------------*/

  int part = 0;
  if (layer == 0) {
    if      (eta < -3.2) part = LArFCAL1C;
    else if (eta < -1.5) part = LArEMECC;
    else if (eta < -1.4) part = LArOverlapC;
    else if (eta <  0.0) part = LArEMBC;
    else if (eta <  1.4) part = LArEMBA;
    else if (eta <  1.5) part = LArOverlapA;
    else if (eta <  3.2) part = LArEMECA;
    else                 part = LArFCAL1A;
  } else {
    if      (eta < -3.2) part = LArFCAL23C;
    else if (eta < -1.5) part = LArHECC;
    else if (eta < -0.9) part = TileEBC;
    else if (eta <  0.0) part = TileLBC;
    else if (eta <  0.9) part = TileLBA;
    else if (eta <  1.5) part = TileEBA;
    else if (eta <  3.2) part = LArHECA;
    else                 part = LArFCAL23A;
  }
  return part;
}

/*---------------------------------------------------------*/
std::string PPrMon::partitionName(int part) {
/*---------------------------------------------------------*/

  std::string name = "";
  switch (part) {
    case LArFCAL1C:   name = "LArFCAL1C";   break;
    case LArEMECC:    name = "LArEMECC";    break;
    case LArOverlapC: name = "LArOverlapC"; break;
    case LArEMBC:     name = "LArEMBC";     break;
    case LArEMBA:     name = "LArEMBA";     break;
    case LArOverlapA: name = "LArOverlapA"; break;
    case LArEMECA:    name = "LArEMECA";    break;
    case LArFCAL1A:   name = "LArFCAL1A";   break;
    case LArFCAL23C:  name = "LArFCAL23C";  break;
    case LArHECC:     name = "LArHECC";     break;
    case TileEBC:     name = "TileEBC";     break;
    case TileLBC:     name = "TileLBC";     break;
    case TileLBA:     name = "TileLBA";     break;
    case TileEBA:     name = "TileEBA";     break;
    case LArHECA:     name = "LArHECA";     break;
    case LArFCAL23A:  name = "LArFCAL23A";  break;
    default:          name = "Unknown";     break;
  }
  return name;
}


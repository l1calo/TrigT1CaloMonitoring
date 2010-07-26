// ********************************************************************
//
// NAME:     PPrSpareMon.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Peter Faulkner
//           
//
// ********************************************************************

#include "LWHists/LWHist.h"
#include "LWHists/TH1F_LW.h"
#include "LWHists/TH2F_LW.h"
#include "LWHists/TH2I_LW.h"
#include "LWHists/TProfile2D_LW.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"
#include "SGTools/StlVectorClids.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1CaloMonitoring/PPrSpareMon.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloLWHistogramTool.h"

#include "TrigT1CaloEvent/TriggerTower_ClassDEF.h"
#include "TrigT1CaloEvent/TriggerTowerCollection.h"
#include "TrigT1CaloUtils/DataError.h"


/*---------------------------------------------------------*/
PPrSpareMon::PPrSpareMon(const std::string & type, const std::string & name,
					           const IInterface* parent)
  : ManagedMonitorToolBase ( type, name, parent ),
    m_SliceNo(15),
    m_errorTool("TrigT1CaloMonErrorTool"),
    m_histTool("TrigT1CaloLWHistogramTool")
/*---------------------------------------------------------*/
{
  declareProperty("BS_TriggerTowerContainer",
                  m_TriggerTowerContainerName = "TriggerTowersSpare");
  declareProperty("ADCHitMap_Thresh",  m_TT_ADC_HitMap_Thresh = 40);

  declareProperty("PathInRootFile",
                  m_PathInRootFile="L1Calo/PPM/SpareChannels") ;
  declareProperty("ErrorPathInRootFile",
                  m_ErrorPathInRootFile="L1Calo/PPM/SpareChannels/Errors") ;
  declareProperty("OnlineTest", m_onlineTest = false,
                  "Test online code when running offline");

}

/*---------------------------------------------------------*/
PPrSpareMon::~PPrSpareMon()
/*---------------------------------------------------------*/
{
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*---------------------------------------------------------*/
StatusCode PPrSpareMon::initialize()
/*---------------------------------------------------------*/
{
  msg(MSG::INFO) << "Initializing " << name() << " - package version "
                 << PACKAGE_VERSION << endreq;

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_errorTool.retrieve();
  if( sc.isFailure() ) {
    msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloMonErrorTool"
                    << endreq;
    return sc;
  }

  sc = m_histTool.retrieve();
  if( sc.isFailure() ) {
    msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloLWHistogramTool"
                    << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPrSpareMon::bookHistograms( bool isNewEventsBlock,
                                        bool isNewLumiBlock, bool isNewRun )
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "in PPrSpareMon::bookHistograms" << endreq;

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  if ( isNewEventsBlock|| isNewLumiBlock) { }

  if( isNewRun ) {

    MonGroup TT_ADC(this, m_PathInRootFile+"/ADC", shift, run);
    MonGroup TT_Error(this, m_ErrorPathInRootFile, shift, run);
    MonGroup TT_ErrorDetail(this, m_ErrorPathInRootFile+"/Detail", expert, run);

    std::string name,title;
    std::stringstream buffer;

    //------------------- ADC Hitmaps for Triggered Timeslice ----------------

    m_histTool->setMonGroup(&TT_ADC);

    buffer.str("");
    buffer << m_TT_ADC_HitMap_Thresh;

    title = "Spare Channels Hit Map of FADC > " + buffer.str()
                                                + " for Triggered Timeslice";
    m_h_TT_HitMap_ADC = m_histTool->bookPPMCrateModuleVsSubmoduleChannel(
                         "ppmspare_2d_tt_adc_HitMap", title, 2, 5);
    title = "Spare Channels Profile Map of FADC for Triggered Timeslice";
    m_p_TT_HitMap_ADC = m_histTool->bookProfilePPMCrateModuleVsSubmoduleChannel(
                         "ppmspare_2d_tt_adc_ProfileMap", title, 2, 5);

    //-------------------------Summary of Errors------------------------------

    m_histTool->setMonGroup(&TT_Error);

    m_h_TT_Error = m_histTool->book1F("ppmspare_1d_ErrorSummary",
                   "Spare Channels Summary of SubStatus Errors", 8, 0., 8.);
    m_histTool->subStatus(m_h_TT_Error);

    m_h_TT_EventNumbers = m_histTool->bookEventNumbers(
      "ppmspare_2d_ErrorEventNumbers",
      "Spare Channels SubStatus Error Event Numbers", 8, 0., 8.);
    m_histTool->subStatus(m_h_TT_EventNumbers, 0, false);
    m_h_TT_ASICEventNumbers = m_histTool->bookEventNumbers(
      "ppmspare_2d_ASICErrorEventNumbers",
      " Spare Channels ASIC Error Field Event Numbers", 8, 0., 8.);
    m_histTool->ppmErrors(m_h_TT_ASICEventNumbers, 0, false);

    //---------------------- SubStatus Word errors ---------------------------
      
    //L1Calo Substatus word
    m_h_TT_error_Crate_25 = m_histTool->bookPPMSubStatusVsCrateModule(
      "ppmspare_2d_Status25",
      "Spare Channels: Errors from TT SubStatus Word (crates 2-5)", 2, 5);

    //error bit field from ASIC data
    m_h_fwPpmError_Crate_25 = m_histTool->bookPPMErrorsVsCrateModule(
      "ppmspare_2d_ErrorField25",
      "Spare Channels: Errors from ASIC error field (crates 2-5)", 2, 5);

    m_histTool->setMonGroup(&TT_ErrorDetail);

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
    for (int error = 0; error < 12; error+=2) {
      for (int crate = 2; crate < 6; crate+=2) {
	buffer.str("");
	buffer << crate;
	name = "ppmspare_2d_"+errNames[error]+errNames[error+1]
	                     +"Crate"+buffer.str();
	title = "ASIC Errors "+errNames[error]+" "
	                      +errNames[error+1]+" for Crates "+buffer.str();
	buffer.str("");
	buffer << (crate+1);
	name += buffer.str();
	title += "-"+buffer.str();
	int nbins = (error != 4) ? 32 : 16;
	TH2F_LW* hist = m_histTool->book2F(name,title,nbins,0,nbins,32,0,32);
	m_histTool->numbers(hist, 0, 15, 2);
	LWHist::LWHistAxis* axis = hist->GetXaxis();
	axis->SetBinLabel(1, errNames[error].c_str());
	if (error != 4) {
	  m_histTool->numbers(hist, 0, 15, 2, 16);
	  axis->SetBinLabel(17, errNames[error+1].c_str());
        }
	axis->SetTitle("MCM");
	m_histTool->ppmCrateModule(hist, crate, crate+1, 0, false);
	m_h_ErrorDetails.push_back(hist);
      }
    }
	  
    //--------------------- number of triggered slice ------------------------
    m_histTool->setMonGroup(&TT_ADC);
    
    m_h_TT_triggeredSlice = m_histTool->book1F(
      "ppmspare_1d_tt_adc_TriggeredSlice",
      "Spare Channels Number of the Triggered Slice", m_SliceNo, 0, m_SliceNo);
    m_histTool->numbers(m_h_TT_triggeredSlice, 0, m_SliceNo-1);
     
    m_histTool->unsetMonGroup();
  }	

  if ( isNewLumiBlock ) { }
    
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode PPrSpareMon::fillHistograms()
/*---------------------------------------------------------*/
{
  const bool debug = msgLvl(MSG::DEBUG);
  if (debug) msg(MSG::DEBUG) << "in fillHistograms()" << endreq;

  // Skip events believed to be corrupt

  if (m_errorTool->corrupt()) {
    if (debug) msg(MSG::DEBUG) << "Skipping corrupt event" << endreq;
    return StatusCode::SUCCESS;
  }

  // Error vector for global overview
  std::vector<int> overview(8);

  //Retrieve TriggerTowers from SG
  StatusCode sc;
  const TriggerTowerCollection* TriggerTowerTES = 0; 
  if (evtStore()->contains<TriggerTowerCollection>(m_TriggerTowerContainerName)) {
    sc = evtStore()->retrieve(TriggerTowerTES, m_TriggerTowerContainerName); 
  } else sc = StatusCode::FAILURE;
  if (sc.isFailure()) {
    if (debug) msg(MSG::DEBUG) << "No TriggerTower found in TES at "
                               << m_TriggerTowerContainerName<< endreq ;
    return StatusCode::SUCCESS;
  }

    
  // =========================================================================
  // ================= Container: TriggerTower ===============================
  // =========================================================================

  TriggerTowerCollection::const_iterator TriggerTowerIterator =
                                                     TriggerTowerTES->begin(); 
  TriggerTowerCollection::const_iterator TriggerTowerIteratorEnd =
                                                     TriggerTowerTES->end(); 
 
  for (; TriggerTowerIterator != TriggerTowerIteratorEnd;
                                                     ++TriggerTowerIterator) {

    //---------------------------- ADC HitMaps -------------------------------

    double crateModule      = (*TriggerTowerIterator)->eta();
    double submoduleChannel = (*TriggerTowerIterator)->phi();
    int icm       = crateModule;
    int isc       = submoduleChannel;
    int crate     = icm/16;
    int module    = icm%16;
    int submodule = isc/4;
    int channel   = isc%4;

    if (crate < 2 || crate > 5) continue;
    
    const int adc = (*TriggerTowerIterator)->emADC()[
                                         (*TriggerTowerIterator)->emADCPeak()];
    if (adc > m_TT_ADC_HitMap_Thresh) {
      m_h_TT_HitMap_ADC->Fill(crateModule-32., submoduleChannel, 1);
    }
    m_p_TT_HitMap_ADC->Fill(crateModule-32., submoduleChannel, adc);

    //------------------------ SubStatus Word errors -------------------------

    using LVL1::DataError;
    DataError error((*TriggerTowerIterator)-> emError());
   
    //Summary

    int ypos = module+(crate-2)*16;

    for (int bit = 0; bit < 8; ++bit) {
      if (error.get(bit + DataError::ChannelDisabled)) {
        m_h_fwPpmError_Crate_25->Fill(bit, ypos);
	m_histTool->fillEventNumber(m_h_TT_ASICEventNumbers, bit);
      }
      if (error.get(bit + DataError::GLinkParity)) {
        m_h_TT_error_Crate_25->Fill(bit, ypos);
	m_h_TT_Error->Fill(bit);
	m_histTool->fillEventNumber(m_h_TT_EventNumbers, bit);
      }
    }

    // Detailed plots by MCM
    ypos = (crate%2)*16+module;
    const int index = (crate-2)/2;
    if (error.get(DataError::ChannelDisabled)) {
      m_h_ErrorDetails[(channel/2)*4+index]->Fill(
                                               (channel%2)*16+submodule, ypos);
    }
    if (error.get(DataError::MCMAbsent)) {
      m_h_ErrorDetails[4+index]->Fill(submodule, ypos);
    }
    if (error.get(DataError::Timeout)) {
      m_h_ErrorDetails[6+index]->Fill(submodule, ypos);
    }
    if (error.get(DataError::ASICFull)) {
      m_h_ErrorDetails[6+index]->Fill(16+submodule, ypos);
    }
    if (error.get(DataError::EventMismatch)) {
      m_h_ErrorDetails[8+index]->Fill(submodule, ypos);
    }
    if (error.get(DataError::BunchMismatch)) {
      m_h_ErrorDetails[8+index]->Fill(16+submodule, ypos);
    }
    if (error.get(DataError::FIFOCorrupt)) {
      m_h_ErrorDetails[10+index]->Fill(submodule, ypos);
    }
    if (error.get(DataError::PinParity)) {
      m_h_ErrorDetails[10+index]->Fill(16+submodule, ypos);
    }

    if (error.get(DataError::ChannelDisabled) ||
        error.get(DataError::MCMAbsent)) overview[crate] |= 1;

    if (error.get(DataError::Timeout)       ||
        error.get(DataError::ASICFull)      ||
        error.get(DataError::EventMismatch) ||
	error.get(DataError::BunchMismatch) ||
        error.get(DataError::FIFOCorrupt)   ||
	error.get(DataError::PinParity)) overview[crate] |= (1 << 1);

    if (error.get(DataError::GLinkParity)   ||
        error.get(DataError::GLinkProtocol) ||
        error.get(DataError::FIFOOverflow)  ||
	error.get(DataError::ModuleError)   ||
        error.get(DataError::GLinkDown)     ||
	error.get(DataError::GLinkTimeout)  ||
        error.get(DataError::BCNMismatch)) overview[crate] |= (1 << 2);
     
     // number of triggered slice
     m_h_TT_triggeredSlice->Fill((*TriggerTowerIterator)->emADCPeak(),1);

  }	     
     
  // Write overview vector to StoreGate
  //std::vector<int>* save = new std::vector<int>(overview);
  //sc = evtStore()->record(save, "L1CaloPPMSpareErrorVector");
  //if (sc != StatusCode::SUCCESS)
  //  {
  //    msg(MSG::ERROR) << "Error recording PPMSpare error vector in TES "
  //                    << endreq;
  //    return sc;
  //  }
  
  return StatusCode::SUCCESS;
}


   
/*---------------------------------------------------------*/
StatusCode PPrSpareMon::procHistograms( bool isEndOfEventsBlock,
                                        bool isEndOfLumiBlock,
					bool isEndOfRun )
/*---------------------------------------------------------*/
{
  msg(MSG::DEBUG) << "in procHistograms" << endreq ;

  if( isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun ) { }
	
  return StatusCode::SUCCESS;
}

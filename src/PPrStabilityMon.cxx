#include <cmath>
#include <utility>
#include <stdlib.h>

#include "TH1F.h"
#include "TString.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"
#include "SGTools/StlVectorClids.h"

#include "LWHists/LWHist.h"
#include "LWHists/TH1F_LW.h"
#include "LWHists/TH2F_LW.h"
#include "LWHists/TH2I_LW.h"
#include "LWHists/TProfile_LW.h"
#include "LWHists/TProfile2D_LW.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

#include "TrigT1CaloMonitoring/PPrStabilityMon.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"
#include "TrigT1CaloMonitoringTools/TrigT1CaloLWHistogramTool.h"
#include "TrigT1CaloToolInterfaces/IL1TriggerTowerTool.h"
#include "TrigT1CaloCalibConditions/L1CaloCoolChannelId.h"
#include "TrigT1CaloCalibTools/L1CaloPprFineTimePlotManager.h"
#include "TrigConfigSvc/ILVL1ConfigSvc.h"

#include "TrigT1CaloEvent/TriggerTowerCollection.h"
#include "TrigT1CaloEvent/TriggerTower_ClassDEF.h"
#include "TrigT1CaloUtils/DataError.h"

using namespace std;

PPrStabilityMon::PPrStabilityMon(const std::string & type, const std::string & name, const IInterface* parent): ManagedMonitorToolBase ( type, name, parent ),
    m_ppmADCMinValue(0),
    m_lumiBlock(0),
    m_storeGate("StoreGateSvc",name),
    m_errorTool("TrigT1CaloMonErrorTool"),
    m_histTool("TrigT1CaloLWHistogramTool"),
    m_ttTool("LVL1::L1TriggerTowerTool/L1TriggerTowerTool"),
    m_plotManager(0),
    m_evtInfo(0)
{
  declareProperty("BS_TriggerTowerContainer",m_TriggerTowerContainerName = "LVL1TriggerTowers");
  declareProperty("ppmADCMinValue", m_ppmADCMinValue=60);
  declareProperty("PathInRootFile", m_PathInRootFile="L1Calo/PPrStabilityMon") ;
}

PPrStabilityMon::~PPrStabilityMon()
{
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*---------------------------------------------------------*/
StatusCode PPrStabilityMon::initialize()
/*---------------------------------------------------------*/
{
  msg(MSG::INFO) << "Initializing " << name() << " - package version "
                 << PACKAGE_VERSION << endreq;

  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  sc = m_errorTool.retrieve();
  if( sc.isFailure() ){msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloMonErrorTool" << endreq;return sc;}

  sc = m_histTool.retrieve();
  if( sc.isFailure() ) {msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloLWHistogramTool" << endreq; return sc;}

  sc = m_ttTool.retrieve();
  if( sc.isFailure() ) {msg(MSG::ERROR) << "Unable to locate Tool L1TriggerTowerTool" << endreq;return sc;}

  sc = m_storeGate.retrieve();
  if( sc.isFailure() ) {msg(MSG::ERROR) << "Unable to locate Tool StoreGateSvcTools "<< endreq; return sc;}

  m_plotManager= new L1CaloPprFineTimePlotManager(this,m_PathInRootFile,m_ppmADCMinValue);

  return StatusCode::SUCCESS;
}

StatusCode PPrStabilityMon::finalize()
{
    delete m_plotManager;
    return StatusCode::SUCCESS;
}

StatusCode PPrStabilityMon::fillHistograms()
{
    StatusCode sc;
    const bool debug = msgLvl(MSG::DEBUG);

    // Skip events believed to be corrupt
    if (m_errorTool->corrupt()){if (debug) msg(MSG::DEBUG) << "Skipping corrupt event" << endreq;return StatusCode::SUCCESS;}

    //Retrieve eventInfo from storeGate;
    m_evtInfo = 0;
    sc = m_storeGate->retrieve(m_evtInfo);
    if( sc.isFailure() ) { msg(MSG::ERROR) <<"Could not retrieve Event Info" <<endreq; return sc;}
   
    //Retrieve TriggerTowers from SG
    const TriggerTowerCollection* trigTwrColl = 0; 

    sc = evtStore()->retrieve(trigTwrColl,m_TriggerTowerContainerName); 
    if (sc.isFailure())
    {
        if (debug) msg(MSG::DEBUG) << "No TriggerTower found at "<< m_TriggerTowerContainerName << endreq ;
        return sc;
    }
    if (debug) msg(MSG::DEBUG)<<"In Fill histograms"<<endreq;

    // ================= Container: TriggerTower ===========================
    
    TriggerTowerCollection::const_iterator TriggerTowerIterator = trigTwrColl->begin(); 
    TriggerTowerCollection::const_iterator TriggerTowerIteratorEnd = trigTwrColl->end(); 
    
    for (; TriggerTowerIterator != TriggerTowerIteratorEnd;++TriggerTowerIterator) 
    {
        const double eta = (*TriggerTowerIterator)->eta();
        const double phi = (*TriggerTowerIterator)->phi();

        const L1CaloCoolChannelId emCoolChannelID = m_ttTool->channelID(eta,phi,0);
        const L1CaloCoolChannelId hadCoolChannelID = m_ttTool->channelID(eta,phi,1);

        unsigned int emcoolID = emCoolChannelID.id();
        unsigned int hadcoolID = hadCoolChannelID.id();

        m_plotManager->SetValues(emcoolID,hadcoolID,m_evtInfo, *TriggerTowerIterator);
    }
    return sc;
}

StatusCode PPrStabilityMon::procHistograms(bool /*isEndofEventsBlock*/, bool /*isEndofLumiBlock*/, bool isEndofRun)
{
    if(isEndofRun){m_plotManager->MakeSummary();}
    return StatusCode::SUCCESS;
}


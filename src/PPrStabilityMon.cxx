
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

#include "TrigT1CaloMonitoring/PPrStabilityMon.h"
#include "TrigT1CaloMonitoring/TrigT1CaloMonErrorTool.h"
#include "TrigT1CaloMonitoringTools/TrigT1CaloLWHistogramTool.h"
#include "TrigT1CaloToolInterfaces/IL1TriggerTowerTool.h"
#include "TrigT1CaloCalibConditions/L1CaloCoolChannelId.h"
#include "TrigT1CaloCalibTools/L1CaloPprFineTimePlotManager.h"
#include "TrigT1CaloCalibTools/L1CaloPprPedestalPlotManager.h"
#include "TrigT1CaloCalibTools/L1CaloPprEtCorrelationPlotManager.h"

#include "TrigT1CaloEvent/TriggerTowerCollection.h"

PPrStabilityMon::PPrStabilityMon(const std::string & type, const std::string & name, const IInterface* parent): ManagedMonitorToolBase ( type, name, parent ),
    m_ppmADCMinValue(0),
    m_lumiBlock(0),
    m_lumiBlockMax(0),
    m_errorTool("TrigT1CaloMonErrorTool"),
    m_histTool("TrigT1CaloLWHistogramTool"),
    m_ttTool("LVL1::L1TriggerTowerTool/L1TriggerTowerTool"),
    m_fineTimePlotManager(0),
    m_pedestalPlotManager(0),
    m_etCorrelationPlotManager(0),
    m_doFineTimeMonitoring(0),
    m_doPedestalMonitoring(0),
    m_doEtCorrelationMonitoring(0),												
    m_evtInfo(0),
    m_fineTimeCut(0),
    m_pedestalMaxWidth(0),
    m_caloCellContainerName(""),
    m_EtMinForEtCorrelation(0)
{
  declareProperty("BS_TriggerTowerContainer",m_TriggerTowerContainerName = "LVL1TriggerTowers");
  declareProperty("ppmADCMinValue", m_ppmADCMinValue=60);
  declareProperty("PathInRootFile", m_PathInRootFile="L1Calo/PPrStabilityMon");
  declareProperty("doFineTimeMonitoring",m_doFineTimeMonitoring = true );
  declareProperty("doPedestalMonitoring",m_doPedestalMonitoring = true );
  declareProperty("doEtCorrelationMonitoring",m_doEtCorrelationMonitoring = true );
  declareProperty("lumiMax", m_lumiBlockMax = 2000);
  declareProperty("fineTimeCut",m_fineTimeCut = 20 );
  declareProperty("pedestalMaxWidth",m_pedestalMaxWidth = 10 );
  declareProperty("caloCellContainerName",m_caloCellContainerName="AllCalo");
  declareProperty("EtMinForEtCorrelation",m_EtMinForEtCorrelation=5);
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
        
    if (m_doFineTimeMonitoring)
    {
        m_fineTimePlotManager = new L1CaloPprFineTimePlotManager(this,
								m_ttTool,
								m_ppmADCMinValue,
								m_lumiBlockMax,
								m_PathInRootFile+"/FineTime");
	m_fineTimePlotManager->SetFineTimeCut(m_fineTimeCut);
    }
    if (m_doPedestalMonitoring)
    {
        m_pedestalPlotManager = new L1CaloPprPedestalPlotManager(this,
								m_ttTool,
								m_lumiBlockMax,
								m_PathInRootFile+"/Pedestal");
	m_pedestalPlotManager->SetPedestalMaxWidth(m_pedestalMaxWidth);
    }
    if (m_doEtCorrelationMonitoring){
        m_etCorrelationPlotManager = new L1CaloPprEtCorrelationPlotManager(this,
								m_ttTool,
								m_lumiBlockMax,
								m_PathInRootFile+"/EtCorrelation");
	m_etCorrelationPlotManager->SetCaloCellContainer(m_caloCellContainerName);
	m_etCorrelationPlotManager->SetEtMin(m_EtMinForEtCorrelation);
    }
    
    return StatusCode::SUCCESS;
    
}

StatusCode PPrStabilityMon::finalize()
{
    if (m_doFineTimeMonitoring) {delete m_fineTimePlotManager;}
    if (m_doPedestalMonitoring) {delete m_pedestalPlotManager;}
    if (m_doEtCorrelationMonitoring) {delete m_etCorrelationPlotManager;}
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
    sc = evtStore()->retrieve(m_evtInfo);
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
    
    if (m_doEtCorrelationMonitoring) {
        sc = m_etCorrelationPlotManager->getCaloCells();
	if (sc.isFailure()) return sc;
    }
    
    // ================= Container: TriggerTower ===========================
    
    TriggerTowerCollection::const_iterator TriggerTowerIterator = trigTwrColl->begin(); 
    TriggerTowerCollection::const_iterator TriggerTowerIteratorEnd = trigTwrColl->end(); 
    
    for (; TriggerTowerIterator != TriggerTowerIteratorEnd; ++TriggerTowerIterator) 
    {
        const double eta = (*TriggerTowerIterator)->eta();
        const double phi = (*TriggerTowerIterator)->phi();

        const L1CaloCoolChannelId emCoolChannelID = m_ttTool->channelID(eta,phi,0);
        const L1CaloCoolChannelId hadCoolChannelID = m_ttTool->channelID(eta,phi,1);

        bool emDead = m_ttTool->disabledChannel(emCoolChannelID);
        bool hadDead= m_ttTool->disabledChannel(hadCoolChannelID);
		
        if (m_doFineTimeMonitoring) {
	    m_fineTimePlotManager->Analyze(m_evtInfo, *TriggerTowerIterator,emDead,hadDead);
	}
	if (m_doPedestalMonitoring) {
	    m_pedestalPlotManager->Analyze(m_evtInfo, *TriggerTowerIterator,emDead,hadDead);
	}
	if (m_doEtCorrelationMonitoring) {
	    m_etCorrelationPlotManager->Analyze(m_evtInfo, *TriggerTowerIterator,emDead,hadDead);
	}
	
    }
    
    return sc;
}

StatusCode PPrStabilityMon::procHistograms(bool /*isEndofEventsBlock*/, bool /*isEndofLumiBlock*/, bool isEndofRun)
{
    //if(isEndofRun){m_fineTimePlotManager->MakeSummary();}
    if(isEndofRun){}
    return StatusCode::SUCCESS;
}


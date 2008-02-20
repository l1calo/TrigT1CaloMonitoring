// ********************************************************************
//
// NAME:     CPMSimBSMon.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Peter Faulkner
//           
//
// ********************************************************************

#include <sstream>
#include <utility>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "StoreGate/StoreGateSvc.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigT1Calo/CPAlgorithm.h"
#include "TrigT1Calo/CMMCPHits.h"
#include "TrigT1Calo/CPMHits.h"
#include "TrigT1Calo/CPMTower.h"
#include "TrigT1Calo/CPMRoI.h"
#include "TrigT1Calo/CoordToHardware.h"
#include "TrigT1Calo/TriggerTower.h"
#include "TrigT1Calo/TriggerTowerKey.h"
#include "TrigT1CaloTools/IL1EmTauTools.h"
#include "TrigT1CaloTools/IL1CPHitsTools.h"
#include "TrigT1Interfaces/Coordinate.h"
#include "TrigT1Interfaces/TrigT1CaloDefs.h"

#include "TrigT1CaloMonitoring/CPMSimBSMon.h"


/*---------------------------------------------------------*/
CPMSimBSMon::CPMSimBSMon(const std::string & type, 
			 const std::string & name,
			 const IInterface* parent)
  : ManagedMonitorToolBase(type, name, parent),
    m_storeGate("StoreGateSvc", name),
    m_emTauTool("LVL1::L1EmTauTools/L1EmTauTools"),
    m_cpHitsTool("LVL1::L1CPHitsTools/L1CPHitsTools"),
    m_monGroup(0)
/*---------------------------------------------------------*/
{
  declareInterface<IMonitorToolBase>(this); 

  declareProperty("CPMTowerLocation",
                 m_cpmTowerLocation  = LVL1::TrigT1CaloDefs::CPMTowerLocation);
  declareProperty("CPMHitsLocation",
                 m_cpmHitsLocation   = LVL1::TrigT1CaloDefs::CPMHitsLocation);
  declareProperty("CMMCPHitsLocation",
                 m_cmmCpHitsLocation = LVL1::TrigT1CaloDefs::CMMCPHitsLocation);
  declareProperty("CPMRoILocation",
                 m_cpmRoiLocation    = LVL1::TrigT1CaloDefs::CPMRoILocation);
  declareProperty("CPMRoILocationRoIB",
                 m_cpmRoiLocationRoib =
		                 LVL1::TrigT1CaloDefs::CPMRoILocation+"RoIB");
  declareProperty("TriggerTowerLocation",
                 m_triggerTowerLocation =
		                 LVL1::TrigT1CaloDefs::TriggerTowerLocation);

  declareProperty("RootDirectory", m_rootDir = "L1Calo");
}

/*---------------------------------------------------------*/
CPMSimBSMon::~CPMSimBSMon()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon:: initialize()
/*---------------------------------------------------------*/
{
  StatusCode sc;

  sc = ManagedMonitorToolBase::initialize();
  if (sc.isFailure()) return sc;

  MsgStream log(msgSvc(), name());
  
  sc = m_storeGate.retrieve();
  if( sc.isFailure() ) {
    log << MSG::ERROR << "Unable to locate Service StoreGateSvc" << endreq;
    return sc;
  }

  sc = m_emTauTool.retrieve();
  if( sc.isFailure() ) {
    log << MSG::ERROR << "Unable to locate Tool L1EmTauTools" << endreq;
    return sc;
  }

  sc = m_cpHitsTool.retrieve();
  if( sc.isFailure() ) {
    log << MSG::ERROR << "Unable to locate Tool L1CPHitsTools" << endreq;
    return sc;
  }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::bookHistograms(bool isNewEventsBlock,
                                           bool isNewLumiBlock, bool isNewRun)
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  if( m_environment == AthenaMonManager::online ) {
    // book histograms that are only made in the online environment...
  }
  	
  if( m_dataType == AthenaMonManager::cosmics ) {
    // book histograms that are only relevant for cosmics data...
  }

  MonGroup monShift ( this, m_rootDir + "/5_CP_Sim", shift, run );

  if ( isNewEventsBlock || isNewLumiBlock ) { }

  if( isNewRun ) {

  //  Error checks

  m_monGroup = &monShift;

  m_h_CPeqSIM = book2D("CPeqSIM",
            "CP Comparison with Simulation - Matches (Events);Crate/Module",
             64, 0, 64, 7, 0, 7);
  m_h_CPeqSIM->SetStats(kFALSE);
  setLabels(m_h_CPeqSIM);
  m_h_CPneSIM = book2D("CPneSIM",
            "CP Comparison with Simulation - Mismatches (Events);Crate/Module",
             64, 0, 64, 7, 0, 7);
  m_h_CPneSIM->SetStats(kFALSE);
  setLabels(m_h_CPneSIM);

  } // end if (isNewRun ...
  
  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::fillHistograms()
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());
  

  //Retrieve Trigger Towers from SG
  const TriggerTowerCollection* triggerTowerTES = 0; 
  StatusCode sc = m_storeGate->retrieve(triggerTowerTES,
                                                     m_triggerTowerLocation); 
  if( sc.isFailure()  ||  !triggerTowerTES ) {
    log << MSG::DEBUG<< "No Trigger Tower container found"<< endreq; 
  }

  //Retrieve CPM Towers from SG
  const CpmTowerCollection* cpmTowerTES = 0; 
  sc = m_storeGate->retrieve(cpmTowerTES, m_cpmTowerLocation); 
  if( sc.isFailure()  ||  !cpmTowerTES ) {
    log << MSG::DEBUG<< "No CPM Tower container found"<< endreq; 
  }
  
  //Retrieve CPM RoIs from SG
  const CpmRoiCollection* cpmRoiTES = 0;
  sc = m_storeGate->retrieve( cpmRoiTES, m_cpmRoiLocation);
  if( sc.isFailure()  ||  !cpmRoiTES ) {
    log << MSG::DEBUG << "No DAQ CPM RoIs container found, trying RoIB"
        << endreq; 
    sc = m_storeGate->retrieve( cpmRoiTES, m_cpmRoiLocationRoib);
    if( sc.isFailure()  ||  !cpmRoiTES ) {
      log << MSG::DEBUG << "No RoIB CPM RoIs container found"<< endreq;
    }
  }
  
  //Retrieve CPM Hits from SG
  const CpmHitsCollection* cpmHitsTES = 0;
  sc = m_storeGate->retrieve( cpmHitsTES, m_cpmHitsLocation);
  if( sc.isFailure()  ||  !cpmHitsTES ) {
    log << MSG::DEBUG << "No CPM Hits container found"<< endreq; 
  }
  
  //Retrieve CMM-CP Hits from SG
  const CmmCpHitsCollection* cmmCpHitsTES = 0;
  sc = m_storeGate->retrieve( cmmCpHitsTES, m_cmmCpHitsLocation);
  if( sc.isFailure()  ||  !cmmCpHitsTES ) {
    log << MSG::DEBUG << "No CMM-CP Hits container found"<< endreq; 
  }

  // Maps to simplify comparisons
  
  TriggerTowerMap ttMap;
  CpmTowerMap     cpMap;
  CpmRoiMap       crMap;
  CpmHitsMap      chMap;
  CmmCpHitsMap    cmMap;
  setupMap(triggerTowerTES, ttMap);
  setupMap(cpmTowerTES, cpMap);
  setupMap(cpmRoiTES, crMap);
  setupMap(cpmHitsTES, chMap);
  setupMap(cmmCpHitsTES, cmMap);

  // Compare Trigger Towers and CPM Towers from data

  const int nCrates = 4;
  const int nCPMs   = 14;
  const int vecsizeCpm = 2 * nCrates * nCPMs;
  ErrorVector errors1(vecsizeCpm);
  compare(ttMap, cpMap, errors1);

  // Compare RoIs simulated from CPM Towers with CPM RoIs from data

  CpmRoiCollection* cpmRoiSIM = 0;
  if (cpmTowerTES) {
    cpmRoiSIM = new CpmRoiCollection;
    simulate(cpMap, cpmRoiSIM);
  }
  CpmRoiMap crSimMap;
  setupMap(cpmRoiSIM, crSimMap);
  ErrorVector errors2(vecsizeCpm);
  compare(crSimMap, crMap, errors2);
  crSimMap.clear();
  delete cpmRoiSIM;

  // Compare CPM Hits simulated from CPM RoIs with CPM Hits from data

  CpmHitsCollection* cpmHitsSIM = 0;
  if (cpmRoiTES) {
    cpmHitsSIM = new CpmHitsCollection;
    simulate(cpmRoiTES, cpmHitsSIM);
  }
  CpmHitsMap chSimMap;
  setupMap(cpmHitsSIM, chSimMap);
  ErrorVector errors3(vecsizeCpm);
  compare(chSimMap, chMap, errors3);
  chSimMap.clear();
  delete cpmHitsSIM;

  // Compare CPM hits with CMM Hits from data

  const int nCMMs = 2;
  const int vecsizeCmm = 2 * nCrates * nCMMs;
  ErrorVector errors4(vecsizeCpm);
  ErrorVector errors5(vecsizeCmm);
  compare(chMap, cmMap, errors4, errors5);

  // Compare Local sums simulated from CMM Hits with Local sums from data

  CmmCpHitsCollection* cmmLocalSIM = 0;
  if (cmmCpHitsTES) {
    cmmLocalSIM = new CmmCpHitsCollection;
    simulate(cmmCpHitsTES, cmmLocalSIM, LVL1::CMMCPHits::LOCAL);
  }
  CmmCpHitsMap cmmLocalSimMap;
  setupMap(cmmLocalSIM, cmmLocalSimMap);
  ErrorVector errors6(vecsizeCmm);
  compare(cmmLocalSimMap, cmMap, errors6, LVL1::CMMCPHits::LOCAL);
  cmmLocalSimMap.clear();
  delete cmmLocalSIM;

  // Compare Local sums with Remote sums from data

  ErrorVector errors7(vecsizeCmm);
  compare(cmMap, cmMap, errors7, LVL1::CMMCPHits::REMOTE_0);

  // Compare Total sums simulated from Remote sums with Total sums from data

  CmmCpHitsCollection* cmmTotalSIM = 0;
  if (cmmCpHitsTES) {
    cmmTotalSIM = new CmmCpHitsCollection;
    simulate(cmmCpHitsTES, cmmTotalSIM, LVL1::CMMCPHits::TOTAL);
  }
  CmmCpHitsMap cmmTotalSimMap;
  setupMap(cmmTotalSIM, cmmTotalSimMap);
  ErrorVector errors8(vecsizeCmm);
  compare(cmmTotalSimMap, cmMap, errors8, LVL1::CMMCPHits::TOTAL);
  cmmTotalSimMap.clear();
  delete cmmTotalSIM;

  // Fill histograms

  const int cpmBins = nCrates * nCPMs;
  const int cmmBins = nCrates * nCMMs;
  for (int crate = 0; crate < nCrates; ++crate) {
    for (int module = 0; module < nCPMs; ++module) {
      int xBin = crate * nCPMs + module;
      int yBin = 0;
      m_h_CPeqSIM->Fill(xBin, yBin++, errors1[xBin]);
      m_h_CPeqSIM->Fill(xBin, yBin++, errors2[xBin]);
      m_h_CPeqSIM->Fill(xBin, yBin++, errors3[xBin]);
      m_h_CPeqSIM->Fill(xBin, yBin,   errors4[xBin]);
      if (module < nCMMs) {
        xBin = crate * nCMMs + module + cpmBins;
        int loc = crate * nCMMs + module;
        m_h_CPeqSIM->Fill(xBin, yBin++, errors5[loc]);
        m_h_CPeqSIM->Fill(xBin, yBin++, errors6[loc]);
        m_h_CPeqSIM->Fill(xBin, yBin++, errors7[loc]);
        m_h_CPeqSIM->Fill(xBin, yBin++, errors8[loc]);
        xBin = crate * nCPMs + module;
      }
      yBin = 0;
      m_h_CPneSIM->Fill(xBin, yBin++, errors1[xBin+cpmBins]);
      m_h_CPneSIM->Fill(xBin, yBin++, errors2[xBin+cpmBins]);
      m_h_CPneSIM->Fill(xBin, yBin++, errors3[xBin+cpmBins]);
      m_h_CPneSIM->Fill(xBin, yBin,   errors4[xBin+cpmBins]);
      if (module < nCMMs) {
        xBin = crate * nCMMs + module + cpmBins;
        int loc = crate * nCMMs + module;
        m_h_CPneSIM->Fill(xBin, yBin++, errors5[loc+cmmBins]);
        m_h_CPneSIM->Fill(xBin, yBin++, errors6[loc+cmmBins]);
        m_h_CPneSIM->Fill(xBin, yBin++, errors7[loc+cmmBins]);
        m_h_CPneSIM->Fill(xBin, yBin++, errors8[loc+cmmBins]);
      }
    }
  }

  return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode CPMSimBSMon::procHistograms(bool isEndOfEventsBlock,
                                  bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
  MsgStream log(msgSvc(), name());

  if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
  }

  return StatusCode::SUCCESS;
}

TH1D* CPMSimBSMon::book1D(std::string name, std::string title,
                                   int nx, double xmin, double xmax)
{
  TH1D *hist = new TH1D(TString(name), TString(title), nx, xmin, xmax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), this->name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< name << endreq;
  }
  
  return hist;
}

TH2D* CPMSimBSMon::book2D(std::string name, std::string title,
                                   int nx, double xmin, double xmax,  
	                           int ny, double ymin, double ymax)
{		
  TH2D *hist = new TH2D(TString(name), TString(title), nx, xmin, xmax,
                                                       ny, ymin, ymax);
  
  if (m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    MsgStream log(msgSvc(), this->name());
    log << MSG::WARNING << "Could not register histogram : " 
	<< name << endreq;
  }
  hist->SetOption("colz");
  
  return hist;
}

//  Compare Trigger Towers and CPM Towers

void CPMSimBSMon::compare(const TriggerTowerMap& ttMap,
                          const CpmTowerMap& cpMap, ErrorVector& errors)
{

  const int maxKey = 0x7fffffff;
  TriggerTowerMap::const_iterator ttMapIter    = ttMap.begin();
  TriggerTowerMap::const_iterator ttMapIterEnd = ttMap.end();
  CpmTowerMap::const_iterator     cpMapIter    = cpMap.begin();
  CpmTowerMap::const_iterator     cpMapIterEnd = cpMap.end();

  while (ttMapIter != ttMapIterEnd || cpMapIter != cpMapIterEnd) {

    int ttKey = maxKey;
    int cpKey = maxKey;
    int ttEm  = 0;
    int ttHad = 0;
    int cpEm  = 0;
    int cpHad = 0;
    double eta = 0.;
    double phi = 0.;

    if (ttMapIter != ttMapIterEnd) ttKey = ttMapIter->first;
    if (cpMapIter != cpMapIterEnd) cpKey = cpMapIter->first;

    if ((cpMapIter == cpMapIterEnd) || (cpKey > ttKey)) {

      // TriggerTower but no CPMTower

      const LVL1::TriggerTower* tt = ttMapIter->second;
      ttEm  = tt->emEnergy();
      ttHad = tt->hadEnergy();
      eta = tt->eta();
      phi = tt->phi();
      ++ttMapIter;

    } else if ((ttMapIter == ttMapIterEnd) || (ttKey > cpKey)) {

      // CPMTower but no TriggerTower

      const LVL1::CPMTower* cp = cpMapIter->second;
      cpEm  = cp->emEnergy();
      cpHad = cp->hadEnergy();
      eta = cp->eta();
      phi = cp->phi();
      ++cpMapIter;

    } else {

      // Have both

      const LVL1::TriggerTower* tt = ttMapIter->second;
      const LVL1::CPMTower*     cp = cpMapIter->second;
      ttEm  = tt->emEnergy();
      ttHad = tt->hadEnergy();
      cpEm  = cp->emEnergy();
      cpHad = cp->hadEnergy();
      eta = tt->eta();
      phi = tt->phi();
      ++ttMapIter;
      ++cpMapIter;
    }
    
    //  Fill in error vector

    const LVL1::Coordinate coord(phi, eta);
    LVL1::CoordToHardware converter;
    const int crate = converter.cpCrate(coord);
    const int cpm   = converter.cpModule(coord);
    const int loc   = crate * 14 + cpm - 1;
    const int cpmBins = 4 * 14;
    if ((ttEm && ttEm == cpEm) || (ttHad && ttHad == cpHad)) errors[loc] = 1;
    if (ttEm != cpEm || ttHad != cpHad) errors[loc+cpmBins] = 1;
  }
}

//  Compare Simulated RoIs with data

void CPMSimBSMon::compare(const CpmRoiMap& roiSimMap, const CpmRoiMap& roiMap,
                                                      ErrorVector& errors)
{
  MsgStream log(msgSvc(), name());

  const int maxKey = 0xffff;
  CpmRoiMap::const_iterator simMapIter    = roiSimMap.begin();
  CpmRoiMap::const_iterator simMapIterEnd = roiSimMap.end();
  CpmRoiMap::const_iterator datMapIter    = roiMap.begin();
  CpmRoiMap::const_iterator datMapIterEnd = roiMap.end();

  while (simMapIter != simMapIterEnd || datMapIter != datMapIterEnd) {

    int simKey = maxKey;
    int datKey = maxKey;
    unsigned int simHits = 0;
    unsigned int datHits = 0;
    int crate = 0;
    int cpm   = 0;

    if (simMapIter != simMapIterEnd) simKey = simMapIter->first;
    if (datMapIter != datMapIterEnd) datKey = datMapIter->first;

    if ((datMapIter == datMapIterEnd) || (datKey > simKey)) {

      // Simulated RoI but no data RoI

      const LVL1::CPMRoI* roi = simMapIter->second;
      simHits = roi->hits();
      crate   = roi->crate();
      cpm     = roi->cpm();
      ++simMapIter;

    } else if ((simMapIter == simMapIterEnd) || (simKey > datKey)) {

      // Data RoI but no simulated RoI

      const LVL1::CPMRoI* roi = datMapIter->second;
      datHits = roi->hits();
      crate   = roi->crate();
      cpm     = roi->cpm();
      ++datMapIter;

    } else {

      // Have both

      const LVL1::CPMRoI* roiS = simMapIter->second;
      const LVL1::CPMRoI* roiD = datMapIter->second;
      simHits = roiS->hits();
      datHits = roiD->hits();
      crate   = roiD->crate();
      cpm     = roiD->cpm();
      ++simMapIter;
      ++datMapIter;
    }
    
    //  Fill in error vector

    const int loc = crate * 14 + cpm - 1;
    const int cpmBins = 4 * 14;
    if (simHits && simHits == datHits) errors[loc] = 1;
    if (simHits != datHits) errors[loc+cpmBins] = 1;
    log << MSG::DEBUG << "DataHits/SimHits: ";
    for (int i = 15; i >= 0; --i) {
      int bit = (datHits >> i) & 0x1;
      log << MSG::DEBUG << bit;
    }
    log << MSG::DEBUG << "/";
    for (int i = 15; i >= 0; --i) {
      int bit = (simHits >> i) & 0x1;
      log << MSG::DEBUG << bit;
    }
    log << MSG::DEBUG << endreq;
  }
}

//  Compare simulated CPM Hits with data

void CPMSimBSMon::compare(const CpmHitsMap& cpmSimMap, const CpmHitsMap& cpmMap,
                                                       ErrorVector& errors)
{

  const int maxKey = 0x7fffffff;
  CpmHitsMap::const_iterator simMapIter    = cpmSimMap.begin();
  CpmHitsMap::const_iterator simMapIterEnd = cpmSimMap.end();
  CpmHitsMap::const_iterator datMapIter    = cpmMap.begin();
  CpmHitsMap::const_iterator datMapIterEnd = cpmMap.end();

  while (simMapIter != simMapIterEnd || datMapIter != datMapIterEnd) {

    int simKey = maxKey;
    int datKey = maxKey;
    unsigned int simHits0 = 0;
    unsigned int simHits1 = 0;
    unsigned int datHits0 = 0;
    unsigned int datHits1 = 0;
    int crate = 0;
    int cpm   = 0;

    if (simMapIter != simMapIterEnd) simKey = simMapIter->first;
    if (datMapIter != datMapIterEnd) datKey = datMapIter->first;

    if ((datMapIter == datMapIterEnd) || (datKey > simKey)) {

      // Simulation Hits but no data Hits

      const LVL1::CPMHits* simh = simMapIter->second;
      simHits0 = simh->HitWord0();
      simHits1 = simh->HitWord1();
      crate    = simh->crate();
      cpm      = simh->module();
      ++simMapIter;

    } else if ((simMapIter == simMapIterEnd) || (simKey > datKey)) {

      // Data Hits but no simulation Hits

      const LVL1::CPMHits* dath = datMapIter->second;
      datHits0 = dath->HitWord0();
      datHits1 = dath->HitWord1();
      crate    = dath->crate();
      cpm      = dath->module();
      ++datMapIter;

    } else {

      // Have both

      const LVL1::CPMHits* simh = simMapIter->second;
      const LVL1::CPMHits* dath = datMapIter->second;
      simHits0 = simh->HitWord0();
      simHits1 = simh->HitWord1();
      datHits0 = dath->HitWord0();
      datHits1 = dath->HitWord1();
      crate    = dath->crate();
      cpm      = dath->module();
      ++simMapIter;
      ++datMapIter;
    }
    
    //  Fill in error vector

    const int loc = crate * 14 + cpm - 1;
    const int cpmBins = 4 * 14;
    if ((simHits0 && simHits0 == datHits0) ||
        (simHits1 && simHits1 == datHits1)) errors[loc] = 1;
    if (simHits0 != datHits0 || simHits1 != datHits1) errors[loc+cpmBins] = 1;
  }
}

//  Compare CPM Hits and CMM Hits

void CPMSimBSMon::compare(const CpmHitsMap& cpmMap, const CmmCpHitsMap& cmmMap,
                          ErrorVector& errors, ErrorVector& errors2)
{

  const int maxKey = 0x7fffffff;
  CpmHitsMap::const_iterator   cpmMapIter    = cpmMap.begin();
  CpmHitsMap::const_iterator   cpmMapIterEnd = cpmMap.end();
  CmmCpHitsMap::const_iterator cmmMapIter    = cmmMap.begin();
  CmmCpHitsMap::const_iterator cmmMapIterEnd = cmmMap.end();

  while (cpmMapIter != cpmMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    int cpmKey = maxKey;
    int cmmKey = maxKey;
    unsigned int cpmHits0 = 0;
    unsigned int cpmHits1 = 0;
    unsigned int cmmHits0 = 0;
    unsigned int cmmHits1 = 0;
    int crate  = 0;
    int cpm    = 0;

    if (cpmMapIter != cpmMapIterEnd) cpmKey = cpmMapIter->first;
    if (cmmMapIter != cmmMapIterEnd) cmmKey = cmmMapIter->first;

    if ((cmmMapIter == cmmMapIterEnd) || (cmmKey > cpmKey)) {

      // CPM Hits but no CMM Hits

      const LVL1::CPMHits* cpmh = cpmMapIter->second;
      cpmHits0 = cpmh->HitWord0();
      cpmHits1 = cpmh->HitWord1();
      crate    = cpmh->crate();
      cpm      = cpmh->module();
      ++cpmMapIter;

    } else if ((cpmMapIter == cpmMapIterEnd) || (cpmKey > cmmKey)) {

      // CMM Hits but no CPM Hits

      const LVL1::CMMCPHits* cmmh = cmmMapIter->second;
      cmmHits0 = cmmh->HitWord0();
      cmmHits1 = cmmh->HitWord1();
      crate    = cmmh->crate();
      cpm      = cmmh->dataID();
      ++cmmMapIter;
      if (cpm > 14) continue;

    } else {

      // Have both

      const LVL1::CPMHits*   cpmh = cpmMapIter->second;
      const LVL1::CMMCPHits* cmmh = cmmMapIter->second;
      cpmHits0 = cpmh->HitWord0();
      cpmHits1 = cpmh->HitWord1();
      cmmHits0 = cmmh->HitWord0();
      cmmHits1 = cmmh->HitWord1();
      crate    = cpmh->crate();
      cpm      = cpmh->module();
      ++cpmMapIter;
      ++cmmMapIter;
    }
    
    //  Fill in error vectors

    const int loc  = crate * 14 + cpm - 1;
    const int loc2 = crate * 2;
    const int cpmBins = 4 * 14;
    const int cmmBins = 4 * 2;
    if ((cpmHits0 && cpmHits0 == cmmHits0) ||
        (cpmHits1 && cpmHits1 == cmmHits1)) errors[loc] = 1;
    if (cpmHits0 != cmmHits0 || cpmHits1 != cmmHits1) errors[loc+cpmBins] = 1;
    if (cpmHits1 && cpmHits1 == cmmHits1) errors2[loc2]   = 1;
    if (cpmHits0 && cpmHits0 == cmmHits0) errors2[loc2+1] = 1;
    if (cpmHits1 != cmmHits1) errors2[loc2+cmmBins]   = 1; // hits1==>cmm 0
    if (cpmHits0 != cmmHits0) errors2[loc2+cmmBins+1] = 1; // hits0==>cmm 1
  }
}

//  Compare Simulated CMM Hits and Data CMM Hits

void CPMSimBSMon::compare(const CmmCpHitsMap& cmmSimMap,
                          const CmmCpHitsMap& cmmMap,
                          ErrorVector& errors, int selection)
{

  const bool local  = (selection == LVL1::CMMCPHits::LOCAL);
  const bool remote = (selection == LVL1::CMMCPHits::REMOTE_0);
  const bool total  = (selection == LVL1::CMMCPHits::TOTAL);
  if (!local && !remote && !total) return;
  std::vector<unsigned int> hits0Sim(3);
  std::vector<unsigned int> hits1Sim(3);
  std::vector<unsigned int> hits0(3);
  std::vector<unsigned int> hits1(3);
  const int maxKey = 0x7fffffff;
  CmmCpHitsMap::const_iterator cmmSimMapIter    = cmmSimMap.begin();
  CmmCpHitsMap::const_iterator cmmSimMapIterEnd = cmmSimMap.end();
  CmmCpHitsMap::const_iterator cmmMapIter       = cmmMap.begin();
  CmmCpHitsMap::const_iterator cmmMapIterEnd    = cmmMap.end();

  while (cmmSimMapIter != cmmSimMapIterEnd || cmmMapIter != cmmMapIterEnd) {

    int cmmSimKey = maxKey;
    int cmmKey    = maxKey;
    unsigned int cmmSimHits0 = 0;
    unsigned int cmmSimHits1 = 0;
    unsigned int cmmHits0 = 0;
    unsigned int cmmHits1 = 0;
    int crate  = 0;
    int dataId = 0;

    if (cmmSimMapIter != cmmSimMapIterEnd) cmmSimKey = cmmSimMapIter->first;
    if (cmmMapIter    != cmmMapIterEnd)    cmmKey    = cmmMapIter->first;

    if ((cmmMapIter == cmmMapIterEnd) || (cmmKey > cmmSimKey)) {

      // Sim CMM Hits but no Data CMM Hits

      const LVL1::CMMCPHits* cmmS = cmmSimMapIter->second;
      ++cmmSimMapIter;
      dataId = cmmS->dataID();
      if (local  && dataId != LVL1::CMMCPHits::LOCAL) continue;
      if (remote && dataId != LVL1::CMMCPHits::LOCAL) continue;
      if (total  && dataId != LVL1::CMMCPHits::TOTAL) continue;
      cmmSimHits0 = cmmS->HitWord0();
      cmmSimHits1 = cmmS->HitWord1();
      crate       = cmmS->crate();

    } else if ((cmmSimMapIter == cmmSimMapIterEnd) || (cmmSimKey > cmmKey)) {

      // Data CMM Hits but no Sim CMM Hits

      const LVL1::CMMCPHits* cmmD = cmmMapIter->second;
      ++cmmMapIter;
      dataId   = cmmD->dataID();
      if (local  && dataId != LVL1::CMMCPHits::LOCAL)    continue;
      if (remote && dataId != LVL1::CMMCPHits::REMOTE_0 &&
                    dataId != LVL1::CMMCPHits::REMOTE_1 &&
		    dataId != LVL1::CMMCPHits::REMOTE_2) continue;
      if (total  && dataId != LVL1::CMMCPHits::TOTAL)    continue;
      cmmHits0 = cmmD->HitWord0();
      cmmHits1 = cmmD->HitWord1();
      crate    = cmmD->crate();

    } else {

      // Have both

      const LVL1::CMMCPHits* cmmS = cmmSimMapIter->second;
      const LVL1::CMMCPHits* cmmD = cmmMapIter->second;
      ++cmmSimMapIter;
      ++cmmMapIter;
      dataId   = cmmS->dataID();
      if (local  && dataId != LVL1::CMMCPHits::LOCAL)    continue;
      if (remote && dataId != LVL1::CMMCPHits::LOCAL    &&
                    dataId != LVL1::CMMCPHits::REMOTE_0 &&
                    dataId != LVL1::CMMCPHits::REMOTE_1 &&
		    dataId != LVL1::CMMCPHits::REMOTE_2) continue;
      if (total  && dataId != LVL1::CMMCPHits::TOTAL)    continue;
      cmmSimHits0 = cmmS->HitWord0();
      cmmSimHits1 = cmmS->HitWord1();
      cmmHits0    = cmmD->HitWord0();
      cmmHits1    = cmmD->HitWord1();
      crate       = cmmS->crate();
    }
    
    //  Fill in error vectors

    if (local || total) {
      const int loc = crate * 2;
      const int cmmBins = 4 * 2;
      if (cmmSimHits1 && cmmSimHits1 == cmmHits1) errors[loc]   = 1;
      if (cmmSimHits0 && cmmSimHits0 == cmmHits0) errors[loc+1] = 1;
      if (cmmSimHits1 != cmmHits1) errors[loc+cmmBins]   = 1; // hits1==>cmm 0
      if (cmmSimHits0 != cmmHits0) errors[loc+cmmBins+1] = 1; // hits0==>cmm 1
    } else {
      if (dataId == LVL1::CMMCPHits::LOCAL) {
        if (crate != 3) {
	  hits0Sim[crate] = cmmSimHits0;
	  hits1Sim[crate] = cmmSimHits1;
        }
      } else if (dataId == LVL1::CMMCPHits::REMOTE_0) {
        hits0[0] = cmmHits0;
        hits1[0] = cmmHits1;
      } else if (dataId == LVL1::CMMCPHits::REMOTE_1) {
        hits0[1] = cmmHits0;
        hits1[1] = cmmHits1;
      } else {
        hits0[2] = cmmHits0;
        hits1[2] = cmmHits1;
      }
    }
  }
  if (remote) {
    for (int crate = 0; crate < 3; ++crate) {
      const int loc = crate * 2;
      const int cmmBins = 4 * 2;
      if (hits1Sim[crate] && hits1Sim[crate] == hits1[crate]) errors[loc]   = 1;
      if (hits0Sim[crate] && hits0Sim[crate] == hits0[crate]) errors[loc+1] = 1;
      if (hits1Sim[crate] != hits1[crate]) errors[loc+cmmBins]   = 1;
      if (hits0Sim[crate] != hits0[crate]) errors[loc+cmmBins+1] = 1;
    }
  }
}

void CPMSimBSMon::setLabels(TH2* hist)
{
  hist->GetXaxis()->SetBinLabel(1, "CPM");
  hist->GetXaxis()->SetBinLabel(8, "0/8");
  hist->GetXaxis()->SetBinLabel(15, "1/1");
  hist->GetXaxis()->SetBinLabel(22, "1/8");
  hist->GetXaxis()->SetBinLabel(29, "2/1");
  hist->GetXaxis()->SetBinLabel(36, "2/8");
  hist->GetXaxis()->SetBinLabel(43, "3/1");
  hist->GetXaxis()->SetBinLabel(50, "3/8");
  hist->GetXaxis()->SetBinLabel(57, "CMM");
  hist->GetXaxis()->SetBinLabel(59, "1/0");
  hist->GetXaxis()->SetBinLabel(61, "2/0");
  hist->GetXaxis()->SetBinLabel(63, "3/0");
  //hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.25);
  int bin = 1;
  hist->GetYaxis()->SetBinLabel(bin++, "Towers");
  hist->GetYaxis()->SetBinLabel(bin++, "RoIs");
  hist->GetYaxis()->SetBinLabel(bin++, "CPMHits");
  hist->GetYaxis()->SetBinLabel(bin++, "CMMHits");
  hist->GetYaxis()->SetBinLabel(bin++, "Local");
  hist->GetYaxis()->SetBinLabel(bin++, "Remote");
  hist->GetYaxis()->SetBinLabel(bin++, "Total");
  hist->GetYaxis()->SetLabelSize(0.045);
}

void CPMSimBSMon::setupMap(const TriggerTowerCollection* coll,
                                 TriggerTowerMap& map)
{
  if (coll) {
    LVL1::TriggerTowerKey towerKey;
    TriggerTowerCollection::const_iterator pos  = coll->begin();
    TriggerTowerCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const double eta = (*pos)->eta();
      if (eta > -2.5 && eta < 2.5) {
        const double phi = (*pos)->phi();
        const int key = towerKey.ttKey(phi, eta);
        map.insert(std::make_pair(key, *pos));
      }
    }
  }
}

void CPMSimBSMon::setupMap(const CpmTowerCollection* coll, CpmTowerMap& map)
{
  if (coll) {
    LVL1::TriggerTowerKey towerKey;
    CpmTowerCollection::const_iterator pos  = coll->begin();
    CpmTowerCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const double eta = (*pos)->eta();
      const double phi = (*pos)->phi();
      const int key = towerKey.ttKey(phi, eta);
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void CPMSimBSMon::setupMap(const CpmRoiCollection* coll, CpmRoiMap& map)
{
  if (coll) {
    CpmRoiCollection::const_iterator pos  = coll->begin();
    CpmRoiCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate = (*pos)->crate();
      const int cpm   = (*pos)->cpm();
      const int chip  = (*pos)->chip();
      const int loc   = (*pos)->location();
      const int key   = (((((crate << 4) | cpm) << 3) | chip) << 3) | loc;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void CPMSimBSMon::setupMap(const CpmHitsCollection* coll, CpmHitsMap& map)
{
  if (coll) {
    CpmHitsCollection::const_iterator pos  = coll->begin();
    CpmHitsCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate = (*pos)->crate();
      const int cpm   = (*pos)->module();
      const int key   = crate * 100 + cpm;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void CPMSimBSMon::setupMap(const CmmCpHitsCollection* coll, CmmCpHitsMap& map)
{
  if (coll) {
    CmmCpHitsCollection::const_iterator pos  = coll->begin();
    CmmCpHitsCollection::const_iterator posE = coll->end();
    for (; pos != posE; ++pos) {
      const int crate  = (*pos)->crate();
      const int dataId = (*pos)->dataID();
      const int key  = crate * 100 + dataId;
      map.insert(std::make_pair(key, *pos));
    }
  }
}

void CPMSimBSMon::simulate(const CpmTowerMap towers,
                                 CpmRoiCollection* rois)
{
  InternalRoiCollection* intRois = new InternalRoiCollection;
  m_emTauTool->findRoIs(&towers, intRois);
  m_cpHitsTool->formCPMRoI(intRois, rois);
  delete intRois;
}

void CPMSimBSMon::simulate(const CpmRoiCollection* rois,
                                 CpmHitsCollection* hits)
{
  m_cpHitsTool->formCPMHits(rois, hits);
}

void CPMSimBSMon::simulate(const CmmCpHitsCollection* hitsIn,
                                 CmmCpHitsCollection* hitsOut, int selection)
{
  if (selection == LVL1::CMMCPHits::LOCAL) {
    m_cpHitsTool->formCMMCPHitsCrate(hitsIn, hitsOut);
  } else if (selection == LVL1::CMMCPHits::TOTAL) {
    m_cpHitsTool->formCMMCPHitsSystem(hitsIn, hitsOut);
  }
}

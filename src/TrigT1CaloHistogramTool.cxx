#include <sstream>

#include "TH1F.h"
#include "TH2F.h"

#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

#include "TrigConfL1Data/CTPConfig.h"
#include "TrigConfL1Data/L1DataDef.h"
#include "TrigConfL1Data/Menu.h"
#include "TrigConfL1Data/TriggerThreshold.h"
#include "TrigT1CaloUtils/DataError.h"

#include "TrigT1CaloMonitoring/TrigT1CaloHistogramTool.h"

// Interface ID

static const InterfaceID IID_ITrigT1CaloHistogramTool(
                                           "TrigT1CaloHistogramTool", 1, 1);

const InterfaceID& TrigT1CaloHistogramTool::interfaceID()
{
  return IID_ITrigT1CaloHistogramTool;
}

// Constructor

TrigT1CaloHistogramTool::TrigT1CaloHistogramTool(const std::string& type,
                                                     const std::string& name,
	    			                     const IInterface*  parent)
        : AthAlgTool(type, name, parent),
	  m_configSvc("TrigConf::TrigConfigSvc/TrigConfigSvc", name),
	  m_monGroup(0), m_phiScaleTT(32./M_PI), m_phiScaleJE(16./M_PI)
{
  declareInterface<TrigT1CaloHistogramTool>(this);

  declareProperty( "LVL1ConfigSvc", m_configSvc, "LVL1 Config Service");
}

// Destructor

TrigT1CaloHistogramTool::~TrigT1CaloHistogramTool()
{
}

// Initialize

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

StatusCode TrigT1CaloHistogramTool::initialize()
{
  msg(MSG::INFO) << "Initializing " << name() << " - package version "
                 << PACKAGE_VERSION << endreq;

  // Connect to the TrigConfigSvc for the trigger configuration:

  StatusCode sc = m_configSvc.retrieve();
  if ( sc.isFailure() ) {
    msg(MSG::ERROR) << "Couldn't connect to " << m_configSvc.typeAndName()
                    << endreq;
    return sc;
  } else msg(MSG::INFO) << "Connected to " << m_configSvc.typeAndName()
                        << endreq;

  return StatusCode::SUCCESS;
}

// Finalize

StatusCode TrigT1CaloHistogramTool::finalize()
{
  return StatusCode::SUCCESS;
}

//===========================================================================
//  Labelling Utilities - General
//===========================================================================

// Get list of threshold names for given type

bool TrigT1CaloHistogramTool::thresholdNames(const std::string& type,
                                             std::vector<std::string>& names)
{
  bool found = false;
  names.clear();
  TrigConf::L1DataDef def;
  int nthresh = 0;
  if (type == def.emType() || type == def.tauType())
                                  nthresh = def.max_EM_Threshold_Number();
  else if (type == def.jetType()) nthresh = def.max_J_Threshold_Number();
  else if (type == def.jfType())  nthresh = def.max_JF_Threshold_Number();
  else if (type == def.jbType())  nthresh = def.max_JB_Threshold_Number();
  else if (type == def.xeType())  nthresh = def.max_XE_Threshold_Number();
  else if (type == def.jeType())  nthresh = def.max_JE_Threshold_Number();
  else if (type == def.teType())  nthresh = def.max_TE_Threshold_Number();
  else return found;

  for (int thr = 0; thr < nthresh; ++thr) {
    std::ostringstream cnum;
    cnum << thr;
    names.push_back(cnum.str());
  }

  const std::vector<TrigConf::TriggerThreshold*>& thresholds(m_configSvc->ctpConfig()->menu()->thresholdVector());
  std::vector<TrigConf::TriggerThreshold*>::const_iterator it;
  for (it = thresholds.begin(); it != thresholds.end(); ++it) {
    const std::string thrType((*it)->type());
    if (type == def.emType() || type == def.tauType()) {
      if (thrType != def.emType() && thrType != def.tauType()) continue;
    } else if (thrType != type) continue;
    const int threshNum = (*it)->thresholdNumber();
    if (threshNum < nthresh) {
      names[threshNum] = (*it)->name();
      found = true;
    }
  }
  return found;
}

// Label bins with threshold names

void TrigT1CaloHistogramTool::thresholds(TH1* hist, const std::string& type,
                                                    int offset, bool xAxis)
{
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  std::vector<std::string> names;
  bool found = thresholdNames(type, names);
  std::vector<std::string>::const_iterator it = names.begin();
  for (int bin = 1; it != names.end(); ++it, ++bin) {
    axis->SetBinLabel(bin+offset, (*it).c_str());
  }
  if (!found) axis->SetTitle("Threshold Number");
}

// Label bins with numbers

void TrigT1CaloHistogramTool::numbers(TH1* hist, int min, int max, int step,
                                                 int offset, bool xAxis)
{
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  int bin = 1 + offset;
  for (int num = min; num <= max; num += step) {
    std::ostringstream cnum;
    cnum << num;
    axis->SetBinLabel(bin, cnum.str().c_str());
    bin += step;
  }
}

// Label bins with number pairs

void TrigT1CaloHistogramTool::numberPairs(TH1* hist, int firstMin, int firstMax,
               int secondMin, int secondMax, int step, int offset, bool xAxis)
{
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  const int numSecond = secondMax - secondMin + 1;
  for (int first = firstMin; first <= firstMax; ++first) {
    int bin = 1 + offset + (first-firstMin)*numSecond;
    for (int second = secondMin; second <= secondMax; second += step) {
      std::ostringstream cnum;
      cnum << first << "/" << second;
      axis->SetBinLabel(bin, cnum.str().c_str());
      bin += step;
    }
  }
}

// Label bins with sub-status error bit names

void TrigT1CaloHistogramTool::subStatus(TH1* hist, bool xAxis)
{
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  const LVL1::DataError err(0);
  for (int bit = 0; bit < 8; ++bit) {
    axis->SetBinLabel(bit + 1,
          (err.bitName(bit + LVL1::DataError::GLinkParity)).c_str());
  }
  //axis->SetLabelSize(0.045);
}

//===========================================================================
//  Labelling Utilities - CPM
//===========================================================================

// Label bins with CPM crate/module

void TrigT1CaloHistogramTool::cpmCrateModule(TH1* hist, bool xAxis)
{
  const int nCrates = 4;
  const int nCPMs = 14;
  numberPairs(hist, 0, nCrates-1, 1, nCPMs, nCPMs/2, 0, xAxis);
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  axis->SetTitle("Crate/Module");
}

// Label bins with CPM RoI threshold names

void TrigT1CaloHistogramTool::cpmThresholds(TH1* hist, int offset, bool xAxis)
{
  thresholds(hist, TrigConf::L1DataDef::emType(), offset, xAxis);
}

//===========================================================================
//  Labelling Utilities - JEM
//===========================================================================

// Label bins with JEM crate/module

void TrigT1CaloHistogramTool::jemCrateModule(TH1* hist, bool xAxis)
{
  const int nCrates = 2;
  const int nJEMs = 16;
  numberPairs(hist, 0, nCrates-1, 0, nJEMs-1, 2, 0, xAxis);
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  axis->SetTitle("Crate/Module");
}

// Label bins with JEM RoI threshold names

void TrigT1CaloHistogramTool::jemThresholds(TH1* hist, int offset, bool xAxis)
{
  int newOffset = offset;
  thresholds(hist, TrigConf::L1DataDef::jetType(), newOffset, xAxis);
  newOffset += TrigConf::L1DataDef::max_J_Threshold_Number();
  thresholds(hist, TrigConf::L1DataDef::jbType(), newOffset, xAxis);
  newOffset += TrigConf::L1DataDef::max_JB_Threshold_Number();
  thresholds(hist, TrigConf::L1DataDef::jfType(), newOffset, xAxis);
}

// Label bins with Main Jet threshold names

void TrigT1CaloHistogramTool::mainJetThresholds(TH1* hist, int offset,
                                                                     bool xAxis)
{
  thresholds(hist, TrigConf::L1DataDef::jetType(), offset, xAxis);
}

// Label bins with Backward Jet threshold names

void TrigT1CaloHistogramTool::backwardJetThresholds(TH1* hist, int offset,
                                                                     bool xAxis)
{
  thresholds(hist, TrigConf::L1DataDef::jbType(), offset, xAxis);
}

// Label bins with Forward Jet threshold names

void TrigT1CaloHistogramTool::forwardJetThresholds(TH1* hist, int offset,
                                                                     bool xAxis)
{
  thresholds(hist, TrigConf::L1DataDef::jfType(), offset, xAxis);
}

// Label bins with JetEt threshold names

void TrigT1CaloHistogramTool::jetEtThresholds(TH1* hist, int offset, bool xAxis)
{
  thresholds(hist, TrigConf::L1DataDef::jeType(), offset, xAxis);
}

// Label bins with MissingEt threshold names

void TrigT1CaloHistogramTool::missingEtThresholds(TH1* hist, int offset,
                                                                     bool xAxis)
{
  thresholds(hist, TrigConf::L1DataDef::xeType(), offset, xAxis);
}

// Label bins with SumEt threshold names

void TrigT1CaloHistogramTool::sumEtThresholds(TH1* hist, int offset, bool xAxis)
{
  thresholds(hist, TrigConf::L1DataDef::teType(), offset, xAxis);
}

//===========================================================================
//  Booking Utilities - General
//===========================================================================

// Book and register a 1D histogram

TH1F* TrigT1CaloHistogramTool::book1F(const std::string& name,
                                      const std::string& title,
                                      int nx, double xmin, double xmax)
{
  TH1F *hist = new TH1F(name.c_str(), title.c_str(), nx, xmin, xmax);
  
  if (m_monGroup && m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    msg(MSG::WARNING) << "Could not register histogram : " 
	              << name << endreq;
  }
  hist->SetStats(kFALSE);
  
  return hist;
}

// Book and register a 2D histogram

TH2F* TrigT1CaloHistogramTool::book2F(const std::string& name,
                                      const std::string& title,
                                      int nx, double xmin, double xmax,  
	                              int ny, double ymin, double ymax)
{		
  TH2F *hist = new TH2F(name.c_str(), title.c_str(), nx, xmin, xmax,
                                                     ny, ymin, ymax);
  
  if (m_monGroup && m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    msg(MSG::WARNING) << "Could not register histogram : " 
	              << name << endreq;
  }
  hist->SetOption("colz");
  hist->SetStats(kFALSE);
  
  return hist;
}

// Book and register a 2D histogram with variable width bins

TH2F* TrigT1CaloHistogramTool::book2F(const std::string& name,
                                      const std::string& title,
                                      int nx, const double* xbins,
                                      int ny, double ymin, double ymax)
{
  TH2F *hist = new TH2F(name.c_str(), title.c_str(), nx, xbins,
                                                     ny, ymin, ymax);

  if (m_monGroup && m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    msg(MSG::WARNING) << "Could not register histogram : "
                      << name << endreq;
  }
  hist->SetOption("colz");
  hist->SetStats(kFALSE);

  return hist;
}

// Book and register a 2D histogram of integers displayed as text

TH2I* TrigT1CaloHistogramTool::book2I(const std::string& name,
                                      const std::string& title,
                                      int nx, double xmin, double xmax,
	                              int ny, double ymin, double ymax)
{
  TH2I *hist = new TH2I(name.c_str(), title.c_str(), nx, xmin, xmax,
                                                     ny, ymin, ymax);

  if (m_monGroup && m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    msg(MSG::WARNING) << "Could not register histogram : "
                      << name << endreq;
  }
  hist->SetOption("text");
  hist->SetStats(kFALSE);

  return hist;
}

//===========================================================================
//  Booking Utilities - CPM
//===========================================================================

// Book CPM crate/module vs thresholds

TH2F* TrigT1CaloHistogramTool::bookCPMCrateModuleVsThreshold(
                const std::string& name, const std::string& title)
{
  TH2F *hist = book2F(name, title, 56, 0, 56, 16, 0, 16);
  cpmCrateModule(hist);
  thresholds(hist, TrigConf::L1DataDef::emType(), 0, false);
  return hist;
}

// Book CPM eta vs phi

TH2F* TrigT1CaloHistogramTool::bookCPMEtaVsPhi(const std::string& name,
                                       const std::string& title, bool isRoi)
{
  // We use phi range 0-64 so tick marks correspond to the bins
  const double phiBin     = M_PI/32.;
  const double halfPhiBin = M_PI/64.;
  const double xmin = (isRoi) ? -2.45 : -2.5;
  const double xmax = (isRoi) ?  2.55 :  2.5;
  TH2F* hist = book2F(name, title, 50, xmin, xmax, 64, 0., 64.);
  hist->SetXTitle("eta");
  for (int chan = 0; chan < 64; chan += 4 ) {
    const double rad = (isRoi) ? (chan + 1)*phiBin : chan*phiBin + halfPhiBin;
    std::ostringstream cnum;
    cnum << chan << "/"
         << std::setiosflags(std::ios::fixed | std::ios::showpoint)
         << std::setprecision(2) << rad;
    hist->GetYaxis()->SetBinLabel(chan+1, cnum.str().c_str());
  }
  hist->GetYaxis()->SetBinLabel(64, "phi");
  return hist;
}

// Book CPM RoI eta vs phi

TH2F* TrigT1CaloHistogramTool::bookCPMRoIEtaVsPhi(const std::string& name,
                                                  const std::string& title)
{
  return bookCPMEtaVsPhi(name, title, true);
}

// Book CPM sub-status errors vs crate/module

TH2F* TrigT1CaloHistogramTool::bookCPMSubStatusVsCrateModule(
                          const std::string& name, const std::string& title)
{
  TH2F* hist = book2F(name, title, 8, 0, 8, 56, 0, 56);
  subStatus(hist);
  cpmCrateModule(hist, false);
  return hist;
}

//===========================================================================
//  Booking Utilities - JEM
//===========================================================================

// Book JEM eta vs phi

TH2F* TrigT1CaloHistogramTool::bookJEMEtaVsPhi(const std::string& name,
                                       const std::string& title, bool isRoi)
{
  // We use phi range 0-32 so tick marks correspond to the bins
  const double phiBin     = M_PI/16.;
  const double halfPhiBin = M_PI/32.;
  const int nxbins = 32;
  const double xbins[nxbins+1] = {-4.9,-3.2,-2.9,-2.7,-2.4,-2.2,-2.0,-1.8,-1.6,
                                  -1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,
				  0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,
				  2.7,2.9,3.2,4.9};
  const double xbinsRoi[nxbins+1] = {-4.0,-3.05,-2.8,-2.55,-2.3,-2.1,-1.9,-1.7,
                                     -1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,
				     0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,
				     2.1,2.3,2.55,2.8,3.05,4.0,4.95};
  TH2F* hist = (isRoi) ? book2F(name, title, nxbins, xbinsRoi, 32, 0., 32.)
                       : book2F(name, title, nxbins, xbins,    32, 0., 32.);
  hist->SetXTitle("eta");
  for (int chan = 0; chan < 32; chan += 2 ) {
    const double rad = (isRoi) ? (chan + 1)*phiBin : chan*phiBin + halfPhiBin;
    std::ostringstream cnum;
    cnum << chan << "/"
         << std::setiosflags(std::ios::fixed | std::ios::showpoint)
	 << std::setprecision(2) << rad;
    hist->GetYaxis()->SetBinLabel(chan+1, cnum.str().c_str());
  }
  hist->GetYaxis()->SetBinLabel(32, "phi");
  return hist;
}

// Book JEM RoI eta vs phi

TH2F* TrigT1CaloHistogramTool::bookJEMRoIEtaVsPhi(const std::string& name,
                                                  const std::string& title)
{
  return bookJEMEtaVsPhi(name, title, true);
}

//===========================================================================
//  Booking Utilities - PPM
//===========================================================================

// Book PPM Em eta vs phi

TH2F* TrigT1CaloHistogramTool::bookPPMEmEtaVsPhi(const std::string name,
                                                 const std::string title)
{
  // We use phi range 0-64 so tick marks correspond to the bins
  const double phiBin     = M_PI/32.;
  const double halfPhiBin = M_PI/64.;
  const int nxbins = 66;
  const double xbins[nxbins+1] = {-4.9,-4.475,-4.050,-3.625,-3.2,-3.1,-2.9,
                                  -2.7,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,
				  -1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,
				  -1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
				  -0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,
				  0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,
				  1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.7,2.9,
				  3.1,3.2,3.625,4.050,4.475,4.9};
  TH2F* hist = book2F(name, title, nxbins, xbins, 64, 0., 64.);
  hist->SetXTitle("eta");
  for (int chan = 0; chan < 64; chan += 4 ) {
    const double rad = chan*phiBin + halfPhiBin;
    std::ostringstream cnum;
    cnum << chan << "/"
         << std::setiosflags(std::ios::fixed | std::ios::showpoint)
         << std::setprecision(2) << rad;
    hist->GetYaxis()->SetBinLabel(chan+1, cnum.str().c_str());
  }
  hist->GetYaxis()->SetBinLabel(64, "phi");
  return hist;

}

// Book PPM Had eta vs phi

TH2F* TrigT1CaloHistogramTool::bookPPMHadEtaVsPhi(const std::string name,
                                                  const std::string title)
{
  return bookPPMEmEtaVsPhi(name, title);

}

//===========================================================================
//  Filling Utilities - CPM
//===========================================================================

// Fill CPM eta vs phi

void TrigT1CaloHistogramTool::fillCPMEtaVsPhi(TH2F* hist, double eta,
                                              double phi, double weight)
{
  const double phiMod = phi * m_phiScaleTT;
  hist->Fill(eta, phiMod, weight);
}

// Fill CPM RoI eta vs phi

void TrigT1CaloHistogramTool::fillCPMRoIEtaVsPhi(TH2F* hist, double eta,
                                                 double phi, double weight)
{
  const double phiMod = phi * m_phiScaleTT - 0.5;
  hist->Fill(eta, phiMod, weight);
}

//===========================================================================
//  Filling Utilities - JEM
//===========================================================================

// Fill JEM eta vs phi

void TrigT1CaloHistogramTool::fillJEMEtaVsPhi(TH2F* hist, double eta,
                                              double phi, double weight)
{
  const double phiMod = phi * m_phiScaleJE;
  if (eta < -3.2 || eta > 3.2) {
    // Fill two bins for FCAL
    hist->Fill(eta, phiMod - 0.5, weight);
    hist->Fill(eta, phiMod + 0.5, weight);
  } else hist->Fill(eta, phiMod, weight);
}

// Fill JEM RoI eta vs phi

void TrigT1CaloHistogramTool::fillJEMRoIEtaVsPhi(TH2F* hist, double eta,
                                                 double phi, double weight)
{
  const double phiMod = phi * m_phiScaleJE - 0.5;
  // JEPRoIDecoder returns eta=3.9 for both of the two forwardmost bins
  if (eta > 3.8 && eta < 4.0) hist->Fill(4.05, phiMod, weight);
  hist->Fill(eta, phiMod, weight);
}

//===========================================================================
//  Filling Utilities - PPM
//===========================================================================

// Fill PPM Em eta vs phi

void TrigT1CaloHistogramTool::fillPPMEmEtaVsPhi(TH2F* hist, double eta,
                                                double phi, double weight)
{
  const double phiMod = phi * m_phiScaleTT;
  double absEta = fabs(eta);
  if (absEta > 3.2) {
    // Fill four bins in phi
    hist->Fill(eta, phiMod + 1.5, weight);
    hist->Fill(eta, phiMod + 0.5, weight);
    hist->Fill(eta, phiMod - 0.5, weight);
    hist->Fill(eta, phiMod - 1.5, weight);
  } else if (absEta > 2.5) {
    // Fill two bins in phi
    hist->Fill(eta, phiMod + 0.5, weight);
    hist->Fill(eta, phiMod - 0.5, weight);
  } else hist->Fill(eta, phiMod, weight);
}

// Fill PPM Had eta vs phi

void TrigT1CaloHistogramTool::fillPPMHadEtaVsPhi(TH2F* hist, double eta,
                                                 double phi, double weight)
{
  // Use EM mapping - puts FCAL2 in left half of bin and FCAL3 in right half
  fillPPMEmEtaVsPhi(hist, eta, phi, weight);
}


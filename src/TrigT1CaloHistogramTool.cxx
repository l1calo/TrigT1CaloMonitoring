#include <sstream>

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "GaudiKernel/IInterface.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"

#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"

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

  declareProperty("EventSamples", m_eventSamples = 10,
                  "Number of Error Event Number Samples");
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

// Label bins with sub-status error bit names

void TrigT1CaloHistogramTool::subStatus(TH1* hist, int offset, bool xAxis)
{
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  const LVL1::DataError err(0);
  for (int bit = 0; bit < 8; ++bit) {
    axis->SetBinLabel(bit + 1 + offset,
          (err.bitName(bit + LVL1::DataError::GLinkParity)).c_str());
  }
}

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

//===========================================================================
//  Labelling Utilities - CPM
//===========================================================================

// Label bins with CPM chip/local coordinate

void TrigT1CaloHistogramTool::cpmChipLocalCoord(TH1* hist, int offset,
                                                                    bool xAxis)
{
  const int nChips = 8;
  const int nLoc   = 8;
  numberPairs(hist, 0, nChips-1, 0, nLoc-1, 4, offset, xAxis);
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  axis->SetTitle("Chip/Local Coord");
}

// Label bins with CPM/CMM crate/module

void TrigT1CaloHistogramTool::cpmCMMCrateModule(TH1* hist, int offset,
                                                                    bool xAxis)
{
  cpmCrateModule(hist, offset, xAxis);
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  axis->SetBinLabel(1+offset,  "CPM");
  axis->SetBinLabel(57+offset, "CMM");
  axis->SetBinLabel(59+offset, "1/0");
  axis->SetBinLabel(61+offset, "2/0");
  axis->SetBinLabel(63+offset, "3/0");
  //axis->SetTitleOffset(1.25);
}

// Label bins with CPM crate/CMM

void TrigT1CaloHistogramTool::cpmCrateCMM(TH1* hist, int offset, bool xAxis)
{
  const int nCrates = 4;
  const int nCMMs   = 2;
  numberPairs(hist, 0, nCrates-1, 0, nCMMs-1, 1, offset, xAxis);
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  axis->SetTitle("Crate/CMM");
}

// Label bins with CPM crate/module

void TrigT1CaloHistogramTool::cpmCrateModule(TH1* hist, int offset, bool xAxis)
{
  const int nCrates = 4;
  const int nCPMs = 14;
  numberPairs(hist, 0, nCrates-1, 1, nCPMs, nCPMs/2, offset, xAxis);
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

void TrigT1CaloHistogramTool::jemCrateModule(TH1* hist, int offset, bool xAxis)
{
  const int nCrates = 2;
  const int nJEMs = 16;
  numberPairs(hist, 0, nCrates-1, 0, nJEMs-1, 2, offset, xAxis);
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
//  Labelling Utilities - PPM
//===========================================================================

// Label bins with PPM crate/module

void TrigT1CaloHistogramTool::ppmCrateModule(TH1* hist, int firstCrate,
                                      int lastCrate, int offset, bool xAxis)
{
  int step = 2;
  if (lastCrate-firstCrate > 1) step = 4;
  numberPairs(hist, firstCrate, lastCrate, 0, 15, step, offset, xAxis);
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  axis->SetTitle("Crate/Module");
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
  registerHist(hist);
  hist->SetStats(kFALSE);
  return hist;
}

// Book and register a 1D histogram with variable width bins

TH1F* TrigT1CaloHistogramTool::book1F(const std::string& name,
                                      const std::string& title,
                                      int nx, const double* xbins)
{
  TH1F *hist = new TH1F(name.c_str(), title.c_str(), nx, xbins);
  registerHist(hist);
  hist->SetStats(kFALSE);
  return hist;
}

// Book and register a 1D profile histogram

TProfile* TrigT1CaloHistogramTool::bookProfile(const std::string& name,
                                               const std::string& title,
                                               int nx, double xmin, double xmax)
{
  TProfile *hist = new TProfile(name.c_str(), title.c_str(), nx, xmin, xmax);
  registerHist(hist);
  hist->SetStats(kFALSE);
  return hist;
}

// Book and register a 1D profile histogram with variable width bins

TProfile* TrigT1CaloHistogramTool::bookProfile(const std::string& name,
                                               const std::string& title,
                                               int nx, const double* xbins)
{
  TProfile *hist = new TProfile(name.c_str(), title.c_str(), nx, xbins);
  registerHist(hist);
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
  registerHist(hist);
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
  registerHist(hist);
  hist->SetOption("colz");
  hist->SetStats(kFALSE);
  return hist;
}

// Book and register a 2D profile histogram

TProfile2D* TrigT1CaloHistogramTool::bookProfile2D(const std::string& name,
                                                   const std::string& title,
                                            int nx, double xmin, double xmax,  
	                                    int ny, double ymin, double ymax)
{		
  TProfile2D *hist = new TProfile2D(name.c_str(), title.c_str(), nx, xmin, xmax,
                                                                 ny, ymin, ymax);
  registerHist(hist);
  hist->SetOption("colz");
  hist->SetStats(kFALSE);
  return hist;
}

// Book and register a 2D profile histogram with variable width bins

TProfile2D* TrigT1CaloHistogramTool::bookProfile2D(const std::string& name,
                                                   const std::string& title,
                                            int nx, const double* xbins,
                                            int ny, double ymin, double ymax)
{
  TProfile2D *hist = new TProfile2D(name.c_str(), title.c_str(), nx, xbins,
                                                                 ny, ymin, ymax);
  registerHist(hist);
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
  registerHist(hist);
  hist->SetOption("text");
  hist->SetStats(kFALSE);
  return hist;
}

// Book and register a 2D histogram containing event numbers as bin contents

TH2I* TrigT1CaloHistogramTool::bookEventNumbers(
                const std::string& name, const std::string& title,
		int ny, double ymin, double ymax)
{
  TH2I* hist = book2I(name, title, m_eventSamples, 0, m_eventSamples,
                                   ny, ymin, ymax);
  if (m_eventSamples <= 10) numbers(hist, 1, m_eventSamples);
  hist->GetXaxis()->SetTitle("Events with Error/Mismatch");
  return hist;
}

// Register a histogram

void TrigT1CaloHistogramTool::registerHist(TH1* hist)
{
  if (m_monGroup && m_monGroup->regHist(hist) != StatusCode::SUCCESS) {
    msg(MSG::WARNING) << "Could not register histogram : "
                      << hist->GetName() << endreq;
  }
}

//===========================================================================
//  Booking Utilities - CPM
//===========================================================================

// Book CPM crate/module vs chip/local coordinate

TH2F* TrigT1CaloHistogramTool::bookCPMCrateModuleVsChipLocalCoord(
                const std::string& name, const std::string& title)
{
  TH2F *hist = book2F(name, title, 56, 0., 56., 64, 0., 64.);
  cpmCrateModule(hist);
  cpmChipLocalCoord(hist, 0, false);
  return hist;
}

// Book CPM crate/module vs FPGA

TH2F* TrigT1CaloHistogramTool::bookCPMCrateModuleVsFPGA(
                const std::string& name, const std::string& title)
{
  TH2F *hist = book2F(name, title, 56, 0., 56., 20, 0., 20.);
  cpmCrateModule(hist);
  numbers(hist, 2, 17, 1, 2, false);
  // colour overlap fpga bins
  TAxis* axis = hist->GetYaxis();
  axis->SetBinLabel(1,  "#color[4]{0}");
  axis->SetBinLabel(2,  "#color[4]{1}");
  axis->SetBinLabel(19, "#color[4]{18}");
  axis->SetBinLabel(20, "#color[4]{19}");
  axis->SetTitle("Serialiser FPGA");
  return hist;
}

// Book CPM crate/module vs thresholds

TH2F* TrigT1CaloHistogramTool::bookCPMCrateModuleVsThreshold(
                const std::string& name, const std::string& title)
{
  TH2F *hist = book2F(name, title, 56, 0., 56., 16, 0., 16.);
  cpmCrateModule(hist);
  cpmThresholds(hist, 0, false);
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

// Book CPM events with error/mismatch vs crate/module

TH2I* TrigT1CaloHistogramTool::bookCPMEventVsCrateModule(
                const std::string& name, const std::string& title)
{
  TH2I* hist = bookEventNumbers(name, title, 56, 0., 56.);
  cpmCrateModule(hist, 0, false);
  return hist;
}

// Book CPM module vs crate

TH2F* TrigT1CaloHistogramTool::bookCPMModuleVsCrate(
                const std::string& name, const std::string& title)
{
  TH2F *hist = book2F(name, title, 14, 1., 15., 4, 0., 4.);
  numbers(hist, 1, 15);
  numbers(hist, 0, 4, 1, 0, false);
  hist->GetXaxis()->SetTitle("Module");
  hist->GetYaxis()->SetTitle("Crate");
  return hist;
}

// Book CPM module vs crate/CMM

TH2F* TrigT1CaloHistogramTool::bookCPMModuleVsCrateCMM(
                const std::string& name, const std::string& title)
{
  TH2F *hist = book2F(name, title, 14, 1., 15., 8, 0., 8.);
  numbers(hist, 1, 15);
  hist->GetXaxis()->SetTitle("Module");
  cpmCrateCMM(hist, 0, false);
  return hist;
}

// Book CPM sub-status errors vs crate/module

TH2F* TrigT1CaloHistogramTool::bookCPMSubStatusVsCrateModule(
                          const std::string& name, const std::string& title)
{
  TH2F* hist = book2F(name, title, 8, 0., 8., 56, 0., 56.);
  subStatus(hist);
  cpmCrateModule(hist, 0, false);
  return hist;
}

// Book CPM Sum/CMM

TH1F* TrigT1CaloHistogramTool::bookCPMSumCMM(const std::string& name,
                                             const std::string& title)
{
  TH1F* hist = book1F(name, title, 16, 0., 16.);
  TAxis* axis = hist->GetXaxis();
  axis->SetBinLabel(1, "L0/0");
  axis->SetBinLabel(2, "L0/1");
  axis->SetBinLabel(3, "L1/0");
  axis->SetBinLabel(4, "L1/1");
  axis->SetBinLabel(5, "L2/0");
  axis->SetBinLabel(6, "L2/1");
  axis->SetBinLabel(7, "L3/0");
  axis->SetBinLabel(8, "L3/1");
  axis->SetBinLabel(9, "R0/0");
  axis->SetBinLabel(10, "R0/1");
  axis->SetBinLabel(11, "R1/0");
  axis->SetBinLabel(12, "R1/1");
  axis->SetBinLabel(13, "R2/0");
  axis->SetBinLabel(14, "R2/1");
  axis->SetBinLabel(15, "T/0");
  axis->SetBinLabel(16, "T/1");
  axis->SetTitle("Sum/CMM");
  return hist;
}

// Book CPM Sum vs Threshold

TH2F* TrigT1CaloHistogramTool::bookCPMSumVsThreshold(const std::string& name,
                                                     const std::string& title)
{
  TH2F* hist = book2F(name, title, 8, 0., 8., 16, 0., 16.);
  TAxis* axis = hist->GetXaxis();
  axis->SetBinLabel(1, "L0");
  axis->SetBinLabel(2, "L1");
  axis->SetBinLabel(3, "L2");
  axis->SetBinLabel(4, "L3");
  axis->SetBinLabel(5, "R0");
  axis->SetBinLabel(6, "R1");
  axis->SetBinLabel(7, "R2");
  axis->SetBinLabel(8, "T");
  axis->SetTitle("Sum (Local/Remote/Total)");
  cpmThresholds(hist, 0, false);
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

// Book PPM Em eta vs phi profile

TProfile2D* TrigT1CaloHistogramTool::bookProfilePPMEmEtaVsPhi(
                              const std::string name, const std::string title)
{
  // todo - remove duplication with above
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
  TProfile2D* hist = bookProfile2D(name, title, nxbins, xbins, 64, 0., 64.);
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

// Book PPM Had eta vs phi profile

TProfile2D* TrigT1CaloHistogramTool::bookProfilePPMHadEtaVsPhi(
                              const std::string name, const std::string title)
{
  return bookProfilePPMEmEtaVsPhi(name, title);

}

// Book PPM events with error/mismatch vs crate/module

TH2I* TrigT1CaloHistogramTool::bookPPMEventVsCrateModule(
                const std::string& name, const std::string& title,
		int firstCrate, int lastCrate)
{
  int nbins = (lastCrate-firstCrate+1)*16;
  TH2I* hist = bookEventNumbers(name, title, nbins, 0, nbins);
  ppmCrateModule(hist, firstCrate, lastCrate, 0, false);
  return hist;
}

//===========================================================================
//  Filling Utilities - General
//===========================================================================

// Fill Error/Mismatch event number

void TrigT1CaloHistogramTool::fillEventNumber(TH2I* hist, double y)
{
  const int biny  = hist->GetYaxis()->FindBin(y);
  const int nbins = hist->GetNbinsX();
  for (int binx = 1; binx <= nbins; ++binx) {
    if (hist->GetBinContent(binx, biny) == 0.) {
      int eventNumber = 0;
      const EventInfo* evInfo = 0;
      StatusCode sc = evtStore()->retrieve(evInfo);
      if (sc.isFailure()) {
        msg(MSG::DEBUG) << "No EventInfo found" << endreq;
      } else {
        const EventID* evID = evInfo->event_ID();
        if (evID) eventNumber = evID->event_number();
      }
      hist->SetBinContent(binx, biny, eventNumber);
      break;
    }
  }
}

//===========================================================================
//  Filling Utilities - CPM
//===========================================================================

// Fill CPM eta vs phi

void TrigT1CaloHistogramTool::fillCPMEtaVsPhi(TH2* hist, double eta,
                                              double phi, double weight)
{
  const double phiMod = phi * m_phiScaleTT;
  hist->Fill(eta, phiMod, weight);
}

// Fill CPM RoI eta vs phi

void TrigT1CaloHistogramTool::fillCPMRoIEtaVsPhi(TH2* hist, double eta,
                                                 double phi, double weight)
{
  const double phiMod = phi * m_phiScaleTT - 0.5;
  hist->Fill(eta, phiMod, weight);
}

//===========================================================================
//  Filling Utilities - JEM
//===========================================================================

// Fill JEM eta vs phi

void TrigT1CaloHistogramTool::fillJEMEtaVsPhi(TH2* hist, double eta,
                                              double phi, double weight)
{
  const double phiMod = phi * m_phiScaleJE;
  if (eta < -3.2 || eta > 3.2) {
    // Fill two bins for FCAL
    hist->Fill(eta, phiMod + 0.5, weight);
    hist->Fill(eta, phiMod - 0.5, weight);
  } else hist->Fill(eta, phiMod, weight);
}

// Fill JEM RoI eta vs phi

void TrigT1CaloHistogramTool::fillJEMRoIEtaVsPhi(TH2* hist, double eta,
                                                 double phi, double weight)
{
  const double phiMod = phi * m_phiScaleJE - 0.5;
  // JEPRoIDecoder returns eta=3.9 for both of the two forwardmost bins
  if (eta > 3.8 && eta < 4.0) hist->Fill(4.05, phiMod, weight);
  hist->Fill(eta, phiMod, weight);
}

// Fill JEM phi allowing for granularity varying with eta

void TrigT1CaloHistogramTool::fillJEMPhi(TH1* hist, double eta, double phi,
                                                                double weight)
{
  const double halfBin = 1./(2.*m_phiScaleJE);
  if (eta < -3.2 || eta > 3.2) {
    // Fill two bins for FCAL
    hist->Fill(phi + halfBin, weight);
    hist->Fill(phi - halfBin, weight);
  } else hist->Fill(phi, weight);
}

//===========================================================================
//  Filling Utilities - PPM
//===========================================================================

// Fill PPM Em eta vs phi

void TrigT1CaloHistogramTool::fillPPMEmEtaVsPhi(TH2* hist, double eta,
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

void TrigT1CaloHistogramTool::fillPPMHadEtaVsPhi(TH2* hist, double eta,
                                                 double phi, double weight)
{
  // Use EM mapping - puts FCAL2 in left half of bin and FCAL3 in right half
  fillPPMEmEtaVsPhi(hist, eta, phi, weight);
}

// Fill PPM phi allowing for granularity varying with eta

void TrigT1CaloHistogramTool::fillPPMPhi(TH1* hist, double eta, double phi,
                                                                double weight)
{
  const double halfBin = 1./(2.*m_phiScaleTT);
  const double absEta  = fabs(eta);
  if (absEta > 3.2) {
    // Fill four bins in phi
    hist->Fill(phi + 3.*halfBin, weight);
    hist->Fill(phi + halfBin,    weight);
    hist->Fill(phi - halfBin,    weight);
    hist->Fill(phi - 3.*halfBin, weight);
  } else if (absEta > 2.5) {
    // Fill two bins in phi
    hist->Fill(phi + halfBin, weight);
    hist->Fill(phi - halfBin, weight);
  } else hist->Fill(phi, weight);
}

// Find bin in Em eta vs phi

int TrigT1CaloHistogramTool::findBinPPMEmEtaVsPhi(TH2* hist, double eta,
                                                             double phi)
{
  double phiMod = phi * m_phiScaleTT;
  if (eta < -2.5 || eta > 2.5) phiMod += 0.5;
  return hist->FindBin(eta, phiMod);
}

// Find bin in Had eta vs phi

int TrigT1CaloHistogramTool::findBinPPMHadEtaVsPhi(TH2* hist, double eta,
                                                              double phi)
{
  return findBinPPMEmEtaVsPhi(hist, eta, phi);
}

// Set bin content and optionally error in Em eta vs phi

void TrigT1CaloHistogramTool::setBinPPMEmEtaVsPhi(TH2* hist, int bin,
                                                  double content, double error)
{
  int binx, biny, binz;
  hist->GetBinXYZ(bin, binx, biny, binz);
  int nbin = 1;
  if      (binx <= 4 || binx >= 63) nbin = 4; // |eta|>3.2
  else if (binx <= 8 || binx >= 59) nbin = 2; // 2.5<|eta|<3.2
  int binyBase = ((biny-1)/nbin)*nbin+1;
  for (int i = 0; i < nbin; ++i) {
    hist->SetBinContent(binx, binyBase+i, content);
    if (error >= 0.) hist->SetBinError(binx, binyBase+i, error);
  }
}

// Set bin content and optionally error in Had eta vs phi

void TrigT1CaloHistogramTool::setBinPPMHadEtaVsPhi(TH2* hist, int bin,
                                                   double content, double error)
{
  setBinPPMEmEtaVsPhi(hist, bin, content, error);
}

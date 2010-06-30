#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>

#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
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
#include "TrigConfigSvc/ILVL1ConfigSvc.h"
#include "TrigT1CaloUtils/DataError.h"
#include "TrigT1CaloUtils/QuadLinear.h"

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
  declareProperty("ShrinkEtaBins", m_shrinkEtaBins = true,
                  "Make all eta bins the same size in eta/phi plots");
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

// Return int/double pair as a string

std::string TrigT1CaloHistogramTool::intDoubleString(int num, double val,
                                                              int precision)
{
  std::ostringstream cnum;
  cnum << num << "/" 
       << std::setiosflags(std::ios::fixed | std::ios::showpoint)
       << std::setprecision(precision) << val;
  return cnum.str();
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

// Label bins with number pairs without skipping bins when stepping

void TrigT1CaloHistogramTool::numberPairs2(TH1* hist,
               int firstMin, int firstMax,
               int secondMin, int secondMax, int step, int offset, bool xAxis)
{
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  int bin = 1 + offset;
  for (int first = firstMin; first <= firstMax; ++first) {
    for (int second = secondMin; second <= secondMax; second += step) {
      std::ostringstream cnum;
      cnum << first << "/" << second;
      axis->SetBinLabel(bin, cnum.str().c_str());
      bin++;
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

// Split long names for Y axis

std::string TrigT1CaloHistogramTool::splitLine(const std::string& word,
                                                                  bool xAxis)
{
  std::string newWord(word);
  if (!xAxis && word.length() > 6) {
    // split at last capital
    std::string::size_type idx =
                               word.find_last_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    if (idx != std::string::npos && idx != 0 && idx != word.length()-1) {
      newWord = "#splitline{" + word.substr(0, idx) + "}{"
                              + word.substr(idx) + "}";
    }
  }
  return newWord;
}

// Label bins with sub-status error bit names

void TrigT1CaloHistogramTool::subStatus(TH1* hist, int offset, bool xAxis)
{
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  const LVL1::DataError err(0);
  for (int bit = 0; bit < 8; ++bit) {
    std::string label(splitLine(err.bitName(bit +
                                LVL1::DataError::GLinkParity), xAxis));
    axis->SetBinLabel(bit + 1 + offset, label.c_str());
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

  if (m_configSvc) {
    const TrigConf::CTPConfig* ctpConf = m_configSvc->ctpConfig();
    if (ctpConf) {
      const TrigConf::Menu* ctpMenu = ctpConf->menu();
      if (ctpMenu) {
        const std::vector<TrigConf::TriggerThreshold*>&
  	                            thresholds(ctpMenu->thresholdVector());
        std::vector<TrigConf::TriggerThreshold*>::const_iterator it;
        for (it = thresholds.begin(); it != thresholds.end(); ++it) {
          const std::string thrType((*it)->type());
          if (type == def.emType() || type == def.tauType()) {
            if (thrType != def.emType() && thrType != def.tauType()) continue;
          } else if (thrType != type) continue;
          const int threshNum = (*it)->thresholdNumber();
          if (threshNum >= 0 && threshNum < nthresh) {
            names[threshNum] = (*it)->name();
            found = true;
          }
        }
      }
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

// Put threshold hit values into a string suitable for printing

std::string TrigT1CaloHistogramTool::thresholdString(int val, int nThresh,
                                                              int nBits)
{
  std::ostringstream cval;
  int mask = (1 << nBits) - 1;
  for (int thr = 0; thr < nThresh; ++thr) {
    int hit = (val >> (nBits*thr)) & mask;
    cval << hit;
    if (thr != nThresh-1) cval << " ";
  }
  return cval.str();
}

// Flag which threshold hit values are non-zero and the same

int TrigT1CaloHistogramTool::thresholdsSame(int val1, int val2,
                                                      int nThresh, int nBits)
{
  int result = 0;
  int mask = (1 << nBits) - 1;
  for (int thr = 0; thr < nThresh; ++thr) {
    int hit1 = (val1 >> (nBits*thr)) & mask;
    int hit2 = (val2 >> (nBits*thr)) & mask;
    if (hit1 && (hit1 == hit2)) result |= (1 << thr);
  }
  return result;
}

// Flag which threshold hit values are different

int TrigT1CaloHistogramTool::thresholdsDiff(int val1, int val2,
                                                      int nThresh, int nBits)
{
  int result = 0;
  int mask = (1 << nBits) - 1;
  for (int thr = 0; thr < nThresh; ++thr) {
    int hit1 = (val1 >> (nBits*thr)) & mask;
    int hit2 = (val2 >> (nBits*thr)) & mask;
    if (hit1 != hit2) result |= (1 << thr);
  }
  return result;
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

// Label bins with JEM/CMM crate/module

void TrigT1CaloHistogramTool::jemCMMCrateModule(TH1* hist, int offset,
                                                                    bool xAxis)
{
  jemCrateModule(hist, offset, xAxis);
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  axis->SetBinLabel(1+offset,  "JEM");
  axis->SetBinLabel(33+offset, "CMM");
  axis->SetBinLabel(35+offset, "1/0");
}

// Label bins with JEM frame/local coord

void TrigT1CaloHistogramTool::jemFrameLoc(TH1* hist, int offset, bool xAxis)
{
  const int nFrame = 8;
  const int nLoc   = 4;
  numberPairs(hist, 0, nFrame-1, 0, nLoc-1, 2, offset, xAxis);
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  axis->SetTitle("Frame/Local Coord");
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

// Label bins with PPM error bit names

void TrigT1CaloHistogramTool::ppmErrors(TH1* hist, int offset, bool xAxis)
{
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  const LVL1::DataError err(0);
  for (int bit = 0; bit < 8; ++bit) {
    std::string label(splitLine(err.bitName(bit +
                                LVL1::DataError::ChannelDisabled), xAxis));
    axis->SetBinLabel(bit + 1 + offset, label.c_str());
  }
}

// Label bins with PPM submodule/channel

void TrigT1CaloHistogramTool::ppmSubmoduleChannel(TH1* hist, int offset,
                                                                 bool xAxis)
{
  numberPairs(hist, 0, 15, 0, 3, 4, offset, xAxis);
  TAxis* axis = (xAxis) ? hist->GetXaxis() : hist->GetYaxis();
  axis->SetTitle("Submodule/Channel");
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
                                               const std::string& title)
{
  TH2F*  hist = 0;
  TAxis* axis = 0;
  if (m_shrinkEtaBins) {
    hist = bookPPMEmEtaVsPhi(name, title);
    axis = hist->GetXaxis();
    for (int bin = 1; bin <= 8; ++bin) {
      axis->SetBinLabel(bin,    "+");
      axis->SetBinLabel(bin+58, "+");
    }
  } else {
    hist = book2F(name, title, 50, -2.5, 2.5, 64, 0., 64.);
    hist->SetXTitle("eta");
    const double phiBin     = M_PI/32.;
    const double halfPhiBin = M_PI/64.;
    axis = hist->GetYaxis();
    for (int chan = 0; chan < 64; chan += 4 ) {
      const double rad = chan*phiBin + halfPhiBin;
      axis->SetBinLabel(chan+1, intDoubleString(chan, rad).c_str());
    }
    axis->SetBinLabel(64, "phi");
  }
  return hist;
}

// Book CPM RoI eta vs phi

TH2F* TrigT1CaloHistogramTool::bookCPMRoIEtaVsPhi(const std::string& name,
                                                  const std::string& title)
{
  TH2F*  hist = 0;
  TAxis* axis = 0;
  if (m_shrinkEtaBins) {
    hist = book2F(name, title, 66, -3.3, 3.3, 64, 0., 64.);
    axis = hist->GetXaxis();
    for (int ch = -24; ch < 26; ch+=4) {
      int chan = ch;
      const double eta = (chan/10.);
      axis->SetBinLabel(chan+33, intDoubleString(chan, eta).c_str());
    }
    for (int bin = 1; bin <= 8; ++bin) {
      axis->SetBinLabel(bin,    "+");
      axis->SetBinLabel(bin+58, "+");
    }
  } else {
    hist = book2F(name, title, 50, -2.45, 2.55, 64, 0., 64.);
    hist->SetXTitle("eta");
  }
  const double phiBin = M_PI/32.;
  axis = hist->GetYaxis();
  for (int chan = 0; chan < 64; chan += 4 ) {
    const double rad = (chan + 1)*phiBin;
    axis->SetBinLabel(chan+1, intDoubleString(chan, rad).c_str());
  }
  if (m_shrinkEtaBins) axis->SetBinLabel(64, "etaVphi");
  else                 axis->SetBinLabel(64, "phi");
  return hist;
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
  numbers(hist, 1, 14);
  numbers(hist, 0, 3, 1, 0, false);
  hist->GetXaxis()->SetTitle("Module");
  hist->GetYaxis()->SetTitle("Crate");
  return hist;
}

// Book CPM module vs crate/CMM

TH2F* TrigT1CaloHistogramTool::bookCPMModuleVsCrateCMM(
                const std::string& name, const std::string& title)
{
  TH2F *hist = book2F(name, title, 14, 1., 15., 8, 0., 8.);
  numbers(hist, 1, 14);
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

// Book JEM crate/module vs EX,Ey,Et

TH2F* TrigT1CaloHistogramTool::bookJEMCrateModuleVsExEyEt(
                             const std::string& name, const std::string& title)
{
  TH2F* hist = book2F(name, title, 32, 0., 32., 3, 0., 3.);
  jemCrateModule(hist);
  TAxis* axis = hist->GetYaxis();
  axis->SetBinLabel(1, "Ex");
  axis->SetBinLabel(2, "Ey");
  axis->SetBinLabel(3, "Et");
  return hist;
}

// Book JEM crate/module vs frame/local coord

TH2F* TrigT1CaloHistogramTool::bookJEMCrateModuleVsFrameLoc(
                             const std::string& name, const std::string& title)
{
  TH2F* hist = book2F(name, title, 32, 0., 32., 32, 0., 32.);
  jemCrateModule(hist);
  jemFrameLoc(hist, 0, false);
  return hist;
}

// Book JEM crate/module vs thresholds

TH2F* TrigT1CaloHistogramTool::bookJEMCrateModuleVsThresholds(
                             const std::string& name, const std::string& title)
{
  TH2F* hist = book2F(name, title, 32, 0., 32., 16, 0., 16.);
  jemCrateModule(hist);
  jemThresholds(hist, 0, false);
  return hist;
}

// Book JEM events with error/mismatch vs crate/module

TH2I* TrigT1CaloHistogramTool::bookJEMEventVsCrateModule(
                const std::string& name, const std::string& title)
{
  TH2I* hist = bookEventNumbers(name, title, 32, 0., 32.);
  jemCrateModule(hist, 0, false);
  return hist;
}

// Book JEM module Vs crate

TH2F* TrigT1CaloHistogramTool::bookJEMModuleVsCrate(const std::string& name,
                                                    const std::string& title)
{
  TH2F* hist = book2F(name, title, 16, 0., 16., 2, 0., 2.);
  numbers(hist, 0, 15);
  numbers(hist, 0, 1, 1, 0, false);
  hist->SetXTitle("Module");
  hist->SetYTitle("Crate");
  return hist;
}

// Book JEM eta

TH1F* TrigT1CaloHistogramTool::bookJEMEta(const std::string& name,
                                          const std::string& title)
{
  const int nxbins = 32;
  const double xbins[nxbins+1] = {-4.9,-3.2,-2.9,-2.7,-2.4,-2.2,-2.0,-1.8,-1.6,
                                  -1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,
				  0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,
				  2.7,2.9,3.2,4.9};
  TH1F* hist = book1F(name, title, nxbins, xbins);
  hist->SetXTitle("eta");
  return hist;
}

// Book JEM eta vs phi

TH2F* TrigT1CaloHistogramTool::bookJEMEtaVsPhi(const std::string& name,
                                               const std::string& title)
{
  TH2F*  hist = 0;
  TAxis* axis = 0;
  if (m_shrinkEtaBins) {
    hist = book2F(name, title, 32, -3.2, 3.2, 32, 0., 32.);
    axis = hist->GetXaxis();
    for (int ch = -11; ch < 12; ch+=2) {
      int chan = ch;
      if (chan >= -1) ++chan;
      const double eta = chan/5. + 0.1;
      axis->SetBinLabel(chan+17, intDoubleString(chan, eta).c_str());
    }
    axis->SetBinLabel(2, "-15/-3.05");
    axis->SetBinLabel(4, "-13/-2.55");
    axis->SetBinLabel(29, "12/2.55");
    axis->SetBinLabel(31, "14/3.05");
  } else {
    const int nxbins = 32;
    const double xbins[nxbins+1] = {-4.9,-3.2,-2.9,-2.7,-2.4,-2.2,-2.0,-1.8,
                                    -1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,
				    0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,
				    2.0,2.2,2.4,2.7,2.9,3.2,4.9};
    hist = book2F(name, title, nxbins, xbins, 32, 0., 32.);
    hist->SetXTitle("eta");
  }
  axis = hist->GetYaxis();
  const double phiBin     = M_PI/16.;
  const double halfPhiBin = M_PI/32.;
  for (int chan = 0; chan < 32; chan += 2 ) {
    const double rad = chan*phiBin + halfPhiBin;
    axis->SetBinLabel(chan+1, intDoubleString(chan, rad).c_str());
  }
  if (m_shrinkEtaBins) axis->SetBinLabel(32, "etaVphi");
  else                 axis->SetBinLabel(32, "phi");
  return hist;
}

// Book JEM RoI eta vs phi

TH2F* TrigT1CaloHistogramTool::bookJEMRoIEtaVsPhi(const std::string& name,
                                                  const std::string& title)
{
  TH2F*  hist = 0;
  TAxis* axis = 0;
  if (m_shrinkEtaBins) {
    hist = book2F(name, title, 32, -3.2, 3.2, 32, 0., 32.);
    axis = hist->GetXaxis();
    for (int ch = -10; ch < 12; ch+=2) {
      int chan = ch;
      const double eta = chan/5.;
      axis->SetBinLabel(chan+16, intDoubleString(chan, eta).c_str());
    }
    axis->SetBinLabel(2, "-14/-2.95");
    axis->SetBinLabel(4, "-12/-2.45");
    axis->SetBinLabel(28, "12/2.45");
    axis->SetBinLabel(30, "14/2.95");
    axis->SetBinLabel(32, "16/4.05");
  } else {
    const int nxbins = 32;
    const double xbins[nxbins+1] = {-4.0,-3.05,-2.8,-2.55,-2.3,-2.1,-1.9,-1.7,
                                    -1.5,-1.3,-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,
  				    0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,
				    2.1,2.3,2.55,2.8,3.05,4.0,4.95};
    hist = book2F(name, title, nxbins, xbins, 32, 0., 32.);
    hist->SetXTitle("eta");
  }
  axis = hist->GetYaxis();
  const double phiBin = M_PI/16.;
  for (int chan = 0; chan < 32; chan += 2 ) {
    const double rad = (chan + 1)*phiBin;
    axis->SetBinLabel(chan+1, intDoubleString(chan, rad).c_str());
  }
  if (m_shrinkEtaBins) axis->SetBinLabel(32, "etaVphi");
  else                 axis->SetBinLabel(32, "phi");
  return hist;
}

// Book JEM energy with bins matching QuadLinear encoding

TH1F* TrigT1CaloHistogramTool::bookJEMQuadLinear(const std::string& name,
                                                 const std::string& title,
						 int scale)
{
  if (scale < 1) scale = 1;
  const int eRange = 256;        //range of encoded value
  const int dRange = 4096*scale; //range of decoded value
  LVL1::QuadLinear expand;
  std::set<int> sorted;
  for (int i = 0; i < eRange; ++i) {
    int val = expand.Expand(i);
    if (val != 0) sorted.insert(val);
  }
  double binedges[eRange+2];
  int nbins = (scale > 1) ? 1 : 0;
  std::set<int>::const_iterator iter  = sorted.begin();
  std::set<int>::const_iterator iterE = sorted.end();
  for (; iter != iterE; ++iter) {
    binedges[nbins] = (*iter)*scale;
    ++nbins;
  }
  binedges[0] = 1;
  //if (scale > 1) binedges[1] = scale;
  binedges[nbins] = dRange;
  TH1F* hist = book1F(name, title, nbins, binedges);
  return hist;
}

// Book JEM main jet thresholds

TH1F* TrigT1CaloHistogramTool::bookMainJetThresholds(const std::string& name,
                                                     const std::string& title)
{
  int nbins = TrigConf::L1DataDef::max_J_Threshold_Number();
  TH1F* hist = book1F(name, title, nbins, 0, nbins);
  mainJetThresholds(hist);
  return hist;
}

// Book JEM backward jet thresholds

TH1F* TrigT1CaloHistogramTool::bookBackwardJetThresholds(
                            const std::string& name, const std::string& title)
{
  int nbins = TrigConf::L1DataDef::max_JB_Threshold_Number();
  TH1F* hist = book1F(name, title, nbins, 0, nbins);
  backwardJetThresholds(hist);
  return hist;
}

// Book JEM forward jet thresholds

TH1F* TrigT1CaloHistogramTool::bookForwardJetThresholds(
                            const std::string& name, const std::string& title)
{
  int nbins = TrigConf::L1DataDef::max_JF_Threshold_Number();
  TH1F* hist = book1F(name, title, nbins, 0, nbins);
  forwardJetThresholds(hist);
  return hist;
}

// Book JEM JetEt thresholds

TH1F* TrigT1CaloHistogramTool::bookJetEtThresholds(const std::string& name,
                                                   const std::string& title)
{
  int nbins = TrigConf::L1DataDef::max_JE_Threshold_Number();
  TH1F* hist = book1F(name, title, nbins, 0, nbins);
  jetEtThresholds(hist);
  return hist;
}

// Book JEM MissingEt thresholds

TH1F* TrigT1CaloHistogramTool::bookMissingEtThresholds(const std::string& name,
                                                       const std::string& title)
{
  int nbins = TrigConf::L1DataDef::max_XE_Threshold_Number();
  TH1F* hist = book1F(name, title, nbins, 0, nbins);
  missingEtThresholds(hist);
  return hist;
}

// Book JEM SumEt thresholds

TH1F* TrigT1CaloHistogramTool::bookSumEtThresholds(const std::string& name,
                                                   const std::string& title)
{
  int nbins = TrigConf::L1DataDef::max_TE_Threshold_Number();
  TH1F* hist = book1F(name, title, nbins, 0, nbins);
  sumEtThresholds(hist);
  return hist;
}

//===========================================================================
//  Booking Utilities - PPM
//===========================================================================

// Book PPM Em eta

TH1F* TrigT1CaloHistogramTool::bookPPMEmEta(const std::string name,
                                            const std::string title)
{
  const int nxbins = 66;
  const double xbins[nxbins+1] = {-4.9,-4.475,-4.050,-3.625,-3.2,-3.1,-2.9,
                                  -2.7,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,
      		                  -1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,
				  -1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
				  -0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,
				  0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,
				  1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.7,2.9,
				  3.1,3.2,3.625,4.050,4.475,4.9};
  TH1F* hist = book1F(name, title, nxbins, xbins);
  hist->SetXTitle("eta");
  return hist;
}

// Book PPM Had eta

TH1F* TrigT1CaloHistogramTool::bookPPMHadEta(const std::string name,
                                             const std::string title)
{
  const int nxbins = 62;
  const double xbins[nxbins+1] = {-4.9,-4.050,-3.2,-3.1,-2.9,
                                  -2.7,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,
      		                  -1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,
				  -1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
				  -0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,
				  0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,
				  1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.7,2.9,
				  3.1,3.2,4.050,4.9};
  TH1F* hist = book1F(name, title, nxbins, xbins);
  hist->SetXTitle("eta");
  return hist;
}
 
// Book PPM Em eta vs phi

TH2F* TrigT1CaloHistogramTool::bookPPMEmEtaVsPhi(const std::string name,
                                                 const std::string title)
{
  TH2F*  hist = 0;
  TAxis* axis = 0;
  if (m_shrinkEtaBins) {
    hist = book2F(name, title, 66, -3.3, 3.3, 64, 0., 64.);
    axis = hist->GetXaxis();
    for (int ch = -25; ch < 25; ch+=4) {
      int chan = ch;
      if (chan >= -1) ++chan;
      const double eta = (chan/10.)+0.05;
      axis->SetBinLabel(chan+34, intDoubleString(chan, eta).c_str());
    }
    axis->SetBinLabel(1, "-49/-4.69");
    axis->SetBinLabel(5, "-32/-3.15");
    axis->SetBinLabel(62, "31/3.15");
    axis->SetBinLabel(66, "44/4.69");
  } else {
    const int nxbins = 66;
    const double xbins[nxbins+1] = {-4.9,-4.475,-4.050,-3.625,-3.2,-3.1,-2.9,
                                    -2.7,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,
  				    -1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,
				    -1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
				    -0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,
				    0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,
				    1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.7,2.9,
				    3.1,3.2,3.625,4.050,4.475,4.9};
    hist = book2F(name, title, nxbins, xbins, 64, 0., 64.);
    hist->SetXTitle("eta");
  }

  axis = hist->GetYaxis();
  const double phiBin     = M_PI/32.;
  const double halfPhiBin = M_PI/64.;
  for (int chan = 0; chan < 64; chan += 4 ) {
    const double rad = chan*phiBin + halfPhiBin;
    axis->SetBinLabel(chan+1, intDoubleString(chan, rad).c_str());
  }
  if (m_shrinkEtaBins) axis->SetBinLabel(64, "etaVphi");
  else                 axis->SetBinLabel(64, "phi");
  return hist;

}

// Book PPM Had eta vs phi

TH2F* TrigT1CaloHistogramTool::bookPPMHadEtaVsPhi(const std::string name,
                                                  const std::string title)
{
  TH2F* hist = bookPPMEmEtaVsPhi(name, title);
  if (m_shrinkEtaBins) {
    TAxis* axis = hist->GetXaxis();
    axis->SetBinLabel(1, "-49/-4.48");
    axis->SetBinLabel(66, "44/4.48");
  }
  return hist;

}

// Book PPM Em eta vs phi profile

TProfile2D* TrigT1CaloHistogramTool::bookProfilePPMEmEtaVsPhi(
                              const std::string name, const std::string title)
{
  // todo - remove duplication with above
  TProfile2D* hist = 0;
  TAxis*      axis = 0;
  if (m_shrinkEtaBins) {
    hist = bookProfile2D(name, title, 66, -3.3, 3.3, 64, 0., 64.);
    axis = hist->GetXaxis();
    for (int ch = -25; ch < 25; ch+=4) {
      int chan = ch;
      if (chan >= -1) ++chan;
      const double eta = (chan/10.)+0.05;
      axis->SetBinLabel(chan+34, intDoubleString(chan, eta).c_str());
    }
    axis->SetBinLabel(1, "-49/-4.69");
    axis->SetBinLabel(5, "-32/-3.15");
    axis->SetBinLabel(62, "31/3.15");
    axis->SetBinLabel(66, "44/4.69");
  } else {
    const int nxbins = 66;
    const double xbins[nxbins+1] = {-4.9,-4.475,-4.050,-3.625,-3.2,-3.1,-2.9,
                                    -2.7,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,
  				    -1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,
				    -1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
				    -0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,
				    0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,
				    1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.7,2.9,
				    3.1,3.2,3.625,4.050,4.475,4.9};
    hist = bookProfile2D(name, title, nxbins, xbins, 64, 0., 64.);
    hist->SetXTitle("eta");
  }

  axis = hist->GetYaxis();
  const double phiBin     = M_PI/32.;
  const double halfPhiBin = M_PI/64.;
  for (int chan = 0; chan < 64; chan += 4 ) {
    const double rad = chan*phiBin + halfPhiBin;
    axis->SetBinLabel(chan+1, intDoubleString(chan, rad).c_str());
  }
  if (m_shrinkEtaBins) axis->SetBinLabel(64, "etaVphi");
  else                 axis->SetBinLabel(64, "phi");
  return hist;

}

// Book PPM Had eta vs phi profile

TProfile2D* TrigT1CaloHistogramTool::bookProfilePPMHadEtaVsPhi(
                              const std::string name, const std::string title)
{
  TProfile2D* hist = bookProfilePPMEmEtaVsPhi(name, title);
  if (m_shrinkEtaBins) {
    TAxis* axis = hist->GetXaxis();
    axis->SetBinLabel(1, "-49/-4.48");
    axis->SetBinLabel(66, "44/4.48");
  }
  return hist;

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

// Book PPM Crate/Module vs Submodule/Channel

TH2F* TrigT1CaloHistogramTool::bookPPMCrateModuleVsSubmoduleChannel(
                const std::string& name, const std::string& title,
                int firstCrate, int lastCrate)
{
  int nbins = (lastCrate-firstCrate+1)*16;
  TH2F* hist = book2F(name, title, nbins, 0., nbins, 64, 0., 64.);
  ppmCrateModule(hist, firstCrate, lastCrate);
  ppmSubmoduleChannel(hist, 0, false);
  return hist;
}

// Book PPM Crate/Module vs Submodule/Channel profile

TProfile2D*
  TrigT1CaloHistogramTool::bookProfilePPMCrateModuleVsSubmoduleChannel(
                const std::string& name, const std::string& title,
                int firstCrate, int lastCrate)
{
  int nbins = (lastCrate-firstCrate+1)*16;
  TProfile2D* hist = bookProfile2D(name, title, nbins, 0., nbins, 64, 0., 64.);
  ppmCrateModule(hist, firstCrate, lastCrate);
  ppmSubmoduleChannel(hist, 0, false);
  return hist;
}

// Book PPM SubStatus vs crate/module

TH2F* TrigT1CaloHistogramTool::bookPPMSubStatusVsCrateModule(
                const std::string& name, const std::string& title,
		int firstCrate, int lastCrate)
{
  int nbins = (lastCrate-firstCrate+1)*16;
  TH2F* hist = book2F(name, title, 8, 0., 8., nbins, 0., nbins);
  subStatus(hist);
  ppmCrateModule(hist, firstCrate, lastCrate, 0, false);
  return hist;
}

// Book PPM ASIC errors vs crate/module

TH2F* TrigT1CaloHistogramTool::bookPPMErrorsVsCrateModule(
                const std::string& name, const std::string& title,
		int firstCrate, int lastCrate)
{
  int nbins = (lastCrate-firstCrate+1)*16;
  TH2F* hist = book2F(name, title, 8, 0., 8., nbins, 0., nbins);
  ppmErrors(hist);
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
  int lastVal = 0;
  for (int binx = 1; binx <= nbins; ++binx) {
    const int val = hist->GetBinContent(binx, biny);
    if (val == 0) {
      int eventNumber = 0;
      const EventInfo* evInfo = 0;
      StatusCode sc = evtStore()->retrieve(evInfo);
      if (sc.isFailure()) {
        msg(MSG::DEBUG) << "No EventInfo found" << endreq;
      } else {
        const EventID* evID = evInfo->event_ID();
        if (evID) eventNumber = evID->event_number();
      }
      if (eventNumber != lastVal) hist->SetBinContent(binx, biny, eventNumber);
      break;
    } else lastVal = val;
  }
}

// Fill weighted thresholds 1D

void TrigT1CaloHistogramTool::fillThresholds(TH1* hist, int val, int nThresh,
                                                      int nBits, int offset)
{
  if (val) {
    int mask = (1 << nBits) - 1;
    for (int thr = 0; thr < nThresh; ++thr) {
      int hit = (val >> (nBits*thr)) & mask;
      if (hit) hist->Fill(thr + offset, hit);
    }
  }
}

// Fill weighted thresholds 2D, X axis

void TrigT1CaloHistogramTool::fillThresholdsVsY(TH2* hist, int val, int y,
                                          int nThresh, int nBits, int offset)
{
  if (val) {
    int mask = (1 << nBits) - 1;
    for (int thr = 0; thr < nThresh; ++thr) {
      int hit = (val >> (nBits*thr)) & mask;
      if (hit) hist->Fill(thr + offset, y, hit);
    }
  }
} 

// Fill weighted thresholds 2D, Y axis

void TrigT1CaloHistogramTool::fillXVsThresholds(TH2* hist, int x, int val,
                                          int nThresh, int nBits, int offset)
{
  if (val) {
    int mask = (1 << nBits) - 1;
    for (int thr = 0; thr < nThresh; ++thr) {
      int hit = (val >> (nBits*thr)) & mask;
      if (hit) hist->Fill(x, thr + offset, hit);
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
  double etaMod = eta;
  if (m_shrinkEtaBins) etaMod -= 0.05;
  hist->Fill(etaMod, phiMod, weight);
}

//===========================================================================
//  Filling Utilities - JEM
//===========================================================================

// Fill JEM eta vs phi

void TrigT1CaloHistogramTool::fillJEMEtaVsPhi(TH2* hist, double eta,
                                              double phi, double weight)
{
  const double phiMod = phi * m_phiScaleJE;
  double etaMod = eta;
  double absEta = fabs(eta);
  if (m_shrinkEtaBins && absEta > 2.4) {
    int offset = 1;
    if      (absEta > 3.2) offset = 4;
    else if (absEta > 2.9) offset = 3;
    else if (absEta > 2.7) offset = 2;
    etaMod = 2.3 + 0.2*offset;
    if (eta < 0.) etaMod = -etaMod;
  }
  if (eta < -3.2 || eta > 3.2) {
    // Fill two bins for FCAL
    hist->Fill(etaMod, phiMod + 0.5, weight);
    hist->Fill(etaMod, phiMod - 0.5, weight);
  } else hist->Fill(etaMod, phiMod, weight);
}

// Fill JEM RoI eta vs phi

void TrigT1CaloHistogramTool::fillJEMRoIEtaVsPhi(TH2* hist, double eta,
                                                 double phi, double weight)
{
  const double phiMod = phi * m_phiScaleJE - 0.5;
  double etaMod = eta;
  double absEta = fabs(eta);
  if (m_shrinkEtaBins && absEta > 2.3) {
    int offset = 1;
    if      (absEta > 4.0)  offset = 5;
    else if (absEta > 3.05) offset = 4;
    else if (absEta > 2.8)  offset = 3;
    else if (absEta > 2.55) offset = 2;
    etaMod = 2.2 + 0.2*offset;
    if (eta < 0.) etaMod = -etaMod;
  }
  const double etaShift = (m_shrinkEtaBins) ? 0.1 : 0.;
  // JEPRoIDecoder returns eta=3.9 for both of the two forwardmost bins
  if (eta > 3.8 && eta < 4.0) {
    const double eta2 = (m_shrinkEtaBins) ? 3.2 : 4.05;
    hist->Fill(eta2 - etaShift, phiMod, weight);
  }
  hist->Fill(etaMod - etaShift, phiMod, weight);
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
  double etaMod = eta;
  double absEta = fabs(eta);
  if (absEta > 3.2) {
    if (m_shrinkEtaBins) {
      etaMod = 2.9 + 0.1*(absEta-3.2)/0.425;
      if (eta < 0.) etaMod = -etaMod;
    }
    // Fill four bins in phi
    hist->Fill(etaMod, phiMod + 1.5, weight);
    hist->Fill(etaMod, phiMod + 0.5, weight);
    hist->Fill(etaMod, phiMod - 0.5, weight);
    hist->Fill(etaMod, phiMod - 1.5, weight);
  } else if (absEta > 2.5) {
    if (m_shrinkEtaBins) {
      etaMod = (absEta > 3.1) ? 2.85 : 2.5 + (absEta-2.5)/2.;
      if (eta < 0.) etaMod = -etaMod;
    }
    // Fill two bins in phi
    hist->Fill(etaMod, phiMod + 0.5, weight);
    hist->Fill(etaMod, phiMod - 0.5, weight);
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
  double etaMod = eta;
  if (m_shrinkEtaBins) {
    double absEta = fabs(eta);
    if (absEta > 3.2) {
      etaMod = 2.9 + 0.1*(absEta-3.2)/0.425;
      if (eta < 0.) etaMod = -etaMod;
    } else if (absEta > 2.5) {
      etaMod = (absEta > 3.1) ? 2.85 : 2.5 + (absEta-2.5)/2.;
      if (eta < 0.) etaMod = -etaMod;
    }
  }
  return hist->FindBin(etaMod, phiMod);
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

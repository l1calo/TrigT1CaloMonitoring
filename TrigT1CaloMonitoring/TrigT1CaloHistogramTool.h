#ifndef TRIGT1CALOBYTESTREAM_TRIGT1CALOHISTOGRAMTOOL_H
#define TRIGT1CALOBYTESTREAM_TRIGT1CALOHISTOGRAMTOOL_H

#include <string>
#include <vector>

#include "GaudiKernel/ServiceHandle.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"

#include "TH1.h"
#include "TH2.h"

#include "TrigConfigSvc/ILVL1ConfigSvc.h"

class TH1;
class TH2;
class TH1F;
class TH2F;
class TProfile;
class TProfile2D;

class IInterface;
class InterfaceID;
class StatusCode;

/** Tool to provide histogramming utilities
 *
 *  @author Peter Faulkner
 */

class TrigT1CaloHistogramTool : public AthAlgTool {

 public:
   TrigT1CaloHistogramTool(const std::string& type, const std::string& name,
                           const IInterface* parent);
   virtual ~TrigT1CaloHistogramTool();

   /// AlgTool InterfaceID
   static const InterfaceID& interfaceID();

   virtual StatusCode initialize();
   virtual StatusCode finalize();

   /// Set current MonGroup
   void setMonGroup(ManagedMonitorToolBase::MonGroup* group)
     { m_monGroup = group; }
   /// Unset current MonGroup
   void unsetMonGroup() { m_monGroup = 0; }

   ////////////////////////////////
   // Labelling Utilities - General
   ////////////////////////////////

   /// Label bins with number pairs
   void numberPairs(TH1* hist, int firstMin, int firstMax,
                               int secondMin, int secondMax,
			       int step = 1, int offset = 0, bool xAxis = true);
   /// Label bins with numbers
   void numbers(TH1* hist, int min, int max, int step = 1, int offset = 0,
                                                           bool xAxis = true);
   /// Label bins with sub-status error bit names
   void subStatus(TH1* hist, int offset = 0, bool xAxis = true);
   /// Get list of threshold names for given type
   bool thresholdNames(const std::string& type,
                             std::vector<std::string>& names);
   /// Label bins with threshold names
   void thresholds(TH1* hist, const std::string& type, int offset = 0,
                                                       bool xAxis = true);

   ////////////////////////////////
   // Labelling Utilities - CPM
   ////////////////////////////////

   /// Label bins with CPM chip/local coordinate
   void cpmChipLocalCoord(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with CPM/CMM crate/module
   void cpmCMMCrateModule(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with CPM crate/CMM
   void cpmCrateCMM(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with CPM crate/module
   void cpmCrateModule(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with CPM RoI threshold names
   void cpmThresholds(TH1* hist, int offset = 0, bool xAxis = true);

   ////////////////////////////////
   // Labelling Utilities - JEM
   ////////////////////////////////

   /// Label bins with JEM crate/module
   void jemCrateModule(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with JEM RoI threshold names
   void jemThresholds(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with Main Jet threshold names
   void mainJetThresholds(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with Backward Jet threshold names
   void backwardJetThresholds(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with Forward Jet threshold names
   void forwardJetThresholds(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with JetEt threshold names
   void jetEtThresholds(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with MissingEt threshold names
   void missingEtThresholds(TH1* hist, int offset = 0, bool xAxis = true);
   /// Label bins with SumEt threshold names
   void sumEtThresholds(TH1* hist, int offset = 0, bool xAxis = true);

   ////////////////////////////////
   // Labelling Utilities - PPM
   ////////////////////////////////

   /// Label bins with PPM crate/module

   void ppmCrateModule(TH1* hist, int firstCrate, int lastCrate,
                                       int offset = 0, bool xAxis = true);

   ////////////////////////////////
   // Booking Utilities - General
   ////////////////////////////////

   /// Book and register a 1D histogram
   TH1F* book1F(const std::string& name, const std::string& title,
                                         int nx, double xmin, double xmax);
   /// Book and register a 1D histogram with variable width bins
   TH1F* book1F(const std::string& name, const std::string& title,
                                         int nx, const double* xbins);
   /// Book and register a 1D profile histogram
   TProfile* bookProfile(const std::string& name, const std::string& title,
                                         int nx, double xmin, double xmax);
   /// Book and register a 1D profile histogram with variable width bins
   TProfile* bookProfile(const std::string& name, const std::string& title,
                                         int nx, const double* xbins);
   /// Book and register a 2D histogram
   TH2F* book2F(const std::string& name, const std::string& title,
                                         int nx, double xmin, double xmax,
		                         int ny, double ymin, double ymax);
   /// Book and register a 2D histogram with variable width bins
   TH2F* book2F(const std::string& name, const std::string& title,
                                         int nx, const double* xbins,
					 int ny, double ymin, double ymax);
   /// Book and register a 2D histogram of integers displayed as text
   TH2I* book2I(const std::string& name, const std::string& title,
                                         int nx, double xmin, double xmax,
					 int ny, double ymin, double ymax);
   /// Book and register a 2D profile histogram
   TProfile2D* bookProfile2D(const std::string& name, const std::string& title,
                                         int nx, double xmin, double xmax,
					 int ny, double ymin, double ymax);
   /// Book and register a 2D profile histogram with variable width bins
   TProfile2D* bookProfile2D(const std::string& name, const std::string& title,
                                         int nx, const double* xbins,
					 int ny, double ymin, double ymax);
   /// Book and register a 2D histogram containing event numbers as bin contents
   TH2I* bookEventNumbers(const std::string& name, const std::string& title,
                                         int ny, double ymin, double ymax);
   /// Register a histogram
   void registerHist(TH1* hist);

   ////////////////////////////////
   // Booking Utilities - CPM
   ////////////////////////////////

   /// Book CPM crate/module vs chip/local coordinate
   TH2F* bookCPMCrateModuleVsChipLocalCoord(const std::string& name,
                                            const std::string& title);
   /// Book CPM crate/module vs FPGA
   TH2F* bookCPMCrateModuleVsFPGA(const std::string& name,
                                  const std::string& title);
   /// Book CPM crate/module vs thresholds
   TH2F* bookCPMCrateModuleVsThreshold(const std::string& name,
                                       const std::string& title);
   /// Book CPM eta vs phi
   TH2F* bookCPMEtaVsPhi(const std::string& name, const std::string& title,
                                                       bool isRoi = false);
   /// Book CPM RoI eta vs phi
   TH2F* bookCPMRoIEtaVsPhi(const std::string& name, const std::string& title);
   /// Book CPM events with error/mismatch vs crate/module
   TH2I* bookCPMEventVsCrateModule(const std::string& name,
                                   const std::string& title);
   /// Book CPM module vs crate
   TH2F* bookCPMModuleVsCrate(const std::string& name,
                              const std::string& title);
   /// Book CPM module vs crate/CMM
   TH2F* bookCPMModuleVsCrateCMM(const std::string& name,
                                 const std::string& title);
   /// Book CPM sub-status errors vs crate/module
   TH2F* bookCPMSubStatusVsCrateModule(const std::string& name,
                                       const std::string& title);
   /// Book CPM Sum/CMM
   TH1F* bookCPMSumCMM(const std::string& name,const std::string& title);
   /// Book CPM Sum vs Threshold
   TH2F* bookCPMSumVsThreshold(const std::string& name,
                               const std::string& title);

   ////////////////////////////////
   // Booking Utilities - JEM
   ////////////////////////////////

   /// Book JEM eta vs phi
   TH2F* bookJEMEtaVsPhi(const std::string& name, const std::string& title,
                                                       bool isRoi = false);
   /// Book JEM RoI eta vs phi
   TH2F* bookJEMRoIEtaVsPhi(const std::string& name, const std::string& title);

   ////////////////////////////////
   // Booking Utilities - PPM
   ////////////////////////////////

   /// Book PPM Em eta vs phi
   TH2F* bookPPMEmEtaVsPhi(const std::string name, const std::string title);
   /// Book PPM Had eta vs phi
   TH2F* bookPPMHadEtaVsPhi(const std::string name, const std::string title);
   /// Book PPM Em eta vs phi profile
   TProfile2D* bookProfilePPMEmEtaVsPhi(const std::string name,
                                        const std::string title);
   /// Book PPM Had eta vs phi profile
   TProfile2D* bookProfilePPMHadEtaVsPhi(const std::string name,
                                         const std::string title);
   /// Book PPM events with error/mismatch vs crate/module
   TH2I* bookPPMEventVsCrateModule(const std::string& name,
                                   const std::string& title,
				   int firstCrate, int lastCrate);

   ////////////////////////////////
   // Filling Utilities - General
   ////////////////////////////////

   /// Fill Error/Mismatch event number
   void fillEventNumber(TH2I* hist, double y);

   ////////////////////////////////
   // Filling Utilities - CPM
   ////////////////////////////////

   /// Fill CPM eta vs phi
   void fillCPMEtaVsPhi(TH2* hist, double eta, double phi, double weight = 1.);
   /// Fill CPM RoI eta vs phi
   void fillCPMRoIEtaVsPhi(TH2* hist, double eta, double phi,
                                                            double weight = 1.);

   ////////////////////////////////
   // Filling Utilities - JEM
   ////////////////////////////////

   ///Fill JEM eta vs phi
   void fillJEMEtaVsPhi(TH2* hist, double eta, double phi, double weight = 1.);
   /// Fill JEM RoI eta vs phi
   void fillJEMRoIEtaVsPhi(TH2* hist, double eta, double phi,
                                                            double weight = 1.);
   /// Fill JEM phi allowing for granularity varying with eta
   void fillJEMPhi(TH1* hist, double eta, double phi, double weight = 1.);

   ////////////////////////////////
   // Filling Utilities - PPM
   ////////////////////////////////

   /// Fill PPM Em eta vs phi
   void fillPPMEmEtaVsPhi(TH2* hist, double eta, double phi,
                                                            double weight = 1.);
   /// Fill PPM Had eta vs phi
   void fillPPMHadEtaVsPhi(TH2* hist, double eta, double phi,
                                                            double weight = 1.);
   /// Fill PPM phi allowing for granularity varying with eta
   void fillPPMPhi(TH1* hist, double eta, double phi, double weight = 1.);
   /// Find bin in Em eta vs phi
   int findBinPPMEmEtaVsPhi(TH2* hist, double eta, double phi);
   /// Find bin in Had eta vs phi
   int findBinPPMHadEtaVsPhi(TH2* hist, double eta, double phi);
   /// Set bin content and optionally error in Em eta vs phi
   void setBinPPMEmEtaVsPhi(TH2* hist, int bin, double content,
                                                double error = -1.);
   /// Set bin content and optionally error in Had eta vs phi
   void setBinPPMHadEtaVsPhi(TH2* hist, int bin, double content,
                                                 double error = -1.);

 private:

   /// Trigger configuration service
   ServiceHandle<TrigConf::ILVL1ConfigSvc> m_configSvc;
   /// Current MonGroup or 0 if not wanted
   ManagedMonitorToolBase::MonGroup* m_monGroup;
   /// Phi scale for trigger tower eta/phi plots
   double m_phiScaleTT;
   /// Phi scale for jet element eta/phi plots
   double m_phiScaleJE;
   /// Number of Error Event Number Samples
   int m_eventSamples;

};

#endif

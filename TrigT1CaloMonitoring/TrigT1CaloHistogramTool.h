#ifndef TRIGT1CALOBYTESTREAM_TRIGT1CALOHISTOGRAMTOOL_H
#define TRIGT1CALOBYTESTREAM_TRIGT1CALOHISTOGRAMTOOL_H

#include <string>
#include <vector>

#include "GaudiKernel/ServiceHandle.h"
#include "AthenaBaseComps/AthAlgTool.h"

#include "TH1.h"
#include "TH2.h"

#include "TrigConfigSvc/ILVL1ConfigSvc.h"

class TH1;
class TH1F;
class TH2F;

class IInterface;
class InterfaceID;
class StatusCode;
class ManagedMonitorToolBase;
class MonGroup;

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

   /// Get list of threshold names for given type
   bool thresholdNames(const std::string& type, std::vector<std::string>& names);
   /// Label bins with threshold names
   void thresholds(TH1* hist, const std::string& type, int offset = 0,
                                                       bool xAxis = true);
   /// Label bins with numbers
   void numbers(TH1* hist, int min, int max, int step = 1, int offset = 0,
                                                           bool xAxis = true);
   /// Label bins with number pairs
   void numberPairs(TH1* hist, int firstMin, int firstMax,
                               int secondMin, int secondMax,
			       int step = 1, int offset = 0, bool xAxis = true);
   /// Label bins with sub-status error bit names
   void subStatus(TH1* hist, bool xAxis = true);

   ////////////////////////////////
   // Labelling Utilities - CPM
   ////////////////////////////////

   /// Label bins with CPM crate/module
   void cpmCrateModule(TH1* hist, bool xAxis = true);
   /// Label bins with CPM RoI threshold names
   void cpmThresholds(TH1* hist, int offset = 0, bool xAxis = true);

   ////////////////////////////////
   // Labelling Utilities - JEM
   ////////////////////////////////

   /// Label bins with JEM crate/module
   void jemCrateModule(TH1* hist, bool xAxis = true);
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
   // Booking Utilities - General
   ////////////////////////////////

   /// Book and register a 1D histogram
   TH1F* book1F(const std::string& name, const std::string& title,
                                         int nx, double xmin, double xmax);
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

   ////////////////////////////////
   // Booking Utilities - CPM
   ////////////////////////////////

   /// Book CPM crate/module vs thresholds
   TH2F* bookCPMCrateModuleVsThreshold(const std::string& name,
                                       const std::string& title);
   /// Book CPM eta vs phi
   TH2F* bookCPMEtaVsPhi(const std::string& name, const std::string& title,
                                                       bool isRoi = false);
   /// Book CPM RoI eta vs phi
   TH2F* bookCPMRoIEtaVsPhi(const std::string& name, const std::string& title);
   /// Book CPM sub-status errors vs crate/module
   TH2F* bookCPMSubStatusVsCrateModule(const std::string& name,
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

   ////////////////////////////////
   // Filling Utilities - CPM
   ////////////////////////////////

   /// Fill CPM eta vs phi
   void fillCPMEtaVsPhi(TH2F* hist, double eta, double phi, double weight = 1.);
   /// Fill CPM RoI eta vs phi
   void fillCPMRoIEtaVsPhi(TH2F* hist, double eta, double phi,
                                                            double weight = 1.);

   ////////////////////////////////
   // Filling Utilities - JEM
   ////////////////////////////////

   ///Fill JEM eta vs phi
   void fillJEMEtaVsPhi(TH2F* hist, double eta, double phi, double weight = 1.);
   /// Fill JEM RoI eta vs phi
   void fillJEMRoIEtaVsPhi(TH2F* hist, double eta, double phi,
                                                            double weight = 1.);

   ////////////////////////////////
   // Filling Utilities - PPM
   ////////////////////////////////

   /// Fill PPM Em eta vs phi
   void fillPPMEmEtaVsPhi(TH2F* hist, double eta, double phi,
                                                            double weight = 1.);
   void fillPPMHadEtaVsPhi(TH2F* hist, double eta, double phi,
                                                            double weight = 1.);

 private:

   /// Trigger configuration service
   ServiceHandle<TrigConf::ILVL1ConfigSvc> m_configSvc;
   /// Current MonGroup or 0 if not wanted
   ManagedMonitorToolBase::MonGroup* m_monGroup;
   /// Phi scale for trigger tower eta/phi plots
   double m_phiScaleTT;
   /// Phi scale for jet element eta/phi plots
   double m_phiScaleJE;

};

#endif

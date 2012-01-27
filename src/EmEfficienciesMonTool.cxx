// ********************************************************************
//
// NAME:     EmEfficienciesMonTool.cxx
// PACKAGE:  TrigT1CaloMonitoring  
//
// AUTHOR:   Hardeep Bansil
//           Adapted for monitoring: Peter Faulkner
//           
//
// ********************************************************************

#include <cmath>
#include <cstdio>

#include "TH1.h"
#include "LWHists/LWHist.h"
#include "LWHists/TH1F_LW.h"
#include "LWHists/TH2F_LW.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/StatusCode.h"
#include "StoreGate/StoreGateSvc.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CoralBase/AttributeList.h"
#include "AthenaPoolUtilities/CondAttrListCollection.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigConfL1Data/L1DataDef.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigT1CaloEvent/TriggerTowerCollection.h"
#include "AnalysisTriggerEvent/LVL1_ROI.h"
#include "AnalysisTriggerEvent/EmTau_ROI.h"
#include "AnalysisTriggerEvent/Jet_ROI.h"
#include "TrigT1CaloCalibToolInterfaces/IL1CaloOfflineTriggerTowerTools.h"
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"
#include "egammaEvent/ElectronContainer.h"
#include "egammaEvent/PhotonContainer.h"
#include "egammaEvent/Electron.h"
#include "egammaEvent/Photon.h"
#include "egammaEvent/egammaPIDdefs.h"
#include "VxVertex/VxContainer.h"
#include "VxVertex/VxTrackAtVertex.h"
#include "CaloEvent/CaloCell.h"
#include "CaloEvent/CaloCluster.h"

#include "TrigT1CaloMonitoring/EmEfficienciesMonTool.h"
#include "TrigT1CaloMonitoringTools/TrigT1CaloLWHistogramTool.h"

/*---------------------------------------------------------*/
EmEfficienciesMonTool::EmEfficienciesMonTool(const std::string & type,
		const std::string & name, const IInterface* parent) 
		  : ManagedMonitorToolBase(type, name, parent),
			m_histTool("TrigT1CaloLWHistogramTool"),
			m_tools("LVL1::L1CaloOfflineTriggerTowerTools/L1CaloOfflineTriggerTowerTools"),
			m_trigger("Trig::TrigDecisionTool/TrigDecisionTool"),
			m_dbPpmDeadChannelsFolder("/TRIGGER/L1Calo/V1/Calibration/PpmDeadChannels"),
			m_triggerTowersLocation("TriggerTowers"),
			m_lvl1RoIsLocation("LVL1_ROI"),
			m_offlineElectronsLocation("ElectronAODCollection"),
			m_offlinePhotonsLocation("PhotonAODCollection"),
			m_primaryVertexLocation("VxPrimaryCandidate"), 
			m_eventInfo(0),
			m_dbPpmDeadChannels(0),
			m_triggerTowers(0), 
			m_lvl1RoIs(0),
			m_offlineElectrons(0), 
			m_offlinePhotons(0), 
			m_primaryVtx(0),
			m_numEvents(0), 
			m_numOffElec(0), 
			m_numOffPhot(0),
			m_numOffElecInContainer(0), 
			m_numOffPhotInContainer(0),
			m_numOffElecPassCuts(0), 
			m_numOffPhotPassCuts(0),
			m_numOffElecTriggered(0),
			m_numOffPhotTriggered(0),
			m_numEmObjPassTrigger(0), 
			m_numEmObjTotal(0),
			m_passed_L1_Jet_Trigger(false),
			m_passed_EF_SingleJet_Trigger(false),
			m_passed_EF_MultiJet_Trigger(false),
			m_passed_EF_egTau_Trigger(false), 
			m_passed_EF_Trigger(false),
			m_emBitMask(0),
			m_h_ClusterRaw_Et_gdEta(0), 
			m_h_ClusterRaw_Et_triggered_gdEta(0),
			m_h_ClusterRaw_Et_triggered_Eff(0),
			m_h_ClusterRaw_10GeV_Eta_vs_Phi(0),
			m_h_ClusterRaw_20GeV_Eta_vs_Phi(0),
			m_h_ClusterRaw_30GeV_Eta_vs_Phi(0), 
			m_h_TrigTower_emBadCalo(0),
			m_h_TrigTower_emDeadChannel(0)

/*---------------------------------------------------------*/
{

	declareProperty("HistogramTool", m_histTool);
	declareProperty("OfflineTriggerTowerTool", m_tools);
	declareProperty("TrigDecisionTool", m_trigger);
	declareProperty("DeadChannelsFolder", m_dbPpmDeadChannelsFolder);
	declareProperty("TriggerTowersLocation", m_triggerTowersLocation);
	declareProperty("RoIsLocation", m_lvl1RoIsLocation);
	declareProperty("OfflineElectronsLocation", m_offlineElectronsLocation);
	declareProperty("OfflinePhotonsLocation", m_offlinePhotonsLocation);
	declareProperty("PrimaryVertexLocation", m_primaryVertexLocation);

	declareProperty("RootDirectory", m_rootDir = "L1Calo");
	declareProperty("UseTrigger", m_useTrigger = true);
	declareProperty("TriggerStrings", m_triggerStrings);
	declareProperty("TestMerging", m_testMerge = false);

	// HSB - python cuts
	declareProperty("useDeltaRMatch", m_useDeltaRMatch = true);
	declareProperty("goodEMDeltaRMatch_Cut", m_goodEMDeltaRMatch_Cut = 0.15);
	declareProperty("goodHadDeltaRMatch_Cut", m_goodHadDeltaRMatch_Cut = 0.2);
	declareProperty("useDeltaEtaPhiMatch", m_useDeltaEtaPhiMatch = false);
	declareProperty("goodEMDeltaPhiMatch_Cut", m_goodEMDeltaPhiMatch_Cut = 0.2);
	declareProperty("goodEMDeltaEtaMatch_Cut", m_goodEMDeltaEtaMatch_Cut = 0.2);
	declareProperty("goodHadDeltaEtaMatch_Cut", m_goodHadDeltaEtaMatch_Cut = 0.3);
	declareProperty("goodHadDeltaPhiMatch_Cut", m_goodHadDeltaPhiMatch_Cut = 0.3);
	declareProperty("deltaRMatchType", m_deltaRMatchType = 1);
	declareProperty("UseEmThresholdsOnly", m_useEmThresholdsOnly = true);
	declareProperty("UseEmTransRegionCut", m_useEmTRcut = true);
	declareProperty("IsEmType", m_isEmType = 32);

	for (int i = 0; i < ROI_BITS; ++i) {
		m_h_ClusterRaw_Et_bitcheck[i] = 0;
		m_h_ClusterRaw_Et_bitcheck_Eff[i] = 0;
		m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig[i] = 0;
		//m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig[i] = 0;
		m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig_Eff[i] = 0;
		//m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig_Eff[i] = 0;
		m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig[i] = 0;
		//m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig[i] = 0;
		m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig_Eff[i] = 0;
		//m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig_Eff[i] = 0;
		m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig[i] = 0;
		//m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig[i] = 0;
		m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig_Eff[i] = 0;
		//m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig_Eff[i] = 0;
	}
}

/*---------------------------------------------------------*/
EmEfficienciesMonTool::~EmEfficienciesMonTool()
/*---------------------------------------------------------*/
{
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*---------------------------------------------------------*/
StatusCode EmEfficienciesMonTool::initialize()
/*---------------------------------------------------------*/
{
	msg(MSG::INFO) << "Initializing " << name() << " - package version " << PACKAGE_VERSION << endreq;

	StatusCode sc;

	sc = ManagedMonitorToolBase::initialize();
	if (sc.isFailure())
		return sc;

	sc = m_histTool.retrieve();
	if (sc.isFailure()) {
		msg(MSG::ERROR) << "Unable to locate Tool TrigT1CaloHistogramTool" << endreq;
		return sc;
	}

	sc = m_tools.retrieve();
	if (sc.isFailure()) {
		msg(MSG::ERROR) << "Cannot retrieve L1CaloOfflineTriggerTowerTools" << endreq;
		return sc;
	}

	// Load in trigger tools
	// This tool gives you access to the physics trigger menu decisions for each event.
	sc = m_trigger.retrieve();
	if (sc.isFailure()) {
		msg(MSG::ERROR) << "Can't get handle on TrigDecisionTool" << endreq;
		return sc;
	}

	std::vector<std::string>::const_iterator iter = m_triggerStrings.begin();
	std::vector<std::string>::const_iterator iterE = m_triggerStrings.end();
	msg(MSG::INFO) << "TriggerStrings:";
	for (; iter != iterE; ++iter)
		msg(MSG::INFO) << " " << *iter;
	msg(MSG::INFO) << endreq;

	return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode EmEfficienciesMonTool::bookHistograms(bool isNewEventsBlock,
		bool isNewLumiBlock, bool isNewRun)
/*---------------------------------------------------------*/
{
	msg(MSG::DEBUG) << "bookHistograms entered" << endreq;

	if (m_environment == AthenaMonManager::online) {
		// book histograms that are only made in the online environment...
	}

	if (m_dataType == AthenaMonManager::cosmics) {
		// book histograms that are only relevant for cosmics data...
	}

	if (isNewEventsBlock || isNewLumiBlock) {
	}

	bool isNewInterval = (m_testMerge) ? isNewLumiBlock : isNewRun;
	if (isNewInterval) {

		std::string dir(m_rootDir + "/Reco/EmEfficiencies");
		ManagedMonitorToolBase::Interval_t interval = (m_testMerge) ? lumiBlock : run;

		MonGroup monEmDead(this, dir + "/DeadOrBadChannels", expert, interval);
		MonGroup monClusterRawNum(this, dir + "/ClusterRaw_Et/numerator", expert, interval);
		MonGroup monClusterRawDen(this, dir + "/ClusterRaw_Et/denominator", expert, interval);
		MonGroup monClusterRawEff(this, dir + "/ClusterRaw_Et", expert, interval, "", "perBinEffPerCent");
		MonGroup monClusterRaw10GeVNum(this, dir + "/ClusterRaw_10GeV_EtaVsPhi/numerator", expert, interval);
		MonGroup monClusterRaw10GeVDen(this, dir + "/ClusterRaw_10GeV_EtaVsPhi/denominator", expert, interval);
		MonGroup monClusterRaw10GeVEff(this, dir + "/ClusterRaw_10GeV_EtaVsPhi", expert, interval, "", "perBinEffPerCent");
		MonGroup monClusterRaw20GeVNum(this, dir + "/ClusterRaw_20GeV_EtaVsPhi/numerator", expert, interval);
		MonGroup monClusterRaw20GeVDen(this, dir + "/ClusterRaw_20GeV_EtaVsPhi/denominator", expert, interval);
		MonGroup monClusterRaw20GeVEff(this, dir + "/ClusterRaw_20GeV_EtaVsPhi", expert, interval, "",  "perBinEffPerCent");
		MonGroup monClusterRaw30GeVNum(this, dir + "/ClusterRaw_30GeV_EtaVsPhi/numerator", expert, interval);
		MonGroup monClusterRaw30GeVDen(this, dir + "/ClusterRaw_30GeV_EtaVsPhi/denominator", expert, interval);
		MonGroup monClusterRaw30GeVEff(this, dir + "/ClusterRaw_30GeV_EtaVsPhi", expert, interval, "", "perBinEffPerCent");

		//Set up EMTAU thresholds array with threshold names
		std::string thrNum[ROI_BITS] = { "0", "1", "2", "3", "4", "5", "6", "7",
		                                 "8", "9", "10", "11", "12", "13", "14", "15" };
		m_emBitMask = m_histTool->thresholdMaskEm();
		TrigConf::L1DataDef def;
		std::vector < std::string > emL1t;
		m_histTool->thresholdNames(def.emType(), emL1t);
		int size = emL1t.size();
		for (int i = 0; i < ROI_BITS; ++i) {
			if (i < size)
				emL1t[i] = "L1_" + emL1t[i];
			else
				emL1t.push_back("L1_??");
		}

		m_histTool->setMonGroup(&monEmDead);

		m_h_TrigTower_emDeadChannel = m_histTool->book2F("TrigTower_emDeadChannel", "EM Trigger Towers with dead channels - #eta V #phi;#eta^{Raw} Cluster;#phi^{Raw} Cluster", 50, -2.5, 2.5, 64, -M_PI, M_PI);

		m_h_TrigTower_emBadCalo = m_histTool->book2F("TrigTower_emBadCalo", "EM Trigger Towers affected by Missing FEBs - #eta V #phi;#eta^{Raw} Cluster;#phi^{Raw} Cluster", 50, -2.5, 2.5, 64, -M_PI, M_PI);

		//Raw Cluster Histograms

		m_histTool->setMonGroup(&monClusterRawDen);

		m_h_ClusterRaw_Et_gdEta = m_histTool->book1F("ClusterRaw_Et_gdEta", "Raw Cluster E_{T};E_{T}^{Raw} Cluster [GeV];Clusters", 100, 0, 100);

		m_histTool->setMonGroup(&monClusterRawNum);

		m_h_ClusterRaw_Et_triggered_gdEta = m_histTool->book1F("ClusterRaw_Et_triggered_gdEta", "Raw Cluster E_{T} (Triggered);E_{T}^{Raw} Cluster [GeV];Clusters", 100, 0, 100);

		std::string name;
		std::string title;
		for (int i = 0; i < ROI_BITS; ++i) {
		    if (emType(i)) {
			name = "ClusterRaw_Et_bitcheck_" + thrNum[i];
			title = "Raw E_{T} for Triggered Clusters passing " + emL1t[i] + ";E_{T}^{raw} Cluster [GeV];Clusters";
			m_h_ClusterRaw_Et_bitcheck[i] = m_histTool->book1F(name, title, 100, 0, 100);
		    } else m_h_ClusterRaw_Et_bitcheck[i] = 0;
		}

		m_histTool->setMonGroup(&monClusterRawEff);

		m_h_ClusterRaw_Et_triggered_Eff	= m_histTool->book1F("ClusterRaw_Et_triggered_Eff", "Raw Cluster E_{T} (Triggered) Efficiency;E_{T}^{Raw} Cluster [GeV];Efficiency %", 100, 0, 100);
		for (int i = 0; i < ROI_BITS; ++i) {
		    if (emType(i)) {
			name = "ClusterRaw_Et_bitcheck_Eff_" + thrNum[i];
			title = "Raw E_{T} for Triggered Clusters passing " + emL1t[i] + " Efficiency;E_{T}^{raw} Cluster [GeV];Efficiency %";
			m_h_ClusterRaw_Et_bitcheck_Eff[i] = m_histTool->book1F(name, title, 100, 0, 100);
		    } else m_h_ClusterRaw_Et_bitcheck_Eff[i] = 0;
		}

		m_histTool->setMonGroup(&monClusterRaw10GeVNum);

		//Et Raw Cluster greater > 10 GeV
		for (int i = 0; i < ROI_BITS; ++i) {
		    if (emType(i)) {
			name = "ClusterRaw_10GeV_Eta_vs_Phi_trig_" + thrNum[i];
			title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 10 GeV);#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig[i] = m_histTool->book2F(name, title, 50, -2.5, 2.5, 64, -M_PI, M_PI);

			//name = "ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig_" + thrNum[i];
			//title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 10 GeV) - Bad Calo and Dead Channel Towers excluded;#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			//m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig[i] = m_histTool->book2F(name,title,50,-2.5,2.5,64,-M_PI,M_PI);
		    } else m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig[i] = 0;
		}

		m_histTool->setMonGroup(&monClusterRaw10GeVDen);

		m_h_ClusterRaw_10GeV_Eta_vs_Phi = m_histTool->book2F("ClusterRaw_10GeV_Eta_vs_Phi",
				"Raw Cluster #eta against #phi (E_{T}^{raw} > 10 GeV);#eta^{Raw} Cluster;#phi^{Raw} Cluster",
				50, -2.5, 2.5, 64, -M_PI, M_PI);

		m_histTool->setMonGroup(&monClusterRaw10GeVEff);

		for (int i = 0; i < ROI_BITS; ++i) {
		    if (emType(i)) {
			name = "ClusterRaw_10GeV_Eta_vs_Phi_trig_Eff_" + thrNum[i];
			title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 10 GeV) Efficiency (%);#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig_Eff[i] = m_histTool->book2F(name, title, 50, -2.5, 2.5, 64, -M_PI, M_PI);

			//name = "ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig_Eff_" + thrNum[i];
			//title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 10 GeV) Efficiency (%) - Bad Calo and Dead Channel Towers excluded;#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			//m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig_Eff[i] = m_histTool->book2F(name,title,50,-2.5,2.5,64,-M_PI,M_PI);
		    } else m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig_Eff[i] = 0;
		}

		m_histTool->setMonGroup(&monClusterRaw20GeVNum);

		//Et Raw Cluster greater > 20 GeV
		for (int i = 0; i < ROI_BITS; ++i) {
		    if (emType(i)) {
			name = "ClusterRaw_20GeV_Eta_vs_Phi_trig_" + thrNum[i];
			title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 20 GeV);#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig[i] = m_histTool->book2F(name, title, 50, -2.5, 2.5, 64, -M_PI, M_PI);

			//name = "ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig_" + thrNum[i];
			//title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 20 GeV) - Bad Calo and Dead Channel Towers excluded;#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			//m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig[i] = m_histTool->book2F(name,title,50,-2.5,2.5,64,-M_PI,M_PI);
		    } else m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig[i] = 0;
		}

		m_histTool->setMonGroup(&monClusterRaw20GeVDen);

		m_h_ClusterRaw_20GeV_Eta_vs_Phi	= m_histTool->book2F("ClusterRaw_20GeV_Eta_vs_Phi", 
				"Raw Cluster #eta against #phi (E_{T}^{raw} > 20 GeV);#eta^{Raw} Cluster;#phi^{Raw} Cluster", 
				50, -2.5, 2.5, 64, -M_PI, M_PI);

		m_histTool->setMonGroup(&monClusterRaw20GeVEff);

		for (int i = 0; i < ROI_BITS; ++i) {
		    if (emType(i)) {
			name = "ClusterRaw_20GeV_Eta_vs_Phi_trig_Eff_" + thrNum[i];
			title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 20 GeV) Efficiency (%);#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig_Eff[i] = m_histTool->book2F(name, title, 50, -2.5, 2.5, 64, -M_PI, M_PI);

			//name = "ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig_Eff_" + thrNum[i];
			//title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 20 GeV) Efficiency (%) - Bad Calo and Dead Channel Towers excluded;#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			//m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig_Eff[i] = m_histTool->book2F(name,title,50,-2.5,2.5,64,-M_PI,M_PI);
		    } else m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig_Eff[i] = 0;
		}

		m_histTool->setMonGroup(&monClusterRaw30GeVNum);

		//Et Raw Cluster greater > 30 GeV
		for (int i = 0; i < ROI_BITS; ++i) {
		    if (emType(i)) {
			name = "ClusterRaw_30GeV_Eta_vs_Phi_trig_" + thrNum[i];
			title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i]	+ " with E_{T}^{raw} > 30 GeV);#eta^{Raw} Cluster;#phi^{Raw} Cluster";
			m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig[i] = m_histTool->book2F(name, title, 50, -2.5, 2.5, 64, -M_PI, M_PI);

			//name = "ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig_" + thrNum[i];
			//title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 30 GeV) - Bad Calo and Dead Channel Towers excluded;#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			//m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig[i] = m_histTool->book2F(name,title,50,-2.5,2.5,64,-M_PI,M_PI);
		    } else m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig[i] = 0;
		}

		m_histTool->setMonGroup(&monClusterRaw30GeVDen);

		m_h_ClusterRaw_30GeV_Eta_vs_Phi = m_histTool->book2F("ClusterRaw_30GeV_Eta_vs_Phi",
				"Raw Cluster #eta V #phi (E_{T}^{raw} > 30 GeV);#eta^{Raw} Cluster;#phi^{Raw} Cluster",
				50, -2.5, 2.5, 64, -M_PI, M_PI);

		m_histTool->setMonGroup(&monClusterRaw30GeVEff);

		for (int i = 0; i < ROI_BITS; ++i) {
		    if (emType(i)) {
			name = "ClusterRaw_30GeV_Eta_vs_Phi_trig_Eff_" + thrNum[i];
			title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 30 GeV) Efficiency (%);#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig_Eff[i] = m_histTool->book2F(name, title, 50, -2.5, 2.5, 64, -M_PI, M_PI);

			//name = "ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig_Eff_" + thrNum[i];
			//title = "Raw Cluster #eta V #phi (Triggered on " + emL1t[i] + " with E_{T}^{raw} > 30 GeV) Efficiency (%) - Bad Calo and Dead Channel Towers excluded;#eta^{Raw} Cluster;#phi^{Raw} Cluster)";
			//m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig_Eff[i] = m_histTool->book2F(name,title,50,-2.5,2.5,64,-M_PI,M_PI);
		    } else m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig_Eff[i] = 0;
		}

		m_histTool->unsetMonGroup();

		// HSB - counters 
		m_numEvents = 0;
		m_numOffElec = 0;
		m_numOffPhot = 0;
		m_numOffElecInContainer = 0;
		m_numOffPhotInContainer = 0;
		m_numOffElecPassCuts = 0;
		m_numOffPhotPassCuts = 0;
		m_numOffElecTriggered = 0;
		m_numOffPhotTriggered = 0;
		m_numEmObjPassTrigger = 0;
		m_numEmObjTotal = 0;
	}

	msg(MSG::DEBUG) << "Leaving bookHistograms" << endreq;

	return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode EmEfficienciesMonTool::fillHistograms()
/*---------------------------------------------------------*/
{
	const bool debug = msgLvl(MSG::DEBUG);
	if (debug) msg(MSG::DEBUG) << "fillHistograms entered" << endreq;

	// Here we can use the trigger menu to decide if we want an event.
	bool useEvent = false;
	if (m_useTrigger == true) {
		typedef std::vector<std::string>::iterator Itr_s;
		for (Itr_s i = m_triggerStrings.begin(); i != m_triggerStrings.end(); ++i) {
			if (m_trigger->isPassed(*i) == true) {
				useEvent = true;
				if (debug) msg(MSG::DEBUG)<< "First requested trigger that fired is : "<< (*i) << " with prescale "<< m_trigger->getPrescale(*i);
				break;
			}
		}

		StatusCode sc = this->triggerChainAnalysis();
		if (sc.isFailure()) {
			if (debug) msg(MSG::DEBUG) << "Problem checking Trigger Chains" << endreq;
			return sc;
		}
	}

	else {
		useEvent = true;
	}

	if (useEvent == true) {
		++m_numEvents;

		StatusCode sc;

		sc = this->loadContainers();
		if (sc.isFailure()) {
			msg(MSG::WARNING) << "Problem loading Athena Containers" << endreq;
			return sc;
		}
		
		if (debug) msg(MSG::DEBUG) << "Run number "<< m_eventInfo->event_ID()->run_number()<< " : Lumi Block "<< m_eventInfo->event_ID()->lumi_block() << " : Event "<< m_eventInfo->event_ID()->event_number() << endreq;
		
		if (m_numEvents == 1) {
			this->triggerTowerAnalysis();
		}

		// Look at vertex requirements
		int numVtx = 0, numTrk = 0;
		if (vertexRequirementsPassed(numVtx, numTrk) == false) {
			if (debug) msg(MSG::DEBUG) << "Event " << m_eventInfo->event_ID()->event_number() << " fails vertex requirements " << endreq;
			return StatusCode::SUCCESS;
		}

		m_numOffElecInContainer = 0;
		m_numOffElecPassCuts = 0;
		m_numOffElecTriggered = 0;
		sc = this->analyseOfflineElectrons();
		if (sc.isFailure()) {
			msg(MSG::WARNING) << "analyseElectrons Failed " << endreq;
			//return sc;
		}
		if (debug) msg(MSG::DEBUG) << "Number of Offline Electrons = "<< m_numOffElecInContainer << " Passing Cuts = "<< m_numOffElecPassCuts << " Triggered = "<< m_numOffElecTriggered << endreq;

		m_numOffPhotInContainer = 0;
		m_numOffPhotPassCuts = 0;
		m_numOffPhotTriggered = 0;
		sc = this->analyseOfflinePhotons();
		if (sc.isFailure()) {
			msg(MSG::WARNING) << "analysePhotons Failed " << endreq;
			//return sc;
		}
		if (debug) msg(MSG::DEBUG) << "Number of Offline Photons = "<< m_numOffPhotInContainer << " Passing Cuts = "<< m_numOffPhotPassCuts << " Triggered = "<< m_numOffPhotTriggered << endreq;

	}

	if (debug)
		msg(MSG::DEBUG) << "Leaving fillHistograms" << endreq;

	return StatusCode::SUCCESS;

}

/*---------------------------------------------------------*/
StatusCode EmEfficienciesMonTool::procHistograms(bool isEndOfEventsBlock,
		bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
	msg(MSG::DEBUG) << "procHistograms entered" << endreq;

	if (isEndOfEventsBlock || isEndOfLumiBlock || isEndOfRun) {
	}

	bool isEndOfInterval = (m_testMerge) ? isEndOfLumiBlock : isEndOfRun;
	if (isEndOfInterval) {
		msg(MSG::DEBUG) << "Number of offline electrons = " << m_numOffElec << endreq;
		msg(MSG::DEBUG) << "Number of offline photons = " << m_numOffPhot << endreq;
		msg(MSG::DEBUG) << "Number of events = " << m_numEvents << endreq;

		efficienciesForMerge(m_h_ClusterRaw_Et_gdEta,
				m_h_ClusterRaw_Et_triggered_gdEta,
				m_h_ClusterRaw_Et_triggered_Eff);

		for (int i = 0; i < ROI_BITS; ++i) {
			if (!emType(i)) continue;
			efficienciesForMerge(m_h_ClusterRaw_Et_gdEta,
					m_h_ClusterRaw_Et_bitcheck[i],
					m_h_ClusterRaw_Et_bitcheck_Eff[i]);
			efficienciesForMerge(m_h_ClusterRaw_10GeV_Eta_vs_Phi,
					m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig[i],
					m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig_Eff[i]);
			//efficienciesForMerge(m_h_ClusterRaw_10GeV_Eta_vs_Phi,
			//m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig[i],
			//m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig_Eff[i]);
			efficienciesForMerge(m_h_ClusterRaw_20GeV_Eta_vs_Phi,
					m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig[i],
					m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig_Eff[i]);
			//efficienciesForMerge(m_h_ClusterRaw_20GeV_Eta_vs_Phi,
			//m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig[i],
			//m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig_Eff[i]);
			efficienciesForMerge(m_h_ClusterRaw_30GeV_Eta_vs_Phi,
					m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig[i],
					m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig_Eff[i]);
			//efficienciesForMerge(m_h_ClusterRaw_30GeV_Eta_vs_Phi,
			//m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig[i],
			//m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig_Eff[i]);
		}
	}

	return StatusCode::SUCCESS;
}

void EmEfficienciesMonTool::efficienciesForMerge(LWHist* lw1, LWHist* lw2,
		LWHist* lw3) {
	TH1* hist1 = lw1->getROOTHistBase();
	TH1* hist2 = lw2->getROOTHistBase();
	TH1* hist3 = lw3->getROOTHistBase();
	// Need errors for Tier0 merge to work correctly
	if (!hist1->GetSumw2()->GetSize()) hist1->Sumw2();
	hist2->Sumw2();
	hist3->Sumw2();
	hist3->Divide(hist2, hist1, 1, 1, "B");
	// Modify errors to suit merging algorithm
	const double OneSigOneSided = 0.159;
	int nbins = hist3->GetNbinsX() + 2;
	if (hist3->GetDimension() == 2)
		nbins *= (hist3->GetNbinsY() + 2);
	for (int bin = 0; bin < nbins; ++bin) {
		double denom = hist1->GetBinContent(bin);
		if (denom == 0.)
			continue;
		double eff = hist3->GetBinContent(bin);
		if (eff == 0. || (float) eff > 0.99) {
			hist3->SetBinError(bin, 1. - std::pow(OneSigOneSided, 1. / denom));
		} else {
			hist3->SetBinError(bin, std::sqrt(eff * (1. - eff) / denom));
		}
	}
	hist3->Scale(100.); // Merge needs % Efficiency
	// Set plot limits - may be ignored by DQ
	if (hist3->GetDimension() == 1) {
		hist3->SetMinimum(0.);
		hist3->SetMaximum(110.);
	}
}

bool EmEfficienciesMonTool::emType(int bitNumber) {
        return (m_emBitMask>>bitNumber)&0x1;
}

/**********************************************/
//Analysis code for offline reconstructed electrons
//Compares electrons with EmTau RoIs
/**********************************************/
StatusCode EmEfficienciesMonTool::analyseOfflineElectrons() {
	//Access all of the offline reconstructed electron candidates in the event
	typedef ElectronContainer::const_iterator Itr_electrons;
	m_numOffElecInContainer = m_offlineElectrons->size();
	m_numEmObjTotal += m_numOffElecInContainer;

	// Create variables for electron properties
	double etaOE = 0.0, phiOE = 0.0, /*EtOE = 0.0, EtCE = 0.0,*/phiCE = 0.0, etaCE = 0.0;
	double EtCEraw = 0.0, phiCEraw = 0.0, etaCEraw = 0.0/*, calRawRatio = 0.0*/;
	// Create variable to determine if selecting the right type of electrons based on criteria in jO
	bool correctType;

	//bool roiValuesFilled = false;

	//Cycle through all of the offline reconstructed electrons
	for (Itr_electrons elItr = m_offlineElectrons->begin(); elItr != m_offlineElectrons->end(); ++elItr) {
		//Keep track of eta, phi and Et as these will be used often
		//----------------------------------------------------------------------
		//EtOE  = (*elItr)->et()/GeV;
		etaOE = (*elItr)->eta();
		phiOE = (*elItr)->phi();
		//----------------------------------------------------------------------
		//EtCE  = (*elItr)->cluster()->et()/GeV;
		etaCE = (*elItr)->cluster()->eta();
		phiCE = (*elItr)->cluster()->phi();
		//----------------------------------------------------------------------
		std::vector<double> rawValues = getRawClusterValuesFromCells(
				const_cast<CaloCluster*> ((*elItr)->cluster()));
		EtCEraw = rawValues.at(0);
		etaCEraw = rawValues.at(1);
		phiCEraw = rawValues.at(2);
		//----------------------------------------------------------------------
		//calRawRatio = (EtCEraw > 0.0) ? EtCE/EtCEraw : -1;

		bool unbiasedTrigger = false;
		if (m_passed_EF_Trigger == true) {
			if (m_passed_EF_SingleJet_Trigger == true) {
				unbiasedTrigger = isolatedEmObjectEF(phiCEraw, etaCEraw);
			} else if (m_passed_EF_egTau_Trigger == true) {
				if (m_passed_EF_SingleJet_Trigger == true) {
					unbiasedTrigger = isolatedEmObjectEF(phiCEraw, etaCEraw);
				} else {
					return StatusCode::SUCCESS;
				}
			} else {
				//unbiasedTrigger = isolatedEmObjectL1(phiCEraw, etaCEraw);
				unbiasedTrigger = true;
			}
		} else {
			unbiasedTrigger = true;
		}

		//If passed the trigger conditions then proceed to start analysis
		if (unbiasedTrigger == true) {

			correctType = false;
			m_numEmObjPassTrigger++;
			m_numOffElec++;

			//Check that the electron matches the IsEm type selected (if any was chosen in jobOptions file)
			correctType = correctIsEmElectron((*elItr));

			// Ask electron which is the highest isEm definition it has passed 
			int isEmCode = 0;
			std::string isEmLevel = isEmLevelElectron((*elItr), isEmCode);

			//Check if the reconstructed electron is reconstructed as a standard egamma object (not forward or softe)
			bool goodAuthor = ((*elItr)->author(egammaParameters::AuthorElectron) == true) ? true : false;
			if (goodAuthor == false) {
				correctType = false;
			}

			//Have the correct type of electron so do some analysis with it   
			if (correctType) {
				//Update counters		
				m_numOffElecPassCuts++;

				if (EtCEraw > 0.0) {

					if (inEgammaGoodEtaRange(etaCEraw, "el")) {
						m_h_ClusterRaw_Et_gdEta->Fill(EtCEraw);

						if (EtCEraw > 10) {
							m_h_ClusterRaw_10GeV_Eta_vs_Phi->Fill(etaCEraw, phiCEraw);
						}
						if (EtCEraw > 20) {
							m_h_ClusterRaw_20GeV_Eta_vs_Phi->Fill(etaCEraw, phiCEraw);
						}
						if (EtCEraw > 30) {
							m_h_ClusterRaw_30GeV_Eta_vs_Phi->Fill(etaCEraw, phiCEraw);
						}
					}
				}

				//Set up useful numbers to keep track of RoI information
				double etaROI = 0.0, phiROI = 0.0/*, EtROI = 0.0*/;
				double dEta = 0.0, dPhi = 0.0, dR = 1001, tempRmin = 0.0;
				double dEtaClus = 0.0, dPhiClus = 0.0, dRClus = 1000, tempRminClus = 0.0;
				double dEtaClRaw = 0.0, dPhiClRaw = 0.0, dRClRaw = 1000, tempRminClRaw = 0.0;
				//double bestEtaROI = 0.0, bestPhiROI = 0.0, bestEtROI = 0.0, bestEtIsol = 0.0;
				//double bestEtResClus = 1000, bestEtResClusRaw = 1000;
				//double bestDeltaPhi = 0.0, bestDeltaEta = 0.0, bestDeltaEt = 0.0;
				//double bestDeltaPhiClus = 0.0, bestDeltaEtaClus = 0.0, bestDeltaEtClus = 0.0;
				//double bestDeltaPhiClRaw = 0.0, bestDeltaEtaClRaw = 0.0, bestDeltaEtClRaw = 0.0;
				uint32_t ROIWord = 0/*, ThrPattern = 0*/;

				//Access the EmTauRoIs
				std::vector < EmTau_ROI > emrois = m_lvl1RoIs->getEmTauROIs();
				typedef std::vector<EmTau_ROI>::const_iterator Itr_emroi;

				//Iterate over all of the EmTauRoIs
				for (Itr_emroi roiItr = emrois.begin(); roiItr != emrois.end(); ++roiItr) {
					bool emThresholdPassed = false;
					if (m_useEmThresholdsOnly == true) {

						std::vector < std::string > thrPassed = (*roiItr).getThresholdNames();
						typedef std::vector<std::string>::iterator Itr_s;
						for (Itr_s i = thrPassed.begin(); i != thrPassed.end(); ++i) {
							if ((*i).find("EM") != std::string::npos) {
								emThresholdPassed = true;
								break;
							}
						}
					} else {
						emThresholdPassed = true;
					}

					if (emThresholdPassed == true) {

						//Get useful values for the EmTauRoI
						etaROI = (*roiItr).getEta();
						phiROI = (*roiItr).getPhi();
						//EtROI  = (*roiItr).getEMClus()/GeV;

						//Calculate the difference in eta and phi between the electron and RoI
						dEta = etaOE - etaROI;
						dPhi = correctDeltaPhi(phiOE - phiROI);
						dEtaClus = etaCE - etaROI;
						dPhiClus = correctDeltaPhi(phiCE - phiROI);
						dEtaClRaw = etaCEraw - etaROI;
						dPhiClRaw = correctDeltaPhi(phiCEraw - phiROI);

						//Calculate delta R
						tempRmin = sqrt(dEta * dEta + dPhi * dPhi);
						tempRminClus = sqrt(dEtaClus * dEtaClus + dPhiClus * dPhiClus);
						tempRminClRaw = sqrt(dEtaClRaw * dEtaClRaw + dPhiClus * dPhiClRaw);

						double tempdR = 0., smallestdRSoFar = 0.;
						switch (m_deltaRMatchType) {
						case 0: // Calibrated Clusters
							tempdR = tempRminClus;
							smallestdRSoFar = dRClus;
							break;
						case 1: // Raw Clusters
							tempdR = tempRminClRaw;
							smallestdRSoFar = dRClRaw;
							break;
						case 2: // Offline Objects
							tempdR = tempRmin;
							smallestdRSoFar = dR;
							break;
						default:
							break;
						}

						//Check if the new delta R is smaller than any previous delta R value.
						//In that case, keep track of the new RoI values
						if (tempdR < smallestdRSoFar) {
							//RoI information
							//bestPhiROI = phiROI;
							//bestEtaROI = etaROI;
							//bestEtROI = EtROI;
							//bestEtIsol = (*roiItr).getEMIsol()/GeV;
							ROIWord = (*roiItr).getROIWord();
							//ThrPattern = (*roiItr).getThrPattern();

							//bestDeltaEta = dEta;
							//bestDeltaPhi = dPhi;
							dR = tempRmin;
							//bestDeltaEt = EtOE-bestEtROI;					   

							//bestDeltaEtaClus = dEtaClus;
							//bestDeltaPhiClus = dPhiClus;
							dRClus = tempRminClus;
							//bestDeltaEtClus = EtCE-bestEtROI;
							//bestEtResClus = (EtCE-bestEtROI)/EtCE;

							//bestDeltaEtaClRaw = dEtaClRaw;
							//bestDeltaPhiClRaw = dPhiClRaw;
							dRClRaw = tempRminClRaw;
							//bestDeltaEtClRaw = EtCEraw-bestEtROI;					     
							//if( EtCEraw > 0.0 ) { bestEtResClusRaw = (EtCEraw-bestEtROI)/EtCEraw; }
						}
					}
				}

				//roiValuesFilled = true;

				//Check to see if there was an RoI to match with an electron cluster
				if (dRClus != 1000) {
					m_numOffElecTriggered++;

					//Check if electron and RoI matched to a very good level (less than cut)
					//if(deltaMatch(dR, bestDeltaEta, bestDeltaPhi, m_goodEMDeltaRMatch_Cut, m_goodEMDeltaEtaMatch_Cut, m_goodEMDeltaPhiMatch_Cut))
					if (dRClus < m_goodEMDeltaRMatch_Cut) {
						if (EtCEraw > 0.0) {
							if (inEgammaGoodEtaRange(etaCEraw, "el")) {
								m_h_ClusterRaw_Et_triggered_gdEta->Fill(EtCEraw);

								//Look at each bit in the RoI word (using bitshift) to see which thresholds were
								//passed and which ones were not
								for (int k = 0; k < ROI_BITS; ++k) {
									if (emType(k) && ((ROIWord >> k) & 1)) {
										m_h_ClusterRaw_Et_bitcheck[k]->Fill(EtCEraw);

										//bool deadBadTower = emObjInDeadBadTower(etaCEraw,phiCEraw);
										if (EtCEraw > 10) {
											m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig[k]->Fill(etaCEraw, phiCEraw);
											//if(deadBadTower == false) { m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig[k]->Fill(etaCEraw,phiCEraw); }

										}
										if (EtCEraw > 20) {
											m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig[k]->Fill(etaCEraw, phiCEraw);
											//if(deadBadTower == false) { m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig[k]->Fill(etaCEraw,phiCEraw); }
										}
										if (EtCEraw > 30) {
											m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig[k]->Fill(etaCEraw, phiCEraw);
											//if(deadBadTower == false) { m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig[k]->Fill(etaCEraw,phiCEraw); }
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return StatusCode::SUCCESS;
}

//----------------------------------------------------------------------
//Analysis code for offline reconstructed photons
//Compares photons with EmTau RoIs
//----------------------------------------------------------------------
StatusCode EmEfficienciesMonTool::analyseOfflinePhotons() {
	// Access all of the offline reconstructed photon candidates in the event
	typedef PhotonContainer::const_iterator Itr_photons;
	m_numOffPhotInContainer = m_offlinePhotons->size();
	m_numEmObjTotal += m_numOffPhotInContainer;

	// Variables for accessing properties of recosntructed photons
	double etaOP = 0.0, phiOP = 0.0, /*EtOP = 0.0, EtCP = 0.0,*/etaCP = 0.0, phiCP = 0.0;
	double EtCPraw = 0.0, etaCPraw = 0.0, phiCPraw = 0.0/*, calRawRatio = 0.0*/;
	// Variable to check if photon is of the right type as defined in the jobOptions
	bool correctType;

	//bool roiValuesFilled = false;

	//Cycle through all of the offline reconstructed photons      
	for (Itr_photons phItr = m_offlinePhotons->begin(); phItr != m_offlinePhotons->end(); ++phItr) {
		//Keep track of eta, phi and Et as these will be used often
		//----------------------------------------------------------------------
		//EtOP  = (*phItr)->et()/GeV;
		etaOP = (*phItr)->eta();
		phiOP = (*phItr)->phi();
		//----------------------------------------------------------------------
		//EtCP  = (*phItr)->cluster()->et()/GeV;
		etaCP = (*phItr)->cluster()->eta();
		phiCP = (*phItr)->cluster()->phi();
		//----------------------------------------------------------------------
		std::vector<double> rawValues = getRawClusterValuesFromCells(const_cast<CaloCluster*> ((*phItr)->cluster()));
		EtCPraw = rawValues.at(0);
		etaCPraw = rawValues.at(1);
		phiCPraw = rawValues.at(2);
		//----------------------------------------------------------------------	
		//calRawRatio = (EtCPraw > 0.0) ? EtCP/EtCPraw : -1;

		//Check, based on trigger chain analysis, that we are using an unbiased trigger
		bool unbiasedTrigger = false;
		//Check that event has passed EF trigger
		if (m_passed_EF_Trigger == true) {
			//Ask if passed EF single jet trigger then check if egamma object near triggered jet
			if (m_passed_EF_SingleJet_Trigger == true) {
				unbiasedTrigger = isolatedEmObjectEF(phiCPraw, etaCPraw);
			} else if (m_passed_EF_egTau_Trigger == true) {
				if (m_passed_EF_SingleJet_Trigger == true) {
					unbiasedTrigger = isolatedEmObjectEF(phiCPraw, etaCPraw);
				} else {
					return StatusCode::SUCCESS;
				}
			} else {
				//unbiasedTrigger = isolatedEmObjectL1(phiCPraw, etaCPraw);
				unbiasedTrigger = true;
			}
		} else {
			unbiasedTrigger = true;
		}

		//If only after the highest energy photon, make sure only select the one that has the right index number
		//Otherwise select all of them	
		if (unbiasedTrigger == true) {
			//Check that the photon matches the IsEm type selected (if any was chosen in jobOptions file)
			correctType = false;
			m_numEmObjPassTrigger++;
			m_numOffPhot++;

			correctType = correctIsEmPhoton((*phItr));

			// Ask photon what is the highest isEm definition that was passed 
			int isEmCode = 0;
			std::string isEmLevel = isEmLevelPhoton((*phItr), isEmCode);

			//Check if the reconstructed photon is reconstructed as a standard egamma object 
			bool goodAuthor = ((*phItr)->author(egammaParameters::AuthorPhoton)
					== true) ? true : false;
			if (goodAuthor == false) {
				correctType = false;
			}

			// Ask is there a conversion - no requirement for photons yet

			//Have the correct type of photon so do some analysis with it    
			if (correctType) {
				//Update counters		
				m_numOffPhotPassCuts++;

				if (EtCPraw > 0.0) {
					if (inEgammaGoodEtaRange(etaCPraw, "ph")) {
						m_h_ClusterRaw_Et_gdEta->Fill(EtCPraw);

						if (EtCPraw > 10) {
							m_h_ClusterRaw_10GeV_Eta_vs_Phi->Fill(etaCPraw, phiCPraw);
						}
						if (EtCPraw > 20) {
							m_h_ClusterRaw_20GeV_Eta_vs_Phi->Fill(etaCPraw, phiCPraw);
						}
						if (EtCPraw > 30) {
							m_h_ClusterRaw_30GeV_Eta_vs_Phi->Fill(etaCPraw, phiCPraw);
						}
					}
				}

				//Set up useful numbers to keep track of RoI information
				double etaROI = 0.0, phiROI = 0.0/*, EtROI = 0.0*/;
				double dEta = 0.0, dPhi = 0.0, dR = 1001, tempRmin = 0.0;
				double dEtaClus = 0.0, dPhiClus = 0.0, dRClus = 1000, tempRminClus = 0.0;
				double dEtaClRaw = 0.0, dPhiClRaw = 0.0, dRClRaw = 1000, tempRminClRaw = 0.0;
				//double bestPhiROI = 0.0, bestEtaROI = 0.0, bestEtROI = 0.0, bestEtIsol = 0.0; 
				//double bestEtResClus = 1000, bestEtResClusRaw = 1000;
				//double bestDeltaPhi = 0.0, bestDeltaEta = 0.0, bestDeltaEt = 0.0;	   
				//double bestDeltaPhiClus = 0.0, bestDeltaEtaClus = 0.0, bestDeltaEtClus = 0.0;	
				//double bestDeltaPhiClRaw = 0.0, bestDeltaEtaClRaw = 0.0, bestDeltaEtClRaw = 0.0;				
				uint32_t ROIWord = 0/*, ThrPattern = 0*/;

				//Access the EmTau RoIs
				std::vector < EmTau_ROI > emrois = m_lvl1RoIs->getEmTauROIs();
				typedef std::vector<EmTau_ROI>::const_iterator Itr_emroi;

				//Iterate over the EmTau RoIs 
				for (Itr_emroi roiItr = emrois.begin(); roiItr != emrois.end(); ++roiItr) {
					bool emThresholdPassed = false;
					if (m_useEmThresholdsOnly == true) {
						std::vector < std::string > thrPassed = (*roiItr).getThresholdNames();
						typedef std::vector<std::string>::iterator Itr_s;
						for (Itr_s i = thrPassed.begin(); i != thrPassed.end(); ++i) {
							if ((*i).find("EM") != std::string::npos) {
								emThresholdPassed = true;
								break;
							}
						}
					} else {
						emThresholdPassed = true;
					}

					if (emThresholdPassed == true) {
						//Get useful values for the EmTauRoI
						etaROI = (*roiItr).getEta();
						phiROI = (*roiItr).getPhi();
						//EtROI  = (*roiItr).getEMClus()/GeV;

						//Calculate the difference in eta and phi between the electron and RoI
						dEta = etaOP - etaROI;
						dPhi = correctDeltaPhi(phiOP - phiROI);
						dEtaClus = etaCP - etaROI;
						dPhiClus = correctDeltaPhi(phiCP - phiROI);
						dEtaClRaw = etaCPraw - etaROI;
						dPhiClRaw = correctDeltaPhi(phiCPraw - phiROI);

						//Calculate deltaR
						tempRmin = sqrt(dEta * dEta + dPhi * dPhi);
						tempRminClus = sqrt(dEtaClus * dEtaClus + dPhiClus * dPhiClus);
						tempRminClRaw = sqrt(dEtaClRaw * dEtaClRaw + dPhiClus * dPhiClRaw);

						double tempdR = 0., smallestdRSoFar = 0.;
						switch (m_deltaRMatchType) {
						case 0:
							tempdR = tempRminClus;
							smallestdRSoFar = dRClus;
							break;
						case 1:
							tempdR = tempRminClRaw;
							smallestdRSoFar = dRClRaw;
							break;
						case 2:
							tempdR = tempRmin;
							smallestdRSoFar = dR;
							break;
						default:
							break;
						}

						//Check if the new deltaR is smaller than any previous deltaR value.
						//In that case, keep track of the new RoI values 
						if (tempdR < smallestdRSoFar) {
							//bestPhiROI = phiROI;
							//bestEtaROI = etaROI;
							//bestEtROI  = EtROI;
							//bestEtIsol = (*roiItr).getEMIsol()/GeV;
							ROIWord = (*roiItr).getROIWord();
							//ThrPattern = (*roiItr).getThrPattern();

							//bestDeltaEta = dEta;
							//bestDeltaPhi = dPhi;
							dR = tempRmin;
							//bestDeltaEt  = EtOP-bestEtROI;

							//bestDeltaEtaClus = dEtaClus;
							//bestDeltaPhiClus = dPhiClus;
							dRClus = tempRminClus;
							//bestDeltaEtClus = EtCP-bestEtROI;
							//bestEtResClus = (EtCP-bestEtROI)/EtCP;

							//bestDeltaEtaClRaw = dEtaClRaw;
							//bestDeltaPhiClRaw = dPhiClRaw;
							dRClRaw = tempRminClRaw;
							//bestDeltaEtClRaw = EtCPraw-bestEtROI;
							//if(EtCPraw > 0.0) { bestEtResClusRaw = (EtCPraw-bestEtROI)/EtCPraw; }
						}
					}
				}

				//roiValuesFilled = true;

				//Check to see if there was an RoI to match with a photon
				if (dRClus != 1000) {
					m_numOffPhotTriggered++;

					//Check if photon and RoI matched to a very good level (less than cut)
					//if(deltaMatch(dR, bestDeltaEta, bestDeltaPhi, m_goodEMDeltaRMatch_Cut, m_goodEMDeltaEtaMatch_Cut, m_goodEMDeltaPhiMatch_Cut))
					if (dRClus < m_goodEMDeltaRMatch_Cut) {
						if (EtCPraw > 0.0) {
							if (inEgammaGoodEtaRange(etaCPraw, "ph")) {
								m_h_ClusterRaw_Et_triggered_gdEta->Fill(EtCPraw);

								//Look at each bit in the RoI word (using bitshift) to see which thresholds were
								//passed and which ones were not
								for (int k = 0; k < ROI_BITS; ++k) {
									if (emType(k) && ((ROIWord >> k) & 1)) {
										m_h_ClusterRaw_Et_bitcheck[k]->Fill(EtCPraw);

										//bool deadBadTower = emObjInDeadBadTower(etaCPraw,phiCPraw);
										if (EtCPraw > 10) {
											m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig[k]->Fill(etaCPraw, phiCPraw);
											//if(deadBadTower == false) { m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig[k]->Fill(etaCPraw,phiCPraw); }
										}
										if (EtCPraw > 20) {
											m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig[k]->Fill(etaCPraw, phiCPraw);
											//if(deadBadTower == false) { m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig[k]->Fill(etaCPraw,phiCPraw); }
										}
										if (EtCPraw > 30) {
											m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig[k]->Fill(etaCPraw, phiCPraw);
											//if(deadBadTower == false) { m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig[k]->Fill(etaCPraw,phiCPraw); }
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return StatusCode::SUCCESS;
}

//------------------------------------------------------------------
// The check to see if an object is triggered w.r.t. an RoI
// It can be done in two ways so allow it to handle either one
//------------------------------------------------------------------
bool EmEfficienciesMonTool::deltaMatch(double dEta, double dPhi, double dR, 
		double goodDEta, double goodDPhi, double goodDR) {
	
	// Calculate if object passes cuts
	bool Rmatch = (dR < goodDR) ? true : false;
	bool EPmatch = ((fabs(dEta) < goodDEta) && (fabs(dPhi) < goodDPhi)) ? true
			: false;

	// First check that both checks are either on or off
	if (m_useDeltaRMatch == true && m_useDeltaEtaPhiMatch == true) {
		if (Rmatch == true && EPmatch == true) {
			return true;
		} else {
			return false;
		}
	} else if (m_useDeltaRMatch == false && m_useDeltaEtaPhiMatch == false) {
		return false;
	}

	// Only one check is being used so match against it
	if (m_useDeltaRMatch == true) {
		if (Rmatch == true) {
			return true;
		} else {
			return false;
		}
	} else if (m_useDeltaEtaPhiMatch == true) {
		if (EPmatch == true) {
			return true;
		} else {
			return false;
		}
	} else {
		return false;
	}
}

//------------------------------------------------------------------
//Asks if there is at least one primary vertex with 3 tracks coming off it
//Useful for selecting collision events
//------------------------------------------------------------------
bool EmEfficienciesMonTool::vertexRequirementsPassed(int &numVtx, int &bestNumTracks) {
	
	bool goodVtx = false;

	//See if we can find the requested primary vertex collection    
	numVtx = m_primaryVtx->size();
	bestNumTracks = 0;

	//See if there are any vertices in the collection
	if (numVtx > 0) {
		//Set up iterators
		VxContainer::const_iterator vertexItr = m_primaryVtx->begin();
		VxContainer::const_iterator vertexItrE = m_primaryVtx->end();

		//Loop over vertices
		for (; vertexItr != vertexItrE; ++vertexItr) {
			//Find out if the vertex has at least 3 tracks coming from it
			std::vector<Trk::VxTrackAtVertex*>* trklist = (*vertexItr)->vxTrackAtVertex();
			int numTracks = trklist->size();

			if (numTracks > bestNumTracks) {
				bestNumTracks = numTracks;
			}

			if (numTracks >= 3) {
				//Met our requirement
				goodVtx = true;
				break;
			}
		}
	}

	return goodVtx;
}

//------------------------------------------------------------------
//Asks if there is at least one primary vertex with 3 tracks coming off it
//Useful for selecting collision events
//------------------------------------------------------------------
bool EmEfficienciesMonTool::emObjInDeadBadTower(double eta, double phi) {
	
	bool deadBadTower = false;
	int etaBin = 0, phiBin = 0;

	phiBin = 1 + (int) ((phi + M_PI) * 64.0 / (2 * M_PI));
	etaBin = 1 + (int) ((eta + 2.5) * 50.0 / (2 * 2.5)); //this will need fixing

	bool badCalo = (m_h_TrigTower_emBadCalo->GetBinContent(etaBin, phiBin) > 0) ? true	: false;
	bool deadChannel = (m_h_TrigTower_emDeadChannel->GetBinContent(etaBin, phiBin) > 0) ? true : false;

	if (badCalo || deadChannel) {
		deadBadTower = true;
	}
	return deadBadTower;
}

//------------------------------------------------------------------
// Get the raw values of the CaloCluster if they are not provided
// directly by the cluster by accessing the cells that make it up
// Uses CaloCells so needs to use ESD files
//------------------------------------------------------------------
std::vector<double> EmEfficienciesMonTool::getRawClusterValuesFromCells(CaloCluster* cc) {
	
	// Add the raw information of the cluster 
	double rawE = 0., rawEta = 0., rawPhi = 0., rawEt = 0.;
	double cellE = 0., cellEta = 0., cellPhi = 0.;

	// Variables for asking if cells cross +/-M_PI boundary
	double phihi = -M_PI, philo = +M_PI, phiDiff = -1;

	CaloCluster::cell_iterator ccIt = cc->cell_begin();
	CaloCluster::cell_iterator ccItE = cc->cell_end();

	// Check if the cells in the cluster cross +/-M_PI boundary
	// by asking what is the biggest phi difference between any of
	// the cells belonging to the cluster
	for (; ccIt != ccItE; ++ccIt) {
		const CaloCell* cell = (*ccIt);
		if (cell) {
			cellPhi = (*ccIt)->phi();
			if (cellPhi > phihi) {
				phihi = cellPhi;
			}
			if (cellPhi < philo) {
				philo = cellPhi;
			}
		}
	}
	// Calculate biggest phi difference
	phiDiff = phihi - philo;

	// Loop over the cells corresponding to the cluster
	for (ccIt = cc->cell_begin(); ccIt != ccItE; ++ccIt) {
		const CaloCell* cell = (*ccIt);
		if (cell) {
			// Get the cell values
			cellE = (*ccIt)->energy();
			cellEta = (*ccIt)->eta();
			cellPhi = (*ccIt)->phi();

			// If the cells cross boundary, shift some so it is continuous
			if (phiDiff > 4 && cellPhi < 0) {
				cellPhi += 2 * M_PI;
			}

			//Add to the total energy sum
			rawE += cellE;
			//Add to the energy weighted eta and phi sum
			rawPhi += cellE * cellPhi;
			rawEta += cellE * cellEta;

		} else {
			msg(MSG::WARNING) << "Problem with Cell within Cluster, check cell pointer: " << cell << endreq;
		}
	}

	std::vector<double> rawV;

	if (rawE > 0) {
		//Unweight the raw eta and phi values
		rawEta /= rawE;
		rawPhi /= rawE;

		//If shifted phi value falls out of range, correct it
		if (phiDiff > 4 && rawPhi > M_PI) {
			rawPhi -= 2 * M_PI;
		}

		//Calculate the raw et from the energy and eta
		rawEt = rawE / (GeV * std::cosh(rawEta));

		rawV.push_back(rawEt);
		rawV.push_back(rawEta);
		rawV.push_back(rawPhi);
	} else {
		rawV.push_back(-10.0);
		rawV.push_back(-10.0);
		rawV.push_back(-10.0);
	}

	return rawV;

}

//------------------------------------------------------------------
//Correct the value of deltaPhi so that it falls in +-M_PI range
//------------------------------------------------------------------
double EmEfficienciesMonTool::correctDeltaPhi(double dPhi) {
	if (fabs(dPhi) > M_PI) {
		dPhi = (dPhi > 0) ? dPhi - 2 * M_PI : dPhi + 2 * M_PI;
	}

	return dPhi;
}

//------------------------------------------------------------------
//Calculate delta R quickly between two objects
//------------------------------------------------------------------
double EmEfficienciesMonTool::calcDeltaR(double eta1, double phi1, double eta2,	double phi2) {
	double dEta = eta1 - eta2;
	double dPhi = correctDeltaPhi(phi1 - phi2);

	double dR = sqrt((dEta * dEta) + (dPhi * dPhi));
	return dR;
}

//------------------------------------------------------------------
//Check for electron that it is of the right isEm type as required from jobOptions
//------------------------------------------------------------------
bool EmEfficienciesMonTool::correctIsEmElectron(const Analysis::Electron* el) {
	bool correctType = false;
	switch (m_isEmType) {
	case 0: //"None"
		correctType = true;
		break;
	case 10: //"ElectronLoose" 
		correctType = (el->isem(egammaPID::ElectronLoose) == 0) ? true : false;
		break;
	case 11: //"ElectronMedium"
		correctType = (el->isem(egammaPID::ElectronMedium) == 0) ? true : false;
		break;
	case 12: //"ElectronTight" 
		correctType = (el->isem(egammaPID::ElectronTight) == 0) ? true : false;
		break;
	case 30: //"ElectronLoose&PhotonLooseCombination" so just ask if electron loose
		correctType = (el->isem(egammaPID::ElectronLoose) == 0) ? true : false;
		break;
	case 31: //"ElectronMedium&PhotonLooseCombination" so just ask if electron medium
		correctType = (el->isem(egammaPID::ElectronMedium) == 0) ? true : false;
		break;
	case 32: //"ElectronTight&PhotonTightCombination" so just ask if electron tight
		correctType = (el->isem(egammaPID::ElectronTight) == 0) ? true : false;
		break;
	default:
		correctType = false;
		break;
	}

	return correctType;
}

//------------------------------------------------------------------
//Check for photon that it is of the right isEm type as required from jobOptions
//------------------------------------------------------------------
bool EmEfficienciesMonTool::correctIsEmPhoton(const Analysis::Photon* ph) {
	bool correctType = false;
	switch (m_isEmType) {
	case 0: //"None"
		correctType = true;
		break;
	case 20: //"PhotonLoose" 
		correctType = (ph->isem(egammaPID::PhotonLoose) == 0) ? true : false;
		break;
	case 21: //"PhotonTight" 
		correctType = (ph->isem(egammaPID::PhotonTight) == 0) ? true : false;
		break;
	case 30: //"PhotonLoose&ElectronLooseCombination" so just ask if photon loose
		correctType = (ph->isem(egammaPID::PhotonLoose) == 0) ? true : false;
		break;
	case 31: //"PhotonLoose&ElectronMediumCombination" so just ask if photon loose
		correctType = (ph->isem(egammaPID::PhotonLoose) == 0) ? true : false;
		break;
	case 32: //"PhotonTight&ElectronTightCombination" so just ask if photon tight
		correctType = (ph->isem(egammaPID::PhotonTight) == 0) ? true : false;
		break;
	default:
		correctType = false;
		break;
	}

	return correctType;
}

//------------------------------------------------------------------
//Check for electron that it is of the right isEm type as required from jobOptions
//------------------------------------------------------------------
std::string EmEfficienciesMonTool::isEmLevelElectron(const Analysis::Electron* el, int &code) {
	std::string isEmLevel = "None";
	code = 0;
	if (el->isem(egammaPID::ElectronLoose) == 0) {
		isEmLevel = "Loose";
		code = 1;
	}
	if (el->isem(egammaPID::ElectronMedium) == 0) {
		isEmLevel = "Medium";
		code = 2;
	}
	if (el->isem(egammaPID::ElectronTight) == 0) {
		isEmLevel = "Tight";
		code = 3;
	}

	return isEmLevel;
}

//------------------------------------------------------------------
//Check for photon that it is of the right isEm type as required from jobOptions
//------------------------------------------------------------------
std::string EmEfficienciesMonTool::isEmLevelPhoton(const Analysis::Photon* ph, int &code) {
	std::string isEmLevel = "None";
	code = 0;
	if (ph->isem(egammaPID::PhotonLoose) == 0) {
		isEmLevel = "Loose";
		code = 1;
	}
	if (ph->isem(egammaPID::PhotonTight) == 0) {
		isEmLevel = "Tight";
		code = 3;
	}

	return isEmLevel;
}

//------------------------------------------------------------------
// Ask if object is within well defined eta range for EM performance
// For egamma trig note this range may change
//------------------------------------------------------------------
bool EmEfficienciesMonTool::inEgammaGoodEtaRange(double eta, std::string /*egType*/) {
	//For performance studies, take all egamma within LAr region
	if (fabs(eta) > 2.50) {
		return false;
	}

	bool inGE = true;
	if (m_useEmTRcut) {
		inGE = (inEMTransR(eta, 0)) ? false : true;
	}
	return inGE;
}

//------------------------------------------------------------------
//Ask if object is within EM barrel
//------------------------------------------------------------------
bool EmEfficienciesMonTool::inEMBarrel(double eta, int sign) {
	bool inEB = (fabs(eta) <= 1.37) ? true : false;

	if (inEB == true) {
		if (sign < 0) {
			inEB = (eta < 0) ? true : false;
		} else if (sign > 0) {
			inEB = (eta > 0) ? true : false;
		}
	}

	return inEB;
}

//------------------------------------------------------------------
//Ask if object is within EM Transition region
//------------------------------------------------------------------
bool EmEfficienciesMonTool::inEMTransR(double eta, int sign) {
	bool inTR = (fabs(eta) > 1.37 && fabs(eta) < 1.52) ? true : false;

	if (inTR == true) {
		if (sign < 0) {
			inTR = (eta < 0) ? true : false;
		} else if (sign > 0) {
			inTR = (eta > 0) ? true : false;
		}
	}

	return inTR;
}

//------------------------------------------------------------------
//Ask if object is within EM endcap
//------------------------------------------------------------------
bool EmEfficienciesMonTool::inEMEndcap(double eta, int sign) {
	bool inEC = (fabs(eta) >= 1.52 && fabs(eta) < 2.50) ? true : false;

	if (inEC == true) {
		if (sign < 0) {
			inEC = (eta < 0) ? true : false;
		} else if (sign > 0) {
			inEC = (eta > 0) ? true : false;
		}
	}

	return inEC;
}

//------------------------------------------------------------------
//Ask if object is has no jets or jet RoIs nearby at L1
//Do this by assuming that the highest ET Jet RoI is the one causing
//the event to trigger at HLT level
//------------------------------------------------------------------
bool EmEfficienciesMonTool::isolatedEmObjectL1(double phi, double eta) {
	bool isolated = true;
	double dREm = 10.0, dR_Max = 10.0;
	double ET_Max = -10.0, EtaMax = -10.0, PhiMax = -10.0;
	double etaROI = 0.0, phiROI = 0.0, ET_ROI = 0.0;

	//Cycle over the rois, get their properties and determine the distance from the object
	const std::vector<Jet_ROI>& jetROIs(m_lvl1RoIs->getJetROIs());
	typedef std::vector<Jet_ROI>::const_iterator Itr_jetroi;
	for (Itr_jetroi roiItr = jetROIs.begin(); roiItr != jetROIs.end(); ++roiItr) {
		etaROI = (*roiItr).getEta();
		phiROI = (*roiItr).getPhi();
		ET_ROI = (*roiItr).getET8x8();
		dREm = calcDeltaR(eta, phi, etaROI, phiROI);

		//If this energy exceeds the current record then store the details
		if (ET_ROI > ET_Max) {
			EtaMax = etaROI;
			PhiMax = phiROI;
			ET_Max = ET_ROI;
			dR_Max = dREm;
		}
	}

	// Check that the object is far away enough from highest ET jet RoI
	if (dR_Max > m_goodHadDeltaRMatch_Cut) {
		isolated = true;
	} else {
		isolated = false;
	}

	return isolated;
}

//------------------------------------------------------------------
//Ask if object is has no jets or jet RoIs nearby at EF
//------------------------------------------------------------------
bool EmEfficienciesMonTool::isolatedEmObjectEF(double phi, double eta) {
	//At present just call L1 function, improve this later
	return isolatedEmObjectL1(phi, eta);
}

//---------------------------------------------------------------
// Trigger Tower Analysis
//---------------------------------------------------------------
void EmEfficienciesMonTool::triggerTowerAnalysis() {
	// We will loop over all TriggerTowers, so we will use an iterator.
	// We use the ++i convention rather than i++. With i++ a copy has to be made and is
	// therefore a lot slower. Don't really matter for ints, certainly does for iterators

	// Look at trigger towers around physics objects

	typedef TriggerTowerCollection::const_iterator Itr_tt;
	typedef std::vector<int>::const_iterator Itr_i;
	for (Itr_tt ttItr = m_triggerTowers->begin(); ttItr	!= m_triggerTowers->end(); ++ttItr) {
		
		// Set up variables to store the Trigger Tower values as they will be required in several places
		double /*ttEtEm  = 0.0, ttEtHad = 0.0,*/ttEta = 0.0, ttPhi = 0.0;

		// Get the values for the electromagnetic and hadronic Et in the towers
		//ttEtEm  = (*ttItr)->emEnergy();
		//ttEtHad = (*ttItr)->hadEnergy();

		// Get the values of eta and phi for the trigger towers
		ttEta = (*ttItr)->eta();
		ttPhi = (*ttItr)->phi();
		double tempTtPhi = (ttPhi > M_PI) ? ttPhi - (2 * M_PI) : ttPhi;

		// Calculate index of trigger tower for reference
		//TriggerTowerCollection::size_type idxVal = (ttItr - m_triggerTowers->begin());

		// Lets look at the database, disabled(dead) channels and for bad Calo Cells
		// We use 4 Database folders:
		// The Calibration folder
		//const coral::AttributeList* emDbCalib = m_tools->emDbAttributes(*ttItr,m_dbPpmChanCalib);
		// The dead channels folder (only has entries for dead channels - no entry = good channel)
		const coral::AttributeList* emDbDead = m_tools->emDbAttributes(*ttItr,
				m_dbPpmDeadChannels);

		int emDead(0), emBadCalo(0);
		if (emDbDead != 0) {
			emDead = m_tools->DeadChannel(emDbDead);
		}
		emBadCalo = m_tools->emBadCalo(*ttItr);

		if (emBadCalo) {
			m_h_TrigTower_emBadCalo->Fill(ttEta, tempTtPhi);
		}
		if (emDead) {
			m_h_TrigTower_emDeadChannel->Fill(ttEta, tempTtPhi);
		}
	}
}

//---------------------------------------------------------------
// Trigger Chain Analysis - find out which trigger chains have 
// passed so that we do not bias our efficiencies
//---------------------------------------------------------------
StatusCode EmEfficienciesMonTool::triggerChainAnalysis() {

	m_passed_L1_Jet_Trigger = false;
	m_passed_EF_SingleJet_Trigger = false;
	m_passed_EF_MultiJet_Trigger = false;
	m_passed_EF_egTau_Trigger = false;
	m_passed_EF_Trigger = false;

	// Get the list of all triggers but do this only once in the event loop
	if (m_configuredChains.size() == 0) {
		m_configuredChains = m_trigger->getListOfTriggers();
	}

	msg(MSG::DEBUG) << "Trigger Analysis: New Event" << endreq;

	//Fill a count in the histogram for every event
	for (std::vector<std::string>::const_iterator it =
			m_configuredChains.begin(); it != m_configuredChains.end(); ++it) {

		//Check that the trigger was passed
		if (m_trigger->isPassed(*it)) {

			//First ask if the event passed the L1 jet trigger items 
			if ((*it).find("L1_J") != std::string::npos) {
				m_passed_L1_Jet_Trigger = true;
				msg(MSG::DEBUG) << "Trigger Analysis: " << *it << endreq;
			}
			//First ask if the event passed the L1 em trigger items as a quick check 
			if ((*it).find("L1_EM") != std::string::npos) {
				msg(MSG::DEBUG) << "Trigger Analysis: " << *it << endreq;
			}

			//Then ask if the event passed the EF trigger chains and determine which ones they were			
			if ((*it).find("EF") != std::string::npos) {

				//Print out the trigger name
				m_passed_EF_Trigger = true;
				msg(MSG::DEBUG) << "Trigger Analysis: " << *it << endreq;

				//Find Event Filter chains corresponding to single jet triggers (keeping it simple)
				if (((*it).find("EF_j") != std::string::npos && (*it).find(
						"EF_je") == std::string::npos) || (*it).find("EF_fj")
						!= std::string::npos || ((*it).find("EF_L1J")
						!= std::string::npos && (*it).find("EMPTY")
						== std::string::npos) || (*it).find("EF_b")
						!= std::string::npos) {
					m_passed_EF_SingleJet_Trigger = true;
				}
				//Find Event Filter chains corresponding to multiple jet triggers (worry about it later)
				if ((*it).find("EF_2j") != std::string::npos || (*it).find(
						"EF_4j") != std::string::npos) {
					m_passed_EF_MultiJet_Trigger = true;
				}
				//Find Event Filter chains corresponding to photons
				if ((*it).find("EF_g") != std::string::npos || (*it).find(
						"EF_2g") != std::string::npos || (*it).find("EF_4g")
						!= std::string::npos) {
					m_passed_EF_egTau_Trigger = true;
				}
				//Find Event Filter chains corresponding to electrons
				if ((*it).find("EF_e") != std::string::npos || (*it).find(
						"EF_2e") != std::string::npos || (*it).find("EF_4e")
						!= std::string::npos) {
					m_passed_EF_egTau_Trigger = true;
				}
				//Find Event Filter chains corresponding to taus
				if ((*it).find("EF_tau") != std::string::npos || (*it).find(
						"EF_2tau") != std::string::npos
						|| (*it).find("EF_4tau") != std::string::npos) {
					m_passed_EF_egTau_Trigger = true;
				}
				//Find Event Filter chains corresponding to missing energy
				//Missing energy could come from electron so it may bias results 
				if ((*it).find("EF_xe") != std::string::npos || (*it).find(
						"EF_xs") != std::string::npos) {
					m_passed_EF_egTau_Trigger = true;
				}

				//Any event filter chains which do not pass any of these checks
				//would be treated as unbiased triggers
			}
		}
	}

	return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------------
// Load important containers
//------------------------------------------------------------------------------------
StatusCode EmEfficienciesMonTool::loadContainers() {
	StatusCode sc;

	m_eventInfo = 0;
	sc = evtStore()->retrieve(m_eventInfo);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load EventInfo" << endreq;
		return sc;
	}

	m_dbPpmDeadChannels = 0;
	sc = detStore()->retrieve(m_dbPpmDeadChannels, m_dbPpmDeadChannelsFolder);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load DB PPM Dead Channels" << endreq;
		return sc;
	}

	m_primaryVtx = 0;
	sc = evtStore()->retrieve(m_primaryVtx, m_primaryVertexLocation);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load Primary Vertices" << endreq;
		return sc;
	}

	if (m_numEvents == 1) {
		m_triggerTowers = 0;
		sc = evtStore()->retrieve(m_triggerTowers, m_triggerTowersLocation);
		if (sc.isFailure()) {
			msg(MSG::WARNING) << "Failed to load Trigger Towers" << endreq;
			return sc;
		}
	}

	m_lvl1RoIs = 0;
	sc = evtStore()->retrieve(m_lvl1RoIs, m_lvl1RoIsLocation);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load LVL1 RoIs" << endreq;
		return sc;
	}

	m_offlineElectrons = 0;
	sc = evtStore()->retrieve(m_offlineElectrons, m_offlineElectronsLocation);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load Offline Electrons" << endreq;
		return sc;
	}

	m_offlinePhotons = 0;
	sc = evtStore()->retrieve(m_offlinePhotons, m_offlinePhotonsLocation);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load Offline Photons" << endreq;
		return sc;
	}

	return sc;
}

// ********************************************************************
//
// NAME:     JetEfficienciesMonTool.cxx
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
#include "AthenaPoolUtilities/CondAttrListCollection.h"
#include "AthenaPoolUtilities/AthenaAttributeList.h"

#include "AthenaMonitoring/AthenaMonManager.h"

#include "TrigConfL1Data/L1DataDef.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigT1CaloEvent/TriggerTowerCollection.h"
#include "AnalysisTriggerEvent/LVL1_ROI.h"
#include "AnalysisTriggerEvent/EmTau_ROI.h"
#include "AnalysisTriggerEvent/Jet_ROI.h"
#include "TrigT1CaloToolInterfaces/IL1TriggerTowerTool.h"
#include "TrigT1CaloCalibToolInterfaces/IL1CaloLArTowerEnergy.h"
#include "EventInfo/EventInfo.h"
#include "EventInfo/EventID.h"
#include "VxVertex/VxContainer.h"
#include "VxVertex/VxTrackAtVertex.h"
#include "JetEvent/JetCollection.h"
#include "JetEvent/Jet.h"
#include "JetUtils/JetCaloQualityUtils.h"
#include "TileEvent/TileTTL1Cell.h"
#include "TrigT1CaloCalibConditions/L1CaloCoolChannelId.h"
#include "Identifier/Identifier.h"

#include "TrigT1CaloMonitoring/JetEfficienciesMonTool.h"
#include "TrigT1CaloMonitoringTools/TrigT1CaloLWHistogramTool.h"

/*---------------------------------------------------------*/
JetEfficienciesMonTool::JetEfficienciesMonTool(const std::string & type,
		const std::string & name, const IInterface* parent) 
		  : ManagedMonitorToolBase(type, name, parent),
			m_histTool("TrigT1CaloLWHistogramTool"),
			m_ttTool("LVL1::L1TriggerTowerTool/L1TriggerTowerTool"),
			m_larEnergy("LVL1::L1CaloLArTowerEnergy/L1CaloLArTowerEnergy"),
			m_trigger("Trig::TrigDecisionTool/TrigDecisionTool"),
			m_dbPpmDeadChannelsFolder("/TRIGGER/L1Calo/V1/Calibration/PpmDeadChannels"),
			m_triggerTowersLocation("TriggerTowers"),
			m_tileTTL1ContainerLocation("TileCellTTL1Container"),
			m_lvl1RoIsLocation("LVL1_ROI"),
			m_offlineJetsLocation("AntiKt4TopoEMJets"),
			m_primaryVertexLocation("VxPrimaryCandidate"), 
			m_eventInfo(0),
			m_dbPpmDeadChannels(0),
			m_triggerTowers(0), 
			m_lvl1RoIs(0),
			m_offlineJets(0),
			m_tileTTL1Container(0),
			m_primaryVertex(0),
			m_numEvents(0), 
			m_numOffJets(0), 
			m_numOffJetsInContainer(0), 
			m_numOffJetsPassCuts(0), 
			m_numOffJetsTriggered(0),
			m_numJetObjPassTrigger(0), 
			m_numJetObjTotal(0),
			m_passed_L1_EM_Trigger(false),
			m_passed_L1_Jet_Trigger(false),
			m_passed_EF_Trigger(false),
			m_passed_EF_SingleJet_Trigger(false),
			m_passed_EF_SingleEgamma_Trigger(false),
			m_passed_EF_MultiJet_Trigger(false),
			m_passed_EF_MultiEgamma_Trigger(false),
			m_passed_EF_Tau_Trigger(false), 
			m_passed_EF_MissingEnergy_Trigger(false),
			m_firstEvent(true),
			m_h_JetEmScale_Et(0),
			m_h_JetEmScale_Et_central(0),
			m_h_JetEmScale_Et_forward(0),
			m_h_JetEmScale_Et_triggered(0),
			m_h_JetEmScale_Et_triggered_central(0),
			m_h_JetEmScale_Et_triggered_forward(0),
			m_h_JetEmScale_Et_triggered_Eff(0),
			m_h_JetEmScale_Et_triggered_central_Eff(0),
			m_h_JetEmScale_Et_triggered_forward_Eff(0),
			m_h_JetEmScale_Eta_vs_Phi(0),
			m_h_JetEmScale_Eta_vs_Phi_triggered(0),
			m_h_JetEmScale_Eta_vs_Phi_triggered_Eff(0),
                        m_h_JetEmScale_50GeV_Eta_vs_Phi(0),
			m_h_JetEmScale_100GeV_Eta_vs_Phi(0),
			m_h_JetEmScale_200GeV_Eta_vs_Phi(0),
			m_h_TrigTower_jetBadCalo(0),
			m_h_TrigTower_jetDeadChannel(0),
			m_h_nPriVtx(0)

/*---------------------------------------------------------*/
{

	declareProperty("HistogramTool", m_histTool);
	declareProperty("TriggerTowerTool", m_ttTool);
	declareProperty("LArTowerEnergyTool", m_larEnergy);
	declareProperty("TrigDecisionTool", m_trigger);
	declareProperty("DeadChannelsFolder", m_dbPpmDeadChannelsFolder);
	declareProperty("TriggerTowersLocation", m_triggerTowersLocation);
	declareProperty("TileCellTTL1Location", m_tileTTL1ContainerLocation);
	declareProperty("RoIsLocation", m_lvl1RoIsLocation);
	declareProperty("OfflineJetsLocation",m_offlineJetsLocation);
	declareProperty("PrimaryVertexLocation", m_primaryVertexLocation);

	declareProperty("RootDirectory", m_rootDir = "L1Calo");
	declareProperty("UseTrigger", m_useTrigger = true);
	declareProperty("TriggerStrings", m_triggerStrings);

	// HSB - python cuts
	declareProperty("goodEMDeltaRMatch_Cut", m_goodEMDeltaRMatch_Cut = 0.5);
	declareProperty("goodHadDeltaRMatch_Cut", m_goodHadDeltaRMatch_Cut = 0.2);
	declareProperty("goodHadDeltaPhiMatch_Cut", m_goodHadDeltaPhiMatch_Cut = 0.3);
	declareProperty("UseEmThresholdsOnly", m_useEmThresholdsOnly = true);
	declareProperty("JetQualityLevel",m_jetQualityLevel = 30);
	declareProperty("NtracksAtPrimaryVertex",m_nTracksAtPrimaryVertex = 4);

	for (int i = 0; i < JET_ROI_BITS; ++i) {
	        m_h_JetEmScale_Et_J_item[i] = 0;
	        m_h_JetEmScale_Et_J_Eff_item[i] = 0;
		m_h_JetEmScale_50GeV_Eta_vs_Phi_J_item[i] = 0;
		m_h_JetEmScale_50GeV_Eta_vs_Phi_J_Eff_item[i] = 0;
		m_h_JetEmScale_100GeV_Eta_vs_Phi_J_item[i] = 0;
		m_h_JetEmScale_100GeV_Eta_vs_Phi_J_Eff_item[i] = 0;
		m_h_JetEmScale_200GeV_Eta_vs_Phi_J_item[i] = 0;
		m_h_JetEmScale_200GeV_Eta_vs_Phi_J_Eff_item[i] = 0;
	}
	for (int i = 0; i < FJET_ROI_BITS; ++i) {
	        m_h_JetEmScale_Et_FJ_J_item[i] = 0;
	        m_h_JetEmScale_Et_FJ_J_Eff_item[i] = 0;
		m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_item[i] = 0;
		m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_Eff_item[i] = 0;
	        m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_item[i] = 0;
	        m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_Eff_item[i] = 0;
		m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_item[i] = 0;
		m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_Eff_item[i] = 0;
		m_linkedHistos[i] = 0;
        }
}

/*---------------------------------------------------------*/
JetEfficienciesMonTool::~JetEfficienciesMonTool()
/*---------------------------------------------------------*/
{
}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*---------------------------------------------------------*/
StatusCode JetEfficienciesMonTool::initialize()
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

	sc = m_ttTool.retrieve();
	if (sc.isFailure()) {
		msg(MSG::ERROR) << "Cannot retrieve L1TriggerTowerTool" << endreq;
		return sc;
	}

	sc = m_larEnergy.retrieve();
	if (sc.isFailure()) {
		msg(MSG::ERROR) << "Cannot retrieve L1CaloLArTowerEnergy" << endreq;
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
StatusCode JetEfficienciesMonTool::bookHistograms(bool isNewEventsBlock,
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

	if (isNewRun) {

		std::string dir(m_rootDir + "/Reco/JetEfficiencies");

		MonGroup monJetDead(this, dir + "/DeadOrBadChannels", expert, run, "", "lowerLB");
		MonGroup monJetEmScaleNum(this, dir + "/JetEmScale_Et/numerator", expert, run);
		MonGroup monJetEmScaleDen(this, dir + "/JetEmScale_Et/denominator", expert, run);
		MonGroup monJetEmScaleEff(this, dir + "/JetEmScale_Et", expert, run, "", "perBinEffPerCent");
		MonGroup monJetEmScaleVtx(this, dir + "/JetEmScale_Et", expert, run);
		MonGroup monJetEmScale50GeVNum(this, dir + "/JetEmScale_50GeV_EtaVsPhi/numerator", expert, run);
		MonGroup monJetEmScale50GeVDen(this, dir + "/JetEmScale_50GeV_EtaVsPhi/denominator", expert, run);
		MonGroup monJetEmScale50GeVEff(this, dir + "/JetEmScale_50GeV_EtaVsPhi", expert, run, "", "perBinEffPerCent");
		MonGroup monJetEmScale100GeVNum(this, dir + "/JetEmScale_100GeV_EtaVsPhi/numerator", expert, run);
		MonGroup monJetEmScale100GeVDen(this, dir + "/JetEmScale_100GeV_EtaVsPhi/denominator", expert, run);
		MonGroup monJetEmScale100GeVEff(this, dir + "/JetEmScale_100GeV_EtaVsPhi", expert, run, "", "perBinEffPerCent");
		MonGroup monJetEmScale200GeVNum(this, dir + "/JetEmScale_200GeV_EtaVsPhi/numerator", expert, run);
		MonGroup monJetEmScale200GeVDen(this, dir + "/JetEmScale_200GeV_EtaVsPhi/denominator", expert, run);
		MonGroup monJetEmScale200GeVEff(this, dir + "/JetEmScale_200GeV_EtaVsPhi", expert, run, "", "perBinEffPerCent");

		//Set up JET thresholds arrays with threshold names
		std::string thrNum[JET_ROI_BITS] = { "0", "1", "2", "3", "4", "5", "6", "7" };
		TrigConf::L1DataDef def;
		std::vector < std::string > jetL1t;
		m_histTool->thresholdNames(def.jetType(), jetL1t);
		int size = jetL1t.size();
		for (int i = 0; i < JET_ROI_BITS; ++i) {
			if (i < size)
				jetL1t[i] = "L1_" + jetL1t[i];
			else
				jetL1t.push_back("L1_??");
		}
		std::vector < std::string > fjetL1t;
		m_histTool->thresholdNames(def.jfType(), fjetL1t);
		size = fjetL1t.size();
		for (int i = 0; i < FJET_ROI_BITS; ++i) {
			if (i < size)
				fjetL1t[i] = "L1_" + fjetL1t[i];
			else
				fjetL1t.push_back("L1_??");
		}

		//work out which L1_Jet histogram connects to which L1_ForwardJet histogram
		int count = 0;
		std::vector < std::string > fjetL1n(FJET_ROI_BITS);
		for(int i = 0; i < FJET_ROI_BITS; i++){
		        m_linkedHistos[i] = 0;
		        std::string strTemp = fjetL1t[i].substr(5);
		        for(int j = 0; j < JET_ROI_BITS; j++){
			        std::string strTemp2 = jetL1t[j].substr(4);
		                if(strTemp == strTemp2) {
		                        m_linkedHistos[i] = j;
					++count;
					break;
		                }
		        }
			fjetL1n[i] = fjetL1t[i] + " || " + jetL1t[m_linkedHistos[i]];
		}
		if (count != FJET_ROI_BITS) {
		        msg(MSG::WARNING) << "Jet ForwardJet mismatch" << endreq;
                }

		m_histTool->setMonGroup(&monJetDead);

		m_h_TrigTower_jetDeadChannel = m_histTool->bookPPMHadEtaVsPhi("TrigTower_jetDeadChannel","jet Trigger Towers with dead channels - #eta against #phi (E_{T} > 5 GeV)");

		m_h_TrigTower_jetBadCalo = m_histTool->bookPPMHadEtaVsPhi("TrigTower_jetBadCalo","jet Trigger Towers - Missing FEBs/Tile Quality - #eta against #phi (E_{T} > 5 GeV)");

		//Raw Jet Histograms

		m_histTool->setMonGroup(&monJetEmScaleDen);

		const int etbins = 200;
		const double etmax = 1000.;
		m_h_JetEmScale_Et = m_histTool->book1F("JetEmScale_Et","Raw Jet E_{T};E_{T} Jet [GeV];Jets", etbins, 0., etmax);
		m_h_JetEmScale_Et_forward = m_histTool->book1F("JetEmScale_Et_forward","Raw Jet E_{T} Forward;E_{T} Jet [GeV];Jets", etbins, 0., etmax);
		m_h_JetEmScale_Et_central = m_histTool->book1F("JetEmScale_Et_central","Raw Jet E_{T} Central;E_{T} Jet [GeV];Jets", etbins, 0., etmax);
		m_h_JetEmScale_Eta_vs_Phi = m_histTool->bookPPMHadEtaVsPhi("JetEmScale_Eta_vs_Phi","Raw Jet #eta v #phi");

		m_histTool->setMonGroup(&monJetEmScaleNum);

		m_h_JetEmScale_Et_triggered = m_histTool->book1F("JetEmScale_Et_triggered","Raw Jet E_{T} (Triggered);E_{T} Jet [GeV];Jets", etbins, 0., etmax);
		m_h_JetEmScale_Et_triggered_forward = m_histTool->book1F("JetEmScale_Et_triggered_forward","Raw Jet E_{T} Forward (Triggered);E_{T} Jet [GeV];Jets", etbins, 0., etmax);
		m_h_JetEmScale_Et_triggered_central = m_histTool->book1F("JetEmScale_Et_triggered_central","Raw Jet E_{T} Central (Triggered);E_{T} Jet [GeV];Jets", etbins, 0., etmax);
		m_h_JetEmScale_Eta_vs_Phi_triggered = m_histTool->bookPPMHadEtaVsPhi("JetEmScale_Eta_vs_Phi_triggered","Raw Jet #eta v #phi");

		std::string name;
		std::string title;
		for (int i = 0; i < JET_ROI_BITS; ++i) {
		        name = "JetEmScale_Et_J_item_" + thrNum[i];
		        title = "Raw E_{T} for Triggered Jets passing " + jetL1t[i] + ";E_{T} Jet [GeV];Jets";
		        m_h_JetEmScale_Et_J_item[i] =  m_histTool->book1F(name, title,  etbins, 0., etmax);
		        if (i < FJET_ROI_BITS) {
		                name = "JetEmScale_Et_FJ_J_item_" + thrNum[i];
			        title = "Raw E_{T} for Triggered Jets passing " + fjetL1n[i] + ";E_{T} Jet [GeV];Jets";
			        m_h_JetEmScale_Et_FJ_J_item[i] = m_histTool->book1F(name, title,  etbins, 0., etmax);
		        }
		}

		m_histTool->setMonGroup(&monJetEmScaleEff);

		m_h_JetEmScale_Et_triggered_Eff = m_histTool->book1F("JetEmScale_Et_triggered_Eff","Raw Jet E_{T} (Triggered) Efficiency;E_{T} Jet [GeV];Efficiency %", etbins, 0., etmax);
		m_h_JetEmScale_Et_triggered_forward_Eff = m_histTool->book1F("JetEmScale_Et_triggered_forward_Eff","Raw Jet E_{T} Forward (Triggered) Efficiency;E_{T} Jet [GeV];Efficiency %", etbins, 0., etmax);
		m_h_JetEmScale_Et_triggered_central_Eff = m_histTool->book1F("JetEmScale_Et_triggered_central_Eff","Raw Jet E_{T} Central (Triggered) Efficiency;E_{T} Jet [GeV];Efficiency %", etbins, 0., etmax);
		m_h_JetEmScale_Eta_vs_Phi_triggered_Eff = m_histTool->bookPPMHadEtaVsPhi("JetEmScale_Eta_vs_Phi_triggered_Eff","Raw Jet #eta v #phi Efficiency (%)");

		for (int i = 0; i < JET_ROI_BITS; ++i) {
		        name = "JetEmScale_Et_J_Eff_item_" + thrNum[i];
		        title = "Raw E_{T} for Triggered Jets passing " + jetL1t[i] + " Efficiency;E_{T} Jet [GeV];Efficiency %";
		        m_h_JetEmScale_Et_J_Eff_item[i] =  m_histTool->book1F(name, title,  etbins, 0., etmax);
		        if (i < FJET_ROI_BITS) {
		                name = "JetEmScale_Et_FJ_J_Eff_item_" + thrNum[i];
			        title = "Raw E_{T} for Triggered Jets passing " + fjetL1n[i] + " Efficiency;E_{T} Jet [GeV]; Efficiency %";
			        m_h_JetEmScale_Et_FJ_J_Eff_item[i] = m_histTool->book1F(name, title,  etbins, 0., etmax);
		        }
		}

		m_histTool->setMonGroup(&monJetEmScaleVtx);

		m_h_nPriVtx = m_histTool->book1F("nPriVtx","Primary Vertex Multiplicity",30,0,30);

		m_histTool->setMonGroup(&monJetEmScale50GeVNum);

		//Et Raw Jet greater > 50 GeV
		for (int i = 0; i < JET_ROI_BITS; ++i) {
		        name = "JetEmScale_50GeV_Eta_vs_Phi_J_item_" + thrNum[i];
		        title = "Raw Jet #eta v #phi (triggered on " + jetL1t[i] + " with E_{T} > 50 GeV)";
		        m_h_JetEmScale_50GeV_Eta_vs_Phi_J_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
		        if (i < FJET_ROI_BITS) {
		                name = "JetEmScale_50GeV_Eta_vs_Phi_FJ_J_item_" + thrNum[i];
			        title = "Raw Jet #eta v #phi (triggered on " + fjetL1n[i] + " with E_{T} > 50 GeV)";
			        m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
                        }
                }

		m_histTool->setMonGroup(&monJetEmScale50GeVDen);

		m_h_JetEmScale_50GeV_Eta_vs_Phi = m_histTool->bookPPMHadEtaVsPhi("JetEmScale_50GeV_Eta_vs_Phi","Raw Jet #eta v #phi (E_{T} > 50 GeV)");

		m_histTool->setMonGroup(&monJetEmScale50GeVEff);

		for (int i = 0; i < JET_ROI_BITS; ++i) {
			name = "JetEmScale_50GeV_Eta_vs_Phi_J_Eff_item_" + thrNum[i];
			title = "Raw Jet #eta v #phi (triggered on " + jetL1t[i] + " with E_{T} > 50 GeV) Efficiency (%)";
			m_h_JetEmScale_50GeV_Eta_vs_Phi_J_Eff_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
			if (i < FJET_ROI_BITS) {
			        name = "JetEmScale_50GeV_Eta_vs_Phi_FJ_J_Eff_item_" + thrNum[i];
				title = "Raw Jet #eta v #phi (triggered on " + fjetL1n[i] + " with E_{T} > 50 GeV) Efficiency (%)";
				m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_Eff_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
		        }
		}

		m_histTool->setMonGroup(&monJetEmScale100GeVNum);

		//Et Raw Jet greater > 100 GeV
		for (int i = 0; i < JET_ROI_BITS; ++i) {
		        name = "JetEmScale_100GeV_Eta_vs_Phi_J_item_" + thrNum[i];
		        title = "Raw Jet #eta v #phi (triggered on " + jetL1t[i] + " with E_{T} > 100 GeV)";
		        m_h_JetEmScale_100GeV_Eta_vs_Phi_J_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
		        if (i < FJET_ROI_BITS) {
		                name = "JetEmScale_100GeV_Eta_vs_Phi_FJ_J_item_" + thrNum[i];
			        title = "Raw Jet #eta v #phi (triggered on " + fjetL1n[i] + " with E_{T} > 100 GeV)";
			        m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
                        }
                }

		m_histTool->setMonGroup(&monJetEmScale100GeVDen);

		m_h_JetEmScale_100GeV_Eta_vs_Phi = m_histTool->bookPPMHadEtaVsPhi("JetEmScale_100GeV_Eta_vs_Phi","Raw Jet #eta v #phi (E_{T} > 100 GeV)");

		m_histTool->setMonGroup(&monJetEmScale100GeVEff);

		for (int i = 0; i < JET_ROI_BITS; ++i) {
			name = "JetEmScale_100GeV_Eta_vs_Phi_J_Eff_item_" + thrNum[i];
			title = "Raw Jet #eta v #phi (triggered on " + jetL1t[i] + " with E_{T} > 100 GeV) Efficiency (%)";
			m_h_JetEmScale_100GeV_Eta_vs_Phi_J_Eff_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
			if (i < FJET_ROI_BITS) {
			        name = "JetEmScale_100GeV_Eta_vs_Phi_FJ_J_Eff_item_" + thrNum[i];
				title = "Raw Jet #eta v #phi (triggered on " + fjetL1n[i] + " with E_{T} > 100 GeV) Efficiency (%)";
				m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_Eff_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
		        }
		}

		m_histTool->setMonGroup(&monJetEmScale200GeVNum);

		//Et Raw Jet greater > 200 GeV
		for (int i = 0; i < JET_ROI_BITS; ++i) {
		        name = "JetEmScale_200GeV_Eta_vs_Phi_J_item_" + thrNum[i];
		        title = "Raw Jet #eta v #phi (triggered on " + jetL1t[i] + " with E_{T} > 200 GeV)";
		        m_h_JetEmScale_200GeV_Eta_vs_Phi_J_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
		        if (i < FJET_ROI_BITS) {
		                name = "JetEmScale_200GeV_Eta_vs_Phi_FJ_J_item_" + thrNum[i];
			        title = "Raw Jet #eta v #phi (triggered on " + fjetL1n[i] + " with E_{T} > 200 GeV)";
			        m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
                        }
                }

		m_histTool->setMonGroup(&monJetEmScale200GeVDen);

		m_h_JetEmScale_200GeV_Eta_vs_Phi = m_histTool->bookPPMHadEtaVsPhi("JetEmScale_200GeV_Eta_vs_Phi","Raw Jet #eta v #phi (E_{T} > 200 GeV)");

		m_histTool->setMonGroup(&monJetEmScale200GeVEff);

		for (int i = 0; i < JET_ROI_BITS; ++i) {
			name = "JetEmScale_200GeV_Eta_vs_Phi_J_Eff_item_" + thrNum[i];
			title = "Raw Jet #eta v #phi (triggered on " + jetL1t[i] + " with E_{T} > 200 GeV) Efficiency (%)";
			m_h_JetEmScale_200GeV_Eta_vs_Phi_J_Eff_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
			if (i < FJET_ROI_BITS) {
			        name = "JetEmScale_200GeV_Eta_vs_Phi_FJ_J_Eff_item_" + thrNum[i];
				title = "Raw Jet #eta v #phi (triggered on " + fjetL1n[i] + " with E_{T} > 200 GeV) Efficiency (%)";
				m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_Eff_item[i] = m_histTool->bookPPMHadEtaVsPhi(name, title);
		        }
		}

		m_histTool->unsetMonGroup();

		// HSB - counters 
		m_numEvents = 0;
		m_numOffJets = 0;
		m_numOffJetsInContainer = 0;
		m_numOffJetsPassCuts = 0;
		m_numOffJetsTriggered = 0;
		m_numJetObjPassTrigger = 0;
		m_numJetObjTotal = 0;
	}

	msg(MSG::DEBUG) << "Leaving bookHistograms" << endreq;

	return StatusCode::SUCCESS;
}

/*---------------------------------------------------------*/
StatusCode JetEfficienciesMonTool::fillHistograms()
/*---------------------------------------------------------*/
{
	const bool debug = msgLvl(MSG::DEBUG);
	if (debug) msg(MSG::DEBUG) << "fillHistograms entered" << endreq;

	// On first event plot disabled channels and bad calo
	if (m_firstEvent) {
		m_firstEvent = false;
		StatusCode sc = this->triggerTowerAnalysis();
		if (sc.isFailure()) {
		        msg(MSG::WARNING) << "Problem analysing Trigger Towers" << endreq;
			return sc;
                }
	}

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

		// Look at vertex requirements
		unsigned int nPriVtx = this->nPrimaryVertex();
		if (nPriVtx < 1) {
			if (debug) msg(MSG::DEBUG) << "Event " << m_eventInfo->event_ID()->event_number() << " fails vertex requirements " << endreq;
			return StatusCode::SUCCESS;
		}

		m_h_nPriVtx->Fill(nPriVtx);

		m_numOffJetsInContainer = 0;
		m_numOffJetsPassCuts = 0;
		m_numOffJetsTriggered = 0;
		sc = this->analyseOfflineJets();
		if (sc.isFailure()) {
			msg(MSG::WARNING) << "analyseJets Failed " << endreq;
			return sc;
		}
		if (debug) msg(MSG::DEBUG) << "Number of Offline Jets = "<< m_numOffJetsInContainer << " Passing Cuts = "<< m_numOffJetsPassCuts << " Triggered = "<< m_numOffJetsTriggered << endreq;

	}

	if (debug)
		msg(MSG::DEBUG) << "Leaving fillHistograms" << endreq;

	return StatusCode::SUCCESS;

}

/*---------------------------------------------------------*/
StatusCode JetEfficienciesMonTool::procHistograms(bool isEndOfEventsBlock,
		bool isEndOfLumiBlock, bool isEndOfRun)
/*---------------------------------------------------------*/
{
	msg(MSG::DEBUG) << "procHistograms entered" << endreq;

	if (isEndOfEventsBlock || isEndOfLumiBlock) {
	}

	if (isEndOfRun) {
		msg(MSG::DEBUG) << "Number of offline jets = " << m_numOffJets << endreq;
		msg(MSG::DEBUG) << "Number of events = " << m_numEvents << endreq;

		m_histTool->efficienciesForMerge(m_h_JetEmScale_Et,
				m_h_JetEmScale_Et_triggered,
				m_h_JetEmScale_Et_triggered_Eff);
		m_histTool->efficienciesForMerge(m_h_JetEmScale_Et_forward,
		                m_h_JetEmScale_Et_triggered_forward,
				m_h_JetEmScale_Et_triggered_forward_Eff);
		m_histTool->efficienciesForMerge(m_h_JetEmScale_Et_central,
		                m_h_JetEmScale_Et_triggered_central,
				m_h_JetEmScale_Et_triggered_central_Eff);
		m_histTool->efficienciesForMerge(m_h_JetEmScale_Eta_vs_Phi,
		                m_h_JetEmScale_Eta_vs_Phi_triggered,
				m_h_JetEmScale_Eta_vs_Phi_triggered_Eff);

		for (int i = 0; i < JET_ROI_BITS; ++i) {
			m_histTool->efficienciesForMerge(m_h_JetEmScale_Et,
			                m_h_JetEmScale_Et_J_item[i],
                                        m_h_JetEmScale_Et_J_Eff_item[i]);
			m_histTool->efficienciesForMerge(m_h_JetEmScale_50GeV_Eta_vs_Phi,
			                m_h_JetEmScale_50GeV_Eta_vs_Phi_J_item[i],
                                        m_h_JetEmScale_50GeV_Eta_vs_Phi_J_Eff_item[i]);
			m_histTool->efficienciesForMerge(m_h_JetEmScale_100GeV_Eta_vs_Phi,
			                m_h_JetEmScale_100GeV_Eta_vs_Phi_J_item[i],
                                        m_h_JetEmScale_100GeV_Eta_vs_Phi_J_Eff_item[i]);
			m_histTool->efficienciesForMerge(m_h_JetEmScale_200GeV_Eta_vs_Phi,
			                m_h_JetEmScale_200GeV_Eta_vs_Phi_J_item[i],
                                        m_h_JetEmScale_200GeV_Eta_vs_Phi_J_Eff_item[i]);
			if (i >= FJET_ROI_BITS) continue;
			// Copy the Jet Trigger information into the Forward Jet histograms (as it should be the logical OR of these two sets)
			m_h_JetEmScale_Et_FJ_J_item[i]->getROOTHistBase()->Add(m_h_JetEmScale_Et_J_item[m_linkedHistos[i]]->getROOTHistBase());
			m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_item[i]->getROOTHistBase()->Add(m_h_JetEmScale_50GeV_Eta_vs_Phi_J_item[m_linkedHistos[i]]->getROOTHistBase());
			m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_item[i]->getROOTHistBase()->Add(m_h_JetEmScale_100GeV_Eta_vs_Phi_J_item[m_linkedHistos[i]]->getROOTHistBase());
			m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_item[i]->getROOTHistBase()->Add(m_h_JetEmScale_200GeV_Eta_vs_Phi_J_item[m_linkedHistos[i]]->getROOTHistBase());
			m_histTool->efficienciesForMerge(m_h_JetEmScale_Et,
			                m_h_JetEmScale_Et_FJ_J_item[i],
                                        m_h_JetEmScale_Et_FJ_J_Eff_item[i]);
			m_histTool->efficienciesForMerge(m_h_JetEmScale_50GeV_Eta_vs_Phi,
			                m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_item[i],
                                        m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_Eff_item[i]);
			m_histTool->efficienciesForMerge(m_h_JetEmScale_100GeV_Eta_vs_Phi,
			                m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_item[i],
                                        m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_Eff_item[i]);
			m_histTool->efficienciesForMerge(m_h_JetEmScale_200GeV_Eta_vs_Phi,
			                m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_item[i],
                                        m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_Eff_item[i]);
		}

	}

	return StatusCode::SUCCESS;
}

//------------------------------------------------------------------
 // The check to see if an object is triggered w.r.t. a RoI
 // It can be done in two ways so allow it to handle either one
//------------------------------------------------------------------
bool JetEfficienciesMonTool::deltaMatch(double etaJet, double dR, double dPhi, bool isForward) {
	
	// If within Tile & HEC then use delta R matching
	// if within FCAL then use delta Phi matching

	double fabsEtaJet = fabs(etaJet); 
	
	if(fabsEtaJet < 2.899) { //within Tile&Hec
		if(dR < m_goodHadDeltaRMatch_Cut) {
			return true;
		} else {
			return false;
		}
	} else if(fabsEtaJet >= 2.899 && fabsEtaJet < 3.2) { //end of HEC so use combined matching between central & forward regions
		if((dR < m_goodHadDeltaRMatch_Cut && !isForward) || (dPhi < m_goodHadDeltaPhiMatch_Cut && isForward)) {
			return true;
		} else {
			return false;
		}
	} else { //within FCAL (3.2 to 4.9)
		if(dPhi < m_goodHadDeltaPhiMatch_Cut) {
			return true;
		} else {
			return false;
		}
	}	
}

//------------------------------------------------------------------
 //Correct the value of deltaPhi so that it falls in +-PI range
//------------------------------------------------------------------
double JetEfficienciesMonTool::correctDeltaPhi(double dPhi) {
	
	if(fabs(dPhi) > M_PI) {
		dPhi = (dPhi > 0) ? dPhi - 2*M_PI : dPhi + 2*M_PI;
	}
	
	return dPhi;
}

//------------------------------------------------------------------
 //Calculate delta R quickly between two objects
//------------------------------------------------------------------
double JetEfficienciesMonTool::calcDeltaR(double eta1, double phi1, double eta2, double phi2) {
	
	double dEta = eta1-eta2;
	double dPhi = correctDeltaPhi(phi1-phi2);
	
	double dR = sqrt((dEta*dEta) + (dPhi*dPhi));
	return dR;
}

//------------------------------------------------------------------
 //Check for electron that it is of the right jet quality type as required from jobOptions
//------------------------------------------------------------------
bool JetEfficienciesMonTool::correctJetQuality(const Jet* jet) {
	
	bool correctType = false;
	switch (m_jetQualityLevel) {
  	    case 0: //"None"
        	    correctType = true; break;
            case 10: //"Jet Loose" 
        	    correctType = (JetCaloQualityUtils::isGood(jet,false)!=0) ? true : false; break;
	    case 20: //"Jet Medium"
        	    correctType = (JetCaloQualityUtils::isGoodMedium(jet,false)!=0) ? true : false; break;
	    case 30: //"Jet Tight"
        	    correctType = (JetCaloQualityUtils::isGoodTight(jet,false)!=0) ? true : false; break;
  	    default:
        	    correctType = false; break;
	}
	return correctType;
}

/**********************************************/
//Analysis code for offline reconstructed jets
//Compares jets with Jet RoIs
/**********************************************/
StatusCode JetEfficienciesMonTool::analyseOfflineJets() {
	    	
	//Access all of the offline reconstructed jet candidates in the event
	typedef JetCollection::const_iterator Itr_jets;
	m_numOffJetsInContainer = m_offlineJets->size();
	m_numJetObjTotal += m_numOffJetsInContainer;
	
	// Create variables for jet properties
	double EtOJ = 0.0, etaOJ = 0.0, fabsEtaOJ = 0.0, phiOJ = 0.0, phiOJ_L1C = 0.0;
	// Create variable to determine if selecting the right type of jets based on criteria in jO
	bool correctType = false;

	bool roiValuesFilled = false;

	//Cycle through all of the offline reconstructed jets
	for (Itr_jets jetItr = m_offlineJets->begin(); jetItr != m_offlineJets->end(); ++jetItr) {
		
		//JetCollection::size_type jIdxVal = (jetItr - m_offlineJets->begin());
		
		//Keep track of eta, phi and Et as these will be used often
		EtOJ = (*jetItr)->et() / GeV;
		etaOJ = (*jetItr)->eta();
		fabsEtaOJ = fabs(etaOJ); 
		phiOJ = (*jetItr)->phi();
		phiOJ_L1C = l1caloPhi(phiOJ);
				
		//Alternate triggers - the ones that will bias the result
		bool altTrigger = m_passed_EF_Tau_Trigger || m_passed_EF_MissingEnergy_Trigger || m_passed_EF_SingleJet_Trigger || m_passed_EF_MultiJet_Trigger;

		//Check that the trigger selection is not biased
		bool unbiasedTrigger = false;
		if (m_passed_EF_Trigger == true) {
			if (m_passed_EF_SingleEgamma_Trigger == true) {
				unbiasedTrigger = isolatedJetObjectEF(phiOJ, etaOJ);
			} else {
				if(altTrigger == true) {
					return StatusCode::SUCCESS;
				} else {
					unbiasedTrigger = true;
				}
			}
		} else if(m_passed_L1_EM_Trigger) {
			unbiasedTrigger = isolatedJetObjectL1(phiOJ, etaOJ);			
		} else {
			return StatusCode::SUCCESS;
		}
		
		//If passed the trigger conditions then proceed to start analysis
		if (unbiasedTrigger == true) {
			
			correctType = false;
			m_numJetObjPassTrigger++;
			
			//Check that the good jet (not bad, ugly) matches the type selected (if any was chosen in jobOptions file)
			correctType = correctJetQuality((*jetItr));

			//Have the correct type of jet so do some analysis with it   
			if (correctType) {
				//Update counters
				m_numOffJets++;
				m_numOffJetsPassCuts++;

				if (EtOJ > 0.0) {
					m_h_JetEmScale_Et->Fill(EtOJ);
					if(fabsEtaOJ >= 2.899) { //treat 2.899 as new forward range due to matching in this region, not 3.2
						m_h_JetEmScale_Et_forward->Fill(EtOJ);
					} else {
						m_h_JetEmScale_Et_central->Fill(EtOJ);
					}
					
					m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_Eta_vs_Phi, etaOJ, phiOJ_L1C);

					if (EtOJ > 50) { 
						m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_50GeV_Eta_vs_Phi, etaOJ, phiOJ_L1C);
					}
					if (EtOJ > 100) { 
						m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_100GeV_Eta_vs_Phi, etaOJ, phiOJ_L1C);
					}
					if (EtOJ > 200) { 
						m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_200GeV_Eta_vs_Phi, etaOJ, phiOJ_L1C);
					}
				}

				//Set up useful numbers to keep track of RoI information
				double etaROI = 0.0, phiROI = 0.0, phiROI_L1C = 0.0;
				double dEta = 0.0, dPhi = 0.0, dR = 1000, temp_dR = 0.0;
				double bestDeltaPhi = 0.0;
				uint32_t ROIWord = 0;
				bool isForward = false, bestIsForward = false;
				
				//Access the Jet RoIs
				typedef std::vector<Jet_ROI>::const_iterator Itr_jetroi;

				//Iterate over all of the Jet RoIs
				for (Itr_jetroi roiItr = m_jetROIs.begin(); roiItr != m_jetROIs.end(); ++roiItr) {
											
					//std::vector<Jet_ROI>::size_type idxVal = (roiItr - m_jetROIs.begin());
					
					//Get useful values for the Jet RoI
					etaROI = (*roiItr).getEta();
					phiROI = (*roiItr).getPhi();
					phiROI_L1C = l1caloPhi(phiROI);
					if(fabs(etaROI) >= 3.2) { isForward = true; }
					else { isForward = false; }
					
					//Calculate the difference in eta and phi between the jet and RoI
					dEta = etaOJ - etaROI;
					dPhi = correctDeltaPhi(phiOJ_L1C - phiROI_L1C);
					
					//Calculate delta R
					temp_dR = sqrt(dEta * dEta + dPhi * dPhi);
					
				
					//Check if the new delta R is smaller than any previous delta R value.
					//In that case, keep track of the new RoI values
					if (temp_dR < dR) {
						
						//RoI information
						ROIWord = (*roiItr).getROIWord();

						bestDeltaPhi = dPhi;
						dR = temp_dR;
						
						bestIsForward = isForward;
					}						
				}

				roiValuesFilled = true;

				//Check to see if there was an RoI to match with an jet
				if (dR != 1000) {
					
					//Check if jet and RoI matched to a very good level (less than cut) - default now 0.2
					if(deltaMatch(etaOJ, dR, bestDeltaPhi, bestIsForward)) {
						
						m_numOffJetsTriggered++;
						
						bool matchToJet[JET_ROI_BITS];
						for (int k = 0; k < JET_ROI_BITS; ++k) { matchToJet[k] = false; }
						
						if (EtOJ > 0.0) {
							
							m_h_JetEmScale_Et_triggered->Fill(EtOJ);
							if(fabsEtaOJ >= 2.899) {
								m_h_JetEmScale_Et_triggered_forward->Fill(EtOJ);							
							} else {
								m_h_JetEmScale_Et_triggered_central->Fill(EtOJ);								
							}
							
							m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_Eta_vs_Phi_triggered, etaOJ, phiOJ_L1C);

							//Look at each bit in the RoI word (using bitshift) to see which thresholds were
							//passed and which ones were not  
							
							// First look at jets only up to 2.899 in mod eta
							if(fabsEtaOJ < 2.899) {
								for (int k = 0; k < JET_ROI_BITS; ++k) {
									matchToJet[k] = false;
									if ((ROIWord >> k) & 1) {
										
										m_h_JetEmScale_Et_J_item[k]->Fill(EtOJ);
										matchToJet[k] = true;

										if (EtOJ > 50) {
											m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_50GeV_Eta_vs_Phi_J_item[k], etaOJ, phiOJ_L1C);
										}
										if (EtOJ > 100) {
											m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_100GeV_Eta_vs_Phi_J_item[k], etaOJ, phiOJ_L1C);
										}
										if (EtOJ > 200) {
											m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_200GeV_Eta_vs_Phi_J_item[k], etaOJ, phiOJ_L1C);
										}
										 
									}
								}
							} else {
								
								if(bestIsForward) { //implies in very forward region (greater than 3.2)
									
									for (int k = 0; k < FJET_ROI_BITS; ++k) {
										if ((ROIWord >> (k+8)) & 1) {
											if(!matchToJet[m_linkedHistos[k]]) {
												m_h_JetEmScale_Et_FJ_J_item[k]->Fill(EtOJ);
											
												if (EtOJ > 50) {
											                m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_item[k], etaOJ, phiOJ_L1C);
												}
												if (EtOJ > 100) {
											                m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_item[k], etaOJ, phiOJ_L1C);
												}
												if (EtOJ > 200) {
											                m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_item[k], etaOJ, phiOJ_L1C);
												}

											}
										}
									}
									
								} else { //between 2.899 and 3.2
									
									for (int k = 0; k < JET_ROI_BITS; ++k) {
										if ((ROIWord >> k) & 1) {
											if(!matchToJet[k]) {
												
												int fj_j_link = -1;											
												for(Int_t l = 0; l < FJET_ROI_BITS; ++l) {
													if(m_linkedHistos[l] == k) { fj_j_link = l;	}
												}
												
												if(fj_j_link != -1) {
													m_h_JetEmScale_Et_FJ_J_item[fj_j_link]->Fill(EtOJ);
																
													if (EtOJ > 50) {
											                        m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_item[fj_j_link], etaOJ, phiOJ_L1C);
													}
													if (EtOJ > 100) {
											                        m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_item[fj_j_link], etaOJ, phiOJ_L1C);
													}
													if (EtOJ > 200) {
											                        m_histTool->fillPPMHadEtaVsPhi(m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_item[fj_j_link], etaOJ, phiOJ_L1C);
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
			}
		}
	}		
		
	return StatusCode::SUCCESS;
}


//------------------------------------------------------------------
 //Ask if object is has no EmTau RoIs nearby at L1
 //Do this by assuming that the highest ET EmTau RoI is the one causing
 //the event to trigger at HLT level
//------------------------------------------------------------------
bool JetEfficienciesMonTool::isolatedJetObjectL1(double phi, double eta) {
	
    bool isolated = true;
	double dREm = 10.0, dR_Max = 10.0; 
	double ET_Max = -10.0;
	double etaROI = 0.0, phiROI = 0.0, ET_ROI = 0.0;

	typedef std::vector<EmTau_ROI>::const_iterator Itr_emTauRoi;

	//Cycle over the rois, get their properties and determine the distance from the object
	for (Itr_emTauRoi roiItr = m_emTauROIs.begin(); roiItr != m_emTauROIs.end(); ++roiItr) {
		bool emThresholdPassed = false;
		if (m_useEmThresholdsOnly == true) {
			std::vector<std::string> thrPassed = (*roiItr).getThresholdNames();
			typedef std::vector<std::string>::iterator Itr_s;				
			for (Itr_s i = thrPassed.begin(); i != thrPassed.end(); ++i) {					
				if ((*i).find("EM") != std::string::npos) {
					emThresholdPassed = true;						
				}
			}
		}
		
		if(emThresholdPassed == true) {			
			etaROI = (*roiItr).getEta();
			phiROI = (*roiItr).getPhi();
			ET_ROI = (*roiItr).getEMClus(); //(*roiItr).getET8x8();			
			dREm = calcDeltaR(eta, phi, etaROI, phiROI);
			
			//If this energy exceeds the current record then store the details
			if(ET_ROI > ET_Max) {
				ET_Max = ET_ROI;
				dR_Max = dREm;
			}
		}
	}
	
	// Check that the object is far away enough from highest ET jet RoI
	if (dR_Max > m_goodEMDeltaRMatch_Cut) { 
		isolated = true;
	} else {
		isolated = false;
	}

	return isolated;
}


//------------------------------------------------------------------
 //Ask if object is has no jets or jet RoIs nearby at EF
//------------------------------------------------------------------
bool JetEfficienciesMonTool::isolatedJetObjectEF(double phi, double eta) {
	
	//At present just call L1 function, change if there are technicalities
    return isolatedJetObjectL1(phi, eta);
}

//---------------------------------------------------------------
 // Trigger Tower Analysis
//---------------------------------------------------------------
StatusCode JetEfficienciesMonTool::triggerTowerAnalysis() {

  m_dbPpmDeadChannels = 0;
  StatusCode sc = detStore()->retrieve(m_dbPpmDeadChannels, m_dbPpmDeadChannelsFolder);
  if (sc.isFailure()) {
    msg(MSG::WARNING) << "Failed to load DB PPM Dead Channels" << endreq;
    return sc;
  }

  m_triggerTowers = 0;
  sc = evtStore()->retrieve(m_triggerTowers, m_triggerTowersLocation);
  if (sc.isFailure()) {
    msg(MSG::WARNING) << "Failed to load Trigger Towers" << endreq;
    return sc;
  }

  m_tileTTL1Container = 0;
  sc = evtStore()->retrieve(m_tileTTL1Container,m_tileTTL1ContainerLocation);
  if (sc.isFailure()) {
    msg(MSG::WARNING) << "Failed to load TileTTL1" << endreq;
    return sc;
  }
  this->preCacheTTL1Cell(m_tileTTL1Container);
       
  typedef TriggerTowerCollection::const_iterator Itr_tt;
  for(Itr_tt ttItr=m_triggerTowers->begin();ttItr!=m_triggerTowers->end();++ttItr) {

    // Get the values of eta and phi for the trigger towers
    double ttEta = (*ttItr)->eta();
    double ttPhi = (*ttItr)->phi();
    const L1CaloCoolChannelId emCoolId(m_ttTool->channelID(ttEta, ttPhi, 0));
    const L1CaloCoolChannelId hadCoolId(m_ttTool->channelID(ttEta, ttPhi, 1));
    const Identifier emIdent(m_ttTool->identifier(ttEta, ttPhi, 0));
    const Identifier hadIdent(m_ttTool->identifier(ttEta, ttPhi, 1));

    // The Disabled Towers DB folder is used for 2011 data (& 2012 hopefully)
    unsigned int emDisabled(0), hadDisabled(0);
	// The dead channels folder (only has entries for dead channels - no entry = good channel)
	CondAttrListCollection::const_iterator itr = m_dbPpmDeadChannels->chanAttrListPair(emCoolId.id());
        if (itr != m_dbPpmDeadChannels->end()) {
                const AthenaAttributeList& attrList(itr->second);
                emDisabled = attrList["ErrorCode"].data<unsigned int>();
        }
	itr = m_dbPpmDeadChannels->chanAttrListPair(hadCoolId.id());
        if (itr != m_dbPpmDeadChannels->end()) {
                const AthenaAttributeList& attrList(itr->second);
                hadDisabled = attrList["ErrorCode"].data<unsigned int>();
        }
    
    int emBadCalo(0),hadBadCalo(0);
	emBadCalo = m_larEnergy->hasMissingFEB(emIdent);
	if (fabs(ttEta) > 1.5) {
	        hadBadCalo = m_larEnergy->hasMissingFEB(hadIdent);
        } else {
		IdTTL1CellMapType::const_iterator ttL1Cell(m_idTTL1CellMap.find(hadIdent));
		IdTTL1CellMapType::const_iterator ttL1Cell_E(m_idTTL1CellMap.end());
		if (ttL1Cell != ttL1Cell_E) {
		    hadBadCalo = (ttL1Cell->second)->qualTower();
		}
        }
    
    if(emBadCalo || hadBadCalo) m_histTool->fillPPMHadEtaVsPhi(m_h_TrigTower_jetBadCalo,ttEta,ttPhi,(double)(emBadCalo+hadBadCalo));
    if(emDisabled || hadDisabled) m_histTool->fillPPMHadEtaVsPhi(m_h_TrigTower_jetDeadChannel,ttEta,ttPhi,(double)(emDisabled+hadDisabled));
  }
  m_idTTL1CellMap.clear();
  
  return StatusCode::SUCCESS;
}

//---------------------------------------------------------------
 // Trigger Chain Analysis - find out which trigger chains have 
 // passed so that we do not bias our efficiencies
//---------------------------------------------------------------
StatusCode JetEfficienciesMonTool::triggerChainAnalysis() {

	// L1 Triggers
	m_passed_L1_Jet_Trigger = false;
	m_passed_L1_EM_Trigger = false;
	// EF Triggers
	m_passed_EF_Trigger = false;
	m_passed_EF_SingleJet_Trigger = false;
	m_passed_EF_MultiJet_Trigger = false;
	m_passed_EF_SingleEgamma_Trigger = false;
	m_passed_EF_MultiEgamma_Trigger = false;
	m_passed_EF_Tau_Trigger = false;
	m_passed_EF_MissingEnergy_Trigger = false;
	
	// Get the list of all triggers but do this only once in the event loop
	if (m_configuredChains.size() == 0) {
		m_configuredChains = m_trigger->getListOfTriggers();
	}

	//std::cout << "Trigger Analysis: New Event ---------------------------------" << std::endl;

	//Fill a count in the histogram for every event
	for (std::vector<std::string>::const_iterator it = m_configuredChains.begin(); it != m_configuredChains.end(); ++it) {
		
		//Check that the trigger was passed
		if (m_trigger->isPassed(*it)) {
			
			//First ask if the event passed the L1 jet trigger items 
			if ((*it).find("L1_J") != std::string::npos && (*it).find("L1_JE") == std::string::npos) {
				m_passed_L1_Jet_Trigger = true;
				//std::cout << "Trigger Analysis: " << *it << std::endl;
			}
			//First ask if the event passed the L1 em trigger items as a quick check 
			if ((*it).find("L1_EM") != std::string::npos) {
				m_passed_L1_EM_Trigger = true;
				//std::cout << "Trigger Analysis: " << *it << std::endl;
			}			
			
			//Then ask if the event passed the EF trigger chains and determine which ones they were			
			if ((*it).find("EF") != std::string::npos) {				
				
				//Print out the trigger name
				m_passed_EF_Trigger = true;				
				//std::cout << "Trigger Analysis: " << *it << std::endl;
				
				//Find Event Filter chains corresponding to single jet triggers (keeping it simple)
				if (((*it).find("EF_j") != std::string::npos && (*it).find("EF_je") == std::string::npos) ||
					 (*it).find("EF_fj") != std::string::npos ||
					 ((*it).find("EF_L1J") != std::string::npos && (*it).find("EMPTY") == std::string::npos) ||
					 (*it).find("EF_l2j") != std::string::npos ||
					 (*it).find("EF_b")	!= std::string::npos) {
					m_passed_EF_SingleJet_Trigger = true;
				}
				//Find Event Filter chains corresponding to multiple jet triggers (worry about it later)
				if ((*it).find("EF_2j") != std::string::npos || 
					(*it).find("EF_4j") != std::string::npos) {
					m_passed_EF_MultiJet_Trigger = true;
				}
				//Find Event Filter chains corresponding to single electrons or photons
				if (((*it).find("EF_e") != std::string::npos || (*it).find("EF_eb") == std::string::npos) ||
					(*it).find("EF_g") != std::string::npos) {					
					m_passed_EF_SingleEgamma_Trigger = true;
				}
				//Find Event Filter chains corresponding to multiple electrons or photons
				if ((*it).find("EF_2e") != std::string::npos || 
					(*it).find("EF_4e")	!= std::string::npos ||
					(*it).find("EF_2g") != std::string::npos || 
					(*it).find("EF_4g")	!= std::string::npos) {
					m_passed_EF_MultiEgamma_Trigger = true;
				}
				//Find Event Filter chains corresponding to taus
				if ((*it).find("EF_tau") != std::string::npos || 
					(*it).find("EF_2tau") != std::string::npos || 
					(*it).find("EF_4tau") != std::string::npos) {
					m_passed_EF_Tau_Trigger = true;
				}
				//Find Event Filter chains corresponding to missing energy
				//Missing energy could come from electron so it may bias results 
				if ((*it).find("EF_xe") != std::string::npos || 
					(*it).find("EF_xs") != std::string::npos) {
					m_passed_EF_MissingEnergy_Trigger = true;
				}
				
				//Any event filter chains which do not pass any of these checks
				//would be treated as unbiased triggers
			}
		}
	}

	return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------------
 // Convert offline analysis phi into value more suitable for L1Calo  
//------------------------------------------------------------------------------------
double JetEfficienciesMonTool::l1caloPhi(const double phi) const{
  if(phi >= 0.0 && phi < M_PI ){
    return phi;
  }
  if(phi < 0.0){
    return phi + (2.0 * M_PI);
  }
  return phi;
}

//------------------------------------------------------------------------------------
 // Return number of primary vertices that have at least a number of tracks (python configurable)  
//------------------------------------------------------------------------------------
unsigned int JetEfficienciesMonTool::nPrimaryVertex(){
  unsigned int nPriVtx(0);
  typedef VxContainer::const_iterator Itr;
  for(Itr i=m_primaryVertex->begin();i!=m_primaryVertex->end();++i){

    if( (*i)->vertexType() == 1 || (*i)->vertexType() == 3 ){
      if( (*i)->vxTrackAtVertex()->size() >= m_nTracksAtPrimaryVertex ){
          nPriVtx++;
      }
    }
  }
  return nPriVtx;
}  

//------------------------------------------------------------------------------------
 // Adjust entries in TileTTL1Cell container to a more suitable format
//------------------------------------------------------------------------------------
void JetEfficienciesMonTool::preCacheTTL1Cell(const TileTTL1CellContainer* cont) {
  if(!cont) return;

  TileTTL1CellContainer::const_iterator it(cont->begin());
  TileTTL1CellContainer::const_iterator itE(cont->end());

  m_idTTL1CellMap.clear();

  for(;it != itE; ++it) {
    m_idTTL1CellMap.insert(IdTTL1CellMapType::value_type((*it)->TTL1_ID(), *it));
  }
}

//------------------------------------------------------------------------------------
// Load important containers
//------------------------------------------------------------------------------------
StatusCode JetEfficienciesMonTool::loadContainers() {
	StatusCode sc;

	m_eventInfo = 0;
	sc = evtStore()->retrieve(m_eventInfo);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load EventInfo" << endreq;
		return sc;
	}

	m_primaryVertex = 0;
	sc = evtStore()->retrieve(m_primaryVertex, m_primaryVertexLocation);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load Primary Vertices" << endreq;
		return sc;
	}

	m_lvl1RoIs = 0;
	sc = evtStore()->retrieve(m_lvl1RoIs, m_lvl1RoIsLocation);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load LVL1 RoIs" << endreq;
		return sc;
	}
	m_emTauROIs.clear();
	m_jetROIs.clear();
	m_emTauROIs = m_lvl1RoIs->getEmTauROIs();
	m_jetROIs = m_lvl1RoIs->getJetROIs();

	m_offlineJets = 0;
	sc = evtStore()->retrieve(m_offlineJets, m_offlineJetsLocation);
	if (sc.isFailure()) {
		msg(MSG::WARNING) << "Failed to load Offline Jets" << endreq;
		return sc;
	}

	return sc;
}

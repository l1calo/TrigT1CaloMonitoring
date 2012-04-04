// ********************************************************************
//
// NAME:     JetEfficienciesMonTool.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Hardeep Bansil
//           Adapted for monitoring: Peter Faulkner
//	     
//
// ********************************************************************
#ifndef TRIGT1CALOMONITORING_JETEFFICIENCIESMONTOOL_H
#define TRIGT1CALOMONITORING_JETEFFICIENCIESMONTOOL_H

#include <string>
#include <vector>
#include <map>

#include "DataModel/DataVector.h"
#include "GaudiKernel/ToolHandle.h"
#include "Identifier/Identifier.h"
#include "TileEvent/TileTTL1CellContainer.h"
#include "AnalysisTriggerEvent/EmTau_ROI.h"
#include "AnalysisTriggerEvent/Jet_ROI.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

class TH1F_LW;
class TH2F_LW;
class StatusCode;
class TrigT1CaloLWHistogramTool;
class EventInfo;
class CondAttrListCollection;
class LVL1_ROI;
class JetCollection;
class Jet;
class TileTTL1Cell;
class VxContainer;

namespace LVL1 {
  class IL1TriggerTowerTool;
  class IL1CaloLArTowerEnergy;
  class TriggerTower;
}
namespace Trig {
  class TrigDecisionTool;
}

class JetEfficienciesMonTool: public ManagedMonitorToolBase
{

public:
  
  JetEfficienciesMonTool(const std::string & type, const std::string & name,
   		         const IInterface* parent);
    

  virtual ~JetEfficienciesMonTool();

  virtual StatusCode initialize();
    
  virtual StatusCode bookHistograms(bool isNewEventsBlock, bool isNewLumiBlock,
                                                           bool isNewRun);
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndOfEventsBlock,
                                    bool isEndOfLumiBlock, bool isEndOfRun);

private:

  //HSB - functions
  //----------------------------------
  double l1caloPhi(const double atlasPhi) const;

  bool deltaMatch(double etaJet, double dR, double dPhi, bool fRoi);
  double correctDeltaPhi(double dPhi);
  double calcDeltaR(double eta1, double phi1, double eta2, double phi2);
	        
  bool correctJetQuality(const Jet* jet);
  
  bool isolatedJetObjectL1(double phi, double eta);
  bool isolatedJetObjectEF(double phi, double eta);    
								        
  StatusCode analyseOfflineJets();
									    
  StatusCode triggerTowerAnalysis();
  StatusCode triggerChainAnalysis();  
  StatusCode loadContainers();
  unsigned int nPrimaryVertex();
  void preCacheTTL1Cell(const TileTTL1CellContainer* cont);

  //----------------------------------    

  static const int JET_ROI_BITS = 8;
  static const int FJET_ROI_BITS = 4;

  ToolHandle<TrigT1CaloLWHistogramTool> m_histTool;
  ToolHandle<LVL1::IL1TriggerTowerTool> m_ttTool;
  ToolHandle<LVL1::IL1CaloLArTowerEnergy> m_larEnergy;
  ToolHandle<Trig::TrigDecisionTool> m_trigger;

  // Configured chains
  std::vector<std::string> m_configuredChains; 
  // Trigger strings
  std::vector<std::string> m_triggerStrings;

  std::string m_dbPpmDeadChannelsFolder;
  std::string m_triggerTowersLocation;
  std::string m_tileTTL1ContainerLocation;
  std::string m_lvl1RoIsLocation;
  std::string m_offlineJetsLocation;
  std::string m_primaryVertexLocation;

  const EventInfo* m_eventInfo;
  const CondAttrListCollection* m_dbPpmDeadChannels;
  const DataVector<LVL1::TriggerTower>* m_triggerTowers;
  const LVL1_ROI* m_lvl1RoIs;
  const JetCollection* m_offlineJets;
  const TileTTL1CellContainer* m_tileTTL1Container;
  const VxContainer* m_primaryVertex;

  std::vector<EmTau_ROI> m_emTauROIs;
  std::vector<Jet_ROI> m_jetROIs;

  // Tile Calorimeter
  typedef std::map<Identifier, const TileTTL1Cell*> IdTTL1CellMapType;
  IdTTL1CellMapType m_idTTL1CellMap;

  /// Root directory
  std::string m_rootDir;

  // Counters for number of events
  unsigned int m_numEvents;
  unsigned int m_numOffJets;
  unsigned int m_numOffJetsInContainer;
  unsigned int m_numOffJetsPassCuts;
  unsigned int m_numOffJetsTriggered;
  unsigned int m_numJetObjPassTrigger;
  unsigned int m_numJetObjTotal;
  
  // Variables for L1 & EF trigger chains 
  bool m_passed_L1_EM_Trigger;
  bool m_passed_L1_Jet_Trigger;
  bool m_passed_EF_Trigger;
  bool m_passed_EF_SingleJet_Trigger;
  bool m_passed_EF_SingleEgamma_Trigger;
  bool m_passed_EF_MultiJet_Trigger;
  bool m_passed_EF_MultiEgamma_Trigger;
  bool m_passed_EF_Tau_Trigger;
  bool m_passed_EF_MissingEnergy_Trigger;

  // Python settable cuts  
  bool m_useEmThresholdsOnly;
            
  double m_goodEMDeltaRMatch_Cut;
  double m_goodHadDeltaRMatch_Cut;
  double m_goodHadDeltaPhiMatch_Cut;
  int m_jetQualityLevel;

  bool m_useTrigger;
  unsigned int m_nTracksAtPrimaryVertex;

  // Flag for disabled channel/bad calo analysis
  bool m_firstEvent;

  //=======================
  //   Histograms
  //=======================

  // Em Scale Jet plots
  TH1F_LW* m_h_JetEmScale_Et;
  TH1F_LW* m_h_JetEmScale_Et_central;
  TH1F_LW* m_h_JetEmScale_Et_forward;
  TH1F_LW* m_h_JetEmScale_Et_triggered;
  TH1F_LW* m_h_JetEmScale_Et_triggered_central; 
  TH1F_LW* m_h_JetEmScale_Et_triggered_forward;
  TH1F_LW* m_h_JetEmScale_Et_triggered_Eff;
  TH1F_LW* m_h_JetEmScale_Et_triggered_central_Eff; 
  TH1F_LW* m_h_JetEmScale_Et_triggered_forward_Eff;
  TH1F_LW* m_h_JetEmScale_Et_J_item[JET_ROI_BITS];
  TH1F_LW* m_h_JetEmScale_Et_FJ_J_item[FJET_ROI_BITS];
  TH1F_LW* m_h_JetEmScale_Et_J_Eff_item[JET_ROI_BITS];
  TH1F_LW* m_h_JetEmScale_Et_FJ_J_Eff_item[FJET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_Eta_vs_Phi;
  TH2F_LW* m_h_JetEmScale_Eta_vs_Phi_triggered;
  TH2F_LW* m_h_JetEmScale_Eta_vs_Phi_triggered_Eff;
  TH2F_LW* m_h_JetEmScale_50GeV_Eta_vs_Phi;
  TH2F_LW* m_h_JetEmScale_50GeV_Eta_vs_Phi_J_item[JET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_item[FJET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_50GeV_Eta_vs_Phi_J_Eff_item[JET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_50GeV_Eta_vs_Phi_FJ_J_Eff_item[FJET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_100GeV_Eta_vs_Phi;
  TH2F_LW* m_h_JetEmScale_100GeV_Eta_vs_Phi_J_item[JET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_item[FJET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_100GeV_Eta_vs_Phi_J_Eff_item[JET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_100GeV_Eta_vs_Phi_FJ_J_Eff_item[FJET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_200GeV_Eta_vs_Phi;
  TH2F_LW* m_h_JetEmScale_200GeV_Eta_vs_Phi_J_item[JET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_item[FJET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_200GeV_Eta_vs_Phi_J_Eff_item[JET_ROI_BITS];
  TH2F_LW* m_h_JetEmScale_200GeV_Eta_vs_Phi_FJ_J_Eff_item[FJET_ROI_BITS];

  //Trigger tower bad/dead
  TH2F_LW* m_h_TrigTower_jetBadCalo;
  TH2F_LW* m_h_TrigTower_jetDeadChannel;

  // Primary Vertex
  TH1F_LW* m_h_nPriVtx;

  int m_linkedHistos[FJET_ROI_BITS];

};

#endif

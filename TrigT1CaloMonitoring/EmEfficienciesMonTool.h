// ********************************************************************
//
// NAME:     EmEfficienciesMonTool.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Hardeep Bansil
//           Adapted for monitoring: Peter Faulkner
//	     
//
// ********************************************************************
#ifndef TRIGT1CALOMONITORING_EMEFFICIENCIESMONTOOL_H
#define TRIGT1CALOMONITORING_EMEFFICIENCIESMONTOOL_H

#include <string>
#include <vector>

#include "DataModel/DataVector.h"
#include "GaudiKernel/ToolHandle.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

class LWHist;
class TH1F_LW;
class TH2F_LW;
class StatusCode;
class CaloCluster;
class TrigT1CaloLWHistogramTool;
class EventInfo;
class CondAttrListCollection;
class LVL1_ROI;
class ElectronContainer;
class PhotonContainer;
class VxContainer;

namespace Analysis {
  class Electron;
  class Photon;
}
namespace LVL1 {
  class IL1TriggerTowerTool;
  class IL1CaloLArTowerEnergy;
  class TriggerTower;
}
namespace Trig {
  class TrigDecisionTool;
}

class EmEfficienciesMonTool: public ManagedMonitorToolBase
{

public:
  
  EmEfficienciesMonTool(const std::string & type, const std::string & name,
  		        const IInterface* parent);
    

  virtual ~EmEfficienciesMonTool();

  virtual StatusCode initialize();
    
  virtual StatusCode bookHistograms(bool isNewEventsBlock, bool isNewLumiBlock,
                                                           bool isNewRun);
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms(bool isEndOfEventsBlock,
                                    bool isEndOfLumiBlock, bool isEndOfRun);

private:

  //HSB - functions
  //----------------------------------
  bool vertexRequirementsPassed(int &numVtx, int &bestNumTracks);  

  bool coordinateMatch(double eta, double phi);
  bool deltaMatch(double dR, double dEta, double dPhi, double gdR, double gdEta, double gdPhi); 
  double correctDeltaPhi(double dPhi);
  double calcDeltaR(double eta1, double phi1, double eta2, double phi2);
	        
  bool emObjInDeadBadTower(double eta, double phi);
			        
  bool inEgammaGoodEtaRange(double eta);    
  bool inEMBarrel(double eta, int sign);
  bool inEMTransR(double eta, int sign);
  bool inEMEndcap(double eta, int sign);
				        
  bool correctIsEmElectron(const Analysis::Electron* el);
  bool correctIsEmPhoton(const Analysis::Photon* ph);
  std::string isEmLevelElectron(const Analysis::Electron* el, int &code);
  std::string isEmLevelPhoton(const Analysis::Photon* ph, int &code);    
				    
  bool isolatedEmObjectL1(double phi, double eta);
  bool isolatedEmObjectEF(double phi, double eta);    
								        
  std::vector<double> getRawClusterValuesFromCells(CaloCluster* cc);
								                
  StatusCode analyseOfflineElectrons();
  StatusCode analyseOfflinePhotons();  
									    
  StatusCode triggerTowerAnalysis();
  StatusCode triggerChainAnalysis();  
  StatusCode loadContainers();

  bool emType(int bitNumber);
  //----------------------------------    

  static const int ROI_BITS = 16;

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
  std::string m_lvl1RoIsLocation;
  std::string m_offlineElectronsLocation;
  std::string m_offlinePhotonsLocation;
  std::string m_primaryVertexLocation;

  const EventInfo* m_eventInfo;
  const CondAttrListCollection* m_dbPpmDeadChannels;
  const DataVector<LVL1::TriggerTower>* m_triggerTowers;
  const LVL1_ROI* m_lvl1RoIs;
  const ElectronContainer* m_offlineElectrons;
  const PhotonContainer* m_offlinePhotons;
  const VxContainer* m_primaryVtx;

  /// Root directory
  std::string m_rootDir;

  // Counters for number of events
  unsigned int m_numEvents;
  unsigned int m_numOffElec;
  unsigned int m_numOffPhot;    
  unsigned int m_numOffElecInContainer;
  unsigned int m_numOffPhotInContainer;    
  unsigned int m_numOffElecPassCuts;
  unsigned int m_numOffPhotPassCuts;
  unsigned int m_numOffElecTriggered;
  unsigned int m_numOffPhotTriggered;
  unsigned int m_numEmObjPassTrigger;
  unsigned int m_numEmObjTotal;
  
  // Variables for L1 & EF trigger chains 
  bool m_passed_L1_Jet_Trigger;
  bool m_passed_EF_SingleJet_Trigger;
  bool m_passed_EF_MultiJet_Trigger;
  bool m_passed_EF_egTau_Trigger;
  bool m_passed_EF_Trigger;

  // Python settable cuts  
  bool m_useEmThresholdsOnly;
            
  double m_goodEMDeltaRMatch_Cut;
  double m_goodEMDeltaEtaMatch_Cut;
  double m_goodEMDeltaPhiMatch_Cut;        
  double m_goodHadDeltaRMatch_Cut;
  double m_goodHadDeltaEtaMatch_Cut;
  double m_goodHadDeltaPhiMatch_Cut;
  bool m_useDeltaRMatch;
  bool m_useDeltaEtaPhiMatch;
  int m_isEmType;
  bool m_removeNoiseBursts;
       
  bool m_useEmTRcut;    
  bool m_useTrigger;
  unsigned int m_emBitMask;

  // disabled channel/bad calo flag
  bool m_firstEvent;

  //=======================
  //   Histograms
  //=======================

  TH1F_LW* m_h_ClusterRaw_Et_gdEta;    
  TH1F_LW* m_h_ClusterRaw_Et_triggered_gdEta;    
  TH1F_LW* m_h_ClusterRaw_Et_triggered_Eff;    
  TH1F_LW* m_h_ClusterRaw_Et_transR;    
  TH1F_LW* m_h_ClusterRaw_Et_triggered_transR;
  TH1F_LW* m_h_ClusterRaw_Et_triggered_transR_Eff;  
  TH1F_LW* m_h_ClusterRaw_Et_bitcheck[ROI_BITS];
  TH1F_LW* m_h_ClusterRaw_Et_bitcheck_Eff[ROI_BITS];
  TH2F_LW* m_h_ClusterRaw_10GeV_Eta_vs_Phi;      
  TH2F_LW* m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig[ROI_BITS];       
  //TH2F_LW* m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig[ROI_BITS]; 
  TH2F_LW* m_h_ClusterRaw_10GeV_Eta_vs_Phi_trig_Eff[ROI_BITS];       
  //TH2F_LW* m_h_ClusterRaw_10GeV_Eta_vs_Phi_noDeadBad_trig_Eff[ROI_BITS]; 
  TH2F_LW* m_h_ClusterRaw_20GeV_Eta_vs_Phi;      
  TH2F_LW* m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig[ROI_BITS];       
  //TH2F_LW* m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig[ROI_BITS];
  TH2F_LW* m_h_ClusterRaw_20GeV_Eta_vs_Phi_trig_Eff[ROI_BITS];       
  //TH2F_LW* m_h_ClusterRaw_20GeV_Eta_vs_Phi_noDeadBad_trig_Eff[ROI_BITS];
  TH2F_LW* m_h_ClusterRaw_30GeV_Eta_vs_Phi;      
  TH2F_LW* m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig[ROI_BITS];
  //TH2F_LW* m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig[ROI_BITS];         
  TH2F_LW* m_h_ClusterRaw_30GeV_Eta_vs_Phi_trig_Eff[ROI_BITS];
  //TH2F_LW* m_h_ClusterRaw_30GeV_Eta_vs_Phi_noDeadBad_trig_Eff[ROI_BITS];         
   
  TH2F_LW* m_h_TrigTower_emBadCalo;   
  TH2F_LW* m_h_TrigTower_emDeadChannel;       
  TH1F_LW* m_h_LAr_Noisy;


};

#endif

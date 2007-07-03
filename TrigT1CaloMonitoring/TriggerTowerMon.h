#// ********************************************************************
//
// NAME:     TriggerTowerMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Ethan-Etienne Woehrling (eew@hep.ph.bham.ac.uk)
//           Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//	     
//
// ********************************************************************
#ifndef TRIGGERTOWERMON_H
#define TRIGGERTOWERMON_H

#include <map>
#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CaloIdentifier/CaloIdManager.h"
#include "CaloIdentifier/CaloLVL1_ID.h"
#include "Identifier/Identifier.h"


#include "TH1.h"
#include "TH2.h"
//class StoreGateSvc;

class TriggerTowerMon: public ManagedMonitorToolBase
{

 public:
  
  TriggerTowerMon(const std::string & type, const std::string & name,
		  const IInterface* parent);
    

  virtual ~TriggerTowerMon();

  virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

private:
  // to get eta/phi position from the offline TT Id fields
  double IDeta(const Identifier& id);
  double IDphi(const Identifier& id);

  int pos_neg_z(double eta) const;
  int region(double eta) const; 

  double etaWidth(double eta) const;
  double phiWidth(double eta) const; 

  int etaIndex(double eta) const; 
  int phiIndex(double eta, double phi) const; 

  std::string m_TriggerTowerContainerName;

  std::string m_DataType;
  std::string m_PathInRootFile;
  
  //=====================================
  //   T1Calo Control Plots
  //=====================================
  
  std::map <Identifier, TH1F*>  m_h_TT_EmADCPeak;
  std::map <Identifier, TH1F*>  m_h_TT_HadADCPeak;
  std::map <Identifier, TH1F*>  m_h_TT_EmLUTPeak;
  std::map <Identifier, TH1F*>  m_h_TT_HadLUTPeak;
  
protected:
   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   // StoreGate service
   StoreGateSvc* m_detStore;;
   // Calorimeter Id manager
   const CaloIdManager* m_caloMgr;
   // CaloLVL1_ID Id helper
   const CaloLVL1_ID* m_lvl1Helper;


   //hitmaps

   TH2F* m_h_TT_EmHitMap_30GeV;
   TH2F* m_h_TT_EmHitMap_10GeV;
   TH2F* m_h_TT_EmHitMap_5GeV;
   TH2F* m_h_TT_EmHitMap_2GeV;
   TH2F* m_h_TT_EmHitMap_1GeV;

   TH2F* m_h_TT_HadHitMap_30GeV;
   TH2F* m_h_TT_HadHitMap_10GeV;
   TH2F* m_h_TT_HadHitMap_5GeV;
   TH2F* m_h_TT_HadHitMap_2GeV;
   TH2F* m_h_TT_HadHitMap_1GeV;

};

#endif

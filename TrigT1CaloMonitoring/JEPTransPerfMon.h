// ********************************************************************
//
// NAME:        JEPTransPerfMon.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of Transmission and Performance of JEP
//
// ********************************************************************

#ifndef JEPTransPerfMon_H
#define JEPTransPerfMon_H

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TH1.h"
#include "TH2.h"

#include "DataModel/DataVector.h"
#include "TrigT1Calo/CMMJetHits.h"
#include "TrigT1Calo/CMMEtSums.h"
#include "TrigT1Calo/JEMHits.h"
#include "TrigT1Calo/JEMEtSums.h"
#include "TrigT1Calo/JetElement.h"
//#include "TrigT1CaloTools/L1JetElementTools.h"

#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"



#include "CaloIdentifier/CaloIdManager.h"
#include "CaloIdentifier/CaloLVL1_ID.h"
#include "Identifier/Identifier.h"
#include "TrigT1CaloCalibTools/L1CaloTTIdTools.h"
#include "CaloTriggerTool/CaloTriggerTowerService.h"



class StoreGateSvc;

namespace LVL1 {
  class CMMEtSums;
  class CMMJetHits;
  class CMMRoI;
  class JEMEtSums;
  class JEMHits;
  class JEMRoI;
  class JetAlgorithm;
  class JetElement;
  class TriggerTower;
  class IL1JEPHitsTools;
  class IL1JetElementTools;
  class IL1JetTools;
  class IL1JEPEtSumsTools;
  class EventInfo;

}

class JEPTransPerfMon : public ManagedMonitorToolBase
{
public:

 	JEPTransPerfMon( const std::string & type, const std::string & name,
	                 const IInterface* parent ); 
	typedef DataVector<LVL1::JetElement> JECollection;
/*         typedef DataVector<const LVL1::CMMJetHits*> pCMMJetHitsCollection; */
        typedef DataVector<LVL1::CMMJetHits> CMMJetHitsCollection;
        typedef DataVector<const LVL1::CMMJetHits> cCMMJetHitsCollection;
        typedef DataVector<LVL1::CMMEtSums> CMMEtSumsCollection;
        typedef DataVector<LVL1::JEMHits> JEMHitsCollection;
	typedef DataVector<LVL1::JEMEtSums> JEMEtSumsCollection;
	typedef DataVector<LVL1::JEMRoI> JemRoiCollection;
	typedef DataVector<LVL1::JetAlgorithm> InternalJemRoi;
	
	

	virtual ~JEPTransPerfMon();
	virtual StatusCode initialize();
	virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
	virtual StatusCode fillHistograms();
	virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );


protected:

 void  TimeSliceMatch(int k, int TT_TS, const JECollection* TT_jetElements, int JE_TS, const JECollection* jetElements, MsgStream::MsgStream* mLog);

	
   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;
   StoreGateSvc* m_eventStore;

   // StoreGate service
   StoreGateSvc* m_detStore;
   // Calorimeter Id manager
   const CaloIdManager* m_caloMgr;
   // CaloLVL1_ID Id helper
   const CaloLVL1_ID* m_lvl1Helper;
   const L1CaloTTIdTools* m_l1CaloTTIdTools;
   
   CaloTriggerTowerService* m_ttSvc;
   // TTOnlineID Id helper
   const TTOnlineID* m_l1ttonlineHelper;
   

   int m_NoEvents;
   
  ToolHandle<LVL1::IL1JEPHitsTools>    m_jepHitsTool;
  ToolHandle<LVL1::IL1JetTools>        m_jetTool;
  ToolHandle<LVL1::IL1JetElementTools> m_jetElementTool;
  ToolHandle<LVL1::IL1JEPEtSumsTools>  m_etSumsTool;
  mutable MsgStream mLog;
  
  
   /** location of data */
   std::string m_BS_JetElementLocation;
   std::string m_BS_TriggerTowerLocation;
   int m_NoLUTSlices;
   bool  m_EventNoInHisto;
   bool m_Offline;

   std::string m_BS_JEMHitsLocation;
   std::string m_Sim_JEMHitsLocation;

   std::string m_BS_JEMEtSumsLocation;   
   std::string m_Sim_JEMEtSumsLocation;   

   std::string m_BS_JEMRoILocation;   
   std::string m_Sim_JEMRoILocation;   

   std::string m_BS_CMMJetHitsLocation;
   std::string m_Sim_CMMJetHitsLocation;

   std::string m_BS_CMMEtSumsLocation;
   std::string m_Sim_CMMEtSumsLocation;


   std::string m_BS_CMMRoILocation;
   std::string m_Sim_CMMRoILocation;

   std::string m_DataType;   
   std::string m_PathInRootFile;   
   std::string m_ErrorPathInRootFile;

   bool m_CompareWithSimulation ;

   /** location of data */


  /** Histos */   
   // JEP 
   TH2F* m_h_SimBSMon_JEP;

   TH2F* m_h_TransCheck_JEP;
   TH2F* m_h_TransCheck_emJetElements;
   TH2F* m_h_TransCheck_hadJetElements;
  
  std::stringstream runNumStr;
  
  int m_evtNum;
  int m_lumiBlock;
  int m_evtBCID;
  int m_runNum;


};


#endif

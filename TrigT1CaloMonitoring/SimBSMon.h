// ********************************************************************
//
// NAME:        SimBSMon.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on CMM level
//
// ********************************************************************

#ifndef SimBSMon_H
#define SimBSMon_H

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TH1.h"
#include "TH2.h"

#include "DataModel/DataVector.h"
#include "TrigT1Calo/CMMJetHits.h"
#include "TrigT1Calo/CMMEtSums.h"
#include "TrigT1Calo/JEMHits.h"
#include "TrigT1Calo/JEMEtSums.h"
#include "TrigT1Calo/JEMRoI.h"


#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"

namespace LVL1 {
  class JEMRoI;
}

class SimBSMon : public ManagedMonitorToolBase
{
public:
        typedef DataVector<LVL1::CMMJetHits> CMMJetHitsCollection;
        typedef DataVector<const LVL1::CMMJetHits> cCMMJetHitsCollection;
        typedef DataVector<LVL1::CMMEtSums> CMMEtSumsCollection;
        typedef DataVector<LVL1::JEMHits> JEMHitsCollection;
	typedef DataVector<LVL1::JEMEtSums> JEMEtSumsCollection;
	typedef DataVector<LVL1::JEMRoI> JEMRoICollection;


	SimBSMon( const std::string & type, const std::string & name,
	                 const IInterface* parent ); 

	virtual ~SimBSMon();

	virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
	virtual StatusCode fillHistograms();
	virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

protected:

   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   /** location of data */
   std::string m_BS_CMMJetHitsLocation;
   std::string m_BS_CMMEtSumsLocation;
   std::string m_BS_CMMRoILocation;

   std::string m_BS_JEMHitsLocation;
   std::string m_BS_JEMEtSumsLocation;   
   std::string m_BS_JEMRoILocation;   

   std::string m_Sim_CMMJetHitsLocation;
   std::string m_Sim_CMMEtSumsLocation;
   std::string m_Sim_CMMRoILocation;

   std::string m_Sim_JEMHitsLocation;
   std::string m_Sim_JEMEtSumsLocation;   
   std::string m_Sim_JEMRoILocation;   

   std::string m_DataType;   
   std::string m_PathInRootFile;   

  /** Histos */   
   // JEM 
   TH1F*  m_h_SimBSMon_JEM_Crate0_Energy;
   TH1F*  m_h_SimBSMon_JEM_Crate1_Energy;

   TH1F*  m_h_SimBSMon_JEM_Crate0_Hits;
   TH1F*  m_h_SimBSMon_JEM_Crate1_Hits;

   TH1F*  m_h_SimBSMon_JEM_Crate0_RoI;
   TH1F*  m_h_SimBSMon_JEM_Crate1_RoI;


};


#endif

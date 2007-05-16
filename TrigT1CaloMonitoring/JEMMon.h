// ********************************************************************
//
// NAME:        JEMMon.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on JEM level
//
// ********************************************************************

#ifndef JEMMon_H
#define JEMMon_H

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TH1.h"
#include "TH2.h"

#include "DataModel/DataVector.h"
#include "TrigT1Calo/JEMHits.h"
#include "TrigT1Calo/JEMEtSums.h"

#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"


class JEMMon : public ManagedMonitorToolBase
{
public:
        typedef DataVector<LVL1::JEMHits> JEMHitsCollection;
	typedef DataVector<LVL1::JEMEtSums> JEMEtSumsCollection;
       

	JEMMon( const std::string & type, const std::string & name,
	                 const IInterface* parent ); 

	virtual ~JEMMon();

	virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
	virtual StatusCode fillHistograms();
	virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

protected:

   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   /** location of data */
   std::string m_BS_JEMHitsLocation;
   std::string m_BS_JEMEtSumsLocation;   

   std::string m_Sim_JEMHitsLocation;
   std::string m_Sim_JEMEtSumsLocation;   

   /** Histos */

   // JEM Hits
   TH1F* m_h_BS_JEMHits_MainHits;
   TH1F* m_h_BS_JEMHits_FwdHitsRight;
   TH1F* m_h_BS_JEMHits_FwdHitsLeft;

   TH1F* m_h_Sim_JEMHits_MainHits;
   TH1F* m_h_Sim_JEMHits_FwdHitsRight;
   TH1F* m_h_Sim_JEMHits_FwdHitsLeft;

   // JEM Et Sums
   TH1F*  m_h_BS_JEMEtSums_Ex;
   TH1F*  m_h_BS_JEMEtSums_Ey;
   TH1F*  m_h_BS_JEMEtSums_Et;

   TH1F*  m_h_Sim_JEMEtSums_Ex;
   TH1F*  m_h_Sim_JEMEtSums_Ey;
   TH1F*  m_h_Sim_JEMEtSums_Et;
 
};


#endif

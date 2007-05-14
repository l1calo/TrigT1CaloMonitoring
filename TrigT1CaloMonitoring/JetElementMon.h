// ********************************************************************
//
// NAME:        JetElementMon.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the inputdata (JetElements) of the JEP
//
// ********************************************************************

#ifndef JetElementMon_H
#define JetElementMon_H

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TH1.h"
#include "TH2.h"

#include "DataModel/DataVector.h"
#include "TrigT1Calo/JetElementMaker.h"

#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"


class JetElementMon : public ManagedMonitorToolBase
{
public:
        typedef DataVector<LVL1::JetElement> JECollection;

	JetElementMon( const std::string & type, const std::string & name,
	                 const IInterface* parent ); 

	virtual ~JetElementMon();

	virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
	virtual StatusCode fillHistograms();
	virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

protected:

   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   /** location of data */
   std::string m_BS_JetElementLocation;
   std::string m_Sim_JetElementLocation;

   /** Histos */
   TH1F* m_h_BS_je_eta;
   TH1F* m_h_BS_je_phi;
   TH1F* m_h_BS_je_emenergy;
   TH1F* m_h_BS_je_hadenergy; 
   TH1F* m_h_BS_je_energy; 

   TH2F* m_h_BS_je_etaphi;
   TH2F* m_h_BS_je_energy_etaphi;

   TH1F* m_h_Sim_je_eta;
   TH1F* m_h_Sim_je_phi;
   TH1F* m_h_Sim_je_emenergy;
   TH1F* m_h_Sim_je_hadenergy; 
   TH1F* m_h_Sim_je_energy; 

   TH2F* m_h_Sim_je_etaphi;
   TH2F* m_h_Sim_je_energy_etaphi;

};


#endif

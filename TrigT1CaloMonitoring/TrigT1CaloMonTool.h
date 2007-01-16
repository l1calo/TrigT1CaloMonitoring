// ********************************************************************
//
// NAME:     TrigT1CaloMonTool.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Ethan-Etienne Woehrling (eew@hep.ph.bham.ac.uk)
//	     
//
// ********************************************************************
#ifndef TRIGT1CALOMONTOOL_H
#define TRIGT1CALOMONTOOL_H

#include <map>
#include "AthenaMonitoring/MonitorToolBase.h"

class TH1D;
class TH2D;
class StoreGateSvc;

class TrigT1CaloMonTool: public MonitorToolBase
{

 public:
  
  TrigT1CaloMonTool(const std::string & type, const std::string & name,
		  const IInterface* parent);
    

  ~TrigT1CaloMonTool();

  StatusCode initialize();
    
  StatusCode bookHists();
  StatusCode fillHists();
  StatusCode finalHists();
  StatusCode checkHists(bool fromFinalize);

private:

  std::string m_towersContName;
  std::string m_TriggerTowerContainerName;
  std::string  m_JetElementContainerName;

  //=====================================
  //   T1Calo Control Plots
  //=====================================

  TH1D* m_h_TT_Em_Et;  
  TH1D* m_h_TT_Had_Et;
  TH1D* m_h_TT_eta;
  TH1D* m_h_TT_phi;

  TH1D* m_h_JE_Em_Et; 
  TH1D* m_h_JE_Had_Et; 
  TH1D* m_h_JE_eta;
  TH1D* m_h_JE_phi;

  TH1D* m_h_TT_Tot_Et;
  TH1D* m_h_JE_Tot_Et; 

  TH1D* m_h_TT_Em10_Et;
  TH1D* m_h_TT_Had10_Et;
  TH1D* m_h_JE_Em10_Et; 
  TH1D* m_h_JE_Had10_Et; 

  TH1D* m_h_TT_key;

  //Calo regions plots

  TH1D* m_h_Barrel_TT_phi;
  TH1D* m_h_Barrel_TT_Em_Et;
  TH1D* m_h_Barrel_TT_Had_Et;

  TH1D* m_h_Barrel_JE_phi;
  TH1D* m_h_Barrel_JE_Em_Et;
  TH1D* m_h_Barrel_JE_Had_Et;

  TH1D* m_h_EC10_TT_Em_Et;
  TH1D* m_h_EC10_TT_Had_Et;
  TH1D* m_h_EC10_JE_Em_Et;
  TH1D* m_h_EC10_JE_Had_Et;

  TH1D* m_h_Barrel10_TT_Em_Et;
  TH1D* m_h_Barrel10_TT_Had_Et;
  TH1D* m_h_Barrel10_JE_Em_Et;
  TH1D* m_h_Barrel10_JE_Had_Et;

  TH1D* m_h_FCAL10_TT_Em_Et;
  TH1D* m_h_FCAL10_TT_Had_Et;
  TH1D* m_h_FCAL10_JE_Em_Et;
  TH1D* m_h_FCAL10_JE_Had_Et;

  TH1D* m_h_EC_TT_phi;
  TH1D* m_h_EC_TT_Em_Et;
  TH1D* m_h_EC_TT_Had_Et;

  TH1D* m_h_EC_JE_phi;
  TH1D* m_h_EC_JE_Em_Et;
  TH1D* m_h_EC_JE_Had_Et;

  TH1D* m_h_FCAL_TT_phi;
  TH1D* m_h_FCAL_TT_Em_Et;
  TH1D* m_h_FCAL_TT_Had_Et;

  TH1D* m_h_FCAL_JE_phi;
  TH1D* m_h_FCAL_JE_Em_Et;
  TH1D* m_h_FCAL_JE_Had_Et;

  TH2D* m_h_TT_etaphi;
  TH2D* m_h_TT_etaphi_hitmap;

  TH2D* m_h_JE_etaphi;
  TH2D* m_h_JE_etaphi_hitmap;

  //Et>10Gev cut
  TH1D* m_h_TT_eta_gt10;
  TH1D* m_h_TT_phi_gt10; 
  TH1D* m_h_TT_EC_phi_gt10; 
  TH1D* m_h_TT_Barrel_phi_gt10; 
  TH1D* m_h_TT_FCAL_phi_gt10; 

  TH1D* m_h_JE_eta_gt10;
  TH1D* m_h_JE_phi_gt10; 
  TH1D* m_h_JE_EC_phi_gt10; 
  TH1D* m_h_JE_Barrel_phi_gt10; 
  TH1D* m_h_JE_FCAL_phi_gt10; 


  //  ==============================
  //  T1Calo comparison with Calo 
  //  ==============================

  // Using CaloTowers from the ESD - CaloTower Plots:
  TH1D* m_h_CaloT_Et;
  TH1D* m_h_CaloT_Et10;
  TH1D* m_h_CaloT_phi;
  TH1D* m_h_CaloT_eta;
  TH1D* m_h_CaloT_phi_gt10;
  TH1D* m_h_CaloT_eta_gt10;
  TH2D* m_h_CaloT_etaphi;
  TH2D* m_h_CaloT_etaphi_hitmap;
  TH1D* m_h_CaloT_key;

  //Trigger-Style-Calo-Towers + TT plots
  //T-S-C-Ts are filled from calocell esd info, but in the style of a TriggerTower

  //TSCT control plots

  TH1D* m_h_Calo_Em_Et;
  TH1D* m_h_Calo_Em_Et10;
  TH1D* m_h_Calo_Had_Et;
  TH1D* m_h_Calo_Had_Et10;
  TH1D* m_h_Calo_phi;
  TH1D* m_h_Calo_eta;
  TH1D* m_h_Calo_phi_gt10;
  TH1D* m_h_Calo_eta_gt10;
  TH2D* m_h_Calo_etaphi;
  TH2D* m_h_Calo_etaphi_hitmap;
  TH1D* m_h_Calo_key;

  //regions
  TH1D* m_h_Barrel_Calo_Em_Et;
  TH1D* m_h_Barrel_Calo_Had_Et;
  TH1D* m_h_Barrel10_Calo_Em_Et;
  TH1D* m_h_Barrel10_Calo_Had_Et;
  TH1D* m_h_Barrel_Calo_phi;

  TH1D* m_h_EC_Calo_Em_Et;
  TH1D* m_h_EC_Calo_Had_Et;
  TH1D* m_h_EC10_Calo_Em_Et;
  TH1D* m_h_EC10_Calo_Had_Et;
  TH1D* m_h_EC_Calo_phi;

  TH1D* m_h_FCAL_Calo_Em_Et;
  TH1D* m_h_FCAL_Calo_Had_Et;
  TH1D* m_h_FCAL10_Calo_Em_Et;
  TH1D* m_h_FCAL10_Calo_Had_Et;
  TH1D* m_h_FCAL_Calo_phi;

  //Compared with TT
  //1D ratio plots
  TH1D* m_h_TT_Calo_eta;
  TH1D* m_h_TT_Calo_Et;

  //2d comparison plots
  TH2D* m_h_TT_Calo_Em_EtTower;
  TH2D* m_h_TT_Calo_Had_EtTower;
  TH2D* m_h_TT_Calo_PhiTower;
  TH2D* m_h_TT_Calo_EtaTower;

  //Discrepancy check plots

  //Plots to look at unmatched TT&Calo towers
  //the D is for Discrepancy

  TH1D* m_h_Ratio_D_Em_Et; 
  TH1D* m_h_Ratio_D_Had_Et; 

  TH1D* m_h_Calo_DEm_under_phi; 
  TH1D* m_h_Calo_DEm_under_eta; 
  TH1D* m_h_Calo_DEm_under_Em_Et; 
  TH1D* m_h_Calo_DEm_under_Had_Et; 
  TH1D* m_h_Calo_DEm_under_TTEm_Et; 
  TH1D* m_h_Calo_DEm_under_TTHad_Et; 

  TH1D* m_h_Calo_DEm_over_phi; 
  TH1D* m_h_Calo_DEm_over_eta; 
  TH1D* m_h_Calo_DEm_over_Em_Et; 
  TH1D* m_h_Calo_DEm_over_Had_Et; 
  TH1D* m_h_Calo_DEm_over_TTEm_Et; 
  TH1D* m_h_Calo_DEm_over_TTHad_Et; 

  TH1D* m_h_Calo_DHad_under_phi; 
  TH1D* m_h_Calo_DHad_under_eta; 
  TH1D* m_h_Calo_DHad_under_Em_Et; 
  TH1D* m_h_Calo_DHad_under_Had_Et; 
  TH1D* m_h_Calo_DHad_under_TTEm_Et; 
  TH1D* m_h_Calo_DHad_under_TTHad_Et; 

  TH1D* m_h_Calo_DHad_over_phi; 
  TH1D* m_h_Calo_DHad_over_eta; 
  TH1D* m_h_Calo_DHad_over_Em_Et; 
  TH1D* m_h_Calo_DHad_over_Had_Et; 
  TH1D* m_h_Calo_DHad_over_TTEm_Et; 
  TH1D* m_h_Calo_DHad_over_TTHad_Et; 

  
  //Calibration plots:

  //Et plotted against Eta
  TH2D* m_h_Calib_TTEM_EtEta;
  TH2D* m_h_Calib_TTHAD_EtEta;
  TH2D* m_h_Calib_CaloT_EtEta;
  TH2D* m_h_Calib_CaloEM_EtEta;
  TH2D* m_h_Calib_CaloHAD_EtEta;

  TH2D* m_h_Calib_EMRatio_ETEta;
  TH2D* m_h_Calib_HADRatio_ETEta;


protected:

  TH1D* book1D(std::string nam, std::string tit, int nx, double xmin, double xmax);
  TH2D* book2D(std::string nam, std::string tit, int nx, double xmin, double xmax, 
	       int ny, double ymin, double ymax);

  std::string m_stem;

  StoreGateSvc* m_StoreGate;


};

#endif

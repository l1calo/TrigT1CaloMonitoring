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

  //Calorimeter plots
  TH1D* m_h_Calo_Et;
  TH1D* m_h_Calo_Et10;
  TH1D* m_h_Calo_phi;
  TH1D* m_h_Calo_eta;
  TH1D* m_h_Calo_phi_gt10;
  TH1D* m_h_Calo_eta_gt10;

  TH2D* m_h_Calo_etaphi;
  TH2D* m_h_Calo_etaphi_hitmap;

  TH1D* m_h_Calo_key;

  //Calo + TT plots

  TH1D* m_h_TT_Calo_eta;
  TH1D* m_h_TT_Calo_Et;
  TH2D* m_h_TT_Calo_EtTower;
  TH2D* m_h_TT_Calo_PhiTower;
  TH2D* m_h_TT_Calo_EtaTower;

  //Discrepancy check plots

  //Plots to look at unmatched TT&Calo towers

  TH1D* m_h_TT_D_under_phi; 
  TH1D* m_h_TT_D_under_eta;
  TH1D* m_h_Calo_D_under_phi;
  TH1D* m_h_Calo_D_under_eta;

  TH1D* m_h_TT_D_over_phi; 
  TH1D* m_h_TT_D_over_eta;
  TH1D* m_h_Calo_D_over_phi;
  TH1D* m_h_Calo_D_over_eta;

  TH1D* m_h_TT_D_over_Et;
  TH1D* m_h_TT_D_over_EmEt;
  TH1D* m_h_TT_D_over_HadEt;
  TH1D* m_h_Calo_D_over_Et;

  TH1D* m_h_TT_D_under_Et;
  TH1D* m_h_TT_D_under_EmEt;
  TH1D* m_h_TT_D_under_HadEt;
  TH1D* m_h_Calo_D_under_Et;

protected:

  TH1D* book1D(std::string nam, std::string tit, int nx, double xmin, double xmax);
  TH2D* book2D(std::string nam, std::string tit, int nx, double xmin, double xmax, 
	       int ny, double ymin, double ymax);

  std::string m_stem;

  StoreGateSvc* m_StoreGate;


};

#endif

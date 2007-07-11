
####################################################################################################################
##################################### Bytestream Optionen ##########################################################
####################################################################################################################

include ("TrigT1CaloMonitoring/TrigT1CaloMonitoring_L1CaloBSConverter.py")

####################################################################################################################
########################################## L1Calo Simulation #######################################################
####################################################################################################################

include ("TrigT1CaloMonitoring/TrigT1CaloMonitoring_L1CaloSimulation.py")

####################################################################################################################
##################################### Monitoring ###################################################################
####################################################################################################################

# *************************
# Application configuration
# *************************

import AthenaCommon.AtlasUnixGeneratorJob

## get a handle on the top sequence of algorithms
from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()

## add an AthenaMonManager algorithm to the list of algorithms to be ran
from AthenaMonitoring.AthenaMonitoringConf import AthenaMonManager
topSequence += AthenaMonManager( "PrimaryManager" )

# ************************
# Monitoring configuration
# ************************

## Setup the output file(s):
THistSvc = Algorithm("THistSvc")

## The string "TestMon" in the argument below is the 'FileKey'
## used by Athena to access the output file internally
THistSvc.Output = ["L1CaloMon DATAFILE='L1CaloMon.root' OPT='RECREATE'"]

## AthenaMonManager is the Algorithm that manages many classes inheriting
## from ManagedMonitorToolBase
monMan = topSequence.PrimaryManager

# Add all the ManagedMonitorToolBase objects
# and configure them
## get a handle on the ToolSvc
from AthenaCommon.AppMgr import ToolSvc as toolSvc

####################### TriggerTower ################################
monMan.AthenaMonTools += [ "TriggerTowerMon/L1TriggerTowerMonTool" ]

toolSvc.L1TriggerTowerMonTool.DataType = "M3"  
toolSvc.L1TriggerTowerMonTool.BS_TriggerTowerContainer = "TriggerTowers"
toolSvc.L1TriggerTowerMonTool.PathInRootFile = "Stats/L1Calo/TriggerTower"
#toolSvc.L1TriggerTowerMonTool.OutputLevel = DEBUG

####################### JetElements ################################
monMan.AthenaMonTools += [ "JetElementMon/BS_L1JetElementMonTool" ]
monMan.AthenaMonTools += [ "JetElementMon/Sim_L1JetElementMonTool" ]

toolSvc.BS_L1JetElementMonTool.DataType = "BS"  #BS or Sim data?
toolSvc.BS_L1JetElementMonTool.JetElementLocation = "JetElements"
toolSvc.BS_L1JetElementMonTool.PathInRootFile = "Stats/L1Calo/JEM/BS/input"
#toolSvc.BS_L1JetElementMonTool.OutputLevel = DEBUG

toolSvc.Sim_L1JetElementMonTool.DataType = "Sim"  #BS or Sim data?
toolSvc.Sim_L1JetElementMonTool.JetElementLocation = "Sim_JetElements"
toolSvc.Sim_L1JetElementMonTool.PathInRootFile = "Stats/L1Calo/JEM/Sim/input"
#toolSvc.Sim_L1JetElementMonTool.OutputLevel = DEBUG

####################### JEMs ################################
monMan.AthenaMonTools += [ "JEMMon/BS_L1JEMMonTool" ]
monMan.AthenaMonTools += [ "JEMMon/Sim_L1JEMMonTool" ]

toolSvc.BS_L1JEMMonTool.DataType = "BS"  #BS or Sim data?
toolSvc.BS_L1JEMMonTool.JEMHitsLocation = "JEMHits"
toolSvc.BS_L1JEMMonTool.JEMEtSumsLocation = "JEMEtSums"
toolSvc.BS_L1JEMMonTool.JEMRoILocation = "JEMRoIs"
toolSvc.BS_L1JEMMonTool.PathInRootFile = "Stats/L1Calo/JEM/BS"
#toolSvc.BS_L1JEMMonTool.OutputLevel = DEBUG

toolSvc.Sim_L1JEMMonTool.DataType = "Sim"  #BS or Sim data?
toolSvc.Sim_L1JEMMonTool.JEMHitsLocation = "Sim_JEMHits"
toolSvc.Sim_L1JEMMonTool.JEMEtSumsLocation = "Sim_JEMEtSums"
toolSvc.Sim_L1JEMMonTool.JEMRoILocation = "Sim_JEMRoIs"
toolSvc.Sim_L1JEMMonTool.PathInRootFile = "Stats/L1Calo/JEM/Sim"
#toolSvc.Sim_L1JEMMonTool.OutputLevel = DEBUG

####################### CMMs ################################
monMan.AthenaMonTools += [ "CMMMon/BS_L1CMMMonTool" ]
monMan.AthenaMonTools += [ "CMMMon/Sim_L1CMMMonTool" ]

toolSvc.BS_L1CMMMonTool.DataType = "BS"  #BS or Sim data?
toolSvc.BS_L1CMMMonTool.CMMJetHitsLocation = "CMMJetHits"
toolSvc.BS_L1CMMMonTool.CMMEtSumsLocation = "CMMEtSums"
toolSvc.BS_L1CMMMonTool.CMMRoILocation = "CMMRoIs"
toolSvc.BS_L1CMMMonTool.JEMHitsLocation = "JEMHits"
toolSvc.BS_L1CMMMonTool.JEMEtSumsLocation = "JEMEtSums"
toolSvc.BS_L1CMMMonTool.PathInRootFile = "Stats/L1Calo/CMM/BS"
#toolSvc.BS_L1CMMMonTool.OutputLevel = DEBUG

toolSvc.Sim_L1CMMMonTool.DataType = "Sim"  #BS or Sim data?
toolSvc.Sim_L1CMMMonTool.CMMJetHitsLocation = "Sim_CMMJetHits"
toolSvc.Sim_L1CMMMonTool.CMMEtSumsLocation = "Sim_CMMEtSums"
toolSvc.Sim_L1CMMMonTool.CMMRoILocation = "Sim_CMMRoIs"
toolSvc.Sim_L1CMMMonTool.JEMHitsLocation = "Sim_JEMHits"
toolSvc.Sim_L1CMMMonTool.JEMEtSumsLocation = "Sim_JEMEtSums"
toolSvc.Sim_L1CMMMonTool.PathInRootFile = "Stats/L1Calo/CMM/Sim"
#toolSvc.Sim_L1CMMMonTool.OutputLevel = DEBUG

####################### Calorimeter ################################
#monMan.AthenaMonTools += [ "TrigT1CaloBSMonTool/L1CaloBSTool" ]
#toolSvc.L1CaloBSTool.DataType = "BS"  #BS or Sim data?
#toolSvc.L1CaloBSTool.BS_TriggerTowerContainer = "TriggerTowers"
#toolSvc.L1CaloBSTool.BS_JetElementContainer = "JetElements"
#toolSvc.L1CaloBSTool.PathInRootFile = "Stats/L1Calo/Calo"
#toolSvc.L1CaloBSTool.OutputLevel = DEBUG

####################### CPMs ################################
monMan.AthenaMonTools += [ "TrigT1CaloCpmMonTool/L1BSCPMMonTool" ]
monMan.AthenaMonTools += [ "TrigT1CaloCpmMonTool/L1SimCPMMonTool" ]

toolSvc.L1BSCPMMonTool.HistogramPrefix = "BS"
toolSvc.L1BSCPMMonTool.TriggerTowerLocation = "TriggerTowers"
toolSvc.L1BSCPMMonTool.CPMTowerLocation = "CPMTowers"
toolSvc.L1BSCPMMonTool.CPMHitsLocation = "CPMHits"
toolSvc.L1BSCPMMonTool.CMMCPHitsLocation = "CMMCPHits"
toolSvc.L1BSCPMMonTool.CPMRoILocation = "CPMRoIs"
#toolSvc.L1BSCPMMonTool.OutputLevel = DEBUG

toolSvc.L1SimCPMMonTool.HistogramPrefix = "Sim"
toolSvc.L1SimCPMMonTool.TriggerTowerLocation = "TriggerTowers"
toolSvc.L1SimCPMMonTool.CPMTowerLocation = "Sim_CPMTowers"
toolSvc.L1SimCPMMonTool.CPMHitsLocation = "Sim_CPMHits"
toolSvc.L1SimCPMMonTool.CMMCPHitsLocation = "Sim_CMMCPHits"
toolSvc.L1SimCPMMonTool.CPMRoILocation = "Sim_CPMRoIs"
#toolSvc.L1SimCPMMonTool.OutputLevel = DEBUG

####################### Simulation vs HW ################################
monMan.AthenaMonTools += [ "SimBSMon/SimBSMonTool" ]

toolSvc.SimBSMonTool.BS_JEMEtSumsLocation = "JEMEtSums"
toolSvc.SimBSMonTool.Sim_JEMEtSumsLocation = "Sim_JEMEtSums"
toolSvc.SimBSMonTool.BS_JEMHitsLocation = "JEMHits"
toolSvc.SimBSMonTool.Sim_JEMHitsLocation = "Sim_JEMHits"
toolSvc.SimBSMonTool.BS_CMMEtSumsLocation = "CMMEtSums"
toolSvc.SimBSMonTool.Sim_CMMEtSumsLocation = "Sim_CMMEtSums"
toolSvc.SimBSMonTool.BS_CMMJetHitsLocation = "CMMJetHits"
toolSvc.SimBSMonTool.Sim_CMMJetHitsLocation = "Sim_CMMJetHits"
toolSvc.SimBSMonTool.PathInRootFile = "Stats/L1Calo/JEM/SimBS"
#toolSvc.SimBSMonTool.OutputLevel = DEBUG




# FileKey must match that given to THistSvc
monMan.FileKey = "L1CaloMon"

# Set global monitoring parameters: see the AthenaMonManager class
# in the Control/AthenaMonitoring package
monMan.ManualDataTypeSetup = True
monMan.DataType            = "collisions"
monMan.Environment         = "tier0"
monMan.ManualRunLBSetup    = True
monMan.Run                 = 1
monMan.LumiBlock           = 1

####################################################################################################################
###################################### Ende ########################################################################
####################################################################################################################

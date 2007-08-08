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
#THistSvc.Output = ["L1CaloMon DATAFILE='L1CaloMon.root' OPT='RECREATE'"]

## AthenaMonManager is the Algorithm that manages many classes inheriting
## from ManagedMonitorToolBase
monMan = topSequence.PrimaryManager

# Add all the ManagedMonitorToolBase objects
# and configure them
## get a handle on the ToolSvc
from AthenaCommon.AppMgr import ToolSvc as toolSvc

####################### TriggerTower ################################
monMan.AthenaMonTools += [ "TriggerTowerMon/L1TriggerTowerMonTool" ]
#toolSvc.L1TriggerTowerMonTool.DataType = "M3"  
toolSvc.L1TriggerTowerMonTool.BS_TriggerTowerContainer = "TriggerTowers"
toolSvc.L1TriggerTowerMonTool.DistPerChannel = False
toolSvc.L1TriggerTowerMonTool.DistPerChannelAndTimeSlice = False
toolSvc.L1TriggerTowerMonTool.LUTHitMap_Thresh0 = 1
toolSvc.L1TriggerTowerMonTool.LUTHitMap_Thresh1 = 3
toolSvc.L1TriggerTowerMonTool.LUTHitMap_Thresh2 = 7
toolSvc.L1TriggerTowerMonTool.ADCHitMap_Thresh = 30
toolSvc.L1TriggerTowerMonTool.PathInRootFile = "L1Calo/PPr"
toolSvc.L1TriggerTowerMonTool.ErrorPathInRootFile = "L1Calo/Errors"
#toolSvc.L1TriggerTowerMonTool.OutputLevel = DEBUG

####################### JetElements ################################
monMan.AthenaMonTools += [ "JetElementMon/BS_L1JetElementMonTool" ]
toolSvc.BS_L1JetElementMonTool.DataType = "BS"  #BS or Sim data?
toolSvc.BS_L1JetElementMonTool.JetElementLocation = "JetElements"
toolSvc.BS_L1JetElementMonTool.PathInRootFile = "L1Calo/JEP/JEM/BS/input"
toolSvc.BS_L1JetElementMonTool.ErrorPathInRootFile = "L1Calo/Errors"
#toolSvc.BS_L1JetElementMonTool.OutputLevel = DEBUG

####################### JEMs ################################
monMan.AthenaMonTools += [ "JEMMon/BS_L1JEMMonTool" ]
toolSvc.BS_L1JEMMonTool.DataType = "BS"  #BS or Sim data?
toolSvc.BS_L1JEMMonTool.JEMHitsLocation = "JEMHits"
toolSvc.BS_L1JEMMonTool.JEMEtSumsLocation = "JEMEtSums"
toolSvc.BS_L1JEMMonTool.JEMRoILocation = "JEMRoIs"
toolSvc.BS_L1JEMMonTool.PathInRootFile = "L1Calo/JEP/JEM/BS"
toolSvc.BS_L1JEMMonTool.ErrorPathInRootFile = "L1Calo/Errors"
#toolSvc.BS_L1JEMMonTool.OutputLevel = DEBUG

####################### CMMs ################################
monMan.AthenaMonTools += [ "CMMMon/BS_L1CMMMonTool" ]
toolSvc.BS_L1CMMMonTool.DataType = "BS"  #BS or Sim data?
toolSvc.BS_L1CMMMonTool.CMMJetHitsLocation = "CMMJetHits"
toolSvc.BS_L1CMMMonTool.CMMEtSumsLocation = "CMMEtSums"
toolSvc.BS_L1CMMMonTool.CMMRoILocation = "CMMRoIs"
toolSvc.BS_L1CMMMonTool.JEMHitsLocation = "JEMHits"
toolSvc.BS_L1CMMMonTool.JEMEtSumsLocation = "JEMEtSums"
toolSvc.BS_L1CMMMonTool.PathInRootFile = "L1Calo/JEP/CMM/BS"
toolSvc.BS_L1CMMMonTool.ErrorPathInRootFile = "L1Calo/Errors"
#toolSvc.BS_L1CMMMonTool.OutputLevel = DEBUG

####################### CPMs ################################
monMan.AthenaMonTools += [ "TrigT1CaloCpmMonTool/L1BSCPMMonTool" ]
toolSvc.L1BSCPMMonTool.HistogramPrefix = "BS"
toolSvc.L1BSCPMMonTool.TriggerTowerLocation = "TriggerTowers"
toolSvc.L1BSCPMMonTool.CPMTowerLocation = "CPMTowers"
toolSvc.L1BSCPMMonTool.CPMHitsLocation = "CPMHits"
toolSvc.L1BSCPMMonTool.CMMCPHitsLocation = "CMMCPHits"
toolSvc.L1BSCPMMonTool.CPMRoILocation = "CPMRoIs"
#toolSvc.L1BSCPMMonTool.OutputLevel = DEBUG


# FileKey must match that given to THistSvc
#monMan.FileKey = "L1CaloMon"
monMan.FileKey = "stat"

# Set global monitoring parameters: see the AthenaMonManager class
# in the Control/AthenaMonitoring package
#monMan.ManualDataTypeSetup = True
#monMan.DataType            = "collisions"
#monMan.Environment         = "tier0"
#monMan.ManualRunLBSetup    = True
#monMan.Run                 = 1
#monMan.LumiBlock           = 1

####################################################################################################################
###################################### Ende ########################################################################
####################################################################################################################

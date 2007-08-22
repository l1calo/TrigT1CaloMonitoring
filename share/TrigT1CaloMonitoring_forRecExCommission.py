####################################################################################################################
##################################### Monitoring ###################################################################
####################################################################################################################

# ************************
# Monitoring configuration
# ************************
topSequence += AthenaMonManager( "L1CaloMonManager" )
L1CaloMan = topSequence.L1CaloMonManager

## get a handle on the ToolSvc
from AthenaCommon.AppMgr import ToolSvc as toolSvc

####################### TriggerTower ################################
L1CaloMan.AthenaMonTools += [ "TriggerTowerMon/L1TriggerTowerMonTool" ]
#toolSvc.L1TriggerTowerMonTool.DataType = "M3"  
toolSvc.L1TriggerTowerMonTool.BS_TriggerTowerContainer = "TriggerTowers"
toolSvc.L1TriggerTowerMonTool.DistPerChannel = True
#toolSvc.L1TriggerTowerMonTool.DistPerChannel = False
toolSvc.L1TriggerTowerMonTool.DistPerChannelAndTimeSlice = False
toolSvc.L1TriggerTowerMonTool.LUTHitMap_Thresh0 = 1
toolSvc.L1TriggerTowerMonTool.LUTHitMap_Thresh1 = 3
toolSvc.L1TriggerTowerMonTool.LUTHitMap_Thresh2 = 7
toolSvc.L1TriggerTowerMonTool.ADCHitMap_Thresh = 30
toolSvc.L1TriggerTowerMonTool.PathInRootFile = "L1Calo/PPr"
toolSvc.L1TriggerTowerMonTool.ErrorPathInRootFile = "L1Calo/Errors/PPr"
#toolSvc.L1TriggerTowerMonTool.OutputLevel = DEBUG


if ATLASCosmicFlags.doLVL1CaloCPJEPMonitoring:
    include ("TrigT1CaloMonitoring/TrigT1CaloMonitoring_L1CaloSimulation.py")

####################### JetElements ################################
    L1CaloMan.AthenaMonTools += [ "JetElementMon/BS_L1JetElementMonTool" ]
    toolSvc.BS_L1JetElementMonTool.DataType = "BS"  #BS or Sim data?
    toolSvc.BS_L1JetElementMonTool.JetElementLocation = "JetElements"
    toolSvc.BS_L1JetElementMonTool.PathInRootFile = "L1Calo/JEP/JEM/BS/input"
    toolSvc.BS_L1JetElementMonTool.ErrorPathInRootFile = "L1Calo/Errors"
    #toolSvc.BS_L1JetElementMonTool.OutputLevel = DEBUG

    L1CaloMan.AthenaMonTools += [ "JetElementMon/Sim_L1JetElementMonTool" ]
    toolSvc.Sim_L1JetElementMonTool.DataType = "Sim"  #BS or Sim data?
    toolSvc.Sim_L1JetElementMonTool.JetElementLocation = "Sim_JetElements"
    toolSvc.Sim_L1JetElementMonTool.PathInRootFile = "L1Calo/JEP/JEM/Sim/input"
    #toolSvc.Sim_L1JetElementMonTool.OutputLevel = DEBUG

####################### JEMs ################################
    L1CaloMan.AthenaMonTools += [ "JEMMon/BS_L1JEMMonTool" ]
    toolSvc.BS_L1JEMMonTool.DataType = "BS"  #BS or Sim data?
    toolSvc.BS_L1JEMMonTool.JEMHitsLocation = "JEMHits"
    toolSvc.BS_L1JEMMonTool.JEMEtSumsLocation = "JEMEtSums"
    toolSvc.BS_L1JEMMonTool.JEMRoILocation = "JEMRoIs"
    toolSvc.BS_L1JEMMonTool.PathInRootFile = "L1Calo/JEP/JEM/BS"
    toolSvc.BS_L1JEMMonTool.ErrorPathInRootFile = "L1Calo/Errors/JEM"
    #toolSvc.BS_L1JEMMonTool.OutputLevel = DEBUG

    L1CaloMan.AthenaMonTools += [ "JEMMon/Sim_L1JEMMonTool" ]
    toolSvc.Sim_L1JEMMonTool.DataType = "Sim"  #BS or Sim data?
    toolSvc.Sim_L1JEMMonTool.JEMHitsLocation = "Sim_JEMHits"
    toolSvc.Sim_L1JEMMonTool.JEMEtSumsLocation = "Sim_JEMEtSums"
    toolSvc.Sim_L1JEMMonTool.JEMRoILocation = "Sim_JEMRoIs"
    toolSvc.Sim_L1JEMMonTool.PathInRootFile = "L1Calo/JEP/JEM/Sim"
    #toolSvc.Sim_L1JEMMonTool.OutputLevel = DEBUG

####################### CMMs ################################
    L1CaloMan.AthenaMonTools += [ "CMMMon/BS_L1CMMMonTool" ]
    toolSvc.BS_L1CMMMonTool.DataType = "BS"  #BS or Sim data?
    toolSvc.BS_L1CMMMonTool.CMMJetHitsLocation = "CMMJetHits"
    toolSvc.BS_L1CMMMonTool.CMMEtSumsLocation = "CMMEtSums"
    toolSvc.BS_L1CMMMonTool.CMMRoILocation = "CMMRoIs"
    toolSvc.BS_L1CMMMonTool.JEMHitsLocation = "JEMHits"
    toolSvc.BS_L1CMMMonTool.JEMEtSumsLocation = "JEMEtSums"
    toolSvc.BS_L1CMMMonTool.PathInRootFile = "L1Calo/JEP/CMM/BS"
    toolSvc.BS_L1CMMMonTool.ErrorPathInRootFile = "L1Calo/Errors/CMM"
    #toolSvc.BS_L1CMMMonTool.OutputLevel = DEBUG

    L1CaloMan.AthenaMonTools += [ "CMMMon/Sim_L1CMMMonTool" ]
    toolSvc.Sim_L1CMMMonTool.DataType = "Sim"  #BS or Sim data?
    toolSvc.Sim_L1CMMMonTool.CMMJetHitsLocation = "Sim_CMMJetHits"
    toolSvc.Sim_L1CMMMonTool.CMMEtSumsLocation = "Sim_CMMEtSums"
    toolSvc.Sim_L1CMMMonTool.CMMRoILocation = "Sim_CMMRoIs"
    toolSvc.Sim_L1CMMMonTool.JEMHitsLocation = "Sim_JEMHits"
    toolSvc.Sim_L1CMMMonTool.JEMEtSumsLocation = "Sim_JEMEtSums"
    toolSvc.Sim_L1CMMMonTool.PathInRootFile = "L1Calo/JEP/CMM/Sim"
    #toolSvc.Sim_L1CMMMonTool.OutputLevel = DEBUG


####################### CPMs ################################
    L1CaloMan.AthenaMonTools += [ "TrigT1CaloCpmMonTool/L1BSCPMMonTool" ]
    toolSvc.L1BSCPMMonTool.HistogramPrefix = "BS"
    toolSvc.L1BSCPMMonTool.TriggerTowerLocation = "TriggerTowers"
    toolSvc.L1BSCPMMonTool.CPMTowerLocation = "CPMTowers"
    toolSvc.L1BSCPMMonTool.CPMHitsLocation = "CPMHits"
    toolSvc.L1BSCPMMonTool.CMMCPHitsLocation = "CMMCPHits"
    toolSvc.L1BSCPMMonTool.CPMRoILocation = "CPMRoIs"
    #toolSvc.L1BSCPMMonTool.OutputLevel = DEBUG

    L1CaloMan.AthenaMonTools += [ "TrigT1CaloCpmMonTool/L1SimCPMMonTool" ]
    toolSvc.L1SimCPMMonTool.HistogramPrefix = "Sim"
    toolSvc.L1SimCPMMonTool.TriggerTowerLocation = "TriggerTowers"
    toolSvc.L1SimCPMMonTool.CPMTowerLocation = "Sim_CPMTowers"
    toolSvc.L1SimCPMMonTool.CPMHitsLocation = "Sim_CPMHits"
    toolSvc.L1SimCPMMonTool.CMMCPHitsLocation = "Sim_CMMCPHits"
    toolSvc.L1SimCPMMonTool.CPMRoILocation = "Sim_CPMRoIs"
    #toolSvc.L1SimCPMMonTool.OutputLevel = DEBUG

####################### Simulation vs HW ################################
    L1CaloMan.AthenaMonTools += [ "SimBSMon/SimBSMonTool" ]
    toolSvc.SimBSMonTool.BS_JEMEtSumsLocation = "JEMEtSums"
    toolSvc.SimBSMonTool.Sim_JEMEtSumsLocation = "Sim_JEMEtSums"
    toolSvc.SimBSMonTool.BS_JEMHitsLocation = "JEMHits"
    toolSvc.SimBSMonTool.Sim_JEMHitsLocation = "Sim_JEMHits"
    toolSvc.SimBSMonTool.BS_CMMEtSumsLocation = "CMMEtSums"
    toolSvc.SimBSMonTool.Sim_CMMEtSumsLocation = "Sim_CMMEtSums"
    toolSvc.SimBSMonTool.BS_CMMJetHitsLocation = "CMMJetHits"
    toolSvc.SimBSMonTool.Sim_CMMJetHitsLocation = "Sim_CMMJetHits"
    toolSvc.SimBSMonTool.PathInRootFile = "L1Calo/SimBS"
    #toolSvc.SimBSMonTool.OutputLevel = DEBUG
    toolSvc.SimBSMonTool.OutputLevel = WARNING


####################### CaloTriggerTower ################################
if ATLASCosmicFlags.doLAr and  ATLASCosmicFlags.doTile:
    from TrigT1Calo.TrigT1CaloConf import LVL1__TriggerTowerMaker
    from AthenaCommon.AlgSequence import AlgSequence
    job = AlgSequence()
    job += LVL1__TriggerTowerMaker( 'TriggerTowerMaker' )
    job.TriggerTowerMaker.CellType = 1
    job.TriggerTowerMaker.TowerNoise = FALSE
    #input (from BS)
    job.TriggerTowerMaker.CaloCellLocation = "AllCalo"
    #output
    job.TriggerTowerMaker.TriggerTowerLocation = "Calo_TriggerTowers"
    #job.TriggerTowerMaker.OutputLevel = DEBUG

    L1CaloMan.AthenaMonTools += [ "CaloTTMon/L1CaloTTMonTool" ]
    toolSvc.L1CaloTTMonTool.CaloTTContainer = "Calo_TriggerTowers"
    toolSvc.L1CaloTTMonTool.Calo_HitMap_Thresh0 = 1
    toolSvc.L1CaloTTMonTool.Calo_HitMap_Thresh1 = 3
    toolSvc.L1CaloTTMonTool.Calo_HitMap_Thresh2 = 7
    toolSvc.L1CaloTTMonTool.PathInRootFile = "L1Calo/CaloTT"
    #toolSvc.L1CaloTTMonTool.OutputLevel = DEBUG


# FileKey must match that given to THistSvc
L1CaloMan.FileKey = "GLOBAL"

####################################################################################################################
###################################### Ende ########################################################################
####################################################################################################################

Offline= not ATLASCosmicFlags.doOnline
CompareWithSimulation=True
MaxEnergyRange = 50





from TrigT1Calo.TrigT1CaloConf import LVL1__JetElementMaker
from AthenaCommon.AlgSequence import AlgSequence
myjob = AlgSequence()
myjob += LVL1__JetElementMaker( 'JetElementMaker' )

# JetElementMaker
#input (from BS)
myjob.JetElementMaker.TriggerTowerLocation ="TriggerTowers"  
#output
myjob.JetElementMaker.JetElementLocation ="Sim_JetElements"
#job.JetElementMaker.OutputLevel = VERBOSE


if CompareWithSimulation:
    include ("TrigT1CaloMonitoring/TrigT1CaloMonitoring_L1CaloSimulation.py")

# ************************
# Monitoring configuration
# ************************
topSequence += AthenaMonManager( "L1CaloMonManager" )
L1CaloMan = topSequence.L1CaloMonManager

## get a handle on the ToolSvc
from AthenaCommon.AppMgr import ToolSvc as toolSvc

#=================================================================================
#================================= PPr ===========================================
#=================================================================================

#--------------------------------------------------------------
# CaloTriggerTowerService
#--------------------------------------------------------------
#required for use of the CaloTriggerTowerService
#include("CaloConditions/CaloConditions_jobOptions.py")

# TT-cell map (maintained by F. Ledroit)
include( "CaloConditions/LArTTCellMap_ATLAS_jobOptions.py")


###################################################################
#
#  Job options file for initializing CaloTriggerTowerService's
#  This jo is called by CaloConditions_jobOptions.py which is the
#  main jo to be used to initialize the CaloTriggerTowerService
#
#  Created 2007-APR-25 by D. Prieur
#==================================================================

#------------------------------------------------------------------
# Add the TTOnOffIdMap & TTOnAttrIdMap Maps COOL folders to IOVDbSvc
#------------------------------------------------------------------du

# Set the connection string
include ( "IOVDbSvc/IOVDbSvc_jobOptions.py" )
IOVDbSvc = Service( "IOVDbSvc" )

# DB connection string
DBConnection = CondDBCool.CALO

# add the DB folders to be opened
FolderNameList  = [ "/CALO/Identifier/CaloTTOnOffIdMapAtlas" ]
FolderNameList += [ "/CALO/Identifier/CaloTTOnAttrIdMapAtlas" ]

# tag to be used
FolderTag = "-0001" # HAD FCAL fix

#------------------------------------------------------------------

for FolderName in FolderNameList:
    if ( FolderTag == "" ):
        TagSpec = "HEAD"
    else:
        from string import *
        TagSpec = join(split(FolderName,'/'),'') + FolderTag

    IOVDbSvc.Folders+=[ FolderName + DBConnection + "<tag>"+TagSpec+"</tag>" ]



include("LArDetDescr/LArDetDescr_joboptions.py")
#include( "LArConditionsCommon/LArIdMap_comm_jobOptions.py" )


#--------------------------------- PPM -------------------------------------------
L1CaloMan.AthenaMonTools += [ "PPrMon/L1PPrMonTool" ]
#toolSvc.L1PPrMonTool.DataType = "M3"  
toolSvc.L1PPrMonTool.BS_TriggerTowerContainer = "TriggerTowers"
#toolSvc.L1PPrMonTool.DistPerChannel = True
toolSvc.L1PPrMonTool.DistPerChannel = False
toolSvc.L1PPrMonTool.DistPerChannelAndTimeSlice = False
toolSvc.L1PPrMonTool.LUTHitMap_Thresh0 = 1
toolSvc.L1PPrMonTool.LUTHitMap_Thresh1 = 3
toolSvc.L1PPrMonTool.LUTHitMap_Thresh2 = 7
toolSvc.L1PPrMonTool.ADCHitMap_Thresh = 30
toolSvc.L1PPrMonTool.MaxEnergyRange = MaxEnergyRange
toolSvc.L1PPrMonTool.PathInRootFile = "L1Calo/1_PPr"
toolSvc.L1PPrMonTool.ErrorPathInRootFile = "L1Calo/01_Errors_PPr"
toolSvc.L1PPrMonTool.EventPathInRootFile = "L1Calo"
toolSvc.L1PPrMonTool.Offline = Offline
#toolSvc.L1PPrMonTool.OutputLevel = DEBUG

#---------------------------- Performance Checks -----------------------------------

#=================================================================================
#=================================== JEP =========================================
#=================================================================================

#------------------------------------ JEM ----------------------------------------
L1CaloMan.AthenaMonTools += [ "JEMMon/BS_L1JEMMonTool" ]
toolSvc.BS_L1JEMMonTool.DataType = "BS"  #BS or Sim data?
toolSvc.BS_L1JEMMonTool.JetElementLocation = "JetElements"
toolSvc.BS_L1JEMMonTool.JEMHitsLocation = "JEMHits"
toolSvc.BS_L1JEMMonTool.JEMEtSumsLocation = "JEMEtSums"
toolSvc.BS_L1JEMMonTool.JEMRoILocation = "JEMRoIs"
toolSvc.BS_L1JEMMonTool.MaxEnergyRange = MaxEnergyRange
toolSvc.BS_L1JEMMonTool.PathInRootFile = "L1Calo/2_JEP_JEM"
toolSvc.BS_L1JEMMonTool.ErrorPathInRootFile = "L1Calo/02_Errors_JEM"
toolSvc.BS_L1JEMMonTool.Offline = Offline
#toolSvc.BS_L1JEMMonTool.OutputLevel = DEBUG

#if CompareWithSimulation:
#    L1CaloMan.AthenaMonTools += [ "JEMMon/Sim_L1JEMMonTool" ]
#    toolSvc.Sim_L1JEMMonTool.DataType = "Sim"  #BS or Sim data?
#    toolSvc.Sim_L1JEMMonTool.JetElementLocation = "Sim_JetElements"
#    toolSvc.Sim_L1JEMMonTool.JEMHitsLocation = "Sim_JEMHits"
#    toolSvc.Sim_L1JEMMonTool.JEMEtSumsLocation = "Sim_JEMEtSums"
#    toolSvc.Sim_L1JEMMonTool.JEMRoILocation = "Sim_JEMRoIs"
#    toolSvc.Sim_L1JEMMonTool.MaxEnergyRange = MaxEnergyRange
#    toolSvc.Sim_L1JEMMonTool.PathInRootFile = "L1Calo/Sim/2_JEP_JEM"
#toolSvc.Sim_L1JEMMonTool.OutputLevel = DEBUG

#----------------------------------- CMM ------------------------------------------
L1CaloMan.AthenaMonTools += [ "CMMMon/BS_L1CMMMonTool" ]
toolSvc.BS_L1CMMMonTool.DataType = "BS"  #BS or Sim data?
toolSvc.BS_L1CMMMonTool.CMMJetHitsLocation = "CMMJetHits"
toolSvc.BS_L1CMMMonTool.CMMEtSumsLocation = "CMMEtSums"
toolSvc.BS_L1CMMMonTool.CMMRoILocation = "CMMRoIs"
toolSvc.BS_L1CMMMonTool.MaxEnergyRange = MaxEnergyRange
toolSvc.BS_L1CMMMonTool.PathInRootFile = "L1Calo/3_JEP_CMM"
toolSvc.BS_L1CMMMonTool.ErrorPathInRootFile = "L1Calo/03_Errors_CMM"
toolSvc.BS_L1CMMMonTool.Offline = Offline
#toolSvc.BS_L1CMMMonTool.OutputLevel = DEBUG

#if CompareWithSimulation:
#    L1CaloMan.AthenaMonTools += [ "CMMMon/Sim_L1CMMMonTool" ]
#    toolSvc.Sim_L1CMMMonTool.DataType = "Sim"  #BS or Sim data?
#    toolSvc.Sim_L1CMMMonTool.CMMJetHitsLocation = "Sim_CMMJetHits"
#    toolSvc.Sim_L1CMMMonTool.CMMEtSumsLocation = "Sim_CMMEtSums"
#   toolSvc.Sim_L1CMMMonTool.CMMRoILocation = "Sim_CMMRoIs"
#   toolSvc.Sim_L1CMMMonTool.MaxEnergyRange = MaxEnergyRange
#   toolSvc.Sim_L1CMMMonTool.PathInRootFile = "L1Calo/Sim/3_JEP_CMM"
#   toolSvc.Sim_L1CMMMonTool.OutputLevel = DEBUG

#--------------------- Transmission and Performance ------------------------------

L1CaloMan.AthenaMonTools += [ "JEPTransPerfMon/JEPTransPerfMonTool" ]
toolSvc.JEPTransPerfMonTool.BS_JetElementLocation = "JetElements"
#toolSvc.JEPTransPerfMonTool.BS_TriggerTowerLocation = "TriggerTowers"
toolSvc.JEPTransPerfMonTool.BS_TriggerTowerLocation = "Sim_JetElements"
toolSvc.JEPTransPerfMonTool.NoLUTSlices=1

toolSvc.JEPTransPerfMonTool.BS_JEMHitsLocation = "JEMHits"
toolSvc.JEPTransPerfMonTool.Sim_JEMHitsLocation = "Sim_JEMHits"
toolSvc.JEPTransPerfMonTool.BS_JEMEtSumsLocation = "JEMEtSums"
toolSvc.JEPTransPerfMonTool.Sim_JEMEtSumsLocation = "Sim_JEMEtSums"
toolSvc.JEPTransPerfMonTool.BS_JEMRoILocation = "JEMRoIs"
toolSvc.JEPTransPerfMonTool.Sim_JEMRoILocation = "Sim_JEMRoIs"

toolSvc.JEPTransPerfMonTool.BS_CMMJetHitsLocation = "CMMJetHits"
toolSvc.JEPTransPerfMonTool.Sim_CMMJetHitsLocation = "Sim_CMMJetHits"
toolSvc.JEPTransPerfMonTool.BS_CMMEtSumsLocation = "CMMEtSums"
toolSvc.JEPTransPerfMonTool.Sim_CMMEtSumsLocation = "Sim_CMMEtSums"
toolSvc.JEPTransPerfMonTool.BS_CMMRoILocation = "CMMRoIs"
toolSvc.JEPTransPerfMonTool.Sim_CMMRoILocation = "Sim_CMMRoIs"

toolSvc.JEPTransPerfMonTool.PathInRootFile = "L1Calo/3_JEP_TransmissionAndPerformance"
toolSvc.JEPTransPerfMonTool.Offline = Offline
toolSvc.JEPTransPerfMonTool.CompareWithSimulation = CompareWithSimulation

#toolSvc.JEPTransPerfMonTool.OutputLevel = VERBOSE
#toolSvc.JEPTransPerfMonTool.OutputLevel = WARNING

#=================================================================================
#===================================== CP ========================================
#=================================================================================
L1CaloMan.AthenaMonTools += [ "TrigT1CaloCpmMonTool/L1BSCPMMonTool" ]
toolSvc.L1BSCPMMonTool.TriggerTowerLocation = "TriggerTowers"
toolSvc.L1BSCPMMonTool.CPMTowerLocation = "CPMTowers"
toolSvc.L1BSCPMMonTool.CPMHitsLocation = "CPMHits"
toolSvc.L1BSCPMMonTool.CMMCPHitsLocation = "CMMCPHits"
toolSvc.L1BSCPMMonTool.CPMRoILocation = "CPMRoIs"
toolSvc.L1BSCPMMonTool.RootDirectory = "L1Calo"
toolSvc.L1BSCPMMonTool.MaxEnergyRange = MaxEnergyRange
toolSvc.L1BSCPMMonTool.SingleDirectory = False
toolSvc.L1BSCPMMonTool.Offline = Offline
#toolSvc.L1BSCPMMonTool.OutputLevel = DEBUG

#if CompareWithSimulation:
#    L1CaloMan.AthenaMonTools += [ "TrigT1CaloCpmMonTool/L1SimCPMMonTool" ]
#    toolSvc.L1SimCPMMonTool.TriggerTowerLocation = "TriggerTowers"
#    toolSvc.L1SimCPMMonTool.CPMTowerLocation = "Sim_CPMTowers"
#    toolSvc.L1SimCPMMonTool.CPMHitsLocation = "Sim_CPMHits"
#    toolSvc.L1SimCPMMonTool.CMMCPHitsLocation = "Sim_CMMCPHits"
#    toolSvc.L1SimCPMMonTool.CPMRoILocation = "Sim_CPMRoIs"
#    toolSvc.L1SimCPMMonTool.RootDirectory = "L1Calo/Sim"
#    toolSvc.L1SimCPMMonTool.MaxEnergyRange = MaxEnergyRange
#    toolSvc.L1SimCPMMonTool.SingleDirectory = False
#toolSvc.L1SimCPMMonTool.OutputLevel = DEBUG





# FileKey must match that given to THistSvc
L1CaloMan.FileKey = "GLOBAL"

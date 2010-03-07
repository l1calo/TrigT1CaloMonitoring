Offline= not athenaCommonFlags.isOnline
if not 'DQMonFlags' in dir():
    print "TrigT1CaloMonitoring_forRecExCommission.py: DQMonFlags not yet imported - I import them now"
    from AthenaMonitoring.DQMonFlags import DQMonFlags

CompareWithSimulation=False
if (globalflags.DataSource() == "data") :
    CompareWithSimulation=True
    
# On Tier0 select monitoring tools according to processing step
if DQMonFlags.monManEnvironment() == 'tier0Raw':
    # Tier0 RAWtoESD step
    l1caloRawMon=True
    l1caloESDMon=False
elif DQMonFlags.monManEnvironment() == 'tier0ESD':
    # Tier0 ESDtoAOD step
    l1caloRawMon=False
    l1caloESDMon=True
else:
    # Anything else
    l1caloRawMon=True
    l1caloESDMon=True


#================================= Monitoring configuration ======================
from AthenaCommon.AlgSequence import AlgSequence
topSequence = AlgSequence()
topSequence += AthenaMonManager( "L1CaloMonManager" )
L1CaloMan = topSequence.L1CaloMonManager

## get a handle on the ToolSvc
from AthenaCommon.AppMgr import ToolSvc

if globalflags.InputFormat() == "bytestream":
    include ("TrigT1CaloByteStream/ReadLVL1CaloBS_jobOptions.py")

if l1caloESDMon:

    #=================================================================================
    #================================= PPr ===========================================
    #=================================================================================
    include("CaloConditions/CaloConditions_jobOptions.py")
    include("LArDetDescr/LArDetDescr_joboptions.py")

    #--------------------------------- PPM -------------------------------------------
    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import PPrMon
    L1PPrMonTool = PPrMon(
        name = "L1PPrMonTool",
        BS_TriggerTowerContainer = "TriggerTowers",
        #DistPerChannel = True,
        DistPerChannel = False,
        DistPerChannelAndTimeSlice = False,
        LUTHitMap_ThreshMax = 10,
        LUTHitMap_LumiBlocks = 10,
        ADCHitMap_Thresh = 50,
        MaxEnergyRange = 256,
        ADCTimingPerChannel = False,
        EMFADCCut = 8,
        HADFADCCut = 8,
        ADCPedestal = 35,
        PathInRootFile = "L1Calo/PPM",
        ErrorPathInRootFile = "L1Calo/PPM/Errors",
        EventPathInRootFile = "L1Calo/Overview",
        Offline = Offline,
        #OutputLevel = VERBOSE,
        #OutputLevel = INFO,
        )
    ToolSvc += L1PPrMonTool
    L1CaloMan.AthenaMonTools += [ L1PPrMonTool ]

    include("TrigT1CaloCalibConditions/L1CaloCalibConditions_jobOptions.py")
    
    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import PPMSimBSMon
    PPMSimBSMonTool = PPMSimBSMon("PPMSimBSMonTool")
    ToolSvc += PPMSimBSMonTool
    L1CaloMan.AthenaMonTools += [ PPMSimBSMonTool ]
    #ToolSvc.PPMSimBSMonTool.OutputLevel = DEBUG
    from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1TriggerTowerTool
    L1TriggerTowerTool = LVL1__L1TriggerTowerTool("L1TriggerTowerTool")
    ToolSvc += L1TriggerTowerTool
    #ToolSvc.L1TriggerTowerTool.OutputLevel = DEBUG
    
if l1caloRawMon:

    #--------------------------------- PPM Spare Channels----------------------------
    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import PPrSpareMon
    L1PPrSpareMonTool = PPrSpareMon(
        name = "L1PPrSpareMonTool",
        BS_TriggerTowerContainer = "TriggerTowersSpare",
        ADCHitMap_Thresh = 40,
        PathInRootFile = "L1Calo/PPM/SpareChannels",
        ErrorPathInRootFile = "L1Calo/PPM/SpareChannels/Errors",
        #OutputLevel = VERBOSE,
        #OutputLevel = INFO,
        )
    ToolSvc += L1PPrSpareMonTool
    L1CaloMan.AthenaMonTools += [ L1PPrSpareMonTool ]

    #---------------------------- Performance Checks -----------------------------------

if l1caloESDMon:

    #=================================================================================
    #=================================== JEP =========================================
    #=================================================================================

    #------------------------------------ JEM ----------------------------------------
    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import JEMMon
    BS_L1JEMMonTool = JEMMon(
        name = "BS_L1JEMMonTool",
        TypeOfData = "BS",  #BS or Sim data?
        JetElementLocation = "JetElements",
        JEMHitsLocation = "JEMHits",
        JEMEtSumsLocation = "JEMEtSums",
        JEMRoILocation = "JEMRoIs",
        MaxEnergyRange = 1024,
        PathInRootFile = "L1Calo/JEM",
        ErrorPathInRootFile = "L1Calo/JEM/Errors/Hardware",
        Offline = Offline,
        #OutputLevel = VERBOSE,
        )
    ToolSvc += BS_L1JEMMonTool
    L1CaloMan.AthenaMonTools += [ BS_L1JEMMonTool ]

    #----------------------------------- CMM ------------------------------------------
    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import CMMMon
    BS_L1CMMMonTool = CMMMon (
        name = "BS_L1CMMMonTool",
        TypeOfData = "BS",  #BS or Sim data?
        CMMJetHitsLocation = "CMMJetHits",
        CMMEtSumsLocation = "CMMEtSums",
        CMMRoILocation = "CMMRoIs",
        MaxEnergyRange = 32767,
        PathInRootFile = "L1Calo/JEM_CMM",
        ErrorPathInRootFile = "L1Calo/JEM_CMM/Errors/Hardware",
        Offline = Offline,
        #OutputLevel = VERBOSE,
        )
    ToolSvc += BS_L1CMMMonTool
    L1CaloMan.AthenaMonTools += [ BS_L1CMMMonTool ]

if l1caloRawMon and globalflags.DataSource() == "data":

    #--------------------- Transmission and Performance ------------------------------
    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import JEPSimBSMon
    JEPSimBSMonTool = JEPSimBSMon("JEPSimBSMonTool",
        JEPHitsTool = "LVL1::L1JEPHitsTools/L1JEPHitsTools_Mon",
        JetTool = "LVL1::L1JetTools/L1JetTools_Mon",
        JEPEtSumsTool = "LVL1::L1JEPEtSumsTools/L1JEPEtSumsTools_Mon",
        CompareWithSimulation = CompareWithSimulation)
    ToolSvc += JEPSimBSMonTool
    L1CaloMan.AthenaMonTools += [ JEPSimBSMonTool ]
    #ToolSvc.JEPSimBSMonTool.OutputLevel = DEBUG

    from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1JEPHitsTools
    L1JEPHitsTools = LVL1__L1JEPHitsTools("L1JEPHitsTools_Mon")
    L1JEPHitsTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
    ToolSvc += L1JEPHitsTools
    from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1JetTools
    L1JetTools = LVL1__L1JetTools("L1JetTools_Mon")
    L1JetTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
    ToolSvc += L1JetTools
    from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1EtTools
    L1EtTools = LVL1__L1EtTools("L1EtTools_Mon")
    L1EtTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
    ToolSvc += L1EtTools
    from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1JEPEtSumsTools
    L1JEPEtSumsTools = LVL1__L1JEPEtSumsTools("L1JEPEtSumsTools_Mon",
                                        EtTool = "LVL1::L1EtTools/L1EtTools_Mon")
    L1JEPEtSumsTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
    ToolSvc += L1JEPEtSumsTools

if l1caloESDMon:

    #=================================================================================
    #===================================== CP ========================================
    #=================================================================================
    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import TrigT1CaloCpmMonTool
    L1BSCPMMonTool = TrigT1CaloCpmMonTool (
        name = "L1BSCPMMonTool",
        TriggerTowerLocation = "TriggerTowers",
        CPMTowerLocation = "CPMTowers",
        CPMHitsLocation = "CPMHits",
        CMMCPHitsLocation = "CMMCPHits",
        CPMRoILocation = "CPMRoIs",
        RootDirectory = "L1Calo",
        MaxEnergyRange = 256,
        #OutputLevel = DEBUG,
        )
    ToolSvc += L1BSCPMMonTool
    L1CaloMan.AthenaMonTools += [ L1BSCPMMonTool ]

if l1caloRawMon and globalflags.DataSource() == "data":

    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import CPMSimBSMon
    CPMSimBSMonTool = CPMSimBSMon("CPMSimBSMonTool",
                              EmTauTool = "LVL1::L1EmTauTools/L1EmTauTools_Mon",
                              CompareWithSimulation = CompareWithSimulation)
    ToolSvc += CPMSimBSMonTool
    L1CaloMan.AthenaMonTools += [ CPMSimBSMonTool ]
    #ToolSvc.CPMSimBSMonTool.OutputLevel = DEBUG

    from TrigT1CaloTools.TrigT1CaloToolsConf import LVL1__L1EmTauTools
    L1EmTauTools = LVL1__L1EmTauTools("L1EmTauTools_Mon")
    L1EmTauTools.LVL1ConfigSvc="TrigConf::TrigConfigSvc/TrigConfigSvc"
    ToolSvc += L1EmTauTools

    #=================================================================================
    #===================================== ROD =======================================
    #=================================================================================
    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import TrigT1CaloRodMonTool
    L1BSRODMonTool = TrigT1CaloRodMonTool (
        name = "L1BSRODMonTool",
        #OutputLevel = DEBUG,
        )
    ToolSvc += L1BSRODMonTool
    L1CaloMan.AthenaMonTools += [ L1BSRODMonTool ]

if globalflags.DataSource() == "data":

    #=================================================================================
    #=============================== Global Overview =================================
    #=================================================================================
    from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import TrigT1CaloGlobalMonTool
    if l1caloRawMon:
        L1GlobalMonTool = TrigT1CaloGlobalMonTool ( name = "L1GlobalMonTool" )
    else:
        L1GlobalMonTool = TrigT1CaloGlobalMonTool ( name = "L1GlobalESDMonTool" )
    ToolSvc += L1GlobalMonTool
    L1CaloMan.AthenaMonTools += [ L1GlobalMonTool ]

#=================================================================================
# FileKey must match that given to THistSvc
L1CaloMan.FileKey             = DQMonFlags.monManFileKey()
L1CaloMan.Environment         = DQMonFlags.monManEnvironment()
L1CaloMan.ManualDataTypeSetup = DQMonFlags.monManManualDataTypeSetup()
L1CaloMan.DataType            = DQMonFlags.monManDataType()

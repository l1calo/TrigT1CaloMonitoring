Offline= not athenaCommonFlags.isOnline
if not 'DQMonFlags' in dir():
    print "TrigT1CaloMonitoring_forRecExCommission.py: DQMonFlags not yet imported - I import them now"
    from AthenaMonitoring.DQMonFlags import DQMonFlags
    
# On Tier0 select monitoring tools according to processing step
if DQMonFlags.monManEnvironment() == 'tier0Raw':
    # Tier0 RAWtoESD step
    l1caloRawMon=True
    l1caloESDMon=True
elif DQMonFlags.monManEnvironment() == 'tier0ESD':
    # Tier0 ESDtoAOD step
    l1caloRawMon=False
    l1caloESDMon=False
else:
    # Anything else
    l1caloRawMon=True
    l1caloESDMon=True
    
if l1caloRawMon or l1caloESDMon:
    
    #================================= Monitoring configuration ======================
    from AthenaCommon.AlgSequence import AlgSequence
    topSequence = AlgSequence()
    L1CaloMan = AthenaMonManager( "L1CaloMonManager" )
    
    ## get a handle on the ToolSvc
    from AthenaCommon.AppMgr import ToolSvc
    
    if globalflags.InputFormat() == "bytestream":
        include ("TrigT1CaloByteStream/ReadLVL1CaloBS_jobOptions.py")
    
    include("CaloConditions/CaloConditions_jobOptions.py")
    if Offline:
        include("LArDetDescr/LArDetDescr_joboptions.py")
    else:
        include("LArConditionsCommon/LArIdMap_comm_jobOptions.py")

    if l1caloESDMon and globalflags.DataSource() == "data":
        include("TrigT1CaloCalibConditions/L1CaloCalibConditionsTier0_jobOptions.py")
    
    from TrigT1CaloMonitoringTools.LVL1CaloMonFlags import LVL1CaloMonFlags
    if LVL1CaloMonFlags.doPPrStabilityMon():
    
        #=================================================================================
        #  Want PPrStabilityMon to run alone
        #=================================================================================
    
        if l1caloESDMon:
    
            from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import PPrStabilityMon
            L1PPrStabilityMonTool = PPrStabilityMon(
                name = "L1PPrStabilityMonTool",
                BS_TriggerTowerContainer = "TriggerTowers",
                ppmADCMinValue = 60,
                lumiMax = 2000,
                fineTimeCut = 20,
                PathInRootFile = "L1Calo/PPrStabilityMon",
                #OutputLevel = DEBUG
                )
            ToolSvc += L1PPrStabilityMonTool
            L1CaloMan.AthenaMonTools += [ L1PPrStabilityMonTool ]
    
    else:
    
        if l1caloESDMon:
    
            #=================================================================================
            #================================= PPr ===========================================
            #=================================================================================
            from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import PPrMon
            L1PPrMonTool = PPrMon(
                name = "L1PPrMonTool",
                BS_TriggerTowerContainer = "TriggerTowers",
                LUTHitMap_ThreshVec=[0,1,2,3,4,5,6,7,10,15,20,33,45,50],
                LUTHitMap_LumiBlocks = 10,
                ADCHitMap_Thresh = 50,
                MaxEnergyRange = 256,
                EMFADCCut = 40,
                HADFADCCut = 40,
                ADCPedestal = 32,
                PathInRootFile = "L1Calo/PPM",
                ErrorPathInRootFile = "L1Calo/PPM/Errors",
                EventPathInRootFile = "L1Calo/Overview",
                #OutputLevel = DEBUG
                )
            ToolSvc += L1PPrMonTool
            L1CaloMan.AthenaMonTools += [ L1PPrMonTool ]
    
        if l1caloESDMon and globalflags.DataSource() == "data":
                
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
                #OutputLevel = DEBUG
                )
            ToolSvc += L1PPrSpareMonTool
            L1CaloMan.AthenaMonTools += [ L1PPrSpareMonTool ]
    
        if l1caloESDMon:
    
            #=================================================================================
            #=================================== JEP =========================================
            #=================================================================================
    
            #------------------------------------ JEM ----------------------------------------
            from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import JEMMon
            L1JEMMonTool = JEMMon(
                name = "L1JEMMonTool",
                JetElementLocation = "JetElements",
                JEMHitsLocation = "JEMHits",
                JEMEtSumsLocation = "JEMEtSums",
                JEMRoILocation = "JEMRoIs",
                MaxEnergyRange = 1024,
                PathInRootFile = "L1Calo/JEM",
                ErrorPathInRootFile = "L1Calo/JEM/Errors/Hardware",
                #OutputLevel = DEBUG
                )
            ToolSvc += L1JEMMonTool
            L1CaloMan.AthenaMonTools += [ L1JEMMonTool ]
    
            #----------------------------------- CMM ------------------------------------------
            from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import CMMMon
            L1CMMMonTool = CMMMon (
                name = "L1CMMMonTool",
                CMMJetHitsLocation = "CMMJetHits",
                CMMEtSumsLocation = "CMMEtSums",
                CMMRoILocation = "CMMRoIs",
                PathInRootFile = "L1Calo/JEM_CMM",
                ErrorPathInRootFile = "L1Calo/JEM_CMM/Errors/Hardware",
                #OutputLevel = DEBUG
                )
            ToolSvc += L1CMMMonTool
            L1CaloMan.AthenaMonTools += [ L1CMMMonTool ]
    
        if l1caloRawMon and globalflags.DataSource() == "data":
    
            #--------------------- Transmission and Performance ------------------------------
            from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import JEPSimBSMon
            JEPSimBSMonTool = JEPSimBSMon("JEPSimBSMonTool",
                JEPHitsTool = "LVL1::L1JEPHitsTools/L1JEPHitsTools_Mon",
                JetTool = "LVL1::L1JetTools/L1JetTools_Mon",
                JEPEtSumsTool = "LVL1::L1JEPEtSumsTools/L1JEPEtSumsTools_Mon",
                )
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
                                      )
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
                prebookThresh = not l1caloESDMon
                L1GlobalMonTool = TrigT1CaloGlobalMonTool ( name = "L1GlobalMonTool",
    	                                                BookCPMThresh = prebookThresh,
    						        BookJEMThresh = prebookThresh,
    						        BookCMMThresh = prebookThresh,
    						        #OutputLevel = DEBUG
                                                          )
            else:
                L1GlobalMonTool = TrigT1CaloGlobalMonTool ( name = "L1GlobalESDMonTool",
    	                                                #OutputLevel = DEBUG
    						  )
            ToolSvc += L1GlobalMonTool
            L1CaloMan.AthenaMonTools += [ L1GlobalMonTool ]
    
        if l1caloRawMon and (globalflags.DataSource() == "data" and Offline
                             and rec.doCalo() and rec.doLArg() and rec.doTile()
                             and (rec.triggerStream() == "JetTauEtmiss"
                               or rec.triggerStream() == "Muons"
                               or rec.triggerStream() == "express")):
    
            #=================================================================================
            #=============================== EM Efficiencies =================================
            #=================================================================================
            trigstring = ['EF_.*']
            from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import EmEfficienciesMonTool
            L1EmEfficienciesMonTool = EmEfficienciesMonTool ( name = "EmEfficienciesMonTool",
                                                                  TriggerStrings = trigstring
                                                            )
            ToolSvc += L1EmEfficienciesMonTool
            L1CaloMan.AthenaMonTools += [ L1EmEfficienciesMonTool ]
            if not hasattr( ToolSvc, "TrigDecisionTool" ):
                from TrigDecisionTool.TrigDecisionToolConf import Trig__TrigDecisionTool
                tdt = Trig__TrigDecisionTool('TrigDecisionTool')
                ToolSvc += tdt
    
        if l1caloRawMon and (globalflags.DataSource() == "data" and Offline
                             and rec.doCalo() and rec.doLArg() and rec.doTile()
                             and (rec.triggerStream() == "Egamma"
                               or rec.triggerStream() == "Muons"
                               or rec.triggerStream() == "express")):
    
            #=================================================================================
            #=============================== Jet Efficiencies ================================
            #=================================================================================
            trigstring = ['EF_.*']
            from TrigT1CaloMonitoring.TrigT1CaloMonitoringConf import JetEfficienciesMonTool
            L1JetEfficienciesMonTool = JetEfficienciesMonTool ( name = "JetEfficienciesMonTool",
                                                                  TriggerStrings = trigstring
                                                              )
            ToolSvc += L1JetEfficienciesMonTool
            L1CaloMan.AthenaMonTools += [ L1JetEfficienciesMonTool ]
            if not hasattr( ToolSvc, "TrigDecisionTool" ):
                from TrigDecisionTool.TrigDecisionToolConf import Trig__TrigDecisionTool
                tdt = Trig__TrigDecisionTool('TrigDecisionTool')
                ToolSvc += tdt
    
    #=================================================================================
    # FileKey must match that given to THistSvc
    L1CaloMan.FileKey             = DQMonFlags.monManFileKey()
    L1CaloMan.Environment         = DQMonFlags.monManEnvironment()
    L1CaloMan.ManualDataTypeSetup = DQMonFlags.monManManualDataTypeSetup()
    L1CaloMan.DataType            = DQMonFlags.monManDataType()
    topSequence += L1CaloMan

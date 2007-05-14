####################################################################################################################
##################################### Konfiguration Simulation #####################################################
####################################################################################################################

AllAlgs = False

# Set flags that determine the input data source
# Set only one flag True:
readZEBRA=False
readPool=False
readBS=  False
readMuonHits=False   
readMuonDigit=False

# These flags can be set in addition:
transientBS=False 
writeBS= False     
doDumpTES = False
doDumpTDS = False

doWriteESD=False # uncomment if do not write ESD
doWriteAOD=FALSE # uncomment if do not write AOD
doWriteTAG=False # uncomment if do not write TAG
doAOD = False
doAODall = False
doHist = False
doTruth = False
useROOTNtuple = False

doTrigger = True
doLArCellDirect = False
doxKalman = False
doiPatRec = False
doEmCluster = False
doCaloCluster = False
doCaloTopoCluster = False
doMoore = False
doMuonboy = False
doConversion = False
doBtagging = False
doEgamma = False
doJetRec = False
doTauRec = False
doMuonIDStandAlone = False
doMuonIDCombined = False      
doStaco = False
doMuTag = False
doTileMuID = False
doMissingET = False
doEFlow = False
doEFlowJet = False
doAODLVL1 = False
doCBNT = False

if DoCaloMon:
    readPool=True
    writeBS= True

if DoBSMon:
    readBS=  True

# Set Flags for the Trigger
include ("TriggerRelease/TriggerFlags.py")

TriggerFlags.doLVL1=False  
TriggerFlags.doLVL2=False
TriggerFlags.doEF=False
TriggerFlags.doMuon=False

TriggerFlags.onlyCalo = True
TriggerFlags.doJetSlice = True
TriggerFlags.doLVL1 = True

doTracking = False # offline track reconstruction

## Set flags for the ntuple Content
TriggerFlags.doNtuple=False
useROOTNtuple=False # False : hbook ntuple

## Set Truth flags
doTruth=False    # needs to be true if creating Fake RoI 

TriggerFlags.readLVL1configFromXML = True
Service("L1Config").thresholdListFileLocation = "LVL1triggerthresholdsCSC-01.xml"
Service("L1Config").triggerMenuFileLocation = "LVL1triggermenuCSC-01.xml"

# DetFlags modifications are best set here (uncomment RecExCommon_flags first)
#include ("RecExCommon/RecExCommon_flags.py")
# switch off ID, calo, or muons
#DetFlags.ID_setOff()
# DetFlags.Calo_setOff()
#DetFlags.Muon_setOff()

if DoCaloMon:
    TriggerFlags.doLVL1=True  

if DoBSMon:
    TriggerFlags.doLVL1=False  


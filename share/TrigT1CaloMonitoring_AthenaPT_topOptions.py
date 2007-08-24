######################################################################
# TopOptions to execute RecExCommissioning in PT
#               (and provide data for online event display)
#
######################################################################

#default flags
include("RecExCommission/RecExCommissionFlags_jobOptions.py")

#important
ATLASCosmicFlags.doOnline=True
ATLASCosmicFlags.doMonitoring=True
ATLASCosmicFlags.doLVL1Calo=True
ATLASCosmicFlags.doLVL1CaloCPJEPMonitoring=True

CBNTAthenaAware = False
ATLASCosmicFlags.doSim=False
ATLASCosmicFlags.JiveXML=False
ATLASCosmicFlags.doInDet=False
ATLASCosmicFlags.doMuons=False
ATLASCosmicFlags.doLAr=False
ATLASCosmicFlags.doTile=False
ATLASCosmicFlags.doCaloTrkMuId=False
ATLASCosmicFlags.doHLT=False
ATLASCosmicFlags.doCTPMon=False
ATLASCosmicFlags.OnlineJiveXML = False
ATLASCosmicFlags.doNtuple = False


#
if ATLASCosmicFlags.doOnline:
    include.block("ByteStreamCnvSvc/BSEventStorageEventSelector_jobOptions.py")

# 
include( "RecExCommission/RecExCommission_topOptions.py" )


##

#MuonCalib_MdtCalibDbAsciiTool.RT_InputFiles = [ "RT_default_comm.dat" ]
#ServiceMgr.ExceptionSvc.Catch="none"
try:
    ByteStreamAddressProviderSvc.TypeNames.remove("LArCellIDC/LArCellIDC")
except ValueError:
    pass
    
MessageSvc.defaultLimit=10000
Service( "ROBDataProviderSvc" ).filterEmptyROB=True

#to avoid crash in SCT monit
#ByteStreamCnvSvc.InitCnvs += [ "SCT_RDO_Container","InDetRawDataCollection<SCT_RDORawData>"]

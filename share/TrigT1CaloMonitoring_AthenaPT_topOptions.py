######################################################################
# TopOptions to execute RecExCommissioning in PT
#               (and provide data for online event display)
#
# 21/08/2007  Jiri.Masik@cern.ch
######################################################################

if (not 'testRun' in dir()):
    testRun=False
    
include( "RecExCommission/RecExCommissionCommonFlags_jobOptions.py" )
rec.doHist = True
rec.doJiveXML=False
rec.doWriteESD=False
rec.doWriteTAG=False
rec.doPerfMon=False
#default flags
include("RecExCommission/RecExCommissionFlags_jobOptions.py")



#now we need also to do

#important
ATLASCosmicFlags.doSim=False
ATLASCosmicFlags.JiveXML=False

if not testRun:
    ATLASCosmicFlags.OnlineJiveXML=False
    
ATLASCosmicFlags.OnlineJiveXML=False
#ATLASCosmicFlags.OnlineJiveXML=False
#athenaCommonFlags.isOnline=True
CBNTAthenaAware = False
ATLASCosmicFlags.doESD=False
ATLASCosmicFlags.doNtuple=False
ATLASCosmicFlags.doFilteredESD=False
#this does not exist yet ut allows to set TRTDummy
ATLASCosmicFlags.CosmicSetup='ATLAS'

#optional switches
ATLASCosmicFlags.CosmicSetup = "M8"
#temporary
globalflags.DetDescrVersion = 'ATLAS-CommNF-07-00-00'
#ATLASCosmicFlags.IOVDbSvcGlobalTagData = 'COMCOND-004-00'
#ATLASCosmicFlags.IOVDbSvcGlobalTagData = 'COMCOND-004-01'
ATLASCosmicFlags.IOVDbSvcGlobalTagData = 'COMCOND-006-00'
ATLASCosmicFlags.doMonitoring=True
ATLASCosmicFlags.doCombinedFit=False # No TRT in M5
#ATLASCosmicFlags.doInDet=False
#InDetCosmicFlags.TRTCalibUseDB=True
#muons disabled for the TR
ATLASCosmicFlags.doInDetMon=False
ATLASCosmicFlags.doTileMon=False
ATLASCosmicFlags.doLArMon=False
ATLASCosmicFlags.doMuonMon=False
ATLASCosmicFlags.doHLTMon=False

rec.doMuons=False
rec.doLArg=False
rec.doTile=False
ATLASCosmicFlags.doLVL1Calo=True
ATLASCosmicFlags.doCaloTrkMuId=False
ATLASCosmicFlags.doHLT=False
rec.doInDet = False
ATLASCosmicFlags.doCTPMon=False
ATLASCosmicFlags.doCaloTrkMuId=False
#ATLASCosmicFlags.CaloDetailedJiveXML = True
ATLASCosmicFlags.CaloDetailedJiveXML = False
ATLASCosmicFlags.useLocalCOOL = False
rec.doDetStatus = False
ATLASCosmicFlags.doDQMonitoring= False
#afs

LArDBConnection=" <dbConnection>impl=cool;techno=sqlite;schema=/det/lar/project/databases/COMCOND200.db;X:COMP200:</dbConnection>"
ATLASCosmicFlags.CondDBCool=LArDBConnection

include("InDetCosmicRecExample/InDetCosmicFlags_jobOptions.py")
InDetCosmicFlags.TRTMCdummyMode=True
InDetCosmicFlags.doMonitoring=False
InDetCosmicFlags.doPixelMonitoring=False
InDetCosmicFlags.doSCTMonitoring=False
InDetCosmicFlags.doTRTMonitoring=False
InDetCosmicFlags.doGlobalMonitoring=False
InDetCosmicFlags.doAlignMonitoring=False
InDetCosmicFlags.maskTRT=False
#
#LArCoolChannelSelection="0,3:92,181:205,231,233,235,237"
#LArCoolChannelSelection="0,92:116,206:230,232,234,236,238,331,338,344,350,1003,1007,1009,1011"
LArCoolChannelSelection="3:238,306,313,319,325,331,338,344,350,1001:1012,1021,1022"

#
if athenaCommonFlags.isOnline:
    include.block("ByteStreamCnvSvc/BSEventStorageEventSelector_jobOptions.py")

# 
include( "RecExCommission/RecExCommission_topOptions.py" )

#include( "ByteStreamCnvSvcBase/BSAddProvSvc_RIO_jobOptions.py" )
#include( "ByteStreamCnvSvcBase/BSAddProvSvc_RDO_jobOptions.py" )

from AthenaCommon.AppMgr import ToolSvc

#
theApp.CreateSvc += [ "ByteStreamCnvSvc" ]
#svcMgr = theApp.serviceMgr()
#svcMgr += ByteStreamCnvSvc()
##

#ToolSvc.MuonCalib_MdtCalibDbAsciiTool.RT_InputFiles = [ "/preseries-scratch/masik/efatl/RT_default_comm.dat" ]
#ServiceMgr.ExceptionSvc.Catch="none"
#try:
 #   ByteStreamAddressProviderSvc.TypeNames.remove("LArCellIDC/LArCellIDC")
    #temporary avoid crash on this container
    #ByteStreamAddressProviderSvc.TypeNames.remove("Muon::TgcPrepDataContainer/TGC_Measurements")
#except ValueError:
 #   pass

try:
    PoolSvc = Service( "PoolSvc" )
    PoolSvc.ReadCatalog.remove("prfile:poolcond/PoolCat_comcond_castor.xml")
except ValueError:
    print 'PoolCat_comcond_castor.xml cannot be removed'

# try:
#     PoolSvc.ReadCatalog.remove("prfile:poolcond/PoolCat_comcond_castor.xml")
# except ValueError:
#     print 'PoolCat_comcond_castor.xml cannot be removed'


MessageSvc.defaultLimit=10000
Service( "ROBDataProviderSvc" ).filterEmptyROB=True
#ToolSvc.LArCablingService.OutputLevel=FATAL
#ToolSvc.calonoisetool.OutputLevel=DEBUG
#TopoMaker.UseCaloNoiseTool=False
#ToolSvc.EventData2XML.OutputLevel=ERROR
#ToolSvc.EventData2XML.OutputLevel=DEBUG

if testRun:
    ToolSvc.EventData2XML.WriteToFile=False

#to avoid crash in SCT monit
#ByteStreamCnvSvc.InitCnvs += [ "SCT_RDO_Container","InDetRawDataCollection<SCT_RDORawData>"]
#MessageSvc.OutputLevel=DEBUG
MessageSvc.OutputLevel = VERBOSE


if rec.doLArg:
    PoolSvc = Service( "PoolSvc" )
    PoolSvc.ReadCatalog +=["xmlcatalog_file:/det/lar/project/databases/oflcond.000003.conditions.simul.pool.v0000/PoolFileCatalog.xml"]
    PoolSvc.ReadCatalog +=["xmlcatalog_file:/det/lar/project/databases/comcond.000006.lar_conditions.recon.pool.v0000/PoolFileCatalog.xml"]   #AK
    PoolSvc.ReadCatalog +=["xmlcatalog_file:/det/lar/project/databases/comcond.000005.lar_conditions.recon.pool.v0000/PoolFileCatalog.xml"]   #AK
    PoolSvc.ReadCatalog +=["xmlcatalog_file:/det/lar/project/databases/comcond.000004.lar_conditions.recon.pool.v0000/PoolFileCatalog.xml"]   #AK
    PoolSvc.ReadCatalog +=["xmlcatalog_file:/det/lar/project/databases/comcond.000003.lar_conditions.recon.pool.v0000/PoolFileCatalog.xml"]   #AK
    PoolSvc.ReadCatalog +=["xmlcatalog_file:/det/lar/project/databases/comcond.000002.lar_conditions.recon.pool.v0000/PoolFileCatalog.xml"]   #AK
    
    #    PoolSvc.ReadCatalog +=["xmlcatalog_file:/det/lar/project/databases/intr130/PoolFileCatalog.xml"]   #AK
    theApp.Dlls += ["LArRawConditions","LArCondAthenaPoolPoolCnv","LArCalibTools","LArCalibUtils"]
    theApp.Dlls += ["MuonByteStream"]



#tmp for
if ATLASCosmicFlags.doLVL1Calo:
    include("JiveXML/DataTypes_Trig.py")
    from TrigJiveXML.TrigJiveXMLConf import JiveXML__CTPDecisionRetriever
    CTPDecRetriever = JiveXML__CTPDecisionRetriever (name = "CTPDecisionRetriever")
    CTPDecRetriever.readCTP = True
                    
Service("PoolSvc").SortReplicas = False

if ATLASCosmicFlags.doLArg and ATLASCosmicFlags.CaloDetailedJiveXML:
    theCaloRetriever = JiveXML__CaloRetriever (name = "CaloRetriever")
    theCaloRetriever.CellThreshold = 150


#print 'jmasik'
if ATLASCosmicFlags.JiveXML:
    print ToolSvc.EventData2XML.DataTypes
    




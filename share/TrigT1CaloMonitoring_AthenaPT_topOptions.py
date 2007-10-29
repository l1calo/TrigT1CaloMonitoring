######################################################################
# TopOptions to execute RecExCommissioning in PT
#               (and provide data for online event display)
#
# 21/08/2007  Jiri.Masik@cern.ch
######################################################################

#default flags
include("RecExCommission/RecExCommissionFlags_jobOptions.py")

#important
ATLASCosmicFlags.doSim=False
ATLASCosmicFlags.JiveXML=False
ATLASCosmicFlags.OnlineJiveXML=False
ATLASCosmicFlags.doOnline=True
CBNTAthenaAware = False

#optional switches
ATLASCosmicFlags.IOVDbSvcGlobalTagData = 'COMCOND-001-00'
ATLASCosmicFlags.doMonitoring=True
ATLASCosmicFlags.doLVL1Calo=True
ATLASCosmicFlags.doInDet=False
ATLASCosmicFlags.doMuons=False
ATLASCosmicFlags.doLAr=False
ATLASCosmicFlags.doTile=False
ATLASCosmicFlags.doCaloTrkMuId=False


#
if ATLASCosmicFlags.doOnline:
    include.block("ByteStreamCnvSvc/BSEventStorageEventSelector_jobOptions.py")

# 
include( "RecExCommission/RecExCommission_topOptions.py" )

#include( "ByteStreamCnvSvcBase/BSAddProvSvc_RIO_jobOptions.py" )
#include( "ByteStreamCnvSvcBase/BSAddProvSvc_RDO_jobOptions.py" )

#
theApp.CreateSvc += [ "ByteStreamCnvSvc" ]
#svcMgr = theApp.serviceMgr()
#svcMgr += ByteStreamCnvSvc()
##

#ToolSvc.MuonCalib_MdtCalibDbAsciiTool.RT_InputFiles = [ "/preseries-scratch/masik/efatl/RT_default_comm.dat" ]
#ServiceMgr.ExceptionSvc.Catch="none"
try:
    ByteStreamAddressProviderSvc.TypeNames.remove("LArCellIDC/LArCellIDC")
    #temporary avoid crash on this container
    ByteStreamAddressProviderSvc.TypeNames.remove("Muon::TgcPrepDataContainer/TGC_Measurements")
except ValueError:
    pass

try:
    PoolSvc.ReadCatalog.remove("prfile:poolcond/PoolCat_comcond_castor.xml")
#    PoolSvc.ReadCatalog.remove("xmlcatalog_file:PoolCat_comcond_castor.xml")
except ValueError:
    pass

MessageSvc.defaultLimit=10000

Service( "ROBDataProviderSvc" ).filterEmptyROB=True
#ToolSvc.LArCablingService.OutputLevel=FATAL
ToolSvc.calonoisetool.OutputLevel=DEBUG
#TopoMaker.UseCaloNoiseTool=False
ToolSvc.EventData2XML.OutputLevel=ERROR
#ToolSvc.OutputLevel=ERROR

#to avoid crash in SCT monit
#ByteStreamCnvSvc.InitCnvs += [ "SCT_RDO_Container","InDetRawDataCollection<SCT_RDORawData>"]
#MessageSvc.OutputLevel=DEBUG


#LArDBConnection=" <dbConnection>impl=cool;techno=sqlite;schema= /det/lar/project/databases/COMCOND200.db;X:COMP200:</dbConnection>"

#PoolSvc.ReadCatalog+=["xmlcatalog_file:/det/lar/project/databases/intr130/PoolFileCatalog.xml"]
#PoolSvc.ReadCatalog+=["xmlcatalog_file:/det/lar/lar/project/databases/comcond.000002.lar_conditions.recon.pool.v0000/PoolFileCatalog.xml"]

theApp.Dlls += ["LArRawConditions","LArCondAthenaPoolPoolCnv","LArCalibTools","LArCalibUtils"]
theApp.Dlls += ["MuonByteStream"]

print 'jmasik'
if ATLASCosmicFlags.JiveXML:
    print ToolSvc.EventData2XML.DataTypes

print theApp.topAlgs
print topSequence


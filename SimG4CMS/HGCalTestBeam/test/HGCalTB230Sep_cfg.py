import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Modifier_hgcaltb_cff import hgcaltb

process = cms.Process('SIM', hgcaltb)

# import of standard configurations
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Geometry.HGCalTBCommonData.testTB230SepXML_cfi')
process.load('Geometry.HGCalCommonData.hgcalEENumberingInitialization_cfi')
process.load('Geometry.HGCalCommonData.hgcalEEParametersInitialization_cfi')
process.load('Geometry.HcalTestBeamData.hcalTB06Parameters_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedFlat_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SimG4CMS.HGCalTestBeam.HGCalTB23Analyzer_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

if 'MessageLogger' in process.__dict__:
    process.MessageLogger.G4cerr=dict()
    process.MessageLogger.G4cout=dict()
    process.MessageLogger.HGCSim=dict()
    process.MessageLogger.CaloSim=dict()
    process.MessageLogger.FlatThetaGun=dict()
    process.MessageLogger.FlatEvtVtx=dict()
    process.MessageLogger.HGCalGeom=dict()
    process.MessageLogger.HGCalGeomX=dict()

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('SingleMuonE200_cfi nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:gensim.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('TBGenSimSep.root')
                                   )

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.generator = cms.EDProducer("FlatRandomEThetaGunProducer",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        MinE = cms.double(199.99),
        MaxE = cms.double(200.01),
        MinTheta = cms.double(0.0),
        MaxTheta = cms.double(0.0),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359),
        PartID = cms.vint32(11)
    ),
    Verbosity = cms.untracked.int32(1),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single electron E 10')
)
process.VtxSmeared.MinZ = -100.0
process.VtxSmeared.MaxZ = -100.0
#process.VtxSmeared.MinX = -1.0
#process.VtxSmeared.MaxX =  1.0
#process.VtxSmeared.MinY = -1.0
#process.VtxSmeared.MaxY =  1.0
process.g4SimHits.OnlySDs = ['HGCalSensitiveDetector', 'HcalTB06BeamDetector']
process.g4SimHits.HGCSD.Detectors = 1
process.g4SimHits.HGCSD.RejectMouseBite = False
process.g4SimHits.HGCSD.RotatedWafer    = False

process.g4SimHits.CaloTrkProcessing.TestBeam = True
process.g4SimHits.HCalSD.ForTBHCAL = True
process.g4SimHits.NonBeamEvent = True
process.g4SimHits.UseMagneticField = False

process.g4SimHits.EventVerbose = 2
process.g4SimHits.SteppingVerbosity = 2
process.g4SimHits.StepVerboseThreshold= 0.1
process.g4SimHits.VerboseEvents = [1]
process.g4SimHits.VertexNumber = []
process.g4SimHits.VerboseTracks =[]

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.analysis_step = cms.Path(process.HGCalTB23Analyzer)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

process.g4SimHits.Physics.type = 'SimG4Core/Physics/FTFP_BERT_EMN'

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,
				process.simulation_step,
#        		        process.analysis_step,
				process.endjob_step,
				process.RAWSIMoutput_step,
				)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 



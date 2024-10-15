
import FWCore.ParameterSet.Config as cms

process = cms.Process("SaveNtuple")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output_ntuplev2.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                               'file:/gpfs/cms/users/jbabbar/work/VertexProd/CMSSW_14_1_0_pre6/work/TTbar_PU/246b9af9-570e-42e4-b188-7e96e65b2832.root'
                            )
)

from RecoMTD.TimingIDTools.mvaTrainingNtuple_cfi import mvaTrainingNtuple

process.mvaTrainingNtuple = mvaTrainingNtuple

process.p = cms.Path(process.mvaTrainingNtuple)


import FWCore.ParameterSet.Config as cms

from RecoMTD.TrackExtender.trackExtenderWithMTD_cfi import *
from RecoMTD.TimingIDTools.mtdTrackQualityMVA_cfi import *

from RecoMTD.TimingIDTools.mtdSoAProducer_cfi import *

mtdSoA = mtdSoAProducer.clone()


fastTimingGlobalRecoTask = cms.Task(trackExtenderWithMTD,mtdSoA,mtdTrackQualityMVA)
fastTimingGlobalReco = cms.Sequence(fastTimingGlobalRecoTask)

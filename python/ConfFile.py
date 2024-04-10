import FWCore.ParameterSet.Config as cms

process = cms.Process("L1Trigger")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('/store/data/Run2024A/Muon1/MINIAOD/PromptReco-v1/000/378/929/00000/ac375b83-e151-467c-979e-bd32d0727f28.root')
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('out_hist.root')
                                   )

process.demo = cms.EDAnalyzer('L1JetAnalyzer',

                              )


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data', '')

process.p = cms.Path(process.demo)

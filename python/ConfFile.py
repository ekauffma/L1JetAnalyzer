import FWCore.ParameterSet.Config as cms

process = cms.Process("L1Trigger")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/378/993/00000/43a9b66c-3aca-4ca1-8ae6-c0593056e33a.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/001/00000/5b29b633-a562-4a23-b4d2-841ec211acb5.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/378/985/00000/5182a2cb-adf0-4cd6-81fb-bb92322442a1.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/003/00000/b4eac488-d3eb-4717-a41d-af243c3d4bcc.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/378/994/00000/72c15c81-42bb-4dc6-a8c7-e89dede14170.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/378/981/00000/857dd3bf-c6e9-4494-96c1-d81b60b60f49.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/378/995/00000/be5c86d0-b977-4c25-8ec8-4ebbe4d1d1ad.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/378/997/00000/f081deca-9095-48c1-98bc-dbd022397291.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/011/00000/f42bf59f-87e0-4512-8875-1a37f8e402dd.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/378/996/00000/d8b7ef10-ab44-40d8-af50-e478891295d9.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/004/00000/d9f2f656-c485-4f3f-849c-6e90bf6d2acf.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/378/998/00000/34ea7a5e-de8b-42d4-84cf-e1a1c127d0c8.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/005/00000/d9ee3d69-0a29-4f1e-bb26-2813e57fdaff.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/006/00000/35da8f34-b6d9-41de-9eeb-cbd42c92ac35.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/378/999/00000/bf15972c-af0e-4cd2-9782-2a71b2e2b30f.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/007/00000/eef8c42b-fa6a-4f83-8f98-1db27bc3bfdd.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/002/00000/4dea65e1-5ef4-4302-8070-a3aef6e740ad.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/008/00000/c8b0109e-2631-44cf-a893-47fd9941a05d.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/000/00000/83d1b547-24ed-4e15-a0a1-7d70fd9f618d.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/009/00000/6b24aa28-b045-46e9-8528-2b836e7542d4.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/010/00000/eea14bea-b5e8-4d19-80c2-8caae01c93ba.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/012/00000/5d05e83c-d0fc-4658-820a-9a010b405f26.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/029/00000/6eef9f8f-99d9-4879-9690-7b0c0b12062d.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/031/00000/a869376f-b67c-4c95-a60b-8538510b5ea4.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/028/00000/beac79dd-e9b5-4ae3-a880-9ac837cd5826.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/058/00000/c49bb1af-540d-4cc0-87c0-01f51b54486d.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/058/00000/9c16e85a-a58c-4e63-a5e3-ec18b3024f53.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/058/00000/82f2bf98-7869-4165-9009-0a7cfe9f71b2.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/058/00000/bbd7a8bd-c450-4a6f-b403-ed62832eb737.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/058/00000/3c1c0c11-a9f4-4f32-9e45-5441b61c03b0.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/058/00000/af91404c-2449-4f65-b822-9c3d3d36bced.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/064/00000/c624b2c3-baca-45b8-b6b8-2091827ad7f0.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/065/00000/719de905-bedd-4f69-9c5b-38864c2813eb.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/142/00000/561305d6-2a96-4c6a-a2e1-0da47fe8ec3f.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/146/00000/91ee022a-1266-4842-8fe9-3633190ba027.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/075/00000/86128edd-937b-4162-b024-814af7f0b73e.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/075/00000/a9d5c81f-e402-4737-9f96-f0ef387d46eb.root',
'/store/data/Run2024B/Muon0/MINIAOD/PromptReco-v1/000/379/075/00000/2c1cf3ef-7fcf-4682-a7d7-5ac83dc73de6.root')
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('out_hist.root')
                                   )

process.demo = cms.EDAnalyzer('L1JetAnalyzer',

                              )


process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data', '')

process.p = cms.Path(process.demo)

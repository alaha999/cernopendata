import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
process = cms.Process("Demo")

#dataset

nEvents = 10

inputFileName = 'list.txt'
outputFileName = 'fileTree.root'

inputEleList = FileUtils.loadListFromFile(inputFileName);

readFilesInputEle = cms.untracked.vstring(*inputEleList)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
        )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

process.source = cms.Source("PoolSource",fileNames = readFilesInputEle)

#Creating Tree
process.ExportTree = cms.EdAnalyzer("Events")

#Calling C++ Script
process.demo = cms.EdAnalyzer('ElectronData',
                              minTracks = cms.untracked.uint32(0),
                              MC=cms.bool(True)
                              )

#Creating outputfile
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(outputFileName)
                                   )

process.p = cms.Path(process.demo)

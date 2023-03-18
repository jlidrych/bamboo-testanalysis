import logging
logger = logging.getLogger("reco base plotter")
from itertools import chain

from bamboo import treefunctions as op
from bamboo import treedecorators as td
from bamboo.root import gbl as ROOT
from bamboo.analysisutils import forceDefine

import definitions as defs
import utils
from basePlotter import basePlotter

class recoBasePlotter(basePlotter):
    """"""
    def __init__(self, args):
        super().__init__(args)

    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument("--btag-sf", choices=['none', 'fixWPMT', 'itFit'], default='fixWPMT', help="Choose b-tag SFs to use (default: %(default)s)")
        parser.add_argument("--decorr-btag", action="store_true", help="Decorrelate b-tagging uncertainties for different working points")
        parser.add_argument("--btag-fix", action="store_true", help="Fix: use preVFP mistag SFs in postVFP")
        parser.add_argument("--no-split-jer", action="store_true", help="Do not produce the split JER variations")
        parser.add_argument("--jes-scheme", choices=["All", "Merged", "Total"], default="Merged", help="JEC uncertainty scheme (default: %(default)s)")

    def prepareTree(self, tree, sample=None, sampleCfg=None, **kwargs):
        tree,noSel,be,lumiArgs = super().prepareTree(tree, sample=sample, sampleCfg=sampleCfg, **kwargs)
        era = sampleCfg["era"]

        # Rochester correction on-the-fly calculation
        if self.args.roc_corr:
            defs.addMuonRocCor(be, tree._Muon, era, sample, self.isMC(sample))

        # configure Jet/Met corrections

        isNotWorker = (self.args.distributed != "worker")
        defs.configureJetMETCorrections(tree, era, isNotWorker, self.isMC(sample), be, sample, splitJER=not self.args.no_split_jer, jesScheme=self.args.jes_scheme)

        return tree,noSel,be,lumiArgs

    def defineBaseSelections(self, tree, noSel, sample, sampleCfg):
        era = sampleCfg["era"]

        if self.isMC(sample):
            noSel = noSel.refine("puWeight", weight=defs.makePUWeight(tree, era, noSel))
        noSel = noSel.refine("flags", cut=defs.flagDef(tree.Flag, era, self.isMC(sample)))

        # calculate (corrected) muon 4-momenta before accessing them
        if self.args.roc_corr:
            forceDefine(tree._Muon.calcProd, noSel)

        twoMuTriggerSel,oneMuTriggerSel = defs.buildMuonSelections(tree, noSel, self.muons, self.vetoMuons, self.electrons, self.vetoElectrons, sample, era, self.isMC(sample))
        twoEleTriggerSel,oneEleTriggerSel = defs.buildElectronSelections(tree, noSel, self.muons, self.vetoMuons, self.electrons, self.vetoElectrons, sample, era, self.isMC(sample))

        # Pre-compute JEC and jetMet variations here, as this is the last step before using jets
        # Call it twice, but the selections are orthogonal, so it will only execute once per event
        # Only for MC since for data tree._Jet does not exist (we just use the default jet collection)
        if self.isMC(sample):
            forceDefine(tree._Jet.calcProd, oneMuTriggerSel)
            forceDefine(tree._Jet.calcProd, oneEleTriggerSel)
            forceDefine(tree._Jet.calcProd, twoMuTriggerSel)
            forceDefine(tree._Jet.calcProd, twoEleTriggerSel)
            if self.bTagWeight:
                forceDefine(self.bTagWeight, oneMuTriggerSel)
                forceDefine(self.bTagWeight, oneEleTriggerSel)
                forceDefine(tree._Jet.calcProd, twoMuTriggerSel)
                forceDefine(tree._Jet.calcProd, twoEleTriggerSel)

        return oneMuTriggerSel,oneEleTriggerSel,twoMuTriggerSel,twoEleTriggerSel

    def defineObjects(self, tree, noSel, sample=None, sampleCfg=None, **kwargs):
        super().defineObjects(tree, noSel, sample, sampleCfg, **kwargs)
        era = sampleCfg["era"]

        self.origMET = tree.MET
        self.corrMET = defs.corrMET(self.origMET, tree.PV, sample, era, self.isMC(sample))

        ###### Add lepton cone pt
        tree.Muon.valueType.conept = td.itemProxy(defs.muonConePt(tree.Muon))
        tree.Electron.valueType.conept = td.itemProxy(defs.eleConePt(tree.Electron))

        ##### Lepton definition and scale factors
        self.muons = op.select(tree.Muon, defs.muonDef(era))
        self.muon = self.muons[0]
        self.electrons = op.select(tree.Electron, defs.eleDef(era))
        self.electron = self.electrons[0]

        self.vetoMuons = op.select(tree.Muon, defs.vetoMuonDef)
        self.vetoElectrons = op.select(tree.Electron, defs.vetoEleDef)

        ##### Jet definition
        self.rawJets = op.select(tree.Jet, defs.jetDef)
        self.cleanedJets = defs.cleanJets(self.rawJets, self.muons, self.electrons)
        self.HT = op.rng_sum(self.cleanedJets, lambda j: j.pt)

        ##### FatJet definition
        self.rawFatJets = op.select(tree.FatJet, defs.fatjetDef)
        self.cleanedFatJets = defs.cleanJets(self.rawFatJets, self.muons, self.electrons)


        ##### B tagging definition
        # order jets by *decreasing* deepFlavour
        self.btagger = "btagDeepFlavB"
        self.cleanedJetsByDeepFlav = op.sort(self.cleanedJets, lambda jet: -getattr(jet, self.btagger))
        # DeepFlavour tagged jets
        self.bJetsM = defs.bTagDef(self.cleanedJets, era, "M", self.btagger)
        self.bJetsT = defs.bTagDef(self.cleanedJets, era, "T", self.btagger)
        self.lightJetsM = defs.lightTagDef(self.cleanedJets, era, "M", self.btagger)

        if self.isMC(sample):
            self.bFlavJets = op.select(self.cleanedJets, lambda j: j.hadronFlavour == 5)
            self.cFlavJets = op.select(self.cleanedJets, lambda j: j.hadronFlavour == 4)
            self.lFlavJets = op.select(self.cleanedJets, lambda j: j.hadronFlavour == 0)
            self.quarkJets = op.select(self.cleanedJets, lambda j: op.AND(j.hadronFlavour == 0, op.in_range(0, op.abs(j.partonFlavour), 4)))
            self.gluonJets = op.select(self.cleanedJets, lambda j: op.AND(j.hadronFlavour == 0, j.partonFlavour == 21))
            self.undefJets = op.select(self.cleanedJets, lambda j: op.AND(j.hadronFlavour == 0, j.partonFlavour == 0))

        # B-tag scale factors
        self.bTagWeight = None
        if self.isMC(sample):
            if self.args.btag_sf == "fixWPMT":
                systMapping = {
                    "pileup": "pileup", # correlated with "pileup" uncertainty
                    "isr": None, # uncorrelated, standalone b-tag uncertainty
                    "fsr": None,
                    "hdamp": None,
                    "qcdscale": None,
                    "topmass": None,
                    "type3": None,
                    "jes": None,
                }
                if self.args.no_split_jer:
                    systMapping["jer"] = "jer"
                else:
                    systMapping["jer0"] = "jer"
                    systMapping["jer1"] = "jer"
                self.bTagWeightFun = defs.makeBtagSFWPs(era, noSel, wps=["M", "T"],
                                                        decorr_wps=self.args.decorr_btag, full_scheme=True, full_scheme_mapping=systMapping,
                                                        prepostVFPfix=self.args.btag_fix)
            elif self.args.btag_sf == "itFit":
                self.bTagWeightFun = defs.makeBtagSFItFit(era, noSel)
            self.bTagWeight = op.rng_product(self.cleanedJets, self.bTagWeightFun)

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None, **kwargs):
        super().postProcess(taskList, config, workdir, resultsdir, **kwargs)

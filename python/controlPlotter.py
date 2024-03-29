import os
import json
import yaml
import glob
import shutil

import logging
logger = logging.getLogger("control plotter")

from bamboo.analysisutils import parseAnalysisConfig
from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo.plots import VariableBinning as VarBin
from bamboo.plots import CutFlowReport
from bamboo import treefunctions as op

import definitions as defs
import utils
import controlPlotDefinition as cp
from recoBasePlotter import recoBasePlotter

class controlPlotter(recoBasePlotter):
    """"""
    def __init__(self, args):
        super(controlPlotter, self).__init__(args)

    def definePlots(s, t, noSel, sample=None, sampleCfg=None):
        super().defineObjects(t, noSel, sample, sampleCfg)
        era = sampleCfg["era"]
        plots = []

        ##### Muon and electron selection, exclusive!
        oneMuTriggerSel, oneEleTriggerSel,twoMuTriggerSel, twoEleTriggerSel, oneEleoneMuTriggerSel  = recoBasePlotter.defineBaseSelections(s, t, noSel, sample, sampleCfg)
        ###### Plots for ==1 lepton (trigger level) ######

        # plots += utils.makeMergedPlots([("1mu", oneMuTriggerSel), ("1ele", oneEleTriggerSel)], "1lep", "nJets", EqBin(10, 0, 10), title="Number of jets", var=op.rng_len(s.cleanedJets))
        # plots += utils.makeMergedPlots([("1mu", oneMuTriggerSel), ("1ele", oneEleTriggerSel)], "1lep", "nVtx", EqBin(100, 0, 100), title="Number of primary vertices", var=t.PV.npvs)

        # for sel,lep,name in [(oneMuTriggerSel, s.muon, "1mu"), (oneEleTriggerSel, s.electron, "1ele")]:
            # plots += cp.makeLeptonPlots(sel, lep, name)
            # plots += cp.makeMETPlots(sel, lep, s.corrMET, name)

            # plots.append(Plot.make1D("{}_MET_dPhi".format(name), op.Phi_mpi_pi(s.corrMET.phi - s.origMET.phi), sel,
                # EqBin(200, -3.1416, 3.1416), title="MET corr phi - MET phi", plotopts=utils.getOpts(name)))
            # plots.append(Plot.make2D("{}_MET_pt_corr_vs_nocorr".format(name), (s.corrMET.pt, s.origMET.pt), sel, (EqBin(20, 0, 500), EqBin(20, 0, 500)), xTitle="corr MET pt", yTitle="orig MET pt", plotopts=utils.getOpts(name)))
            # plots.append(Plot.make2D("{}_MET_phi_corr_vs_nocorr".format(name), (s.corrMET.phi, s.origMET.phi), sel, (EqBin(20, -3.1416, 3.1416), EqBin(20, -3.1416, 3.1416)), xTitle="corr MET phi", yTitle="orig MET phi", plotopts=utils.getOpts(name)))

        ###### Cutflow report

        yields = CutFlowReport(s.yieldsPrefix)
        plots.append(yields)
        ###### Plots for ==1 lepton, >=0 jets ######
        oneMu0JetSel = oneMuTriggerSel.refine("1muon_0jets",cut=True)
        oneEle0JetSel = oneEleTriggerSel.refine("1ele_0jets",cut=True)
        twoMu0JetSel = twoMuTriggerSel.refine("2muon_0jets",cut=True)
        twoEle0JetSel = twoEleTriggerSel.refine("2ele_0jets",cut=True)
        oneEleoneMu0JetSel = oneEleoneMuTriggerSel.refine("elemuon_0jets",cut=True)
        
        yields.add(oneMu0JetSel, "1mu0j")
        yields.add(oneEle0JetSel, "1ele0j")
        yields.add(twoMu0JetSel, "2mu0j")
        yields.add(twoEle0JetSel, "2ele0j")
        yields.add(oneEleoneMu0JetSel, "1ele1mu0j")
        ###### Plots for ==1 lepton, >=1 jets ######
        oneJetSel = op.rng_len(s.cleanedJets) >= 1
        oneMu1JetSel = oneMuTriggerSel.refine("1muon_1jets", cut=oneJetSel)
        oneEle1JetSel = oneEleTriggerSel.refine("1ele_1jets", cut=oneJetSel)
        yields.add(oneMu1JetSel, "1mu1j")
        yields.add(oneMu1JetSel, "1ele1j")
        oneFatJetSel = op.rng_len(s.cleanedFatJets) >=1
        oneMu1FatJetSel = oneMuTriggerSel.refine("1muon_1fatjets",cut=oneFatJetSel)
        oneEle1FatJetSel = oneEleTriggerSel.refine("1ele_1fatjets",cut=oneFatJetSel)
        yields.add(oneMu1FatJetSel, "1mu1fj")
        yields.add(oneMu1FatJetSel, "1ele1fj")
        
        for sel,lep,name in [(oneMu0JetSel, s.muon, "1mu_0j"), (oneEle0JetSel, s.electron, "1ele_0j")]:
             plots += cp.makeLeptonPlots(sel, lep, name)
             plots += cp.makeJetPlots(sel, s.cleanedJets, name)
             plots += cp.makeBJetPlots(sel, s.cleanedJetsByDeepFlav, name + "_byDeepFlav")
             plots += cp.makeMETPlots(sel, lep, s.corrMET, name)

        for sel,lep,name in [(oneMu1FatJetSel, s.muon, "1mu_1fj"), (oneEle1FatJetSel, s.electron, "1ele_1fj")]:
             plots += cp.makeFatJetPlots(sel, s.cleanedFatJets, name,1)

        for sel,lep,name in [(twoMu0JetSel, s.muons, "2mu_0j"), (twoEle0JetSel, s.electrons, "2ele_0j")]:
            plots += cp.makeDileptonPlots(sel, lep, name)
            plots += cp.makeJetPlots(sel, s.cleanedJets, name)

        for sel,ele,mu,name in [(oneEleoneMu0JetSel,s.electrons,s.muons,"elemu_0j")]:
            plots += cp.makeEleMuPlots(sel, ele,mu,name)
            plots += cp.makeJetPlots(sel, s.cleanedJets, name)
            plots += cp.makeFatJetPlots(sel, s.cleanedFatJets, name,1)
            plots += cp.makeBJetPlots(sel, s.cleanedJetsByDeepFlav, name + "_byDeepFlav")
            plots += cp.makeMissingETPlots(sel, s.corrMET, name)


        plots += utils.makeMergedPlots(
            [(f"2mu_0j", twoMu0JetSel), (f"2ele_0j", twoEle0JetSel)],
            f"2lep_0j", "MET_pt", EqBin(60, 0, 600), title="MET p_{T} (GeV)", var=s.corrMET.pt)

        plots += utils.makeMergedPlots(
            [(f"2mu_0j", twoMu0JetSel), (f"2ele_0j", twoEle0JetSel)],
            f"2lep_0j", "MET_phi", EqBin(60, -3.1416, 3.1416), title="MET #phi", var=s.corrMET.phi)

        plots += utils.makeMergedPlots(
            [(f"2mu_0j", twoMu0JetSel), (f"2ele_0j", twoEle0JetSel)],
            f"2lep_0j", "nJets", EqBin(9, 0.5, 9.5), title="Number of jets", var=op.rng_len(s.cleanedJets))

        plots += utils.makeMergedPlots(
            [(f"2mu_0j", twoMu0JetSel), (f"2ele_0j", twoEle0JetSel)],
            f"2lep_0j", "nBDeepFlavM", EqBin(7, -0.5, 6.5), title="Number of medium b-tagged jets", var=op.rng_len(s.bJetsM))

        ###### Plots for ==1 lepton, >=2 jets ######
        twoJetSel = op.rng_len(s.cleanedJets) >= 2
        oneMu2JetSel = oneMuTriggerSel.refine("1muon_2jets", cut=twoJetSel)
        oneEle2JetSel = oneEleTriggerSel.refine("1ele_2jets", cut=twoJetSel)
        yields.add(oneMu2JetSel, "1mu2j")
        yields.add(oneMu2JetSel, "1ele2j")
        for sel,lep,name in [(oneMu2JetSel, s.muon, "1mu_2j"), (oneEle2JetSel, s.electron, "1ele_2j")]:
            plots += cp.makeDijetPlots(sel, s.cleanedJets, name)
        ##### Plots for ==1 lepton, >=2 jets, >= 1 b jets ######

        oneMu2Jet1BSel = oneMu2JetSel.refine("1muon_2jets_1b", cut=op.rng_len(s.bJetsM) >= 1, weight=s.bTagWeight)
        oneEle2Jet1BSel = oneEle2JetSel.refine("1ele_2jets_1b", cut=op.rng_len(s.bJetsM) >= 1, weight=s.bTagWeight)
        yields.add(oneMu2Jet1BSel, "1mu2j1b")
        yields.add(oneEle2Jet1BSel, "1ele2j1b")

        plots += utils.makeMergedPlots(
            [(f"1mu_2j_1b", oneMu2Jet1BSel), (f"1ele_2j_1b", oneEle2Jet1BSel)],
            f"1lep_2j_1b", "nVtx", EqBin(100, 0, 100), title="Number of primary vertices", var=t.PV.npvs)

        plots += utils.makeMergedPlots(
            [(f"1mu_2j_1b", oneMu2Jet1BSel), (f"1ele_2j_1b", oneEle2Jet1BSel)],
            f"1lep_2j_1b", "nJets", EqBin(8, 1.5, 9.5), title="Number of jets", var=op.rng_len(s.cleanedJets))

        plots += utils.makeMergedPlots(
            [(f"1mu_2j_1b", oneMu2Jet1BSel), (f"1ele_2j_1b", oneEle2Jet1BSel)],
            f"1lep_2j_1b", "nBDeepFlavM", EqBin(6, 0.5, 6.5), title="Number of medium b-tagged jets", var=op.rng_len(s.bJetsM))

        for sel,lep,name in [(oneMu2Jet1BSel, s.muon, f"1mu_2j_1b"), (oneEle2Jet1BSel, s.electron, f"1ele_2j_1b")]:
            plots += cp.makeLeptonPlots(sel, lep, name, binScaling=2)
            # plots += cp.makeJetPlots(sel, s.cleanedJets, name, binScaling=2, allJets=True)
            # plots += cp.makeBJetPlots(sel, s.cleanedJetsByDeepFlav, name + "_byDeepFlav")
            # plots += cp.makeMETPlots(sel, lep, s.corrMET, name, binScaling=2)
            # plots += cp.makeHEMPlots(sel, lep, s.cleanedJets, name)

        return plots

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super(controlPlotter, self).postProcess(taskList, config, workdir, resultsdir,
             makeBU=True,
             removeBatch=True,
             createEnvelope=True,
             moveSystHists=True,
             removeSignalOverlap=False)
        self.runPlotIt(taskList, config, workdir, resultsdir)

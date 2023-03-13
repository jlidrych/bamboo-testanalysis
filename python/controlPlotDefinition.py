from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

import definitions as defs
import utils

def makeLeptonPlots(sel, lepton, uname, binScaling=1):
    plots = []

    if "mu" in uname:
        flav = "Muon"
    if "ele" in uname:
        flav = "Electron"

    plots.append(Plot.make1D(f"{uname}_lep_pt", lepton.pt, sel,
            EqBin(60 // binScaling, 30., 530.), title="%s p_{T} (GeV)" % flav,
            plotopts=utils.getOpts(uname)))
    plots.append(Plot.make1D(f"{uname}_lep_eta", lepton.eta, sel,
            EqBin(50 // binScaling, -2.4, 2.4), title="%s eta" % flav,
            plotopts=utils.getOpts(uname, **{"log-y": False})))
    plots.append(Plot.make1D(f"{uname}_lep_phi", lepton.phi, sel,
            EqBin(50 // binScaling, -3.1416, 3.1416), title="%s phi" % flav,
            plotopts=utils.getOpts(uname, **{"log-y": False})))

    return plots

def makeJetPlots(sel, jets, uname, maxJet=4, binScaling=1, allJets=False):
    plots = []

    if allJets:
        plots.append(Plot.make1D(f"{uname}_jets_pt", op.map(jets, lambda j: j.pt), sel,
                EqBin(70 // binScaling, 30., 730.), title="Jets p_{T} (GeV)",
                plotopts=utils.getOpts(uname)))
        plots.append(Plot.make1D(f"{uname}_jets_eta", op.map(jets, lambda j: j.eta), sel,
                EqBin(48 // binScaling, -2.4, 2.4), title=f"Jets #eta",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
        plots.append(Plot.make1D(f"{uname}_jets_phi", op.map(jets, lambda j: j.phi), sel,
                EqBin(50 // binScaling, -3.1416, 3.1416), title="Jets #phi",
                plotopts=utils.getOpts(uname, **{"log-y": False})))

    for i in range(maxJet):
        plots.append(Plot.make1D(f"{uname}_jet{i+1}_pt", jets[i].pt, sel,
                EqBin(60 // binScaling, 30., 730. - min(4, i) * 100), title=f"{utils.getCounter(i+1)} jet p_{{T}} (GeV)",
                plotopts=utils.getOpts(uname)))
        plots.append(Plot.make1D(f"{uname}_jet{i+1}_eta", jets[i].eta, sel,
                EqBin(50 // binScaling, -2.4, 2.4), title=f"{utils.getCounter(i+1)} jet #eta",
                plotopts=utils.getOpts(uname, **{"log-y": False})))

    return plots

def makeFatJetPlots(sel, fatjets, uname, maxJet=4, binScaling=1, allJets=False):
    plots = []

    if allJets:
        plots.append(Plot.make1D(f"{uname}_fatjets_pt", op.map(fatjets, lambda fj: fj.pt), sel,
                EqBin(70 // binScaling, 30., 730.), title="Fatjets p_{T} (GeV)",
                plotopts=utils.getOpts(uname)))
        plots.append(Plot.make1D(f"{uname}_fatjets_eta", op.map(fatjets, lambda fj: fj.eta), sel,
                EqBin(48 // binScaling, -2.4, 2.4), title=f"Fatjets #eta",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
        plots.append(Plot.make1D(f"{uname}_fatjets_phi", op.map(fatjets, lambda fj: fj.phi), sel,
                EqBin(50 // binScaling, -3.1416, 3.1416), title="Fatjets #phi",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
        plots.append(Plot.make1D(f"{uname}_fatjets_mass", op.map(fatjets, lambda fj: fj.mass), sel,
                EqBin(50 // binScaling, 0., 250.), title="Fatjets mass (GeV)",
                plotopts=utils.getOpts(uname)))
        plots.append(Plot.make1D(f"{uname}_fatjets_sdmass", op.map(fatjets, lambda fj: fj.msoftdrop), sel,
                EqBin(50 // binScaling, 0., 250.), title="Fatjets SD mass (GeV)",
                plotopts=utils.getOpts(uname)))

    for i in range(maxJet):
        plots.append(Plot.make1D(f"{uname}_fatjet{i+1}_pt", fatjets[i].pt, sel,
                EqBin(80 // binScaling, 200., 1000. - min(4, i) * 100), title=f"{utils.getCounter(i+1)} fatjet p_{{T}} (GeV)",
                plotopts=utils.getOpts(uname)))
        plots.append(Plot.make1D(f"{uname}_fatjets{i+1}_mass", fatjets[i].mass, sel,
                EqBin(50 // binScaling, 0., 250.), title=f"{utils.getCounter(i+1)} fatjets mass (GeV)",
                plotopts=utils.getOpts(uname)))
        plots.append(Plot.make1D(f"{uname}_fatjets{i+1}_sdmass", fatjets[i].msoftdrop, sel,
                EqBin(50 // binScaling, 0., 250.), title=f"{utils.getCounter(i+1)} fatjets SD mass (GeV)",
                plotopts=utils.getOpts(uname)))

    return plots

def makeExtraJetPlots(sel, jets, uname, binScaling=2, pTthresh=30):
    plots = []

    for i in range(2):
        plots.append(Plot.make1D(f"{uname}_extra_jet{i+1}_pt", jets[i].pt, sel,
                EqBin(60 // binScaling, pTthresh, 600. + pTthresh - i*200), title=f"{utils.getCounter(i+1)} extra jet p_{{T}} (GeV)",
                plotopts=utils.getOpts(uname)))
        plots.append(Plot.make1D(f"{uname}_extra_jet{i+1}_eta", jets[i].eta, sel,
                EqBin(50 // binScaling, -2.4, 2.4), title=f"{utils.getCounter(i+1)} extra jet #eta",
                plotopts=utils.getOpts(uname, **{"log-y": False})))
    plots.append(Plot.make1D(f"{uname}_extra_jet_DR", op.deltaR(jets[0].p4, jets[1].p4),
            sel, EqBin(60 // binScaling, 0.4, 3.), title="Extra jets #Delta R(bb)",
            plotopts=utils.getOpts(uname, **{"log-y": False})))
    plots.append(Plot.make1D(f"{uname}_extra_jet_M", op.invariant_mass(jets[0].p4, jets[1].p4),
            sel, EqBin(60 // binScaling, 15., 435.), title="Extra jets M(bb) (GeV)",
            plotopts=utils.getOpts(uname)))
    plots.append(Plot.make1D(f"{uname}_extra_jet_pT", (jets[0].p4 + jets[1].p4).Pt(),
            sel, EqBin(60 // binScaling, 20., 620.), title="Extra jets p_{T}(bb) (GeV)",
            plotopts=utils.getOpts(uname)))

    return plots

def makeBJetPlots(sel, jets, uname):
    plots = []

    for i in range(4):
        plots.append(Plot.make1D(f"{uname}_jet{i+1}_deepFlav", jets[i].btagDeepFlavB,
                sel, EqBin(60, 0., 1.), title=f"{utils.getCounter(i+1)}-highest jet deepFlavour", plotopts=utils.getOpts(uname)))

    return plots

def makeMETPlots(sel, lepton, met, uname, binScaling=1):
    plots = []

    plots.append(Plot.make1D(f"{uname}_MET_pt", met.pt, sel,
            EqBin(60 // binScaling, 0., 600.), title="MET p_{T} (GeV)",
            plotopts=utils.getOpts(uname)))
    plots.append(Plot.make1D(f"{uname}_MET_phi", met.phi, sel,
            EqBin(60 // binScaling, -3.1416, 3.1416), title="MET #phi",
            plotopts=utils.getOpts(uname, **{"log-y": False})))
    # plots.append(Plot.make1D(f"{uname}_MET_lep_deltaPhi",
    #         op.Phi_mpi_pi(lepton.phi - met.phi), sel, EqBin(60 // binScaling, -3.1416, 3.1416),
    #         title="#Delta #phi (lepton, MET)", plotopts=utils.getOpts(uname, **{"log-y": False})))

    MT = op.sqrt( 2. * met.pt * lepton.p4.Pt() * (1. - op.cos(op.Phi_mpi_pi(met.phi - lepton.p4.Phi()))) )
    plots.append(Plot.make1D(f"{uname}_MET_MT", MT, sel,
            EqBin(60 // binScaling, 0., 600.), title="Lepton M_{T} (GeV)",
            plotopts=utils.getOpts(uname)))

    return plots

def makeDijetPlots(sel,jet, uname, binScaling=1):
    plots = []

    mjj = op.invariant_mass(jet[0].p4,jet[1].p4)
    deltaeta = op.abs(jet[0].p4.Eta() - jet[1].p4.Eta())
    deltaphi = op.Phi_mpi_pi(jet[0].p4.Phi() - jet[1].p4.Phi())

    plots.append(Plot.make1D(f"{uname}_dijet_mjj", mjj, sel,
                             EqBin(200 // binScaling, 0., 1000.), title="dijet M_{jj} (GeV)",
                             plotopts=utils.getOpts(uname)))
    plots.append(Plot.make1D(f"{uname}_dijet_abseta", deltaeta, sel,
                             EqBin(60 // binScaling, 0., 5.), title="dijet #Delta#eta_{jj}",
                             plotopts=utils.getOpts(uname)))
    plots.append(Plot.make1D(f"{uname}_dijet_deltaphi", deltaphi, sel,
                             EqBin(60 // binScaling, -3.14, 3.14), title="dijet #Delta#phi_{jj}",
                             plotopts=utils.getOpts(uname)))

    return plots

def makeHEMPlots(sel, lepton, jets, uname):
    plots = []

    jets_eta = op.map(jets, lambda jet: jet.eta)
    jets_phi = op.map(jets, lambda jet: jet.phi)
    jets_pt = op.map(jets, lambda jet: jet.pt)

    jetsHEM = op.select(jets, lambda jet: op.AND(jet.eta < -1.3, jet.phi > -1.57, jet.phi < -0.87))
    jetsHEM_pt = op.map(jetsHEM, lambda jet: jet.pt)
    jetsHEM_sel = sel.refine(uname + "HEM", cut=op.rng_len(jetsHEM) > 0)

    jetsNoHEM = op.select(jets, lambda jet: op.OR(op.AND(jet.eta < -1.3, op.OR(jet.phi < -1.57, jet.phi > -0.87)), jet.eta > 1.3))
    jetsNoHEM_pt = op.map(jetsNoHEM, lambda jet: jet.pt)

    jetsNoHEM_sel = sel.refine(uname + "NoHEM", cut=op.rng_len(jetsNoHEM) > 0)

    plots.append(Plot.make1D(f"{uname}_jet_pt", jets_pt, sel,
            EqBin(30, 30., 730.), title="Jet p_{T} (GeV)", plotopts=utils.getOpts(uname)))
    plots.append(Plot.make1D(f"{uname}_HEM_jet_pt", jetsHEM_pt, jetsHEM_sel,
            EqBin(30, 30., 730.), title="Jet p_{T} inside HEM failure (GeV)",
            plotopts=utils.getOpts(uname)))
    plots.append(Plot.make1D(f"{uname}_NoHEM_jet_pt", jetsNoHEM_pt, jetsNoHEM_sel,
            EqBin(30, 30., 730.), title="Jet p_{T} outside HEM failure (GeV)",
            plotopts=utils.getOpts(uname)))
    plots.append(Plot.make2D(f"{uname}_jet_eta_vs_phi", (jets_eta, jets_phi), sel,
            (EqBin(60, -2.4, 2.4), EqBin(60, -3.1416, 3.1416)), xTitle="Jet eta", yTitle="Jet phi", plotopts=utils.getOpts(uname)))
    plots.append(Plot.make2D(f"{uname}_lep_eta_vs_phi", (lepton.eta, lepton.phi), sel,
            (EqBin(60, -2.4, 2.4), EqBin(60, -3.1416, 3.1416)), xTitle="Lepton eta", yTitle="Lepton phi", plotopts=utils.getOpts(uname)))

    return plots

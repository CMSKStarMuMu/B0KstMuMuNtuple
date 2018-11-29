import ROOT
import numpy

from ROOT import TGraph, TFile, TMultiGraph
import math, array

n_points  = 15
n_samples = 11
scale     = 0.00225 ##1/440

infile = TFile('outcome_wp_finding_massFlatten_onlyRT_NoB0PsiCut.root','read')

graphs_vy = []
for i in range(n_samples):
    graphs_vy_sig_tmp = infile.Get('graph_signal_sample%s'%i).GetY()
    graphs_vy_bkg_tmp = infile.Get('graph_background_sample%s'%i).GetY()

    graphs_vy_sig_tmp.SetSize(infile.Get('graph_signal_sample%s'%i).GetN())
    graphs_vy_bkg_tmp.SetSize(infile.Get('graph_background_sample%s'%i).GetN())

    ## now convert to usual python objects
    vy = []
    vy_sig = array.array('f',graphs_vy_sig_tmp)
    vy_bkg = array.array('f',graphs_vy_bkg_tmp)
    for j in range(len(vy_sig)):
        vy.append(scale*vy_sig[j] / math.sqrt( scale*vy_sig[j] + vy_bkg[j]))
    
    graphs_vy_tmp2 = numpy.asarray(vy )
    graphs_vy.append(graphs_vy_tmp2)

    if i==0: 
        graphs_vx_tmp = infile.Get('graph_signal_sample%s'%i).GetX()
        graphs_vx_tmp.SetSize(infile.Get('graph_signal_sample%s'%i).GetN())
        vx = array.array('f',graphs_vx_tmp)
        graphs_vx = numpy.asarray(vx )


norm_graphs = []
for i in range(n_samples):
    graphs_vy_array = numpy.asarray(graphs_vy[i] )
    
    norm_graphs.append( ROOT.TGraph(len(graphs_vx),graphs_vx, graphs_vy_array )) # convert lists to arrays
    norm_graphs[i].GetXaxis().SetTitle('bdt cut')
    norm_graphs[i].GetYaxis().SetTitle('S/#sqrt{S+B} a.u.')
    norm_graphs[i].SetTitle('')
    norm_graphs[i].SetMarkerStyle(8)
    norm_graphs[i].SetName('graph_sample%s'%i)
    norm_graphs[i].SetMarkerColor(ROOT.kViolet+i)
    norm_graphs[i].SetMarkerStyle(20+i)


mg = TMultiGraph()
for i in range(n_samples):
  mg.Add(norm_graphs[i])

canv = ROOT.TCanvas()
mg.Draw('AP')
mg.GetXaxis().SetTitle('BDT score')
mg.GetYaxis().SetTitle('S/#sqrt{S+B}')

leg = ROOT.TLegend(0.4,0.16,0.7,0.46)
for i in range(11):
    leg.AddEntry(norm_graphs[i], 'BDT %s'%i, 'p')
leg.SetBorderSize(0)
leg.Draw()
canv.SaveAs('bdt_wps_MassFlatten_onlyRT_noB0PsiCut_scale%s.pdf'%scale)


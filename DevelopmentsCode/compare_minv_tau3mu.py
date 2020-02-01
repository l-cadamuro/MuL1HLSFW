import ROOT
import sys
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

toplot = [
    'true',
    'm0',
    'vivado_hls'

]

plots = {
    'true' : 'minv_ta3mu_true.txt',
    'm0'   : 'minv_ta3mu_m0_true.txt',
    'vivado_hls' : '/home/zynq/luca/MuL1HLSFW/MuonAlgorithms/tau3mu_inv_mass_hw_results.txt',
}

colors = {
    'true'       : ROOT.kGray+1,
    'm0'         : ROOT.kBlack,
    'vivado_hls' : ROOT.kRed,
}

legend = {
    'true'       : "True m_{#mu#mu} from L1 TT",
    'm0'         : "True m_{#mu#mu} from L1 TT w/ m(#mu) = 0 approx.",
    'vivado_hls' : "Vidado HLS impl."
}

data = {}

for pl in toplot:
    print '... doing ', pl
    fin = open(plots[pl])
    data[pl] = [float(x.strip()) for x in fin]
    fin.close()

xmmins = [min(x) for x in data.values()]
xmmaxs = [max(x) for x in data.values()]

# xmin = min(xmmins)
# xmax = max(xmmaxs)

xmin = 1.0
xmax = 2.

delta = xmax-xmin

histos = {}
for tp in toplot:
    histos[tp] = ROOT.TH1F('h_' + tp, ';m_{#mu#mu};a.u.', 100, xmin - 0.1*delta , xmax + 0.1*delta)
    for d in data[tp]:
        histos[tp].Fill(d)

frame = ROOT.TH1F('h_framte', ';m_{#mu#mu};a.u.', 100, xmin - 0.1*delta , xmax + 0.1*delta)
ymmins = [h.GetMinimum() for h in histos.values()]
ymmaxs = [h.GetMaximum() for h in histos.values()]
ymin = 0 # min(ymmins)
ymax = max(ymmaxs)

frame.SetMinimum(ymin)
frame.SetMaximum(1.15*ymax)

# frame.SetMinimum(1.e-5)
# frame.SetMaximum(1.15*ymax)

c1 = ROOT.TCanvas('c1', 'c1', 600, 600)
c1.SetFrameLineWidth(3)
# c1.SetLogy(True)
frame.Draw()

# now plot
for tp in toplot:
    h = histos[tp]
    if tp in colors:
        h.SetLineColor(colors[tp])
    h.Draw('same')

# legend
leg = ROOT.TLegend(0.1, 0.6, 0.5, 0.88)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
for tp in toplot:
    name = legend[tp] if tp in legend else tp
    leg.AddEntry(histos[tp], name, 'l')
leg.Draw()

c1.Update()
c1.Print('minv_comparison_tau3mu.pdf', 'pdf')
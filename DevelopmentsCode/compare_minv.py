import ROOT
import sys
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

toplot = [
    # 'true',
    # 'm0',
    'hw_float',
    # 'hw_dphi_4Kx9',
    'hw_dphi_2Kx18', ## the best
    # 'hw_dphi_1Kx36',
    # 'hw_full_ADDR6-3_OUT9-7', ## into a single LUT indexing th1, th2
    # 'hw_full_ADDR6-3_OUT18-10', ## into two LUT indexing th1, th2, each in output has half value of the function
    # 'hw_full_ADDR6-3_OUT32-14', ## just to show that more precision in output does not matter
    # 'hw_full_thetaCos12-3_9-1_floatDiv',
    # 'hw_full_thetaCos11-3_18-1_floatDiv', ## the best
    # 'hw_full_thetaCos10-3_36-1_floatDiv',
    # 'prova',
    # 'hw_full_bestphitheta_hwDiv_2Kx18',
    # 'full_hw_HLS_div',
    # 'uGTalgo_hw_cosh1Kx36',
    'uGTalgo_hw_cosh2Kx18', ## best results
    # 'uGTalgo_hw_cosh4Kx9',
    'vivado_hls'

]

plots = {
    'true'                               : 'minv_JPsiMuMu_true.txt',
    'm0'                                 : 'minv_JPsiMuMu_massless.txt',
    'hw_float'                           : 'minv_JPsiMuMu_floatprecision.txt', 
    'hw_dphi_4Kx9'                       : 'minv_JPsiMuMu_hwdphi_CosLUT_4Kx9.txt', 
    'hw_dphi_2Kx18'                      : 'minv_JPsiMuMu_hwdphi_CosLUT_2Kx18.txt', 
    'hw_dphi_1Kx36'                      : 'minv_JPsiMuMu_hwdphi_CosLUT_1Kx36.txt',
    'hw_full_ADDR6-3_OUT9-7'             : 'minv_JPsiMuMu_hwfull_ThetaLUT_4Kx9_ADDR6-3_OUT9-7.txt',
    'hw_full_ADDR6-3_OUT18-10'           : 'minv_JPsiMuMu_hwfull_ThetaLUT_4Kx9_ADDR6-3_OUT18-10.txt',
    'hw_full_ADDR6-3_OUT32-14'           : 'minv_JPsiMuMu_hwfull_ThetaLUT_4Kx9_ADDR6-3_OUT32-14.txt',
    'hw_full_thetaCos12-3_9-1_floatDiv'  : 'minv_hw_full_thetaCos12-3_9-1_floatDiv.txt', 
    'hw_full_thetaCos11-3_18-1_floatDiv' : 'minv_hw_full_thetaCos11-3_18-1_floatDiv.txt', 
    'hw_full_thetaCos10-3_36-1_floatDiv' : 'minv_hw_full_thetaCos10-3_36-1_floatDiv.txt', 
    'hw_full_bestphitheta_hwDiv_2Kx18'   : 'minv_hw_full_bestphitheta_hwDiv_2Kx18.txt',
    'full_hw_HLS_div'                    : 'minv_full_hw_HLS_div.txt',
    'uGTalgo_hw_floatprecision'          : 'minv_uGTalgo_hw_floatprecision.txt',
    # 'prova'                            : 'minv_JPsiMuMu_hwdphi_CosLUT_1Kx36_v2.txt', 
    # 'prova'                            : 'prova.txt',
    'uGTalgo_hw_cosh1Kx36'               : 'minv_uGTalgo_hw_cosh1Kx36.txt',
    'uGTalgo_hw_cosh2Kx18'               : 'minv_uGTalgo_hw_cosh2Kx18.txt',
    'uGTalgo_hw_cosh4Kx9'                : 'minv_uGTalgo_hw_cosh4Kx9.txt',
    'vivado_hls'                         : '/home/zynq/luca/MuL1HLSFW/MuonAlgorithms/inv_mass_hw_results.txt'
}

colors = {
    'true'                                : ROOT.kGray+1,
    'm0'                                  : ROOT.kBlack,
    'hw_float'                            : ROOT.kRed,
    'hw_dphi_4Kx9'                        : ROOT.kBlue,
    'hw_dphi_2Kx18'                       : ROOT.kGreen,
    'hw_dphi_1Kx36'                       : ROOT.kOrange,
    'hw_full_ADDR6-3_OUT9-7'              : ROOT.kViolet,
    'hw_full_ADDR6-3_OUT18-10'            : ROOT.kCyan,
    'hw_full_ADDR6-3_OUT32-14'            : ROOT.kGreen+1,
    'hw_full_thetaCos12-3_9-1_floatDiv'   : ROOT.kOrange,
    'hw_full_thetaCos11-3_18-1_floatDiv'  : ROOT.kBlue,
    'hw_full_thetaCos10-3_36-1_floatDiv'  : ROOT.kAzure,
    'prova'                               : ROOT.kBlack,
    'hw_full_bestphitheta_hwDiv_2Kx18'    : ROOT.kBlack,
    'full_hw_HLS_div'                     : ROOT.kOrange,
    'uGTalgo_hw_floatprecision'           : ROOT.kOrange,
    'uGTalgo_hw_cosh1Kx36'                : ROOT.kRed,
    'uGTalgo_hw_cosh2Kx18'                : ROOT.kBlack,
    'uGTalgo_hw_cosh4Kx9'                 : ROOT.kGreen,
    'vivado_hls'                          : ROOT.kBlue,

}

legend = {
    'true'                               : "True m_{#mu#mu}",
    'm0'                                 : "True m_{#mu#mu} w/ m(#mu) = 0 approx.",
    'hw_float'                           : "Floating pt. m_{#mu#mu} w/ digital inputs",
    'hw_dphi'                            : "HW cos(#Delta#varphi)",
    'hw_full_ADDR6-3_OUT9-7'             : 'hw full 6/3 #rightarrow 9/7',
    'hw_full_ADDR6-3_OUT18-10'           : 'hw full 6/3 #rightarrow 18/10',
    'hw_full_ADDR6-3_OUT32-14'           : 'hw full 6/3 #rightarrow 32/14',
    'hw_full_thetaCos12-3_9-1_floatDiv'  : "hw full cos#theta 12/3 #rightarrow 9/1+ float div",
    'hw_full_thetaCos11-3_18-1_floatDiv' : "hw full cos#theta 11/3 #rightarrow 18/1+ float div",
    'hw_full_thetaCos10-3_36-1_floatDiv' : "hw full cos#theta 10/3 #rightarrow 36/1+ float div",
    'hw_full_bestphitheta_hwDiv_2Kx18'   : "hw full hw div 2Kx18",

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

xmin = 3.0
xmax = 3.15

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
c1.Print('minv_comparison.pdf', 'pdf')
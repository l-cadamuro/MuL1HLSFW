import ROOT
import sys
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(True)

nin  = 'minv_JPsiMuMu.txt'

if len(sys.argv) > 1:
    nin = sys.argv[1]

opdf  = nin.replace('.txt', '.pdf')
oroot = nin.replace('.txt', '.root')

print '... running on file : ', nin
print '... out pdf name    : ', opdf
print '... out root name   : ', oroot

# fin = open('minv_Bsmumu.txt')
fin = open(nin)
data = [float(x.strip()) for x in fin]
# print data

xmin = min(data)
xmax = max(data)
delta = xmax-xmin

print xmin, xmax

h = ROOT.TH1D ('h', ';m_{#mu#mu};a.u.', 100, xmin - 0.1*delta , xmax + 0.1*delta)
for d in data:
    h.Fill(d)
c1 = ROOT.TCanvas('c1', 'c1', 600, 600)
h.Draw()
# c1.SetLogy()
c1.Print(opdf, 'pdf')

fOut = ROOT.TFile(oroot, 'recreate')
h.Write()
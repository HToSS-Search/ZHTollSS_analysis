# Standard importts
import os,sys,socket,argparse
import os
import ROOT
import math
from array import array
import numpy as np
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad
from ROOT import kBlack, kBlue, kRed
#import tdrstyle
#import CMS_lumi
import yaml

ROOT.gROOT.SetBatch(True)
ReducedBinning = False
correctforMCPU = False
#tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
# RooFit
ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libRooFitCore.so")
#ROOT.gROOT.SetStyle("Plain") # Not sure this is needed
ROOT.gSystem.SetIncludePath( "-I$ROOFITSYS/include/" )
# colors = [ROOT.kBlack,  ROOT.kRed, ROOT.kBlue,ROOT.kGreen+2,ROOT.kMagenta+1, ROOT.kOrange+1, ROOT.kTeal-1,ROOT.kRed-3, ROOT.kCyan+2]
colors = [45, 85]
markers = [20, 21, 22, 33, 47]
lumi_scale = {'UL2016_APV': 19500, 'UL2016': 16800,'UL2017': 41480, '2017F': 13540,'UL2018': 59830 } #in pb-1


def gaussFit(hist, masspoint, ctau, output, suffix=""):
    tmp_sigma, tmp_mean = hist.GetRMS(), hist.GetMean()
    xlabel = "m_{hh} (GeV)"
    mass = ROOT.RooRealVar("mass", xlabel, 0., 3.);
    dh_hist = ROOT.RooDataHist("dh_hist", "dh_hist", ROOT.RooArgList(mass),
                                    ROOT.RooFit.Import(hist));

    # plot the data hist with error from sum of weighted events
    frame = mass.frame(ROOT.RooFit.Title("Mass"))
    dh_hist.plotOn(frame, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.MarkerSize(2), ROOT.RooFit.MarkerColor(ROOT.kGray+3),ROOT.RooFit.LineColor(ROOT.kGray+3), ROOT.RooFit.Name("hist"))
    
    # create a simple gaussian pdf
    gauss_mean = ROOT.RooRealVar("mean", "mean", tmp_mean, tmp_mean-0.5, tmp_mean+0.5)
    gauss_sigma = ROOT.RooRealVar("sigma jer", "sigma gauss", 0.004, 0.004, 0.03)
    gauss = ROOT.RooGaussian("gauss", "gauss", mass, gauss_mean, gauss_sigma)

    alp=0.3

    minX, maxX = masspoint - 0.03, masspoint + 0.03 #Consider fit around the masspoint 1.2 of range 0.03 on either side
    print("Min,Max:",minX,maxX)

    gauss.fitTo(dh_hist,ROOT.RooFit.Save(),ROOT.RooFit.SumW2Error(True) ,ROOT.RooFit.Range(minX,maxX))
    gauss.plotOn(frame, ROOT.RooFit.Name("fit"))


    frame.SetMaximum(frame.GetMaximum() * 1.25)
    frame.SetMinimum(0)
    frame.getAttMarker().SetMarkerColor(ROOT.kBlue)
    frame.getAttMarker().SetMarkerSize(1.5)
    frame.getAttLine().SetLineColor(ROOT.kBlue)


    # add chi2 info
    chi2_text = ROOT.TPaveText(0.2, 0.6, 0.2, 0.88, "NBNDC")
    chi2_text.SetTextAlign(11)
    chi2_text.AddText("#chi^{2} fit = %s" %round(frame.chiSquare(),2))
    chi2_text.AddText("#sigma "+"= {} #pm {}".format(round(gauss_sigma.getVal(),3), round(gauss_sigma.getError(),3)))
    chi2_text.AddText("#mu "+"= {} #pm {}".format(round(gauss_mean.getVal(),3), round(gauss_mean.getError(),3)))

    chi2_text.SetTextSize(0.03)
    chi2_text.SetTextColor(2)
    chi2_text.SetShadowColor(0)
    chi2_text.SetFillColor(0)
    chi2_text.SetLineColor(0)
    frame.addObject(chi2_text)
    
    #Canvas Drawing
    cfit, pad1, pad2 = createCanvasPads("cfit",boundary=0.3)
    pad1.cd()
    pad1.SetLogx(False)
    frame.GetYaxis().SetTitle('Events/'+str(hist.GetBinWidth(1)))
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetLabelSize(0.)
    frame.GetXaxis().SetTitleOffset(1.2)
    frame.GetYaxis().SetTitleOffset(1.6)
    frame.GetXaxis().SetRangeUser(tmp_mean-0.1,tmp_mean+0.1)
    frame.Draw()

    ymax=frame.GetMaximum()
    minXt,maxXt=float(masspoint)-2.5*gauss_sigma.getVal(),float(masspoint)+2.5*gauss_sigma.getVal()
    l1 = ROOT.TLine(minXt,0,minXt,ymax)
    l2 = ROOT.TLine(maxXt,0,maxXt,ymax)
    l1.SetLineWidth(2)
    l2.SetLineWidth(2)
    l1.Draw("same")
    l2.Draw("same")
    

    legxlow, legxhigh = 0.75,0.9
    legend = ROOT.TLegend(legxlow,0.75,legxhigh,0.9)
    legend.AddEntry(frame.findObject("hist"),"ZH","lep")
    legend.AddEntry(frame.findObject("fit"),"gauss","l")
    legend.AddEntry(l1,"bounds","l")
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetTextSize(0.03)
    legend.Draw("same")

    pad2.cd()
    pad2.SetLogx(False)
    hpull = frame.pullHist("hist","fit")
    frame3 = mass.frame(ROOT.RooFit.Title("Pull Distribution")) ;
    hpull.SetMarkerSize(1.5)
    frame3.addPlotable(hpull,"P") ;
    frame3.GetYaxis().SetTitle("pull")
    frame3.GetXaxis().SetTitle(xlabel)
    frame3.GetYaxis().SetTitleSize(0.1)
    frame3.GetXaxis().SetTitleSize(0.1)
    frame3.GetYaxis().SetLabelSize(0.1)
    frame3.GetXaxis().SetLabelSize(0.1)
    frame3.GetXaxis().SetTitleOffset(1.2)
    frame3.GetYaxis().SetTitleOffset(0.6)
    frame3.GetXaxis().SetRangeUser(tmp_mean-0.1,tmp_mean+0.1)
    frame3.Draw()
    fit_filename = "fit_MS" + str(masspoint).replace(".","p") + "_ctau" + ctau + suffix
    # if not os.path.exists(fit_plot_directory): os.makedirs(fit_plot_directory)
    #cfit.SaveAs(os.path.join(output, fit_filename + ".pdf"))
    cfit.SaveAs(os.path.join(output, fit_filename + ".png"))
    del cfit

    return round(gauss_sigma.getVal(),5), round(gauss_sigma.getError(),5), round(gauss_mean.getVal(),5), round(gauss_mean.getError(),5)

def createCanvasPads(savename, boundary=0.25):
    # tdrStyle.SetPadTopMargin(0.05)
    # tdrStyle.SetPadBottomMargin(0.13)
    # tdrStyle.SetPadLeftMargin(0.16)
    # tdrStyle.SetPadRightMargin(0.02)
    ylength_c = 2400
    c = ROOT.TCanvas(savename, savename, 2200, ylength_c)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, boundary, 1, 1.0)
    pad1.SetTopMargin(0.07)
    pad1.SetBottomMargin(0.05)  # joins upper and lower plot
    pad1.SetLeftMargin(0.16)
    pad1.SetRightMargin(0.02)
    #pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, boundary+0.01)
    pad2.SetTopMargin(0.)  # joins upper and lower plot
    pad2.SetBottomMargin(0.13/boundary)
    pad2.SetLeftMargin(0.16)
    pad2.SetRightMargin(0.02)

    #pad2.SetGridx()
    pad2.Draw()

    return c, pad1, pad2


def main():
    parser = argparse.ArgumentParser(description='Plot stacked histogram')
    parser.add_argument("--input",  dest="fin",   help="input folder 1", type=str)
    parser.add_argument("-o","--output", dest="out", help="Output file name", type=str)
    parser.add_argument("-y", "--year",   dest="year",   help="data year", type=str)
    parser.add_argument("--hname", dest = "hname", help = "name of the histogram", type=str)

    args = parser.parse_args()
    lumi_factor = lumi_scale[args.year]

    cwd = os.getcwd()
    indir = cwd + '/'+ args.fin 
    filepath = indir+ "/ZH/"+ "output_" + 'ZHTollSS_MH125_MS1p2_ctauS0'+ ".root"
    print(filepath)
    
    fin = ROOT.TFile(filepath, "READ")
    histname = args.hname
    
    h_1 = fin.Get(histname)
    print(h_1.GetName())
    h_1.Scale(lumi_factor)
    mass = 1.2 #GeV
    ctau = '0'

    sigma,sigma_err,mean,mean_err = gaussFit(h_1 ,mass, ctau, args.out, "_"+histname)
    print(sigma, sigma_err)
    minX, maxX = float(mass)-2.5*sigma, float(mass)+2.5*sigma
    width = maxX - minX
    print("Check selected interval:",minX, maxX, width)
    
    return 0

if __name__ == '__main__':
    main()






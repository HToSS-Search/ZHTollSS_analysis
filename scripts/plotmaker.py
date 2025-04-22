# Standard importts
import os,sys,socket,argparse
import shutil
import ROOT
import math
from array import array
import numpy as np
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad
from ROOT import kBlack, kBlue, kRed
#import tdrstyle
#import CMS_lumi
import gc

# Collect previous objects and disable circular garbage collector (reference counter still works)
gc.collect()
gc.disable()

ROOT.gStyle.SetOptStat(111111)
ROOT.gROOT.SetBatch(True)
#tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gSystem.SetIncludePath( "-I$ROOFITSYS/include/" )

colors = [ROOT.kBlack,  ROOT.kRed, ROOT.kBlue,ROOT.kGreen+2,ROOT.kMagenta+1, ROOT.kOrange+1, ROOT.kTeal-1,ROOT.kRed-3, ROOT.kCyan+2]
markers = [20, 21, 22, 33, 47]
m_size = 2


def createCanvasPads(savename,boundary=0.25):
    ylength_c=2400
    c = TCanvas(savename, savename, 2200, ylength_c)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, boundary, 1, 1.0)
    pad1.SetTopMargin(0.07)
    pad1.SetBottomMargin(0.05)  # joins upper and lower plot
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.04)
    #pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0, 1, boundary+0.01)
    pad2.SetTopMargin(0.05)  # joins upper and lower plot
    pad2.SetBottomMargin(0.13/boundary)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.04)
    # pad2.SetBottomMargin(0.1)
    #pad2.SetGridx()
    pad2.Draw()

    return c, pad1, pad2

def createRatio(h1, h2):
    # npts = h1.GetNbinsX()
    # sfx=[None]*npts
    # sfx_err=[None]*npts
    # sfy=[None]*npts
    # sfy_err=[None]*npts

    # for i in range(npts):
    #     # mcEff, dataEff = 0., 0.
    #     # h2.GetPoint(i,mcEff,sfx[i])
    #     # h1.GetPoint(i,dataEff,sfx[i])
    #     # mcErr, dataErr = h2.GetErrorY(i), h1.GetErrorY(i)
    #     # sfx[i] = h2.GetPointX(i)
    #     mcEff, mcErr, dataEff, dataErr = h2.GetY()[i], h2.GetErrorY(i), h1.GetY()[i], h1.GetErrorY(i)
    #     sfx[i] = h2.GetX()[i]
    #     sfx_err[i] = h2.GetErrorX(i)
    #     sfy[i] = dataEff/mcEff if mcEff else 0.0
    #     sfy_err[i] = 0.0
    #     if dataEff and mcEff:
    #         sfy_err[i] = sfy[i] * ((dataErr / dataEff)**2 + (mcErr / mcEff)**2)**0.5
    # h3 = ROOT.TGraphErrors(npts,array('d',sfx),array('d',sfy),array('d',sfx_err),array('d',sfy_err))
    h3 = h1.Clone()
    h3.Divide(h2)
    h3.SetLineColor(kBlack)
    h3.SetMarkerStyle(21)
    h3.SetTitle("Data/MC")
    h3.SetMinimum(0.)
    h3.SetMaximum(2.)
    # Set up plot for markers and errors
    #h3.Sumw2()
    #h3.SetStats(0)
    #for i in range(npts):
    #    h3.SetBinContent(i+1,sfy[i])
    #    h3.SetBinError(i+1,sfy_err[i])
    # Adjust y-axis settings
    y = h3.GetYaxis()
    y.SetTitle("Data/MC")
    y.SetNdivisions(505)
    y.SetTitleSize(0.15)
    y.SetTitleFont(42)
    y.SetTitleOffset(0.5)
    y.SetLabelFont(42)
    #y.SetLabelOffset(0.007)
    y.SetLabelSize(0.14)

    # Adjust x-axis settings
    x = h3.GetXaxis()
    x.SetTitleSize(40)
    x.SetTitleFont(42)
    x.SetTitleOffset(4.0)
    x.SetLabelFont(42)
    x.SetLabelSize(0.)
    #x.SetLimits(0.,105.)

    return h3

def main():
    ROOT.TH1.AddDirectory(ROOT.kFALSE)
    parser = argparse.ArgumentParser(description='Plot stacked histogram')
    parser.add_argument("-y", "--year",   dest="year",   help="data year", type=str)
    # parser.add_argument("-s", "--signal",   dest="sig",   help="HtoSS_MS2_ctauS0 or HtoSS_MS2_ctauS1 etc", type=str)
    parser.add_argument("-o","--output", dest="out", help="Output file name", type=str)
    parser.add_argument("-i","--input", dest="input", help="Input directory name", type=str)
    #parser.add_argument("--i2","--input2", dest="input2", help="Input directory name pion", type=str)
    # parser.add_argument("-m","--multiply", dest="mult", default=1,help="Multiply signal by a factor", type=int)
    # parser.add_argument("-c","--category", dest="category",help="prompt/displacedmumu/displacedhh/displaced", type=str)
    # parser.add_argument("--tf", dest="tf", default=1,help="Transfer factor depending on CR - Check BkgEst.py", type=float)
    parser.add_argument("--yhigh", dest="yhigh", default=500,help="y-axis multiplicative factor for ymax", type=float)
    parser.add_argument("--log", dest="log", help="true for plotting with logY, false by default", action="store_true")
    parser.add_argument("--norm", dest="norm", help="true for plotting with normalisation, false by default", action="store_true")
    parser.add_argument("--cat", dest="cat", help="true for plotting the categorized plots, false by default", action="store_true")
    # parser.add_argument("--mass", dest="mass", help="mass of scalar", type=str)
    # parser.add_argument("--analysis", dest="analysis", help="true for plotting after analysis, false by default (after skim)", action="store_true")

    # parser.add_argument("--noratio", dest="noratio", help="true for not plotting with ratio, false by default", action="store_true")
    # add an option to plot just one plot accessible name in histo_dict; change savename accordingly
    args = parser.parse_args()
    # create required parts
    cwd = os.getcwd()


    #Set up the histogram dictionaries in same way as defined in RDataFrames_ZH_Analyzer.py
    histo_dict = {
        'CutFlow':{'hname':'h_CutFlow','label': "Cuts", 'xlow':-0.5,'xhigh':8.5,'hrebin':1, '2d':False},
        #'ZBosonMass': {'hname':'h_ZBosonMass','label': "m_{Z} (GeV)", 'xlow':70,'xhigh':110,'hrebin':1, '2d':False},
        'ChargedHadron1relIso': {'hname':'h_ChargedHadron1relIso','label': "leading h^{#pm}(s1) PF Rel. isolation", 'xlow':0.,'xhigh':10.,'hrebin':1, '2d':False},
        'ChargedHadron2relIso': {'hname':'h_ChargedHadron2relIso','label': "subleading h^{#pm}(s1) PF Rel. isolation", 'xlow':0.,'xhigh':10.,'hrebin':1, '2d':False},
        'ChargedHadron3relIso': {'hname':'h_ChargedHadron3relIso','label': "leading h^{#pm}(s2) PF Rel. isolation", 'xlow':0.,'xhigh':10.,'hrebin':1, '2d':False},
        'ChargedHadron4relIso': {'hname':'h_ChargedHadron4relIso','label': "subleading h^{#pm}(s2) PF Rel. isolation", 'xlow':0.,'xhigh':10.,'hrebin':1, '2d':False},
        'Scalar1Mass': {'hname':'h_Scalar1Mass','label': "m_{s1} (GeV)", 'xlow':0.8,'xhigh':3.6,'hrebin':1, '2d':False},
        'Scalar2Mass': {'hname':'h_Scalar2Mass','label': "m_{s2} (GeV)", 'xlow':0.8,'xhigh':3.6,'hrebin':1, '2d':False},
        'Scalar12Mass': {'hname':'h_Scalar12Mass','label': "m_{s2} (GeV)", 'xlow':0.8,'xhigh':3.6,'hrebin':1, '2d':False},
        #'Scalar1Lxy': {'hname':'h_Scalar1Lxy','label': "L_{xy}^{s1} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':False},
        #'Scalar1LxySig': {'hname':'h_Scalar1LxySig','label': "L_{xy}^{s1}/#Delta L_{xy}^{s1}", 'xlow':0.,'xhigh':200.,'hrebin':1, '2d':False},
        #'Scalar2Lxy': {'hname':'h_Scalar2Lxy','label': "L_{xy}^{s2} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':False},
        #'Scalar2LxySig': {'hname':'h_Scalar2LxySig','label': "L_{xy}^{s2}/#Delta L_{xy}^{s2}", 'xlow':0.,'xhigh':200.,'hrebin':1, '2d':False},
        #'HiggsBosonMass_Loose_Loose': {'hname':'h_HiggsBosonMass_Loose','label': "m_{H} (GeV)", 'xlow':110,'xhigh':140,'hrebin':1, '2d':False},
        'ChargedHadron1relIso_Loose': {'hname':'h_ChargedHadron1relIso_Loose','label': "leading h^{#pm}(s1) PF Rel. isolation", 'xlow':0.,'xhigh':10.,'hrebin':1, '2d':False},
        'ChargedHadron2relIso_Loose': {'hname':'h_ChargedHadron2relIso_Loose','label': "subleading h^{#pm}(s1) PF Rel. isolation", 'xlow':0.,'xhigh':10.,'hrebin':1, '2d':False},
        'ChargedHadron3relIso_Loose': {'hname':'h_ChargedHadron3relIso_Loose','label': "leading h^{#pm}(s2) PF Rel. isolation", 'xlow':0.,'xhigh':10.,'hrebin':1, '2d':False},
        'ChargedHadron4relIso_Loose': {'hname':'h_ChargedHadron4relIso_Loose','label': "subleading h^{#pm}(s2) PF Rel. isolation", 'xlow':0.,'xhigh':10.,'hrebin':1, '2d':False},
        'Scalar1Mass_Loose': {'hname':'h_Scalar1Mass_Loose','label': "m_{s1} (GeV)", 'xlow':0.8,'xhigh':3.6,'hrebin':1, '2d':False},
        'Scalar2Mass_Loose': {'hname':'h_Scalar2Mass_Loose','label': "m_{s2} (GeV)", 'xlow':0.8,'xhigh':3.6,'hrebin':1, '2d':False},
        'Scalar12Mass_Loose': {'hname':'h_Scalar12Mass_Loose','label': "m_{s2} (GeV)", 'xlow':0.8,'xhigh':3.6,'hrebin':1, '2d':False},
        #'Scalar1Lxy_Loose': {'hname':'h_Scalar1Lxy_Loose','label': "L_{xy}^{s1} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':False},
        #'Scalar1LxySig_Loose': {'hname':'h_Scalar1LxySig_Loose','label': "L_{xy}^{s1}/#Delta L_{xy}^{s1}", 'xlow':0.,'xhigh':200.,'hrebin':1, '2d':False},
        #'Scalar2Lxy_Loose': {'hname':'h_Scalar2Lxy_Loose','label': "L_{xy}^{s2} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':False},
        #'Scalar2LxySig_Loose': {'hname':'h_Scalar2LxySig_Loose','label': "L_{xy}^{s2}/#Delta L_{xy}^{s2}", 'xlow':0.,'xhigh':200.,'hrebin':1, '2d':False},
        'ScalarMass': {'hname':'h_ScalarMass','labelx': "m_{s1} (GeV)",'labely': "m_{s2} (GeV)", 'xlow':0.8,'xhigh':3.6,'hrebin':1, '2d':True},
        'ScalarMass_Loose': {'hname':'h_ScalarMass_Loose','labelx': "m_{s1} (GeV)",'labely': "m_{s2} (GeV)", 'xlow':0.8,'xhigh':3.6,'hrebin':1, '2d':True},
        #'ScalarLxy': {'hname':'h_ScalarLxy','labelx': "L_{xy}^{s1} (cm)" ,'labely': "L_{xy}^{s2} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':True},
        'ScalarLxySig': {'hname':'h_ScalarLxySig','labelx': "L_{xy}^{s1}/#Delta L_{xy}^{s1}" ,'labely': "L_{xy}^{s2}/#Delta L_{xy}^{s2}", 'xlow':0.,'xhigh':200.,'hrebin':1, '2d':True},
        #'ScalarLxy_Loose': {'hname':'h_ScalarLxy_Loose','labelx': "L_{xy}^{s1} (cm)", 'labely': "L_{xy}^{s2} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':True},
        'ScalarLxySig_Loose': {'hname':'h_ScalarLxySig_Loose','labelx': "L_{xy}^{s1}/#Delta L_{xy}^{s1}", 'labely': "L_{xy}^{s2}/#Delta L_{xy}^{s2}", 'xlow':0.,'xhigh':200.,'hrebin':1, '2d':True},
        #'Scalar1MassLxy': {'hname':'h_Scalar1MassLxy','labelx': "m_{s1} (GeV)" ,'labely': "L_{xy}^{s1} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':True},
        #'Scalar2MassLxy': {'hname':'h_Scalar1MassLxy','labelx': "m_{s2} (GeV)" ,'labely': "L_{xy}^{s2} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':True},
        #'Scalar1MassLxy_Loose': {'hname':'h_Scalar1MassLxy_Loose','labelx': "m_{s1} (GeV)" ,'labely': "L_{xy}^{s1} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':True},
        #'Scalar2MassLxy_Loose': {'hname':'h_Scalar1MassLxy_Loose','labelx': "m_{s2} (GeV)" ,'labely': "L_{xy}^{s2} (cm)", 'xlow':0.,'xhigh':20.,'hrebin':1, '2d':True},
    }


    cat_dict = {
        'prompt': {'id_': '_prompt', 'label': 'prompt'}, 
        'displaceds1': {'id_': '_displaceds1', 'label': 'displaced s1'},
        'displaceds2': {'id_': '_displaceds2', 'label': 'displaced s2'},
        'displaced': {'id_': '_displaced', 'label': 'displaced'}
    }

    datasets_dict = {
        'ZHTollSS_MH125_MS1p2_ctauS0':{'type':'signal','label':"#splitline{m_{S}=1.2 GeV,}{c#tau = 0.1mm}",'color':ROOT.kOrange+1,'integral':-1},
        'SingleMuonRun2017F':{'type':'data','label':"SingleMuonRun2017F",'color':ROOT.kGreen+3,'integral':-1}
    } 
    
    ctau_style_map = {'ctauS0':1,'ctauS1':5,'ctauS10':9,'ctauS100':10, 'Data':3}

    lumi_scale = {'UL2016_APV': 19500, 'UL2016': 16800,'UL2017':  41480,'UL2018': 59830, '2017F': 13540} #in pb-1
    lumi_factor = lumi_scale[args.year]


    indir = args.input
    for dname in datasets_dict:
        print(dname)
        if 'Run' in dname:
            print('here')
            tag = '_onlyData'
            fin = ROOT.TFile(indir+"/Data/"+ "output_" + dname + ".root", "READ")
        elif 'ZHTollSS' in dname:
            print('here1')
            tag = '_onlySignal'
            fin = ROOT.TFile(indir+"/ZH/"+ "output_" + dname + ".root", "READ")
        fin.ls()
        if args.cat:
            print('CATEGORY PLOTS')
            for cat_key in cat_dict:
                catName = cat_dict[cat_key]['id_']
                catLabel = cat_dict[cat_key]['label']
                for key in histo_dict:
                    if 'CutFlow' in key: 
                        continue
                    savename=key + catName + tag
                    hprop = histo_dict[key]
                    histname = hprop['hname'] + catName
                    print(histname)
                    h_1 = fin.Get(histname)
                    print(h_1.GetName())
                    leg = ROOT.TLegend(0.5, 0.65, 0.95, 0.92)
                    leg.SetBorderSize(0)
                    leg.SetFillStyle(0)
                    leg.SetNColumns(2)
                    leg.SetTextSize(0.022)
                    boundary_percent = 0.35
                    ylength_c = int(2400*(1-boundary_percent+0.15))
                    c1 = TCanvas(savename, savename, 2200, ylength_c)
                    pad1 = ROOT.TPad("pad1", "pad1", 0, 0, 1, 1)
                    pad1.SetTopMargin(0.06)
                    pad1.SetBottomMargin(0.15)
                    pad1.SetLeftMargin(0.16)
                    pad1.SetRightMargin(0.16)
                    #pad1.SetGridx()
                    pad1.Draw()
                    # Lower ratio plot is pad2
                    c1.cd()
                    pad1.cd()
                    
                    if ('ZHTollSS' in dname) and ('CutFlow' not in key):
                        h_1.Scale(lumi_factor)
                    #if args.norm:
                    #    h_1.Scale(1/h_1.Integral())
                    #datasets_dict[dname]['integral']=h_1.Integral()
                    if args.log:
                        # c1.SetLogy()
                        pad1.SetLogy()
                        if hprop['2d']: pad1.SetLogx()
                        # pad1.SetLogx()
                    if hprop['2d']:
                        h_1.GetXaxis().SetTitle(hprop['labelx'])
                        h_1.GetYaxis().SetTitle(hprop['labely'])
                        ROOT.gStyle.SetStatX(0.12)
                        ROOT.gStyle.SetStatW(0.12)
                        ROOT.gStyle.SetStatY(0.5)
                        h_1.Draw('COLZ')
                    else:
                        color_ = datasets_dict[dname]['color']
                        h_1.SetFillColor(color_)
                        h_1.SetFillStyle(0) # hollow hist
                        h_1.SetMarkerColor(color_)
                        h_1.SetLineColor(color_)
                        h_1.SetMarkerStyle(8)
                        #h_1.SetLineStyle(ctau_style_map[dname.split('_')[-1]])
                        h_1.SetLineWidth(4)
                        h_1.SetMarkerSize(0.8) if ('CutFlow' in key) else h_1.SetMarkerSize(m_size)
                        if (not 'CutFlow' in key):
                            h_1.GetXaxis().SetRangeUser(hprop['xlow'],hprop['xhigh'])
                            h_1.GetXaxis().SetTitle(hprop['label'])
                            h_1.GetXaxis().SetTitleSize(0.05)
                            h_1.GetYaxis().SetTitleSize(0.05)
                            h_1.GetXaxis().SetLabelSize(0.045)
                            h_1.GetYaxis().SetLabelSize(0.045)
                            h_1.GetXaxis().SetTitleOffset(1.1)
                            h_1.GetYaxis().SetTitleOffset(1.4)
                        drawTag = 'hist text' if ('CutFlow' in key) else 'hist' 
                        h_1.Draw(drawTag)
                    leg.AddEntry(h_1, datasets_dict[dname]['label'], "l")


                    leg.Draw()
                    pad1.Modified()
                    pad1.Update()           
                    #CMS_lumi.cmsText = 'CMS'
                    #CMS_lumi.writeExtraText = True
                    #CMS_lumi.extraText = 'Work in Progress'
                    #CMS_lumi.lumi_13TeV = args.year+" MC"
                    #CMS_lumi.lumiTextSize = 0.5
                    #CMS_lumi.cmsTextSize=1.
                    #CMS_lumi.CMS_lumi(pad1, 4, 11)

                    pad1.Modified()
                    pad1.Update()
                    c1.Modified()
                    c1.Update()
                    c1.SaveAs(args.out+'/'+str(cat_key)+'/'+savename+'.png')
        else:
            for key in histo_dict:
                leg = ROOT.TLegend(0.5, 0.65, 0.95, 0.92)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.SetNColumns(2)
                leg.SetTextSize(0.022)
                boundary_percent = 0.35
                ylength_c = int(2400*(1-boundary_percent+0.15))
                c1 = TCanvas(savename, savename, 2200, ylength_c)
                pad1 = ROOT.TPad("pad1", "pad1", 0, 0, 1, 1)
                pad1.SetTopMargin(0.06)
                pad1.SetBottomMargin(0.15)
                pad1.SetLeftMargin(0.16)
                pad1.SetRightMargin(0.16)
                #pad1.SetGridx()
                pad1.Draw()
                # Lower ratio plot is pad2
                c1.cd()
                pad1.cd()
                hprop = histo_dict[key]
                histname = hprop['hname']
                h_1 = fin.Get(histname)
                print(h_1.GetName())

                #h_1.Rebin(hprop['hrebin'])
                if ('ZHTollSS' in dname) and ('CutFlow' not in key):
                    h_1.Scale(lumi_factor)
                #if args.norm:
                #    h_1.Scale(1/h_1.Integral())
                #datasets_dict[dname]['integral']=h_1.Integral()
                if args.log:
                    # c1.SetLogy()
                    pad1.SetLogy()
                    if hprop['2d']: pad1.SetLogx()
                    # pad1.SetLogx()
                if hprop['2d']:
                    h_1.GetXaxis().SetTitle(hprop['labelx'])
                    h_1.GetYaxis().SetTitle(hprop['labely'])
                    ROOT.gStyle.SetStatX(0.12)
                    ROOT.gStyle.SetStatW(0.12)
                    ROOT.gStyle.SetStatY(0.5)
                    h_1.Draw('COLZ')
                else:
                    color_ = datasets_dict[dname]['color']
                    h_1.SetFillColor(color_)
                    h_1.SetFillStyle(0) # hollow hist
                    h_1.SetMarkerColor(color_)
                    h_1.SetLineColor(color_)
                    h_1.SetMarkerStyle(8)
                    #h_1.SetLineStyle(ctau_style_map[dname.split('_')[-1]])
                    h_1.SetLineWidth(4)
                    h_1.SetMarkerSize(0.8) if ('CutFlow' in key) else h_1.SetMarkerSize(m_size)
                    if (not 'CutFlow' in key):
                        h_1.GetXaxis().SetRangeUser(hprop['xlow'],hprop['xhigh'])
                        h_1.GetXaxis().SetTitle(hprop['label'])
                        h_1.GetXaxis().SetTitleSize(0.05)
                        h_1.GetYaxis().SetTitleSize(0.05)
                        h_1.GetXaxis().SetLabelSize(0.045)
                        h_1.GetYaxis().SetLabelSize(0.045)
                        h_1.GetXaxis().SetTitleOffset(1.1)
                        h_1.GetYaxis().SetTitleOffset(1.4)
                    drawTag = 'hist text' if ('CutFlow' in key) else 'hist' 
                    h_1.Draw(drawTag)
                leg.AddEntry(h_1, datasets_dict[dname]['label'], "l")


                leg.Draw()
                pad1.Modified()
                pad1.Update()           
                #CMS_lumi.cmsText = 'CMS'
                #CMS_lumi.writeExtraText = True
                #CMS_lumi.extraText = 'Work in Progress'
                #CMS_lumi.lumi_13TeV = args.year+" MC"
                #CMS_lumi.lumiTextSize = 0.5
                #CMS_lumi.cmsTextSize=1.
                #CMS_lumi.CMS_lumi(pad1, 4, 11)

                pad1.Modified()
                pad1.Update()
                c1.Modified()
                c1.Update()
                c1.SaveAs(args.out+'/'+savename+'.png')
                #c1.SaveAs(args.out+'/'+savename+'.pdf')
            fin.Close()
    gc.enable()

if __name__ == '__main__':
    main()

import ROOT
import os,argparse

# Define the exponential function
def exponential_fit(x, par):
    return par[0] * ROOT.TMath.Exp(-par[1] * x[0])



def transfer_factor(category, output, hist, xlow_CR1=110, xhigh_CR1 =122.5, xlow_SR=122.5, xhigh_SR=127.5, xlow_CR2=127.5, xhigh_CR2=140):
    hist.Rebin(2)
    total = hist.Integral()

    #Define the mass ranges
    mass = ROOT.RooRealVar("mass", "m(hh#mu#mu)", 125, 110, 140)
    mass.setRange("loM", 110, 122.5 )
    mass.setRange("hiM", 127.5, 140 )
    mass.setRange("loSB", 120, 122.5 )
    mass.setRange("hiSB", 127.5, 130 )
    mass.setRange("peak", 122.5, 127.5 )
    mass.setRange("peakcore", 124.2, 125.8 )
    mass.setRange("full", 110, 140 )
    fit_range = "loM,hiM"
    
    #Set up the dataHist over the full 110 to 140 range
    data_all = ROOT.RooDataHist("data_all","data_all", ROOT.RooArgList(mass),ROOT.RooFit.Import(hist) )
    #Define the background model slope parameter
    alpha_all = ROOT.RooRealVar("alpha_all", "alpha_all", -0.1, -0.2, 0.)
    # alpha_all = ROOT.RooRealVar("alpha_all", "alpha_all", -0.01, -1, -0.005)
    
    #Set up the background model (exponential pdf)
    model_bkg_all = ROOT.RooExponential("model_bkg_all", "model_bkg_all", mass, alpha_all )
    #This is the normalisation parameter
    N_all=ROOT.RooRealVar("N_all", "N_all", 0, 5*total);
    
    #Create and extended pdf, making the normalisation a free parameter as well
    extmodel_bkg_all=ROOT.RooExtendPdf("extmodel_bkg_all", "extmodel_bkg_all", model_bkg_all, N_all, "full");

    #Only fit the exponential model only to the sidebands CR (110 to 122.5) and (127.5 to 140)
    extmodel_bkg_all.fitTo(data_all, ROOT.RooFit.Range(fit_range))
    
    canvas = ROOT.TCanvas(hist.GetName(), "Exponential Fit Canvas", 800, 600)
    plot = mass.frame()
    data_all.plotOn( plot,ROOT.RooFit.Name("dist"), ROOT.RooFit.MarkerColor(0), ROOT.RooFit.LineColor(ROOT.kGreen) )
    extmodel_bkg_all.plotOn(plot,ROOT.RooFit.NormRange(fit_range),ROOT.RooFit.Range("full"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("extfit"), ROOT.RooFit.LineStyle(2))
    
    # add chi2 info
    chi2_text = ROOT.TPaveText(0.2,0.7,0.2,0.9,"NBNDC")
    chi2_text.SetTextAlign(11)
    chi2_text.AddText("#chi^{2} fit = %s" %round(plot.chiSquare("extfit","dist"),2))
    chi2_text.AddText("#alpha "+"= {} #pm {}".format(round(alpha_all.getVal(),3), round(alpha_all.getError(),3)) )
    chi2_text.AddText("N "+"= {} #pm {}".format(round(N_all.getVal(),3), round(N_all.getError(),3)) )

    chi2_text.SetTextSize(0.03)
    chi2_text.SetTextColor(2)
    chi2_text.SetShadowColor(0)
    chi2_text.SetFillColor(0)
    chi2_text.SetLineColor(0)
    chi2_text.SetName("chi2_text")
    plot.addObject(chi2_text)

    # Plot the histogram and the fit function
    hist.GetXaxis().SetRangeUser(110,140)
    hist.GetXaxis().SetTitle("m(hh#mu#mu)")
    hist.SetMarkerStyle(8)
    hist.SetMarkerSize(2)
    hist.Draw("pe")

    plot.Draw("same")
    chi2_text.Draw("same")
    

    #TF Calculation
    bkg_SR = extmodel_bkg_all.createIntegral(ROOT.RooArgSet(mass), ROOT.RooFit.Range("peak")).getVal()
    bkg_CR = extmodel_bkg_all.createIntegral(ROOT.RooArgSet(mass), ROOT.RooFit.Range(fit_range)).getVal()
    print(bkg_SR, bkg_CR)
    tf = bkg_SR / bkg_CR

    # Show the canvas
    canvas.Update()
    canvas.SaveAs(os.path.join(output, category+"_TF_Fit_"+hist.GetName()+".png"))
    #mass.Delete();data_all.Delete();alpha_all.Delete();model_bkg_all.Delete();N_all.Delete();plot.Delete();extmodel_bkg_all.Delete();
    del mass;del data_all;del alpha_all;del model_bkg_all;del N_all;del plot;del extmodel_bkg_all;
    return tf

def main():
    parser = argparse.ArgumentParser(description='Plot stacked histogram')
    parser.add_argument("--input",  dest="fin",   help="input folder 1", type=str)
    parser.add_argument("-o","--output", dest="out", help="Output file name", type=str)
    #parser.add_argument("-y", "--year",   dest="year",   help="data year", type=str)
    parser.add_argument("--cat", dest="cat",   help="category", type=str)
    parser.add_argument("--hname", dest = "hname", help = "name of the histogram without category handle", type=str)

    args = parser.parse_args()
    #lumi_factor = lumi_scale[args.year]

    cwd = os.getcwd()
    indir = cwd + '/'+ args.fin
    filepath = indir+ "/Data/"+ "output_" + 'SingleMuonRun2017F'+ ".root"
    print(filepath)
    fin = ROOT.TFile(filepath, "READ")
    
    #Fit on either prompt Higgs mass spectrum or the combined displaced categories
    if args.cat == 'prompt':
        histname = args.hname + '_prompt'
        h_1 = fin.Get(histname)
    else:
        h1 = fin.Get(args.hname + '_displaceds1')
        h2 = fin.Get(args.hname + '_displaceds2')
        h3 = fin.Get(args.hname + '_displaced')
        
        h_1 = h1.Clone()
        h_1.Add(h2)
        h_1.Add(h3)
        

    transferFactor = transfer_factor(args.cat, args.out, h_1)
    print(transferFactor)
    return 0

if __name__ == '__main__':
    main()


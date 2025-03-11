def createCanvasPads(savename,boundary=0.25):
    ylength_c=2400
    c = TCanvas(savename, savename, 2200, ylength_c)
    # Upper histogram plot is pad1
    # pad1 = TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
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


def createRatio(h1, h2, errflag=False): #h1 is hdata, h2 is hmc; h1/h2; h1 is MC, h2 is MCTruth
    npts_data = h1.GetN()
    npts_mc = h2.GetNbinsX()
    npts = min(npts_mc,npts_data)
    print(npts_data,npts_mc)
    # npts = npts_mc
    sfx=[0]*npts
    sfx_err=[0]*npts
    sfx_err_hi=[0]*npts
    sfx_err_lo=[0]*npts
    sfy=[0]*npts
    sfy_err=[0]*npts
    sfy_err_lo=[0]*npts
    sfy_err_hi=[0]*npts

    for i in range(npts):
        if errflag:
            mcEff, mcErr, dataEff, dataErr = h2.GetBinContent(i), 0., h1.GetBinContent(i), h1.GetBinError(i)
            sfx_err_hi[i], sfx_err_lo[i] = h1.GetErrorXhigh(i), h1.GetErrorXlow(i)
            sfx[i] = h1.GetX()[i]
        else:
            mcEff, mcErr, dataEff, dataErr = h2.GetBinContent(i), h2.GetBinError(i), h1.GetY()[i], h1.GetErrorY(i)
            sfx_err_hi[i], sfx_err_lo[i] = h2.GetBinWidth(i), h2.GetBinWidth(i)
            sfx[i] = h2.GetBinCenter(i)
        sfy[i] = dataEff/mcEff if mcEff else 0.0
        sfy_err[i] = 0.0
        if dataEff and mcEff:
            sfy_err[i] = sfy[i] * ((dataErr / dataEff)**2 + (mcErr / mcEff)**2)**0.5
        sfy_err_lo[i],sfy_err_hi[i] = sfy_err[i], sfy_err[i]
        sfx_err[i] = sfx_err_lo[i], sfy_err_hi[i]
        print (mcEff,mcErr,dataEff, dataErr,h1.GetX()[i], sfy[i],sfy_err[i],sfx_err[i]) #get last point in graph and print error
    h3 = ROOT.TGraphAsymmErrors(npts,array('d',sfx),array('d',sfy),array('d',sfx_err_lo),array('d',sfx_err_hi),array('d',sfy_err_lo), array('d',sfy_err_hi))
    h3.SetLineColor(kBlue)
    h3.SetMarkerStyle(8)
    return h3

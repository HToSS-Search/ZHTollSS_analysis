import ROOT
import numpy as np


outFile = ROOT.TFile("plots/UL2017/MC.root", "READ")

i = 0
canvasNameList = ["canvas_turnon_HLT_IsoMu27", "canvas_turnon_HLT_IsoMu24", "canvas_turnon_HLT_IsoTkMu24", "canvas_pt", "canvas_mu_eta", "canvas_mu_phi", "canvas_deltaR", "canvas_invmass"]
print("Following canvasses will be displayed")
print(canvasNameList)
canvasList = []
for canvasName in canvasNameList:
	c = outFile.Get(canvasName)
	c.Draw()
	c.Update()

outFile.Close()
print("\n Program running... Press Enter to stop.")
input()  # Waits for the Enter key to be pressed
print("Program stopped.")



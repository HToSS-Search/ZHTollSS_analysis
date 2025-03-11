import ROOT
import numpy as uu
import os, argparse, re

ROOT.gStyle.SetOptStat(111110)

#Custom object which contains name of object and object to be displayed
class displayObj:
    def __init__(self, objStr, obj):
        self.name = objStr
        self.obj = obj

#Helper functions for canvas generation
def displayTurnOn(turnon, triggerName, dataName):
    canvas = ROOT.TCanvas('canvas_turnon_' + triggerName + '_' + dataName, "Turn On Curve " + triggerName, 800, 600)
    canvas.SetGrid()
    turnon.Draw("AP")
    xpos = 0
    if triggerName[-2:].isdigit():
        xpos = int(triggerName[-2:])
        vLine = ROOT.TLine(xpos, canvas.GetUymin(), xpos, canvas.GetUymax())
        vLine.SetLineStyle(2)
        vLine.SetLineWidth(2)
        vLine.SetLineColor(ROOT.kBlue)
        vLine.Draw("same")
    else:
        None

    return canvas
    
def displayHist(hist, varName, dataName, particleName):
    canvas = ROOT.TCanvas('canvas_hist_' + particleName + '_' + varName + '_' + dataName, particleName +' ' + varName + ' Distribution', 800, 600)
    hist.Draw('HIST')

    return canvas

def display2DHist(hist, name, dataName):
    canvas = ROOT.TCanvas('canvas_hist2D_' + name + '_' + dataName, name + ' 2D Distribution', 800, 600)
    hist.Draw('COLZ')
    ROOT.gStyle.SetStatX(0.3)  # Move left (default is 0.9)
    ROOT.gStyle.SetStatY(0.9)
    return canvas

def displayPTHist(listOfPt, dataName, particleName):
    canvas = ROOT.TCanvas('canvas_hist_' + particleName + '_leadingsubleadingpt_' + dataName, particleName + ' pT Distribution', 800, 600)
    i = 0
    maxlist = [hist.GetMaximum() for hist in listOfPt]
    max_ = max(maxlist) * 1.1
    for hist in listOfPt:
        hist.SetMaximum(max_)
        if (i==1):
            hist.Draw("HIST same")
            hist.SetLineColor(ROOT.kRed)
        else:
            hist.Draw("HIST ")
        i+=1
    max_ = max(maxlist) * 1.1
    return canvas

def displayCompare(compareList, varName, year):
    canvas = ROOT.TCanvas('canvas_compare2_' + varName, 'Data vs MC ' + varName, 800, 600)
    
    if 'Data' in compareList[1].name:
        temp_displayObj = compareList[0]
        compareList[0] = compareList[1]
        compareList[1] = temp_displayObj
    else:
        None
    compareList[0].obj.SetLineColor(ROOT.kRed)
    
    # Get the maximum values of the histograms
    max_data = compareList[0].obj.GetMaximum()
    max_mc = compareList[1].obj.GetMaximum()
                
    # Determine the overall maximum and add padding
    max_value = max(max_data, max_mc) * 1.1  # Add 10% padding
    
    # Set the maximum for the upper pad
    compareList[0].obj.SetMaximum(max_value)
    compareList[1].obj.SetMaximum(max_value)
    
    ROOT.SetOwnership(compareList[0].obj, False)
    ROOT.SetOwnership(compareList[1].obj, False)
    ratioPlot = ROOT.TRatioPlot(compareList[0].obj, compareList[1].obj, "divsym")
    ROOT.SetOwnership(ratioPlot, False)

    
    ratioPlot.SetH1DrawOpt("HIST E")
    ratioPlot.SetH2DrawOpt("HIST E")
    
    ratioPlot.Draw()
    
    ratioPlot.GetUpperRefObject().SetTitle('Comparison ' + varName + ' ' + year)
    ratioPlot.GetLowerRefYaxis().SetTitle('Data/MC')
    ratioPlot.GetLowerRefGraph().SetLineColor(ROOT.kBlack)
    ratioPlot.GetXaxis().SetTitle(varName)

    ratioMax = 1.5
    ratioMin = 0.5
    
    ratioPlot.GetLowerRefGraph().GetYaxis().SetRangeUser(ratioMin, ratioMax)
    
    ratioPlot.GetUpperPad().cd()
    legend = ROOT.TLegend(0.8, 0.55, 0.95, 0.65)
    legend.AddEntry(compareList[0].obj, 'Data', 'l')
    legend.AddEntry(compareList[1].obj, 'MC', 'l')
    legend.SetTextSize(0.03)
    legend.Draw("SAME")
    ratioPlot.GetUpperPad().GetListOfPrimitives().Add(legend)
    ratioPlot.GetUpperPad().Update()
    ratioPlot.GetUpperPad().Draw()
    return canvas


def displayObjects(listToDisplay, compareTag):
    canvasList = []
    leadingSubleading = []
    if compareTag:
        varNameUsed = []
        compareDuo = []
        for i, displayObj in enumerate(listToDisplay):
            idx_ = displayObj.name.find('_')
            idx2_ = displayObj.name.find("_", idx_ + 1)
            idx3_= displayObj.name.find('_', idx2_ + 1)
            idx4_ = displayObj.name.find('_', idx3_ + 1)
            idx5_ = displayObj.name.find('_', idx4_ + 1)
            varName = displayObj.name[idx2_+1:idx3_]
            year = displayObj.name[idx4_+1:]
            if (varName not in varNameUsed):
                varNameUsed.append(varName)
                for j in range(i+1, len(listToDisplay)):
                    displayObjToCompare = listToDisplay[j]
                    if varName in displayObjToCompare.name:
                        print(f"match: {displayObj.name} vs {displayObjToCompare.name}")
                        compareDuo = [displayObj, displayObjToCompare]
                        compareDuoPrint = [displayObj.name, displayObjToCompare.name]
                        print(compareDuoPrint)
                        canvasList.append(displayCompare(compareDuo, varName, year))
                    else:
                        None
                        #print(f"No match: {displayObj.name} vs {displayObjToCompare.name}")
            else:
                None
                #print(f"{varName} already in:", varNameUsed)

    else:
        #Only display the selected objects by the user
        for displayObj in listToDisplay:
        #Different implementations needed for different types of plots
            if 'turnon' in displayObj.name:
                idx_ = displayObj.name.find('_')
                idxlast_ = displayObj.name.rfind('_')
                triggerName = displayObj.name[idx_:idxlast_]
                dataName = displayObj.name[idxlast_:]
                canvasList.append(displayTurnOn(displayObj.obj, triggerName, dataName))
            elif ('hist2D' in displayObj.name):
                idx_ = displayObj.name.find('_')
                idxlast_ = displayObj.name.rfind('_')
                name = displayObj.name[idx_ + 1:idxlast_]
                dataName = displayObj.name[idxlast_:]
                canvasList.append(display2DHist(displayObj.obj, name, dataName))
            elif ( ('hist' in displayObj.name) and not ('2D' in displayObj.name) ):
                idx_ = displayObj.name.find('_')
                idx2_ = displayObj.name.find("_", idx_ + 1)
                idx3_= displayObj.name.find('_', idx2_ + 1)
                idxlast_ = displayObj.name.rfind('_')
                particleName = displayObj.name[idx_ + 1:idx2_]
                varName = displayObj.name[idx2_+1:idx3_]
                dataName = displayObj.name[idxlast_:]
                if ('pT' in varName) and ('leading' in displayObj.name):
                    leadingSubleading.append(displayObj.obj)
                    print(displayObj.name)
                    #Only if there are two histograms (leading and subleading pTs) draw them on canvas
                    if (len(leadingSubleading) == 2):
                        canvasList.append(displayPTHist(leadingSubleading, dataName, particleName))
                        #Remove the already drawn pair of leading and subleading histograms from list ST other leading/subleading can be drawn together in the correct pairs
                        leadingSubleading = []
                    else: None
                else: 
                    canvasList.append(displayHist(displayObj.obj, varName, dataName, particleName))
            #ADD OTHER ELIFS LATER FOR OTHER PLOT IMPLEMENTATION
            else: None
    return canvasList       

parser = argparse.ArgumentParser(description='Display analysis plots')
parser.add_argument('-o', '--outputfile', dest='outputfile', help = 'Provide the .root file from which plots are read', type=str)
args = parser.parse_args()

#Open the ROOT File and print out content for user
outFile = ROOT.TFile(args.outputfile, "READ")
keys = outFile.GetListOfKeys()
keynames = [key.GetName() for key in keys]
print("\nThe objects inside the .root are the following:\n" + str(keynames) + "\n")

#Ask for user input which plots want to be seen
user_input1 = input(
    f"Please specify all objects to be displayed (comma-separated).\nType 'all' to display all loaded objects.\nType 'Data' or 'MC' for associated plots.\nType 'Compare' if you want to compare data and MC\n"
    )

#Generate and add to canvasList accordingly
inputDisplayList = []
listToDisplay = []
listOfData = []
compareTag = False
if (user_input1.lower() == 'all'):
    #if all selected use the complete list
    #A list of custom display objects is generated which takes the obj name and obj itself from the .root file
    i = 0
    for key in keys:
        objName = key.GetName()
        listToDisplay.append(displayObj(objName, outFile.Get(objName)))    
        i+=1
elif (user_input1.lower() == 'data') or (user_input1.lower() == 'mc'):
    user_input2 = input("Please select a year\n")
    for key in keys:
        objName = key.GetName()
        if (user_input2 in objName) and (user_input1 in objName) :
            print(objName)
            listToDisplay.append(displayObj(objName, outFile.Get(objName))) 
elif (user_input1.lower() == 'compare'):
    user_input2 = input("Please select a year to compare or type 'all' to add all data and compare to it\n")
    if user_input2.lower() == 'all':
        #implement summation of data hists and use the MC_Scaled2017Combined hist (need to generate this still in analysis script)
        for key in keys:
            objName = key.GetName()
            if ('Data' in objName) and ('turnon' not in objName) and ('subleading' not in objName):
                print(objName)
                listOfData.append(displayObj(objName, outFile.Get(objName)))
            if ('MC_TOTAL' in objName) and ('subleading' not in objName):
                listToDisplay.append(displayObj(objName, outFile.Get(objName)))
        compareTag = True
        print('___________________________________')
    else:
        for key in keys:
            objName = key.GetName()
            if (user_input2 in objName) and ('turnon' not in objName) and ('subleading' not in objName):
                print(objName)
                listToDisplay.append(displayObj(objName, outFile.Get(objName)))
                compareTag = True
else:
    inputDisplayList = [obj.strip() for obj in user_input.split(",")]
    #Based on the list of selected names of objects by user take objects from .root accordingly
    i = 0
    for inputObjName in inputDisplayList:
        listToDisplay.append(displayObj(inputObjName, outFile.Get(inputObjName)))
        i += 1

#Implementation of the combined comparison

if len(listOfData) != 0:
    listOfCombinedData = []
    varNameList = []
    for displayObj_data in listOfData:
        idx_ = displayObj_data.name.find('_')
        idx2_ = displayObj_data.name.find("_", idx_ + 1)
        idx3_= displayObj_data.name.find('_', idx2_ + 1)
        idxlast_ = displayObj_data.name.rfind('_')
        varName = displayObj_data.name[idx2_+1:idx3_]
        particleName = displayObj_data.name[idx_ + 1:idx2_]
        if varName not in varNameList:
            varNameList.append(varName)
    for varName in varNameList:
        listOfDataVarName = []
        for displayObj_data in listOfData:
            if varName in displayObj_data.name:
                listOfDataVarName.append(displayObj_data)
        print('GROUPED BY VAR')
        histDataCombinedVar = listOfDataVarName[0].obj.Clone()
        histDataCombinedVar.SetName('hist_'+particleName+'_'+varName+'_Data_TOTAL')
        histDataCombinedVar.SetTitle('Distribution ' + varName + ' Data combined') 
        for i in range(1, len(listOfDataVarName)):
            print(listOfDataVarName[i].name)
            histDataCombinedVar.Add(listOfDataVarName[i].obj)
        listOfCombinedData.append(histDataCombinedVar)
    for hist in listOfCombinedData:
        listToDisplay.append(displayObj(hist.GetName(), hist)) 
            

for display in listToDisplay:
    print(display.name)

canvasList = displayObjects(listToDisplay, compareTag)

#Draw and update the canvasses
for canvas in canvasList:
    canvas.Update()
    canvas.Draw()


outFile.Close()
print("\n Program running... Press Enter to stop.")
input()  # Waits for the Enter key to be pressed
print("Program stopped.")



import ROOT
import numpy as np


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
    
def displayHist(hist, varName, dataName):
    canvas = ROOT.TCanvas('canvas_hist_' + varName + '_' + dataName, 'Muon ' + varName + ' Distribution', 800, 600)
    hist.Draw('HIST E')

    return canvas

def displayPTHist(listOfPt, dataName):
    canvas = ROOT.TCanvas('canvas_hist_leadingsubleadingpt_' + dataName, 'Muon pT Distribution', 800, 600)
    i = 0
    for hist in listOfPt:
        if (i==1):
            hist.Draw("HIST E same")
        else:
            hist.Draw("HIST E")
        i+=1
    return canvas

def displayCompare(compareList, varName):
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

    legend = ROOT.TLegend(0.2, 0.2, 0.7, 0.7)
    legend.AddEntry(compareList[0].obj, 'Data', 'l')
    legend.AddEntry(compareList[1].obj, 'MC', 'l')
    legend.Draw()

    ratioPlot.GetLowerRefYaxis().SetTitle('Data/MC')
    ratioPlot.GetLowerRefGraph().SetLineColor(ROOT.kBlack)
    ratioPlot.GetXaxis().SetTitle(varName)


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
            varName = displayObj.name[idx2_+1:idx3_]
            if (varName not in varNameUsed):
                varNameUsed.append(varName)
                for j in range(i+1, len(listToDisplay)):
                    displayObjToCompare = listToDisplay[j]
                    if varName in displayObjToCompare.name:
                        print(f"match: {displayObj.name} vs {displayObjToCompare.name}")
                        compareDuo = [displayObj, displayObjToCompare]
                        compareDuoPrint = [displayObj.name, displayObjToCompare.name]
                        print(compareDuoPrint)
                        canvasList.append(displayCompare(compareDuo, varName))
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
            elif ('hist' in displayObj.name):
                idx_ = displayObj.name.find('_')
                idx2_ = displayObj.name.find("_", idx_ + 1)
                idx3_= displayObj.name.find('_', idx2_ + 1)
                idxlast_ = displayObj.name.rfind('_')
                particleName = displayObj.name[idx_ + 1:idx2_]
                varName = displayObj.name[idx2_+1:idx3_]
                dataName = displayObj.name[idxlast_:]
                if (varName == 'pT') and ('leading' in displayObj.name):
                    leadingSubleading.append(displayObj.obj)
                    #Only if there are two histograms (leading and subleading pTs) draw them on canvas
                    if (len(leadingSubleading) == 2):
                        canvasList.append(displayPTHist(leadingSubleading, dataName))
                        #Remove the already drawn pair of leading and subleading histograms from list ST other leading/subleading can be drawn together in the correct pairs
                        leadingSubleading = []
                    else: None
                else: 
                    canvasList.append(displayHist(displayObj.obj, varName, dataName))
            #ADD OTHER ELIFS LATER FOR OTHER PLOT IMPLEMENTATION
            else: None
    return canvasList       

#Open the ROOT File and print out content for user
outFile = ROOT.TFile("plots/UL2017/MC.root", "READ")
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
    user_input2 = input("Please select a year to compare or type add to add all data and compare to it\n")
    if user_input2.lower() == 'all':
        #implement summation of data hists and use the MC_Scaled2017Combined hist (need to generate this still in analysis script)
        None
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

canvasList = displayObjects(listToDisplay, compareTag)

#Draw and update the canvasses
for canvas in canvasList:
    canvas.Update()
    canvas.Draw()


outFile.Close()
print("\n Program running... Press Enter to stop.")
input()  # Waits for the Enter key to be pressed
print("Program stopped.")



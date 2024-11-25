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
    hist.Draw()

    return canvas

def displayPTHist(listOfPt, dataName):
    canvas = ROOT.TCanvas('canvas_hist_leadingsubleadingpt_' + dataName, 'Muon pT Distribution', 800, 600)
    i = 0
    for hist in listOfPt:
        if (i==1):
            hist.Draw("same")
        else:
            hist.Draw()
        i+=1
    return canvas

def displayObjects(listToDisplay):
    canvasList = []
    leadingSubleading = []
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
    f"Please specify all objects to be displayed (comma-separated).\nType 'all' to display all loaded objects\nType 'year' to select a year\n"
    )

#Generate and add to canvasList accordingly
inputDisplayList = []
listToDisplay = []
if (user_input1.lower() == 'all'):
    #if all selected use the complete list
    #A list of custom display objects is generated which takes the obj name and obj itself from the .root file
    i = 0
    for key in keys:
        objName = key.GetName()
        listToDisplay.append(displayObj(objName, outFile.Get(objName)))    
        i+=1
elif (user_input1.lower() == 'year'):
    user_input2 = input("Please select a year or type 'MC' to view the simuation plots\n")
    for key in keys:
        objName = key.GetName()
        if user_input2 in objName:
            print(objName)
            listToDisplay.append(displayObj(objName, outFile.Get(objName))) 

else:
    inputDisplayList = [obj.strip() for obj in user_input.split(",")]
    #Based on the list of selected names of objects by user take objects from .root accordingly
    i = 0
    for inputObjName in inputDisplayList:
        listToDisplay.append(displayObj(inputObjName, outFile.Get(inputObjName)))
        i += 1

canvasList = displayObjects(listToDisplay)

#Draw and update the canvasses
for canvas in canvasList:
    canvas.Update()
    canvas.Draw()


outFile.Close()
print("\n Program running... Press Enter to stop.")
input()  # Waits for the Enter key to be pressed
print("Program stopped.")



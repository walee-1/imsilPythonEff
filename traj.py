#! /usr/bin/env python

import numpy as np

import os
import math

import re
import argparse as ap

def binUpFunc(num,binSize):
    return math.ceil(num/binSize)*binSize

def cceFunc(depth,gammaVal,tauVal):
    return(1-gammaVal*math.exp(-depth/tauVal))

def binmids(binList):
    arr=[0]*(len(binList)-1)
    for i in range(len(binList)-1):
        arr[i]=(binList[i]+binList[i+1])/2
    return arr


def cceCalc(binList,enLossTab,gamma,tau):
    cceArray=[0]*(len(enLossTab))
    cceArray=np.asarray(cceArray, dtype=np.float32)
    for i in range(len(enLossTab)):
        cceArray[i]=enLossTab[i]/3.6*(cceFunc(binList[i]/10,gamma,tau)+cceFunc(binList[i+1]/10,gamma,tau))/2
    return(cceArray)

def trajSorterFunc():
    global noIons
    global singleIonData
    global window
    global gammaVals
    global tauVals
    global writeFilecceBest
    global writeFilecceWorse
    global writeFileTraj
    global writeFileHisto
    noIons+=1
    binSize=100
    flag=False
    localArray=np.array(singleIonData)
    singleIonData=[]
    maxDepth=np.amax(localArray,axis=0)[2]
    if maxDepth>window:
        maxBin=binUpFunc(maxDepth-window,binSize)
        if maxBin==maxDepth-window:
            binList=np.arange(0,maxBin+2*binSize,binSize)
        else:
            binList=np.arange(0,maxBin+binSize,binSize)
        enLossTab=[0]*(len(binList)-1)
        enLossTab=np.asarray(enLossTab, dtype=np.float32)
        for i in range(len(localArray)):
            if localArray[i][2]>window:
                tempIndex=math.floor((localArray[i][2]-window)/binSize)
                enLossTab[tempIndex]+=localArray[i][6]
        negPos=np.where(enLossTab<0)
        if len(negPos)>0 and len(negPos[0])>0:
            flag=True
        while flag:
                flag=False

                for i in range(len(negPos[0])):
                    index=negPos[0][i]
                    if index==0:
                        enLossTab[index]=0
                    else:
                        enLossTab[index-1]=enLossTab[index-1]+enLossTab[index]
                        enLossTab[index]=0
                negPos=np.where(enLossTab<0)
                if len(negPos)>0 and len(negPos[0])>0:
                   flag=True
    else:
        enLossTab=np.array([0],dtype=np.float32)
        binList=np.arange(0,binSize+binSize,binSize)

    cceWorse=cceCalc(binList,enLossTab,gammaVals[1],tauVals[1])
    cceBest=cceCalc(binList,enLossTab,gammaVals[0],tauVals[0])
    bins=binmids(binList)
    cceBestWrite=np.vstack((bins,cceBest))
    cceWorseWrite=np.vstack((bins,cceWorse))
    enTabWrite=np.vstack((bins,enLossTab))
    with open(writeFilecceBest,"a") as f:
        np.savetxt(f,cceBestWrite, fmt='%1.3f')
    with open(writeFilecceWorse,"a") as f:
        np.savetxt(f,cceWorseWrite, fmt='%1.3f')
    with open(writeFileTraj,"a") as f:
        np.savetxt(f,enTabWrite, fmt='%1.3f')
    totalCceBest=np.sum(cceBest)
    totalCceWorse=np.sum(cceWorse)
    totalEnLoss=np.sum(enLossTab)
    histoWrite=np.array([totalEnLoss,totalCceBest,totalCceWorse],dtype=np.float32).T
    with open(writeFileHisto,"a") as f:
        np.savetxt(f,histoWrite[np.newaxis], fmt='%1.3f')


def fileCheck(fileName):
    if os.path.exists(fileName):
        id=1
        fileNameWoTXT=re.sub('.txt','',fileName)
        fileName2=''.join([fileNameWoTXT,'_Copy',str(id),'.txt'])
        while os.path.exists(fileName2):
            fileName2=re.sub('_Copy[\d+]\.txt','',fileName2)
            id+=1
            fileName2=''.join([fileName2,'_Copy',str(id),'.txt'])
        os.rename(fileName,fileName2)

parser=ap.ArgumentParser(
    description='''This program goes track by track and writes the output to 4 files and one Help File (Traj_OutputHelp).''',
    epilog='''Please read the output help file for info on output files'''
)
parser.add_argument('--removeFlag',type=bool,help='Set to True if you want to remove the .tra files after processing')
parser.add_argument('--fileName',type=str,required=True,help='Input file name, should be of the format text_runNo_enEN,angANG.tra')
parser.add_argument('--window',type=int,required=True,help='Thickness of entrance window to ignore, please provide in Angstorms')
parser.add_argument('--outputPath',type=str,help='This is just if you need a specific path, otherwise the path is ~/plgad_results/')
parser.add_argument('--windowName',type=str,help='this is for the output file name apart from HONALSI')
args=parser.parse_args()



if args.outputPath:
    outPath=os.path.join(os.path.expanduser('~'),args.outputPath)
else:
    outPath=os.path.join(os.path.expanduser('~'),'plgad_results_python')
inpFile=os.path.splitext(os.path.basename(args.fileName))[0]

if args.removeFlag:
    rmFlag=args.removeFlag
else:
    rmFlag=False

window=args.window

fileSplit=inpFile.split("_")
runNo=fileSplit[1]


en=fileSplit[2]
en=en.strip("en")


ang=fileSplit[3]
ang=ang.strip("ang")

noIons=0
singleIonData=[]

ehPair={}
gammaVals=np.array([0.6,0.9],dtype=np.float32)
tauVals=np.array([50,100],dtype=np.float32)
if args.windowName:
    winName=args.windowName
else:
    winName=fileSplit[0]



outFileName=[winName,'_en',str(en),'_ang',str(ang),'_',str(runNo)]
writeFile = ''.join(outFileName)
writeFilecceBest = os.path.join(outPath,''.join(['cceBest_',writeFile,'.txt']))
writeFilecceWorse = os.path.join(outPath,''.join(['cceWorse_',writeFile,'.txt']))
writeFileTraj = os.path.join(outPath,''.join(['Traj_',writeFile,'.txt']))
writeFileHisto=os.path.join(outPath,''.join(['ToHisto_',writeFile,'.txt']))
fileCheck(writeFilecceBest)
fileCheck(writeFilecceWorse)
fileCheck(writeFileTraj)
fileCheck(writeFileHisto)
helpFile=os.path.join(outPath,"Traj_OutputHelp.txt")


with open(helpFile,'w') as fobj:
    fobj.write("The following files are written: cceBest, cceWorse, Traj and finally ToHisto. ToHisto is the file that is used by efficiency.py to calculate efficiencies\n")
    fobj.write("if a file already exists with the same name (apart from the helpfile), the file is renamed to sameName_Copy(INT).txt and the newer file is written with the desired name\n")
    fobj.write("for cceBest, Worse and Traj, the information for each line is written by using two lines so for 5 ions, 10 lines will be written\n")
    fobj.write("the first line is the middle of the histogramming bins and the second line is the weighted cce e/hpairs or for the case of Traj file the raw energy lost.\n")
    fobj.write("You need to divide by 3.6 for silicon to get the raw e/h pairs generated\n")
    fobj.write("The last and the most important file for calculation of efficiency is: ToHisto...txt, the format is: totalEnLost cceBest cceWorse\n")

prevLine=""

with open(args.fileName,'r') as fobj:
    for line in fobj:
        numbers=[float(num) for num in line.split()]
        if numbers[7]==1.0 and numbers[8]==1.0 and numbers[9]==1.0 and prevLine!="":
            trajSorterFunc()

        prevLine=line
        numbers2=[float(num) for num in prevLine.split()]
        singleIonData.append(numbers2)
trajSorterFunc()

if rmFlag:
    os.remove(args.fileName)
# print(singleIonData)
# print(noIons)
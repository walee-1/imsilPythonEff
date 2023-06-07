#! /usr/bin/env python


import numpy as np

import os
import math
import re

import argparse as ap




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


def binDef(list):
    global binSize
    maxBin=math.ceil(np.max(list)/100)*100
    return(np.arange(0,maxBin,binSize))

def effCalc(list):
    global totalIons
    return(np.cumsum(list)/totalIons)

def binmids(binList):
    arr=[0]*(len(binList)-1)
    for i in range(len(binList)-1):
        arr[i]=(binList[i]+binList[i+1])/2
    return arr

def histoFileWrite(list,fileName):
    global totalIons
    bins=binDef(list)
    histo=np.histogram(list,bins)
    #wPedestal=np.cumsum(histo[0])/totalIons
    woPedestal=[sum(histo[0][i:len(histo[0])])/totalIons for i in range(1,len(histo[0]))]
    woPedestal=np.array(woPedestal)
    #woPedestal=np.cumsum(histo[0][1:])
    DetectedEff=(1-woPedestal)*100
    binMids=binmids(bins)

    pedWrite=np.vstack([binMids[1:],woPedestal])
    detWrite=np.vstack([binMids[1:],DetectedEff])

    with open(fileName,"a") as f:
       np.savetxt(f,detWrite,fmt='%1.3f')
       np.savetxt(f,pedWrite, fmt='%1.3f')
       #np.savetxt(f,histoWr,fmt='%1.3f')
       f.write(str(totalIons))
        



#####MAIN STARTS HERE!~
### Arugments Processing

parser=ap.ArgumentParser(    
    description='''This program goes through the toHisto file and writes 3 output files and one Help File (Track_OutputHelp).''',
    epilog='''Please read the output help file for info on output files'''
    )

parser.add_argument('--path',type=str,help="input Path for the toHisto file, if not given it is ~/plgad_results_python/")
parser.add_argument('--fileName',type=str,required=True,help="input File name, ToHisto_... without runNo and.txt")
parser.add_argument('--noRuns',type=int,required=True,help="Total number of runs to combine")
parser.add_argument('--outputPath',type=str,help="output path, by default:~/plgad_results_python/")
parser.add_argument('--binSize',type=float,help="binSize to histogram data by, by default it is set to 10")
args=parser.parse_args()

if args.path:
    inPath=args.path
else:
    inPath=os.path.join(os.path.expanduser('~'),'plgad_results_python')

if args.binSize:
    binSize=args.binSize
else:
    binSize=10

if args.outputPath:
    outPath=os.path.join(os.path.expanduser('~'),args.outputPath)
else:
    outPath=os.path.join(os.path.expanduser('~'),'plgad_results_python')

inpFile=os.path.basename(args.fileName)
runNos=args.noRuns

outFileName=inpFile.replace("ToHisto_","")

writeFileCceBest= os.path.join(outPath,''.join(['HistocceBest_',outFileName,'.txt']))
writeFileCceWorse= os.path.join(outPath,''.join(['HistocceWorse_',outFileName,'.txt']))
writeFileTotal= os.path.join(outPath,''.join(['HistoTotalEn_',outFileName,'.txt']))
fileCheck(writeFileCceBest)
fileCheck(writeFileCceWorse)
fileCheck(writeFileTotal)

##Main Core of reading file and doing processing
totalEn=[]
cceBest=[]
cceWorse=[]


for i in range(runNos):
    f=os.path.join(outPath,args.fileName+"_"+str(i)+".txt")
    with open(f,'r') as fobj:
        for line in fobj:
            numbers=[float(num) for num in line.split()]
            totalEn.append(numbers[0])
            cceBest.append(numbers[1])
            cceWorse.append(numbers[2])

totalIons=len(totalEn)
helpFile=os.path.join(outPath,"Efficiency_OutputHelp.txt")


with open(helpFile,'w') as fobj:
    fobj.write("The following files are written: cceBest, cceWorse, TotalEnLost\n")
    fobj.write("if a file already exists with the same name (apart from the helpfile), the file is renamed to sameName_Copy(INT).txt and the newer file is written with the desired name\n")
    fobj.write("for cceBest, Worse and EnLost, the information for the histograms is written in 2lines, first is the bins and the second what they contain\n")
    fobj.write("First 2 lines are for Lost efficiency and the second two are for efficiency and the last line is the total number of Ions Total.\n")

histoFileWrite(cceBest,writeFileCceBest)
histoFileWrite(cceWorse,writeFileCceWorse)
histoFileWrite(totalEn,writeFileTotal)





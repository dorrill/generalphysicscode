#This Python script reads in text-based data of ADC-counts from a CAEN digitizer.
#In particular, it is used to analyze data from PMTs connected to the digitizer.
#Both SiPMs and traditional HV PMTs were used. The script can plot events (
#Voltage vs. time) and fit for the compton edges of various isotopes is they were
#used in testing.

#Usage: python3 plotCAENData.py filename

#Imported libraries below
import matplotlib.pyplot as plt
import numpy as np
import sys


#   #   #   #   #   #   #   #   #   # Functions #   #   #   #   #   #   #   #   #   #

def analyzeWaveforms(dataFile):
    #Creates histograms of Max pulse heights per event (max V), pulse width, and the pulse integral
    allADCCounts,allMaxes=findVoltages(dataFile)
    #print(allMaxes)
    newWidths=findWaveformWidths(allADCCounts,allMaxes)
    #print(newWidths)
    pulseintegrals=integrateWaveform(allADCCounts)
    #print(pulseintegrals)
    generalHistogram(allMaxes,26,"Max Pulse Heights (ADC counts)","Pulse Heights","red")
    generalHistogram(newWidths,26,"Pulse Widths (ns)","Pulse Widths @50% Max","cyan")
    generalHistogram(pulseintegrals,25,"Pulse Integrals (ADC-ns)","Pulse Integrals","green")

def findVoltages(dataFile):
    print("Finding ADC Counts, Voltages")
    maxVoltages = []    #stores max value of voltage for each channel of each event
    channelADCCounts = []    #empty array to put voltages in. Refreshes here with each new line / channel
    allChannelADCCounts = []
    badMaxV = 0 #count how many weird values we get close to 0
    for line in dataFile:
        if line[0] != "#" and len(line) < 8:      #skip file header
            line = line.strip() #strip away no numeric characters from data
            zeroed_V = int(line)-baselineoffset #noise floor, baseline is roughly at 9814
            channelADCCounts.append(zeroed_V)

        elif len(line) > 10 and len(channelADCCounts)>0: #when you reach end of event, at header for next one
            allChannelADCCounts.append(channelADCCounts)
            maxV = int(np.amax(channelADCCounts))
            maxVoltages.append(maxV)
            channelADCCounts = []
    return allChannelADCCounts,maxVoltages

def findWaveformWidths(allChannelVoltages,maxVoltages):
    #Finds the FWHM (full-width, at half maximum of the voltage for a pulse - i.e. region of interest)
    print("Finding Widths in ns")
    widths=[]
    sampletime = 4 #4ns
    x1=-1          #used to calculate width
    x2=-1
    for i in range(len(allChannelVoltages)):
        channelVoltages=allChannelVoltages[i]
        max=maxVoltages[i]
        leftSide="yes"
        for j in range(111,190):
            V=channelVoltages[j]
            if V>0.4*max and leftSide=="yes":
                m=(V-channelVoltages[j-1])/1 #slope = deltaV/deltax but delta x = 1
                if m == 0:
                    x1 = j
                else:
                    b=V-(m*j)
                    x1=(V-b)/m
                leftSide="no"
            elif V<0.4*max and leftSide=="no":
                m=(V-channelVoltages[j-1])/1 #slope = deltaV/deltax but delta x = 1
                if m == 0:
                    x2 = j
                else:
                    b=V-(m*j)
                    x2=(V-b)/m
                width = (x2-x1)*sampletime
                #print("x1,x2, width: (",x1,",",x2,",",width,")")
                widths.append(width)
                break
    return widths

def generalHistogram(dataPoints,numBins,axislabel,title,color1): #Will use max voltage of each event as a proxy for energy, then histogram them
    (binVals, bins, patches1) = plt.hist(dataPoints,numBins,color=color1,histtype='bar', ec='black')
    plt.title(title)
    plt.xlabel(axislabel)
    plt.ylabel("#Events per bin")
    plotName = "./"+fileNameStr
    #plt.savefig()
    plt.show()

def generateTimes(numSamples):
    #Generates times for each ADC count value (i.e. time for each Voltage)
    #For Caen digitizer samples taken at 250 MHz --> 4ns between samples, 260 total samples
    times=[]
    sampletime = 4 #4ns
    for i in range(0,numSamples):
        times.append(sampletime*i)
    return times


def histogramEnergiesDownwardPulses(dataFile): #Will use max voltage of each event as a proxy for energy, then histogram them. For downward going, negative Voltages from a HV PMT
    minVoltages = []    #stores max value of voltage for each channel of each event
    channelADCCounts = []    #empty array to put voltages in. Refreshes here with each new line / channel
    badMinV = 0 #count how many weird values we get close to 0
    for line in dataFile:
        if line[0] != "#" and len(line) < 8:      #skip file header
            line = line.strip() #strip away no numeric characters from data
            zeroed_V = 9814-int(line) #noise floor, baseline is roughly at 9814
            channelADCCounts.append(zeroed_V)
        elif len(line) > 10 and len(channelADCCounts)>0: #when you reach end of event, at header for next one
            maxV = int(np.amax(channelADCCounts))
            if maxV < 8000:
                minVoltages.append(maxV)
            else:
                badMinV += 1
            channelADCCounts = []

    print("Length min voltages: ",len(minVoltages))
    print("Bad minV Values (<100): ", badMinV)
    #plt.hist(minVoltages,80,(2500,10000))   #histogram the maximum voltages with some # bins
    (binVals, bins, patches1) = plt.hist(minVoltages,100,(0,4000),color='red')
    #binVals, bins = np.histogram(minVoltages,80,(0,4000),color='red')
    #print("binVals: \n",binVals, "\nBins: \n",bins,"/n"
    #fitNa22Edges(binVals,bins)
    #fitCs137Edges(binVals,bins)
    #fitCo60Edges(binVals,bins)
    #dataTypeString = fileNameStr[3].split("_")
    dataTypeString = "SiPM LED"
    plt.title("%s Events, ADC Counts"%(dataTypeString))
    plt.xlabel("ADC Counts (downward pulses, noise at 9800)")
    plt.ylabel("#Events per bin (100 bins)")
    plotName = "./"+fileNameStr
    #plt.savefig()
    plt.show()

def histogramEnergies(dataFile): #Will use max voltage of each event as a proxy for energy, then histogram them. For positive (upward going) signals
    maxVoltages = []    #stores max value of voltage for each channel of each event
    channelADCCounts = []    #empty array to put voltages in. Refreshes here with each new line / channel
    badMaxV = 0 #count how many weird values we get close to 0
    for line in dataFile:
        if line[0] != "#" and len(line) < 8:      #skip file header
            line = line.strip() #strip away no numeric characters from data
            zeroed_V = int(line)-baselineoffset #noise floor, baseline is roughly at 9814
            #print("Zeroed_V: ",zeroed_V)
            channelADCCounts.append(zeroed_V)
        elif len(line) > 10 and len(channelADCCounts)>0: #when you reach end of event, at header for next one
            #print("Length ADC count list: ",len(channelADCCounts))
            maxV = int(np.amax(channelADCCounts))
            maxVoltages.append(maxV)
            #print minV
            #if maxV < 8000:
            #    minVoltages.append(maxV)
                #print("minV: ", minV)
            #else:
            #    badMinV += 1
            channelADCCounts = []

    print("Length max voltages: ",len(maxVoltages))
    print("Bad maxV Values (<100): ", badMaxV)
    #plt.hist(minVoltages,80,(2500,10000))   #histogram the maximum voltages with some # bins
    (binVals, bins, patches1) = plt.hist(maxVoltages,61,(100,1000),color='red')
    #binVals, bins = np.histogram(minVoltages,80,(0,4000),color='red')
    #print("binVals: \n",binVals, "\nBins: \n",bins,"/n"
    #fitNa22Edges(binVals,bins)
    #fitCs137Edges(binVals,bins)
    #fitCo60Edges(binVals,bins)
    #dataTypeString = fileNameStr[3].split("_")
    dataTypeString = "SiPM LED"
    plt.title("%s Events, ADC Counts"%(dataTypeString))
    plt.xlabel("ADC Counts (minus noise at 9800)")
    plt.ylabel("#Events per bin (100 bins)")
    plotName = "./"+fileNameStr
    #plt.savefig()
    plt.show()

def integrateWaveform(allADCs):
    #"Integrates" by finding the area under the curve from data from all evets from the CAEN digitizer. It's slow.
    #Caen digitizer samples at 250 MHz --> 4ns between samples,260 total samples
    #Integrating by trapezoidal rule 0.5(v1+v2)*(dt)
    print("Integrating Waveforms")
    integrals=[]
    sampletime = 4.0 #4ns
    integral=0
    for i in range(len(allADCs)):
        volts=allADCs[i]
        integral=0
        for j in range(111,190):
            integral+=0.5*sampletime*(volts[j]+volts[j+1])
        integrals.append(integral)
    return integrals

def plotEvents(dataFile):
    #Plot ADC-counts (basically Voltages) from each event on one canvas
    eventNum = 0        #first event will be 0
    for line in dataFile:
        if line[0] == "R":  #R is the first character in a new event, so plot the previous event here
            if eventNum > 0:
                sampletimes=generateTimes(len(channelVoltages))
                plt.plot(sampletimes,channelVoltages)
                plt.title("SiPM Pulses from TestStand (ADC counts)")
                plt.xlabel("ns (250 MHz sampling)")
                plt.ylabel("ADC Counts (minus noise at 9800)")
                plt.show()
                answer = "y"#input("Continue? ") #break plotting
                if answer == "n" or answer == "no" or answer=="No": #allow user to quit
                    break
            eventNum+=1
            channelVoltages = []
        elif len(line)>3 and len(line)<6:      #skip file header
            voltage=int(line)
            channelVoltages.append(voltage-baselineoffset)

#   #   #   #   #   #   #   #   #   # Main #   #   #   #   #   #   #   #   #   #

baselineoffset=7980
file = sys.argv[1]
fileNameStr = str(file)
#fileNameStr = str(file).split('.')
#fileNameStr = str(fileNameStr[1]).split('/')
#prefix = str(fileNameStr[3])
print("Compiling max voltages from ", fileNameStr)
with open(file) as f:
    #histogramEnergies(f)    #Try out histogram out
    #f.seek(0)
    #print("Beginning to plot events!")
    #plotEvents(f)           #Try to plot all 8 channels of an event on one plot
    analyzeWaveforms(f)

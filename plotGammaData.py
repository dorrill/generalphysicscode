#$This Python Script Can Take Generic Text-based data of voltages from a CAEN readout system
#and fit for the compton edges of various isotopes. Although it's tailored to specific data,
#the techniques here can apply to similar data from other systems.

#Imported libraries below
import matplotlib.pyplot as plt
import numpy as np
import sys


#   #   #   #   #   #   #   #   #   # Functions #   #   #   #   #   #   #   #   #   #

#Simple histogram of random x-y points
def testHistogram2d():
    # Generate a normal distribution with center at x=0, y = 1 and plot it.
    N_points = 16000
    x = np.random.randn(N_points)
    y = np.random.randn(N_points) + 8   #y points are offset by 8 from the origin
    plt.hist2d(x, y, bins=(50, 50), cmap=plt.cm.jet)
    plt.show()

    N_points = 320000   #try again with 20x more data points
    x = np.random.randn(N_points)
    y = np.random.randn(N_points) + 1
    plt.hist2d(x, y, bins=(50, 50), cmap=plt.cm.nipy_spectral)  #cmap is the color-map type. Try gist_earth or ocean
    plt.show()

def histogramEnergies(dataFile): #Will use max voltage of each event as a proxy for energy, then histogram them
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

    print "Length min voltages: ",len(minVoltages)
    print "Bad minV Values (<100): ", badMinV
    #plt.hist(minVoltages,80,(2500,10000))   #histogram the maximum voltages with some # bins
    (binVals, bins, patches1) = plt.hist(minVoltages,100,(0,4000),color='red')
    #binVals, bins = np.histogram(minVoltages,80,(0,4000),color='red')
    #print "binVals: \n",binVals, "\nBins: \n",bins,"/n"
    fitNa22Edges(binVals,bins)
    #fitCs137Edges(binVals,bins)
    #fitCo60Edges(binVals,bins)
    dataTypeString = fileNameStr[3].split("_")
    dataTypeString = dataTypeString[1]
    plt.title("%s Events, ADC Counts"%(dataTypeString))
    plt.xlabel("ADC Counts (downward pulses, noise at 9800)")
    plt.ylabel("#Events per bin (100 bins)")
    plotName = "./IIT_lab_data/"+fileNameStr[2]+"/"+prefix+"_hist"
    plt.savefig(plotName)
    plt.show()

def fitNa22Edges(binVals,bins): #works right for 100 bins
    vals511 = []
    vals1275 = []
    peakCts511 = -1
    maxBin511 = -1
    for i in range(23,38):
        vals511.append(binVals[i])
        if binVals[i] > peakCts511:
            peakCts511 = binVals[i]
            maxBin511 = bins[i]
    vals1275 = []
    peakCts1275 = -1
    maxBin1275 = -1
    for i in range(60,75):
        vals1275.append(binVals[i])
        if binVals[i] > peakCts1275:
            peakCts1275 = binVals[i]
            maxBin1275 = bins[i]

    print "511 packets, peak, MeV/ADC: ", peakCts511," ",maxBin511, " ",(maxBin511/0.511)
    print "1275 packets, peak: ", peakCts1275," ",maxBin1275, " ",(maxBin1275/1.275)

def fitCs137Edges(binVals,bins): #works right for 100 bins
    vals662 = []
    peakCts662 = -1
    maxBin662 = -1
    for i in range(23,38):
        vals662.append(binVals[i])
        if binVals[i] > peakCts662:
            peakCts662 = binVals[i]
            maxBin662 = bins[i]

    print "0.662 peakcts, peak, MeV/ADC: ", peakCts662," ",maxBin662, " ",(maxBin662/0.4220)

def fitCs137Edges(binVals,bins): #works right for 100 bins
    vals662 = []
    peakCts662 = -1
    maxBin662 = -1
    for i in range(23,38):
        vals662.append(binVals[i])
        if binVals[i] > peakCts662:
            peakCts662 = binVals[i]
            maxBin662 = bins[i]

    print "0.662 peakcts, peak, MeV/ADC: ", peakCts662," ",maxBin662, " ",(maxBin662/0.4220)

def fitMn54Edges(binVals,bins): #works right for 100 bins
    vals = []
    peakCts = -1
    maxBin = -1
    for i in range(23,38):
        vals.append(binVals[i])
        if binVals[i] > peakCts:
            peakCts = binVals[i]
            maxBin = bins[i]

    print "Peakcts, peak, MeV/ADC: ", peakCts662," ",maxBin662, " ",(maxBin662/0.4220)



def plotEvents(dataFile):
    eventNum = 0        #first event will be 0
    for line in dataFile:
        channelVoltages = []    #empty array to put voltages in. Refreshes here with each new line / channel
        if line[0] != "#":      #skip file header
            line = line.strip() #strip away no numeric characters from data
            line = line.split(" ")  #split data by spaces (each Voltage in the file is separate by a space)

            i = 7               #data column where voltages start
            while i < len(line):    #add everything after column 7 from data file to voltages
                channelVoltages.append(float(line[i]))
                i += 1

            if int(line[1]) == eventNum: #If still the same event, add new channel voltages to the plot
                plt.plot(channelVoltages) #add this channel to voltage plot, but don't show plot yet
            else: #new event, so let's plot it and starts a new plot afterward
                plt.show() #show the previous channels from the previous event
                plt.plot(channelVoltages) #Begin new plot by adding new voltage
                eventNum = int(line[1]) #set eventNumber to new event's number. Must use "int" or it will be characters
                answer = raw_input("Continue? ") #break plotting
                if answer == "n" or answer == "no": #allow user to quit
                    break


#   #   #   #   #   #   #   #   #   # Main #   #   #   #   #   #   #   #   #   #
file = sys.argv[1]
fileNameStr = str(file).split('.')
fileNameStr = str(fileNameStr[1]).split('/')
prefix = str(fileNameStr[3])
print "Compiling max voltages from ", prefix
with open(file) as f:
    histogramEnergies(f)    #Try out histogram out
    #print "Beginning to plot events!"
    #plotEvents(f)           #Try to plot all 8 channels of an event on one plot

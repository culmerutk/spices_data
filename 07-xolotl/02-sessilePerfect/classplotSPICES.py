LATTICE_CONSTANT=0.36
PI = 3.14159265359
VISIBILITY = 2.0
DOSE_RATE = 5.0E-4
DENSITY_CONVERSION = 1.0E27
LABEL_SIZE = 12

class clusterData():
    def __init__(self,fileName,typeName,clusterType,firstSize,
                 plotName,plotStyle,plotColor):
        self.fileName=fileName
        self.typeName=typeName
        self.clusterType=clusterType
        self.firstSize=firstSize
        self.plotName=plotName
        self.plotStyle=plotStyle
        self.plotColor=plotColor
        self.readData()
        #self.deleteFirst()
        #self.negativeNegative()
        self.scaleData()
        self.allAverage()
        self.visibleAverage()
        return

    def readData(self):
        with open(self.fileName, "r") as f:
            rawdata = f.readlines()
        rawdata = [line.split() for line in rawdata]
        timestepcounter = -1
        timestep = []
        time = []
        conc = []
        tempconc = []
        for line in rawdata:
            if line[0] == "time":
                timestepcounter = timestepcounter + 1
                timestep.append(int(timestepcounter))
                time.append(float(line[1]))
                if (timestepcounter > 0):
                    tempconc = sorted(tempconc, key = lambda x: (x[0]))
                    tempconc = [(float(x[1])/100.0) for x in tempconc]
                    conc.append(tempconc)
                    tempconc = []
            else:
                value = [abs(int(line[0])),float(line[2])]
                tempconc.append(value)
        tempconc = sorted(tempconc, key = lambda x: (x[0]))
        tempconc = [(float(x[1])/100.0) for x in tempconc]
        conc.append(tempconc)

        #data = [[float(number) for number in line.split()] for line in rawdata]
        #timestep = [line[0] for line in data]
        #time = [line[1] for line in data]
        #conc = data
        #for line in conc:
        #    del line[1]
        #    del line[0]
        self.timestep=timestep
        self.time=time
        self.concentration=conc
        return

    def scaleData(self):
        self.time = [value*DOSE_RATE for value in self.time]
        self.concentration = [[value*DENSITY_CONVERSION for value in line] for line in self.concentration]
        return

    def deleteFirst(self):
        del self.time[0]
        del self.timestep[0]
        del self.concentration[0]
        return

    def negativeNegative(self):
        for line in self.concentration:
            for value in line:
                if value < 0:
                    value = -value
        return

    def getRadius(self,size):
        if self.typeName == "vType" or self.typeName == "voidType" or self.typeName == "iType":
            return (LATTICE_CONSTANT * (float(size)*3.0/(16.0*PI))**(1.0/3.0))
        elif self.typeName == "perfectType":
            return (LATTICE_CONSTANT * (float(size)*2.0**0.5/(4.0*PI))**0.5)
        else:
            return (LATTICE_CONSTANT * (float(size)*3.0**0.5/(4.0*PI))**0.5)

    def getDiameter(self,size):
        return (2.0 * self.getRadius(size))

    def allAverage(self):
        density = []
        average_size = []
        for line in self.concentration:
            temp_density = 0.0
            temp_average_size = 0.0
            size = self.firstSize
            for value in line:
                temp_density = temp_density + value
                temp_average_size = temp_average_size + value * self.getDiameter(size)
                size = size + 1
            density.append(temp_density)
            if temp_density == 0.0:
                average_size.append(0.0)
            else:
                average_size.append((temp_average_size/temp_density))
        self.averageDensity=density
        self.averageDiameter=average_size
        return

    def visibleAverage(self):
        density = []
        average_size = []
        for line in self.concentration:
            temp_density = 0.0
            temp_average_size = 0.0
            size = self.firstSize
            for value in line:
                diam = self.getDiameter(size)
                if diam > VISIBILITY:
                    temp_density = temp_density + value
                    temp_average_size = temp_average_size + value * diam
                size = size + 1
            density.append(temp_density)
            if temp_density == 0.0:
                average_size.append(0.0)
            else:
                average_size.append((temp_average_size/temp_density))
        self.visibleDensity = density
        self.visibleDiameter = average_size
        return (density,average_size)

def main():
    # Read all cluster data
    data = []
    data.append(clusterData("iType.dat","iType","point",1,
                            "INTERSTITIAL",'-','black'))
    data.append(clusterData("frankType.dat","frankType","loop",5,
                            "FRANK LOOP",'--','black'))
    data.append(clusterData("perfectType.dat","perfectType","loop",5,
                            "PERFECT LOOP",':','black'))
    data.append(clusterData("vType.dat","vType","point",1,
                            "VACANCY",'-','red'))
    data.append(clusterData("faultedType.dat","faultedType","loop",6,
                            "FAULTED LOOP",'--','red'))
    data.append(clusterData("voidType.dat","voidType","void",6,
                            "VOID",':','red'))

    # Prepare for plotting
    import matplotlib.pyplot as plt

    # Plot total concentration
    plt.subplot(231)
    for cluster in data:
        plt.plot(cluster.time,cluster.averageDensity,label=cluster.plotName,
                 linestyle=cluster.plotStyle,color=cluster.plotColor)
    plt.xlabel('DOSE (dpa)',fontsize=LABEL_SIZE)
    plt.ylabel('DENSITY (per cubic meter)',fontsize=LABEL_SIZE)
    plt.yscale('log')
    #plt.xlim(0.0,9.0)
    plt.ylim(1E12, 1E26)
    lgd = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc=3,ncol=1,
                     mode="expand",borderaxespad=0.,fontsize=(LABEL_SIZE-2))

    # Plot average size of all clusters
    plt.subplot(234)
    for cluster in data:
        plt.plot(cluster.time,cluster.averageDiameter,label=cluster.plotName,
                 linestyle=cluster.plotStyle,color=cluster.plotColor)
    plt.xlabel('DOSE (dpa)',fontsize=LABEL_SIZE)
    plt.ylabel('AVERAGE DIAMETER (nm)',fontsize=LABEL_SIZE)
    #plt.xlim(0.0,9.0)

    # Get visible loop data
    loops=[]
    for cluster in data:
        if cluster.clusterType == "loop":
            loops.append(cluster)
    loopDensity=[sum(values) for values in
                  zip(*[cluster.visibleDensity for cluster in loops])]
    loopTime=loops[0].time
    loopDiameter=[((numerator/denominator) if denominator != 0.0 else 0.0) for numerator,denominator in zip([sum(x) for x in zip(*[[(size*density) for size,density in zip(cluster.visibleDiameter,cluster.visibleDensity)]for cluster in loops])],loopDensity)]

    # Plot loop density
    plt.subplot(232)
    plt.plot(loopTime,loopDensity,label="LOOP DENSITY",
             linestyle='-',color='black')
    plt.xlabel('DOSE (dpa)',fontsize=LABEL_SIZE)
    plt.ylabel('DENSITY (per cubic meter)',fontsize=LABEL_SIZE)
    #plt.xlim(0.0,9.0)

    #Plot loop diameter
    plt.subplot(235)
    plt.plot(loopTime,loopDiameter,label="LOOP DIAMETER",
             linestyle='-',color='black')
    plt.xlabel('DOSE (dpa)',fontsize=LABEL_SIZE)
    plt.ylabel('AVERAGE DIAMETER (nm)',fontsize=LABEL_SIZE)
    #plt.xlim(0.0,9.0)

    # Direct calculation to check visible loop diameter
    #i = 0
    #size = len(loops[0].time)
    #dCheck=[]
    #rhoCheck=[]
    #while i < size:
    #    rhoCheck.append(data[1].visibleDensity[i]+data[2].visibleDensity[i]+data[4].visibleDensity[i])
    #    dCheck.append(data[1].visibleDensity[i]*data[1].visibleDiameter[i] +
    #                  data[2].visibleDensity[i]*data[2].visibleDiameter[i] +
    #                  data[4].visibleDensity[i]*data[4].visibleDiameter[i])
    #    dCheck[i] = dCheck[i] / rhoCheck[i] if rhoCheck[i] != 0.0 else 0.0
    #    i=i+1

    # Plot void density
    plt.subplot(233)
    plt.plot(data[5].time,data[5].visibleDensity,label="VOID DENSITY",
             linestyle='-',color='black')
    plt.xlabel('DOSE (dpa)',fontsize=LABEL_SIZE)
    plt.ylabel('DENSITY (per cubic meter)',fontsize=LABEL_SIZE)
    #plt.xlim(0.0,9.0)

    #Plot void diameter
    plt.subplot(236)
    plt.plot(data[5].time,data[5].visibleDiameter,label="VOID DIAMETER",
             linestyle='-',color='black')
    plt.xlabel('DOSE (dpa)',fontsize=LABEL_SIZE)
    plt.ylabel('AVERAGE DIAMETER (nm)',fontsize=LABEL_SIZE)
    #plt.xlim(0.0,9.0)

    plt.tight_layout()
    plt.savefig("SPICES_PLOT.png", format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    #plt.show()
    return

if __name__ == '__main__':
    main()

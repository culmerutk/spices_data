file = open("density.txt",'r')
lines = file.readlines()
size = len(lines)
splitlines=[]
time_density=[]
int_density=[]
frank_density=[]
perfect_density=[]
vac_density=[]
faulted_density=[]
void_density=[]
for index in range(size):
    splitlines.append(lines[index].split())

index = 1
while index < size:
    time_density.append(float(splitlines[index][0])*5E-4)
    int_density.append(float(splitlines[index][1]))
    frank_density.append(float(splitlines[index][2]))
    perfect_density.append(float(splitlines[index][3]))
    vac_density.append(float(splitlines[index][4]))
    faulted_density.append(float(splitlines[index][5]))
    void_density.append(float(splitlines[index][6]))
    index = index + 1

file = open("diameter.txt",'r')
lines = file.readlines()
size = len(lines)
splitlines=[]
time_diameter=[]
frank_diameter=[]
perfect_diameter=[]
faulted_diameter=[]
void_diameter=[]
for index in range(size):
    splitlines.append(lines[index].split())

index = 1
while index < size:
    time_diameter.append(float(splitlines[index][0])*5E-4)
    frank_diameter.append(float(splitlines[index][1]))
    perfect_diameter.append(float(splitlines[index][2]))
    faulted_diameter.append(float(splitlines[index][3]))
    void_diameter.append(float(splitlines[index][4]))
    index = index + 1

file = open("loop.txt",'r')
lines = file.readlines()
size = len(lines)
splitlines=[]
loop_time=[]
loop_density=[]
loop_diameter=[]
for index in range(size):
    splitlines.append(lines[index].split())

index = 1
while index < size:
    loop_time.append(float(splitlines[index][0])*5E-4)
    loop_density.append(float(splitlines[index][1]))
    loop_diameter.append(float(splitlines[index][2]))
    index = index + 1

file = open("void.txt",'r')
lines = file.readlines()
size = len(lines)
splitlines=[]
visible_void_time=[]
visible_void_density=[]
visible_void_diameter=[]
for index in range(size):
    splitlines.append(lines[index].split())

index = 1
while index < size:
    visible_void_time.append(float(splitlines[index][0])*5E-4)
    visible_void_density.append(float(splitlines[index][1]))
    visible_void_diameter.append(float(splitlines[index][2]))
    index = index + 1

import numpy as np
import matplotlib.pyplot as plt

plt.plot(time_density,int_density,label="INTERSTITIAL DENSITY",linewidth=2.0,linestyle='-',color='black')
plt.plot(time_density,frank_density,label="FRANK LOOP DENSITY",linewidth=2.0,linestyle='--',color='black')
plt.plot(time_density,perfect_density,label="PERFECT LOOP DENSITY",linewidth=2.0,linestyle=':',color='black')
plt.plot(time_density,vac_density,label="VACANCY DENSITY",linewidth=2.0,linestyle='-',color='red')
plt.plot(time_density,faulted_density,label="FAULTED LOOP DENSITY",linewidth=2.0,linestyle='--',color='red')
plt.plot(time_density,void_density,label="VOID DENSITY",linewidth=2.0,linestyle=':',color='red')

plt.xlabel('DOSE (dpa)',fontsize=16)
plt.ylabel('DENSITY (per cubic meter)',fontsize=16)
plt.yscale('log')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=14)
plt.savefig("density.pdf", format='pdf', bbox_inches='tight')
plt.clf()

plt.plot(time_diameter,frank_diameter,label="FRANK LOOP DIAMETER",linewidth=2.0,linestyle='--',color='black')
plt.plot(time_diameter,perfect_diameter,label="PERFECT LOOP DIAMETER",linewidth=2.0,linestyle=':',color='black')
plt.plot(time_diameter,faulted_diameter,label="FAULTED LOOP DIAMETER",linewidth=2.0,linestyle='--',color='red')
plt.plot(time_diameter,void_diameter,label="VOID DIAMETER",linewidth=2.0,linestyle=':',color='red')

plt.xlabel('DOSE (dpa)',fontsize=16)
plt.ylabel('DIAMETER (nm)',fontsize=16)
#plt.yscale('log')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=14)
plt.savefig("diameter.pdf", format='pdf', bbox_inches='tight')
plt.clf()

plt.plot(loop_time,loop_density,label="LOOP DENISTY",linewidth=2.0,linestyle='-',color='black')
plt.xlabel('DOSE (dpa)',fontsize=16)
plt.ylabel('DENSITY (per cubic meter)',fontsize=16)
plt.savefig("loop_density.pdf", format='pdf', bbox_inches='tight')
plt.clf()
plt.plot(loop_time,loop_diameter,label="LOOP DIAMETER (nm)",linewidth=2.0,linestyle='-',color='black')
plt.xlabel('DOSE (dpa)',fontsize=16)
plt.ylabel('DIAMETER (nm)',fontsize=16)
plt.savefig("loop_diameter.pdf", format='pdf', bbox_inches='tight')
plt.clf()

plt.plot(visible_void_time,visible_void_density,label="LOOP DENISTY",linewidth=2.0,linestyle='-',color='black')
plt.xlabel('DOSE (dpa)',fontsize=16)
plt.ylabel('DENSITY (per cubic meter)',fontsize=16)
plt.savefig("void_density.pdf", format='pdf', bbox_inches='tight')
plt.clf()
plt.plot(visible_void_time,visible_void_diameter,label="LOOP DIAMETER (nm)",linewidth=2.0,linestyle='-',color='black')
plt.xlabel('DOSE (dpa)',fontsize=16)
plt.ylabel('DIAMETER (nm)',fontsize=16)
plt.savefig("void_diameter.pdf", format='pdf', bbox_inches='tight')

plt.show()

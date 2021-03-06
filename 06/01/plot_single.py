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
    int_density.append(float(splitlines[index][1])*1E25)
    frank_density.append(float(splitlines[index][2])*1E25)
    perfect_density.append(float(splitlines[index][3])*1E25)
    vac_density.append(float(splitlines[index][4])*1E25)
    faulted_density.append(float(splitlines[index][5])*1E25)
    void_density.append(float(splitlines[index][6])*1E25)
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
    frank_diameter.append(float(splitlines[index][1])*2.0)
    perfect_diameter.append(float(splitlines[index][2])*2.0)
    faulted_diameter.append(float(splitlines[index][3])*2.0)
    void_diameter.append(float(splitlines[index][4])*2.0)
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
    loop_density.append(float(splitlines[index][1])*1E25)
    loop_diameter.append(float(splitlines[index][2])*2.0)
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
    visible_void_density.append(float(splitlines[index][1])*1E25)
    visible_void_diameter.append(float(splitlines[index][2])*2.0)
    index = index + 1

fs = 12

import numpy as np
import matplotlib.pyplot as plt

plt.subplot(231)

plt.plot(time_density,int_density,label="INTERSTITIAL",linewidth=2.0,linestyle='-',color='black')
plt.plot(time_density,frank_density,label="FRANK LOOP",linewidth=2.0,linestyle='--',color='black')
plt.plot(time_density,perfect_density,label="PERFECT LOOP",linewidth=2.0,linestyle=':',color='black')
plt.plot(time_density,vac_density,label="VACANCY",linewidth=2.0,linestyle='-',color='red')
plt.plot(time_density,faulted_density,label="FAULTED LOOP",linewidth=2.0,linestyle='--',color='red')
plt.plot(time_density,void_density,label="VOID",linewidth=2.0,linestyle=':',color='red')

plt.xlabel('DOSE (dpa)',fontsize=fs)
plt.ylabel('DENSITY (per cubic meter)',fontsize=fs)
plt.yscale('log')
#lgd = plt.legend(bbox_to_anchor=(0., -0.5, -0.5, 1.), loc=3, ncol=1, mode="expand", borderaxespad=0., fontsize=10)
lgd = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=1, mode="expand", borderaxespad=0., fontsize=10)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=14)
#plt.savefig("density.png", format='png', bbox_inches='tight')

plt.subplot(234)

plt.plot(time_diameter,frank_diameter,label="FRANK LOOP DIAMETER",linewidth=2.0,linestyle='--',color='black')
plt.plot(time_diameter,perfect_diameter,label="PERFECT LOOP DIAMETER",linewidth=2.0,linestyle=':',color='black')
plt.plot(time_diameter,faulted_diameter,label="FAULTED LOOP DIAMETER",linewidth=2.0,linestyle='--',color='red')
plt.plot(time_diameter,void_diameter,label="VOID DIAMETER",linewidth=2.0,linestyle=':',color='red')

plt.xlabel('DOSE (dpa)',fontsize=fs)
plt.ylabel('DIAMETER (nm)',fontsize=fs)
#plt.yscale('log')
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0., fontsize=12)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, fontsize=14)
#plt.savefig("diameter.png", format='png', bbox_inches='tight')

plt.subplot(232)

plt.plot(loop_time,loop_density,label="LOOP DENISTY",linewidth=2.0,linestyle='-',color='black')
plt.xlabel('DOSE (dpa)',fontsize=fs)
plt.ylabel('DENSITY (per cubic meter)',fontsize=fs)
#plt.savefig("loop_density.png", format='png', bbox_inches='tight')
plt.subplot(235)
plt.plot(loop_time,loop_diameter,label="LOOP DIAMETER (nm)",linewidth=2.0,linestyle='-',color='black')
plt.xlabel('DOSE (dpa)',fontsize=fs)
plt.ylabel('DIAMETER (nm)',fontsize=fs)
#plt.savefig("loop_diameter.png", format='png', bbox_inches='tight')

plt.subplot(233)

plt.plot(visible_void_time,visible_void_density,label="LOOP DENISTY",linewidth=2.0,linestyle='-',color='black')
plt.xlabel('DOSE (dpa)',fontsize=fs)
plt.ylabel('DENSITY (per cubic meter)',fontsize=fs)
#plt.savefig("void_density.png", format='png', bbox_inches='tight')
plt.subplot(236)
plt.plot(visible_void_time,visible_void_diameter,label="LOOP DIAMETER (nm)",linewidth=2.0,linestyle='-',color='black')
plt.xlabel('DOSE (dpa)',fontsize=fs)
plt.ylabel('DIAMETER (nm)',fontsize=fs)
#plt.savefig("all.png", format='png', bbox_inches='tight')

plt.tight_layout()

plt.savefig("all.png", format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
#plt.show()

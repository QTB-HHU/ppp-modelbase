'H. A. Berthon, W. A. Bubb, and P. W. Kuchel. 13C n.m.r.                      \
isotopomer and computer-simulation studies of the non-oxidative pen-          \
tose phosphate pathway of human erythrocytes. Biochemical Journal,            \
296(2):379–387, Dec. 1993.'                                                   


import modelbase
import numpy as np
import matplotlib.pyplot as plt
import pppmodel


###############################################################################
'--------------------------------Fig. 2B--------------------------------------'
###############################################################################

# start concentrations

y0d ={'Ru5P':1e-7,
'Xu5P':1e-7,
'Rib5P':1e-7,
'G3P':1e-7,
'Sed7P':1e-7,
'Ery4P':1e-7,
'Fru6P':1e-7,
'DHAP':1e-7,
'Sed1_7P2':1e-7,
'Fru1_6P2':30e-3,
'Man6P':1e-7,
'Glc6P':20e-3,
'Glc1P':1e-7,
'UDPGlc':1e-7,
'Glc1P':1e-7,
'Gal1P':1e-7,
'UDPGal1P':1e-7}

# model call and generation of isotopomer start vector

PPP_model = pppmodel.PPP()
Test = pppmodel.LabelPPP()

y0 = Test.set_initconc_cpd_labelpos(y0d)

y0[Test.cpdIdDict['Glc6P100000']] = 30e-3

Test.par.update(PPP_model.par2) # Berthon1993 parameters, see PPPmodel

# time simulation


begin=0
end=300*60
steps=10000

T = np.linspace(begin,end,steps)

s = modelbase.Simulator(Test)
s.integrator.linear_solver = 'SPGMR'
s.integrator.atol = 1e-15
s.integrator.rtol = 1e-15


import time 
start_time = time.time()
s.timeCourse(T,y0)
print("--- %s seconds ---" % 
(time.time() - start_time))

# plot

plt.plot(s.getT()/60,s.getVarByName('Glc6P100000')*1000, label = 'Glc6P-C1')
plt.plot(s.getT()/60,s.getVarByName('Glc6P001000')*5*1000, label = 'Glc6P-C3')
plt.plot(s.getT()/60,s.getVarByName('Ru5P10000')*5*1000, label ='Rub5P-C1')
plt.xlim(0,17500/60)
plt.xlabel('Time [min]')
plt.ylabel('Concentrations [mM]')
plt.legend()
plt.show()


################################################################################
#'--------------------------------Fig. 3--------------------------------------'
################################################################################
#
# start concentrations

y0d ={'Ru5P':1e-7,
'Xu5P':1e-7,
'Rib5P':1e-7,
'G3P':1e-7,
'Sed7P':1e-7,
'Ery4P':1e-7,
'Fru6P':5e-6,
'DHAP':5e-6,
'Sed1_7P2':1e-7,
'Fru1_6P2':28.5e-3,
'Man6P':1e-7,
'Glc6P':19e-3,
'Glc1P':1e-7}



# model call and generation of isotopomer start vector

PPP_model = pppmodel.PPP()
Test = pppmodel.LabelPPP()

y0 = Test.set_initconc_cpd_labelpos(y0d)

y0[Test.cpdIdDict['Glc6P111111']] = 28.5e-3


Test.par.update(PPP_model.par2) # Berthon1993 parameters, see PPPmodel


# time simulation

begin=0
end=50000
steps=10000

T = np.linspace(begin,end,steps)

s = modelbase.Simulator(Test)
s.integrator.linear_solver = 'SPGMR'
s.integrator.atol = 1e-15
s.integrator.rtol = 1e-15


import time 
start_time = time.time()
s.timeCourse(T,y0)
print("--- %s seconds ---" % 
(time.time() - start_time))


# plot 3a

plt.plot(s.getT()/60,s.getVarByName('Glc6P111111')*1000, label = '[1,2,3,4,5,6-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P000000')*1000, label = 'Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P111000')*1000, label ='[1,2,3-$^{13}C$]-Glc6P')
plt.xlim(0,50000/60)
plt.xlabel('Time [min]')
plt.ylabel('Concentration [mM]')
plt.legend()
plt.show()


#plot 3b

plt.plot(s.getT()/60,s.getVarByName('Glc6P110000')*1000, label = '[1,2-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P001000')*1000, label = '[3-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P001111')*1000, label ='[3,4,5,6-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P000111')*1000, label = '[4,5,6-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P110111')*1000, label ='[1,2,4,5,6-$^{13}C$]-Glc6P')
plt.xlim(0,50000/60)
plt.xlabel('Time [min]')
plt.ylabel('Concentrations [mM]')
plt.legend()
plt.show()


# plot 3c

plt.plot(s.getT()/60,s.getVarByName('Glc6P100000')*1000+s.getVarByName('Glc6P010000')*1000, label = '[1-$^{13}C$]-Glc6P+[2-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P011000')*1000+s.getVarByName('Glc6P101000')*1000, label = '[2,3-$^{13}C$]-Glc6P+[1,3-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P011111')*1000+s.getVarByName('Glc6P101111')*1000, label ='[2,3,4,5,6-$^{13}C$]-Glc6P+[1,3,4,5,6-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P100111')*1000+s.getVarByName('Glc6P010111')*1000, label ='[1,4,5,6-$^{13}C$]-Glc6P+[2,4,5,6-$^{13}C$]-Glc6P')
plt.xlabel('Time [min]')
plt.ylabel('Concentrations [mM]')
plt.xlim(0,50000/60)
plt.legend()
plt.show()

###############################################################################
'--------------------------------Fig. 6--------------------------------------'
###############################################################################

#TA DEFICIT
# start concentration, start vector and parameters, see simulation of Fig. 3
# time simulation

begin=0
end=50000
steps=10000

T = np.linspace(begin,end,steps)

Test.par.update({'TA_Et':0})

s = modelbase.Simulator(Test)
s.integrator.linear_solver = 'SPGMR'
s.integrator.atol = 1e-15
s.integrator.rtol = 1e-15


import time 
start_time = time.time()
s.timeCourse(T,y0)
print("--- %s seconds ---" % 
(time.time() - start_time))

# plot 6a

plt.plot(s.getT()/60,s.getVarByName('Glc6P111111')*1000, label = '[1,2,3,4,5,6-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P000000')*1000, label = 'Glc6')
plt.plot(s.getT()/60,s.getVarByName('Glc6P111000')*1000, label = '[1,2,3-$^{13}C$]-Glc6P')
plt.xlim(0,50000/60)
plt.xlabel('Time [min]')
plt.ylabel('Concentrations [mM]')
plt.legend()
plt.show()

# plot 6b

plt.plot(s.getT()/60,s.getVarByName('Glc6P110000')*1000, label = '[1,2-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P001000')*1000, label = '[3-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P001111')*1000, label = '[3,4,5,6-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P000111')*1000, label = '[4,5,6-$^{13}C$]-Glc6P')
plt.plot(s.getT()/60,s.getVarByName('Glc6P110111')*1000, label = '[1,2,4,5,6-$^{13}C$]-Glc6P')
plt.xlim(0,50000/60)
plt.xlabel('Time [min]')
plt.ylabel('Concentrations [mM]')
plt.legend()
plt.show()

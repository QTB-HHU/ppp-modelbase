###############################################################################
#   author: Tim Nies
#   license: GPL3 (you should have obtained a copy by downloading this package)   
#   name: simulate.py
#
#   methods for simulating the basic model contained in PPPmodel (without modelbase).
#   The modelbase was used to simulate the expanded model see LabelPPPmodel.py
#   and examples/Berthon1993.py
###############################################################################




import numpy as np
import scipy.integrate as scp
import pppmodel

class Simulate(pppmodel.PPP):
        #FUNCTION FOR TIMEINTEGRATION
    def timeintegrator(self,x0,T):
        
        time = T
        results0 = [x0]

        integrator = scp.ode(self.model).set_integrator('lsoda').set_initial_value(x0,0)   
        cnt = 1  


        while cnt < len(time):
            zahl = integrator.integrate(time[cnt])
            results0.append(zahl) 
            cnt += 1
        
        return(results0)
        
        
    #CALCULATES EQUILIBRIUM CONSTANTS 
    @staticmethod  
    def calcKeq(steadystate):
        names=['Ru5PE','RIB5PI','TK1','TK2','TA','TIM','AL1','AL2','MI','GPI','PM','UDPGAL1PT']
        basis=steadystate
        res=[]
        res.append(basis[1]/basis[0]) #Ru5PE
        res.append(basis[0]/basis[2]) #RIB5PI
        res.append((basis[3]*basis[4])/(basis[1]*basis[2])) #TK1
        res.append((basis[3]*basis[6])/(basis[1]*basis[5])) #TK2
        res.append((basis[5]*basis[6])/(basis[4]*basis[3])) #TA
        res.append(basis[7]/basis[3]) #TIM
        res.append((basis[5]*basis[7])/basis[8])#Al1
        res.append((basis[7]*basis[3])/basis[9])#AL2
        res.append(basis[6]/basis[10]) #MI
        res.append(basis[11]/basis[6]) #GPI
        res.append(basis[11]/basis[12])
        
        print('---------------------------------------------'+'\n'+'\n'+'equilibrium constants:'+'\n')
    
        for i in range(len(res)):
            print('Keq_'+ names[i]+': ', res[i])

            
            
            
    #CALCULATES THE STEADY STATE      
    def steadystate(self,x0,toleranz=1e-25):
        labels = ['Ru5P','Xu5P','Rib5P','G3P','Sed7P','Ery4P','Fru6P','DHAP','Sed1_7P2','Fru1_6P2','Man6P','Glc6P','Glc1P','UDPGlc','Gal1P','UDPGal1P']
        time=range(1000000000)
        results0=[x0]
        error=np.linalg.norm(results0[0],ord=2)
        integrator=scp.ode(self.model).set_integrator('lsoda').set_initial_value(x0,0)
        cnt=1
        
        while cnt < len(time) and error > toleranz:
            zahl = integrator.integrate(time[cnt])
            error=np.linalg.norm(zahl-results0[-1],ord=2)
            results0.append(zahl) 
            cnt += 1
            
        print('time point of the steady_state:','\t',integrator.t)
        print('\n','steady state:')
        
        for i in range(len(results0[-1])):
            print('\n','-------------------------','\n',labels[i],':','\t',results0[-1][i])
        
        
        return results0[-1]
        

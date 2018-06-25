#BASIC MODEL FROM MCINTYRE ET AL. (1989)

class PPP(object):
    
    def __init__(self):
        
        #Parameters from McIntyre et al. 1989
        
        self.par= {'RU5PE_k1':3.91e6, 'RU5PE_k2': 438,'RU5PE_k3':305,'RU5PE_k4':1.49e6,'RU5PE_Et':4.22e-6,
                   'RIB5PI_k1':21600,'RIB5PI_k2':14.2,'RIB5PI_k3':33.3,'RIB5PI_k4':60900,'RIB5PI_Et':14.2e-6,
                   'TK1_k1':2.16e5,'TK1_k2':38,'TK1_k3':34,'TK1_k4':1.56e5,'TK1_k5':3.29e5,'TK1_k6':175,'TK1_k7':40,'TK1_k8':4.48e4,'TK1_Et':0.33e-6,
                   'TK2_k1':2.16e5,'TK2_k2':38,'TK2_k3':34,'TK2_k4':1.56e5,'TK2_k5':2.24e6,'TK2_k6':175,'TK2_k7':40,'TK2_k8':2.13e4,'TK2_Et':0.33e-6,
                   'TA_k1':5.8e5,'TA_k2':45.3,'TA_k3':16.3,'TA_k4':1.01e6,'TA_k5':4.9e5,'TA_k6':60,'TA_k7':17,'TA_k8':7.9e4,'TA_Et':0.69e-6,
                   'TIM_k1':3.7e7,'TIM_k2':1320,'TIM_k3':1.46e4,'TIM_k4':1.9e7,'TIM_Et':1.14e-6,
                   'AL2_k1':1.07e7, 'AL2_k2':233, 'AL2_k3':1900, 'AL2_k4':1.12e7, 'AL2_k5':70, 'AL2_k6':6.41e6, 'AL2_Et':0.37e-6,
                   'AL1_k1':8.47e6,'AL1_k2':151,'AL1_k3':117,'AL1_k4':8.43e5,'AL1_k5':70,'AL1_k6':6.41e6,'AL1_Et':0.37e-6,
                   'MI_k1':1.19e6,'MI_k2':800,'MI_k3':800,'MI_k4':1.23e6,'MI_Et':0.628e-9,
                   'GPI_k1':3.98e7,'GPI_k2':1290,'GPI_k3':1550,'GPI_k4':1.56e7,'GPI_Et':0.2e-6,
                   'PM_k1':4.2e6,'PM_k2':80.9,'PM_k3':242.6,'PM_k4':7.2e5,'PM_k5':1.3e7,'PM_k6':10,'PM_Et':0.178e-6,
                   'UDPGal1PT_k1':7.58e5, 'UDPGal1PT_k2':100, 'UDPGal1PT_k3':150, 'UDPGal1PT_k4':8.75e5, 'UDPGal1PT_k5':3.96e5, 'UDPGal1PT_k6':110, 'UDPGal1PT_k7':113, 'UDPGal1PT_k8':8.73e5, 'UDPGal1PT_Et':0.048e-6}
                   
                   
        #Parameters from Berthon et al. 1993
                   
        self.par2= {'RU5PE_k1':3.91e6, 'RU5PE_k2': 438,'RU5PE_k3':305,'RU5PE_k4':1.49e6,'RU5PE_Et':2.96e-6,
                   'RIB5PI_k1':7.24e4,'RIB5PI_k2':3.138e1,'RIB5PI_k3':33.3,'RIB5PI_k4':8.29e4,'RIB5PI_Et':6.9e-5,
                   'TK1_k1':2.16e5,'TK1_k2':5.5e1,'TK1_k3':7.7e1,'TK1_k4':9.0e6,'TK1_k5':3.95e5,'TK1_k6':2.53e2,'TK1_k7':7.7e1,'TK1_k8':1.95e3,'TK1_Et':1.7e-7,
                   'TK2_k1':2.16e5,'TK2_k2':5.5e1,'TK2_k3':7.7e1,'TK2_k4':9.0e6,'TK2_k5':8.96e7,'TK2_k6':2.53e2,'TK2_k7':7.7e1,'TK2_k8':3.08e4,'TK2_Et':1.7e-7,
                   'TA_k1':6.6e5,'TA_k2':45.3,'TA_k3':16.3,'TA_k4':5.0e6,'TA_k5':4.58e6,'TA_k6':60,'TA_k7':17,'TA_k8':1.66e5,'TA_Et':9.1e-7,
                   'TIM_k1':3.7e7,'TIM_k2':1320,'TIM_k3':1.46e4,'TIM_k4':1.9e7,'TIM_Et':7.96e-7,
                   'AL2_k1':1.07e7, 'AL2_k2':233, 'AL2_k3':1900, 'AL2_k4':2.8e7, 'AL2_k5':7.0e1, 'AL2_k6':6.41e6, 'AL2_Et':3.4e-7,
                   'AL1_k1':8.47e6,'AL1_k2':151,'AL1_k3':117,'AL1_k4':4.05e6,'AL1_k5':7.0e1,'AL1_k6':6.41e6,'AL1_Et':3.4e-7,
                   'MI_k1':1.19e6,'MI_k2':800,'MI_k3':800,'MI_k4':1.23e6,'MI_Et':1.7e-10,
                   'GPI_k1':3.98e7,'GPI_k2':1290,'GPI_k3':1550,'GPI_k4':1.56e7,'GPI_Et':1.4e-7,
                   'PM_k1':80.9,'PM_k2':4.2e6,'PM_k3':7.2e5,'PM_k4':2.426e2,'PM_k5':1.0e1,'PM_k6':1.3e7,'PM_Et':1.0e-8,
                   'UDPGal1PT_k1':7.58e5, 'UDPGal1PT_k2':100, 'UDPGal1PT_k3':150, 'UDPGal1PT_k4':8.75e5, 'UDPGal1PT_k5':3.96e5, 'UDPGal1PT_k6':110, 'UDPGal1PT_k7':113, 'UDPGal1PT_k8':8.73e5, 'UDPGal1PT_Et':3.36e-8}


    #rate laws derived by using the King-Altman procedure, see E.L King and C. Altman 1956

    def ratelawA(self,S,P,k1,k2,k3,k4,Et):
        return Et*((k1*k3*S-k2*k4*P)/(k2+k3+k1*S+k4*P))
        
    def ratelawB(self,A,B,P,Q,k1,k2,k3,k4,k5,k6,k7,k8,Et):
        return Et*((k1*k3*k5*k7*A*B-k2*k4*k6*k8*P*Q)/(k1*k3*(k6+k7)*A+k5*k7*(k2+k3)*B+k2*k4*(k6+k7)*P+k6*k8*(k2+k3)*Q+k1*k5*(k3+k7)*A*B+k4*k8*(k2+k6)*P*Q+k5*k8*(k2+k3)*B*Q+k1*k4*(k6+k7)*A*P))
        
    def ratelawC(self,S,P,Q,k1,k2,k3,k4,k5,k6,Et):
        return Et*((k1*k3*k5*S-k2*k4*k6*P*Q)/((k2+k3)*k5+k1*(k3+k5)*S+k2*k4*P+(k2+k3)*k6*Q+k1*k4*S*P+k4*k6*P*Q))
    
    
    #specific rate laws
    
    def Ru5PE(self,S,P):
        return self.ratelawA(S,P,self.par['RU5PE_k1'],self.par['RU5PE_k2'],self.par['RU5PE_k3'],self.par['RU5PE_k4'],self.par['RU5PE_Et'])
    
    def RIB5PI(self,S,P):
        return self.ratelawA(S,P,self.par['RIB5PI_k1'],self.par['RIB5PI_k2'],self.par['RIB5PI_k3'],self.par['RIB5PI_k4'],self.par['RIB5PI_Et'])
    
    def TK1(self,A,B,P,Q):
        return self.ratelawB(A,B,P,Q,self.par['TK1_k1'],self.par['TK1_k2'],self.par['TK1_k3'],self.par['TK1_k4'],self.par['TK1_k5'],self.par['TK1_k6'],self.par['TK1_k7'],self.par['TK1_k8'],self.par['TK1_Et'])
    
    def TK2(self,A,B,P,Q):
        return self.ratelawB(A,B,P,Q,self.par['TK2_k1'],self.par['TK2_k2'],self.par['TK2_k3'],self.par['TK2_k4'],self.par['TK2_k5'],self.par['TK2_k6'],self.par['TK2_k7'],self.par['TK2_k8'],self.par['TK2_Et'])
    
    def TA(self,A,B,P,Q):
        return self.ratelawB(A,B,P,Q,self.par['TA_k1'],self.par['TA_k2'],self.par['TA_k3'],self.par['TA_k4'],self.par['TA_k5'],self.par['TA_k6'],self.par['TA_k7'],self.par['TA_k8'],self.par['TA_Et'])
    
    def TIM(self,S,P):
        return self.ratelawA(S,P,self.par['TIM_k1'],self.par['TIM_k2'],self.par['TIM_k3'],self.par['TIM_k4'],self.par['TIM_Et'])
    
    def AL1(self,S,P,Q):
        return self.ratelawC(S,P,Q,self.par['AL1_k1'],self.par['AL1_k2'],self.par['AL1_k3'],self.par['AL1_k4'],self.par['AL1_k5'],self.par['AL1_k6'],self.par['AL1_Et'])
    
    def AL2(self,S,P,Q):
        return self.ratelawC(S,P,Q,self.par['AL2_k1'],self.par['AL2_k2'],self.par['AL2_k3'],self.par['AL2_k4'],self.par['AL2_k5'],self.par['AL2_k6'],self.par['AL2_Et'])
    
    def MI(self,S,P):
        return self.ratelawA(S,P,self.par['MI_k1'],self.par['MI_k2'],self.par['MI_k3'], self.par['MI_k4'],self.par['MI_Et'])
    
    def GPI(self,S,P):
        return self.ratelawA(S,P,self.par['GPI_k1'],self.par['GPI_k2'],self.par['GPI_k3'], self.par['GPI_k4'],self.par['GPI_Et'])
        
    def PM(self,S,P):
        return self.ratelawA(S,P,self.par['PM_k1'],self.par['PM_k2'],self.par['PM_k3'],self.par['PM_k4'],self.par['PM_Et'])
        
    def UDPGal1PT(self,A,B,P,Q):
        return self.ratelawB(A,B,P,Q,self.par['UDPGal1PT_k1'],self.par['UDPGal1PT_k2'],self.par['UDPGal1PT_k3'],self.par['UDPGal1PT_k4'],self.par['UDPGal1PT_k5'],self.par['UDPGal1PT_k6'],self.par['UDPGal1PT_k7'],self.par['UDPGal1PT_k8'],self.par['UDPGal1PT_Et'])
    
    
    #Model
    
    def model(self,T,x0):
        
        Ru5P=x0[0] 
        Xu5P=x0[1]
        Rib5P=x0[2]
        G3P=x0[3]
        Sed7P=x0[4]
        Ery4P=x0[5]
        Fru6P=x0[6]
        DHAP=x0[7]
        Sed1_7P2=x0[8]
        Fru1_6P2=x0[9]
        Man6P=x0[10]
        Glc6P=x0[11]
        Glc1P=x0[12]
        UDPGlc = x0[13]
        Gal1P =x0[14]
        UDPGal1P = x0[15]
        
        dRu5P=self.RIB5PI(Rib5P,Ru5P)-self.Ru5PE(Ru5P,Xu5P)
        
        dXu5P=self.Ru5PE(Ru5P,Xu5P)-self.TK1(Xu5P,Rib5P,G3P,Sed7P)-self.TK2(Xu5P,Ery4P,G3P,Fru6P)
        
        dRib5P=-self.RIB5PI(Rib5P,Ru5P)-self.TK1(Xu5P,Rib5P,G3P,Sed7P)
        
        dG3P=self.TK1(Xu5P,Rib5P,G3P,Sed7P)+self.AL2(Fru1_6P2,G3P,DHAP)-self.TA(Sed7P,G3P,Ery4P,Fru6P)-self.TIM(G3P,DHAP)+self.TK2(Xu5P,Ery4P,G3P,Fru6P)
        
        dSed7P=self.TK1(Xu5P,Rib5P,G3P,Sed7P)-self.TA(Sed7P,G3P,Ery4P,Fru6P)
        
        dEry4P=self.TA(Sed7P,G3P,Ery4P,Fru6P)+self.AL1(Sed1_7P2,Ery4P,DHAP)-self.TK2(Xu5P,Ery4P,G3P,Fru6P)
        
        dFru6P=self.TA(Sed7P,G3P,Ery4P,Fru6P)+self.MI(Man6P,Fru6P)+self.TK2(Xu5P,Ery4P,G3P,Fru6P)-self.GPI(Fru6P,Glc6P)
        
        dDHAP=self.AL1(Sed1_7P2,Ery4P,DHAP)+self.AL2(Fru1_6P2,G3P,DHAP)+self.TIM(G3P,DHAP)
        
        dSed1_7P2=-self.AL1(Sed1_7P2,Ery4P,DHAP)
        
        dFru1_6P2=-self.AL2(Fru1_6P2,G3P,DHAP)
        
        dMan6P=-self.MI(Man6P,Fru6P)
        
        dGlc6P=self.GPI(Fru6P,Glc6P)+self.PM(Glc1P,Glc6P)
        
        dGlc1P=-self.PM(Glc1P,Glc6P)-self.UDPGal1PT(UDPGlc,Glc1P,Gal1P,UDPGal1P)
        
        dUDPGlc = -self.UDPGal1PT(UDPGlc,Glc1P,Gal1P,UDPGal1P)
        
        dGal1P = self.UDPGal1PT(UDPGlc,Glc1P,Gal1P,UDPGal1P)
        
        dUDPGal1P = self.UDPGal1PT(UDPGlc,Glc1P,Gal1P,UDPGal1P)
        
        return [dRu5P,dXu5P,dRib5P,dG3P,dSed7P,dEry4P,dFru6P,dDHAP,dSed1_7P2,dFru1_6P2,dMan6P,dGlc6P,dGlc1P,dUDPGlc,dGal1P,dUDPGal1P]
    
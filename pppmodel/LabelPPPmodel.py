###############################################################################
#   author: Tim Nies
#   license: GPL3 (you should have obtained a copy by downloading this package)   
#   name: LabelPPPmodel.py
#
#   Reimplemention and expansion of the basic model contained in PPPmodel.py to
#   a system that is able to simulate labelings in metabolic pathways using the
#   modelbase
###############################################################################


import modelbase
import pppmodel



#BEGIN OF THE LABELLING SIMULATIONS USING THE PACKAGE modelbase.     
        
class LabelPPP(modelbase.LabelModel):
    
    def __init__(self):
    
        super().__init__()
        
        #GET THE ORIGINAL PARAMETER        
        basic_PPPmodel = pppmodel.PPP()
        self.origPars = basic_PPPmodel.par.copy()
        
        #UPDATE PARAMETER OF THE LABELMODEL
        self.par.update(self.origPars)
        
        # ADD BASIC COMPOUNDS
        self.add_base_cpd('Ru5P',5)
        self.add_base_cpd('Xu5P',5)
        self.add_base_cpd('Rib5P',5)
        self.add_base_cpd('G3P',3)
        self.add_base_cpd('Sed7P',7)
        self.add_base_cpd('Ery4P',4)
        self.add_base_cpd('Fru6P',6)
        self.add_base_cpd('DHAP',3)
        self.add_base_cpd('Sed1_7P2',7)
        self.add_base_cpd('Fru1_6P2',6)
        self.add_base_cpd('Man6P',6)
        self.add_base_cpd('Glc6P',6)
        self.add_base_cpd('Glc1P',6)

        
        
        
        # GENERAL RATE EQUATIONS
        def ratelawAf(S,Stot,Ptot,k1,k2,k3,k4,Et):
            return Et*((k1*k3*S)/(k2+k3+k1*Stot+k4*Ptot))
            
        def ratelawAr(P,Stot,Ptot,k1,k2,k3,k4,Et):
            return Et*((k2*k4*P)/(k2+k3+k1*Stot+k4*Ptot))
        
        def ratelawBf(A,B,Atot,Btot,Ptot,Qtot,k1,k2,k3,k4,k5,k6,k7,k8,Et):
            return Et*((k1*k3*k5*k7*A*B)/(k1*k3*(k6+k7)*Atot+k5*k7*(k2+k3)*Btot+k2*k4*(k6+k7)*Ptot+k6*k8*(k2+k3)*Qtot+k1*k5*(k3+k7)*Atot*Btot+k4*k8*(k2+k6)*Ptot*Qtot+k5*k8*(k2+k3)*Btot*Qtot+k1*k4*(k6+k7)*Atot*Ptot))
            
        def ratelawBr(P,Q,Atot,Btot,Ptot,Qtot,k1,k2,k3,k4,k5,k6,k7,k8,Et):
            return Et*((k2*k4*k6*k8*P*Q)/(k1*k3*(k6+k7)*Atot+k5*k7*(k2+k3)*Btot+k2*k4*(k6+k7)*Ptot+k6*k8*(k2+k3)*Qtot+k1*k5*(k3+k7)*Atot*Btot+k4*k8*(k2+k6)*Ptot*Qtot+k5*k8*(k2+k3)*Btot*Qtot+k1*k4*(k6+k7)*Atot*Ptot))
        
        def ratelawCf(S,Stot,Ptot,Qtot,k1,k2,k3,k4,k5,k6,Et):
            return Et*((k1*k3*k5*S)/((k2+k3)*k5+k1*(k3+k5)*Stot+k2*k4*Ptot+(k2+k3)*k6*Qtot+k1*k4*Stot*Ptot+k4*k6*Ptot*Qtot))
            
        def ratelawCr(P,Q,Stot,Ptot,Qtot,k1,k2,k3,k4,k5,k6,Et):
            return Et*((k2*k4*k6*P*Q)/((k2+k3)*k5+k1*(k3+k5)*Stot+k2*k4*Ptot+(k2+k3)*k6*Qtot+k1*k4*Stot*Ptot+k4*k6*Ptot*Qtot))
    
    
    
    
    
        # SPECIFIC RATE EQUATIONS
    
    
        def Ru5PEf(p,S,Stot,Ptot):
            return ratelawAf(S,Stot,Ptot,p.RU5PE_k1,p.RU5PE_k2,p.RU5PE_k3,p.RU5PE_k4,p.RU5PE_Et)
            
        self.add_carbonmap_reaction('Ru5PEf',Ru5PEf,[0,1,2,3,4],['Ru5P'],['Xu5P'],'Ru5P','Ru5P_total','Xu5P_total')
        
        
        
        def Ru5PEr(p,P,Stot,Ptot):
            return ratelawAr(P,Stot,Ptot,p.RU5PE_k1,p.RU5PE_k2,p.RU5PE_k3,p.RU5PE_k4,p.RU5PE_Et)
            
        self.add_carbonmap_reaction('Ru5PEr',Ru5PEr,[0,1,2,3,4],['Xu5P'],['Ru5P'],'Xu5P','Ru5P_total','Xu5P_total')
                
        
        
        

        def RIB5PIf(p,S,Stot,Ptot):
            return ratelawAf(S,Stot,Ptot,p.RIB5PI_k1,p.RIB5PI_k2,p.RIB5PI_k3,p.RIB5PI_k4,p.RIB5PI_Et)
            
        self.add_carbonmap_reaction('RIB5PIf',RIB5PIf,[0,1,2,3,4],['Rib5P'],['Ru5P'],'Rib5P','Rib5P_total','Ru5P_total')
        
        
        
        def RIB5PIr(p,P,Stot,Ptot):
            return ratelawAr(P,Stot,Ptot,p.RIB5PI_k1,p.RIB5PI_k2,p.RIB5PI_k3,p.RIB5PI_k4,p.RIB5PI_Et)
            
        self.add_carbonmap_reaction('RIB5PIr',RIB5PIr,[0,1,2,3,4],['Ru5P'],['Rib5P'],'Ru5P','Rib5P_total','Ru5P_total')
        
        
        
        
        
        def TK1f(p,A,B,Atot,Btot,Ptot,Qtot):
            return ratelawBf(A,B,Atot,Btot,Ptot,Qtot,p.TK1_k1,p.TK1_k2,p.TK1_k3,p.TK1_k4,p.TK1_k5,p.TK1_k6,p.TK1_k7,p.TK1_k8,p.TK1_Et)
            
        self.add_carbonmap_reaction('TK1f',TK1f,[2,3,4,0,1,5,6,7,8,9],['Xu5P','Rib5P'],['G3P','Sed7P'],'Xu5P','Rib5P','Xu5P_total','Rib5P_total','G3P_total','Sed7P_total')
        
        
        
        def TK1r(p,P,Q,Atot,Btot,Ptot,Qtot):
            return ratelawBr(P,Q,Atot,Btot,Ptot,Qtot,p.TK1_k1,p.TK1_k2,p.TK1_k3,p.TK1_k4,p.TK1_k5,p.TK1_k6,p.TK1_k7,p.TK1_k8,p.TK1_Et)
            
        self.add_carbonmap_reaction('TK1r',TK1r,[3,4,0,1,2,5,6,7,8,9],['G3P','Sed7P'],['Xu5P','Rib5P'],'G3P','Sed7P','Xu5P_total','Rib5P_total','G3P_total','Sed7P_total')
        
        
        
        
        
        def TK2f(p,A,B,Atot,Btot,Ptot,Qtot):
            return ratelawBf(A,B,Atot,Btot,Ptot,Qtot,p.TK2_k1,p.TK2_k2,p.TK2_k3,p.TK2_k4,p.TK2_k5,p.TK2_k6,p.TK2_k7,p.TK2_k8,p.TK2_Et)
            
        self.add_carbonmap_reaction('TK2f',TK2f,[2,3,4,0,1,5,6,7,8],['Xu5P','Ery4P'],['G3P','Fru6P'],'Xu5P','Ery4P','Xu5P_total','Ery4P_total','G3P_total','Fru6P_total')
        
        
        
        def TK2r(p,P,Q,Atot,Btot,Ptot,Qtot):
            return ratelawBr(P,Q,Atot,Btot,Ptot,Qtot,p.TK2_k1,p.TK2_k2,p.TK2_k3,p.TK2_k4,p.TK2_k5,p.TK2_k6,p.TK2_k7,p.TK2_k8,p.TK2_Et)
            
        self.add_carbonmap_reaction('TK2r',TK2r,[3,4,0,1,2,5,6,7,8],['G3P','Fru6P'],['Xu5P','Ery4P'],'G3P','Fru6P','Xu5P_total','Ery4P_total','G3P_total','Fru6P_total')

        
        


        def TAf(p,A,B,Atot,Btot,Ptot,Qtot):
            return ratelawBf(A,B,Atot,Btot,Ptot,Qtot,p.TA_k1,p.TA_k2,p.TA_k3,p.TA_k4,p.TA_k5,p.TA_k6,p.TA_k7,p.TA_k8,p.TA_Et)
        
        self.add_carbonmap_reaction('TAf',TAf,[3,4,5,6,0,1,2,7,8,9],['Sed7P','G3P'],['Ery4P','Fru6P'],'Sed7P','G3P','Sed7P_total','G3P_total','Ery4P_total','Fru6P_total')    
        
        
        
        def TAr(p,P,Q,Atot,Btot,Ptot,Qtot):
            return ratelawBr(P,Q,Atot,Btot,Ptot,Qtot,p.TA_k1,p.TA_k2,p.TA_k3,p.TA_k4,p.TA_k5,p.TA_k6,p.TA_k7,p.TA_k8,p.TA_Et)
        
        self.add_carbonmap_reaction('TAr',TAr,[4,5,6,0,1,2,3,7,8,9],['Ery4P','Fru6P'],['Sed7P','G3P'],'Ery4P','Fru6P','Sed7P_total','G3P_total','Ery4P_total','Fru6P_total')    
       
        
        
        
        
        def TIMf(p,S,Stot,Ptot):
            return ratelawAf(S,Stot,Ptot,p.TIM_k1,p.TIM_k2,p.TIM_k3,p.TIM_k4,p.TIM_Et)
            
        self.add_carbonmap_reaction('TIMf',TIMf,[2,1,0],['G3P'],['DHAP'],'G3P','G3P_total','DHAP_total')
        
        
        
        def TIMr(p,P,Stot,Ptot):
            return ratelawAr(P,Stot,Ptot,p.TIM_k1,p.TIM_k2,p.TIM_k3,p.TIM_k4,p.TIM_Et)
            
        self.add_carbonmap_reaction('TIMr',TIMr,[2,1,0],['DHAP'],['G3P'],'DHAP','G3P_total','DHAP_total')
        
        
        
        
        
        def AL1f(p,S,Stot,Ptot,Qtot):
            return ratelawCf(S,Stot,Ptot,Qtot,p.AL1_k1,p.AL1_k2,p.AL1_k3,p.AL1_k4,p.AL1_k5,p.AL1_k6,p.AL1_Et)
        
        self.add_carbonmap_reaction('AL1f',AL1f,[3,4,5,6,0,1,2],['Sed1_7P2'],['Ery4P','DHAP'],'Sed1_7P2','Sed1_7P2_total','Ery4P_total','DHAP_total')    
        
        
        
        def AL1r(p,P,Q,Stot,Ptot,Qtot):
            return ratelawCr(P,Q,Stot,Ptot,Qtot,p.AL1_k1,p.AL1_k2,p.AL1_k3,p.AL1_k4,p.AL1_k5,p.AL1_k6,p.AL1_Et)
        
        self.add_carbonmap_reaction('AL1r',AL1r,[4,5,6,0,1,2,3],['Ery4P','DHAP'],['Sed1_7P2'],'Ery4P','DHAP','Sed1_7P2_total','Ery4P_total','DHAP_total')  
        
        
        
        
        
        def AL2f(p,S,Stot,Ptot,Qtot):
            return ratelawCf(S,Stot,Ptot,Qtot,p.AL2_k1,p.AL2_k2,p.AL2_k3,p.AL2_k4,p.AL2_k5,p.AL2_k6,p.AL2_Et)
        
        self.add_carbonmap_reaction('AL2f',AL2f,[3,4,5,0,1,2],['Fru1_6P2'],['G3P','DHAP'],'Fru1_6P2','Fru1_6P2_total','G3P_total','DHAP_total')
        
        
        
        def AL2r(p,P,Q,Stot,Ptot,Qtot):
            return ratelawCr(P,Q,Stot,Ptot,Qtot,p.AL2_k1,p.AL2_k2,p.AL2_k3,p.AL2_k4,p.AL2_k5,p.AL2_k6,p.AL2_Et)
        
        self.add_carbonmap_reaction('AL2r',AL2r,[3,4,5,0,1,2],['G3P','DHAP'],['Fru1_6P2'],'G3P','DHAP','Fru1_6P2_total','G3P_total','DHAP_total')
          
          
          
          
          
        def MIf(p,S,Stot,Ptot):
            return ratelawAf(S,Stot,Ptot,p.MI_k1,p.MI_k2,p.MI_k3, p.MI_k4,p.MI_Et)
            
        self.add_carbonmap_reaction('MIf',MIf,[0,1,2,3,4,5],['Man6P'],['Fru6P'],'Man6P','Man6P_total','Fru6P_total')
        
        
        
        def MIr(p,P,Stot,Ptot):
            return ratelawAr(P,Stot,Ptot,p.MI_k1,p.MI_k2,p.MI_k3, p.MI_k4,p.MI_Et)
            
        self.add_carbonmap_reaction('MIr',MIr,[0,1,2,3,4,5],['Fru6P'],['Man6P'],'Fru6P','Man6P_total','Fru6P_total')
        
        
        
        
        
        def GPIf(p,S,Stot,Ptot):
            return ratelawAf(S,Stot,Ptot,p.GPI_k1,p.GPI_k2,p.GPI_k3, p.GPI_k4,p.GPI_Et)
            
        self.add_carbonmap_reaction('GPIf',GPIf,[0,1,2,3,4,5],['Fru6P'],['Glc6P'],'Fru6P','Fru6P_total','Glc6P_total')
        
        
                
        def GPIr(p,P,Stot,Ptot):
            return ratelawAr(P,Stot,Ptot,p.GPI_k1,p.GPI_k2,p.GPI_k3, p.GPI_k4,p.GPI_Et)
            
        self.add_carbonmap_reaction('GPIr',GPIr,[0,1,2,3,4,5],['Glc6P'],['Fru6P'],'Glc6P','Fru6P_total','Glc6P_total')
        
        
        
            
                
        def PMf(p,S,Stot,Ptot):
            return ratelawAf(S,Stot,Ptot,p.PM_k1,p.PM_k2,p.PM_k3,p.PM_k4,p.PM_Et)
            
        self.add_carbonmap_reaction('PMf',PMf,[0,1,2,3,4,5],['Glc1P'],['Glc6P'],'Glc1P','Glc1P_total','Glc6P_total')
        
        
        
        def PMr(p,P,Stot,Ptot):
            return ratelawAr(P,Stot,Ptot,p.PM_k1,p.PM_k2,p.PM_k3,p.PM_k4,p.PM_Et)
            
        self.add_carbonmap_reaction('PMr',PMr,[0,1,2,3,4,5],['Glc6P'],['Glc1P'],'Glc6P','Glc1P_total','Glc6P_total')


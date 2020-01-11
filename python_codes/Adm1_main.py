from  scipy  import*
from  scipy.integrate import odeint
import numpy as np
import load_parameters
import Adm1_ode
d=load_parameters.load_parameters()
locals().update(d)


t= linspace(0,200 ,2001) # ( t0= 0, tf=200, pas=0.1)

def main (t):

     Xinit=np.array([S_su, S_aa, S_fa, S_va,S_bu, S_pro, S_ac, S_h2,
                    S_ch4, S_IC, S_IN, S_I,X_xc, X_ch, X_pr, X_li,
                    X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I,
                    S_cat, S_an, S_hva, S_hbu,S_hpro, S_hac, S_hco3,
                    S_nh3,S_gas_h2, S_gas_ch4, S_gas_co2, Q_D,T_D,
                    S_D1_D, S_D2_D, S_D3_D, X_D4_D, X_D5_D])

     x= odeint (Adm1_ode.test_ode, Xinit , t)
     return(x)
x=main(t)

     
#--- Post computation:--------------------------------------------------------------------------------------------------

K_w = (pow(10,(-14)))*np.exp((55900/(R*100))*(1/T_base-1/T_op))  # 2.08e-14
p_gas_h2o = 0.0313*np.exp(5290*(1/T_base-1/T_op))  #//0.0557 at 35°C
pH=zeros(2001)
Stot=zeros(2001)
Xbact=zeros(2001)
p_gas_h2=zeros(2001)
p_gas_ch4=zeros(2001)
p_gas_co2=zeros(2001)
P_gas=zeros(2001)
q_gas=zeros(2001)
T_COD=zeros(2001)
VFA=zeros(2001)
IN_N=zeros(2001)
ch4_per=zeros(2001)

#----------Computation of PH:--------------------------------------------------------------------------------------------
Phi = x[:,24] + (x[:,10]-x[:,31]) - x[:,30] - x[:,29]/64 - x[:,28]/112 - x[:,27]/160 - x[:,26]/208 - x[:,25]
x[:,40] = - Phi/2 + sqrt(pow(Phi, 2) + 4*K_w)/2
pH[:]=-log10(x[:,40])
x[:,41]=pH[:]


              

#----------Computations of transfer rates, gas flow rate for the figures : ----------------------------------------------

Stot[:]=x[:,0]+x[:,1]+x[:,2]+x[:,3]+x[:,4]+x[:,5]+x[:,6]+x[:,7]+x[:,8]+x[:,9]+x[:,10]+x[:,11] #substrats sum
Xbact[:]=x[:,16]+x[:,17]+x[:,18]+x[:,19]+x[:,20]+x[:,21]+x[:,22]  #Micro-organisms sum
T_COD[:]=Stot[:]+x[:,12]+x[:,13]+x[:,14]+x[:,15]+x[:,23]+Xbact[:]
VFA[:]=x[:,3]+x[:,4]+x[:,5]+x[:,6]
IN_N[:]=x[:,10]

            

#----------GAS PHASE : --------------------------------------------------------------------------------------------------
            
p_gas_h2[:]= x[:,32]*R*T_op/16
x[:,35]=p_gas_h2[:]
p_gas_ch4[:] = x[:,33]*R*T_op/64 
x[:,36]=p_gas_ch4[:]
p_gas_co2[:] = x[:,34]*R*T_op
x[:,37]=p_gas_co2[:]

            
P_gas[:] = p_gas_h2[:] + p_gas_ch4[:] + p_gas_co2[:] + p_gas_h2o 
x[:,38]=P_gas[:]

#----% CH4 : ------------------------------------------------------------------------------------------------------------

ch4_per[:]=(x[:,36]/x[:,38])*100

#---- Simplified gas calculation Batstone 20002 : ------------------------------------------------------------------------

q_gas[:] = k_P*(P_gas[:]-P_atm)*P_gas[:]/P_atm

#q_gas=max(q_gas,zeros(size(q_gas)))  
x[:,39]= q_gas[:]

# output_variables :----------------------------------------------

output_var=np.array([x[2000,0],x[2000,1],x[2000,2],x[2000,3],x[2000,4],\
                     x[2000,5],x[2000,6],x[2000,7],x[2000,8],x[2000,9],\
                     x[2000,10],x[2000,11],x[2000,12],x[2000,13],x[2000,14],\
                     x[2000,15],x[2000,16],x[2000,17],x[2000,18],x[2000,19],\
                     x[2000,20],x[2000,21],x[2000,22],x[2000,23],x[2000,24],\
                     x[2000,25],x[2000,26],x[2000,27],x[2000,28],x[2000,29],\
                     x[2000,30],x[2000,31],x[2000,32],x[2000,33],x[2000,34],\
                     x[2000,35],x[2000,36],x[2000,37],x[2000,38],x[2000,39],\
                     x[2000,40],x[2000,41]])


#--- Data saving:---------------
np.savetxt('X.txt',x)        
np.savetxt('data_to_plot.txt', (pH, T_COD, VFA, IN_N, Xbact, ch4_per, q_gas))
np.savetxt('output_variables.txt',output_var)

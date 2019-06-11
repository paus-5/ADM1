from  scipy  import*
from  scipy.integrate import odeint
import numpy as np
import  matplotlib.pyplot  as  plt

def ADM1(t):



    #// stoechimetric parameter
    f_sI_xc = 0 #0.1;
    f_xI_xc = 0.53#0.2;
    f_ch_xc = 0.41#0.2;
    f_pr_xc = 0.06#0.2;
    f_li_xc = 0#0.3; 

    N_xc = 0.0376/14
    N_I = 0.06/14
    N_aa = 0.007
    C_xc = 0.02786
    C_sI = 0.03
    C_ch = 0.0313
    C_pr = 0.03
    C_li = 0.022
    C_xI = 0.03
    C_su = 0.0313
    C_aa = 0.03
    f_fa_li = 0.95
    C_fa = 0.0217
    f_h2_su = 0.19
    f_bu_su = 0.13
    f_pro_su = 0.27
    f_ac_su = 0.41
    N_bac = 0.08/14
    C_bu = 0.025
    C_pro = 0.0268
    C_ac = 0.0313
    C_bac = 0.0313
    Y_su = 0.1
    f_h2_aa = 0.06
    f_va_aa = 0.23
    f_bu_aa = 0.26
    f_pro_aa = 0.05
    f_ac_aa = 0.40
    C_va = 0.024
    Y_aa = 0.08
    Y_fa = 0.06
    Y_c4 = 0.06
    Y_pro = 0.04
    C_ch4 = 0.0156
    Y_ac = 0.05
    Y_h2 = 0.06
    
    #// Biochemical parameter:
    k_dis = 0.5
    k_hyd_ch = 10
    k_hyd_pr = 10
    k_hyd_li = 10
    K_S_IN = 1e-4
    k_m_su = 30
    K_S_su = 0.5
    pH_UL_aa = 5.5
    pH_LL_aa = 4
    k_m_aa = 50
    K_S_aa = 0.3
    k_m_fa = 6
    K_S_fa = 0.4
    K_Ih2_fa = 5e-6
    k_m_c4 = 20
    K_S_c4 = 0.2
    K_Ih2_c4 = 1e-5
    k_m_pro = 13
    K_S_pro = 0.1
    K_Ih2_pro = 3.5e-6
    k_m_ac = 8
    K_S_ac = 0.15
    K_I_nh3 = 0.0018
    pH_UL_ac = 7
    pH_LL_ac = 6
    k_m_h2 = 35
    K_S_h2 = 7e-6
    pH_UL_h2 = 6
    pH_LL_h2 = 5
    k_dec_Xsu = 0.02
    k_dec_Xaa = 0.02
    k_dec_Xfa = 0.02
    k_dec_Xc4 = 0.02
    k_dec_Xpro = 0.02
    k_dec_Xac = 0.02
    k_dec_Xh2 = 0.02
    
    #//temperature
    T_base = 298.15
    T_op = 308.15

    #// physiochemical parameter, corrected in temperature
    R = 0.083145
    K_w = (pow(10,(-14)))*np.exp((55900/(R*100))*(1/T_base-1/T_op))  #// 2.08e-14
    K_a_va = pow(10,-4.86)
    K_a_bu = pow(10,-4.82)
    K_a_pro = pow(10,-4.88)
    K_a_ac =  pow(10,-4.76)
    K_a_co2 = pow(10,(-6.35))*np.exp((7646/(R*100))*(1/T_base-1/T_op))  #// 4.94e-7 at 35°C;
    K_a_IN = pow(10,(-9.25))*np.exp((51965/(R*100))*(1/T_base-1/T_op))  #//1.11e-9 at 35°C;
    k_A_Bva = 1e10    #//1e8; according to STR 
    k_A_Bbu = 1e10     #//1e8; according to STR  
    k_A_Bpro = 1e10    #//1e8; according to STR 
    k_A_Bac = 1e10     #//1e8; according to STR 
    k_A_Bco2 = 1e10    #//1e8; according to STR 
    k_A_BIN = 1e10     #//1e8; according to STR 
    P_atm = 1.013
    p_gas_h2o = 0.0313*np.exp(5290*(1/T_base-1/T_op))  #//0.0557 at 35°C
    kLa = 200
    K_H_co2 = 0.035*np.exp((-19410/(R*100))*(1/T_base-1/T_op)) #//0.0271 at 35°C;
    K_H_ch4 = 0.0014*np.exp((-14240/(R*100))*(1/T_base-1/T_op)) #//0.00116 at 35°C;
    K_H_h2 = 7.8e-4*np.exp((-4180/(R*100))*(1/T_base-1/T_op)) #//7.38e-4 at 35°C;

    #// Physical parameter
    V_liq =3400 #1400; //3400
    V_gas =300  #100;  // 300
    k_P = 5e4

    #///////////////////////////////////////////////////////////////////////
    #///////////////////// Input characterisation //////////////////////////
    #///////////////////////////////////////////////////////////////////////
    u=np.array([0,0,0,0,0,0,0,0,0,0.04,0.01,0,248,0,0,0,0,0,0,0,0,0,0,0,0.04,0.02,170,35])
    q_in=u[26]

    #//////////////////////////////////////////////////////////////////////
    #//////////////////////// initial conditions
    #//////////////////////////////////////////////////////////////////////
    S_su =  0.024309
    S_aa =  0.010808
    S_fa =  0.29533
    S_va =  0.02329
    S_bu = 0.031123
    S_pro = 0.043974
    S_ac =  0.50765
    S_h2 =  4.9652e-007
    S_ch4 = 0.055598
    S_IC = 0.10258
    S_IN = 0.10373
    S_I = 3.2327
    X_xc = 7.5567
    X_ch = 0.074679
    X_pr = 0.074679
    X_li = 0.11202
    X_su = 0.57565
    X_aa = 0.43307
    X_fa = 0.44433
    X_c4 = 0.18404
    X_pro = 0.087261
    X_ac = 0.57682
    X_h2 =  0.28774
    X_I =  18.6685
    S_cat =  3.3531e-042
    S_an =  1.5293e-041
    S_hva = 0.023204
    S_hbu = 0.031017
    S_hpro =  0.043803
    S_hac =  0.50616
    S_hco3 = 0.092928
    S_nh3 = 0.0021958
    S_gas_h2 = 1.9096e-005
    S_gas_ch4 = 1.5103
    S_gas_co2 = 0.013766
    Q_D = 166.5395
    T_D = 35
    S_D1_D = 1
    S_D2_D = 1
    S_D3_D = 1
    X_D4_D = 106.0099
    X_D5_D = 106.0099
    S_H_ion = 5.3469e-008

    Xinit=np.array([S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN, S_I, X_xc, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I, S_cat,S_an, S_hva, S_hbu, S_hpro, S_hac, S_hco3, S_nh3, S_gas_h2, S_gas_ch4, S_gas_co2, Q_D, T_D, S_D1_D, S_D2_D, S_D3_D, X_D4_D, X_D5_D])
    
    #//////////////////////////////////////////////////////////////////////
    #/////////////////////////// Solver ODE///////////////////////////////
    #//////////////////////////////////////////////////////////////////////

    import ODE_adm1
    x=odeint(ODE_adm1,Xinit,t)
    print(x)


    
    #//////////////////////////////////////////////////////////////////////// 
    #//////////////pH reconstruction from electoneutrality
    #////////////////////////////////////////////////////////////////////////

    stot=zeros(2001)
    Xbact=zeros(2001)
    p_gas_h2=zeros(2001)
    p_gas_ch4=zeros(2001)
    p_gas_co2=zeros(2001)
    P_gas=zeros(2001)

    for i in range (0,2001):

        Phi = x[i][24] + (x[i][10]-x[i][31]) - x[i][30] - x[i][29]/64 - x[i][28]/112 - x[i][27]/160 - x[i][26]/208 - x[i][25]
        x[i][40] = - Phi/2 + sqrt(pow(Phi, 2) + 4*K_w)/2
        pH=-log10(x[i][40])
        print("pH ==== ",pH)
        x[i][41]=pH #% calculate pH from S_H+
          
          
        #------------ Computations of transfer rates, gas flow rate for the figures : 

        Stot[i]=x[i][0]+x[i][1]+x[i][2]+x[i][3]+x[i][4]+x[i][5]+x[i][6]+x[i][7]+x[i][8]+x[i][9]+x[i][10]+x[i][11] #substrats sum
        Xbact[i]=x[i][16]+x[i][17]+x[i][18]+x[i][19]+x[i][20]+x[i][21]+x[i][22]  #Micro-organisms sum 

        #----------GAS PHASE :
        
        p_gas_h2[i]= x[i][32]*R*T_op/16
        x[i][35]=p_gas_h2[i]
        p_gas_ch4[i] = x[i][33]*R*T_op/64 
        x[i][36]=p_gas_ch4[i]
        p_gas_co2[i] = x[i][34]*R*T_op
        x[i][37]=p_gas_co2[i]

        
        P_gas[i] = p_gas_h2[i] + p_gas_ch4[i] + p_gas_co2[i] + p_gas_h2o 
        x[i][38]=P_gas[i]
            
        ##---// simplified gas calculation Batstone 2002 :

        q_gas = k_P*(P_gas[i]-P_atm)*P_gas[i]/P_atm
        q_gas=np.max(q_gas,zeros(1,size((q_gas),2)))  
        x[i][39]= q_gas

        plt . plot (t , x[t][41])
        plt . show ()



    return 0

t=np.linspace(0,200,2001)
x= ADM1(t)





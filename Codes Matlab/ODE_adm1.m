%//______________________________________________________________________________
%// Modelling the anaerobic digestion with ADM1
%//
%// ODE Implementation of ADM1.
%// Copyright Christian Rosen, Ulf Jeppsson, IEA, Lund University, Lund, Sweden
%// www.iea.lth.se
%
%// Translate into scilab by Cyrille Charnier,2014, 
%//BioEnTech, 74 avenue Paul Sabatier, 11100 Narbonne, France % and % INRA-LBE, Avenue des Etangs, 11100 Narbonne, France
%// cyrille.charnier@bioentech.eu % cyrille.charnier@supagro.inra.fr
%// laurent.lardon@bioentech.eu   % eric.latrille@supagro.inra.fr
%//______________________________________________________________________________
function dx=ODE_adm1(x,t); 

%// filter for negative value
for i= 1:42
   if x(i) < 0 
      x(i) = 0;
end
end
%%AJOUTS:
R = 0.083145;
T_op = 308.15;
T_base = 298.15;
p_gas_h2o = 0.0313*exp(5290*(1/T_base-1/T_op));  %0.0557 at 35°C;
k_P = 5e4;
P_atm = 1.013;
K_w = (10^(-14))*exp((55900/(R*100))*(1/T_base-1/T_op));  % 2.08e-14;
pH_UL_aa = 5.5;
pH_LL_aa = 4;
pH_UL_ac = 7;
pH_LL_ac = 6;
pH_UL_h2 = 6;
pH_LL_h2 = 5;
K_S_IN = 1e-4;
K_Ih2_fa = 5e-6;
K_Ih2_c4 = 1e-5;
K_Ih2_pro = 3.5e-6;
K_I_nh3 = 0.0018;
k_dis = 0.5;
k_hyd_ch = 10;
k_hyd_pr = 10;
k_hyd_li = 10;
k_m_su = 30;
K_S_su = 0.5;
k_m_aa = 50;
K_S_aa = 0.3;
k_m_fa = 6;
K_S_fa = 0.4;
k_m_c4 = 20;
K_S_c4 = 0.2;
k_m_pro = 13;
K_S_pro = 0.1;
k_m_ac = 8;
K_S_ac = 0.15;
k_m_h2 = 35;
K_S_h2 = 7e-6;
k_dec_Xsu = 0.02;
k_dec_Xaa = 0.02;
k_dec_Xfa = 0.02;
k_dec_Xc4 = 0.02;
k_dec_Xpro = 0.02;
k_dec_Xac = 0.02;
k_dec_Xh2 = 0.02;
K_a_va = 10^-4.86;
K_a_bu = 10^-4.82;
K_a_pro = 10^-4.88;
K_a_ac =  10^-4.76;
K_a_co2 = 10^(-6.35)*exp((7646/(R*100))*(1/T_base-1/T_op));  % 4.94e-7 at 35°C;
K_a_IN = 10^(-9.25)*exp((51965/(R*100))*(1/T_base-1/T_op));  % 1.11e-9 at 35°C;
k_A_Bva = 1e10;     %1e8; according to STR ;
k_A_Bbu = 1e10;     %1e8; according to STR  
k_A_Bpro = 1e10;    %1e8; according to STR 
k_A_Bac = 1e10;     %1e8; according to STR 
k_A_Bco2 = 1e10;    %1e8; according to STR 
k_A_BIN = 1e10;     %1e8; according to STR 
kLa = 200;
K_H_co2 = 0.035*exp((-19410/(R*100))*(1/T_base-1/T_op)); %0.0271 at 35°C;
K_H_ch4 = 0.0014*exp((-14240/(R*100))*(1/T_base-1/T_op)); %0.00116 at 35°C;
K_H_h2 = 7.8e-4*exp((-4180/(R*100))*(1/T_base-1/T_op)); %7.38e-4 at 35°C;
C_xc = 0.02786;
C_sI = 0.03;
C_ch = 0.0313;
C_pr = 0.03;
C_li = 0.022;
C_xI = 0.03;
C_su = 0.0313;
C_aa = 0.03;
f_fa_li = 0.95;
C_fa = 0.0217;
f_h2_su = 0.19;
f_bu_su = 0.13;
f_pro_su = 0.27;
f_ac_su = 0.41;
N_bac = 0.08/14;
C_bu = 0.025;
C_pro = 0.0268;
C_ac = 0.0313;
C_bac = 0.0313;
Y_su = 0.1;
f_h2_aa = 0.06;
f_va_aa = 0.23;
f_bu_aa = 0.26;
f_pro_aa = 0.05;
f_ac_aa = 0.40;
C_va = 0.024;
Y_aa = 0.08;
Y_fa = 0.06;
Y_c4 = 0.06;
Y_pro = 0.04;
C_ch4 = 0.0156;
Y_ac = 0.05;
Y_h2 = 0.06;
f_sI_xc = 0;%0.1;
f_xI_xc = 0.53;%0.2;
f_ch_xc = 0.41;%0.2;
f_pr_xc = 0.06;%0.2;
f_li_xc = 0;%0.3; 
N_xc = 0.0376/14;
N_I = 0.06/14;
N_aa = 0.007;
V_liq =3400; %1400; //3400
V_gas =300;% 100;  // 300
u=[0;0;0;0;0;0;0;0;0;0.04;0.01;0;248;0;0;0;0;0;0;0;0;0;0;0;0.04;0.02;170;35];




%// pression
p_gas_h2 = x(33)*R*T_op/16;
p_gas_ch4 = x(34)*R*T_op/64;
p_gas_co2 = x(35)*R*T_op;
  P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o;
%  // simplified gas flow calculation, Batstone 2002
  q_gas = k_P*(P_gas-P_atm);  %*P_gas/P_atm;


%// pH algebric calculation
  phi = x(25)+(x(11)-x(32))-x(31)-x(30)/64-x(29)/112-x(28)/160-x(27)/208-x(26);
  S_H_ion = -phi*0.5+0.5*sqrt(phi*phi+4*K_w);
  if S_H_ion<=0 then
      S_H_ion=10^-12;
  end
  pH_op = -log10(S_H_ion); %// pH


%//// STRs function => ADM1 officiel
%//if (pH_op < pH_UL_aa) then
%//   I_pH_aa = exp(-3.0*(pH_op-pH_UL_aa)*(pH_op-pH_UL_aa)/((pH_UL_aa-pH_LL_aa)*(pH_UL_aa-pH_LL_aa)));
%//else
%//   I_pH_aa = 1.0;
%//   end
%//if (pH_op < pH_UL_ac) then
%//   I_pH_ac = exp(-3.0*(pH_op-pH_UL_ac)*(pH_op-pH_UL_ac)/((pH_UL_ac-pH_LL_ac)*(pH_UL_ac-pH_LL_ac)));
%//else
%//   I_pH_ac = 1.0;
%//   end
%//if (pH_op < pH_UL_h2) then
%//   I_pH_h2 = exp(-3.0*(pH_op-pH_UL_h2)*(pH_op-pH_UL_h2)/((pH_UL_h2-pH_LL_h2)*(pH_UL_h2-pH_LL_h2)));
%//else
%//   I_pH_h2 = 1.0;
%//   end
%
%
%// Hill function on pH
%//pHLim_aa = 10^-((pH_UL_aa + pH_LL_aa)/2);
%//pHLim_ac = 10^-((pH_UL_ac + pH_LL_ac)/2);
%//pHLim_h2 = 10^-((pH_UL_h2 + pH_LL_h2)/2);
%//I_pH_aa = (pH_op^25)/((pH_op^25)+(pHLim_aa^25));
%//I_pH_ac = (pH_op^60)/((pH_op^60)+(pHLim_ac^60));
%//I_pH_h2 = (pH_op^60)/((pH_op^60)+(pHLim_h2^60));

eps=0.000001; %// ADM1 2002;
%// MNPs function 
%//a_aa = 25/(pH_UL_aa-pH_LL_aa+eps);
%//a_ac = 25/(pH_UL_ac-pH_LL_ac+eps);
%//a_h2 = 25/(pH_UL_h2-pH_LL_h2+eps);
%//
%//I_pH_aa = 0.5*(1+tanh(a_aa*(pH_op/pHLim_aa - 1)));
%//I_pH_ac = 0.5*(1+tanh(a_ac*(pH_op/pHLim_ac - 1)));
%//I_pH_h2 = 0.5*(1+tanh(a_h2*(pH_op/pHLim_h2 - 1)));
%//
%
%// Hill function on SH+ => ulf version 2008
pHLim_aa = 10^(-(pH_UL_aa + pH_LL_aa)/2);
pHLim_ac = 10^(-(pH_UL_ac + pH_LL_ac)/2);
pHLim_h2 = 10^(-(pH_UL_h2 + pH_LL_h2)/2);
n_aa=3/(pH_UL_aa-pH_LL_aa);
n_ac=3/(pH_UL_ac-pH_LL_ac);
n_h2=3/(pH_UL_h2-pH_LL_h2);
I_pH_aa = (pHLim_aa^n_aa)/((S_H_ion^n_aa)+(pHLim_aa^n_aa));
I_pH_ac = (pHLim_ac^n_ac)/((S_H_ion^n_ac)+(pHLim_ac^n_ac));
I_pH_h2 = (pHLim_h2^n_h2)/((S_H_ion^n_h2)+(pHLim_h2^n_h2));


%// Inhibition
I_IN_lim = 1/(1+K_S_IN/x(11));
I_h2_fa = 1/(1+x(8)/K_Ih2_fa);
I_h2_c4 = 1/(1+x(8)/K_Ih2_c4);
I_h2_pro = 1/(1+x(8)/K_Ih2_pro);
I_nh3 = 1/(1+x(32)/K_I_nh3);

inhib(1) = I_pH_aa*I_IN_lim;
inhib(2) = inhib(1)*I_h2_fa;
inhib(3) = inhib(1)*I_h2_c4;
inhib(4) = inhib(1)*I_h2_pro;
inhib(5) = I_pH_ac*I_IN_lim*I_nh3;
inhib(6) = I_pH_h2*I_IN_lim;

%// Biochemical Process rates
proc1 = k_dis*x(13);
proc2 = k_hyd_ch*x(14);
proc3 = k_hyd_pr*x(15);
proc4 = k_hyd_li*x(16);
proc5 = (k_m_su*x(1)/(K_S_su+x(1)))*x(17)*inhib(1);
proc6 = (k_m_aa*x(2)/(K_S_aa+x(2)))*x(18)*inhib(1);
proc7 = (k_m_fa*x(3)/(K_S_fa+x(3)))*x(19)*inhib(2);
proc8 = (k_m_c4*x(4)/(K_S_c4+x(4)))*x(20)*(x(4)/(x(4)+x(5)+eps))*inhib(3);
proc9 = (k_m_c4*x(5)/(K_S_c4+x(5)))*x(20)*(x(5)/(x(4)+x(5)+eps))*inhib(3);
proc10 = (k_m_pro*x(6)/(K_S_pro+x(6)))*x(21)*inhib(4);
proc11 = (k_m_ac*x(7)/(K_S_ac+x(7)))*x(22)*inhib(5);
proc12 = (k_m_h2*x(8)/(K_S_h2+x(8)))*x(23)*inhib(6);
proc13 = k_dec_Xsu*x(17);
proc14 = k_dec_Xaa*x(18);
proc15 = k_dec_Xfa*x(19);
proc16 = k_dec_Xc4*x(20);
proc17 = k_dec_Xpro*x(21);
proc18 = k_dec_Xac*x(22);
proc19 = k_dec_Xh2*x(23);

%// Acid-base rates
procA4 = k_A_Bva*(x(27)*(K_a_va+S_H_ion)-K_a_va*x(4));
procA5 = k_A_Bbu*(x(28)*(K_a_bu+S_H_ion)-K_a_bu*x(5));
procA6 = k_A_Bpro*(x(29)*(K_a_pro+S_H_ion)-K_a_pro*x(6));
procA7 = k_A_Bac*(x(30)*(K_a_ac+S_H_ion)-K_a_ac*x(7));
procA10 = k_A_Bco2*(x(31)*(K_a_co2+S_H_ion)-K_a_co2*x(10));
procA11 = k_A_BIN*(x(32)*(K_a_IN+S_H_ion)-K_a_IN*x(11));

%// Gas transfer rates
procT8 = kLa*(x(8)-16*K_H_h2*p_gas_h2);
procT9 = kLa*(x(9)-64*K_H_ch4*p_gas_ch4);
procT10 = kLa*((x(10)-x(31))-K_H_co2*p_gas_co2);

%// reaction parameters
stoich1 = -C_xc + f_sI_xc*C_sI + f_ch_xc*C_ch + f_pr_xc*C_pr + f_li_xc*C_li + f_xI_xc*C_xI;
stoich2 = -C_ch+C_su;
stoich3 = -C_pr+C_aa;
stoich4 = -C_li + (1-f_fa_li)*C_su + f_fa_li*C_fa;
stoich5 = -C_su + (1-Y_su)*(f_bu_su*C_bu+f_pro_su*C_pro+f_ac_su*C_ac) + Y_su*C_bac;
stoich6 = -C_aa + (1-Y_aa)*(f_va_aa*C_va+f_bu_aa*C_bu+f_pro_aa*C_pro+f_ac_aa*C_ac) + Y_aa*C_bac;
stoich7 = -C_fa + (1-Y_fa)*0.7*C_ac + Y_fa*C_bac;
stoich8 = -C_va + (1-Y_c4)*0.54*C_pro+(1-Y_c4)*0.31*C_ac+Y_c4*C_bac;
stoich9 = -C_bu + (1-Y_c4)*0.8*C_ac+Y_c4*C_bac;
stoich10 = -C_pro + (1-Y_pro)*0.57*C_ac+Y_pro*C_bac;
stoich11 = -C_ac + (1-Y_ac)*C_ch4+Y_ac*C_bac;
stoich12 = (1-Y_h2)*C_ch4 + Y_h2*C_bac;
stoich13 = -C_bac + C_xc;

reac1 = proc2+(1-f_fa_li)*proc4-proc5;
reac2 = proc3-proc6;
reac3 = f_fa_li*proc4-proc7;
reac4 = (1-Y_aa)*f_va_aa*proc6-proc8;
reac5 = (1-Y_su)*f_bu_su*proc5+(1-Y_aa)*f_bu_aa*proc6-proc9;
reac6 = (1-Y_su)*f_pro_su*proc5+(1-Y_aa)*f_pro_aa*proc6+(1-Y_c4)*0.54*proc8-proc10;
reac7 = (1-Y_su)*f_ac_su*proc5+(1-Y_aa)*f_ac_aa*proc6+(1-Y_fa)*0.7*proc7+(1-Y_c4)*0.31*proc8+(1-Y_c4)*0.8*proc9+(1-Y_pro)*0.57*proc10-proc11;
reac8 = (1-Y_su)*f_h2_su*proc5+(1-Y_aa)*f_h2_aa*proc6+(1-Y_fa)*0.3*proc7+(1-Y_c4)*0.15*proc8+(1-Y_c4)*0.2*proc9+(1-Y_pro)*0.43*proc10-proc12-procT8;
reac9 = (1-Y_ac)*proc11+(1-Y_h2)*proc12-procT9;
reac10 = - stoich1*proc1 - stoich2*proc2 - stoich3*proc3 - stoich4*proc4 - stoich5*proc5 - stoich6*proc6 - stoich7*proc7 - stoich8*proc8 - stoich9*proc9 - stoich10*proc10 - stoich11*proc11 - stoich12*proc12 - stoich13*proc13 - stoich13*proc14 - stoich13*proc15 - stoich13*proc16 - stoich13*proc17 - stoich13*proc18 - stoich13*proc19 - procT10;
reac11 = (N_xc-f_xI_xc*N_I-f_sI_xc*N_I-f_pr_xc*N_aa)*proc1 - Y_su*N_bac*proc5 + (N_aa-Y_aa*N_bac)*proc6 - Y_fa*N_bac*proc7 - Y_c4*N_bac*proc8 - Y_c4*N_bac*proc9 - Y_pro*N_bac*proc10 - Y_ac*N_bac*proc11 - Y_h2*N_bac*proc12 + (N_bac-N_xc)*(proc13+proc14+proc15+proc16+proc17+proc18+proc19);
reac12 = f_sI_xc*proc1;
reac13 = -proc1+proc13+proc14+proc15+proc16+proc17+proc18+proc19;
reac14 = f_ch_xc*proc1-proc2;
reac15 = f_pr_xc*proc1-proc3;
reac16 = f_li_xc*proc1-proc4;
reac17 = Y_su*proc5-proc13;
reac18 = Y_aa*proc6-proc14;
reac19 = Y_fa*proc7-proc15;
reac20 = Y_c4*proc8+Y_c4*proc9-proc16;
reac21 = Y_pro*proc10-proc17;
reac22 = Y_ac*proc11-proc18;
reac23 = Y_h2*proc12-proc19;
reac24 = f_xI_xc*proc1;

%// Gas flow calculation
%//q_gas = R*T_op*V_liq*(procT8/16+procT9/64+procT10)/(P_atm-p_gas_h2o);
  if q_gas < 0 then
    q_gas = 0.0
end

%////////// differential equation ///////////////////////////////////////////////
   
dx(1) = 1/V_liq*(u(27)*(u(1)-x(1)))+reac1;  %// Ssu
dx(2) = 1/V_liq*(u(27)*(u(2)-x(2)))+reac2;  %// Saa
dx(3) = 1/V_liq*(u(27)*(u(3)-x(3)))+reac3;  %// Sfa
dx(4) = 1/V_liq*(u(27)*(u(4)-x(4)))+reac4; %// Sva */
dx(5) = 1/V_liq*(u(27)*(u(5)-x(5)))+reac5; %// Sbu */
dx(6) = 1/V_liq*(u(27)*(u(6)-x(6)))+reac6;  %// Spro */
dx(7) = 1/V_liq*(u(27)*(u(7)-x(7)))+reac7; %// Sac */
dx(8) = 1/V_liq*(u(27)*(u(8)-x(8)))+reac8; %// Sh2
dx(9) = 1/V_liq*(u(27)*(u(9)-x(9)))+reac9; %// Sch4
dx(10) = 1/V_liq*(u(27)*(u(10)-x(10)))+reac10; %// SIC */
dx(11) = 1/V_liq*(u(27)*(u(11)-x(11)))+reac11; %// SIN */
dx(12) = 1/V_liq*(u(27)*(u(12)-x(12)))+reac12; %// SI
dx(13) = 1/V_liq*(u(27)*(u(13)-x(13)))+reac13; %//xc
dx(14) = 1/V_liq*(u(27)*(u(14)-x(14)))+reac14; %//Xch
dx(15) = 1/V_liq*(u(27)*(u(15)-x(15)))+reac15; %//Xpr
dx(16) = 1/V_liq*(u(27)*(u(16)-x(16)))+reac16; %//Xli
dx(17) = 1/V_liq*(u(27)*(u(17)-x(17)))+reac17; %//Xsu
dx(18) = 1/V_liq*(u(27)*(u(18)-x(18)))+reac18; %//Xaa
dx(19) = 1/V_liq*(u(27)*(u(19)-x(19)))+reac19; %//Xfa
dx(20) = 1/V_liq*(u(27)*(u(20)-x(20)))+reac20; %//Xc4
dx(21) = 1/V_liq*(u(27)*(u(21)-x(21)))+reac21; %//Xpro
dx(22) = 1/V_liq*(u(27)*(u(22)-x(22)))+reac22; %//Xac
dx(23) = 1/V_liq*(u(27)*(u(23)-x(23)))+reac23; %//Xh2
dx(24) = 1/V_liq*(u(27)*(u(24)-x(24)))+reac24; %//XI

dx(25) = 1.0/V_liq*(u(27)*(u(25)-x(25))); %// Scat+ */
dx(26) = 1.0/V_liq*(u(27)*(u(26)-x(26))); %// San- */

dx(27) = -procA4;  %// Sva- */
dx(28) = -procA5;  %// Sbu- */
dx(29) = -procA6;  %// Spro- */
dx(30) = -procA7;  %// Sac- */
dx(31) = -procA10; %// SHCO3- */
dx(32) = -procA11; %// SNH3 */

dx(33) = -x(33)*q_gas/V_gas+procT8*V_liq/V_gas; %// Sgas H2
dx(34) = -x(34)*q_gas/V_gas+procT9*V_liq/V_gas; %// Sgas CH4
dx(35) = -x(35)*q_gas/V_gas+procT10*V_liq/V_gas; %//Sgas CO2


%// Dummy states*/
dx(36) = 0; 
dx(37) = 0; 
dx(38) = 0; 
dx(39) = 0;
dx(40) = 0;
dx(41) = 0;
dx(42) = 0;


endfunction

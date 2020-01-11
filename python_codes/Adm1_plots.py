import  matplotlib.pyplot  as  plt
import numpy as np
from  scipy  import*


t=linspace(0,200,2001)


#--- PH :--------------------------------------

plt.tick_params(axis = 'both', labelsize = 6)
M = np.loadtxt('data_to_plot.txt')
plt.xlabel(" time (d)", size=7) 
plt.plot(t,M[0,:], '-')
plt.title('pH',size=7)
plt.savefig("Plot_pH.png")
plt.close()

#--- DCO :--------------------------------------

plt.tick_params(axis = 'both', labelsize = 6)
plt.xlabel(" time (d)", size=7)
plt.ylabel(" kg.m$^{-3}$", size=7) 
plt.plot(t,M[1,:], '-')
plt.title('Total COD (kg.m$^{-3}$)', size=7)
plt.savefig("Total COD.png")
plt.close()

#--- VFA :--------------------------------------

plt.tick_params(axis = 'both', labelsize = 6)
plt.xlabel(" time (d)", size=7)
plt.ylabel(" kg COD.m$^{-3}$", size=7) 
plt.plot(t,M[2,:], '-')
plt.title('Volatil Fatty Acid (kg COD.m$^{-3}$) ', size=7)
plt.savefig("Volatil Fatty Acid.png")
plt.close()

#--- IN_N :--------------------------------------

plt.tick_params(axis = 'both', labelsize = 6)
plt.xlabel(" time (d)", size=7)
plt.ylabel(" kmole N.m$^{-3}$", size=7) 
plt.plot(t,M[3,:], '-')
plt.title('Inorganic Nitrogen (kmole N.m$^{-3}$)', size=7)
plt.savefig("Inorganic Nitrogen.png")
plt.close()

#--- Bacterial_pop :--------------------------------------

plt.tick_params(axis = 'both', labelsize = 6)
plt.xlabel(" time (d)", size=7)
plt.ylabel(" kg COD.m$^{-3}$", size=7) 
plt.plot(t,M[4,:], '-')
plt.title('Bacterial_population (kg COD.m$^{-3}$)', size=7)
plt.savefig("X_bact.png")
plt.close()

#--- %CH4:--------------------------------------

plt.tick_params(axis = 'both', labelsize = 6)
plt.xlabel(" time (d)", size=7) 
plt.plot(t,M[5,:], '-')
plt.title('% CH4', size=7)
plt.savefig("% CH4.png")
plt.close()

#--- Gas flow:--------------------------------------

plt.tick_params(axis = 'both', labelsize = 6)
plt.xlabel(" time (d)", size=7)
plt.ylabel(" Nm$^{3}$.d$^{-1}$", size=7) 
plt.plot(t,M[6,:], '-')
plt.title('Gas_flow_rate (Nm$^{3}$.d$^{-1}$) ', size=7)
plt.savefig("Gas_flow_rate.png")



'''
*Si vous souhaitez afficher tous les graphes sur la meme figure
 --------------------------------------------------------------

plt.figure

#--- Total_COD:----------------------------------------------------------------
plt.gcf().subplots_adjust(left = 0.2, bottom = 0.2, right = 0.9, top = 0.9, wspace = 0.7, hspace = 0.9)
plt.subplot(3,2,1)
plt.tick_params(axis = 'both', labelsize = 7)
plt.plot(t,M[1,:],color='red')
#plt.ylabel("T_COD")
plt.title('Total COD (kg.m$^{-3}$)', size=7)

#--- VFA : ---------------------------------------------
plt.subplot(3,2,2)
plt.tick_params(axis = 'both', labelsize = 7)
plt.plot(t,M[2,:],color='green')
plt.title('Volatil Fatty Acid (kg COD.m$^{-3}$) ', size=7)

#---- Inorganic_nitrogen : ------------------------------
plt.subplot(3,2,4)
plt.tick_params(axis = 'both', labelsize = 7)
plt.plot(t,M[3,:],color='yellow')
plt.title('Inorganic Nitrogen (kmole N.m$^{-3}$)', size=7)

#---- Gas flow rate :------------------------------------
plt.subplot(3,2,7)
plt.tick_params(axis = 'both', labelsize = 7)
plt.plot(t,M[6,:],color='black')
plt.xlabel('time (d)', size=7)
plt.title('Gas_flow_rate (Nm$^{3}$.d$^{-1}$) ', size=7)

#--- % CH_4 :--------------------------------------------
plt.subplot(3,2,6)
plt.tick_params(axis = 'both', labelsize = 7)
plt.plot(t,M[7,:],color='brown')
plt.xlabel('time (d)',size=7)
plt.title('% CH_4', size=7)

#--- Bacterial population :------------------------------

plt.subplot(3,2,3)
plt.tick_params(axis = 'both', labelsize = 7)
plt.plot(t,M[4,:],color='blue')
plt.title('Bacterial_population (kg COD.m$^{-3}$)', size=7)

plt.savefig('FIGURE2.png')
plt.close()'''

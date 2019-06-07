from  scipy  import*
from  scipy.integrate import odeint
import numpy as np
import  matplotlib.pyplot  as  plt

a=-2
t= linspace(0,10 ,30)

def equadiff_resolution (x,t) :
    return a*x

x= odeint (equadiff_resolution, [0.2] , t)
plt . plot (t , x)
plt . show ()

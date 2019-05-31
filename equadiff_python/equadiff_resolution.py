from  scipy  import*
from  scipy.integrate  import  odeint
import  matplotlib . pyplot  as  plt
def equadiff_resolution(x,t):
    return -2*x(t)
    t=linspace(0,1,12)
    x0 = 0.1
    dx=-2*x
    solution= odeint ( equadiff_resolution, x0 , t)
    x= solution[ : , 0 ]
    pypl.plot(t,dx)
    plt.show()
    return 0
    

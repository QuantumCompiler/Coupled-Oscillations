############### Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

############### Equations
def model(IC,t):
    x, dx, o1, do1, o2, do2 = IC
    mp=0.08
    ms=2.47
    mw=0.38
    l=0.357
    r=0.04
    g=9.8066
    I=0.0005472
    a=(2*((I/r**2+mp)/mp)+ms+0.5*mw)
    b=((do1)**2*np.sin(o1)*np.cos(o1))
    B=((do2)**2*np.sin(o2)*np.cos(o2))
    c=(g*np.sin(o1)/(l)*a)
    d=(g*np.sin(o2)/(l)*a)
    f=((do1)**2*np.sin(o1)*np.cos(o2))
    g=((do2)**2*np.sin(o2)*np.cos(o1))
    h=(a-np.cos(o1)**2)
    j=(a-np.cos(o2)**2)
    k=(np.cos(o1)*np.cos(o2))
    l=(np.cos(o1)**2*np.cos(o2)**2)
    d2dx=(l*((do1)**2*np.sin(o1)+(do2)**2*np.sin(o2))+
          0.5*g*(np.sin(2*o1)+np.sin(2*o2)))/(a+((np.cos(o1))**2)+(np.cos(o2))**2)
    d2do1=(b+B-c-k+(B+f-d)/j)/((1-(l/(h*j)))*h)
    d2do2=(B+b-d-k+(b+g-c)/h)/((1-(l/(j*h)))*j)
    f=[dx,d2dx,do1,d2do1,do2,d2do2]
    return f

############### Initial Conditions
x0=0.0
dx0=0.0
o10=-0.0434
do10=0.0
o20=0.0378
do20=0.0
ic=[x0,dx0,o10,do10,o20,do20]

############### Time Space
t=np.linspace(0,120,2400)

############### Solution to ODE
y=odeint(model,ic,t)

############### Plots
############### x(t)
'''plt.plot(t,y[:,0])
plt.title('X(t)')
plt.xlabel('Time in Seconds')
plt.ylabel('Position in Meters')
plt.show()'''

############### dx(t)
'''plt.plot(t,y[:,1])
plt.title('dX(t)')
plt.xlabel('Time in Seconds')
plt.ylabel('Velocity in m/s')
plt.show()'''

############### o1(t)
'''plt.plot(t,y[:,2])
plt.title('o1(t)')
plt.xlabel('Time in Seconds')
plt.ylabel('Angular Position in Radians')
plt.show()'''

'''############### do1(t)
plt.plot(t,y[:,3])
plt.title('do1(t)')
plt.xlabel('Time in Seconds')
plt.ylabel('Angular Velocity in Radians/Sec')
plt.show()'''

############### o2(t)
'''plt.plot(t,y[:,4])
plt.title('o2(t)')
plt.xlabel('Time in Seconds')
plt.ylabel('Angular Position in Radians')
plt.show()'''

############### do2(t)
'''plt.plot(t,y[:,5])
plt.title('do2(t)')
plt.xlabel('Time in Seconds')
plt.ylabel('Angular Velocity in Radians/Sec')
plt.show()'''

############### Dual
plt.plot(t,y[:,2],label='o1')
plt.plot(t,y[:,4],label='o2')
plt.title('')
plt.xlabel('Time in Seconds')
plt.ylabel('Angular Position in Radians')
plt.legend(loc='upper left')
plt.show()

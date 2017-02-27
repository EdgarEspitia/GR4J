import matplotlib.pyplot as plt
import numpy as np
import GR4J as hm
import criteria as cri
from scipy.optimize import minimize

data = np.genfromtxt("input_data.csv", delimiter = ",")
obs = data[:,3]

def fun_obj(x):

    '''
    Nash-Suttclife function
    Objetive function for search the parameter of model GR4J
    '''

    X =[0]
    X.extend(x)  
    my_basin = hm.GR4J(data, 260.0, X, 100.0, 100.0)
    sim =  my_basin.run()
    return np.sum((obs-sim)**2) / (np.sum((obs-obs.mean())**2))


# X0 is a estimation where the minimum can be located
x0 = [320, 2.42, 60, 1.3]
res = minimize(fun_obj, x0, method = 'nelder-mead',
               options = {'xtol': 1e-8, 'disp': True})

 
print (res.x)


X = [0]
X.extend(res.x)
my_basin = hm.GR4J(data, 260.0, X, 192.06, 48.73 )
Q_fore = my_basin.run()
print("Model GR4J")
print("parameters Fore. Flows caculate by  Nelder-Mead algoritm " )
print(X)
print( "NSE = {0:.2f} " .format( cri.nash(data[:, 3], Q_fore) ) )


X = [0, 320.10 , 2.420,  69.62, 1.38 ]
my_basin = hm.GR4J(data, 260.0, X, 192.06, 48.73 )
Q_fore_opt = my_basin.run()
print("parameters Fore_opt. Flows  " )
print(X)
print( "NSE = {0:.2f} " .format( cri.nash(data[:, 3], Q_fore_opt) ) )





fig, ax = plt.subplots()
line1, = ax.plot(data[:,0], data[:, 3]     , label='Obs. Flows ')
line2, = ax.plot(data[:, 0], Q_fore ,  label='Fore. Flows')
line3, = ax.plot(data[:, 0], Q_fore_opt ,  label='Fore_opt. Flows')
ax.legend(loc='upper right')
plt.show()

import numpy as np
import matplotlib.pyplot as plt
# importin the GR4J model
import GR4J as hm

# Read the data file from csv file with delimiter ","
data = np.genfromtxt("input_data.csv", delimiter = ",")

# Parameters value
#X = [0, 320.11, 2.42,  69.63, 1.39 ]
X = [0, 320.107309225224 , 2.42084329704772,  69.6276180318729, 1.38913446625445 ]

#   GR4J(PEQ, A, X, S0, R0):
# A the area of basin
# S0 and R0 inicial values 
my_basin = hm.GR4J(data, 260.0, X, 192.064385535134, 48.7393326223111 )
#for i in range(1,1):
Q_fore = my_basin.run()


for i in range(1, len(Q_fore)):
    print(Q_fore[i])


fig, ax = plt.subplots()
line1, = ax.plot(data[:,0], data[:, 3]     , label='Obs. Flows ')
line2, = ax.plot(data[:, 0], Q_fore ,  label='Fore.Flows')
ax.legend(loc='upper right')
plt.show()








    




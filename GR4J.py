import numpy as np

'''
This class modelate a proccess rainfall-runoff

Perrin, C., C. Michel, et al. (2003). "Improvement of a parsimonious model for streamflow simulation." Journal of Hydrology 279(1-4): 275-289 

http://www.cemagref.fr/webgr/Modelesgb/gr4j/fonctionnement_gr4jgb.htm 
'''

class GR4J:
    'Class hydrological model GR4J'

    def __init__(self, PEQ, A, X, S0, R0):
        self.set_PEQ(PEQ)
        self.a = np.copy(A) # Area of basin
        n = PEQ.shape[0] # Get the lenght of data series
        self.X =  np.copy(X)  # Parameters of model
        self.Pn = np.zeros((n))
        self.En = np.zeros((n))
        self.R_ = np.copy(R0)
        self.S_ = np.copy(S0)
       
    # Define the imputs 
    def set_PEQ(self, PEQ):
        self.P = np.copy(PEQ[:,1])
        self.E = np.copy(PEQ[:,2])
        self.Q = np.copy(PEQ[:,3])

        

    # Definition of precipitation P in mm/day
    def set_P(self, P):
        self.P = np.copy(P)

    
    # Definito of E the potential evapotranspiration (PE) in mm/day
    def set_E(self, E):
        self.E = np.copy(E)

    # Definition of flow in m^3/s 
    def set_Q(self, Q):
        self.Q = np.copy(Q)

    def set_S_0(self, S_0):
        self.S_ = np.copy(S_0)

    def set_R_0(self, R_0):
        self.R_ = np.copy(R_0)

    #Definition the area of basin
    def set_a(self,a):
        self.a = np.copy(a)

    def set_parameters(self, X):
        self.X = np.copy(X)

        
        
    # Print status of variables    

    def gen_unity_hydrogram(self):

        # n can be arbitrary but in order to reduce computing time is better chose an size the parameter X[4] + 2, it is because the hydragram univtary only get values for geater integer value of X[] 
        n = int(np.ceil(self.X[4]))+ 2
        #n = 10
        self.HU1 = np.zeros(n)
        self.SH1 = np.zeros(n)
        self.HU1_09 = np.zeros(n)
        m = 2 * n

        self.HU2 = np.zeros(m)
        self.SH2 = np.zeros(m)
        self.HU2_01 = np.zeros(m)
        
        for i in range(1,n,1):
            if i < self.X[4] : 
                self.SH1[i] = np.power( float(i) /self.X[4]  , 5.0/2.0)
            else:
                self.SH1[i] = 1.0

        for i in range(1,m,1):
            if i < self.X[4] : 
                self.SH2[i] = 1.0/2.0 * np.power( float(i) /self.X[4]  , 5.0/2.0)
            elif i > self.X[4] and i < 2.0*self.X[4]  :
                self.SH2[i] = 1.0 - 1.0/2.0 * np.power( 2 - float(i)/self.X[4] , 5.0/2.0)
            elif i > 2.0*self.X[4]:
                self.SH2[i] = 1.0
    
        for i in range(1,n-1,1):
            self.HU1[i] = self.SH1[i] - self.SH1[i-1]

        for i in range(1,m-1,1):
            self.HU2[i] = self.SH2[i] - self.SH2[i-1]
    
        
        
    def run(self):

        # Creating hydrograms
        self.gen_unity_hydrogram()
        
     
        self.Es_ = 0.0    
        self.F_ = 0.0
        self.Pr_ = 0.0
        self.Perc_ = 0.0
        self.Qr_  = 0.0
        self.Qd_ = 0.0
        self.Ps_ = 0
        
        
        
        # N Lenght of data series
        N = self.P.size 
        n = len(self.HU1) # n Lenght of hydrogram 1
        m = len(self.HU2) # n Lenght of hydrogram 2
        
        #Determination of net rainfall and PE
        for i in range(0, N ):
            if self.P[i] >= self.E[i]:
                # Equation 1
                self.Pn[i] = self.P[i] - self.E[i]
                self.En[i] = 0
            else:
                #Equation 2
                self.Pn[i] = 0
                self.En[i] = self.E[i] - self.P[i]


            #Production (SMA) store
            if self.Pn[i] > 0:
                #Equation 3
                self.Ps_ = ( self.X[1]*(1.0 - (self.S_  / self.X[1])**2.0  ) * np.tanh( self.Pn[i]/self.X[1] ) )       / ( 1 + self.S_  /self.X[1] * np.tanh(self.Pn[i]/self.X[1] ))
            else:
                self.Ps_ = 0

            if self.En[i] > 0:
                #Equation 4
                self.Es_ = ((2.0 - self.S_ /self.X[1])*self.S_ *np.tanh(self.En[i]/self.X[1]) ) / (1 + (1-self.S_ /self.X[1])*np.tanh(self.En[i]/self.X[1]))
            else:
                self.Es_ = 0


            #Equation 5

            self.S_ = self.S_ - self.Es_ + self.Ps_

            #Equation 6
            self.Perc_ = self.S_ * (  1.0 - np.power( 1.0 + np.power(4*self.S_ /(9.0*self.X[1] ) ,4.0)  , -1.0/4.0)  )

            #Equation 7
            self.S_ = self.S_ - self.Perc_

            #Equation 8
            self.Pr_ = self.Perc_ + (self.Pn[i] - self.Ps_)

            #Equation 9 to 17 in  gen_unity_hydrogram
            if i > 0:
                for j in range(1, n-1):
                    if j == n-1:
                        self.HU1_09[j] = self.Pr_ *  self.HU1[j] * 0.9
                    else:
                        self.HU1_09[j] = self.HU1_09[j+1] + self.Pr_ *  self.HU1[j] * 0.9
                
                for j in range(1, m -1):
                    if j == m-1:
                        self.HU2_01[j] = self.Pr_ *  self.HU2[j] * 0.1
                    else:
                        self.HU2_01[j] = self.HU2_01[j+1] + self.Pr_ *  self.HU2[j] * 0.1

            else:
                self.HU1_09 = self.Pr_ *  self.HU1 * 0.9
                self.HU2_01 = self.Pr_ *  self.HU2 * 0.1

            #Equation 18
            self.F_ =  self.X[2] * np.power( self.R_ /self.X[3]  , 7.0/2.0)
            #Equation 19
            self.R_ = max( 0.0, self.R_ + self.HU1_09[1] + self.F_ )
            
            #Equation 20
            self.Qr_ = self.R_ * ( 1.0 - np.power(1 + np.power( self.R_ /self.X[3] ,4.0)  , - 1.0/4.0 )   )

            #Equation 21
            self.R_ = self.R_ - self.Qr_

            #Equation 22
            self.Qd_ = max(0,self.HU2_01[1] + self.F_)
                
            #Equatiion 23
            self.Q[i] = self.Qr_ + self.Qd_

        
        #Retrun the unitary flow
        #return self.Q
#        for i in range(0, len(self.Q)):
 #           print(self.Q[i])
        # Return the flow as m3/s
        return self.Q * self.a * 1000.0/86400.0 

    



    

    

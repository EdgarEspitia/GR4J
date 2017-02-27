import numpy as np

def nash(obs, sim):
    '''
    Nash suttclife criteria
    '''
    if ( len( obs ) != len( sim  )):
        print "ERROR INPUT DATA!"
        return -1

    return 1 - np.sqrt( np.sum((obs-sim)**2) / (np.sum((obs-obs.mean())**2)))


def RMSE( obs , sim ):
    '''
    Root mean square error
    '''
    if ( len( obs ) != len( sim  )):
        print "ERROR INPUT DATA!"
        return -1           
    return np.sqrt( sum((obs-sim )**2  )  / len(sim) )



import numpy as np
import matplotlib.pyplot as plt

def moving_average(x, w):
    """ 
    x: input signal
    w: window
    """
    return np.convolve(x, np.ones(w), 'valid') / w

plt.show()

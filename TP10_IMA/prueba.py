import numpy as np
import matplotlib.pyplot as plt
from scipy import stats 

asd = np.load('asd.npy')
a = stats.mode(asd)

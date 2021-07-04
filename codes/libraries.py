import numpy as np
import time
import sys
import scipy
import pickle
from scipy import interpolate
from scipy.integrate import odeint
from scipy.interpolate import RegularGridInterpolator
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rc
from astropy import units as u
from astropy import constants as const
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

import dsharp_opac as opacity

import warnings
warnings.filterwarnings("ignore")


#sudo apt-get install texlive-latex-extra texlive-fonts-recommended dvipng   #para intentar labels en LaTex:    usetex=True



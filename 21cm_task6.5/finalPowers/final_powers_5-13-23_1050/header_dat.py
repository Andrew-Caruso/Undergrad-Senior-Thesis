import pandas as pd
#import cosmolopy.distance as cd
from astropy.cosmology import WMAP9 as cosmo

#import nbodykit as nbdk
#from nbodykit.lab import ArrayMesh
#from nbodykit.source.catalog import ArrayCatalog

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import pickle

import sys, glob
import multiprocessing as mp
import seaborn as sns
from colossus.cosmology import cosmology
from colossus.lss import mass_function
import scipy.ndimage as ndi


SNAPNUM_FLAG = False
NPROCS = 96


#Cosmology
omegam = 0.305147
omegab = 0.0482266
sig8 = 0.82033
hh = 0.68
ns = 0.9667
fbaryon = omegab/omegam
rho_c_cgs = 1.87890e-29

c = 3E10
MPC2KM = 3.086E19
KM2CM = 1E5
MPC2CM = KM2CM*MPC2KM


MPROTON = 1.6726e-24 # in grams
KBOLTZ = 1.38066e-16 # in erg/K
MYRTOSEC = 3.15569e13
KPCTOCM = 3.08560e21 # cm/kpc
PCTOCM = 3.08568025e18 # cm per parsec 
MPCTOCM = 3.08560e24 # cm/Mpc
TCMB = 2.726 #CMB temperature at z=0
MUB = 1.22 #Mean molecular weight, for primordial
LIGHTSPEED = 3.0e5 #km s^-1
LIGHTSPEED_CGS = 3e10 #cm/s 
mH_CGS = 1.67223e-24 # Mass of Hydrogen in grams
me_CGS = 9.10939e-28 # Mass of electron (g) 
e_c = 4.803e-10 #electron charge (esu) 
mHe_CGS = 4*mH_CGS 
H0_CGS   = 3.24086e-18
Tcmb0 = 2.7255 # //Present day CMB temp in K
GAMMA_IDEAL = 5.0/3.0 

OMEGAM = 0.305147
OMEGAB = 0.0482266
HUBBLE_h = 0.68
YHe = 0.2453
XHy = 1.0 - YHe
FBARYON = OMEGAB/OMEGAM



EION_HEII = 54.4 #ionization potential of HeII in eV
SIGMA_HEII = 6.3e-18/4.*1.21 #cm^2 -- photoionization cross section for HeII
SIGMA_HI = 6.3463e-18# cm^2 -- photoionization cross section for HI
EVTOERG = 1.60217646e-12
H0_CGS   = 3.24086e-18 #in 1/s.  Does not include h.
K2eV = KBOLTZ/EVTOERG
eV2K = 1.0/K2eV
HI_eV = 13.6
HI_K = HI_eV*eV2K
Planck_h = 6.626075540e-27 #Planck's constant in erg/s



GAL_ALPHA = 2.0
GAMMA_RH = 1.0E-13

#Initialize cosmology
params = {'flat': True, 'H0': hh*100., 'Om0': omegam, 'Ob0': omegab, 'sigma8': sig8, 'ns': ns}
cosmology.addCosmology('myCosmo', params)
cosmo = cosmology.setCosmology('myCosmo')

#Physical Constants
rhoc_ast = 2.743e11 # 2.77550e11 #present day crit density in Msun/h / (Mpc/h)^3
rho0 = omegam*rhoc_ast # in Msun/h / (Mpc/h)^3
rho0 = rho0*(hh**3/hh)*(1.989E33)/(3.086E24**3) # comoving cgs units


cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'omega_k_0': 0.0, 'h' : 0.68}


# Hydro Simulation Properties 
HYDRO_DIR = "/expanse/lustre/scratch/ccain002/temp_project/for_garett/N1024_run/"
UGAS_BASE = "u_gas.z="
VGAS_BASE = "v_gas.z="
UGAS_NELEM = 5
VGAS_NELEM = 3
LBOX_HYDRO = 200.
NGRID_HYDRO = 1024
DELTAX_HYDRO = float(LBOX_HYDRO)/float(NGRID_HYDRO) #Mpc/h
Z_HYDRO = '07.0000'

# Box Properties 
DATAFIELDS_DIRECTORY = "" #"/expanse/lustre/scratch/ccain002/temp_project/naidu_project/" # "/expanse/lustre/scratch/ccain002/temp_project/for_anson/caseA_models/FINAL/"
LBOX = 200.
NGRID = 200
DELTAX = LBOX/float(NGRID)
FLEXRT_HYDRO_PROP = float(NGRID)/float(NGRID_HYDRO)

# Sightlines Properties
SKEW_LENGTH = 100.
SKEW_NGRID = 100
SKEW_RES = DELTAX #SKEW_LENGTH/float(SKEW_NGRID)
SKEW_HYDRO_NGRID = int(SKEW_NGRID/FLEXRT_HYDRO_PROP+0.5)
SKEW_HYDRO_DELTAX = float(SKEW_LENGTH)/float(SKEW_HYDRO_NGRID)
NSIGHTLINES = 100


# LAE Catalog Properites
CATALOG_DIRECTORY = "/expanse/lustre/scratch/ngang002/temp_project/LAE_Bubble_Project_Scratch/LAE_Catalog_M1e9_nodc/"
NDIM = 23


# Galaxy Overdensity Properties
DELTA_DIRECTORY = "/expanse/lustre/scratch/ngang002/temp_project/LAE_Bubble_Project_Scratch/Bright_Galaxy_Sightlines"

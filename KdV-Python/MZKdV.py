# Author:  Jake Price
# Date:    February 8, 2021
# Purpose: Collection of custom functions for running CMA simulations of KdV,
#          analysis functions, and renormalization scripts. Generates some images
#          from past papers, but not all (didn't want / need to duplicate everything)
#          Translation of code from UW PhD in Matlab.

# import libraries
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import animation
import glob
import re



def fftnorm(u_full):
    """Computes normalized FFT (such that FFT and IFFT are symmetrically normalized)
    
    Parameters
    ----------
    u_full : 1D Numpy Array (N,)
        The vector whose discrete FFT is to be computed

    Returns
    -------
    normalizedFFT : 1D Numpy Array (N,)
        The transformed version of that vector
    """
    
    N = u_full.shape[0]
    normalizedFFT = np.fft.fft(u_full)*1/N
    return normalizedFFT

def ifftnorm(u_full):
    """Computes normalized IFFT (such that FFT and IFFT are symmetrically normalized)
    
    Parameters
    ----------
    u_full : 1D Numpy Array (N,)
        The vector whose discrete IFFT is to be computed

    Returns
    -------
    normalizedIFFT : 1D Numpy Array (N,)
        The transformed version of that vector
    """
    
    N = u_full.shape[0]
    normalizedIFFT = np.real(np.fft.ifft(u_full)*N)
    return normalizedIFFT

def convolutionSumKdV(u,v,alpha):
    """Computes convolution sum associated with RHS of KdV ODE
    
    C_k(u,v) = -(alpha * 1i * k) / 2 * sum_{i+j = k} u_i v_j
    
    Computed in real space to avoid loops and then converted back
    to Fourier space.
    
    Parameters
    ----------
    u : 1D Numpy Array (N,)
        One of the two vectors being convolved
        
    v : 1D Numpy Array (N,)
        One of the two vectors being convolved
        
    alpha : float
        Degree of nonlinearity in KdV

    Returns
    -------
    convo : 1D Numpy Array (N,)
        Convolution of the two vectors
    """
    
    # generate array of wavenumbers
    L = u.shape[0]
    k = np.concatenate([np.arange(0,L/2),np.arange(-L/2,0)])
    if v.shape[0]!=L:
        raise NameError('u and v must be the same length.')
    
    # compute double sum in real space, then apply scalar multiplier
    convo = fftnorm(ifftnorm(u)*ifftnorm(v))
    convo = -alpha/2*1j*k*convo
    return convo


# RHS: Right hand side functions for CMA and non-renormalized KdV

def markovKdV(u,M,alpha):
    """Computes nonlinear part of Markov term in KdV
    
    C_k(u,v) = -(alpha * 1i * k) / 2 * sum_{i+j = k} u_i v_j
    
    where the sum of i and j is over a "full" system with M positive modes (user specified)
    
    Computed in real space to avoid loops and then converted back
    to Fourier space.
    
    Parameters
    ----------
    u : 1D Numpy Array (N,)
        Positive modes of state vector whose RHS is being computed
        
    M : int
        Number of positive modes in "full" model for intermediary calculations
        
    alpha : float
        Degree of nonlinearity in KdV

    Returns
    -------
    nonlin0 : 1D Numpy Array (2*M,)
        Nonlinear part of Markov term for given state vector
        
    u_full : 1D Numpy array (2*M,)
        "full" state vector for use in later computations
    """
    
    # construct full Fourier vector from only the positive modes
    u_full = np.zeros(2*M) +1j*np.zeros(2*M)
    u_full[0:u.shape[0]] = u
    u_full[2*M-u.shape[0]+1:] = np.conj(np.flip(u[1:]))
    
    # compute the convolution sum
    nonlin0 = convolutionSumKdV(u_full,u_full,alpha)
    return nonlin0,u_full

def tModelKdV(u_full,nonlin0,alpha,F_modes):
    """Computes t-model term in KdV
    
    C_k(u,v) = -(alpha * 1i * k) / 2 * sum_{i+j = k, i and j in F} u_i v_j
    
    where the sum of i and j is over a "full" system with M positive modes (user specified)
    
    Computed in real space to avoid loops and then converted back
    to Fourier space.
    
    Parameters
    ----------
    u_full : Numpy array (2M,1)
             Current state of u in full form
            
    nonlin0 : Numpy array (2M,1)
              Markov term (for convolving)
             
    alpha : float
        Degree of nonlinearity in KdV
        
    F_modes : Numpy array
              Set of resolved modes (and aliasing modes) to zero out
              
    Returns
    -------
    nonlin1 : 1D Numpy Array (2*M,)
        t-model term
        
    uuStar : 1D Numpy array (2*M,)
        unresolved modes of state vector convolved with itself
    """
    
    uuStar = np.copy(nonlin0)
    uuStar[F_modes] = 0
    
    nonlin1 = 2*convolutionSumKdV(u_full, uuStar, alpha)
    
    return nonlin1,uuStar

def t2ModelKdV(u_full,nonlin0,uuStar,alpha,F_modes,G_modes,k,epsilon):
    """Computes second order ROM term in KdV
    
    *see paper / symbolic notebook for expression*
    
    Computed in real space to avoid loops and then converted back
    to Fourier space.
    
    Parameters
    ----------
    u_full : Numpy array (2M,1)
             Current state of u in full form
            
    nonlin0 : Numpy array (2M,1)
              Markov term (for convolving)
              
    uuStar : 1D Numpy array (2*M,)
        Unresolved modes of state vector convolved with itself
             
    alpha : float
        Degree of nonlinearity in KdV
        
    F_modes : Numpy array
              Set of resolved modes (and aliasing modes) to zero out
              
    G_modes : Numpy array
              Set of unresolved modes (and aliasing modes) to zero out
              
    k  :  Numpy array (2M,1)
          Array of wavenumbers
          
    epsilon : float
        Size of linear term (stiffness)
        
    Returns
    -------
    nonlin2 : 1D Numpy Array (2*M,)
        t2-model term
        
    uk3 : 1D Numpy array (2*M,)
        Resolved modes of state vector multiplied by k^3
        
    uu : 1D Numpy array (2*M,)
        Resolved modes of state vector convolved with itself
        
    A, AStar, B, BStar, C, CStar, D, DStar : 1D Numpy arrays (2*M,)
        Specific convolutions used as inner terms in future terms
    """
    
    # compute inner convolutions
    uu = np.copy(nonlin0)
    uu[G_modes] = 0
    
    uk3 = k**3*u_full
    
    A = k**3*uu
    AStar = k**3*uuStar
    
    B = convolutionSumKdV(1j*epsilon**2*uk3+uu,u_full,alpha)
    BStar = np.copy(B)
    B[G_modes] = 0
    BStar[F_modes] = 0
    
    C = convolutionSumKdV(uuStar,u_full,alpha)
    CStar = np.copy(C)
    C[G_modes] = 0
    CStar[F_modes] = 0
    
    D = convolutionSumKdV(uuStar,uuStar,alpha)
    DStar = np.copy(D)
    D[G_modes] = 0
    DStar[F_modes] = 0
    
    # compute actual term
    nonlin2 = -2*convolutionSumKdV(u_full,1j*epsilon**2*AStar - 2*BStar + 2*CStar,alpha) - 2*D
    
    return nonlin2,uk3,uu,A,AStar,B,BStar,C,CStar,D,DStar

def t3ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,A,AStar,B,BStar,C,CStar,DStar):
    """Computes third order ROM term in KdV
    
    *see paper / symbolic notebook for expression*
    
    Computed in real space to avoid loops and then converted back
    to Fourier space.
    
    Parameters
    ----------
    alpha : float
        Degree of nonlinearity in KdV
        
    F_modes : Numpy array
              Set of resolved modes (and aliasing modes) to zero out
              
    G_modes : Numpy array
              Set of unresolved modes (and aliasing modes) to zero out
              
    k  :  Numpy array (2M,1)
          Array of wavenumbers
          
    epsilon : float
        Size of linear term (stiffness)
    
    u_full : Numpy array (2M,1)
             Current state of u in full form
             
    uu : 1D Numpy array (2*M,)
        Resolved modes of state vector convolved with itself
            
    uuStar : 1D Numpy array (2*M,)
        Unresolved modes of state vector convolved with itself
        
    uk3 : 1D Numpy array (2*M,)
        Resolved modes of state vector multiplied by k^3
             
    A, AStar, B, BStar, C, CStar, DStar : 1D Numpy arrays (2*M,)
        Specific convolutions used as inner terms in future terms

        
    Returns
    -------
    nonlin3 : 1D Numpy Array (2*M,)
        t3-model term
        
    uk6 : 1D Numpy array (2*M,)
        Resolved modes of state vector multiplied by k^6
        nonlin3,uk6,E,EStar,F,FStar
        
    E, EStar, F, FStar : 1D Numpy arrays (2*M,)
        Specific convolutions used as inner terms in future terms
    """
    
    # compute internal convolutions
    uk6 = k**3*uk3
    
    E = convolutionSumKdV(1j*epsilon**2*uk3+uu,1j*epsilon**2*uk3+uu,alpha)
    EStar = np.copy(E)
    E[G_modes] = 0
    EStar[F_modes] = 0
    
    F = convolutionSumKdV(uuStar,1j*epsilon**2*uk3+uu,alpha)
    FStar = np.copy(F)
    F[G_modes] = 0
    FStar[F_modes] = 0
    
    int1 = -2*BStar+CStar
    int2 = (convolutionSumKdV(u_full,
                             -epsilon**4*uk6
                             +1j*epsilon**2*(A+AStar)
                             +2*(B-2*C)
                             +2*(CStar-2*BStar),
                             alpha))
    int2[F_modes] = 0
    int3 = EStar-FStar
    int4 = np.copy(DStar)
    int5 = CStar-BStar
    
    # compute actual 3rd order term
    nonlin3 = (2*convolutionSumKdV(u_full,-k**3*epsilon**4*AStar
                                  +2*1j*epsilon**2*k**3*int1
                                  +2*int2+2*int3+2*int4,alpha) 
              +6*convolutionSumKdV(uuStar,1j*epsilon**2*AStar + 2*int5,alpha))
    
    return nonlin3,uk6,E,EStar,F,FStar

def t4ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,uk6,A,AStar,B,BStar,C,CStar,D,DStar,E,EStar,F,FStar):
    """Computes fourth order ROM term in KdV
    
    *see paper / symbolic notebook for expression*
    
    Computed in real space to avoid loops and then converted back
    to Fourier space.
    
    Parameters
    ----------
    alpha : float
        Degree of nonlinearity in KdV
        
    F_modes : Numpy array
              Set of resolved modes (and aliasing modes) to zero out
              
    G_modes : Numpy array
              Set of unresolved modes (and aliasing modes) to zero out
              
    k  :  Numpy array (2M,1)
          Array of wavenumbers
          
    epsilon : float
        Size of linear term (stiffness)
    
    u_full : Numpy array (2M,1)
             Current state of u in full form
             
    uu : 1D Numpy array (2*M,)
        Resolved modes of state vector convolved with itself
            
    uuStar : 1D Numpy array (2*M,)
        Unresolved modes of state vector convolved with itself
        
    uk3 : 1D Numpy array (2*M,)
        Resolved modes of state vector multiplied by k^3
        
    uk6 : 1D Numpy array (2*M,)
        Resolved modes of state vector multiplied by k^6
             
    A, AStar, B, BStar, C, CStar, DStar, E, EStar, F, FStar : 1D Numpy arrays (2*M,)
        Specific convolutions used as inner terms in future terms

        
    Returns
    -------
    nonlin4 : 1D Numpy Array (2*M,)
        t4-model term
    """
    
    # compute internal convolutions
    internal1 = (convolutionSumKdV(u_full,-epsilon**4*uk6+1j*epsilon**2*(A+AStar)
                                   +2*B-4*C-4*BStar+2*CStar,alpha))
    
    internal1[F_modes] = 0
    
    internal2 = (1j*epsilon**2*k**3*convolutionSumKdV(u_full,-3*epsilon**4*uk6
                                                      +1j*epsilon**2*(3*A+AStar)
                                                      -2*(-3*B+5*C)
                                                      +2*(-3*BStar+CStar),alpha))
    internal2[F_modes] = 0
    
    auxiliary1 = 2*convolutionSumKdV(u_full,epsilon**4*uk6-1j*epsilon**2*(A+3*AStar)
                                     +2*(3*C-B)+2*(5*BStar-3*CStar),alpha)
    auxiliary1[G_modes] = 0
    
    auxiliary2 = 2*convolutionSumKdV(u_full,-3*epsilon**4*uk6+1j*epsilon**2*(3*A+AStar)
                                     +2*(3*B-5*C)+2*(-3*BStar+CStar),alpha)
    auxiliary2[F_modes] = 0
    
    internal3 = convolutionSumKdV(u_full,1j*k**3*uk6*epsilon**6
                                  +k**3*epsilon**4*(A-AStar)
                                  +2*1j*epsilon**2*k**3*(3*C-B)
                                  +2*1j*epsilon**2*k**3*(-3*BStar+CStar)
                                  +auxiliary1+auxiliary2
                                  -2*(E-2*F)
                                  +2*(3*EStar-2*FStar)
                                  -6*D+2*DStar,alpha)
    internal3[F_modes]= 0
    
    internal4 = convolutionSumKdV(1j*epsilon**2*uk3+uu,3*epsilon**4*uk6-1j*epsilon**2*(3*A+AStar)
                                  +2*(-3*B+5*C)+2*(3*BStar-CStar),alpha)
    internal4[F_modes] = 0
    
    internal5 = convolutionSumKdV(uuStar,-epsilon**4*uk6+1j*epsilon**2*(A+3*AStar)
                                  +2*B-6*C-10*BStar+6*CStar,alpha)
    internal5[F_modes] = 0
    
    # compute actual fourth order term
    nonlin4 = (2*convolutionSumKdV(u_full,-1j*epsilon**6*k**6*AStar
                                  +2*k**6*epsilon**4*(3*BStar-CStar)
                                  +2*internal2
                                  +2*internal3
                                  +2*internal4
                                  -2*k**3*1j*epsilon**2*(2*FStar-3*EStar)
                                  +2*k**3*1j*epsilon**2*DStar
                                  +2*internal5,alpha)
               +8*convolutionSumKdV(uuStar,-k**3*epsilon**4*AStar
                                    +2*1j*epsilon**2*k**3*(-2*BStar+CStar)
                                    +2*internal1
                                    +2*(EStar-FStar)
                                    +2*DStar,alpha)
               -48*convolutionSumKdV(BStar,1j*epsilon**2*AStar+2*CStar,alpha)
               +6*convolutionSumKdV(1j*epsilon**2*AStar+2*(BStar+CStar),
                                    1j*epsilon**2*AStar+2*(BStar+CStar),alpha)
               )
    
    nonlin4 = -nonlin4
    return nonlin4



def RHSKdV(t,u,params):
    """
    Computes the RHS for a full KdV or ROM simulation. For use in solver.
    
    Parameters
    ----------
    t : float
        Current time
        
    u : Numpy array (N,)
        Current state vector
              
    params : Dictionary
             Dictionary of relevant parameters (see below)
        N : float, number of positive modes in simulation
        M : float, number of positive modes in "full" intermediate compuation
        alpha : float, degree of nonlinearity in KdV
        epsilon : float, size of linear term (stiffness)
        tau : float, time decay modifier
        coeffs : Numpy array, renormalization coefficients for ROM (None if no ROM)

        
    Returns
    -------
    RHS : 1D Numpy array (N,)
          Derivative of each positive mode in state vector
    """
    
    # extract parameters from dictionary
    N = params['N']
    M = params['M']
    alpha = params['alpha']
    epsilon = params['epsilon']
    tau = params['tau']
    coeffs = params['coeffs']
    
    # construct wavenumber array
    k = np.concatenate([np.arange(0,M),np.arange(-M,0)])
    
    
    # Linear and Markov term
    nonlin0,u_full = markovKdV(u,M,alpha)
    RHS = 1j*k[0:N]**3*epsilon**2*u + nonlin0[0:N]
    
    if (np.any(coeffs == None)):
        order = 0
    else:
        order = coeffs.shape[0]
    
    if (order >= 1):
        # compute t-model term
        
        # define which modes are resolved / unresolved in full array
        F_modes = np.concatenate([np.arange(0,N),np.arange(2*N-1,M+N+2),np.arange(2*M-N+1,2*M)])
        G_modes = np.arange(N,2*M-N+1)
    
        # compute t-model term
        nonlin1,uuStar = tModelKdV(u_full,nonlin0,alpha,F_modes)
        RHS = RHS + coeffs[0]*nonlin1[0:N]*t**(1-tau)
        
        order = coeffs.shape[0]
    
    if (order >= 2):
        # compute t2-model term
        nonlin2,uk3,uu,A,AStar,B,BStar,C,CStar,D,DStar = t2ModelKdV(u_full,nonlin0,uuStar,alpha,F_modes,G_modes,k,epsilon)
        RHS = RHS + coeffs[1]*nonlin2[0:N]*t**(2*(1-tau))
    
    if (order >= 3):
        # compute t3-model term
        nonlin3,uk6,E,EStar,F,FStar = t3ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,A,AStar,B,BStar,C,CStar,DStar)
        RHS = RHS + coeffs[2]*nonlin3[0:N]*t**(3*(1-tau))
    
    if (order == 4):
        # compute t4-model term
        nonlin4 = t4ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,uk6,A,AStar,B,BStar,C,CStar,D,DStar,E,EStar,F,FStar)
        RHS = RHS + coeffs[3]*nonlin4[0:N]*t**(4*(1-tau))

    return RHS


def getMass(u,N):
    """Computes mass in first N modes for all timesteps from solution array u
    
    Parameters
    ----------
    u : 2D Numpy Array (M,tList)
        Positive modes of state vector for all timesteps
        
    N : int
        Number of positive modes to include in mass measurement
        
    Returns
    -------
    mass : 1D Numpy Array (tList,)
        Energy in first N modes at all timesteps
    """

    mass = np.sum(2*(abs(u[0:N,]))**2,0)
    return mass




def runSim(params):
    """
    Runs an actual ROM or non-ROM simulation of KdV
    
    Parameters
    ----------
    params : Dictionary
             Dictionary of relevant parameters (see below)
        N : float, number of positive modes in simulation
        M : float, number of positive modes in "full" intermediate compuation
        alpha : float, degree of nonlinearity in KdV
        epsilon : float, size of linear term (stiffness)
        tau : float, time decay modifier
        coeffs : Numpy array, renormalization coefficients for ROM (None if no ROM)
        IC : function handle, initial condition of simulation
        endtime : float, final time to simulate to
        timesteps: Numpy array, specific timesteps for which to save solution

        
    Returns
    -------
    uSim : ODE solver output
           Output solution from sp.integrate.solve_ivp (includes state vector at all timesteps, time vector, etc.)
    """
    
    # unpack parameters from dictionary
    N = params['N']
    IC = params['IC']
    endtime = params['endtime']
    timesteps = params['timesteps']
    
    # generate initial condition
    x = np.linspace(0,2*np.pi-2*np.pi/(2*N),2*N)
    y = IC(x)
    uFull = fftnorm(y)
    u = uFull[0:N]
    
    # define RHS in form appropriate for solve_ivp
    def myRHS(t,y):
        out = RHSKdV(t,y,params)
        return out
    
    # solve the IVP
    uSim = sp.integrate.solve_ivp(fun = myRHS, t_span = [0,endtime], y0 = u,method = "BDF", t_eval = timesteps)
    return uSim


def makeRealSpace(u,N):
    """Takes a completed simulation and finds the real space solution at all timesteps for a chosen subset of modes
    
    Parameters
    ----------
    u : Numpy array (M,t)
        Output of simulation giving energy in first M positive modes for all timesteps t
    
    N : int
        Number of positive modes to use in real space
        
        
    Returns
    -------
    x     : Numpy vector (2xN,1)
            x-grid for plotting purposes
        
    uReal : Numpy array (2xN,t)
            Real space solution at all times
    """
    
    # identify shapes of arrays
    uShape = u.shape
    numTimes = uShape[1]
    
    # drop modes we don't wish to keep
    uNew = u[0:N,:]
    
    # generate full vector (with negative modes)
    uFull = np.zeros((2*N,numTimes)) + 1j*0
    uFull[0:N,:] = uNew
    uFull[2*N-N+1:,:] = np.conj(np.flip(uNew[1:,:],0))
    
    # initialize output
    uReal = np.zeros(uFull.shape)
    
    # take inverse transform for each timestep
    # NOTE: is there a vectorized way to do this?
    for i in np.arange(0,numTimes):
        uReal[:,i] = ifftnorm(uFull[:,i])
    
    return uReal
    
    



def makeAnimations(uList,t,legendList):
    """
    Creates an animation from a list of simulations
    
    Parameters
    ----------
    uList : List of Numpy arrays of size (N,T)
            Set of state vector evolutions to animate
            
    t : Numpy array (T,)
        Timesteps associated with simulations (must all be the same)
        
    legendList : List of strings
                 Labels for each simulation
        
    Returns
    -------
    anim : animation object
           output from animation.FuncAnimation
    """
    # identify the resolution to use for plots and generate x grid
    N = min([x.shape[0] for x in uList])
    xgrid = np.linspace(0,2*np.pi*(2*N-1)/(2*N),2*N)
    
    # generate real space solutions
    realSols = [makeRealSpace(x,N) for x in uList]
    
    # initialize figure
    myFig = plt.figure()
    ax = plt.subplot()
    ax.axis(xmin = 0,xmax = 2*np.pi-np.pi/N,ymin = -2, ymax = 4)
    
    # create empty list of lines to populate each iteration
    lineList = [ax.plot([],[]) for i in range(len(uList))]

    # define function to draw each frame
    def makeFrame(n):
        for i in range(len(uList)):
            lineList[i][0].set_data(xgrid,realSols[i][:,n])
        plt.title('t = '+str(round(t[n],1)))
        plt.legend(legendList, loc = "upper right")
        return lineList

    # generate animation
    anim = animation.FuncAnimation(fig = myFig,func = makeFrame,frames = t.shape[0])
    return anim

def renormalize(fullM, endtime, Nlist, Mlist, epsilon, alpha, tau, timesteps, IC = np.sin, plots = False):
    """
    Finds renormalization coefficients based on a single simulation. If the
    simulation doesn't yet exist, it creates it
    
    Parameters
    ----------
    fullM : int
            Size of full simulation to base fits on
            
    endtime : int
        Endtime of full simulation
        
    Nlist : list of ints
            List of resolutions for which to find coefficients
            
    Mlist : list of ints
            List of intermediary "full" simulations to use for ROMs
            
    epsilon : float
        size of linear term (stiffness)
        
    alpha : float
        degree of nonlinearity in KdV
        
    tau : float
        time decay modifier
        
    timesteps : Numpy array
        specific timesteps for which to save solution
        
    IC : function handle
        initial condition of simulation (default np.sin)
        
    plots : boolean
        Indicates whether to generate plots (default: False)
        
    Returns
    -------
    
    coeeffsArray1 : Numpy array (length(Nlist),1)
        Renormalization coefficients for t-model only
        
    coeffsArray2 : Numpy array (length(Nlist),2)
        Renormalization coefficients for t-model and t2-model only
        
    coeffsArray3 : Numpy array (length(Nlist),3)
        Renormalization coefficients for t1-t3-models
        
    coeffsArray4 : Numpy array (length(Nlist),4)
        Renormalization coefficients for t1-t4-models
        
    coeffsArray2only : Numpy array (length(Nlist),1)
        Renormalization coefficients for t2-model only
        
    coeffsArray24only : Numpy array (length(Nlist),2)
        Renormalization coefficients for t2-model and t4-model only
        
    fitLines : Dict
        Contains scaling law fits for each ROM coefficients
        of form   c = -b * N^a
        Terms given are a, b, and r (correlation coefficient of fit)
        
    err : Dict
        Contains least-squares error for each fit for each model and resolution
    """
    
    # Check if full simulation has already been constructed
    #   if so, load it, if not, generate it
    try:
        uFull = np.load("u" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+".npy")
        tFull = np.load("t" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+".npy")
    except:
        fullParams = {
            'N': fullM,
            'M': int(3/2*fullM),
            'alpha': 1,
            'epsilon': epsilon,
            'tau': 1,
            'coeffs': None,
            'IC': IC,
            'endtime': endtime,
            'timesteps': timesteps
            }

        uSimFull = runSim(fullParams)
        uFull = uSimFull.y
        tFull = uSimFull.t
        np.save( "u" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p'),uFull)
        np.save( "t" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p'),tFull)
     
    # initialize output arrays
    coeffsArray1 = np.zeros((Nlist.shape[0],1))
    coeffsArray2 = np.zeros((Nlist.shape[0],2))
    coeffsArray3 = np.zeros((Nlist.shape[0],3))
    coeffsArray4 = np.zeros((Nlist.shape[0],4))
    coeffsArray2only = np.zeros((Nlist.shape[0],1))
    coeffsArray24only = np.zeros((Nlist.shape[0],2))
    
    # recover number of timesteps
    numSteps = tFull.shape[0]
    
    # initialize least squares error output
    err = {"t-model" : np.zeros((Nlist.shape[0],1)),
           "t2-model" : np.zeros((Nlist.shape[0],1)),
           "t3-model" : np.zeros((Nlist.shape[0],1)),
           "t4-model" : np.zeros((Nlist.shape[0],1)),
           "t2-model only" : np.zeros((Nlist.shape[0],1)),
           "t2- and t4-models" : np.zeros((Nlist.shape[0],1))}

    # loop through all resolutions
    for j in np.arange(0,Nlist.shape[0]):
    
        # Find number of positive terms in ROM, in intermediate calculations, and wavenumber array
        N = Nlist[j]
        M = Mlist[j]
        k = np.concatenate([np.arange(0,M),np.arange(-M,0)])
    
        # Gather first derivative data for fitting purposes
        exactEnergy = np.zeros((N,numSteps))
        R0Energy = np.zeros((N,numSteps))
        R1Energy = np.zeros((N,numSteps))
        R2Energy = np.zeros((N,numSteps))
        R3Energy = np.zeros((N,numSteps))
        R4Energy = np.zeros((N,numSteps))
    
        # plug exact solution into exact RHS and all ROM terms and find energy contribution of each
        for i in np.arange(0,numSteps):
            
            # exact RHS
            exactRHS,dummyU = markovKdV(uFull[:,i],int(fullM*3/2),alpha)
            exactEnergy[:,i] = np.real(exactRHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(exactRHS[0:N])*uFull[0:N,i])
        
            # Markov RHS
            nonlin0,u_full = markovKdV(uFull[0:N,i],M,alpha)
            R0RHS = nonlin0
            R0Energy[:,i] = np.real(R0RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R0RHS[0:N])*uFull[0:N,i])
        
            # First order RHS term
            F_modes = np.concatenate([np.arange(0,N),np.arange(2*N-1,M+N+2),np.arange(2*M-N+1,2*M)])
            G_modes = np.arange(N,2*M-N+1)
            nonlin1,uuStar = tModelKdV(u_full,nonlin0,alpha,F_modes)
            R1RHS = nonlin1*tFull[i]**(1-tau)
            R1Energy[:,i] = np.real(R1RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R1RHS[0:N])*uFull[0:N,i])
        
            # Second order RHS term
            nonlin2,uk3,uu,A,AStar,B,BStar,C,CStar,D,DStar = t2ModelKdV(u_full,nonlin0,uuStar,alpha,F_modes,G_modes,k,epsilon)
            R2RHS = nonlin2*tFull[i]**(2*(1-tau))
            R2Energy[:,i] = np.real(R2RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R2RHS[0:N])*uFull[0:N,i])
        
            # Third order RHS term
            nonlin3,uk6,E,EStar,F,FStar = t3ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,A,AStar,B,BStar,C,CStar,DStar)
            R3RHS = nonlin3*tFull[i]**(3*(1-tau))
            R3Energy[:,i] = np.real(R3RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R3RHS[0:N])*uFull[0:N,i])
            
            # Fourth order RHS term
            nonlin4 = t4ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,uk6,A,AStar,B,BStar,C,CStar,D,DStar,E,EStar,F,FStar)
            R4RHS = nonlin4*tFull[i]**(4*(1-tau))
            R4Energy[:,i] = np.real(R4RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R4RHS[0:N])*uFull[0:N,i])
        
        if j == 0:
            R0Energy0 = np.copy(R0Energy)
            R1Energy0 = np.copy(R1Energy)
            R2Energy0 = np.copy(R2Energy)
            R3Energy0 = np.copy(R3Energy)
            R4Energy0 = np.copy(R4Energy)
        
        ##################################################
        # Use least-squares fit to identify coefficients #
        ##################################################
        
        # t-model coefficient
        coeffsArray1[j,:] = np.sum((exactEnergy - R0Energy)*R1Energy)/np.sum(R1Energy*R1Energy)
        err["t-model"][j] = np.sum((exactEnergy - R0Energy - coeffsArray1[j,0]*R1Energy)**2)
        
        # t2-model coefficient
        LSMatrix = (np.array([[np.sum(R1Energy*R1Energy),np.sum(R1Energy*R2Energy)],
                             [np.sum(R2Energy*R1Energy),np.sum(R2Energy*R2Energy)]]))
        LSb = (np.array([np.sum(R1Energy*(exactEnergy-R0Energy)),np.sum(R2Energy*(exactEnergy-R0Energy))]))
        coeffsArray2[j,:] = np.linalg.solve(LSMatrix,LSb)
        err["t2-model"][j] = np.sum((exactEnergy - R0Energy - coeffsArray2[j,0]*R1Energy - coeffsArray2[j,1]*R2Energy)**2)
        
        # t3-model coefficient
        LSMatrix = (np.array([[np.sum(R1Energy*R1Energy),np.sum(R1Energy*R2Energy),np.sum(R1Energy*R3Energy)],
                             [np.sum(R2Energy*R1Energy),np.sum(R2Energy*R2Energy),np.sum(R2Energy*R3Energy)],
                             [np.sum(R3Energy*R1Energy),np.sum(R3Energy*R2Energy),np.sum(R3Energy*R3Energy)]]))
        LSb = (np.array([np.sum(R1Energy*(exactEnergy-R0Energy)),np.sum(R2Energy*(exactEnergy-R0Energy)),np.sum(R3Energy*(exactEnergy-R0Energy))]))
        coeffsArray3[j,:] = np.linalg.solve(LSMatrix,LSb)
        err["t3-model"][j] = np.sum((exactEnergy - R0Energy - coeffsArray3[j,0]*R1Energy - coeffsArray3[j,1]*R2Energy - coeffsArray3[j,2]*R3Energy)**2)
        
        # t4-model coefficient
        LSMatrix = (np.array([[np.sum(R1Energy*R1Energy),np.sum(R1Energy*R2Energy),np.sum(R1Energy*R3Energy),np.sum(R1Energy*R4Energy)],
                             [np.sum(R2Energy*R1Energy),np.sum(R2Energy*R2Energy),np.sum(R2Energy*R3Energy),np.sum(R2Energy*R4Energy)],
                             [np.sum(R3Energy*R1Energy),np.sum(R3Energy*R2Energy),np.sum(R3Energy*R3Energy),np.sum(R3Energy*R4Energy)],
                             [np.sum(R4Energy*R1Energy),np.sum(R4Energy*R2Energy),np.sum(R4Energy*R3Energy),np.sum(R4Energy*R4Energy)]]))
        LSb = (np.array([np.sum(R1Energy*(exactEnergy-R0Energy)),np.sum(R2Energy*(exactEnergy-R0Energy)),np.sum(R3Energy*(exactEnergy-R0Energy)),np.sum(R4Energy*(exactEnergy-R0Energy))]))
        coeffsArray4[j,:] = np.linalg.solve(LSMatrix,LSb)
        err["t4-model"][j] = np.sum((exactEnergy - R0Energy - coeffsArray4[j,0]*R1Energy - coeffsArray4[j,1]*R2Energy - coeffsArray4[j,2]*R3Energy - coeffsArray4[j,3]*R4Energy)**2)
        
        # t2-model with *no* t-model
        coeffsArray2only[j,:] = np.sum((exactEnergy - R0Energy)*R2Energy)/np.sum(R2Energy*R2Energy)
        err["t2-model only"][j] = np.sum((exactEnergy - R0Energy - coeffsArray2only[j,0]*R2Energy)**2)
        
        # t2-model and t4-model with *no* t-model or t3-model
        LSMatrix = (np.array([[np.sum(R2Energy*R2Energy),np.sum(R2Energy*R4Energy)],
                             [np.sum(R4Energy*R2Energy),np.sum(R4Energy*R4Energy)]]))
        LSb = (np.array([np.sum(R2Energy*(exactEnergy-R0Energy)),np.sum(R4Energy*(exactEnergy-R0Energy))]))
        coeffsArray24only[j,:] = np.linalg.solve(LSMatrix,LSb)
        err["t2- and t4-models"][j] = np.sum((exactEnergy - R0Energy - coeffsArray24only[j,0]*R2Energy - coeffsArray24only[j,1]*R4Energy)**2)
    
    # Generate plots if desired
    if plots:
        
        # Plot 1: Qualitative comparison of each term contributing to energy movement
        N = Nlist[0]
        fig1, ax1 = plt.subplots(3,2)
        ax1[0,0].plot(tFull,np.sum(exactEnergy[0:N,:],0))
        ax1[0,0].set_title("Exact Energy Decay")
        ax1[0,1].plot(tFull,np.sum(R0Energy0[0:N,:],0))
        ax1[0,1].set_title("Markov Energy Decay")
        ax1[1,0].plot(tFull,np.sum(R2Energy0[0:N,:],0))
        ax1[1,0].set_title("R2 Energy Decay")
        ax1[1,1].plot(tFull,np.sum(R1Energy0[0:N,:],0))
        ax1[1,1].set_title("R1 Energy Decay")
        ax1[2,0].plot(tFull,np.sum(R4Energy0[0:N,:],0))
        ax1[2,0].set_title("R4 Energy Decay")
        ax1[2,1].plot(tFull,np.sum(R3Energy0[0:N,:],0))
        ax1[2,1].set_title("R3 Energy Decay")
    
        fig1.suptitle("N = "+str(N)+" Energy Decays")
        plt.tight_layout()

        # remove axis labels to not crowd plots (since only qualitative comparisons desired)
        for i in range(0,3):
            for j in range(0,2):
                #ax1[i,j].tick_params(labelbottom=False,labelleft=False)
                ax1[i,j].tick_params(labelleft=False)
    
        # compute best fit lines for coefficients in log-log space
        fitLines = {"t-model" : np.zeros((1,3)),
                    "t2-model" : np.zeros((2,3)),
                    "t3-model" : np.zeros((3,3)),
                    "t4-model" : np.zeros((4,3)),
                    "t2-model only" : np.zeros((1,3)),
                    "t2- and t4-models" : np.zeros((2,3))}
    
            
        fig2, ax2 = plt.subplots(2,2)
        # t-model
        ax2[0,0].scatter(np.log(Nlist),np.log(abs(coeffsArray1[:,0])))
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray1[:,0])))
        ax2[0,0].plot(np.log(Nlist),intercept + slope*np.log(Nlist))
        fitLines["t-model"][:] = np.array([slope,np.exp(intercept),r_value])
        
        # t2-model
        ax2[0,0].scatter(np.log(Nlist),np.log(abs(coeffsArray2[:,0])),color="red")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray2[:,0])))
        ax2[0,0].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "red")
        fitLines["t2-model"][0,:] = np.array([slope,np.exp(intercept),r_value])
        
        ax2[0,1].scatter(np.log(Nlist),np.log(abs(coeffsArray2[:,1])),color="red")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray2[:,1])))
        ax2[0,1].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "red")
        fitLines["t2-model"][1,:] = np.array([slope,np.exp(intercept),r_value])
        
        
        # t3-model
        ax2[0,0].scatter(np.log(Nlist),np.log(abs(coeffsArray3[:,0])),color="green")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray3[:,0])))
        ax2[0,0].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "green")
        fitLines["t3-model"][0,:] = np.array([slope,np.exp(intercept),r_value])
        
        ax2[0,1].scatter(np.log(Nlist),np.log(abs(coeffsArray3[:,1])),color="green")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray3[:,1])))
        ax2[0,1].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "green")
        fitLines["t3-model"][1,:] = np.array([slope,np.exp(intercept),r_value])
        
        ax2[1,0].scatter(np.log(Nlist),np.log(abs(coeffsArray3[:,2])),color="green")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray3[:,2])))
        ax2[1,0].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "green")
        fitLines["t3-model"][2,:] = np.array([slope,np.exp(intercept),r_value])
        
        
        # t4-model
        ax2[0,0].scatter(np.log(Nlist),np.log(abs(coeffsArray4[:,0])),color="purple")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray4[:,0])))
        ax2[0,0].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "purple")
        fitLines["t4-model"][0,:] = np.array([slope,np.exp(intercept),r_value])
        
        ax2[0,1].scatter(np.log(Nlist),np.log(abs(coeffsArray4[:,1])),color="purple")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray4[:,1])))
        ax2[0,1].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "purple")
        fitLines["t4-model"][1,:] = np.array([slope,np.exp(intercept),r_value])
        
        ax2[1,0].scatter(np.log(Nlist),np.log(abs(coeffsArray4[:,2])),color="purple")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray4[:,2])))
        ax2[1,0].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "purple")
        fitLines["t4-model"][2,:] = np.array([slope,np.exp(intercept),r_value])
        
        ax2[1,1].scatter(np.log(Nlist),np.log(abs(coeffsArray4[:,3])),color="purple")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray4[:,3])))
        ax2[1,1].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "purple")
        fitLines["t4-model"][3,:] = np.array([slope,np.exp(intercept),r_value])
        
        
        # t2-model alone
        ax2[0,1].scatter(np.log(Nlist),np.log(abs(coeffsArray2only[:,0])),color="cyan")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray2only[:,0])))
        ax2[0,1].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "cyan")
        fitLines["t2-model only"][:] = np.array([slope,np.exp(intercept),r_value])
        
        
        # t2- and t4-model alone
        ax2[0,1].scatter(np.log(Nlist),np.log(abs(coeffsArray24only[:,0])),color="black")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray24only[:,0])))
        ax2[0,1].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "black")
        fitLines["t2- and t4-models"][0,:] = np.array([slope,np.exp(intercept),r_value])
        
        ax2[1,1].scatter(np.log(Nlist),np.log(abs(coeffsArray24only[:,1])),color="black")
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray24only[:,1])))
        ax2[1,1].plot(np.log(Nlist),intercept + slope*np.log(Nlist), color = "black")
        fitLines["t2- and t4-models"][1,:] = np.array([slope,np.exp(intercept),r_value])
    
    
        ax2[0,0].set_title("t-model")
        ax2[0,1].set_title("t2-model")
        ax2[1,0].set_title("t3-model")
        ax2[1,1].set_title("t4-model")
    
        customLines = [plt.Line2D([0],[0], color = "blue"),
                       plt.Line2D([0],[0], color = "red"),
                       plt.Line2D([0],[0], color = "green"),
                       plt.Line2D([0],[0], color = "purple"),
                       plt.Line2D([0],[0], color = "cyan"),
                       plt.Line2D([0],[0], color = "black")]
    
        ax2[0,1].legend(customLines,["First Order Model","Second Order Model",
                                     "Third Order Model","Fourth Order Model",
                                     "Only Second Order","Second and Fourth Order"],
                        prop = {"size":5})
    
        fig2.suptitle("Renormalization Coefficients (log(a) vs log(N))")
        plt.subplots_adjust(right=0.7)
        plt.tight_layout()
        
    # calculate best fit lines if plotting didn't occur
    else:
        fitLines = {"t-model" : np.zeros((1,3)),
                "t2-model" : np.zeros((2,3)),
                "t3-model" : np.zeros((3,3)),
                "t4-model" : np.zeros((4,3)),
                "t2-model only" : np.zeros((1,3)),
                "t2- and t4-models" : np.zeros((2,3))}
        
        # t-model
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray1[:,0])))
        fitLines["t-model"][:] = np.array([slope,np.exp(intercept),r_value])
        
        # second order ROM
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray2[:,0])))
        fitLines["t2-model"][0,:] = np.array([slope,np.exp(intercept),r_value])
        
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray2[:,1])))
        fitLines["t2-model"][1,:] = np.array([slope,np.exp(intercept),r_value])
        
        # third order ROM
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray3[:,0])))
        fitLines["t3-model"][0,:] = np.array([slope,np.exp(intercept),r_value])
        
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray3[:,1])))
        fitLines["t3-model"][1,:] = np.array([slope,np.exp(intercept),r_value])
        
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray3[:,2])))
        fitLines["t3-model"][2,:] = np.array([slope,np.exp(intercept),r_value])
        
           # fourth order ROM 
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray4[:,0])))
        fitLines["t4-model"][0,:] = np.array([slope,np.exp(intercept),r_value])
        
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray4[:,1])))
        fitLines["t4-model"][1,:] = np.array([slope,np.exp(intercept),r_value])
        
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray4[:,2])))
        fitLines["t4-model"][2,:] = np.array([slope,np.exp(intercept),r_value])
        
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray4[:,3])))
        fitLines["t4-model"][3,:] = np.array([slope,np.exp(intercept),r_value])
        
        # only t2-model
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray2only[:,0])))
        fitLines["t2-model only"][:] = np.array([slope,np.exp(intercept),r_value])
        
        # only t2- and t4-models
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray24only[:,0])))
        fitLines["t2- and t4-models"][0,:] = np.array([slope,np.exp(intercept),r_value])
        
        slope,intercept,r_value,p_value,std_err = sp.stats.linregress(np.log(Nlist), np.log(abs(coeffsArray24only[:,1])))
        fitLines["t2- and t4-models"][1,:] = np.array([slope,np.exp(intercept),r_value])
    
    return coeffsArray1,coeffsArray2,coeffsArray3,coeffsArray4,coeffsArray2only,coeffsArray24only,fitLines,err


def scalingLaws(fullM, endtime, Nlist, Mlist, epsilonList, alpha, tau, timesteps, IC = np.sin, plots = False):
    """
    Finds renormalization coefficients based on a simulations with a range of
    epsilon values.
    
    Parameters
    ----------
    fullM : int
            Size of full simulation to base fits on
            
    endtime : int
        Endtime of full simulation
        
    Nlist : list of ints
            List of resolutions for which to find coefficients
            
    Mlist : list of ints
            List of intermediary "full" simulations to use for ROMs
            
    epsilonList : list of floats
        size of linear term (stiffness)
        
    alpha : float
        degree of nonlinearity in KdV
        
    tau : float
        time decay modifier
        
    timesteps : Numpy array
        specific timesteps for which to save solution
        
    IC : function handle
        initial condition of simulation (default np.sin)
        
    plots : boolean
        Indicates whether to generate plots (default: False)
        
    Returns
    -------
    
    coeeffsArray1 : Numpy array (length(Nlist),1)
        Renormalization coefficients for t-model only
        
    coeffsArray2 : Numpy array (length(Nlist),2)
        Renormalization coefficients for t-model and t2-model only
        
    coeffsArray3 : Numpy array (length(Nlist),3)
        Renormalization coefficients for t1-t3-models
        
    coeffsArray4 : Numpy array (length(Nlist),4)
        Renormalization coefficients for t1-t4-models
        
    coeffsArray2only : Numpy array (length(Nlist),1)
        Renormalization coefficients for t2-model only
        
    coeffsArray24only : Numpy array (length(Nlist),2)
        Renormalization coefficients for t2-model and t4-model only
        
    fitLines : Dict
        Contains scaling law fits for each ROM coefficients
        of form   c = -b * N^a
        Terms given are a, b, and r (correlation coefficient of fit)
    """
    
    # initialize output arrays
    c1 = np.zeros((len(Nlist),1,len(epsilonList)))
    c2 = np.zeros((len(Nlist),2,len(epsilonList)))
    c3 = np.zeros((len(Nlist),3,len(epsilonList)))
    c4 = np.zeros((len(Nlist),4,len(epsilonList)))
    c2only = np.zeros((len(Nlist),1,len(epsilonList)))
    c24only = np.zeros((len(Nlist),2,len(epsilonList)))
    
    # loop through all epsilon values
    for i in np.arange(0,len(epsilonList)):
        
        # renormalize for given epsilon value and save results
        coeffsArray1,coeffsArray2,coeffsArray3,coeffsArray4,coeffsArray2only,coeffsArray24only,fitLines,err = renormalize(fullM = fullM, endtime = endtime, Nlist = Nlist, Mlist = Mlist, epsilon = epsilonList[i], alpha = alpha, tau = tau, timesteps = timesteps, IC = IC, plots = False)
        c1[:,:,i] = coeffsArray1
        c2[:,:,i] = coeffsArray2
        c3[:,:,i] = coeffsArray3
        c4[:,:,i] = coeffsArray4
        c2only[:,:,i] = coeffsArray2only
        c24only[:,:,i] = coeffsArray24only
    
    # pack results into dictionary for output    
    coefficients = {"t-model" : c1,
                    "t2-model" : c2,
                    "t3-model" : c3,
                    "t4-model" : c4,
                    "t2-model only" : c2only,
                    "t2- and t4-models" : c24only}
    
    
    # initialize output with best fit scaling laws
    fitLines = {"t-model" : np.zeros((1,3)),
                "t2-model" : np.zeros((2,3)),
                "t3-model" : np.zeros((3,3)),
                "t4-model" : np.zeros((4,3)),
                "t2-model only" : np.zeros((1,3)),
                "t2- and t4-models" : np.zeros((2,3))}
    
    # find the scaling laws for each coefficient
    
    # t-model coefficient
    fitLines["t-model"][0,:] = epsilonNscalingLaw(c1[:,0,:],Nlist,epsilonList)
    
    # Second order model coefficients
    fitLines["t2-model"][0,:] = epsilonNscalingLaw(c2[:,0,:],Nlist,epsilonList)
    fitLines["t2-model"][1,:] = epsilonNscalingLaw(c2[:,1,:],Nlist,epsilonList)
    
    # Third order model coefficients
    fitLines["t3-model"][0,:] = epsilonNscalingLaw(c3[:,0,:],Nlist,epsilonList)
    fitLines["t3-model"][1,:] = epsilonNscalingLaw(c3[:,1,:],Nlist,epsilonList)
    fitLines["t3-model"][2,:] = epsilonNscalingLaw(c3[:,2,:],Nlist,epsilonList)
    
    # Fourth order model coefficients
    fitLines["t4-model"][0,:] = epsilonNscalingLaw(c4[:,0,:],Nlist,epsilonList)
    fitLines["t4-model"][1,:] = epsilonNscalingLaw(c4[:,1,:],Nlist,epsilonList)
    fitLines["t4-model"][2,:] = epsilonNscalingLaw(c4[:,2,:],Nlist,epsilonList)
    fitLines["t4-model"][3,:] = epsilonNscalingLaw(c4[:,3,:],Nlist,epsilonList)
    
    # Only t2-model coefficient
    fitLines["t2-model only"][0,:] = epsilonNscalingLaw(c2only[:,0,:],Nlist,epsilonList)
    
    # Only t2- and t4-models coefficients
    fitLines["t2- and t4-models"][0,:] = epsilonNscalingLaw(c24only[:,0,:],Nlist,epsilonList)
    fitLines["t2- and t4-models"][1,:] = epsilonNscalingLaw(c24only[:,1,:],Nlist,epsilonList)
    
    
    # make plots
    fig1,ax1 = plt.subplots(1,2)
    fig2,ax2 = plt.subplots(2,2)
    fig3,ax3 = plt.subplots(3,2)
    fig4,ax4 = plt.subplots(4,2)
    fig5,ax5 = plt.subplots(1,2)
    fig6,ax6 = plt.subplots(2,2)

    # loop through epsilon values
    for i in np.arange(len(epsilonList)):
        
        # t-model coefficient
        ax1[0].scatter(np.log(Nlist),np.log(-c1[:,0,i]))

        # Second order model coefficients
        ax2[0,0].scatter(np.log(Nlist),np.log(-c2[:,0,i]))
        ax2[1,0].scatter(np.log(Nlist),np.log(-c2[:,1,i]))
        
        # Third order model coefficients
        ax3[0,0].scatter(np.log(Nlist),np.log(-c3[:,0,i]))
        ax3[1,0].scatter(np.log(Nlist),np.log(-c3[:,1,i]))
        ax3[2,0].scatter(np.log(Nlist),np.log(-c3[:,2,i]))
        
        # Fourth order model coefficients
        ax4[0,0].scatter(np.log(Nlist),np.log(-c4[:,0,i]))
        ax4[1,0].scatter(np.log(Nlist),np.log(-c4[:,1,i]))
        ax4[2,0].scatter(np.log(Nlist),np.log(-c4[:,2,i]))
        ax4[3,0].scatter(np.log(Nlist),np.log(-c4[:,3,i]))
        
        # Only t2-model
        ax5[0].scatter(np.log(Nlist),np.log(-c2only[:,0,i]))
        
        # Only t2- and t4-models
        ax6[0,0].scatter(np.log(Nlist),np.log(-c24only[:,0,i]))
        ax6[1,0].scatter(np.log(Nlist),np.log(-c24only[:,1,i]))
    
        
        # plot best fit lines
        myEps = epsilonList[i]
        
        myFit = fitLines["t-model"][0,:]
        ax1[0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        
        myFit = fitLines["t2-model"][0,:]
        ax2[0,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        myFit = fitLines["t2-model"][1,:]
        ax2[1,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))


        myFit = fitLines["t3-model"][0,:]
        ax3[0,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        myFit = fitLines["t3-model"][1,:]
        ax3[1,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        myFit = fitLines["t3-model"][2,:]
        ax3[2,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        
        myFit = fitLines["t4-model"][0,:]
        ax4[0,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        myFit = fitLines["t4-model"][1,:]
        ax4[1,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        myFit = fitLines["t4-model"][2,:]
        ax4[2,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        myFit = fitLines["t4-model"][3,:]
        ax4[3,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
    
        myFit = fitLines["t2-model only"][0,:]
        ax5[0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        
        myFit = fitLines["t2- and t4-models"][0,:]
        ax6[0,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))
        myFit = fitLines["t2- and t4-models"][1,:]
        ax6[1,0].plot(np.log(Nlist),np.log(-myFit[0])+myFit[1]*np.log(Nlist)+myFit[2]*np.log(myEps))

    

    
    # loop through epsilon values
    for j in np.arange(len(Nlist)):
        
        # t-model coefficient
        ax1[1].scatter(np.log(epsilonList),np.log(-c1[j,0,:]))
        
        # Second order model coefficients
        ax2[0,1].scatter(np.log(epsilonList),np.log(-c2[j,0,:]))
        ax2[1,1].scatter(np.log(epsilonList),np.log(-c2[j,1,:]))
            
        # Third order model coefficients
        ax3[0,1].scatter(np.log(epsilonList),np.log(-c3[j,0,:]))
        ax3[1,1].scatter(np.log(epsilonList),np.log(-c3[j,1,:]))
        ax3[2,1].scatter(np.log(epsilonList),np.log(-c3[j,2,:]))
            
        # Fourth order model coefficients
        ax4[0,1].scatter(np.log(epsilonList),np.log(-c4[j,0,:]))
        ax4[1,1].scatter(np.log(epsilonList),np.log(-c4[j,1,:]))
        ax4[2,1].scatter(np.log(epsilonList),np.log(-c4[j,2,:]))
        ax4[3,1].scatter(np.log(epsilonList),np.log(-c4[j,3,:]))
        
        # Only t2-model
        ax5[1].scatter(np.log(epsilonList),np.log(-c2only[j,0,:]))
        
        # Only t2- and t4-models
        ax6[0,1].scatter(np.log(epsilonList),np.log(-c24only[j,0,:]))
        ax6[1,1].scatter(np.log(epsilonList),np.log(-c24only[j,1,:]))
        
        
        # plot best fit lines
        myN = Nlist[j]
            
        myFit = fitLines["t-model"][0,:]
        ax1[1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        
        myFit = fitLines["t2-model"][0,:]
        ax2[0,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        myFit = fitLines["t2-model"][1,:]
        ax2[1,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        
        
        myFit = fitLines["t3-model"][0,:]
        ax3[0,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        myFit = fitLines["t3-model"][1,:]
        ax3[1,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        myFit = fitLines["t3-model"][2,:]
        ax3[2,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        
        myFit = fitLines["t4-model"][0,:]
        ax4[0,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        myFit = fitLines["t4-model"][1,:]
        ax4[1,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        myFit = fitLines["t4-model"][2,:]
        ax4[2,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        myFit = fitLines["t4-model"][3,:]
        ax4[3,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        
        myFit = fitLines["t2-model only"][0,:]
        ax5[1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        
        myFit = fitLines["t2- and t4-models"][0,:]
        ax6[0,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        myFit = fitLines["t2- and t4-models"][1,:]
        ax6[1,1].plot(np.log(epsilonList),np.log(-myFit[0])+myFit[1]*np.log(myN)+myFit[2]*np.log(epsilonList))
        
        
    # label all plots
    fig1.suptitle("t-model")
    ax1[0].set_title("log(a1) vs log(N)")
    ax1[0].legend(["epsilon = "+str(round(epsilonList[i]),2) for i in range(len(epsilonList))],prop = {"size":5})
    ax1[1].set_title("log(a1) vs log(epsilon)")
    ax1[1].legend(["N = "+str(Nlist[i]) for i in range(len(Nlist))],prop = {"size":5})
    plt.tight_layout()
    
    fig2.suptitle("Second Order Renormalization")
    ax2[0,0].set_title("log(a1) vs log(N)")
    ax2[1,0].set_title("log(a2) vs log(N)")
    ax2[0,0].legend(["epsilon = "+str(round(epsilonList[i]),2) for i in range(len(epsilonList))],prop = {"size":5})
    ax2[0,1].set_title("log(a1) vs log(epsilon)")
    ax2[1,1].set_title("log(a1) vs log(epsilon)")
    ax2[0,1].legend(["N = "+str(Nlist[i]) for i in range(len(Nlist))],prop = {"size":5})
    plt.tight_layout()
    
    fig3.suptitle("Third Order Renormalization")
    ax3[0,0].set_title("log(a1) vs log(N)")
    ax3[1,0].set_title("log(a2) vs log(N)")
    ax3[2,0].set_title("log(a3) vs log(N)")
    ax3[0,0].legend(["epsilon = "+str(round(epsilonList[i]),2) for i in range(len(epsilonList))],prop = {"size":5})
    ax3[0,1].set_title("log(a1) vs log(epsilon)")
    ax3[1,1].set_title("log(a2) vs log(epsilon)")
    ax3[2,1].set_title("log(a3) vs log(epsilon)")
    ax3[0,1].legend(["N = "+str(Nlist[i]) for i in range(len(Nlist))],prop = {"size":5})
    plt.tight_layout()
    
    fig4.suptitle("Fourth Order Renormalization")
    ax4[0,0].set_title("log(a1) vs log(N)")
    ax4[1,0].set_title("log(a2) vs log(N)")
    ax4[2,0].set_title("log(a3) vs log(N)")
    ax4[3,0].set_title("log(a4) vs log(N)")
    ax4[0,0].legend(["epsilon = "+str(round(epsilonList[i]),2) for i in range(len(epsilonList))],prop = {"size":5})
    ax4[0,1].set_title("log(a1) vs log(epsilon)")
    ax4[1,1].set_title("log(a2) vs log(epsilon)")
    ax4[2,1].set_title("log(a3) vs log(epsilon)")
    ax4[3,1].set_title("log(a4) vs log(epsilon)")
    ax4[0,1].legend(["N = "+str(Nlist[i]) for i in range(len(Nlist))],prop = {"size":5})
    plt.tight_layout()
        
    fig5.suptitle("Only t2-Model Renormalization")
    ax5[0].set_title("log(a2) vs log(N)")
    ax5[0].legend(["epsilon = "+str(round(epsilonList[i]),2) for i in range(len(epsilonList))],prop = {"size":5})
    ax5[1].set_title("log(a2) vs log(epsilon)")
    ax5[1].legend(["N = "+str(Nlist[i]) for i in range(len(Nlist))],prop = {"size":5})
    plt.tight_layout()
        
    fig6.suptitle("Second and Fourth Order Renormalization")
    ax6[0,0].set_title("log(a2) vs log(N)")
    ax6[1,0].set_title("log(a4) vs log(N)")
    ax6[0,0].legend(["epsilon = "+str(round(epsilonList[i]),2) for i in range(len(epsilonList))],prop = {"size":5})
    ax6[0,1].set_title("log(a2) vs log(epsilon)")
    ax6[1,1].set_title("log(a4) vs log(epsilon)")
    ax6[0,1].legend(["N = "+str(Nlist[i]) for i in range(len(Nlist))],prop = {"size":5})
    plt.tight_layout()
    
    return coefficients,fitLines

def epsilonNscalingLaw(coeffArray,Nlist,epsilonList):
    
    numEps = len(epsilonList)
    numN = len(Nlist)
    epsilonTile = np.tile(epsilonList,(numN,1))
    Ntile = np.transpose(np.tile(Nlist,(numEps,1)))

    LSMatrix = (np.array([[numEps*numN,np.sum(np.log(Ntile)),np.sum(np.log(epsilonTile))],
                     [np.sum(np.log(Ntile)),np.sum(np.log(Ntile)**2),np.sum(np.log(Ntile)*np.log(epsilonTile))],
                     [np.sum(np.log(epsilonTile)),np.sum(np.log(Ntile)*np.log(epsilonTile)),np.sum(np.log(epsilonTile)**2)]])
                )

    LSb = np.array([np.sum(np.log(np.abs(coeffArray))),np.sum(np.log(np.abs(coeffArray))*np.log(Ntile)),np.sum(np.log(np.abs(coeffArray))*np.log(epsilonTile))])

    sol = np.linalg.solve(LSMatrix,LSb)
    sol[0] = -np.exp(sol[0])
    return sol


def findError(compareList,exact,t):
    """
    Finds the two norm of the error between a list of ROMs and an exact solution.
    
    Parameters
    ----------
    compareList : List of Numpy arrays of size (N,T)
            Set of state vector evolutions to find errors from
            
    exact : Numpy array of size (N,T)
        Exact solution for the same timesteps
            
    t : Numpy array (T,)
        Timesteps associated with simulations (must all be the same)
        
    Returns
    -------
    errList : List of Numpy arrays of size (T,1)
        Arrays with the two-norm of the error at all timesteps for each ROM
    """
    
    # find the ROM size
    N = compareList[0].shape[0]
    
    # generate real space solutions
    realSols = [makeRealSpace(x,N) for x in compareList]
    exactSol = makeRealSpace(exact,N)
    
    # compute two norm of error at all times
    errList =[np.sum((i - exactSol)**2,0) for i in realSols]
    
    return errList


def renormalizeRobust(fullM, endtime, Nlist, Mlist, epsilon, alpha, tau, timesteps, IC = np.sin, plots = False):
    """
    Finds renormalization coefficients based on a single simulation. If the
    simulation doesn't yet exist, it creates it
    
    Parameters
    ----------
    fullM : int
            Size of full simulation to base fits on
            
    endtime : int
        Endtime of full simulation
        
    Nlist : list of ints
            List of resolutions for which to find coefficients
            
    Mlist : list of ints
            List of intermediary "full" simulations to use for ROMs
            
    epsilon : float
        size of linear term (stiffness)
        
    alpha : float
        degree of nonlinearity in KdV
        
    tau : float
        time decay modifier
        
    timesteps : Numpy array
        specific timesteps for which to save solution
        
    IC : function handle
        initial condition of simulation (default np.sin)
        
    plots : boolean
        Indicates whether to generate plots (default: False)
        
    Returns
    -------
    
    coeeffsArray1 : Numpy array (length(Nlist),1)
        Renormalization coefficients for t-model only
        
    coeffsArray2 : Numpy array (length(Nlist),2)
        Renormalization coefficients for t-model and t2-model only
        
    coeffsArray3 : Numpy array (length(Nlist),3)
        Renormalization coefficients for t1-t3-models
        
    coeffsArray4 : Numpy array (length(Nlist),4)
        Renormalization coefficients for t1-t4-models
        
    coeffsArray2only : Numpy array (length(Nlist),1)
        Renormalization coefficients for t2-model only
        
    coeffsArray24only : Numpy array (length(Nlist),2)
        Renormalization coefficients for t2-model and t4-model only
        
    fitLines : Dict
        Contains scaling law fits for each ROM coefficients
        of form   c = -b * N^a
        Terms given are a, b, and r (correlation coefficient of fit)
    """
    
    # Check if full simulation has already been constructed
    #   if so, load it, if not, generate it
    try:
        uFull = np.load("u" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+IC.__name__+".npy")
        tFull = np.load("t" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+IC.__name__+".npy")
    except:
        fullParams = {
            'N': fullM,
            'M': int(3/2*fullM),
            'alpha': 1,
            'epsilon': epsilon,
            'tau': 1,
            'coeffs': None,
            'IC': IC,
            'endtime': endtime,
            'timesteps': timesteps
            }

        uSimFull = runSim(fullParams)
        uFull = uSimFull.y
        tFull = uSimFull.t
        np.save( "u" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+IC.__name__,uFull)
        np.save( "t" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+IC.__name__,tFull)

    
    # recover number of timesteps
    numSteps = tFull.shape[0]
    
    # initialize output arrays
    coeffsArray1 = np.zeros((Nlist.shape[0],numSteps - 30,1))
    coeffsArray2 = np.zeros((Nlist.shape[0],numSteps - 30,2))
    coeffsArray3 = np.zeros((Nlist.shape[0],numSteps - 30,3))
    coeffsArray4 = np.zeros((Nlist.shape[0],numSteps - 30,4))
    
    coeffsArray2only = np.zeros((Nlist.shape[0],numSteps - 30,1))
    coeffsArray24only = np.zeros((Nlist.shape[0],numSteps - 30,2))
     

    # loop through all resolutions
    for j in np.arange(0,Nlist.shape[0]):
    
        # Find number of positive terms in ROM, in intermediate calculations, and wavenumber array
        N = Nlist[j]
        M = Mlist[j]
        k = np.concatenate([np.arange(0,M),np.arange(-M,0)])
    
        # Gather first derivative data for fitting purposes
        exactEnergy = np.zeros((N,numSteps))
        R0Energy = np.zeros((N,numSteps))
        R1Energy = np.zeros((N,numSteps))
        R2Energy = np.zeros((N,numSteps))
        R3Energy = np.zeros((N,numSteps))
        R4Energy = np.zeros((N,numSteps))
    
        # plug exact solution into exact RHS and all ROM terms and find energy contribution of each
        for i in np.arange(0,numSteps):
            
            # exact RHS
            exactRHS,dummyU = markovKdV(uFull[:,i],int(fullM*3/2),alpha)
            exactEnergy[:,i] = np.real(exactRHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(exactRHS[0:N])*uFull[0:N,i])
        
            # Markov RHS
            nonlin0,u_full = markovKdV(uFull[0:N,i],M,alpha)
            R0RHS = nonlin0
            R0Energy[:,i] = np.real(R0RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R0RHS[0:N])*uFull[0:N,i])
        
            # First order RHS term
            F_modes = np.concatenate([np.arange(0,N),np.arange(2*N-1,M+N+2),np.arange(2*M-N+1,2*M)])
            G_modes = np.arange(N,2*M-N+1)
            nonlin1,uuStar = tModelKdV(u_full,nonlin0,alpha,F_modes)
            R1RHS = nonlin1*tFull[i]**(1-tau)
            R1Energy[:,i] = np.real(R1RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R1RHS[0:N])*uFull[0:N,i])
        
            # Second order RHS term
            nonlin2,uk3,uu,A,AStar,B,BStar,C,CStar,D,DStar = t2ModelKdV(u_full,nonlin0,uuStar,alpha,F_modes,G_modes,k,epsilon)
            R2RHS = nonlin2*tFull[i]**(2*(1-tau))
            R2Energy[:,i] = np.real(R2RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R2RHS[0:N])*uFull[0:N,i])
        
            # Third order RHS term
            nonlin3,uk6,E,EStar,F,FStar = t3ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,A,AStar,B,BStar,C,CStar,DStar)
            R3RHS = nonlin3*tFull[i]**(3*(1-tau))
            R3Energy[:,i] = np.real(R3RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R3RHS[0:N])*uFull[0:N,i])
            
            # Fourth order RHS term
            nonlin4 = t4ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,uk6,A,AStar,B,BStar,C,CStar,D,DStar,E,EStar,F,FStar)
            R4RHS = nonlin4*tFull[i]**(4*(1-tau))
            R4Energy[:,i] = np.real(R4RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R4RHS[0:N])*uFull[0:N,i])

        
        ##################################################
        # Use least-squares fit to identify coefficients #
        ##################################################
        
        for i in np.arange(30,numSteps):
            
            exactEnergySnip = exactEnergy[:,0:i]
            R0EnergySnip = R0Energy[:,0:i]
            R1EnergySnip = R1Energy[:,0:i]
            R2EnergySnip = R2Energy[:,0:i]
            R3EnergySnip = R3Energy[:,0:i]
            R4EnergySnip = R4Energy[:,0:i]
            
            # t-model coefficient
            coeffsArray1[j,i-30,:] = np.sum((exactEnergySnip - R0EnergySnip)*R1EnergySnip)/np.sum(R1EnergySnip*R1EnergySnip)
        
            # t2-model coefficient
            LSMatrix = (np.array([[np.sum(R1EnergySnip*R1EnergySnip),np.sum(R1EnergySnip*R2EnergySnip)],
                                  [np.sum(R2EnergySnip*R1EnergySnip),np.sum(R2EnergySnip*R2EnergySnip)]]))
            LSb = (np.array([np.sum(R1EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R2EnergySnip*(exactEnergySnip-R0EnergySnip))]))
            coeffsArray2[j,i-30,:] = np.linalg.solve(LSMatrix,LSb)
        
            # t3-model coefficient
            LSMatrix = (np.array([[np.sum(R1EnergySnip*R1EnergySnip),np.sum(R1EnergySnip*R2EnergySnip),np.sum(R1EnergySnip*R3EnergySnip)],
                                  [np.sum(R2EnergySnip*R1EnergySnip),np.sum(R2EnergySnip*R2EnergySnip),np.sum(R2EnergySnip*R3EnergySnip)],
                                  [np.sum(R3EnergySnip*R1EnergySnip),np.sum(R3EnergySnip*R2EnergySnip),np.sum(R3EnergySnip*R3EnergySnip)]]))
            LSb = (np.array([np.sum(R1EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R2EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R3EnergySnip*(exactEnergySnip-R0EnergySnip))]))
            coeffsArray3[j,i-30,:] = np.linalg.solve(LSMatrix,LSb)
        
            # t4-model coefficient
            LSMatrix = (np.array([[np.sum(R1EnergySnip*R1EnergySnip),np.sum(R1EnergySnip*R2EnergySnip),np.sum(R1EnergySnip*R3EnergySnip),np.sum(R1EnergySnip*R4EnergySnip)],
                                  [np.sum(R2EnergySnip*R1EnergySnip),np.sum(R2EnergySnip*R2EnergySnip),np.sum(R2EnergySnip*R3EnergySnip),np.sum(R2EnergySnip*R4EnergySnip)],
                                  [np.sum(R3EnergySnip*R1EnergySnip),np.sum(R3EnergySnip*R2EnergySnip),np.sum(R3EnergySnip*R3EnergySnip),np.sum(R3EnergySnip*R4EnergySnip)],
                                  [np.sum(R4EnergySnip*R1EnergySnip),np.sum(R4EnergySnip*R2EnergySnip),np.sum(R4EnergySnip*R3EnergySnip),np.sum(R4EnergySnip*R4EnergySnip)]]))
            LSb = (np.array([np.sum(R1EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R2EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R3EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R4EnergySnip*(exactEnergySnip-R0EnergySnip))]))
            coeffsArray4[j,i-30,:] = np.linalg.solve(LSMatrix,LSb)
        
            # t2-model with *no* t-model
            coeffsArray2only[j,i-30,:] = np.sum((exactEnergySnip - R0EnergySnip)*R2EnergySnip)/np.sum(R2EnergySnip*R2EnergySnip)
        
            # t2-model and t4-model with *no* t-model or t3-model
            LSMatrix = (np.array([[np.sum(R2EnergySnip*R2EnergySnip),np.sum(R2EnergySnip*R4EnergySnip)],
                                  [np.sum(R4EnergySnip*R2EnergySnip),np.sum(R4EnergySnip*R4EnergySnip)]]))
            LSb = (np.array([np.sum(R2EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R4EnergySnip*(exactEnergySnip-R0EnergySnip))]))
            coeffsArray24only[j,i-30,:] = np.linalg.solve(LSMatrix,LSb)
    
    for ind in np.arange(Nlist.shape[0]):
        fig1,ax1 = plt.subplots(2,2)
        fig1.suptitle("N = "+str(Nlist[ind]))
        ax1[0,0].plot(timesteps[30:],coeffsArray1[ind,:,0],color = "blue")
        ax1[0,0].plot(timesteps[30:],coeffsArray2[ind,:,0],color = "red")
        ax1[0,0].plot(timesteps[30:],coeffsArray3[ind,:,0],color = "green")
        ax1[0,0].plot(timesteps[30:],coeffsArray4[ind,:,0],color = "black")
        ax1[0,0].set_title("t-model")

        ax1[0,1].plot([],[],color = "blue")
        ax1[0,1].plot(timesteps[30:],coeffsArray2[ind,:,1], color = "red")
        ax1[0,1].plot(timesteps[30:],coeffsArray3[ind,:,1], color = "green")
        ax1[0,1].plot(timesteps[30:],coeffsArray4[ind,:,1], color = "black")
        ax1[0,1].plot(timesteps[30:],coeffsArray2only[ind,:,0],color = "cyan")
        ax1[0,1].plot(timesteps[30:],coeffsArray24only[ind,:,0], color = "magenta")
        ax1[0,1].set_title("t2-model")
        ax1[0,1].legend(["First order","Second order","Third order","Fourth order","Only t2","t2 and t4"],prop = {"size":5})

        ax1[1,0].plot(timesteps[30:],coeffsArray3[ind,:,2], color = "green")
        ax1[1,0].plot(timesteps[30:],coeffsArray4[ind,:,2],color = "black")
        ax1[1,0].set_title("t3-model")
    
        ax1[1,1].plot(timesteps[30:],coeffsArray4[ind,:,3], color = "black")
        ax1[1,1].plot(timesteps[30:],coeffsArray24only[ind,:,1], color = "magenta")
        ax1[1,1].set_title("t4-model")
        plt.tight_layout()
    
    return coeffsArray1,coeffsArray2,coeffsArray3,coeffsArray4,coeffsArray2only,coeffsArray24only

def renormalizeWindow(fullM, endtime, width, Nlist, Mlist, epsilon, alpha, tau, timesteps, IC = np.sin, plots = False):
    """
    Finds renormalization coefficients using sliding window least squares.
    
    Parameters
    ----------
    fullM : int
            Size of full simulation to base fits on
            
    endtime : int
        Endtime of full simulation
        
    width : float
        Size of sliding window to use in fitting
        
    Nlist : list of ints
            List of resolutions for which to find coefficients
            
    Mlist : list of ints
            List of intermediary "full" simulations to use for ROMs
            
    epsilon : float
        size of linear term (stiffness)
        
    alpha : float
        degree of nonlinearity in KdV
        
    tau : float
        time decay modifier
        
    timesteps : Numpy array
        specific timesteps for which to save solution
        
    IC : function handle
        initial condition of simulation (default np.sin)
        
    plots : boolean
        Indicates whether to generate plots (default: False)
        
    Returns
    -------
    
    coeeffsArray1 : Numpy array (length(Nlist),1)
        Renormalization coefficients for t-model only
        
    coeffsArray2 : Numpy array (length(Nlist),2)
        Renormalization coefficients for t-model and t2-model only
        
    coeffsArray3 : Numpy array (length(Nlist),3)
        Renormalization coefficients for t1-t3-models
        
    coeffsArray4 : Numpy array (length(Nlist),4)
        Renormalization coefficients for t1-t4-models
        
    coeffsArray2only : Numpy array (length(Nlist),1)
        Renormalization coefficients for t2-model only
        
    coeffsArray24only : Numpy array (length(Nlist),2)
        Renormalization coefficients for t2-model and t4-model only
        
    fitLines : Dict
        Contains scaling law fits for each ROM coefficients
        of form   c = -b * N^a
        Terms given are a, b, and r (correlation coefficient of fit)
    """
    
    # Check if full simulation has already been constructed
    #   if so, load it, if not, generate it
    try:
        uFull = np.load("u" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+IC.__name__+".npy")
        tFull = np.load("t" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+IC.__name__+".npy")
    except:
        fullParams = {
            'N': fullM,
            'M': int(3/2*fullM),
            'alpha': 1,
            'epsilon': epsilon,
            'tau': 1,
            'coeffs': None,
            'IC': IC,
            'endtime': endtime,
            'timesteps': timesteps
            }

        uSimFull = runSim(fullParams)
        uFull = uSimFull.y
        tFull = uSimFull.t
        np.save( "u" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+IC.__name__,uFull)
        np.save( "t" + str(fullM) + "t" + str(endtime)+"e"+str(round(epsilon,2)).replace('.','p')+IC.__name__,tFull)

    
    # recover number of timesteps
    numSteps = tFull.shape[0]
    
    widthSteps = round(width/(tFull[1]-tFull[0]))
    
    # initialize output arrays
    coeffsArray1 = np.zeros((Nlist.shape[0],numSteps - widthSteps+1,1))
    coeffsArray2 = np.zeros((Nlist.shape[0],numSteps - widthSteps+1,2))
    coeffsArray3 = np.zeros((Nlist.shape[0],numSteps - widthSteps+1,3))
    coeffsArray4 = np.zeros((Nlist.shape[0],numSteps - widthSteps+1,4))
    
    coeffsArray2only = np.zeros((Nlist.shape[0],numSteps - widthSteps+1,1))
    coeffsArray24only = np.zeros((Nlist.shape[0],numSteps - widthSteps+1,2))
    
    exact1 = np.zeros((Nlist.shape[0],1))
    exact2 = np.zeros((Nlist.shape[0],2))
    exact3 = np.zeros((Nlist.shape[0],3))
    exact4 = np.zeros((Nlist.shape[0],4))
    exact2o = np.zeros((Nlist.shape[0],1))
    exact24o = np.zeros((Nlist.shape[0],2))
     

    # loop through all resolutions
    for j in np.arange(0,Nlist.shape[0]):
    
        # Find number of positive terms in ROM, in intermediate calculations, and wavenumber array
        N = Nlist[j]
        M = Mlist[j]
        k = np.concatenate([np.arange(0,M),np.arange(-M,0)])
    
        # Gather first derivative data for fitting purposes
        exactEnergy = np.zeros((N,numSteps))
        R0Energy = np.zeros((N,numSteps))
        R1Energy = np.zeros((N,numSteps))
        R2Energy = np.zeros((N,numSteps))
        R3Energy = np.zeros((N,numSteps))
        R4Energy = np.zeros((N,numSteps))
    
        # plug exact solution into exact RHS and all ROM terms and find energy contribution of each
        for i in np.arange(0,numSteps):
            
            # exact RHS
            exactRHS,dummyU = markovKdV(uFull[:,i],int(fullM*3/2),alpha)
            exactEnergy[:,i] = np.real(exactRHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(exactRHS[0:N])*uFull[0:N,i])
        
            # Markov RHS
            nonlin0,u_full = markovKdV(uFull[0:N,i],M,alpha)
            R0RHS = nonlin0
            R0Energy[:,i] = np.real(R0RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R0RHS[0:N])*uFull[0:N,i])
        
            # First order RHS term
            F_modes = np.concatenate([np.arange(0,N),np.arange(2*N-1,M+N+2),np.arange(2*M-N+1,2*M)])
            G_modes = np.arange(N,2*M-N+1)
            nonlin1,uuStar = tModelKdV(u_full,nonlin0,alpha,F_modes)
            if tFull[i] == 0:
                R1RHS = nonlin1*0
            else:
                R1RHS = nonlin1*tFull[i]**(1-tau)
            R1Energy[:,i] = np.real(R1RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R1RHS[0:N])*uFull[0:N,i])
        
            # Second order RHS term
            nonlin2,uk3,uu,A,AStar,B,BStar,C,CStar,D,DStar = t2ModelKdV(u_full,nonlin0,uuStar,alpha,F_modes,G_modes,k,epsilon)
            R2RHS = nonlin2*tFull[i]**(2*(1-tau))
            R2Energy[:,i] = np.real(R2RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R2RHS[0:N])*uFull[0:N,i])
        
            # Third order RHS term
            nonlin3,uk6,E,EStar,F,FStar = t3ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,A,AStar,B,BStar,C,CStar,DStar)
            R3RHS = nonlin3*tFull[i]**(3*(1-tau))
            R3Energy[:,i] = np.real(R3RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R3RHS[0:N])*uFull[0:N,i])
            
            # Fourth order RHS term
            nonlin4 = t4ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,uk6,A,AStar,B,BStar,C,CStar,D,DStar,E,EStar,F,FStar)
            R4RHS = nonlin4*tFull[i]**(4*(1-tau))
            R4Energy[:,i] = np.real(R4RHS[0:N]*np.conj(uFull[0:N,i]) + np.conj(R4RHS[0:N])*uFull[0:N,i])

        
        exact1[j,:] = np.sum((exactEnergy - R0Energy)*R1Energy)/np.sum(R1Energy*R1Energy)

        LSMatrix = (np.array([[np.sum(R1Energy*R1Energy),np.sum(R1Energy*R2Energy)],
                                  [np.sum(R2Energy*R1Energy),np.sum(R2Energy*R2Energy)]]))
        LSb = (np.array([np.sum(R1Energy*(exactEnergy-R0Energy)),np.sum(R2Energy*(exactEnergy-R0Energy))]))
        exact2[j,:] = np.linalg.solve(LSMatrix,LSb)
        
        LSMatrix = (np.array([[np.sum(R1Energy*R1Energy),np.sum(R1Energy*R2Energy),np.sum(R1Energy*R3Energy)],
                                  [np.sum(R2Energy*R1Energy),np.sum(R2Energy*R2Energy),np.sum(R2Energy*R3Energy)],
                                  [np.sum(R3Energy*R1Energy),np.sum(R3Energy*R2Energy),np.sum(R3Energy*R3Energy)]]))
        LSb = (np.array([np.sum(R1Energy*(exactEnergy-R0Energy)),np.sum(R2Energy*(exactEnergy-R0Energy)),np.sum(R3Energy*(exactEnergy-R0Energy))]))
        exact3[j,:] = np.linalg.solve(LSMatrix,LSb)
        
        LSMatrix = (np.array([[np.sum(R1Energy*R1Energy),np.sum(R1Energy*R2Energy),np.sum(R1Energy*R3Energy),np.sum(R1Energy*R4Energy)],
                            [np.sum(R2Energy*R1Energy),np.sum(R2Energy*R2Energy),np.sum(R2Energy*R3Energy),np.sum(R2Energy*R4Energy)],
                            [np.sum(R3Energy*R1Energy),np.sum(R3Energy*R2Energy),np.sum(R3Energy*R3Energy),np.sum(R3Energy*R4Energy)],
                            [np.sum(R4Energy*R1Energy),np.sum(R4Energy*R2Energy),np.sum(R4Energy*R3Energy),np.sum(R4Energy*R4Energy)]]))
        LSb = (np.array([np.sum(R1Energy*(exactEnergy-R0Energy)),np.sum(R2Energy*(exactEnergy-R0Energy)),np.sum(R3Energy*(exactEnergy-R0Energy)),np.sum(R4Energy*(exactEnergy-R0Energy))]))
        exact4[j,:] = np.linalg.solve(LSMatrix,LSb)
        

        exact2o[j,:] = np.sum((exactEnergy - R0Energy)*R2Energy)/np.sum(R2Energy*R2Energy)
        
          
        LSMatrix = (np.array([[np.sum(R2Energy*R2Energy),np.sum(R2Energy*R4Energy)],
                             [np.sum(R4Energy*R2Energy),np.sum(R4Energy*R4Energy)]]))
        LSb = (np.array([np.sum(R2Energy*(exactEnergy-R0Energy)),np.sum(R4Energy*(exactEnergy-R0Energy))]))
        exact24o[j,:] = np.linalg.solve(LSMatrix,LSb)
            
        ##################################################
        # Use least-squares fit to identify coefficients #
        ##################################################
        
        for i in np.arange(0,numSteps-widthSteps+1):
            
            exactEnergySnip = exactEnergy[:,i:i+widthSteps]
            R0EnergySnip = R0Energy[:,i:i+widthSteps]
            R1EnergySnip = R1Energy[:,i:i+widthSteps]
            R2EnergySnip = R2Energy[:,i:i+widthSteps]
            R3EnergySnip = R3Energy[:,i:i+widthSteps]
            R4EnergySnip = R4Energy[:,i:i+widthSteps]
            
            
            # t-model coefficient
            coeffsArray1[j,i,:] = np.sum((exactEnergySnip - R0EnergySnip)*R1EnergySnip)/np.sum(R1EnergySnip*R1EnergySnip)
        
            # t2-model coefficient
            LSMatrix = (np.array([[np.sum(R1EnergySnip*R1EnergySnip),np.sum(R1EnergySnip*R2EnergySnip)],
                                  [np.sum(R2EnergySnip*R1EnergySnip),np.sum(R2EnergySnip*R2EnergySnip)]]))
            LSb = (np.array([np.sum(R1EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R2EnergySnip*(exactEnergySnip-R0EnergySnip))]))
            coeffsArray2[j,i,:] = np.linalg.solve(LSMatrix,LSb)
        
            # t3-model coefficient
            LSMatrix = (np.array([[np.sum(R1EnergySnip*R1EnergySnip),np.sum(R1EnergySnip*R2EnergySnip),np.sum(R1EnergySnip*R3EnergySnip)],
                                  [np.sum(R2EnergySnip*R1EnergySnip),np.sum(R2EnergySnip*R2EnergySnip),np.sum(R2EnergySnip*R3EnergySnip)],
                                  [np.sum(R3EnergySnip*R1EnergySnip),np.sum(R3EnergySnip*R2EnergySnip),np.sum(R3EnergySnip*R3EnergySnip)]]))
            LSb = (np.array([np.sum(R1EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R2EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R3EnergySnip*(exactEnergySnip-R0EnergySnip))]))
            coeffsArray3[j,i,:] = np.linalg.solve(LSMatrix,LSb)
        
            # t4-model coefficient
            LSMatrix = (np.array([[np.sum(R1EnergySnip*R1EnergySnip),np.sum(R1EnergySnip*R2EnergySnip),np.sum(R1EnergySnip*R3EnergySnip),np.sum(R1EnergySnip*R4EnergySnip)],
                                  [np.sum(R2EnergySnip*R1EnergySnip),np.sum(R2EnergySnip*R2EnergySnip),np.sum(R2EnergySnip*R3EnergySnip),np.sum(R2EnergySnip*R4EnergySnip)],
                                  [np.sum(R3EnergySnip*R1EnergySnip),np.sum(R3EnergySnip*R2EnergySnip),np.sum(R3EnergySnip*R3EnergySnip),np.sum(R3EnergySnip*R4EnergySnip)],
                                  [np.sum(R4EnergySnip*R1EnergySnip),np.sum(R4EnergySnip*R2EnergySnip),np.sum(R4EnergySnip*R3EnergySnip),np.sum(R4EnergySnip*R4EnergySnip)]]))
            LSb = (np.array([np.sum(R1EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R2EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R3EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R4EnergySnip*(exactEnergySnip-R0EnergySnip))]))
            coeffsArray4[j,i,:] = np.linalg.solve(LSMatrix,LSb)
        
            # t2-model with *no* t-model
            coeffsArray2only[j,i,:] = np.sum((exactEnergySnip - R0EnergySnip)*R2EnergySnip)/np.sum(R2EnergySnip*R2EnergySnip)
        
            # t2-model and t4-model with *no* t-model or t3-model
            LSMatrix = (np.array([[np.sum(R2EnergySnip*R2EnergySnip),np.sum(R2EnergySnip*R4EnergySnip)],
                                  [np.sum(R4EnergySnip*R2EnergySnip),np.sum(R4EnergySnip*R4EnergySnip)]]))
            LSb = (np.array([np.sum(R2EnergySnip*(exactEnergySnip-R0EnergySnip)),np.sum(R4EnergySnip*(exactEnergySnip-R0EnergySnip))]))
            coeffsArray24only[j,i,:] = np.linalg.solve(LSMatrix,LSb)
    
    for ind in np.arange(Nlist.shape[0]):
        fig1,ax1 = plt.subplots(2,2)
        fig1.suptitle("N = "+str(Nlist[ind]))
        ax1[0,0].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray1[ind,:,0],color = "blue")
        ax1[0,0].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray2[ind,:,0],color = "red")
        ax1[0,0].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray3[ind,:,0],color = "green")
        ax1[0,0].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray4[ind,:,0],color = "black")
        ax1[0,0].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact1[ind,0],exact1[ind,0]],color="blue",linestyle=":")
        ax1[0,0].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact2[ind,0],exact2[ind,0]],color="red",linestyle=":")
        ax1[0,0].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact3[ind,0],exact3[ind,0]],color="green",linestyle=":")
        ax1[0,0].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact4[ind,0],exact4[ind,0]],color="black",linestyle=":")
        ax1[0,0].set_title("t-model")

        ax1[0,1].plot([],[],color = "blue")
        ax1[0,1].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray2[ind,:,1], color = "red")
        ax1[0,1].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray3[ind,:,1], color = "green")
        ax1[0,1].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray4[ind,:,1], color = "black")
        ax1[0,1].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray2only[ind,:,0],color = "cyan")
        ax1[0,1].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray24only[ind,:,0], color = "magenta")
        ax1[0,1].set_title("t2-model")
        ax1[0,1].legend(["First order","Second order","Third order","Fourth order","Only t2","t2 and t4"],prop = {"size":5})
        ax1[0,1].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact2[ind,1],exact2[ind,1]],color="red",linestyle=":")
        ax1[0,1].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact3[ind,1],exact3[ind,1]],color="green",linestyle=":")
        ax1[0,1].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact4[ind,1],exact4[ind,1]],color="black",linestyle=":")
        ax1[0,1].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact2o[ind,0],exact2o[ind,0]],color="cyan",linestyle=":")
        ax1[0,1].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact24o[ind,0],exact24o[ind,0]],color="magenta",linestyle=":")

        ax1[1,0].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray3[ind,:,2], color = "green")
        ax1[1,0].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray4[ind,:,2],color = "black")
        ax1[1,0].set_title("t3-model")
        ax1[1,0].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact3[ind,2],exact3[ind,2]],color="green",linestyle=":")
        ax1[1,0].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact4[ind,2],exact4[ind,2]],color="black",linestyle=":")
    
        ax1[1,1].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray4[ind,:,3], color = "black")
        ax1[1,1].plot(timesteps[0:numSteps-widthSteps+1],coeffsArray24only[ind,:,1], color = "magenta")
        ax1[1,1].set_title("t4-model")
        ax1[1,1].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact4[ind,3],exact4[ind,3]],color="black",linestyle=":")
        ax1[1,1].plot([timesteps[0],timesteps[numSteps-widthSteps+1]],[exact24o[ind,1],exact24o[ind,1]],color="magenta",linestyle=":")
        plt.tight_layout()
    
    return coeffsArray1,coeffsArray2,coeffsArray3,coeffsArray4,coeffsArray2only,coeffsArray24only


def renormalizeTau(fullM, endtime, Nlist, Mlist, epsilon, alpha, tauList, timesteps, IC = np.sin):
    """
    Tests a range of tau values for fitting coefficients.
    
    Parameters
    ----------
    fullM : int
        Resolution of full model to use for fitting
        
    endtime : float
        Final time to use for fitting
        
    epsilon : float
        size of linear term (stiffness)
        
    alpha : float
        degree of nonlinearity in KdV
        
    timesteps : Numpy array
        specific timesteps for which to save solution
        
    tauList : float
        Grid of tau values to test (default 0:0.05:1)
        
    IC : function handle
        initial condition of simulation (default np.sin)
     
        
    Returns
    -------
    out : dict
        Contains optimal coefficients for each model for each value of tau.
        't-model' is len(tauList) x len(Nlist) x 1,
        't2-model' is len(tauList) x len(Nlist) x 2, etc.
        
    err : 
    
    
    
    """
    
    out = {"t-model" : np.zeros((tauList.shape[0],Nlist.shape[0],1)),
           "t2-model" : np.zeros((tauList.shape[0],Nlist.shape[0],2)),
           "t3-model" : np.zeros((tauList.shape[0],Nlist.shape[0],3)),
           "t4-model" : np.zeros((tauList.shape[0],Nlist.shape[0],4)),
           "t2-model only" : np.zeros((tauList.shape[0],Nlist.shape[0],1)),
           "t2- and t4-models" : np.zeros((tauList.shape[0],Nlist.shape[0],2))
           }
    
    errList = []
    
    for i in np.arange(tauList.shape[0]):
        result = renormalize(fullM, endtime, Nlist, Mlist, epsilon, alpha, tauList[i], timesteps, IC = IC, plots = False)
        out["t-model"][i,:,:] = result[0]
        out["t2-model"][i,:,:] = result[1]
        out["t3-model"][i,:,:] = result[2]
        out["t4-model"][i,:,:] = result[3]
        out["t2-model only"][i,:,:] = result[4]
        out["t2- and t4-models"][i,:,:] = result[5]
        
        errList.append(result[7])
        
    return out,errList

def automatedROM(N,alpha,epsilon,timesteps,fitTime = 10,tauTests = np.arange(0,1.05,0.05),IC = np.sin,tol = 1e-3):
    """
    Automatically finds optimal tau and coefficients for an ROM and runs the ROM.
    Also produces reference exact solution
    
    Parameters
    ----------
    N : int
        Resolution of ROM
        
    alpha : float
        degree of nonlinearity in KdV
        
    epsilon : float
        size of linear term (stiffness)
        
    timesteps : Numpy array
        specific timesteps for which to save solution
            
    fitTime : float
        Fits are made over window of exact solution from 0 to fitTime (default 10)
        
    tauTests : float
        Grid of tau values to test (default 0:0.05:1)
        
    IC : function handle
        initial condition of simulation (default np.sin)
        
    tol : float
        Tolerance for declaring full model "resolved"
        There must be less than this relative error in the first half of the 
        full modes up to the end time (default 10^-3)
     
        
    Returns
    -------
    
    simMarkov : SciPy integration object
        Simulation up to end time of Markov model
        
    sim2 : SciPy integration object
        Simulation up to end time of 2nd order model model
        
    sim4 : SciPy integration object
        Simulation up to end time of fourth order model
        
    coefficients : 
    
    errors
    
    """
    
    endtime = timesteps[-1]
    
    M = 16
    unresolved = True
    #print("Constructing reference exact solution...")
    
    try:
        fileList = glob.glob("u" + '[0-9]*' + "t" + str(int(endtime))+"e"+str(round(epsilon,2)).replace('.','p')+".npy")
        myFile = fileList[0]
        uFull = np.load(myFile)
        
        fileList = glob.glob("t" + '[0-9]*' + "t" + str(int(endtime))+"e"+str(round(epsilon,2)).replace('.','p')+".npy")
        myFile = fileList[0]
        tFull = np.load(myFile)
        
        #print("Success! (it was already saved)\n")
        
        M = int(re.findall(r'\d+', fileList[0])[0])
    except:
        # find resolved simulation    
        while unresolved:
        
            M = 2*M
            #print("\nCurrently testing M = "+str(M))
            
            fullParams = {
                'N': M,
                'M': int(3/2*M),
                'alpha': 1,
                'epsilon': epsilon,
                'tau': 1,
                'coeffs': None,
                'IC': IC,
                'endtime': endtime,
                'timesteps': timesteps
                }

            uSimFull = runSim(fullParams)
            uFull = uSimFull.y
            tFull = uSimFull.t
            energyCheck = getMass(uFull,int(M/2))
            
            print("Maximum mass deviation in first "+str(M/2)+" modes: "+str(max(abs(energyCheck - energyCheck[0]))))
            
            if max(abs(energyCheck - energyCheck[0])) < tol:
                print("Success!\n")
                unresolved = False
                np.save("u" + str(M) + "t" + str(int(endtime))+"e"+str(round(epsilon,2)).replace('.','p')+".npy",uFull)
                np.save("t" + str(M) + "t" + str(int(endtime))+"e"+str(round(epsilon,2)).replace('.','p')+".npy",tFull)
        
    uFit = uFull[:,tFull<=fitTime]
    tFit = tFull[tFull<=fitTime]
    np.save("u" + str(M) + "t" + str(int(fitTime))+"e"+str(round(epsilon,2)).replace('.','p')+".npy",uFit)
    np.save("t" + str(M) + "t" + str(int(fitTime))+"e"+str(round(epsilon,2)).replace('.','p')+".npy",tFit)
    
    # Find specific fitting window
    
    #print("Finding coefficients for range of tau values...")
    coefficients,errors = renormalizeTau(fullM = M,
                                         endtime = fitTime,
                                         Nlist = np.array([N]),
                                         Mlist = np.array([N*3]),
                                         epsilon = epsilon,
                                         alpha = alpha,
                                         tauList = tauTests,
                                         timesteps = tFit,
                                         IC = IC)
    
    err2 = [fit["t2-model only"][0][0] for fit in errors]
    c2 = coefficients["t2-model only"][err2.index(min(err2))][0]
    
    err4 = [fit["t2- and t4-models"][0][0] for fit in errors]
    c4 = coefficients["t2- and t4-models"][err4.index(min(err4))][0]

    
    #print("Coefficients found!\n")
    
    paramsMarkov = {
                'N': N,
                'M': int(3*N),
                'alpha': 1,
                'epsilon': epsilon,
                'tau': 0,
                'coeffs': None,
                'IC': IC,
                'endtime': endtime,
                'timesteps': timesteps
                }
    
    params2 = {
                'N': N,
                'M': int(3*N),
                'alpha': 1,
                'epsilon': epsilon,
                'tau': tauTests[err2.index(min(err2))],
                'coeffs': np.array([0,c2[0]]),
                'IC': IC,
                'endtime': endtime,
                'timesteps': timesteps
                }
    
    params4 = {
                'N': N,
                'M': int(3*N),
                'alpha': 1,
                'epsilon': epsilon,
                'tau': tauTests[err4.index(min(err4))],
                'coeffs': np.array([0,c4[0],0,c4[1]]),
                'IC': IC,
                'endtime': endtime,
                'timesteps': timesteps
                }
    
    params2scaling = {
                'N': N,
                'M': int(3*N),
                'alpha': 1,
                'epsilon': epsilon,
                'tau': 1,
                'coeffs': np.array([0,-0.7615*N**-5.8081*epsilon**-3.7681]), # NOTE GET THE LATEST SCALING LAWS HERE
                'IC': IC,
                'endtime': endtime,
                'timesteps': timesteps
                }
    
    params4scaling = {
                'N': N,
                'M': int(3*N),
                'alpha': 1,
                'epsilon': epsilon,
                'tau': 1,
                'coeffs': np.array([0,-1.2473*N**-5.7356*epsilon**-3.6910,0,-0.3675*N**-11.4719*epsilon**-7.3881]), # NOTE GET THE LATEST SCALING LAWS HERE
                'IC': IC,
                'endtime': endtime,
                'timesteps': timesteps
                }
    
    #print("Running Markov simulation...\n")
    simMarkov = runSim(paramsMarkov)
    
    #print("Running second order simulation with tau = "+str(tauTests[err2.index(min(err2))])+"...\n")
    sim2 = runSim(params2)
    
    #print("Running fourth order simulation with tau = "+str(tauTests[err4.index(min(err4))])+"...")
    sim4 = runSim(params4)
    
    sim2scale = runSim(params2scaling)
    
    sim4scale = runSim(params4scaling)
    
    return simMarkov,sim2,sim4,sim2scale,sim4scale,coefficients,errors


# coeffsArray1,coeffsArray2,coeffsArray3,coeffsArray4,coeffsArray2only,coeffsArray24only,fitLines

# old functions that should be defunct

# def RHSFullKdV(u,M,alpha,epsilon):
#     """Computes RHS of KdV with no reduced order model and dealiasing
    
#     Parameters
#     ----------
#     u : 1D Numpy Array (M,)
#         Positive modes of state vector whose RHS is being computed
        
#     M : int
#         Number of positive modes in "full" model for intermediary calculations
        
#     alpha : float
#         Degree of nonlinearity in KdV
        
#     epsilon : float
#         Size of linear term (stiffness)

#     Returns
#     -------
#     RHS : 1D Numpy Array (M,)
#         RHS of ODE for state vector u
#     """
    
#     # compute nonlinear part of RHS using "full" model with M positive modes
#     # the factor of 3/2 is for dealiasing
#     a,b = markovKdV(u,int(3/2*M),alpha)
    
#     # compute the linear part of the RHS and add it to the nonlinear part
#     L = u.shape[0]
#     k = np.arange(0,L)
#     RHS = 1j*k**3*epsilon**2*u + a[0:L]
#     return RHS

# def RHSMarkovKdV(u,N,M,alpha,epsilon):
#     """Computes RHS of KdV with no reduced order model and dealiasing
    
#     Parameters
#     ----------
#     u : 1D Numpy Array (N,)
#         Positive modes of state vector whose RHS is being computed
        
#     N : int
#         Number of positive modes in ROM
        
#     M : int
#         Number of positive modes in "full" model for intermediary calculations
        
#     alpha : float
#         Degree of nonlinearity in KdV
        
#     epsilon : float
#         Size of linear term (stiffness)

#     Returns
#     -------
#     RHS : 1D Numpy Array (N,)
#         RHS of ODE for state vector u
#     """
    
#     # compute nonlinear part of RHS using "full" model with M positive modes
#     a,b = markovKdV(u,M,alpha)
    
#     # compute the linear part of the RHS and add it to the nonlinear part
#     L = u.shape[0]
#     k = np.arange(0,L)
#     RHS = 1j*k**3*epsilon**2*u + a[0:L]
#     return RHS

# def RHStModelKdV(u,N,M,alpha,epsilon,t):
#     """Computes RHS of KdV for the t-model with coefficient one
    
#     Parameters
#     ----------
#     u : 1D Numpy Array (N,)
#         Positive modes of state vector whose RHS is being computed
        
#     N : int
#         Number of positive modes in ROM
        
#     M : int
#         Number of positive modes in "full" model for intermediary calculations
        
#     alpha : float
#         Degree of nonlinearity in KdV
        
#     epsilon : float
#         Size of linear term (stiffness)
        
#     t  : float
#          Current timestep

#     Returns
#     -------
#     RHS : 1D Numpy Array (N,)
#         RHS of ODE for state vector u
#     """
    
#     # compute nonlinear part of RHS using "full" model with M positive modes
#     nonlin0,u_full = markovKdV(u,M,alpha)
    
#     # define which modes are resolved in full array
#     F_modes = np.concatenate([np.arange(0,N),np.arange(2*N-1,M+N+2),np.arange(2*M-N+1,2*M)])
    
#     # compute t-model term
#     nonlin1,uuStar = tModelKdV(u_full,nonlin0,alpha,F_modes)
    
#     # compute the linear part of the RHS and add it to the nonlinear part
#     L = u.shape[0]
#     k = np.arange(0,L)
    
#     # combine linear, Markov, and t-model terms
#     RHS = 1j*k**3*epsilon**2*u + nonlin0[0:L] + nonlin1[0:L]*t
#     return RHS

# def RHSt2ModelKdV(u,N,M,alpha,epsilon,t,tau,coeffs):
#     """Computes RHS of KdV for the t2-model with variable coefficients
    
#     Parameters
#     ----------
#     u : 1D Numpy Array (N,)
#         Positive modes of state vector whose RHS is being computed
        
#     N : int
#         Number of positive modes in ROM
        
#     M : int
#         Number of positive modes in "full" model for intermediary calculations
        
#     alpha : float
#         Degree of nonlinearity in KdV
        
#     epsilon : float
#         Size of linear term (stiffness)
        
#     t  : float
#          Current timestep
         
#     tau : float
#           Time decay modifier (default 0 - time dependence of memory fully retained)
         
#     coeffs : 1D Numpy Array (2,1)
#         Renormalization coefficients for t-model and t2-model

#     Returns
#     -------
#     RHS : 1D Numpy Array (N,)
#         RHS of ODE for state vector u
#     """
    
#     # compute nonlinear part of RHS using "full" model with M positive modes
#     nonlin0,u_full = markovKdV(u,M,alpha)
    
#     # define which modes are resolved / unresolved in full array
#     F_modes = np.concatenate([np.arange(0,N),np.arange(2*N-1,M+N+2),np.arange(2*M-N+1,2*M)])
#     G_modes = np.arange(N,2*M-N+1)
    
#     # compute t-model term
#     nonlin1,uuStar = tModelKdV(u_full,nonlin0,alpha,F_modes)
    
#     k = np.concatenate([np.arange(0,M),np.arange(-M,0)])
    
#     # compute t2-model term
#     nonlin2,uk3,uu,A,AStar,B,BStar,C,CStar,D,DStar = t2ModelKdV(u_full,nonlin0,uuStar,alpha,F_modes,G_modes,k,epsilon)
    
#     # combine linear, Markov, t-model, and t2-model terms
#     RHS = 1j*k[0:N]**3*epsilon**2*u + nonlin0[0:N] + coeffs[0]*nonlin1[0:N]*t**(1-tau) + coeffs[1]*nonlin2[0:N]*t**(2*(1-tau))
#     return RHS

# def RHSt4ModelKdV(u,N,M,alpha,epsilon,t,tau,coeffs):
#     """Computes RHS of KdV for the t4-model with variable coefficients
    
#     Parameters
#     ----------
#     u : 1D Numpy Array (N,)
#         Positive modes of state vector whose RHS is being computed
        
#     N : int
#         Number of positive modes in ROM
        
#     M : int
#         Number of positive modes in "full" model for intermediary calculations
        
#     alpha : float
#         Degree of nonlinearity in KdV
        
#     epsilon : float
#         Size of linear term (stiffness)
        
#     t  : float
#          Current timestep
         
#     tau : float
#           Time decay modifier (default 0 - time dependence of memory fully retained)
         
#     coeffs : 1D Numpy Array (2,1)
#         Renormalization coefficients for t-model and t2-model

#     Returns
#     -------
#     RHS : 1D Numpy Array (N,)
#         RHS of ODE for state vector u
#     """
    
#     # compute nonlinear part of RHS using "full" model with M positive modes
#     nonlin0,u_full = markovKdV(u,M,alpha)
    
#     # define which modes are resolved / unresolved in full array
#     F_modes = np.concatenate([np.arange(0,N),np.arange(2*N-1,M+N+2),np.arange(2*M-N+1,2*M)])
#     G_modes = np.arange(N,2*M-N+1)
    
#     # compute t-model term
#     nonlin1,uuStar = tModelKdV(u_full,nonlin0,alpha,F_modes)
    
#     k = np.concatenate([np.arange(0,M),np.arange(-M,0)])
    
#     # compute t2-model term
#     nonlin2,uk3,uu,A,AStar,B,BStar,C,CStar,D,DStar = t2ModelKdV(u_full,nonlin0,uuStar,alpha,F_modes,G_modes,k,epsilon)
    
#     # compute t3-model term
#     nonlin3,uk6,E,EStar,F,FStar = t3ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,A,AStar,B,BStar,C,CStar,DStar)
    
#     # compute t4-model term
#     nonlin4 = t4ModelKdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uuStar,uk3,uk6,A,AStar,B,BStar,C,CStar,D,DStar,E,EStar,F,FStar)
    
#     # combine linear, Markov, t-model, and t2-model terms
#     RHS = (1j*k[0:N]**3*epsilon**2*u + nonlin0[0:N] + coeffs[0]*nonlin1[0:N]*t**(1-tau)
#     + coeffs[1]*nonlin2[0:N]*t**(2*(1-tau))
#     + coeffs[2]*nonlin3[0:N]*t**(3*(1-tau))
#     + coeffs[3]*nonlin4[0:N]*t**(4*(1-tau)))
#     return RHS


# def nonRenormSim(N, endtime = 10, alpha = 1, epsilon = 0.1, IC = np.sin, timesteps = None):
#     """Runs a non-renormalized simulation of KdV - use large N so it is resolved
    
#     Parameters
#     ----------
#     N : int
#         Number of positive modes in simulation
        
#     endtime : float
#         Final time to run simulation to (default 10)
        
#     alpha : float
#         Degree of nonlinearity in KdV (default 1)
        
#     epsilon : float
#         Size of linear term (stiffness) (default 0.1)
        
#     IC : function handle
#         Initial condition (default sin(x))
        
        
#     Returns
#     -------
#     uSim : ODE Solution object
#         Data from completed simulation
#     """
    
#     # generate initial condition
#     x = np.linspace(0,2*np.pi-2*np.pi/(2*N),2*N)
#     y = IC(x)
#     uFull = fftnorm(y)
#     u = uFull[0:N]
    
#     # define RHS for integration
#     # (RHS FullKdV uses "full" size with 3/2*N positive modes to dealias results)
#     def myRHS(t,y):
#         out = RHSFullKdV(y,N,alpha,epsilon)
#         return out
    
#     # run simulation using stiff solver
#     uSim = sp.integrate.solve_ivp(fun = myRHS, t_span = [0,endtime], y0 = u,method = "BDF", t_eval = timesteps)
#     return uSim



# def markovSim(N, M = None, endtime = 10, alpha = 1, epsilon = 0.1, IC = np.sin, timesteps = None):
#     """Runs a non-renormalized simulation of KdV - use large N so it is resolved
    
#     Parameters
#     ----------
#     N : int
#         Number of positive modes in simulation
        
#     M : int
#         Number of positive modes in "full" model for intermediary calculations (default 2N)
        
#     endtime : float
#         Final time to run simulation to (default 10)
        
#     alpha : float
#         Degree of nonlinearity in KdV (default 1)
        
#     epsilon : float
#         Size of linear term (stiffness) (default 0.1)
        
#     IC : function handle
#         Initial condition (default sin(x))
        
        
#     Returns
#     -------
#     uSim : ODE Solution object
#         Data from completed simulation
#     """
    
#     # default value of M is double N
#     if M is None:
#         M = 2*N
    
#     # generate initial condition
#     x = np.linspace(0,2*np.pi-2*np.pi/(2*N),2*N)
#     y = IC(x)
#     uFull = fftnorm(y)
#     u = uFull[0:N]
    
#     # define RHS for integration
#     # (RHS FullKdV uses "full" size with 3/2*N positive modes to dealias results)
#     def myRHS(t,y):
#         out = RHSMarkovKdV(y,N,M,alpha,epsilon)
#         return out
    
#     # run simulation using stiff solver
#     uSim = sp.integrate.solve_ivp(fun = myRHS, t_span = [0,endtime], y0 = u,method = "BDF", t_eval = timesteps)
#     return uSim

# def tModelSim(N, M = None, endtime = 10, alpha = 1, epsilon = 0.1, IC = np.sin, timesteps = None):
#     """Runs a basic t-model simulation of KdV
    
#     Parameters
#     ----------
#     N : int
#         Number of positive modes in simulation
        
#     M : int
#         Number of positive modes in "full" model for intermediary calculations (default 3N)
        
#     endtime : float
#         Final time to run simulation to (default 10)
        
#     alpha : float
#         Degree of nonlinearity in KdV (default 1)
        
#     epsilon : float
#         Size of linear term (stiffness) (default 0.1)
        
#     IC : function handle
#         Initial condition (default sin(x))
        
#     timesteps : Numpy array (t,1)
#         Specific timesteps to save
        
        
#     Returns
#     -------
#     uSim : ODE Solution object
#         Data from completed simulation
#     """
    
#     # default value of M is triple N
#     if M is None:
#         M = 3*N
    
#     # generate initial condition
#     x = np.linspace(0,2*np.pi-2*np.pi/(2*N),2*N)
#     y = IC(x)
#     uFull = fftnorm(y)
#     u = uFull[0:N]
    
#     # define RHS for integration
#     # (RHS FullKdV uses "full" size with 3/2*N positive modes to dealias results)
#     def myRHS(t,y):
#         out = RHStModelKdV(y,N,M,alpha,epsilon,t)
#         return out
    
#     # run simulation using stiff solver
#     uSim = sp.integrate.solve_ivp(fun = myRHS, t_span = [0,endtime], y0 = u,method = "BDF", t_eval = timesteps)
#     return uSim

# def t2ModelSim(N, M = None, endtime = 10, alpha = 1, epsilon = 0.1, tau = 0, IC = np.sin, timesteps = None, coeffs = np.array([1,-0.5])):
#     """Runs a t2-model simulation of KdV
    
#     Parameters
#     ----------
#     N : int
#         Number of positive modes in simulation
        
#     M : int
#         Number of positive modes in "full" model for intermediary calculations (default 3N)
        
#     endtime : float
#         Final time to run simulation to (default 10)
        
#     alpha : float
#         Degree of nonlinearity in KdV (default 1)
        
#     epsilon : float
#         Size of linear term (stiffness) (default 0.1)
        
#     tau : float
#         Time decay modifier (default 0 - time dependence of memory fully retained)
        
#     IC : function handle
#         Initial condition (default sin(x))
        
#     timesteps : Numpy array (t,1)
#         Specific timesteps to save
        
#     coeffs : Numpy array (2,1)
#         Renormalization coefficients on t-model and t2-model terms (default [1,-0.5])
        
        
#     Returns
#     -------
#     uSim : ODE Solution object
#         Data from completed simulation
#     """
    
#     # default value of M is triple N
#     if M is None:
#         M = 3*N
    
#     # generate initial condition
#     x = np.linspace(0,2*np.pi-2*np.pi/(2*N),2*N)
#     y = IC(x)
#     uFull = fftnorm(y)
#     u = uFull[0:N]
    
#     # define RHS for integration
#     # (RHS FullKdV uses "full" size with 3/2*N positive modes to dealias results)
#     def myRHS(t,y):
#         out = RHSt2ModelKdV(y,N,M,alpha,epsilon,t,tau,coeffs)
#         return out
    
#     # run simulation using stiff solver
#     uSim = sp.integrate.solve_ivp(fun = myRHS, t_span = [0,endtime], y0 = u,method = "BDF", t_eval = timesteps)
#     return uSim

# def t4ModelSim(N, M = None, endtime = 10, alpha = 1, epsilon = 0.1, tau = 0, IC = np.sin, timesteps = None, coeffs = np.array([1,-0.5,1/6,-1/24])):
#     """Runs a t2-model simulation of KdV
    
#     Parameters
#     ----------
#     N : int
#         Number of positive modes in simulation
        
#     M : int
#         Number of positive modes in "full" model for intermediary calculations (default 3N)
        
#     endtime : float
#         Final time to run simulation to (default 10)
        
#     alpha : float
#         Degree of nonlinearity in KdV (default 1)
        
#     epsilon : float
#         Size of linear term (stiffness) (default 0.1)
        
#     tau : float
#         Time decay modifier (default 0 - time dependence of memory fully retained)
        
#     IC : function handle
#         Initial condition (default sin(x))
        
#     timesteps : Numpy array (t,1)
#         Specific timesteps to save
        
#     coeffs : Numpy array (4,1)
#         Renormalization coefficients on ROM terms (default [1,-0.5,1/6,-1/24])
        
        
#     Returns
#     -------
#     uSim : ODE Solution object
#         Data from completed simulation
#     """
    
#     # default value of M is triple N
#     if M is None:
#         M = 3*N
    
#     # generate initial condition
#     x = np.linspace(0,2*np.pi-2*np.pi/(2*N),2*N)
#     y = IC(x)
#     uFull = fftnorm(y)
#     u = uFull[0:N]
    
#     # define RHS for integration
#     # (RHS FullKdV uses "full" size with 3/2*N positive modes to dealias results)
#     def myRHS(t,y):
#         out = RHSt4ModelKdV(y,N,M,alpha,epsilon,t,tau,coeffs)
#         return out
    
#     # run simulation using stiff solver
#     uSim = sp.integrate.solve_ivp(fun = myRHS, t_span = [0,endtime], y0 = u,method = "BDF", t_eval = timesteps)
#     return uSim

# def makeAnimation(u,t,N):
#     """Takes a completed simulation and creates a real space animation of a subset of modes
    
#     Parameters
#     ----------
#     u : Numpy array (M,t)
#         Output of simulation giving energy in first M positive modes for all timesteps t
    
#     t : Numpy array (t,1)
#         Timesteps associated with the simulation
    
#     N : int
#         Number of positive modes to use in real space
        
        
#     Returns
#     -------
#     anim : Animation object
#             Real space solution animated
#     """
    
#     # generate real space solutions
#     xgrid,yarray = makeRealSpace(u,N)
    
#     # initialize figure
#     myFig = plt.figure()
#     ax = plt.subplot(111)
#     ax.axis(xmin = 0,xmax = 2*np.pi,ymin = -2, ymax = 3)
#     myLine, = ax.plot([],[], 'b')

#     # define function to draw each frame
#     def makeFrame(n):
#         myLine.set_data(xgrid,yarray[:,n])
#         plt.title('t = '+str(round(t[n],1)))
#         return myLine

#     # generate animation
#     anim = animation.FuncAnimation(fig = myFig,func = makeFrame,frames = t.shape[0])
    
#     return anim



# def makeAnimationFullMarkov(uFull,uMarkov,t):
#     """Takes a completed full and Markov simulation and creates a real space animation comparing them
    
#     Parameters
#     ----------
#     uFull : Numpy array (M,t)
#         Output of full simulation giving energy in first M positive modes for all timesteps t
    
#     uMarkov : Numpy array (M,t)
#         Output of Markov ROM simulation of size N for all timesteps t
    
#     t : Numpy array (t,1)
#         Timesteps associated with the simulation
        
        
#     Returns
#     -------
#     anim : Animation object
#             Real space solution animated
#     """
    
#     N = uMarkov.shape[0]
    
#     # generate real space solutions
#     xgrid,yarrayFull = makeRealSpace(uFull,N)
#     xgrid,yarrayMarkov = makeRealSpace(uMarkov,N)
    
#     # initialize figure
#     myFig = plt.figure()
#     ax = plt.subplot(111)
#     ax.axis(xmin = 0,xmax = 2*np.pi-np.pi/N,ymin = -2, ymax = 3)
#     myFullLine, = ax.plot([],[], 'b')
#     myMarkovLine, = ax.plot([],[], 'r')

#     # define function to draw each frame
#     def makeFrame(n):
#         myFullLine.set_data(xgrid,yarrayFull[:,n])
#         myMarkovLine.set_data(xgrid,yarrayMarkov[:,n])
#         plt.title('t = '+str(round(t[n],1)))
#         plt.legend((myFullLine,myMarkovLine),("Exact Solution","Markov Model"))
#         return myFullLine,myMarkovLine

#     # generate animation
#     anim = animation.FuncAnimation(fig = myFig,func = makeFrame,frames = t.shape[0])
    
#     return anim

# def makeAnimationtModel(uFull,uMarkov,utModel,t):
#     """Takes a completed full, Markov, and t-model simulation and creates a real space animation comparing them
    
#     Parameters
#     ----------
#     uFull : Numpy array (M,t)
#         Output of full simulation giving energy in first M positive modes for all timesteps t
    
#     uMarkov : Numpy array (M,t)
#         Output of Markov ROM simulation of size N for all timesteps t
        
#     utmodel : Numpy array (M,t)
#         Output of t-model ROM simulation of size N for all timesteps t
    
#     t : Numpy array (t,1)
#         Timesteps associated with the simulation
        
        
#     Returns
#     -------
#     anim : Animation object
#             Real space solution animated
#     """
    
#     N = uMarkov.shape[0]
    
#     # generate real space solutions
#     xgrid,yarrayFull = makeRealSpace(uFull,N)
#     xgrid,yarrayMarkov = makeRealSpace(uMarkov,N)
#     xgrid,yarraytModel = makeRealSpace(utModel,N)
    
#     # initialize figure
#     myFig = plt.figure()
#     ax = plt.subplot(111)
#     ax.axis(xmin = 0,xmax = 2*np.pi-np.pi/N,ymin = -2, ymax = 3)
#     myFullLine, = ax.plot([],[], 'b')
#     myMarkovLine, = ax.plot([],[], 'r')
#     mytModelLine, = ax.plot([],[], 'g')

#     # define function to draw each frame
#     def makeFrame(n):
#         myFullLine.set_data(xgrid,yarrayFull[:,n])
#         myMarkovLine.set_data(xgrid,yarrayMarkov[:,n])
#         mytModelLine.set_data(xgrid,yarraytModel[:,n])
#         plt.title('t = '+str(round(t[n],1)))
#         plt.legend((myFullLine,myMarkovLine,mytModelLine),("Exact Solution","Markov Model","t-Model"))
#         return myFullLine,myMarkovLine

#     # generate animation
#     anim = animation.FuncAnimation(fig = myFig,func = makeFrame,frames = t.shape[0])
    
#     return anim

# def makeAnimationROMS(uFull,uMarkov,utModel,ut2Model,ut4Model,t):
#     """
#     """
    
#     N = uMarkov.shape[0]
    
    
    
#     # generate real space solutions
#     xgrid,yarrayFull = makeRealSpace(uFull,N)
#     xgrid,yarrayMarkov = makeRealSpace(uMarkov,N)
#     xgrid,yarraytModel = makeRealSpace(utModel,N)
#     xgrid,yarrayt2Model = makeRealSpace(ut2Model,N)
#     xgrid,yarrayt4Model = makeRealSpace(ut4Model,N)
    
#     # initialize figure
#     myFig = plt.figure()
#     ax = plt.subplot(111)
#     ax.axis(xmin = 0,xmax = 2*np.pi-np.pi/N,ymin = -2, ymax = 3)
#     myFullLine, = ax.plot([],[], 'blue')
#     myMarkovLine, = ax.plot([],[], 'red')
#     mytModelLine, = ax.plot([],[], 'green')
#     myt2ModelLine, = ax.plot([],[], 'purple')
#     myt4ModelLine, = ax.plot([],[], 'black')

#     # define function to draw each frame
#     def makeFrame(n):
#         myFullLine.set_data(xgrid,yarrayFull[:,n])
#         myMarkovLine.set_data(xgrid,yarrayMarkov[:,n])
#         mytModelLine.set_data(xgrid,yarraytModel[:,n])
#         myt2ModelLine.set_data(xgrid,yarrayt2Model[:,n])
#         myt4ModelLine.set_data(xgrid,yarrayt4Model[:,n])
#         plt.title('t = '+str(round(t[n],1)))
#         plt.legend((myFullLine,myMarkovLine,mytModelLine,myt2ModelLine,myt4ModelLine),("Exact Solution","Markov Model","t-Model","2nd Order CMA","4th Order CMA"))
#         return myFullLine,myMarkovLine

#     # generate animation
#     anim = animation.FuncAnimation(fig = myFig,func = makeFrame,frames = t.shape[0])
    
#     return anim
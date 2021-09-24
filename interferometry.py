"""
interferometry.py is a Python routine that computes the response and
sensitivity functions of space-based GW detectors, e.g., LISA and Taiji,
to the detection of stochastic gravitational wave backgrounds.
"""

import astropy.constants as const
import astropy.units as u
import numpy as np

def read_sens(SNR=10, T=4, interf='LISA', Xi=False):

    """
    Function that reads the sensitivity Omega (expressed as a GW energy
    density spectrum), and the power law sensitivity (PLS)
    of the interferometer chosen (options are 'LISA', 'Taiji', and the
    LISA-Taiji network 'comb').

    Arguments:
        SNR -- signal-to-noise ratio (SNR) of the resulting PLS (default 10)
        T -- duration of the mission (in years) of the resulting PLS
             (default 4)
        interf -- option to chose the interferometer (default 'LISA')
        Xi -- option to return the helical sensitivity and PLS (default False)

    Returns:
        fs -- array of frequencies
        interf_Om -- sensitivity of the interferometer 'interf'
                     expressed as a GW energy density spectrum
        interf_OmPLS --  PLS of the interferometer 'interf'
        interf_Xi -- helical sensitivity of the interferometer 'interf'
        interf_XiPLS -- helical PLS of the interferometer 'interf'
    """

    import pandas as pd

    fact = SNR/np.sqrt(T)
    dir = 'detector_sensitivity'

    if interf=='LISA':

        fs, LISA_Om = read_csv(dir, 'LISA_Omega')
        fs, LISA_OmPLS = read_csv(dir, 'LISA_OmegaPLS')
        LISA_OmPLS *= fact
        if Xi:
            fs, LISA_Xi = read_csv(dir, 'LISA_Xi', b='Xi')
            fs, LISA_XiPLS = read_csv(dir, 'LISA_XiPLS', b='Xi')
            LISA_XiPLS *= fact
            return fs, LISA_Om, LISA_OmPLS, LISA_Xi, LISA_XiPLS
        else: return fs, LISA_Om, LISA_OmPLS

    if interf=='Taiji':

        fs, Taiji_Om = read_csv(dir, 'Taiji_Omega')
        fs, Taiji_OmPLS = read_csv(dir, 'Taiji_OmegaPLS')
        Taiji_OmPLS *= fact
        if Xi:
            fs, Taiji_Xi = read_csv(dir, 'Taiji_Xi', b='Xi')
            fs, Taiji_XiPLS = read_csv(dir, 'Taiji_XiPLS', b='Xi')
            Taiji_XiPLS *= fact
            return fs, Taiji_Om, Taiji_OmPLS, Taiji_Xi, Taiji_XiPLS
        else: return fs, Taiji_Om, Taiji_OmPLS

    if interf=='comb':

        fs, LISA_Taiji_Xi = read_csv(dir, 'LISA_Taiji_Xi', b='Xi')
        fs, LISA_Taiji_XiPLS = read_csv(dir, 'LISA_Taiji_XiPLS', b='Xi')
        LISA_Taiji_XiPLS *= fact
        return fs, LISA_Taiji_Xi, LISA_Taiji_XiPLS

def read_csv(dir, file, a='f', b='Omega'):

    """
    Function that reads a csv file with two arrays and returns them.

    Arguments:
        dir -- directory that contains the file
        file -- name of the csv file

    Returns:
        x -- first array of the file
        y -- second array of the file
        a -- identifier in pandas dataframe of first array (default 'f')
        b -- identifier in pandas dataframe of second array (default 'Omega')
    """

    import pandas as pd

    df = pd.read_csv(dir + '/' + file + '.csv')
    x = np.array(df[a])
    y = np.array(df[b])

    return x, y

# not used at the moment (replaced by read_csv)
def read_pic(dir, file):

    """
    Function that reads a pickle variable file and returns it.

    Arguments:
        dir -- directory that contains the file
        file -- name of the pickle file

    Returns:
        x -- pickle variable contained in the file
    """

    import pickle

    f = open(dir + '/' + file, 'rb')
    x = pickle.load(f)
    f.close()

    return x

def Poms_f(f, P=15, L=0):

    """
    Function that computes the power spectral density (PSD) of the optical
    metrology system noise.

    Arguments:
        f -- frequency array (should be in units of Hz)
        P -- noise parameter (default 15 for LISA)
             value for Taiji is P = 8
        L -- length of the interferometer arm (default is L = 2.5e6 km for LISA)

    Returns:
        Poms -- oms PSD noise
    """

    f_mHz = f.to(u.mHz)
    if L == 0: L = 2.5e6*u.km
    L_pm = L.to(u.pm)
    Poms = P**2/L_pm.value**2*(1+(2/f_mHz.value)**4)/u.Hz

    return Poms

def Pacc_f(f, A=3, L=0):

    """
    Function that computes the power spectral density (PSD) of the mass
    acceleration noise.

    Arguments:
        f -- frequency array (should be in units of Hz)
        A -- noise parameter (default 3 for LISA)
             value for Taiji is A = 3
        L -- length of the interferometer arm (default is L = 2.5e6 km for LISA)

    Returns:
        Pacc -- mass acceleration PSD noise
    """

    f_mHz = f.to(u.mHz)
    if L == 0: L = 2.5e6*u.km
    L_fm = L.to(u.fm)
    c = const.c
    Loc = L/c
    Loc = Loc.to(u.s)

    fsinv = (c/2/np.pi/f/L)
    fsinv = fsinv.to(1)

    Pacc = A**2*Loc.value**4/L_fm.value**2*(1 + (.4/f_mHz.value)**2)*(1 + \
           (f_mHz.value/8)**4)*fsinv.value**4/u.Hz

    return Pacc

def Pn_f(f, P=15, A=3, L=0):

    """
    Function that computes the noise power spectral density (PSD) of an
    interferometer channel X, Pn(f), and of the cross-correlation of two
    different interferometer channels XY, Pn_cross(f).

    Arguments:
        f -- frequency array (should be in units of Hz)
        P -- noise parameter (default 15 for LISA)
        A -- noise parameter (default 3 for LISA)
        L -- length of the interferometer arm (default is L = 2.5e6 km for LISA,
             otherwise should be in units of km)

    Returns:
        Pn -- noise PSD
        Pn_cross -- cross-correlation noise PSD
    """

    Poms = Poms_f(f, P=P, L=L)
    Pacc = Pacc_f(f, A=A, L=L)
    if L == 0: L = 2.5e6*u.km
    c = const.c
    f0 = c/2/np.pi/L
    f_f0 = f.to(u.Hz)/f0.to(u.Hz)

    Pn = Poms + (3 + np.cos(2*f_f0.value))*Pacc
    Pn_cross = -.5*np.cos(f_f0.value)*(Poms + 4*Pacc)

    return Pn, Pn_cross

def Pn_TDI(f, P=15, A=3, L=0):

    """
    Function that computes the noise power spectral density (PSD) of the
    interferometer channels A and T (Sagnac channel) of the TDI combination
    of channels.

    Arguments:
        f -- frequency array (should be in units of Hz)
        P -- noise parameter (default 15 for LISA)
        A -- noise parameter (default 3 for LISA)
        L -- length of the interferometer arm (default is L = 2.5e6 km for LISA,
             otherwise should be in units of km)

    Returns:
        PnA -- noise PSD of the TDI channel A
        PnT -- noise PSD of the TDI channel T
    """

    PnX, PnXY = Pn_f(f, P=P, A=A, L=L)
    PnA = 2*(PnX - PnXY)/3
    PnT = (PnX + 2*PnXY)/3

    return PnA, PnT

def R_f(f, L=0):

    """
    Function that computes the analytical fit of the response function.

    Arguments:
        f -- frequency array (should be in units of Hz)
        L -- length of the interferometer arm (default is L = 2.5e6 km for LISA,
             otherwise should be in units of km)

    Returns:
        Rf -- response function (using analytical fit)
    """

    if L == 0: L = 2.5e6*u.km
    c = const.c
    f0 = c/2/np.pi/L
    f_f0 = f.to(u.Hz)/f0.to(u.Hz)
    Rf = .3/(1 + .6*f_f0.value**2)

    return Rf

def Sn_f(f, P=15, A=3, L=0):

    """
    Function that computes the strain sensitivity using the analytical fit
    for an interferometer channel X.

    Arguments:
        f -- frequency array (should be in units of Hz)
        P -- noise parameter (default 15 for LISA)
        A -- noise parameter (default 3 for LISA)
        L -- length of the interferometer arm (default is L = 2.5e6 km for LISA,
             otherwise should be in units of km)

    Returns:
        Rf -- response function (using analytical fit)
    """

    Pn = Pn_f(f, P=P, A=A, L=L)
    Rf = R_f(f, L=L)

    return Pn/Rf

def delta(a, b):

    """
    Function that returns the Kronecker delta delta(a, b).
    Returns 1 if a = b and 0 otherwise.
    """

    if a==b: return 1
    else: return 0

def compute_interferometry(f, L=0):

    """
    Function that computes numerically the monopole and dipole response
    functions of an interferometer channel of a space-based GW detector.

    Arguments:
        f -- frequency array (should be in units of Hz)
        L -- length of the interferometer arm (default is L = 2.5e6 km for LISA,
             otherwise should be in units of km)

    Returns:
        MAA -- monopole response function of the TDI channel A
        MEE -- monopole response function of the TDI channel E (same as A)
        MAE -- monopole response function of the TDI correlation of the
               channels A and E
        MTT -- monopole response function of the TDI channel T (Sagnac channel)
        DAE -- dipole response function of the TDI correlation of the
               channels A and E
    """

    if L == 0: L = 2.5e6*u.km
    c = const.c

    # integration over sky directions (theta, phi)
    theta = np.linspace(0, np.pi, 50)
    phi = np.linspace(0, 2*np.pi, 50)
    # array of wave numbers
    k = 2*np.pi*f/c
    kL = L*k
    kL = kL.to(1)

    kLij, th, ph = np.meshgrid(kL, theta, phi, indexing='ij')

    kx1 = 0
    kx2 = np.cos(th)
    kx3 = .5*(np.sqrt(3)*np.cos(ph)*np.sin(th) + np.cos(th))

    kU1 = kx2 - kx1
    kU2 = kx3 - kx2
    kU3 = kx1 - kx3

    TkU1 = np.exp(-1j*kLij*(1+kU1)/2)*np.sinc(kLij.value*(1 - kU1)/2/np.pi)
    TkU1 += np.exp(1j*kLij*(1-kU1)/2)*np.sinc(kLij.value*(1 + kU1)/2/np.pi)
    TkU2 = np.exp(-1j*kLij*(1+kU2)/2)*np.sinc(kLij.value*(1 - kU2)/2/np.pi)
    TkU2 += np.exp(1j*kLij*(1-kU2)/2)*np.sinc(kLij.value*(1 + kU2)/2/np.pi)
    TkU3 = np.exp(-1j*kLij*(1+kU3)/2)*np.sinc(kLij.value*(1 - kU3)/2/np.pi)
    TkU3 += np.exp(1j*kLij*(1-kU3)/2)*np.sinc(kLij.value*(1 + kU3)/2/np.pi)
    TkmU1 = np.exp(-1j*kLij*(1-kU1)/2)*np.sinc(kLij.value*(1 + kU1)/2/np.pi)
    TkmU1 += np.exp(1j*kLij*(1+kU1)/2)*np.sinc(kLij.value*(1 - kU1)/2/np.pi)
    TkmU2 = np.exp(-1j*kLij*(1-kU2)/2)*np.sinc(kLij.value*(1 + kU2)/2/np.pi)
    TkmU2 += np.exp(1j*kLij*(1+kU2)/2)*np.sinc(kLij.value*(1 - kU2)/2/np.pi)
    TkmU3 = np.exp(-1j*kLij*(1-kU3)/2)*np.sinc(kLij.value*(1 + kU3)/2/np.pi)
    TkmU3 += np.exp(1j*kLij*(1+kU3)/2)*np.sinc(kLij.value*(1 - kU3)/2/np.pi)

    U1 = np.array([0, 0, 1])
    U2 = .5*np.array([np.sqrt(3), 0, -1])
    U3 = -.5*np.array([np.sqrt(3), 0, 1])

    c = np.matrix([[2, -1, -1], [0, -np.sqrt(3), np.sqrt(3)],
                  [1,1,1]])/3

    QA = np.zeros((3, 3, len(kL), len(theta), len(phi)))*1j
    QE = np.zeros((3, 3, len(kL), len(theta), len(phi)))*1j
    QT = np.zeros((3, 3, len(kL), len(theta), len(phi)))*1j
    for i in range(0, 3):
        for j in range(0,3):
            Q1 = .25*np.exp(-1j*kLij*kx1)*(TkU1*U1[i]*U1[j] - TkmU3*U3[i]*U3[j])
            Q2 = .25*np.exp(-1j*kLij*kx2)*(TkU2*U2[i]*U2[j] - TkmU1*U1[i]*U1[j])
            Q3 = .25*np.exp(-1j*kLij*kx3)*(TkU3*U3[i]*U3[j] - TkmU2*U2[i]*U2[j])
            QA[i,j,:,:,:] = Q1*c[0,0] + Q2*c[0,1] + Q3*c[0,2]
            QE[i,j,:,:,:] = Q1*c[1,0] + Q2*c[1,1] + Q3*c[1,2]
            QT[i,j,:,:,:] = Q1*c[2,0] + Q2*c[2,1] + Q3*c[2,2]

    k1 = np.cos(ph)*np.sin(th)
    k2 = np.sin(ph)*np.sin(th)
    k3 = np.cos(th)

    e1ab = np.zeros((3, 3, len(kL), len(theta), len(phi)))*1j
    for i in range(0, 3):
        if i==0: ki=k1
        elif i==1: ki=k2
        else: ki=k3
        for j in range(0, 3):
            if j==0: kj=k1
            elif j==1: kj=k2
            else: kj=k3
            e1ab[i,j,:,:,:] = delta(i,j) - ki*kj
            if i==0:
                if j==1: e1ab[i,j,:,:,:] += -1j*k3
                elif j==2: e1ab[i,j,:,:,:] += 1j*k2
            elif i==1:
                if j==0: e1ab[i,j,:,:,:] += 1j*k3
                elif j==2: e1ab[i,j,:,:,:] += -1j*k1
            else:
                if j==0: e1ab[i,j,:,:,:] += -1j*k2
                elif j==1: e1ab[i,j,:,:,:] += 1j*k1

    FAA = np.zeros((len(kL), len(theta), len(phi)))*1j
    FAE = np.zeros((len(kL), len(theta), len(phi)))*1j
    FEE = np.zeros((len(kL), len(theta), len(phi)))*1j
    FTT = np.zeros((len(kL), len(theta), len(phi)))*1j

    for a in range(0, 3):
        for b in range(0, 3):
            for c in range(0, 3):
                for d in range(0, 3):
                    eabcd = .25*e1ab[a,c,:,:,:]*e1ab[b,d,:,:,:]
                    FAA += eabcd*QA[a,b,:,:,:]*np.conjugate(QA[c,d,:,:,:])
                    FAE += eabcd*QA[a,b,:,:,:]*np.conjugate(QE[c,d,:,:,:])
                    FEE += eabcd*QE[a,b,:,:,:]*np.conjugate(QE[c,d,:,:,:])
                    FTT += eabcd*QT[a,b,:,:,:]*np.conjugate(QT[c,d,:,:,:])

    # Monopole and dipole response functions of LISA
    MAA1 = np.trapz(FAA*np.sin(th), th, axis=1)
    MAA = np.trapz(MAA1, phi, axis=1)/np.pi
    MAE1 = np.trapz(FAE*np.sin(th), th, axis=1)
    MAE = np.trapz(MAE1, phi, axis=1)/np.pi
    MEE1 = np.trapz(FEE*np.sin(th), th, axis=1)
    MEE = np.trapz(MEE1, phi, axis=1)/np.pi
    MTT1 = np.trapz(FTT*np.sin(th), th, axis=1)
    MTT = np.trapz(MTT1, phi, axis=1)/np.pi
    DAE1 = 1j*np.trapz(FAE*np.sin(th)**2*np.sin(ph), th, axis=1)
    DAE = np.trapz(DAE1, phi, axis=1)/np.pi

    return MAA, MAE, MEE, MTT, DAE

def refine_M(f, M, A=.3, exp=0):

    """
    Function that refines the response function by appending
    a A*f^exp in the low-frequency regime.

    Arguments:
        f -- frequency array (should be in units of Hz)
        M -- response function to be refined at low frequencies
        A -- amplitude of the response function at low frequencies
             (default 0.3 for the LISA monopole response function)
        exp -- exponent of the response function at low frequencies (default 0
               for the LISA monopole response function)

    Returns:
        fs -- refined array of frequencies
        Ms -- refined response function
    """

    ff0 = np.logspace(-6, np.log10(f[0].value), 1000)*u.Hz
    fs = np.append(ff0, f)
    Ms = np.append(A*ff0.value**exp, np.real(M))
    return fs, Ms

def compute_response_LISA_Taiji(save=True):

    import pandas as pd

    f = np.logspace(-4, 0, 5000)*u.Hz

    # LISA response functions
    L = 2.5e6*u.km
    MAA, MAE, MEE, MTT, DAE = compute_interferometry(f, L=L)
    # Taiji response functions
    L = 3e6*u.km
    MAA_Tai, MAE_Tai, MEE_Tai, MTT_Tai, DAE_Tai = \
            compute_interferometry(f, L=L)

    # refine response functions at low frequencies
    fs, MAs = refine_M(f, MAA)
    fs, MAs_Tai = refine_M(f, MAA_Tai)
    fs, MTs = refine_M(f, MTT, A=1.709840e6, exp=6)
    fs, MTs_Tai = refine_M(f, MTT_Tai, A=5.105546e6, exp=6)
    fs, MAEs = refine_M(f, MAE, A=-2.7870505059493287e-5)
    fs, MAEs_Tai = refine_M(f, MAE_Tai, A=-2.7870505059493287e-5)
    fs, DAEs = refine_M(f, DAE, A=.2)
    fs, DAEs_Tai = refine_M(f, DAE_Tai, A=.2)

    # Write response functions in csv files
    if save:
        df = pd.DataFrame({'frequency': fs, 'MA': MAs, 'MT': MTs,
                           'MAE': MAEs, 'DAE': DAEs})
        df.to_csv('detector_sensitivity/LISA_response_f.csv')
        df_Tai = pd.DataFrame({'frequency': fs, 'MA': MAs_Tai, 'MT': MTs_Tai,
                               'MAE': MAEs_Tai, 'DAE': DAEs_Tai})
        df_Tai.to_csv('detector_sensitivity/Taiji_response_f.csv')

def read_response_LISA_Taiji():

    import pandas as pd

    df = pd.read_csv('detector_sensitivity/LISA_response_f.csv')
    df_Tai = pd.read_csv('detector_sensitivity/Taiji_response_f.csv')
    fs = np.array(df['frequency'])
    MAs = np.array(df['MA'])
    MAEs = np.array(df['MAE'])
    MTs = np.array(df['MT'])
    DAEs = np.array(df['DAE'])
    MAs_Tai = np.array(df_Tai['MA'])
    MAEs_Tai = np.array(df_Tai['MAE'])
    MTs_Tai = np.array(df_Tai['MT'])
    DAEs_Tai = np.array(df_Tai['DAE'])
    fs = fs*u.Hz

    return fs, MAs, MAEs, MTs, DAEs, MAs_Tai, MAEs_Tai, MTs_Tai, DAEs_Tai

def read_MAC(M='MAC', V='V'):

    """
    Function that reads the V response functions of the cross-correlated
    channels of the LISA-Taiji network.

    Argument:
        M -- string of the channels to be read (options are default 'MAC',
             'MAD', 'MEC', 'MED')
        V -- can be changed to read the 'I' response function (default 'V')
    """

    import pandas as pd
    dir = 'detector_sensitivity/LISA_Taiji/'
    df = pd.read_csv(dir + M + '_' + V + '.csv')
    f = np.array(df['f'])
    MAC = np.array(df['M'])
    inds = np.argsort(f)
    f = f[inds]
    MAC = MAC[inds]
    f = f*u.Hz
    return f, MAC

def Oms(f, S):

    """
    Function that returns the sensitivity in terms of the GW energy density.

    Arguments:
        f -- frequency array (should be in units of Hz)
        S -- strain sensitivity function

    Returns:
        Omega -- GW energy density sensitivity
    """

    H0 = 100*u.km/u.s/u.Mpc
    H0 = H0.to(u.Hz)
    A = 8*np.pi**2/3/H0**2
    Omega = S*A*f**3
    return Omega

def OmPLS(f, Oms, beta, SNR=1, T=0, Xi=0):

    """
    Function that computes the power law sensitivity (PLS).

    Arguments:
        f -- frequency array (should be in units of Hz)
        Oms -- GW energy density sensitivity
        SNR -- threshold signal-to-noise ratio (SNR) to compute the PLS
               (default 1)
        T -- duration of the observation (in units of year, default 1)
        Xi -- allows to compute PLS for polarization signals using the dipole
              response function setting Xi = 0.25 (default 0)

    Returns:
        Omega -- GW energy density power law sensitivity (PLS)
    """

    Cbeta = np.zeros(len(beta))
    if T == 0: T = 1*u.yr
    T = T.to(u.s)
    for i in range(0, len(beta)):
        aux = f.value**(2*beta[i])
        aa = abs(1 - Xi*beta[i])
        Cbeta[i] = SNR*.5/aa/np.sqrt(np.trapz(aux/Oms**2, f.value)*T.value)

    funcs = np.zeros((len(f), len(beta)))
    for i in range(0, len(beta)): funcs[:,i] = f.value**beta[i]*Cbeta[i]
    Omega = np.zeros(len(f))
    for j in range(0, len(f)): Omega[j] = np.max(funcs[j,:])

    return Omega

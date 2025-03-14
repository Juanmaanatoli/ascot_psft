"""
Generate input from EQDSK files.

"""
import numpy as np
import numpy.matlib as matlib
import scipy.interpolate as interp
from ascot_psft.preprocessing import ReadEQDSK as readeqdsk
from a5py.ascot5io import B_2DS


def generate(fn_in, cocosin=0, fn_out=None, sigma_ip=0, sigma_b0=0, desc=None):
    """
    Function to create the ascot5 equilibrium field starting from an eqdsk
    
    Parameters
        fn_in    (str): filename of the eqdsk to read
        cocosin  (int): default=0. cocos of the eqdsk. If 0, it will ask you for it
        fn_out   (str): default=None. Name of ascot5 h5 file in output
        sigma_ip (int): default=0. Sign of current desired in output (compliant with COCOS)
        sigma_b0 (int): default=0. Sign of Bfield desired in output (compliant with COCOS)
	desc (str): default = None. Description. 
    Arguments
        B2D (dict): dict with the variables to be written to ascot5
    """

    _eqd = readeqdsk.ReadEQDSK(fn_in)
    if cocosin==0:
        cocosin = input('Please insert input cocos of eqdsk:')
        cocosin = int(cocosin)
    cocos_check(_eqd, cocosin)

    if cocosin==3:
        print("Lucky! Eqdsk already in cocos for ASCOT5! (number 3)")
        eqd=_eqd
    else:
        eqd = cocos_transform(_eqd, cocosin, sigma_ip, sigma_b0)

    B2D = {}    
    
    B2D["Rmin"] = eqd.rboxleft
    B2D["Rmax"] = eqd.rboxleft+eqd.rboxlength
    B2D["nR"]   = eqd.nrbox
    
    B2D["zmin"] = -eqd.zboxlength/2.
    B2D["zmax"] = eqd.zboxlength/2.
    B2D["nz"]   = eqd.nzbox
    
    B2D["axisR"] = eqd.Raxis
    B2D["axisz"] = eqd.Zaxis
    
    B2D["psiRz"] = eqd.psi
    B2D["psiaxis"] = eqd.psiaxis
    B2D["psisepx"] = eqd.psiedge
        
    B2D["B_R"] = eqd.psi*0
    B2D["B_z"] = eqd.psi*0
    
    # Toroidal component is more complicated.
    # It can be evaluated from Btor = F/R but
    # we need to map F(psi) to F(R,z) first.
    # However, F(psi) is given only inside the plasma.
    fspl = interp.interp1d(eqd.psi_grid, eqd.T, kind='linear', bounds_error = False, fill_value=eqd.T[-1])
    fpolrz = fspl(eqd.psi)
    Rvec = np.linspace(eqd.rboxleft, eqd.rboxleft+eqd.rboxlength, eqd.nrbox)
    zvec = np.linspace(-eqd.zboxlength/2., eqd.zboxlength/2., eqd.nzbox)
    Rgrid = matlib.repmat(Rvec,eqd.nzbox,1)  
    
    B2D["B_phi"] = fpolrz/Rgrid
    
    # Check psi value
    fpsi = interp.RectBivariateSpline(Rvec, zvec, eqd.psi.T)
    psimaxis = fpsi(eqd.Raxis,eqd.Zaxis)
    print(f'\npsiaxis is: {B2D["psiaxis"]}')
    print(f'Interpolation gives: {psimaxis[0]}')
    if np.abs(B2D["psiaxis"]-psimaxis[0])<0.01:
        print('These more or less agree!')
    else:
        print('Please check that these agree!')
    print(f'\npsisepx is: {B2D["psisepx"]}')
    
    if fn_out != None:
        B_2DS.write_hdf5(fn_out, B2D["Rmin"], B2D["Rmax"], B2D["nR"], B2D["zmin"], B2D["zmax"], B2D["nz"], \
                   B2D["axisR"], B2D["axisz"], B2D["psiRz"].T, B2D["psiaxis"], B2D["psisepx"], \
                   B2D["B_R"].T, B2D["B_phi"].T, B2D["B_z"].T, desc=desc)
    
    return B2D

def cocos_check(eqd, COCOS):
    """cocos check
    This function does not trust the user and check the COCOS of the input eqdsk
    
    Parameters:
        eqd :  eqdsk structure read with ReadEQDSK
        COCOS (int): input cocos.
    Attributes:
        None
    """
    cocos_dict = fill_cocosdict(COCOS)
    sign_q = np.sign(eqd.q[0])
    sign_ip = np.sign(eqd.Ip)
    sign_b0 = np.sign(eqd.B0EXP)
    if sign_q*cocos_dict['sigma_rhothetaphi']*sign_ip*sign_b0<0:
        print(f"The sign of q is not consistent. sign(q)={sign_q}, it should be sigma_rhothetaphi*sign_ip*sign_b0={cocos_dict['sigma_rhothetaphi']}*{sign_ip}*{sign_b0}={cocos_dict['sigma_rhothetaphi']*sign_ip*sign_b0}")
    if np.sign(np.mean(eqd.T))*sign_b0<0:
        print(f"The sign of F and B0 do not correspond")
    psi_increase = np.sign(eqd.psiedge-eqd.psiaxis) # 1 if increasing, -1 if decreasing

    sigma_bp = psi_increase/sign_ip
    sigma_rhothetaphi = sign_q/(sign_ip*sign_b0)

    # Finding if psi has 2pi value
    a_minor = max(eqd.R)-min(eqd.R)
    R0 = (max(eqd.R)+min(eqd.R))/2.
    deltastar_axis = 4.*(eqd.psi_grid[-1]-eqd.psi_grid[0])/a_minor**2
    mu0 = 4*np.pi*1e-7
    rjphi_nopi = - mu0*eqd.pprime[0]*R0**2 - eqd.TTprime[0]
    rjphi_4pi2 = rjphi_nopi* 4* np.pi**2
    cocos_gt_10= False if abs(rjphi_nopi-deltastar_axis) < abs(rjphi_4pi2-deltastar_axis) else True

    cocos_read = _sign_to_cocos(sigma_bp, sigma_rhothetaphi, cocos_gt_10)

    if COCOS not in cocos_read:
        error_msg = "=== \n"
        error_msg += f"You said cocos {COCOS}, I read cocos {cocos_read[0:2]} (depending on direction of phi)"
        error_msg += " \nWe strongly suggest to check the cocos you have.\nThis is fundamental for correct B field creation.\n"
        error_msg += "=== \n"
        print(error_msg)
        error_msg = f" COCOS in eqdsk {cocos_read} is different from the one in input [{COCOS}]"
        raise Exception(error_msg)

    print("Good! Your cocos matches!")

def _sign_to_cocos(sigma_bp, sigma_rhothetaphi, cocos_gt_10):
    """
    Associating the sign with the correct cocos
    
    Parameters:
        sigma_bp (int): sigma_bp as in cocos definition
        sigma_rhothetaphi (int): sigma_rhothetaphi as in cocos definition
        cocos_gt_10 (bool): true if the cocos should be >10
    Attributes:
        cocos_read (array): cocos inferred by input eqdsk
    """
    cocos_read=[0,0]
    if sigma_bp>0:
        if sigma_rhothetaphi>0:
            cocos_read = [1, 2, 11, 12]
        else:
            cocos_read = [5, 6, 15, 16]
    else:
        if sigma_rhothetaphi>0:
            cocos_read = [7, 8, 17, 18]
        else:
            cocos_read = [3, 4, 13, 14]
    cocos_read = np.array(cocos_read)
    #if cocos_gt_10:
    #    cocos_read = cocos_read[cocos_read>10]
    #else:
    #    cocos_read = cocos_read[cocos_read<10]
    return cocos_read

def cocos_transform(eqd, COCOS, sigma_ip=0, sigma_b0=0):
    """ cocos transformation
    This function converts the magnetic input from their starting cocos to cocos 3 (needed by ascot5)

    https://www.sciencedirect.com/science/article/pii/S0010465512002962

    Table 4 is not updated. EFIT has usually cocos=7, but if people try to keep q>0, then cocos changes.

    Parameters:
        COCOS (int): input cocos
        sigma_ip (int): output sign of ip requested. If 0, it will keep the ones in the eqdsk
        sigma_b0 (int): output sign of b0 requested. If 0, it will keep the ones in the eqdsk
        
    Attributes:
        None
    """
    print("COCOS tranformation from "+str(COCOS)+" to 3")
    eqdout = eqd

    cocosin = fill_cocosdict(COCOS)
    cocosin['sigma_ip'] = np.sign(eqd.Ip)
    cocosin['sigma_b0'] = np.sign(eqd.B0EXP)

    #These cocos are for COCOS 3
    cocosout = fill_cocosdict(3)
    #Checking the signs of the current and field desired as output
    cocosout['sigma_ip'] = np.sign(eqd.Ip)    if sigma_ip == 0 else sigma_ip
    cocosout['sigma_b0'] = np.sign(eqd.B0EXP) if sigma_b0 == 0 else sigma_b0

    # Define effective variables: sigma_Ip_eff, sigma_B0_eff, sigma_Bp_eff, exp_Bp_eff as in Appendix C
    sigma_Bp_eff = cocosin['sigma_Bp'] * cocosout['sigma_Bp']
    exp_Bp_eff   = cocosout['exp_Bp']  - cocosin['exp_Bp']
    sigma_rhothetaphi_eff = cocosin['sigma_rhothetaphi'] * cocosout['sigma_rhothetaphi']
    if 'sigma_ip' in cocosout.keys() and 'sigma_b0' in cocosout.keys():
        print(f"Requested sign(Ip)= {cocosout['sigma_ip']}, sign(b0)= {cocosout['sigma_b0']}")
        sigma_Ip_eff = cocosin['sigma_ip']*cocosout['sigma_ip']
        sigma_B0_eff = cocosin['sigma_b0']*cocosout['sigma_b0']
    else:
        print('No sign of Ip nor B0 requested')
        sigma_Ip_eff = cocosin['sigma_RphiZ']*cocosout['sigma_RphiZ']
        sigma_B0_eff = cocosin['sigma_RphiZ']*cocosout['sigma_RphiZ']

    # Define input
    F_in       = eqd.T
    FFprime_in = eqd.TTprime
    pprime_in  = eqd.pprime

    psirz_in   = eqd.psi
    psiaxis_in = eqd.psiaxis
    psiedge_in = eqd.psiedge
    psigrid_in = eqd.psi_grid

    q_in  = eqd.q
    b0_in = eqd.B0EXP
    ip_in = eqd.Ip
    
    # Transform
    F         = F_in       * sigma_B0_eff
    FFprime   = FFprime_in * sigma_Ip_eff * sigma_Bp_eff / (2*np.pi)**exp_Bp_eff
    pprime    = pprime_in  * sigma_Ip_eff * sigma_Bp_eff / (2*np.pi)**exp_Bp_eff
    
    _fact_psi = sigma_Ip_eff * sigma_Bp_eff * (2*np.pi)**exp_Bp_eff
    psirz     = psirz_in   * _fact_psi 
    psi_grid  = psigrid_in * _fact_psi
    psiaxis   = psiaxis_in * _fact_psi
    psiedge   = psiedge_in * _fact_psi
    
    q  = q_in  * sigma_Ip_eff * sigma_B0_eff * sigma_rhothetaphi_eff
    b0 = b0_in * sigma_B0_eff
    ip = ip_in * sigma_Ip_eff

    # Define output
    eqdout.T       = F
    eqdout.TTprime = FFprime
    eqdout.pprime  = pprime
    
    eqdout.psi      = psirz
    eqdout.psi_grid = psi_grid
    eqdout.psiaxis  = psiaxis
    eqdout.psiedge  = psiedge
    
    eqdout.q     = q
    eqdout.B0EXP = b0
    eqdout.Ip    = ip
    return eqdout

def fill_cocosdict(COCOS):
    """
    Function to fill the dictionary with the COCOS variables

    Parameters:
        COCOS (int): input cocos
    Args:
        cocosdict (dict): dictionary with cocos variables
    """
    cocos_keys = ['sigma_Bp', 'sigma_RphiZ', 'sigma_rhothetaphi',\
          'sign_q_pos', 'sign_pprime_pos', 'exp_Bp']
    cocosdict = dict.fromkeys(cocos_keys)

    cocosdict['exp_Bp'] = 1 if COCOS > 10 else 0 # if COCOS>=10, this should be 1

    if COCOS==1 or COCOS==11:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==2 or COCOS==12:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==3 or COCOS==13:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = +1
    elif COCOS==4 or COCOS==14:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = +1
    elif COCOS==5 or COCOS==15:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==6 or COCOS==16:
        cocosdict['sigma_Bp']          = +1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = -1
        cocosdict['sign_q_pos']        = -1
        cocosdict['sign_pprime_pos']   = -1
    elif COCOS==7 or COCOS==17:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = +1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = +1
    elif COCOS==8 or COCOS==18:
        cocosdict['sigma_Bp']          = -1
        cocosdict['sigma_RphiZ']       = -1
        cocosdict['sigma_rhothetaphi'] = +1
        cocosdict['sign_q_pos']        = +1
        cocosdict['sign_pprime_pos']   = +1
    else:
        raise ValueError(f"COCOS {COCOS} does not exists \n")

    return cocosdict

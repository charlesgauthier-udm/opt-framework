import numpy as np
import fenv
import time
from numba import njit
import lib_fixed_turbation.soilcprocessesf90.turbation

@njit
def resveg(litrplsz,soilplsz,envltr,envsol,reduceatdepth,bsratelt,bsratesc,fcancmx,ncell):
    ltresveg = np.zeros((20,ncell,11))
    scresveg = np.zeros((20,ncell,11))
    for i in range(0,20):
        for j in range(0,11):
            for k in range(0,ncell):
                if fcancmx[j,k] > 0:
                    ltresveg[i,k,j] = bsratelt[j] * litrplsz[j,i,k] * envltr[i,k] * reduceatdepth[i,k]
                    scresveg[i,k,j] = bsratesc[j] * soilplsz[j,i,k] * envsol[i,k] * reduceatdepth[i,k]
                else:
                    ltresveg[i,k,j] = 0.0
                    scresveg[i,k,j] = 0.0
    return ltresveg, scresveg

@njit
def update_pool(Cinput, ltresveg, humicfac,spinfast,scresveg,litrplsz,soilplsz,delt,ncell):

    dx1 = np.zeros((11,20,ncell))
    dx2 = np.zeros((11,20,ncell))
    for i in range(0,20):
        for j in range(0,11):
            for k in range(0,ncell):
                dx1[j,i,k] = Cinput[j,i,k] - (ltresveg[i,k,j] * (1+humicfac[j]))
                dx2[j,i,k] = spinfast * ((ltresveg[i,k,j] * humicfac[j]) - scresveg[i,k,j])
                litrplsz[j,i,k] += dx1[j,i,k] * delt
                soilplsz[j,i,k] += dx2[j,i,k] * delt
    return litrplsz,soilplsz

def model(param,ncell,iccp2,ignd,litrmass,soilcmas,tf,delt,ndays, thice, thpor, thliq, psisat, tbar, b,
          itermax,spinup,zbotw,actlyr,SAND,fcancmx,delzw,Cinput,rank):
    """
    Function that simulates
    :param param: PFT specific parameters t obe optimized
    :param litrmass: Mass of carbon in litter pool, per layer, per PFT [kg C m^-2]
    :param soilcmas: Mass of carbon in soil pool, per layer, per PFT [kg C m^-2
    :param ignd: Number of soil layers
    :param iccp2: Number of CTEM PFTs + bare + LUC
    :param itermax: Number of iterations in the simulation [days]
    :param tf: Length of the forcing data [yr]
    :param delt: Length of a timestep [yr]
    :param spinup: Number of spinup iterations
    :param ndays: Number of days of data
    :param thice: Frozen water content [m^3 m^-3]
    :param thpor: Total soil porority [m^3 m^-3]
    :param thliq: Liquid water content [m^3 m^-3]
    :param psisat: Soil moisture suction at saturation [m]
    :param isand: Flag for soil type -4 = ice, -3 = rock
    :param tbar: Temperature at each soil layer [K]
    :param b: Clapp and Hornberger empirical param []
    :param zbotw: Bottom of permeable soil layers [m]
    :param actlyr: Max annual active layer depth [m]
    :param SAND: Sand content of soils [%]
    :param delzw: Thickness of permeable part of soil layers [m]
    :param Cinput: Carbon input to litter pool, dim(ndays, iccp2, ignd, ncell) [kg C m^-2 yr-1]
    :param fcancmx: PFT coverage fraction of gridcells + bare + LUC
    :param ncell: Number of gridcell to simulate
    :return: litter and soil pool size per layer [kg C m^-2]
    """
    # Variables needed for run
    # Array of pft specific variables, 9 CTEM pfts + Bare + LUC

    apfts = np.array([[param['bsratelt_NdlEvgTr'], param['bsratelt_NdlDcdTr'], param['bsratelt_BdlEvgTr'],
                       param['bsratelt_BdlDCoTr'], param['bsratelt_BdlDDrTr'], 0.6,0.6, param['bsratelt_GrassC3'], param['bsratelt_GrassC4'],
                       0.5605, 0.5605],
         [param['bsratesc_NdlEvgTr'], param['bsratesc_NdlDcdTr'], param['bsratesc_BdlEvgTr'],
          param['bsratesc_BdlDCoTr'], param['bsratesc_BdlDDrTr'], 0.035, 0.035,
          param['bsratesc_GrassC3'], param['bsratesc_GrassC4'], 0.02258, 0.02258],
         [param['humicfac_NdlEvgTr'], param['humicfac_NdlDcdTr'], param['humicfac_BdlEvgTr'],
          param['humicfac_BdlDCoTr'], param['humicfac_BdlDDrTr'], 0.10, 0.10,
          param['humicfac_GrassC3'], param['humicfac_GrassC4'], 0.45, 0.45]])
    threshold = 0.001  # equilibrium threshold

    # Initializing soil C content vector
    litrplsz = litrmass[:,:iccp2,:].transpose(1,0,2)   # initial carbon mass of litter pool [kg C / m^2]
    soilplsz = soilcmas[:, :iccp2,:].transpose(1,0,2)  # initial carbon mass of soil pool [kg C / m^2]
    del(litrmass)
    del(soilcmas)

    # arrays for outputs
    litter = np.empty([ignd, ncell, 20440])
    soil = np.zeros([ignd, ncell, 20440])
    resp = np.zeros([ncell, 20440])

    count = 0  # iteration count
    output_count = 0  # number of output counts
    eqflag = 0  # Equilibrium flag, raised once equilibrium is reached
    totsoil_prev = 0
    totsoil_prev2 = np.zeros(7300)
    extra_loop = 0  # number of output counts

    # Assigning parameter values from algorithm to their corresponding variable
    bsratelt = apfts[0, :]  # turnover rate of litter pool [kg C/kg C.yr]
    bsratesc = apfts[1, :]  # turnover rate of soil pool [kg C/kg C.yr]
    humicfac = apfts[2, :]  # Humification factor
    kterm = 3.0 #param['kterm']  # Constant that determines the depth at which turbation stops
    biodiffus = param['biodiffus'] # Diffusivity coefficient for bioturbation [m^2/day]
    cryodiffus = param['cryodiffus'] # Diffusivity coefficient for cryoturbation [m^2/day]
    tanha = param['tanha'] # Constant a for tanh formulation of respiration Q10 determination
    tanhb = param['tanhb'] # Constant b for tanh formulation of respiration Q10 determination
    tanhc = param['tanhc'] # Constant c for tanh formulation of respiration Q10 determination
    tanhd = param['tanhd'] # Constant d for tanh formulation of respiration Q10 determination
    r_depthredu = param['r_depthredu']
    tcrit = param['t_crit']
    frozered = 0.1#param['frozered']
    mois_a = 10000 #param['moisa']
    mois_b = 6 #param['moisb']
    mois_c = param['moisc']

    reduceatdepth = np.exp(-zbotw / r_depthredu)

    # Computing environmental modifiers for the whole run
    envmodltr, envmodsol = fenv.fenv_numba(ndays, thice, thpor, thliq, psisat, ignd, SAND, tbar, b, ncell, tanha, tanhb,
                                     tanhc, tanhd, tcrit, frozered, mois_a, mois_b, mois_c)
    del(thice)
    del(thliq)
    del(thpor)
    del(tbar)
    del(b)
    del(psisat)
    sample = int(20 / delt)
    # Time loop
    #print('STARTING MODEL RUN')
    while count < itermax:
        if not eqflag:
            t = count % sample  # cycling on the length of imput files

        # Modifying transfert coeff if spinup=True
        if eqflag:
            spinfast = 1  # transfer coeff from litter pool to soil
        else:
            spinfast = 10

        # Computing respiration over time step for both pool
        ltresveg,scresveg = resveg(litrplsz,soilplsz,envmodltr[t],envmodsol[t],reduceatdepth,bsratelt, bsratesc, fcancmx,ncell)
        # Updating pools
        litrplsz, soilplsz = update_pool(Cinput[t],ltresveg,humicfac,spinfast,scresveg, litrplsz, soilplsz, delt, ncell)
        # Calling turbation subroutine
        litrplsz_fortran,soilplsz_fortran = lib_fixed_turbation.soilcprocessesf90.turbation.soilcprocesses(int(1), ncell, int(10), cryodiffus,biodiffus,
                                                    kterm, zbotw.transpose(), np.round(SAND).astype('int').transpose(),
                                                    actlyr[t], spinfast, int(0), litrplsz.transpose(2,0,1), soilplsz.transpose(2,0,1))
        litrplsz = litrplsz_fortran.transpose(1,2,0)
        soilplsz = soilplsz_fortran.transpose(1,2,0)

        if eqflag==False:  # Check for equilibrium every 20 years
            totsoil_prev2[t] = np.sum(soilplsz) # Summing on all gridcells for soil pool
        count += 1  # Iteration count
        # Equilibrium check
        if count % 7300 == 0 and eqflag==False:  # Check for equilibrium every 20 years
            totsol = np.nanmean(totsoil_prev2)  # Summing on all gridcells for soil pool
            diff = np.abs(totsoil_prev - totsol)
            if diff/totsol <= threshold:  # If threshold is crossed, raise equilibrium flag
                eqflag = 1
            totsoil_prev = totsol

        if eqflag:  # Once equilibrium flag is up, start outputting
            t = count % sample
            if extra_loop >= 7300:
                gridltr = np.sum(litrplsz.transpose(1, 0, 2) * fcancmx, axis=1)  # Litter pool size at grid level
                gridsol = np.sum(soilplsz.transpose(1, 0, 2) * fcancmx, axis=1)  # SoilC pool size at grid level
                gridrsp = np.sum(((ltresveg + scresveg) * fcancmx.transpose()),
                                 axis=(0, 2))  # Respiration at grid level
                t = output_count + 22265
                # Adding to output arrays
                litter[:, :, output_count] = gridltr
                soil[:, :, output_count] = gridsol
                resp[:, output_count] = gridrsp
                output_count += 1  # Updating number of outputs
            extra_loop +=1

        if output_count >= 20440:  # Once equilibrium is reached we output tf * 365 more iterations and stop the simul #TODO have to generalize that
            break
    return litter, soil, resp

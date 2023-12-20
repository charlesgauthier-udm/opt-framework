import numpy as np
from numba import jit,njit

@njit
def custom_clip(val):
    """
    Similar to np.clip(val,0.2,1.0), but can be use with numba
    :param val: input value
    :return: clipped value
    """
    if val > 1.0:
        return 1.0
    elif val < 0.2:
        return 0.2
    else:
        return val


def fenv(ndays, thice, thpor, thliq, psisat, ignd, SAND, tbar, b, ncell,tanha,tanhb,tanhc,tanhd,tcrit,frozered,mois_a,mois_b,mois_c):
    """
    Function that computes environmental modifiers i.e. f_15(Q_10) and f(psi) described in:
    https://cccma.gitlab.io/classic/namespaceheterotrophicrespirationmod.html, following:
    https://gitlab.com/jormelton/classic/-/blob/develop/src/heterotrophicRespirationMod.f90
    :param ndays: Number of days of data
    :param thice: Frozen water content [m^3 m^-3]
    :param thpor: Total soil porority [m^3 m^-3]
    :param thliq: Liquid water content [m^3 m^-3]
    :param psisat: Soil moisture suction at saturation [m]
    :param ignd: Number of soil layers
    :param isand: Flag for soil type, -4 = ice, -3 = rock
    :param tbar: Temperature at each soil layer [K]
    :param b: Clapp and Hornberger empirical param []
    :param ncell: Number of gridcells simulated at once
    :param tanha: Constant a for tanh formulation of respiration Q10 determination
    :param tanhb: Constant b for tanh formulation of respiration Q10 determination
    :param tanhc: Constant c for tanh formulation of respiration Q10 determination
    :param tanhd: Constant d for tanh formulation of respiration Q10 determination
    :param tcrit: # Temperature below which respiration is inhibited [C]
    :param frozered:
    :param mois_a:
    :param mois_b:
    :param mois_c:
    :return:
    """
    # Defining variables
    isand = np.where((SAND == -4) | (SAND == -3))       # Flag for soil type -4 = ice, -3 = rock
    psisat = np.tile(psisat, (ndays,1,1))
   # tcrit = -1.0  # Temperature below which respiration is inhibited [C]
    tfrez = 273.16  # Freezing point of water [K]
   # frozered = 0.1  # Factor to reduce resp. by for temp below tcrit
    # Find the matric potential
    ipore = np.where(thice > thpor)  # Where the whole pore space is ice
    psi = psisat * (thliq / (thpor - thice))**(-b)
    psi[ipore] = mois_a    # pore space is ice, suction is very high

    # Find moisture scalar
    ltrmoscl = np.zeros([ndays, ignd, ncell])  # Litter pool moisture scalar
    socmoscl = np.zeros([ndays, ignd, ncell])  # Soil pool moisture scalar
    # Different conditions on psi
    cvrdry = np.where(psi >= mois_a)    # Soil is very dry
    cdrmst = np.where((mois_b < psi) & (psi < mois_a)) # Soil is dry to moist
    coptim = np.where((mois_c <= psi) & (psi <= mois_b))  # Soil is optimal
    ctoowt = np.where((psisat < psi) & (psi < mois_c)) # Soil is too wet
    csatur = np.where(psi <= psisat) # Soil is saturated
    # Computing the moisture scalars
    ltrmoscl[cvrdry] = 0.2
    socmoscl[cvrdry] = 0.2

    ltrmoscl[cdrmst] = 1.0 - 0.8 * ((np.log10(psi[cdrmst]) - np.log10(mois_b)) / (np.log10(mois_a) - np.log10(mois_b)))
    socmoscl[cdrmst] = ltrmoscl[cdrmst]

    ltrmoscl[coptim] = 1
    socmoscl[coptim] = 1

    socmoscl[ctoowt] = 1.0 - 0.5 * ((np.log10(4.0) - np.log10(psi[ctoowt])) / (np.log10(4.0) - np.log10(psisat[ctoowt])))
    ltrmoscl[ctoowt] = socmoscl[ctoowt]
    if 0 in np.unique(ctoowt):
        ltrmoscl[0] = 1

    socmoscl[csatur] = 0.5
    ltrmoscl[csatur] = socmoscl[csatur]
    if 0 in np.unique(csatur):
        ltrmoscl[0] = 1
    for i in range(0,ndays):
        ltrmoscl[i][isand] = 0.2
        socmoscl[i][isand] = 0.2

    ltrmoscl.clip(0.2, 1.0, out=ltrmoscl)
    socmoscl.clip(0.2, 1.0, out=socmoscl)

    ff_p = np.array([ltrmoscl, socmoscl])

    # Computing q10 response function
    litrq10 = tanha + tanhb * np.tanh(tanhc * (tanhd - (tbar - tfrez)))
    solcq10 = litrq10

    q10funcLitr = litrq10**(0.1 * (tbar - tfrez - 15.0))
    q10funcSoilC= solcq10**(0.1 * (tbar - tfrez - 15.0))
    # Accounting for frozen soil
    ifroz = np.where((tbar - tfrez) <= tcrit)
    q10funcLitr[ifroz] *= frozered
    q10funcSoilC[ifroz] *= frozered
    ff_T = np.array([q10funcLitr,q10funcSoilC])
    envmod = ff_T * ff_p
    return envmod

@njit
def fenv_numba(ndays, thice, thpor, thliq, psisat, ignd, SAND, tbar, b, ncell,tanha,tanhb,tanhc,tanhd,tcrit,frozered,mois_a,mois_b,mois_c):
    """
    Function that computes environmental modifiers i.e. f_15(Q_10) and f(psi) described in:
    https://cccma.gitlab.io/classic/namespaceheterotrophicrespirationmod.html, following:
    https://gitlab.com/jormelton/classic/-/blob/develop/src/heterotrophicRespirationMod.f90
    :param ndays: Number of days of data
    :param thice: Frozen water content [m^3 m^-3]
    :param thpor: Total soil porority [m^3 m^-3]
    :param thliq: Liquid water content [m^3 m^-3]
    :param psisat: Soil moisture suction at saturation [m]
    :param ignd: Number of soil layers
    :param isand: Flag for soil type, -4 = ice, -3 = rock
    :param tbar: Temperature at each soil layer [K]
    :param b: Clapp and Hornberger empirical param []
    :param ncell: Number of gridcells simulated at once
    :param tanha: Constant a for tanh formulation of respiration Q10 determination
    :param tanhb: Constant b for tanh formulation of respiration Q10 determination
    :param tanhc: Constant c for tanh formulation of respiration Q10 determination
    :param tanhd: Constant d for tanh formulation of respiration Q10 determination
    :param tcrit: # Temperature below which respiration is inhibited [C]
    :param frozered:
    :param mois_a:
    :param mois_b:
    :param mois_c:
    :return:
    """
    tfrez = 273.16  # Freezing point of water [K]
    # Initializing arrays
    psi = np.zeros((ndays,ignd,ncell))
    ltrmoscl = np.zeros((ndays,ignd,ncell))
    socmoscl = np.zeros((ndays,ignd, ncell))
    q10funcLitr = np.zeros((ndays,ignd,ncell))
    q10funcSoilc = np.zeros((ndays,ignd,ncell))
    # First we find the matric potential
    for k in range(0,ndays):
        for j in range(0,ignd):
            for i in range(0,ncell):
                if (SAND[j,i] != - 3) & (SAND[j,i] != -4): # Soil isnt ice or bedrock
                    # No lower limit on psi
                    if thice[k,j,i] <= thpor[j,i]: # There is some remaining pore space
                        psi[k,j,i] = psisat[j,i] * (thliq[k,j,i] / (thpor[j,i] - thice[k,j,i])) ** -b[j,i]
                    else:
                        psi[k,j,i] = mois_a

                    # Now find the moisture scalar
                    if psi[k,j,i] >= mois_a: # Very dry
                        ltrmoscl[k,j,i] = 0.2
                        socmoscl[k,j,i] = 0.2
                    elif (psi[k,j,i] < mois_a) & (psi[k,j,i] > mois_b): # dry to moist
                        ltrmoscl[k,j,i] = 1.0 - 0.8 * ((np.log10(psi[k,j,i] - np.log10(mois_b))) /
                                                       (np.log10(mois_a) - np.log10(mois_b)))
                        socmoscl[k,j,i] = ltrmoscl[k,j,i]
                    elif (psi[k,j,i] <= mois_b) & (psi[k,j,i] >= mois_c): # Optimal
                        ltrmoscl[k,j,i] = 1.0
                        socmoscl[k,j,i] = 1.0
                    elif (psi[k,j,i] < mois_c) & (psi[k,j,i] > psisat[j,i]): # Too wet (for soil C)
                        socmoscl[k,j,i] = 1.0 - 0.5 * ((np.log10(mois_c) - np.log10(psi[k,j,i]))/
                                                       (np.log10(mois_c) - np.log10(psisat[j,i])))
                        if j==0: # Only the litter at the surface are not impeded by too moist conditions
                            ltrmoscl[k,j,i] = 1.0
                        else:
                            ltrmoscl[k,j,i] = socmoscl[k,j,i]

                    elif psi[k,j,i] <= psisat[j,i]: # Saturated conditions (impeded respiration)
                        socmoscl[k,j,i] = 0.5
                        if j == 0: # Only litter at surface is not impeded at all
                            ltrmoscl[k,j,i] = 1.0
                        else: # Otherwise assume same environment as soil C
                            ltrmoscl[k,j,i] = socmoscl[k,j,i]
                else: # Soil is bedrok or ice
                    socmoscl[k,j,i] = 0.2
                    ltrmoscl[k,j,i] = 0.2

                ltrmoscl[k,j,i] = custom_clip(ltrmoscl[k,j,i])
                socmoscl[k,j,i] = custom_clip(socmoscl[k,j,i])
                #ltrmoscl[k,j,i] = ltrmoscl[k,j,i].clip(0.2, 1.0)
                #socmoscl[k,j,i] = socmoscl[k,j,i].clip(0.2, 1.0)

                # Next we find the q10 response function
                litrq10 = tanha + tanhb * (np.tanh(tanhc * (tanhd - (tbar[k,j,i] - tfrez))))

                # Soil has the same values
                soilcq10 = litrq10

                # Apply reduction when soil layer freezes
                if (tbar[k,j,i] - tfrez) > tcrit:
                    q10funcLitr[k,j,i] = litrq10**(0.1 * (tbar[k,j,i] - tfrez - 15.0))
                    q10funcSoilc[k,j,i]= soilcq10**(0.1* (tbar[k,j,i] - tfrez - 15.0))
                else: # frozen/freezing soils
                    q10funcLitr[k, j, i] = litrq10 ** (0.1 * (tbar[k, j, i] - tfrez - 15.0)) * frozered
                    q10funcSoilc[k, j, i] = soilcq10 ** (0.1 * (tbar[k, j, i] - tfrez - 15.0)) * frozered
    envmod_litr = ltrmoscl * q10funcLitr
    envmod_soil = socmoscl * q10funcSoilc


    return envmod_litr, envmod_soil

def fenv2(ndays, thice, thpor, thliq, psisat, ignd,isand,tbar,b,nruns,tanha,tanhb,tanhc,tanhd,tcrit,frozered,mois_a,mois_b,mois_c):
    """
    Computes temperature and moisture modifier, outputs modifier per pool per pft
    :param i: current iteration
    :return: environmental modifier dim npool x npfts
    """
    # Defining variables
    psisat = np.tile(psisat, (ndays,1,1))
    tfrez = 273.16  # Freezing point of water [K]
    # Find the matric potential
    ipore = np.where(thice > thpor)  # Where the whole pore space is ice
    psi = psisat * (thliq / (thpor - thice))**(-b)
    psi[ipore] = 10000.0    # pore space is ice, suction is very high
    psi = np.tile(psi,(nruns))
    psisat = np.tile(psisat,(nruns))
    # Find moisture scalar
    ltrmoscl = np.zeros([ndays, ignd, nruns])  # Litter pool moisture scalar
    socmoscl = np.zeros([ndays, ignd, nruns])  # Soil pool moisture scalar
    # Different conditions on psi
    cvrdry = np.where(psi >= mois_a)    # Soil is very dry
    cdrmst = np.where((mois_b < psi) & (psi < mois_a)) # Soil is dry to moist
    coptim = np.where((mois_c <= psi) & (psi <= mois_b))  # Soil is optimal
    ctoowt = np.where((psisat < psi) & (psi < mois_c)) # Soil is too wet
    csatur = np.where(psi <= psisat) # Soil is saturated
    # Computing the moisture scalars
    ltrmoscl[cvrdry] = 0.2
    socmoscl[cvrdry] = 0.2

    ltrmoscl[cdrmst] = 1.0 - 0.8 * ((np.log10(psi[cdrmst]) - np.log10(mois_b[cdrmst[2]])) / (np.log10(mois_a[cdrmst[2]]) - np.log10(mois_b[cdrmst[2]])))
    socmoscl[cdrmst] = ltrmoscl[cdrmst]

    ltrmoscl[coptim] = 1
    socmoscl[coptim] = 1

    socmoscl[ctoowt] = 1.0 - 0.5 * ((np.log10(mois_c[ctoowt[2]]) - np.log10(psi[ctoowt])) / (np.log10(mois_c[ctoowt[2]]) - np.log10(psisat[ctoowt])))
    ltrmoscl[ctoowt] = socmoscl[ctoowt]
    if 0 in np.unique(ctoowt):
        ltrmoscl[0] = 1

    socmoscl[csatur] = 0.5
    ltrmoscl[csatur] = socmoscl[csatur]
    if 0 in np.unique(csatur):
        ltrmoscl[0] = 1
    for i in range(0,ndays):
        ltrmoscl[i][isand] = 0.2
        socmoscl[i][isand] = 0.2

    ltrmoscl.clip(0.2, 1.0, out=ltrmoscl)
    socmoscl.clip(0.2, 1.0, out=socmoscl)

    ff_p = np.array([ltrmoscl, socmoscl])

    # Computing q10 response function
    litrq10 = tanha + tanhb * np.tanh(tanhc * (tanhd - (tbar - tfrez)))
    solcq10 = litrq10

    q10funcLitr = litrq10**(0.1 * (tbar - tfrez - 15.0))
    q10funcSoilC = solcq10**(0.1 * (tbar - tfrez - 15.0))
    # Accounting for frozen soil
    atbar = (tbar-tfrez) - tcrit
    ifroz = np.where(atbar <= 0.0)
    temp = np.zeros_like(atbar)
    temp[ifroz] = True
    temp *= frozered
    temp[np.where(temp==0)] = 1.

    q10funcLitr *= temp
    q10funcSoilC *= temp
    ff_T = np.array([q10funcLitr,q10funcSoilC])
    envmod = ff_T * ff_p
    return envmod

def fenv3(ndays, thice, thpor, thliq, psisat, ignd,isand,tbar,b,nruns,tanha,tanhb,tanhc,tanhd,tcrit,frozered,mois_a,mois_b,mois_c):
    """
    Computes temperature and moisture modifier, outputs modifier per pool per pft
    :param i: current iteration
    :return: environmental modifier dim npool x npfts
    """
    # Defining variables
    psisat = np.tile(psisat, (ndays,1,1))
    tfrez = 273.16  # Freezing point of water [K]
    # Find the matric potential
    ipore = np.where(thice > thpor)  # Where the whole pore space is ice
    psi = psisat * (thliq / (thpor - thice))**(-b)
    psi[ipore] = 10000.0    # pore space is ice, suction is very high
    psi = np.tile(psi,(nruns))
    psisat = np.tile(psisat,(nruns))
    # Find moisture scalar
    ltrmoscl = np.zeros([ndays, ignd, nruns])  # Litter pool moisture scalar
    socmoscl = np.zeros([ndays, ignd, nruns])  # Soil pool moisture scalar
    # Different conditions on psi
    cvrdry = np.where(psi >= mois_a)    # Soil is very dry
    cdrmst = np.where((mois_b < psi) & (psi < mois_a)) # Soil is dry to moist
    coptim = np.where((mois_c <= psi) & (psi <= mois_b))  # Soil is optimal
    ctoowt = np.where((psisat < psi) & (psi < mois_c)) # Soil is too wet
    csatur = np.where(psi <= psisat) # Soil is saturated
    # Computing the moisture scalars
    ltrmoscl[cvrdry] = 0.2
    socmoscl[cvrdry] = 0.2

    ltrmoscl[cdrmst] = 1.0 - 0.8 * ((np.log10(psi[cdrmst]) - np.log10(mois_b)) / (np.log10(mois_a) - np.log10(mois_b)))
    socmoscl[cdrmst] = ltrmoscl[cdrmst]

    ltrmoscl[coptim] = 1
    socmoscl[coptim] = 1

    socmoscl[ctoowt] = 1.0 - 0.5 * ((np.log10(mois_c) - np.log10(psi)) / (np.log10(mois_c) - np.log10(psisat)))
    ltrmoscl[ctoowt] = socmoscl[ctoowt]
    if 0 in np.unique(ctoowt):
        ltrmoscl[0] = 1

    socmoscl[csatur] = 0.5
    ltrmoscl[csatur] = socmoscl[csatur]
    if 0 in np.unique(csatur):
        ltrmoscl[0] = 1
    for i in range(0,ndays):
        ltrmoscl[i][isand] = 0.2
        socmoscl[i][isand] = 0.2

    ltrmoscl.clip(0.2, 1.0, out=ltrmoscl)
    socmoscl.clip(0.2, 1.0, out=socmoscl)

    ff_p = np.array([ltrmoscl, socmoscl])

    # Computing q10 response function
    litrq10 = tanha + tanhb * np.tanh(tanhc * (tanhd - (tbar - tfrez)))
    solcq10 = litrq10

    q10funcLitr = litrq10**(0.1 * (tbar - tfrez - 15.0))
    q10funcSoilC = solcq10**(0.1 * (tbar - tfrez - 15.0))
    # Accounting for frozen soil
    atbar = (tbar-tfrez) - tcrit
    ifroz = np.where(atbar <= 0.0)
    temp = np.zeros_like(atbar)
    temp[ifroz] = True
    temp *= frozered
    temp[np.where(temp==0)] = 1.

    q10funcLitr *= temp
    q10funcSoilC *= temp
    ff_T = np.array([q10funcLitr,q10funcSoilC])
    envmod = ff_T * ff_p
    return envmod


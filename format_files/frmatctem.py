# rmatctem function script
import numpy as np
from numba import jit

@jit
def unfrozenroot(maxannualactlyr, zbotw, rmatctem):
    iccp2 = 11
    ignd = 20
    ndays = rmatctem.shape[0]
    ncell = rmatctem.shape[3]
    unfrozenroot = np.zeros((ndays, iccp2,ignd,ncell))
    botlyr = 1
    for i in range(0,ncell):
        for h in range(0, ndays):
            for k in range(0,ignd):
                if maxannualactlyr[h,i] < zbotw[k,i]:
                    break
                botlyr = k

            for j in range(0, iccp2):
                if botlyr == ignd:
                    unfrozenroot[h,j,:,i] = rmatctem[h,j,:,i]
                else:
                    frznrtlit = np.sum(rmatctem[h,j,botlyr + 0:ignd,i])
                    for k in range(0, botlyr):
                        unfrozenroot[h,j,k,i] = rmatctem[h,j,k,i] + rmatctem[h,j,k,i] / (1.0 - frznrtlit) *frznrtlit

    return unfrozenroot


@jit
def rmatctemf2(zbotw, alpha, rootmass, avertmas, abar, soildpth, maxannualactlyr, mxrtdpth):
    """
    Function that computes rthe rmatctem variable using numba on-the-fly compilator
    :param zbotw: Bottom of soil layers [m] dim(ignd)
    :param alpha: Parameter determining how roots grow dim(icc)
    :param rootmass: Root mass [kg C m^-2] dim(ndays,icc)
    :param avertmas: Average root biomass [kg C m^-2] dim(icc)
    :param abar: Parameter determining average root profile dim(icc)
    :param soildpth: Soil permeable depth [m] dim(scalar)
    :param maxannualactlyr: Active layer thickness maximum over the e-folding period specified by parameter eftime [m] dim(scalar)
    :param mxrtdpth: mxrtdpth: Maximum rooting depth [m] dim(icc)
    :return: rmatctem, fraction of live roots in each soil layer for each pft dim(icc, ignd)
    """

    iccp2 = alpha.shape[0]
    ignd = zbotw.shape[0]
    ndays = rootmass.shape[0]
    ncell = zbotw.shape[1]
    abszero = 1e-12
    # estimate parameter b of variable root profile parameterization
    b = abar * (avertmas ** alpha)
    # Estimate 99% rooting depth
    useb = b
    usealpha = alpha

    a = np.zeros((ndays,iccp2,ncell))
    totala = np.zeros((ndays,iccp2,ncell))
    rmatctem = np.zeros((ndays,iccp2, ignd,ncell))
    etmp = np.zeros((ndays, iccp2, ignd,ncell))
    tab = np.zeros((ndays,iccp2,ncell))
    for j in range(0,ncell):
        rootdpth = (4.605 * (rootmass[:,:,j] ** alpha)) / b
        for i in range(0,iccp2):
            for k in range(0,ndays):
                if rootdpth[k,i] > min(soildpth[j],maxannualactlyr[k,j],zbotw[ignd-1,j],mxrtdpth[i]):
                    rootdpth[k,i] = min(soildpth[j],maxannualactlyr[k,j],zbotw[ignd-1,j],mxrtdpth[i])
                    if rootdpth[k,i] <= abszero:
                        a[k,i,j] = 100.0
                    else:
                        a[k,i,j] = 4.605 / rootdpth[k,i]
                else:
                    if rootmass[k,i,j] <= abszero:
                        a[k,i,j] = 100.0
                    else:
                        a[k,i,j] = useb[i] / (rootmass[k,i,j] ** usealpha[i])
        for i in range(0,iccp2):
            for k in range(0,ndays):

                kend = 9999 # Initialize at dummy value

                # Using parameter 'a' we can find fraction of roots in each soil layer
                zroot = rootdpth[k,i]
                totala[k,i,j] = 1.0 - np.exp(-a[k,i,j] * zroot)

                # If rootdepth is shallower than the bottom of the first layer
                if zroot <= zbotw[0,j]:
                    rmatctem[k,i,0,j] = 1.0
                    rmatctem[k,i,1:,j] = 0.0
                    kend = 0
                else:
                    for tempk in range(1,ignd):
                        if (zroot <= zbotw[tempk,j]) & (zroot > zbotw[tempk-1,j]):
                            kend = tempk
                    if kend == 9999:
                        print('ERROR KEND IS NOT ASSIGNED')
                        kend = 1
                    etmp[k,i,0,j] = np.exp(-a[k,i,j] * zbotw[0,j])
                    rmatctem[k,i,0,j] = (1.0 - etmp[k,i,0,j]) / totala[k,i,j]
                    if kend == 1:
                        etmp[k,i,kend,j] = np.exp(-a[k,i,j] * zroot)
                        rmatctem[k,i,kend,j] = (etmp[k,i,kend-1,j] - etmp[k,i,kend,j]) / totala[k,i,j]
                    elif kend > 1:
                        for tempk in range(1,kend):
                            etmp[k,i,tempk,j] = np.exp(-a[k,i,j] * zbotw[tempk,j])
                            rmatctem[k,i,tempk,j] = (etmp[k,i,tempk-1,j] - etmp[k,i,tempk,j]) / totala[k,i,j]

                        etmp[k,i,kend,j] = np.exp(-a[k,i,j] * zroot)
                        rmatctem[k,i,kend,j] = (etmp[k,i,kend-1,j] - etmp[k,i,kend,j]) / totala[k,i,j]

                tab[k,i,j] = kend
    rmatctem = unfrozenroot(maxannualactlyr, zbotw, rmatctem)
    return rmatctem
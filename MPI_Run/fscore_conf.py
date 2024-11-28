import numpy as np
import pickle
import xarray as xr
from scipy.interpolate import interp1d
import numpy as np


def custom_score_conf(true_csoil, true_resp, sim_csoil, sim_resp, typicall_std):
    # Average csoil pool to get mean state
    av_sim_csoil = np.nanmean(sim_csoil,axis=2) # Time average
    av_sim_resp = np.nanmean(sim_resp,axis=1)

    av_true_csoil = np.nanmean(true_csoil,axis=0) # Time average
    av_true_resp = np.nanmean(true_resp,axis=0)

    # Normalized arrays
    norm_true_csoil = av_true_csoil / av_true_csoil
    norm_sim_csoil = av_sim_csoil / av_true_csoil

    norm_true_resp = av_true_resp / av_true_resp
    norm_sim_resp = av_sim_resp / av_true_resp

    #======================================================
    # SOIL C SCORES
    #======================================================
    # Error oriented score for soil C
    eo_csoil = typicall_std * (norm_sim_csoil - norm_true_csoil) ** 4
    eo_csoil = np.clip(eo_csoil, 0, 100) # Limit upper bound to 100
    eo_csoil = np.nanmean(eo_csoil) # Average over gridcells and soil layers

    # Mean oriented score for soil C
    mo_csoil = 1 - np.exp(-np.abs(norm_sim_csoil - norm_true_csoil))
    mo_csoil = np.nanmean(mo_csoil) # Average over gridcells and soil layers

    #======================================================
    # RESP SCORES
    #======================================================
    # Error oriented score for resp
    eo_resp = 1.0 * (norm_sim_resp - norm_true_resp) ** 4
    eo_resp = np.clip(eo_resp, 0, 100)
    eo_resp = np.nanmean(eo_resp)

    # Mean oriented score for resp
    mo_resp = 1 - np.exp(-np.abs(norm_sim_resp - norm_true_resp))
    mo_resp = np.nanmean(mo_resp)
    return eo_csoil, mo_csoil, eo_resp, mo_resp


def custom_score_wosis(wosis, classic, classic_zbot, typicall_std):
    """
    Function that compares CLASSIC's simulated carbon content with the WOSIS dataset and computes different scores based on the comparison
    :param wosis: Wosis dataset carbon content formated on the classic grid
    :param classic: simulated carbon content from CLASSIC
    :param classic_zbot: soil layer thickness of CLASSIC
    :param typicall_std: Typical normalized standard deviation from the wosis dataset
    :return: rmse: Root mean square error
             eo_score: Error oriented score
             mo_score: Mean oriented score
    """
    # Empty arrays to store formatted data
    av_tab = []  # Array of average value from Wosis at a certain depth
    model_tab = []#np.empty(0)  # Model value
    classic_orc_tot = np.mean(classic,axis=2)
    for h in range(0,len(wosis)):  # Looping on gridcell containing Wosis data
        # Observation data
        depth = np.array(wosis[h][1][0])  # List of all layer depths in the gridcell
        orgc = np.array(wosis[h][0][0])  # Array of wosis soil C value for all layers in the gridcell

        # Model output
        classic_orgc = classic_orc_tot[:,h]  # Model simulated soil C
        temp = np.zeros(21)  # Formatting to add layer at surface with 0 soil C
        temp[1:] = classic_orgc
        classic_orgc = temp
        values = np.unique(depth)  # All possible values of depth in the gridcell

        for k in range(0,len(values)):  # Looping on all possible depth value in the grid cell
            current_depth = values[k]  # Current depth that is being used
            ind = np.where(depth == values[k])  # Index of all layers that are at the current depth

            av = np.mean(orgc[ind])  # Average value of the layers at current depth

            f = interp1d(classic_zbot,classic_orgc)  # Interpolation function from model's layer depth and soil C values
            classic_inter_orgc = f(current_depth)  # Interpolating CLASSIC soil C at current depth

            #classic_inter_orgc = interpolate(classic_zbot,classic_orgc,current_depth)

            # Adding to empty arrays
            av_tab.append(av)
            model_tab.append(classic_inter_orgc)

    av_tab = np.array(av_tab)
    model_tab = np.array(model_tab)
    rmse = np.sum((model_tab - av_tab) ** 2) / len(av_tab)  # Computing root mean square error
    tab_eo_score = []  # Empty array for error oriented score
    tab_mo_score = []  # Empty array for mean oriented score
    for i in range(0,len(av_tab)):  # Looping on comparable values
        current_wosis = av_tab[i]  # Wosis value
        current_classic = model_tab[i]  # Classic value
        ncurrent_wosis = current_wosis / current_wosis  # Normalized wosis value
        ncurrent_classic = current_classic / current_wosis  # Normalized CLASSIC value
        current_eo_score = typicall_std * (ncurrent_classic - ncurrent_wosis) ** 4  # Error oriented score
        tab_eo_score.append(current_eo_score)  # Adding to array
        current_mo_score = 1 - np.exp(-np.abs(ncurrent_classic - ncurrent_wosis))  # Mean oriented score
        tab_mo_score.append(current_mo_score)  # Adding to array

    ind_eo = np.where(np.array(tab_eo_score) < 100)
    eo_corr = np.array(tab_eo_score)[ind_eo]
    eo_score = np.nansum(np.array(eo_corr)) / len(eo_corr)  # Total error oriented score
    mo_score = np.nansum(np.array(tab_mo_score)) / len(tab_mo_score)  # Total mean oriente score
    return eo_score, mo_score


def custom_score_srdb(srdb,classic,rmrveg):
    """
    Function that compares
    :param srdb: Annual respiration values from the SRDB V5 dataset
    :param classic: Simulated respiration values from the CLASSIC run
    :return: rmse: root mean square error
             eo_score: error oriented score
             mo_score: mean oriented score
    """
    srdb_resp = srdb[:, 2]  # Soil respiration data from SRDB
    srdb_begin_yr = srdb[:, 4]  # Start year of observation data
    srdb_stop_yr = srdb[:, 5]  # Stop year

    # Creating empty arrays to store data for comparison
    classic_resp = np.array([])  # Simulated soil respiration
    observd_resp = np.array([])  # Observational soil respiration data

    # Getting every datapoint for comparison
    for i in range(0, len(srdb)):  # Looping on each gridcell containing SRDB data
        # Converting start/stop year to index of the model output
        ibegin = (
            np.floor((365 * srdb_begin_yr[i][0][:, 0] - 1900 * 365 - 22265 + 1)).astype(int))  # array idex for begin yr

        istop = (np.floor(365 * srdb_stop_yr[i][0][:, 0] - 1900 * 365 - 22265 + 1)).astype(
            int)  # array index for stop yr

        for j in range(0, len(ibegin)):  # Looping on all data in a single gridcell
            # Adding datapoints to empty arrays
            current_resp = np.mean(classic[i, ibegin[j]:istop[j]]) # CLASSIC's simulated heterotrophic resp
            istopyear = istop[j] // 365     # Indice of the year of the last measurment
            auto_resp = rmrveg[istopyear,i] # Autotrophic respiration corresponding to this year
            current_resp += auto_resp       # Soil respiration (Auto + hetero)
            classic_resp = np.append(classic_resp, current_resp)

            observd_resp = np.append(observd_resp, srdb_resp[i][0][j, 0])

    tab_eo_score = [] # Empty array for error oriented score
    tab_mo_score = [] # Empty array for mean oriented score
    for i in range(0, len(observd_resp)): # Looping on all comparable values
        current_srdb = observd_resp[i]      # Srdb value
        current_classic = classic_resp[i]   # Classic value
        ncurrent_wosis = current_srdb / current_srdb        # Normalized srdb value
        ncurrent_classic = current_classic / current_srdb   # Normalized CLASSIC value
        current_eo_score = 1 * (ncurrent_classic - ncurrent_wosis) ** 4  # Error oriented score with typicall std of 1
        tab_eo_score.append(current_eo_score) # Adding to empty array
        current_mo_score = 1 - np.exp(-np.abs(ncurrent_classic - ncurrent_wosis))  # Mean oriented score
        tab_mo_score.append(current_mo_score) # Adding to empty array
        dummy = 1
    ind = np.where(np.array(tab_eo_score) < 100)
    eo_corr = np.array(tab_eo_score)[ind]
    eo_score = np.nansum(np.array(eo_corr)) / len(eo_corr)  # Total error oriented score
    mo_score = np.nansum(np.array(tab_mo_score)) / len(tab_mo_score) # Total mean oriented score
    return eo_score, mo_score

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
NAME:
Created By: Logan Wright
Created On: Dec 19, 2018
Contact: logan.wright(at)colorado.edu

Code was adopted from original MATLAB Code

Description of the atmospheric correction methods can be found in:
[1] L. A. Wright, “Hyperspectral Observations for Atmospheric Remote Sensing:
        Instrumentation, Atmospheric Correction, and Spectral Unmixing,”
        University of Colorado Boulder, 2018.
[2] L. A. Wright, B. C. Kindel, P. Pilewskie, N. Leisso, T. U. Kampe, and
        K. S. Schmidt, “Below-Cloud Atmospheric Correction of Airborne
        Hyperspectral Imagery Using Simultaneous Solar Spectral Irradiance
        Observations,” IN REVIEW, 2018.


DESCRIPTION:



'''
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import os.path

import modtran_tools
from load_nis import load_nis
from super_resample import super_resample

def standard_correction(obs,L0,I0,trans,mu):
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:

    Outputs:

    '''

    transmittance = super_resample((trans['Ts'] + trans['ts']) * (trans['T'] + trans['t']),trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
    spherical_albedo = super_resample(trans['sph'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])

    n_spectra =obs['spectra'].shape
    n_wvl = obs['resp_func']['wvl'].shape
    R = np.zeros((n_spectra[0],n_wvl[0]))

    for n in range(n_spectra[0]):
#        spectrum = super_resample(obs['spectra'][n,:],obs['resp_func']['wvl0'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
        spectrum = obs['spectra'][n,0:152]
        R[n,:] = np.pi * (spectrum - L0) / (spherical_albedo * np.pi * (spectrum - L0) + mu * I0 * ( transmittance ) )

    return R

def irrad_correction(obs,L0,If_dn,trans,mu = None):
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:

    Outputs:

    '''

    # This transmittance calculation is more correct!
    transmittance = super_resample((trans['Ts'] + trans['ts']) * (trans['T'] + trans['t']),trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])

    # Mirror Matlab Setup
    t_t = super_resample(trans['t'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
    t_Tso = super_resample(trans['Tso'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
    transmittance = t_t**2 + 2*t_t*np.sqrt(t_Tso) + t_Tso

    spherical_albedo = super_resample(trans['sph'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])

    n_spectra =obs['spectra'].shape
    n_wvl = obs['resp_func']['wvl'].shape
    R = np.zeros((n_spectra[0],n_wvl[0]))

    # CHECK IF Number of Spectra == Number of irradiance
    if obs['spectra'].shape[0] != If_dn.shape[0]:
        print('WARNING: Non-matching number of irradiance spectra')

    if mu is not None:
        for n in range(n_spectra[0]):
            # NEEDS TO BE ADJUSTED
#            spectrum = super_resample(obs['spectra'][n,:],obs['resp_func']['wvl0'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
            spectrum = obs['spectra'][n,0:152]
            R[n,:] = np.pi * (spectrum - L0[n,:]) / (spherical_albedo * np.pi * (spectrum - L0[n,:]) + If_dn[n,:] * ( transmittance ) )
    else:
        for n in range(n_spectra[0]):
#            spectrum = super_resample(obs['spectra'][n,:],obs['resp_func']['wvl0'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
            spectrum = obs['spectra'][n,0:152]
            R[n,:] = np.pi * (spectrum - L0[n,:]) / (spherical_albedo * np.pi * (spectrum - L0[n,:]) + If_dn[n,:] * ( transmittance ) )

    return R

def albedo(obs,Idn):
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:

    Outputs:

    '''
    n_spectra =obs['spectra'].shape
    n_wvl = obs['resp_func']['wvl'].shape
    A = np.zeros((n_spectra[0],n_wvl[0]))

    for n in range(n_spectra[0]):
#        spectrum = super_resample(obs['spectra'][n,:],obs['resp_func']['wvl0'],obs['resp_func']['wvl'],obs['resp_func']['fwhm'])
        spectrum = obs['spectra'][n,0:152]
        A[n,:] = spectrum * np.pi / Idn[n,:]

    return A

def adjacency_correction(obs,Idn,Iup,L0,Iup0,trans,mu = None, type = 'Standard'):
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:

    Outputs:

    '''
    n_spectra =obs['spectra'].shape
    n_wvl = obs['resp_func']['wvl'].shape
    R = np.zeros((n_spectra[0],n_wvl[0]))

    if type == 'Standard':
        if mu is None:
            print('ERROR! mu value required for standard adjacency correction!')
            return None
        # Correct Transmittance Calculations
#        transmittance = super_resample((trans['Ts'] + trans['ts']) * (trans['T'] + trans['t']),trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
        spherical_albedo = super_resample(trans['sph'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
#
#        transmittance1 = super_resample((trans['ts'] + trans['Ts']) * trans['t'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
#        transmittance2 = super_resample((trans['ts'] + trans['Ts']) * trans['T'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])

#        # Mirror Matlab Setup
        t_t = super_resample(trans['t'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
        t_Tso = super_resample(trans['Tso'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
        transmittance = t_t**2 + 2*t_t*np.sqrt(t_Tso) + t_Tso

        transmittance1 = super_resample((trans['ts'] + trans['Ts']) * trans['t'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
        transmittance2 = super_resample((trans['ts'] + trans['Ts']) * trans['T'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])

        for n in range(n_spectra[0]):
#            spectrum = super_resample(obs['spectra'][n,:],obs['resp_func']['wvl0'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
           spectrum = obs['spectra'][n,0:152]
           Rh_bar = (Iup[n,:] - Iup0) / ((mu * Idn * transmittance) + (Iup[n,:] - Iup0) * spherical_albedo)
           R[n,:] = (spectrum - L0 - ((Idn * transmittance1) / np.pi) * (Rh_bar / (1 - Rh_bar * spherical_albedo))) * ((np.pi * (1 - Rh_bar * spherical_albedo)) / (Idn * transmittance2))

    elif type == 'Irradiance':
        # Correct Transmittance Calculations
#        transmittance = super_resample(trans['ts'] + trans['Ts'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
        spherical_albedo = super_resample(trans['sph'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
#
#        transmittance1 = super_resample((trans['ts'] + trans['Ts']) * trans['t'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
#        transmittance2 = super_resample((trans['ts'] + trans['Ts']) * trans['T'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])

        # Mirror Matlab Setup
        t_t = super_resample(trans['t'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
        t_Tso = super_resample(trans['Tso'],trans['wvl'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
        transmittance = (t_t + np.sqrt(t_Tso))**2

        transmittance1 = (t_t + np.sqrt(t_Tso))*t_t
        transmittance2 = (t_t + np.sqrt(t_Tso))*np.sqrt(t_Tso)

        for n in range(n_spectra[0]):
#            spectrum = super_resample(obs['spectra'][n,:],obs['resp_func']['wvl0'], obs['resp_func']['wvl'], obs['resp_func']['fwhm'])
#            R[n,:] = np.pi * (spectrum - L0) / (spherical_albedo * np.pi * (spectrum - L0) + Idn * ( transmittance ) )

            spectrum = obs['spectra'][n,0:152]
            Rh_bar = (Iup[n,:] - Iup0[n,:]) / ((Idn[n,:] * transmittance) + (Iup[n,:] - Iup0[n,:]) * spherical_albedo)
            R[n,:] = (spectrum - L0[n,:] - ((Idn[n,:] * transmittance1) / np.pi) * (Rh_bar / (1 - Rh_bar * spherical_albedo))) * ((np.pi * (1 - Rh_bar * spherical_albedo)) / (Idn[n,:] * transmittance2))

    else:
        print('Invalid Type: Type "',type,'" is not a valid selection')

    return R

def mu_correction():
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:

    Outputs:

    '''
    # Mu Correction
    R = pi * (Rob - rad0) / (sph_alb_1km * pi *(Rob - rad0) + Ifdn_obs_trim * (t_1km**2 + (1 + mu) * t_1km * sqrt(Tso_1km) + mu*Tso_1km))

    rho_bar_ret_mu = (Ifup_obs_trim - Ifup0) / ((Ifdn_obs_trim *(t_1km**2 + (1 + mu) * t_1km * sqrt(Tso_1km) + mu*Tso_1km) + (Ifup_obs_trim - Ifup0) * sph_alb_1km))
    R = (Rob - rad0 - (Ifdn_obs_trim * (t_1km + sqrt(Tso_1km)) * t_1km) / pi * (rho_bar_ret_mu / (1 - rho_bar_ret_mu * sph_alb_1km))) * (pi / (Ifdn_obs_trim * (t_1km + sqrt(Tso_1km)) * mu * sqrt(Tso_1km))) * (1 - rho_bar_ret_mu * sph_alb_1km)

#    disp('NoAdj Mu:');disp(sqrt(sum(((mean(veg_ret_noadj_mu,1)-veg_r)*100)**2)/numel(veg_r)))
#    disp('Ret Mu:');disp(sqrt(sum(((mean(veg_ret_mu,1)-veg_r)*100)**2)/numel(veg_r)))


def load_ASD(filepath):
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:

    Outputs:

    '''

    print('NOT YET IMPLEMENTED')
    return

def plot_results(*args, title = 'Reflectances', yrng = (0,1), save = False):
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:

    Outputs:

    '''

    fig1 = plt.figure(figsize = (4,4))
    ax1 = plt.subplot()
    plt.xlabel('Wavelength [nm]',fontsize = 12)
    plt.ylabel('Reflectance',fontsize = 12)
    ax1.axis([350,1025,0,0.1])
    ax1.set_xticks([400,600,800,1000])
    plt.tight_layout(pad=2.0)

    for arg in args:
        if np.size(arg['spectra'].shape) == 2:
            mean_spec = np.mean(arg['spectra'],axis = 0)
        else:
            mean_spec = arg['spectra']
        ax1.plot(arg['wvl'],mean_spec,color = arg['color'],LineWidth = 1.5)
        plt.ylim(yrng)

    if save == True:
        for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
            item.set_fontsize(20)
        fig1.savefig( title + '.eps',dpi = 300, format = 'eps')

def calc_rmse(baseline, wvl, *args):
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:
        RMSE should be in format of list (name (string),RMSE, RMSE, RMSE, RMSE, ...)

    Outputs:

    '''
    out = dict()

#    wv_bands = np.logical_xor(wvl < 928, wvl > 989)
#    wv_bands = np.logical_xor(wv_bands,wvl > 1095)
#    wv_bands = np.where(wv_bands == True)

    for arg in args:
        RMSE = np.sqrt(np.sum(((np.mean(arg['spectra'], axis = 0) - baseline['spectra'])*100)**2)/len(baseline['spectra']))
#        RMSE_nowv = np.sqrt(np.sum(((np.mean(arg['spectra'][:,wv_bands], axis = 0) - baseline['spectra'][wv_bands])*100)**2)/len(baseline['spectra'][wv_bands]))

        out[arg['name']] = RMSE

    return out

def sampling_issue(I0_orig,SunEllipticFactor,zwvl,ssim_zfwhm,neon_wvl,neon_fwhm,t_Ts,t_ts):
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:

    Outputs:

    '''
    # Demonstrate Sampling Issue
    I0_ssim = SunEllipticFactor*super_resample(I0_orig[:,1],I0_orig[:,0],zwvl,ssim_zfwhm)
    Ig_ssim = I0_ssim*super_resample(t_Ts+t_ts,wvl,zwvl,ssim_zfwhm)
    Ig_NIS = I0*super_resample(t_Ts+t_ts,wvl,neon_wvl,neon_fwhm)
    Ig_ssim_resampled = super_resample(I0_ssim[17:end]*super_resample(t_Ts+t_ts,wvl,zwvl[17:end],ssim_zfwhm[17:end]),zwvl[17:end],neon_wvl,neon_fwhm)
    plt.figure()
    # plot(I0_orig(:,1),I0_orig(:,2),'k')
    plt.plot(neon_wvl,Ig_NIS,'Color',[0,0,0])
    plt.plot(zwvl,Ig_ssim,'Color',[33,102,172]/255)
    plt.plot(neon_wvl,Ig_ssim_resampled,'Color',[178,24,43]/255)

    plt.legend('NIS Sampling','SSIM Sampling','SSIM Resampled to NIS')
    plt.ylabel('Spectral Irradiance [W m^{-2} nm^{-1}]')
    plt.xlabel('Wavelength [nm]')
    plt.title('Modelled Downwelling Irradiance')
    plt.set(gca,'FontSize',18)
    plt.axis([720,800,0.7,1.3])
    plt.grid()
    plt.savefig('IrradianceDiffRespFunc_O2.svg','svg')
    plt.axis([350,1025,0,1.8]);
    plt.savefig('IrradianceDiffRespFunc_Full.svg','svg')

    plt.figure;
    plt.plot(neon_wvl,Ig_NIS-Ig_ssim_resampled,'k')
    plt.ylabel('Spectral Irradiance [W m^{-2} nm^{-1}]');
    plt.xlabel('Wavelength [nm]')
    plt.title('Difference Between NIS and SSIM Sampling')
    plt.axis([720,800,-0.15,0.1])
    plt.grid()
    plt.set(gca,'FontSize',18)
    plt.savefig('DiffRespFunc_O2.svg','svg')

def empirical_wvl_fit():
    '''
    Created by: Logan Wright
    Created On: December 20 2018

    Translated from original MATLAB Code

    DESCRIPTION:

    Inputs:

    Outputs:

    '''
    import numpy.matlib
    # save('Emp_Fit_SSIM_Jun16_Tarp48.mat','neon_wvl','tarp48_ret','tarp48_ret_TOA')

    # Test "Empirical Smoothing"
    # "Best" fit NIS offset is 1.28
    dlambda_NIS = np.arange(-1.3,-1.2,0.01);#1.28;#1.2:0.01:1.4;
    neon_wvl2 = np.matlib.repmat(neon_wvl,1,dlambda_NIS.size) + np.matlib.repmat(dlambda_NIS,152,1)
    dlambda_SSIM_zen = np.arange(-3,3,0.1)
    zwvl2 = np.matlib.repmat(zwvl,1,dlambda_SSIM_zen.size) + np.matlib.repmat(dlambda_SSIM_zen,439,1)
    for x in range(dlambda_NIS.size):
        rad0 = super_resample(r0,flx_wvl,neon_wvl2[:,x],neon_fwhm);
        T_wvl = super_resample(t_T_1km,wvl,neon_wvl2[:,x],neon_fwhm);
        t_wvl = super_resample(t_t_1km,wvl,neon_wvl2[:,x],neon_fwhm);
        Ts_wvl = super_resample(t_Ts_1km,wvl,neon_wvl2[:,x],neon_fwhm);
        ts_wvl = super_resample(t_ts_1km,wvl,neon_wvl2[:,x],neon_fwhm);
        Tso_wvl = super_resample(t_Tso_1km,wvl,neon_wvl2[:,x],neon_fwhm);
        sph_alb_layer = super_resample(t_sph_alb_1km,wvl,neon_wvl2[:,x],neon_fwhm);
        Tsts_Tt = super_resample((t_Ts + t_ts)*(t_T+t_t),wvl,neon_wvl2[:,x],neon_fwhm);

        tarp48_ret_TOA_WVLS[:,x] = pi*(Rob-rad0)/(sph_alb*pi*(Rob-rad0)+mu*I0*(Tsts_Tt));

        for y in range(dlambda_SSIM_zen.size):
            Ifdn_obs_trim = super_resample(zspect[ssirtime_ind,:].T,zwvl2[:,y],neon_wvl2[:,x],neon_fwhm)
            Ifup_obs_trim = super_resample(nspect[ssirtime_ind,:].T,nwvl,neon_wvl2[:,x],neon_fwhm)
            Ifup0 = super_resample(rho_a_I,flx_wvl,neon_wvl2[:,x],neon_fwhm)*Ifdn_obs_trim

            tarp48_ret_noadj_WVLS[:,x,y] = pi * (Rob - rad0) / (sph_alb_layer * pi * ...
                (Rob - rad0) + Ifdn_obs_trim * (T_wvl + t_wvl)**2)
            tarp48_ret_noadj_WVLS2[:,x,y] = pi * (Rob - rad0) / (sph_alb_layer * pi * ...
                (Rob - rad0) + Ifdn_obs_trim * (t_wvl**2 + 2 * t_wvl * sqrt(Tso_wvl) + Tso_wvl))




    dlambda_SSIM_zen = np.arange(-5,3,0.1)
    zwvl2 = np.matlib.repmat(zwvl,1,length(dlambda_SSIM_zen)) + np.matlib.repmat(dlambda_SSIM_zen,439,1)
    for x in range(dlambda_SSIM_zen.size):
        Ifdn_obs_trim = super_resample(zspect[ssirtime_ind,:].T,zwvl2[:,y],neon_wvl,neon_fwhm)
        tarp48_ret_TEST[:,x] = pi*(Rob-rad0)/(sph_alb_layer*pi*(Rob-rad0) + Ifdn_obs_trim*(Tt2_d[1:152]))

    dlambda_SSIM_nad = np.arange(-3,3,0.1)
    # dfwhm = -3:0.1:3
    neon_wvl2 = np.matlib.repmat(neon_wvl,1,dlambda_NIS.size) + np.matlib.repmat(dlambda_NIS,152,1)
    nwvl2 = np.matlib.repmat(nwvl,1,dlambda_SSIM_nad.size) + np.matlib.repmat(dlambda_SSIM_nad,256,1)
    zwvl2 = np.matlib.repmat(zwvl,1,dlambda_SSIM_zen.size) + np.matlib.repmat(dlambda_SSIM_zen,439,1)
    neon_fwhm2 = np.matlib.repmat(neon_fwhm,1,dfwhm.size) + np.matlib.repmat(dfwhm,152,1)
    Rob = np.mean(tarp48_spc[:,0:152],axis = 1).T
    for x in range(dlambda_SSIM_zen.size):
    #     rad_resample = super_resample(data_rad(:,2)/100,data_rad(:,1),neon_wvl2(:,x),neon_fwhm)
    #     sph_alb_resample(:,1) = super_resample(sph_alb_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm)
    #     sph_alb_resample(:,2) = super_resample(sph_alb_orig(:,2),wvl,neon_wvl2(:,x),neon_fwhm)
    #     Ts_resample = super_resample(Ts_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm)
    #     ts_resample = super_resample(ts_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm)
    #     T_resample = super_resample(T_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm)
    #     t_resample = super_resample(t_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm)
    #     sph_alb_layer_resample = (sph_alb_resample(:,1)-sph_alb_resample(:,2)*(t_resample(:,1)+T_resample(:,1)))/(1-sph_alb_resample(:,2)*(t_resample(:,1)+T_resample(:,1)))
        for y in range(dlambda_SSIM_nad.size):
            Ifdn_obs_trim_NIS_SSIM = super_resample(Ifdn_obs.T,zwvl2[:,x],neon_wvl,neon_fwhm)
            Ifup_obs_trim_NIS_SSIM = super_resample(Ifup_obs.T,nwvl2[:,y],neon_wvl,neon_fwhm)
            Ifup0_resample = rho_a_I*Ifdn_obs_trim_NIS_SSIM
            rad0_resample = rho_a_R*Ifdn_obs_trim_NIS_SSIM
            rho_bar_ret_NIS_SSIM = (Ifup_obs_trim_NIS_SSIM-Ifup0_resample)/((Ifdn_obs_trim_NIS_SSIM*(T[:,1]+t[:,1])**2)+(Ifup_obs_trim_NIS_SSIM-Ifup0_resample)*sph_alb_layer)

            tarp48_ret_SSIM_TEST[:,x,y] = (Rob - rad0_resample - (Ifdn_obs_trim_NIS_SSIM*(T[:,1]+t[:,1])*t[:,1])/pi*(rho_bar_ret_NIS_SSIM/(1-rho_bar_ret_NIS_SSIM*sph_alb_layer)))*(pi/(Ifdn_obs_trim_NIS_SSIM*(T[:,1]+t[:,1])*T[:,1]))*(1-rho_bar_ret_NIS_SSIM*sph_alb_layer)

if __name__ == '__main__':
    # Below Cloud Surface Reflectance Retrieval

    # Set Date
    date = '20150608'

    # This section includes filenames and constants for each day/flight line
    if date == '20150608':
        NISfile = os.path.abspath('../ATMO_CORR_Code/NIS01_20150608_165842_rdn')
        tarp3_coord = np.array(((303,307),(363,369)))
        tarp48_coord = np.array(((312,319),(364,370)))
        veg_coord = np.array(((323,328),(379,420)))
        EWroad_coord = np.array(((309,348),(353,356)))
        NSroad_coord = np.array(((333,336),(316,354)))
        filename = [os.path.abspath('../ATMO_CORR_Code/NEON_20150608_Baseline'),
            os.path.abspath('../ATMO_CORR_Code/NEON_20150608_1kmAtmo_40H2O')]
        SunEllipticFactor = 0.97083    # Determined by Day of Year
        mu = np.cos(np.deg2rad(31))                 # Determined by Time of Day
    #     SSIRfile = '/Users/wrightad/Documents/Data/NEON/FlightData/20150608_SSIR.mat'
    #     SSIRfile = '/Users/wrightad/Documents/Data/NEON/FlightData/20150608_SSIR_CALIB_0602.mat'
        SSIRfile = os.path.abspath('../ATMO_CORR_Code/20150608_SSIR_CALIBSPECT_Jul2017.mat')
        asdfile = os.path.abspath('../ATMO_CORR_Code/20150608_ASD_REFL.mat')

    elif date == '20150616':
        NISfile = os.path.abspath('../ATMO_CORR_Code/NIS01_20150616_151613_rdn')
        tarp3_coord = np.array(((139,145),(231,237)))
        tarp48_coord = np.array(((149,155)(231,237)))
        veg_coord = np.array(((154,160),(252,280)))
        EWroad_coord = np.array(((147,185),(219,224)))
        NSroad_coord = np.array(((167,172),(181,218)))
        filename = (os.path.abspath('../ATMO_CORR_Code/NEON_20150616_Baseline'),
            os.path.abspath('../ATMO_CORR_Code/NEON_20150616_259mAtmo_30H2O'))
        SunEllipticFactor = 0.96950    # Determined by Day of Year
        mu = np.cos(np.deg2rad(50.2))             # Determined by Time of Day
        SSIRfile = os.path.abspath('../ATMO_CORR_Code/20150616_SSIR_CALIBSPECT_20150616FIELDB.mat')
        asdfile = os.path.abspath('../ATMO_CORR_Code/20150616_ASD_REFL.mat')
    #     asdfile = '/Users/wrightad/Documents/Data/NEON/GroundData/ASD_Refl/20150608_ASD_REFL.mat'

    elif date == '20150617':
        NISfile = os.path.abspath('../ATMO_CORR_Code/NIS01_20150617_203142_rdn')
        tarp3_coord = np.array(((173,180),(199,205)))
        tarp48_coord = np.array(((182,190),(199,207)))
        veg_coord = np.array(((195,199),(219,245)))
        EWroad_coord = np.array(((172,226),(189,192)))
        NSroad_coord = np.array(((207,211),(158,190)))
        filename = (os.path.abspath('../ATMO_CORR_Code/NEON_20150617_Baseline'),
            os.path.abspath('../ATMO_CORR_Code/NEON_20150617_1kmAtmo_40H2O'))
        SunEllipticFactor = 0.96933    # Determined by Day of Year
        mu = np.cos(np.deg2rad(25.0))                # Determined by Time of Day
        SSIRfile = os.path.abspath('../ATMO_CORR_Code/20150617_SSIR.mat')
        asdfile = os.path.abspath('../ATMO_CORR_Code/20150617_ASD_REFL.mat')

    N = len(filename)

    conv = 10000   # Scale Factor to convert [W cm^-2 nm^-1] to [W m^-2 nm^-2], used for .flx MODTRAN output file

    # Load NIS Data
    nis_datacube = load_nis(NISfile)
#    metadata = load_nis(NISfile+'_obs_ort')
    nistime = nis_datacube.nav['gps_time']

    # Modify NEON Instrument Response Function
    neon_wvl0 = nis_datacube.resp_func['wvl'] + 1.28     # Add Empirically Determined Offset
    neon_fwhm0 = nis_datacube.resp_func['fwhm']

    # Load NEON Instrument Response Function # Try Other Response Function
    respfunc = np.genfromtxt('../ATMO_CORR_Code/nis2013_fmed2_2013.dat')
    neon_wvl0 = respfunc[:,1]+1.28
    neon_fwhm0 = respfunc[:,2]

    # Set Wavelength Ranges
    # SSIR Index (25) to Neon Wvl (364) - Zenith
    # SSIR Index (25) to NEON Wvl (152) - Nadir
    neon_wvl_zen = neon_wvl0[0:364]
    neon_fwhm_zen = neon_fwhm0[0:364]
    neon_wvl_nad = neon_wvl0[0:152]
    neon_fwhm_nad = neon_fwhm0[0:152]
    # neon_wvl_nad = neon_wvl0(1:125)
    #neon_fwhm_nad = neon_fwhm0(1:125)
    # This truncates at 1000 nm (near Si/InGaAs joining)
    neon_wvl = neon_wvl_nad
    neon_fwhm = neon_fwhm_nad
    # L = len(neon_wvl0)

    # Load MODTRAN Irradiance
    data = np.genfromtxt(os.path.abspath('../ATMO_CORR_Code/SUN01med2irradwnNormt.dat'),dtype = float, skip_header = 2)
    data = data[1:50000,:]
    #     Convert to wavelength increments
    data[:,0] = 1e7/data[:,0]
    data[:,1] = data[:,1]*(1/data[:,0]**2)*1e11
    I0 = data[data[:,0]<3000,:]

    I0_orig = I0

    #     Convolve to NEON Wavelength
    I0 = SunEllipticFactor*super_resample(I0[:,1],I0[:,0],neon_wvl,neon_fwhm)

    # Load in MODTRAN Output
    baseline_acd = modtran_tools.load_acd(filename[0]+'.acd')

    # Load MODTRAN Irradiance
    baseline_flx = modtran_tools.load_flx(filename[0]+'.flx')
    flight_ind = 10; # Index of flight altitude (nominally 2.7 km on 6/8 and 6/17)

    # flx_fup = conv*super_resample(flx_data.data(:,8),flx_data.data(:,1),neon_wvl,neon_fwhm)
    # flx_fdn = flx_data.data[:,33]+flx_data.data[:,34]
    # flx_fup = flx_data.data[:,32]
    Iup0 = np.flip(baseline_flx['upwelling'][flight_ind])*10000 # factor of 10,000 to convert from W*cm-2*nm-1 to W*m-2*nm-1

    # Load MODTRAN Radiance
    baseline_7sc = modtran_tools.load_7sc(filename[0]+'.7sc')
    wvl_7sc = baseline_7sc['wvl']
    r0 = baseline_7sc['path_radiance']/100 # factor of 100 to convert from uW*cm-2*sr-1*nm-1 to W*m2*sr-1*nm-1
#    Iup0 = baseline_7sc['path_irradiance']/100 # factor of 100 to convert from uW*cm-2*nm-1 to W*m2*nm-1

#    r0 = r0/100

    # r0 = r0[0:,:]/100
    # Iup = Iup[0,:]

    # Calculate Atmospheric Reflectance

    rho_a_I = np.flip(baseline_flx['upwelling'][flight_ind]/(baseline_flx['downwelling'][flight_ind]))
    rho_a_R = r0/(1e4*np.flipud(baseline_flx['downwelling'][flight_ind]))

    rad0_neon = super_resample(r0,wvl_7sc,neon_wvl,neon_fwhm)
    Iup0_neon = super_resample(Iup0,wvl_7sc,neon_wvl,neon_fwhm)

    path_reflectances = dict([('wvl',wvl_7sc),('rad_path_refl',rho_a_R),('irrad_path_refl',rho_a_I)])

    flight_acd = modtran_tools.load_acd(filename[1] + '.acd')

    # Load SSIR Data
    ssim = sio.loadmat(SSIRfile)
    zwvl = ssim['zwvl']-2.5
    nwvl = ssim['nwvl']

    # ssim_zwvl = zwvl
    zfwhm0 = np.array(((365.25, 7.67),
                       (404.84, 8.11),
                       (435.83, 7.25),
                       (546.07, 7.88),
                       (912.30, 10.86),
                       (965.78, 15.97),
                       (1244.30, 12.70),
                       (1694.06, 9.75),
                       (1791.46, 10.79),
                       (2061.62, 17.58)))
    ssim['zfwhm'] = np.zeros((439,1))
    ssim['zfwhm'][0:192] = np.interp(zwvl[0:192],zfwhm0[0:5,0],zfwhm0[0:5,1])#,'linear','extrap')
    ssim['zfwhm'][192:] = np.interp(zwvl[192:],zfwhm0[5:,0],zfwhm0[5:,1])#,'linear','extrap')
    # ssim_nwvl = nwvl
    nfwhm0 = np.array(((365.25, 7.50),
                            (404.84, 7.85),
                            (435.83, 7.28),
                            (546.07, 8.16),
                            (912.30, 10.75)))
    ssim_nfwhm = np.interp(nwvl,nfwhm0[:,0],nfwhm0[:,1])#'linear','extrap')

    ssim['wvl'] = zwvl[24:225]



    # Load Ground Reflectances
    asd_20150617 = sio.loadmat(os.path.abspath('../ATMO_CORR_Code/20150617_ASD_REFL.mat'))
    # Average Ground Reflectance Spectra
    tarp3_asd17 = super_resample(np.mean(asd_20150617['tarp03'],axis = 0),np.squeeze(asd_20150617['asd_wvl']),neon_wvl,neon_fwhm)/100
    tarp48_asd17 = super_resample(np.mean(asd_20150617['tarp48'],axis = 0),np.squeeze(asd_20150617['asd_wvl']),neon_wvl,neon_fwhm)/100
    veg_asd17 = super_resample(np.mean(asd_20150617['veg'],axis = 0),np.squeeze(asd_20150617['asd_wvl']),neon_wvl,neon_fwhm)/100
    ewroad_asd17 = super_resample(np.mean(asd_20150617['ewroad'],axis = 0),np.squeeze(asd_20150617['asd_wvl']),neon_wvl,neon_fwhm)/100
    nsroad_asd17 = super_resample(np.mean(asd_20150617['nsroad'],axis = 0),np.squeeze(asd_20150617['asd_wvl']),neon_wvl,neon_fwhm)/100

    asd_20150616 = sio.loadmat(os.path.abspath('../ATMO_CORR_Code/20150616_ASD_REFL.mat'))
    # Average Ground Reflectance Spectra
    tarp3_asd16 = super_resample(np.mean(asd_20150616['tarp03'],axis = 0),np.squeeze(asd_20150616['asd_wvl']),neon_wvl,neon_fwhm)/100
    tarp48_asd16 = super_resample(np.mean(asd_20150616['tarp48'],axis = 0),np.squeeze(asd_20150616['asd_wvl']),neon_wvl,neon_fwhm)/100
    veg_asd16 = super_resample(np.mean(asd_20150616['veg'],axis = 0),np.squeeze(asd_20150616['asd_wvl']),neon_wvl,neon_fwhm)/100
    ewroad_asd16 = super_resample(np.mean(asd_20150616['ewroad'],axis = 0),np.squeeze(asd_20150616['asd_wvl']),neon_wvl,neon_fwhm)/100
    nsroad_asd16 = super_resample(np.mean(asd_20150616['nsroad'],axis = 0),np.squeeze(asd_20150616['asd_wvl']),neon_wvl,neon_fwhm)/100


    asd_20150608 = sio.loadmat(os.path.abspath('../ATMO_CORR_Code/20150608_ASD_REFL.mat'))
    # Average Ground Reflectance Spectra
    tarp3_asd08 = super_resample(np.mean(asd_20150608['tarp03'],axis = 0),np.squeeze(asd_20150608['asd_wvl']),neon_wvl,neon_fwhm)/100
    tarp48_asd08 = super_resample(np.mean(asd_20150608['tarp48'],axis = 0),np.squeeze(asd_20150608['asd_wvl']),neon_wvl,neon_fwhm)/100
    veg_asd08 = super_resample(np.mean(asd_20150608['veg'],axis = 0),np.squeeze(asd_20150608['asd_wvl']),neon_wvl,neon_fwhm)/100
    ewroad_asd08 = super_resample(np.mean(asd_20150608['ewroad'],axis = 0),np.squeeze(asd_20150608['asd_wvl']),neon_wvl,neon_fwhm)/100
    nsroad_asd08 = super_resample(np.mean(asd_20150608['nsroad'],axis = 0),np.squeeze(asd_20150608['asd_wvl']),neon_wvl,neon_fwhm)/100

    plt.figure
    plt.subplot(2,2,1)
    plt.plot(neon_wvl,tarp3_asd17,'k--')
    plt.plot(neon_wvl,tarp3_asd16,'k:')
    plt.plot(neon_wvl,tarp3_asd08,'k')

    plt.subplot(2,2,2)
    plt.plot(neon_wvl,tarp48_asd17,'b--')
    plt.plot(neon_wvl,tarp48_asd16,'b:')
    plt.plot(neon_wvl,tarp48_asd08,'b')

    plt.subplot(2,2,3)
    plt.plot(neon_wvl,veg_asd17,color = [34/255,139/255,34/255],linestyle = '--')
    plt.plot(neon_wvl,veg_asd16,color = [34/255,139/255,34/255],linestyle = ':')
    plt.plot(neon_wvl,veg_asd08,color = [34/255,139/255,34/255])

    plt.subplot(2,2,4)
    plt.plot(neon_wvl,ewroad_asd17,color = [156/255,102/255,31/255],linestyle = '--')
    plt.plot(neon_wvl,ewroad_asd16,color = [156/255,102/255,31/255],linestyle = ':')
    plt.plot(neon_wvl,ewroad_asd08,color = [156/255,102/255,31/255])
    plt.legend(['Jun 17','Jun 16','Jun 8'])

    asdtemp = sio.loadmat(asdfile)
    asd_wvl = np.squeeze(asdtemp['asd_wvl'])
    # asd_wvl = asd_wvl'
#    if date == '20150608':
#        asd_wvl = asd_wvl.T

    # Set Color Tables
    enhanced = [178/255,24/255,43/255] # Dark Red
    light_enhanced = [253/255,219/255,199/255] # Light Red
    intermediate = [211/255,84/255,0/255]   # Dark Orange
    light_intermediate = [235/255,152/255,78/255]   # Light Orange
    standard = [33/255,102/255,172/255] # Blue
    light_standard = [146/255,197/255,222/255]  # Light Blue

    standard_adj = [152/255,78/255,163/255] # Purple

    albedo_color = [153/255,163/255,164/255]    # Gray
    asd_color = [0,0,0]

    target_list = [dict([('name','3% Tarp'),('coord',tarp3_coord),
                         ('fname',date+'_tarp03'),('ref',np.mean(asdtemp['tarp03'],axis = 0)),
                         ('range',(0,0.1))]),
                   dict([('name','48% Tarp'),('coord',tarp48_coord),
                         ('fname',date+'_tarp48'),('ref',np.mean(asdtemp['tarp48'],axis = 0)),
                         ('range',(0,1))]),
                   dict([('name','Vegetation'),('coord',veg_coord),
                         ('fname',date+'_veg'),('ref',np.mean(asdtemp['veg'],axis = 0)),
                         ('range',(0,0.4))]),
                   dict([('name','E-W Road'),('coord',EWroad_coord),
                         ('fname',date+'_ewroad'),('ref',np.mean(asdtemp['ewroad'],axis = 0)),
                         ('range',(0,0.6))]),
                   dict([('name','North-South Road'),('coord',NSroad_coord),
                         ('fname',date+'_nsroad'),('ref',np.mean(asdtemp['nsroad'],axis = 0)),
                         ('range',(0,0.6))])]

    # rad0 = super_resample(r0,wvl_7sc,neon_wvl,neon_fwhm)
    # Ifup0 = super_resample(rho_a_I*Ifdn_obs_trim_NIS_SSIM,flx_wvl,neon_wvl2[:,x],neon_fwhm)*Ifdn_obs_trim
    rmse_vals = dict()

    for target in target_list:
        print(target['name'])

        target['spectra'] = np.reshape(nis_datacube.data_cube[target['coord'][1,0]:target['coord'][1,1],target['coord'][0,0]:target['coord'][0,1],:],(-1,426))
        target['resp_func'] = dict([('wvl',neon_wvl),('fwhm',neon_fwhm),('wvl0',neon_wvl0)])

        # Time Dependent Solver
        obstime =  nistime[target['coord'][1,0]:target['coord'][1,1],target['coord'][0,0]:target['coord'][0,1]].flatten()
        temp_Ifdn = list()
        temp_Ifup = list()

        Iup = np.zeros((target['spectra'].shape[0],target['resp_func']['wvl'].shape[0]))
        Idn = np.zeros((target['spectra'].shape[0],target['resp_func']['wvl'].shape[0]))

        for i in range(len(obstime)):
#            obs_time = np.mean(nistime[target['coord'][1,0]:target['coord'][1,1],target['coord'][0,0]:target['coord'][0,1]])
            ssirtime_ind = np.argmin(np.abs(ssim['ssirtime']-obstime[i]))
            temp_Ifdn.append(ssim['zspect'][ssirtime_ind,:])
            temp_Ifup.append(ssim['nspect'][ssirtime_ind,:])

            # For R_stand_adj:
            Iup[i,:] = super_resample(ssim['nspect'][ssirtime_ind,:],nwvl,target['resp_func']['wvl'], target['resp_func']['fwhm'])

            # For R_enhan:
            Idn[i,:] = super_resample(ssim['zspect'][ssirtime_ind,:],zwvl,target['resp_func']['wvl'], target['resp_func']['fwhm'])

#        Ifdn_obs = dict([('wvl',ssim['zwvl']),('If_dn',np.array(temp_Ifdn))])
#        Ifup_obs = dict([('wvl',ssim['nwvl']),('If_up',np.array(temp_Ifup))])


#        # Note - I should rewrite to find the closest SSIR observation to each pixel's time, not the mean time
#        obs_time = np.mean(nistime[target['coord'][1,0]:target['coord'][1,1],target['coord'][0,0]:target['coord'][0,1]])
#        ssirtime_ind = np.argmin(np.abs(ssim['ssirtime']-obs_time))
#        Ifdn_obs = dict([('wvl',ssim['zwvl']),('If_dn',ssim['zspect'][ssirtime_ind,:])])
#        Ifup_obs = dict([('wvl',ssim['nwvl']),('If_up',ssim['nspect'][ssirtime_ind,:])])
#        target['spectra'] = np.reshape(nis_datacube.data_cube[target['coord'][1,0]:target['coord'][1,1],target['coord'][0,0]:target['coord'][0,1],:],(-1,426))
#        target['resp_func'] = dict([('wvl',neon_wvl),('fwhm',neon_fwhm),('wvl0',neon_wvl0)])

        # Standard Atmospheric Correction using Whole Atmosphere MODTRAN
        R_stand = standard_correction(target, rad0_neon, I0, baseline_acd, mu)
        R_stand = dict([('spectra',R_stand),('wvl',neon_wvl),('color',standard),('name','Standard')])
        # Whole Atmosphere MODTRAN with adjacency correction using Upwelling Irradiance
#        Iup = super_resample(Ifup_obs['If_up'],Ifup_obs['wvl'],target['resp_func']['wvl'], target['resp_func']['fwhm'])
        R_stand_adj = adjacency_correction(target, I0, Iup, rad0_neon, Iup0_neon, baseline_acd, mu = mu, type = 'Standard')
        R_stand_adj = dict([('spectra',R_stand_adj),('wvl',neon_wvl),('color',standard_adj),('name','Standard with Adjacency')])

        # Enhanced Atmospheric Correction using Downwelling Irradiance
#        Idn = super_resample(Ifdn_obs['If_dn'],Ifdn_obs['wvl'],target['resp_func']['wvl'], target['resp_func']['fwhm'])
        rad0_obs = super_resample(rho_a_R,wvl_7sc,neon_wvl,neon_fwhm) * Idn
        R_enhan = irrad_correction(target, rad0_obs, Idn, flight_acd)
        R_enhan = dict([('spectra',R_enhan),('wvl',neon_wvl),('color',intermediate),('name','Enhanced')])
        # Enhanced Atmospheric Correction Using Up- and Downwelling Irradiances
        Iup0_obs = super_resample(rho_a_I,wvl_7sc,neon_wvl,neon_fwhm) * Idn
        R_enhan_adj = adjacency_correction(target, Idn, Iup, rad0_obs, Iup0_obs, flight_acd, type = 'Irradiance')
        R_enhan_adj = dict([('spectra',R_enhan_adj),('wvl',neon_wvl),('color',enhanced),('name','Enhanced with Adjacency')])

        R_albedo = albedo(target,Idn)
        R_albedo = dict([('spectra',R_albedo),('wvl',neon_wvl),('color',albedo_color),('name','Flight Level Albedo')])

        asd = dict([('spectra',target['ref']/100),('wvl',asd_wvl),('color',asd_color),('name','ASD-GroundTruth')])
        asd_neon = dict([('spectra',super_resample(target['ref'],asd_wvl,neon_wvl,neon_fwhm)/100),('wvl',asd_wvl),('color',asd_color),('name','ASD-GroundTruth')])

        plot_results(asd,R_albedo,R_stand, R_stand_adj, R_enhan, R_enhan_adj, title = target['fname'], yrng = target['range'], save = True)
        rmse_vals[target['name']] = calc_rmse(asd_neon, neon_wvl, R_albedo, R_stand, R_stand_adj, R_enhan, R_enhan_adj)

    plt.figure()
    plt.plot(neon_wvl0,super_resample(baseline_flx['downwelling'][flight_ind]/mu,baseline_flx['wvl'],neon_wvl0,neon_fwhm0)*10000,'r')
    plt.plot(zwvl,ssim['zspect'][ssirtime_ind,:],'k')

    model_Ifdn = super_resample(baseline_flx['downwelling'][flight_ind]/mu,baseline_flx['wvl'],neon_wvl0,neon_fwhm0)*10000
    obs_Ifdn = ssim['zspect'][ssirtime_ind,:]

    sio.savemat(date.format('%s')+'_downirradcomp.mat',{'model_Ifdn':model_Ifdn,'obs_Ifdn':obs_Ifdn,'neon_wvl0':neon_wvl0,'zwvl':zwvl})


    keys = tuple(rmse_vals[target['name']].keys())
    fid = open('RMSE_Output.txt','w')
    fid.write('atmo_corr.py RMSE Outputs\nTarget\t{0[0]!s}\t{0[1]!s}\t{0[2]!s}\t{0[3]!s}\t{0[4]!s}\n'.format(keys))

    for line in rmse_vals:
        fid.write('{0!s}\t{1[0]:f}\t{1[1]:f}\t{1[2]:f}\t{1[3]:f}\t{1[4]:f}\n'.format(line,list(rmse_vals[line].values())))

    fid.close()

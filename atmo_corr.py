#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
NAME:
Created By: Logan Wright
Created On: Dec 19, 2018

Code was adopted from original MATLAB Code

DESCRIPTION:


'''


import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

def standard_correction(Lobs,L0,I0,transmittance,spherical_albedo,mu):
    R = np.pi * (Lobs - L0) / (spherical_albedo * np.pi * (Lobs - L0) + mu * I0 * ( (Ts + ts) * (T + t) ) )
    return R

def irrad_correction(Lobs,L0,Ifup0,Ifdn,transmittance,spherical_albedo,mu = None):
    R = np.pi * (Lobs - L0) / (spherical_albedo * np.pi * (Lobs - L0) + Ifdn .* (T + t)**2)
    return R

def albedo(Lobs,Idn):
    A = Rob * np.pi / Ifdn_obs_trim
    return A

def adjacency_correction(Lobs,L0,Ifup0,Ifdn,Ifup,transmittance,spherical_albedo,mu = None):
    Rh_bar = (Ifup - Ifup0) / ((Ifdn * (T + t)**2) + (Ifup - Ifup0) .* sph_alb_1km))

    R = (Lobs - L0 - ((Ifdn * (t + T) * t) / np.pi) * (Rh_bar / (1 - Rh_bar * spherical_albedo))) *
            ((np.pi .* (1 - Rh_bar * spherical_albedo)) / (Ifdn * (t + T) * T)) ;

    return R

def mu_correction():
    # Mu Correction
    tarp48_ret_noadj_mu(j,:) = pi .* (Rob - rad0) ./ (sph_alb_1km .* pi .*(Rob - rad0) + Ifdn_obs_trim .* (t_1km.^2 + (1 + mu) * t_1km .* sqrt(Tso_1km) + mu.*Tso_1km));

    rho_bar_ret_mu = (Ifup_obs_trim - Ifup0) ./ ((Ifdn_obs_trim .*(t_1km.^2 + (1 + mu) * t_1km .* sqrt(Tso_1km) + mu.*Tso_1km) + (Ifup_obs_trim - Ifup0) .* sph_alb_1km));
    R = (Rob - rad0 - (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* t_1km) ./ pi .* (rho_bar_ret_mu ./ (1 - rho_bar_ret_mu .* sph_alb_1km))) .* (pi ./ (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* mu .* sqrt(Tso_1km))) .* (1 - rho_bar_ret_mu .* sph_alb_1km);

def load_ASD(filepath):

    return (wvl,R)

def atmcorr(wvl,radiance,model,irrad_up = None,irrad_down = None):
    n = np.size(radiance)
    # Observation Retrieval
    [~,ssirtime_ind] = np.min(np.abs(ssirtime-obs_time[j]));

    Ifdn_obs_trim = super_resample(zspect(ssirtime_ind,:)',zwvl,neon_wvl,neon_fwhm);
    Ifup_obs_trim = super_resample(nspect(ssirtime_ind,:)',nwvl,neon_wvl,neon_fwhm);

    Ifup0 = super_resample(rho_a_I,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;
    rad0 = super_resample(rho_a_R,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;

    tarp3_ret_noadj(j,:) = pi .* (Rob - rad0) ./ (sph_alb_1km .* pi .* ...
            (Rob - rad0) + Ifdn_obs_trim .* (t_1km.^2 + 2 * t_1km .* sqrt(Tso_1km) + Tso_1km));

    rho_bar_ret = (Ifup_obs_trim - Ifup0) ./ ((Ifdn_obs_trim .* ...
            ((t_1km + sqrt(Tso_1km)).^2) + (Ifup_obs_trim - Ifup0) .* sph_alb_1km));
    tarp3_ret(j,:) = (Rob - rad0 - (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* t_1km) ./ ...
            pi .* (rho_bar_ret ./ (1 - rho_bar_ret .* sph_alb_1km))) .* ...
            (pi ./ (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* sqrt(Tso_1km))) .* (1 - rho_bar_ret .* sph_alb_1km);

    flight_alb(j,:) = Rob*pi./Ifdn_obs_trim;

    plt.figure;hold on
    # plot(neon_wvl,tarp3_ret_TOA,'Color',light_standard)
    # plot(neon_wvl,tarp3_ret_noadj,'Color',light_intermediate) # Intermediate
    # plot(neon_wvl,tarp3_ret,'Color',light_enhanced)
    plt.plot(asd_wvl,mean(tarp03,1)'./100,'k','LineWidth',1.5)
    # plot(neon_wvl,mean(tarp3_ret_TOA,1),'Color',standard,'LineWidth',1.5)
    plt.plot(neon_wvl,mean(flight_alb,1),'Color',[153,163,164]./255,'LineWidth',1.5)
    plt.plot(neon_wvl,mean(tarp3_ret_TOA_B,1),'-','Color',standard,'LineWidth',1.5)
    plt.plot(neon_wvl,mean(tarp3_ret_noadj,1),'Color',intermediate,'LineWidth',1.5) # Intermediate
    plt.plot(neon_wvl,mean(tarp3_ret,1),'Color',enhanced,'LineWidth',1.5)
    plt.ylabel('Reflectance','FontSize',12);xlabel('Wavelength [nm]','FontSize',12);
    plt.axis([350,1025,0,0.1])
    plt.set(gca,'FontSize',12,'XColor','k','YColor','k');box on;
    fig = plt.gcf;
    plt.set(fig,'PaperUnits','inches','PaperPosition',[0 0 4 4]);
    plt.savefig(fig,'-depsc','-r300','TestTarp3Retrieval.eps')

    wv_bands = xor(neon_wvl < 928, neon_wvl > 989);
    wv_bands = xor(wv_bands,neon_wvl > 1095);

    tarp03_r = super_resample(mean(tarp03,1),asd_wvl,neon_wvl,neon_fwhm)'./100;

    plt.figure;hold on;plot(neon_wvl,(mean(tarp3_ret,1)-tarp03_r)*100,'Color',enhanced)
    plt.plot(neon_wvl,(mean(tarp3_ret_noadj,1)-tarp03_r)*100,'Color',intermediate)
    plt.plot(neon_wvl,(mean(tarp3_ret_TOA_B,1)-tarp03_r)*100,'Color',standard);grid on;
    RMSE_tarp3 = [sqrt(sum(((mean(tarp3_ret,1)-tarp03_r)*100).^2)./numel(tarp03_r));
        sqrt(sum(((mean(tarp3_ret_noadj,1)-tarp03_r)*100).^2)./numel(tarp03_r));
        sqrt(sum(((mean(tarp3_ret_TOA_B,1)-tarp03_r)*100).^2)./numel(tarp03_r))];
    RMSE_tarp3_no_wv = [sqrt(sum(((mean(tarp3_ret(:,wv_bands),1)-tarp03_r(wv_bands))*100).^2)./numel(tarp03_r(wv_bands)));
        sqrt(sum(((mean(tarp3_ret_noadj(:,wv_bands),1)-tarp03_r(wv_bands))*100).^2)./numel(tarp03_r(wv_bands)));
        sqrt(sum(((mean(tarp3_ret_TOA_B(:,wv_bands),1)-tarp03_r(wv_bands))*100).^2)./numel(tarp03_r(wv_bands)))];
    plt.title('Tarp 3#')




if __name__ == '__main__':
    # Below Cloud Surface Reflectance Retrieval

    ## Section 1: Load Data and Set constants
    # Set Date
    date = '20150608'

    # This section includes filenames and constants for each day/flight line
    if date == '20150608':
        NISfile = '/Users/wrightad/Documents/Data/NEON/Flight_ROIs/20150608/NIS01_20150608_165842_rdn'
        tarp3_coord = np.array(((304,307),(364,369)))
        tarp48_coord = np.array(((313,319),(365,370)))
        veg_coord = np.array(((324,328),(380,420)))
        EWroad_coord = np.array(((310,348),(354,356)))
        NSroad_coord = np.array(((334,336),(317,354)))
        filename = ['/Users/wrightad/Documents/Modtran/NEON_ATMCORR/NEON_20150608_Baseline',
            '/Users/wrightad/Documents/Modtran/NEON_ATMCORR/NEON_20150608_1kmAtmo_40H2O']
        SunEllipticFactor = 0.97083    # Determined by Day of Year
        mu = np.cosd(31)                  # Determined by Time of Day
    #     SSIRfile = '/Users/wrightad/Documents/Data/NEON/FlightData/20150608_SSIR.mat'
    #     SSIRfile = '/Users/wrightad/Documents/Data/NEON/FlightData/20150608_SSIR_CALIB_0602.mat'
        SSIRfile = '/Users/wrightad/Documents/Data/NEON/FlightData/20150608_SSIR_CALIBSPECT_Jul2017.mat'
        asdfile = '/Users/wrightad/Documents/Data/NEON/GroundData/ASD_Refl/20150608_ASD_REFL.mat'

    elif date == '20150616':
        NISfile = '/Users/wrightad/Documents/Data/NEON/Flight_ROIs/20150616/NIS01_20150616_151613_rdn'
        tarp3_coord = np.array(((140,145),(232,237)))
        tarp48_coord = np.array(((150,155)(232,237)))
        veg_coord = np.array(((155,160),(253,280)))
        EWroad_coord = np.array(((148,185),(220,224)))
        NSroad_coord = np.array(168,172),(182,218)))
        filename = ('/Users/wrightad/Documents/Modtran/NEON_ATMCORR/NEON_20150616_Baseline',
            '/Users/wrightad/Documents/Modtran/NEON_ATMCORR/NEON_20150616_259mAtmo_30H2O')
        SunEllipticFactor = 0.96950    # Determined by Day of Year
        mu = np.cosd(50.2)                # Determined by Time of Day
        SSIRfile = '/Users/wrightad/Documents/Data/NEON/FlightData/20150616_SSIR_CALIBSPECT_20150616FIELDB.mat'
        asdfile = '/Users/wrightad/Documents/Data/NEON/GroundData/ASD_Refl/20150616_ASD_REFL.mat'
    #     asdfile = '/Users/wrightad/Documents/Data/NEON/GroundData/ASD_Refl/20150608_ASD_REFL.mat'

    elif date == '20150617':
        NISfile = '/Users/wrightad/Documents/Data/NEON/Flight_ROIs/20150617/NIS01_20150617_203142_rdn'
        tarp3_coord = np.array(((174,180),(200,205)))
        tarp48_coord = np.array(((183,190),(200,207)))
        veg_coord = np.array(((196,199),(220,245)))
        EWroad_coord = np.array(((173,226),(190,192)))
        NSroad_coord = np.array(((208,211),(159,190)))
        filename = ('/Users/wrightad/Documents/Modtran/NEON_ATMCORR/NEON_20150617_Baseline',
            '/Users/wrightad/Documents/Modtran/NEON_ATMCORR/NEON_20150617_1kmAtmo_40H2O')
        SunEllipticFactor = 0.96933    # Determined by Day of Year
        mu = cosd(25.0)                # Determined by Time of Day
        SSIRfile = '/Users/wrightad/Documents/Data/NEON/FlightData/20150617_SSIR.mat'
        asdfile = '/Users/wrightad/Documents/Data/NEON/GroundData/ASD_Refl/20150617_ASD_REFL.mat'

    N = len(filename)

    conv = 10000;   # Scale Factor to convert [W cm^-2 nm^-1] to [W m^-2 nm^-2], used for .flx MODTRAN output file

    # Load NEON Instrument Response Function
    # respfunc = np.genfromtxt('/Users/wrightad/Documents/Data/NEON/nis2013_fmed2_2013.dat');
    # neon_wvl0 = respfunc[:,2]+1.28;
    # neon_fwhm0 = respfunc[:,3];

    # Set Wavelength Ranges
    # SSIR Index (25) to Neon Wvl (364) - Zenith
    # SSIR Index (25) to NEON Wvl (152) - Nadir
    neon_wvl_zen = neon_wvl0[1:364];neon_fwhm_zen = neon_fwhm0[1:364];
    neon_wvl_nad = neon_wvl0[1:152];neon_fwhm_nad = neon_fwhm0[1:152];
    # neon_wvl_nad = neon_wvl0(1:125);neon_fwhm_nad = neon_fwhm0(1:125); # This
    # Truncates at 1000 nm (near Si/InGaAs joining)
    neon_wvl = neon_wvl_nad;
    neon_fwhm = neon_fwhm_nad;
    L = len(neon_wvl0);

    # Load MODTRAN Irradiance
    data = np.genfromtxt('/Users/wrightad/Documents/Modtran/DATA/SUN01med2irradwnNormt.dat',' ',2);
    data = data.data(2:end,:);
    #     Convert to wavelength increments
    I0[:,1] = 1e7./data[:,1]
    I0[:,2] = data[:,2].*(1./I0[:,1].^2).*1e11
    I0 = I0[I0[]:,1]<3000,:]

    I0_orig = I0

    #     Convolve to NEON Wavelength
    I0 = SunEllipticFactor.*super_resample(I0[:,2],I0[:,1],neon_wvl,neon_fwhm)



    # Load in MODTRAN Output
    [t_freq, t_ts, t_Tso, t_t, t_T, t_sph_alb] = load_acd(filename[0]+'.acd')
    t_Ts = t_Tso./t_T
    t_Ts(np.isnan(t_Ts)) = 0
    t_Ts(np.isinf(t_Ts)) = 0
    wvl = 1e7./(t_freq)
    T = super_resample(t_T,wvl,neon_wvl,neon_fwhm)
    t = super_resample(t_t,wvl,neon_wvl,neon_fwhm)
    Ts = super_resample(t_Ts,wvl,neon_wvl,neon_fwhm)
    ts = super_resample(t_ts,wvl,neon_wvl,neon_fwhm)
    Tso = super_resample(t_Tso,wvl,neon_wvl,neon_fwhm)
    sph_alb = super_resample(t_sph_alb,wvl,neon_wvl,neon_fwhm)

    # Load MODTRAN Irradiance
    flx_data = np.genfromtxt(filename[0]+'.flx'))
    flx_wvl = flx_data.data[:,1]
    # flx_fup = conv.*super_resample(flx_data.data(:,8),flx_data.data(:,1),neon_wvl,neon_fwhm)
    flx_fdn = flx_data.data[:,33]+flx_data.data[:,34]
    flx_fup = flx_data.data[:,32]

    # Load MODTRAN Radiance
    data_7sc = np.genfromtxt(filename[0]+'.7sc',' ',11)
    wvl_7sc = data_7sc.data[0:-1,0]
    r0 = data_7sc.data[:,9]
    Iup = data_7sc.data[:,10]

    r0 = r0[0:-1,:]./100
    Iup = Iup[0:-1,:]

    # Calculate Atmospheric Reflectance
    rho_a_I = flx_fup./flx_fdn
    rho_a_R = np.flip(r0)./(1e4*flx_fdn)

    rad0 = super_resample(r0,wvl_7sc,neon_wvl,neon_fwhm)
    Iup0_neon = super_resample(Iup,wvl_7sc,neon_wvl,neon_fwhm)

    [~, t_ts_1km, t_Tso_1km, t_t_1km, t_T_1km, t_sph_alb_1km] = load_acd(strcat(filename{2},'.acd'))
    t_Ts_1km = t_Tso_1km./t_T_1km
    t_Ts_1km[np.isnan(t_Ts_1km)] = 0
    t_Ts_1km[np.isinf(t_Ts_1km)] = 0
    T_1km = super_resample(t_T_1km,wvl,neon_wvl,neon_fwhm)
    t_1km = super_resample(t_t_1km,wvl,neon_wvl,neon_fwhm)
    Ts_1km = super_resample(t_Ts_1km,wvl,neon_wvl,neon_fwhm)
    ts_1km = super_resample(t_ts_1km,wvl,neon_wvl,neon_fwhm)
    Tso_1km = super_resample(t_Tso_1km,wvl,neon_wvl,neon_fwhm)
    sph_alb_1km = super_resample(t_sph_alb_1km,wvl,neon_wvl,neon_fwhm)

    flx_data_1km = importdata(filename[1]+'.flx')
    t_T_derived_1km = flx_data_1km.data[:,4]./flx_data_1km.data[:,34]
    t_Ts_derived_1km = flx_data_1km.data[:,4]./flx_data_1km.data[:,-1]
    t_T_derived = flx_data.data[:,4]./flx_data.data[:,34]
    t_Ts_derived = flx_data.data[:,4]./flx_data.data[:,-1]

    T_derived_1km = super_resample(flx_data_1km.data[:,4]./flx_data_1km.data[:,34],flx_data_1km.data[:,1],neon_wvl,neon_fwhm)
    Ts_derived_1km = super_resample(flx_data_1km.data[:,4]./flx_data_1km.data[:,-1],flx_data_1km.data[:,1],neon_wvl,neon_fwhm)
    T_derived = super_resample(flx_data.data[:,4]./flx_data.data[:,34],flx_data.data[:,1],neon_wvl,neon_fwhm)
    Ts_derived = super_resample(flx_data.data[:,4]./flx_data.data[:,end],flx_data.data[:,1],neon_wvl,neon_fwhm)

    # Load NIS Data
    data = Load_NIS(NISfile)
    data = data./100   # Divide by 100 to convert to W/sr/m^2/nm
    metadata = Load_NIS(NISfile+'_obs_ort')
    nistime = metadata[:,:,10]

    # Load SSIR Data
    sio.loadmat(SSIRfile)
    zwvl = zwvl-2.5
    nwvl = nwvl

    ssim_zwvl = zwvl
    ssim_zfwhm0 = [365.25, 7.67;
                   404.84, 8.11;
                   435.83, 7.25;
                   546.07, 7.88;
                   912.30, 10.86;
                   965.78, 15.97;
                   1244.30, 12.70;
                   1694.06, 9.75;
                   1791.46, 10.79;
                   2061.62, 17.58];
    ssim_zfwhm = np.zeros(439,1);
    ssim_zfwhm(1:192) = np.interp(ssim_zfwhm0(1:5,1),ssim_zfwhm0(1:5,2),ssim_zwvl(1:192),'linear','extrap');
    ssim_zfwhm(193:end) = np.interp(ssim_zfwhm0(6:end,1),ssim_zfwhm0(6:end,2),ssim_zwvl(193:end),'linear','extrap');
    ssim_nwvl = nwvl;
    ssim_nfwhm0 = [365.25, 7.50;
                   404.84, 7.85;
                   435.83, 7.28;
                   546.07, 8.16;
                   912.30, 10.75];
    ssim_nfwhm = interp1(ssim_nfwhm0(:,1),ssim_nfwhm0(:,2),ssim_nwvl,'linear','extrap');

    ssim_wvl = ssim_zwvl(26:225);

    # Demonstrate Sampling Issue
    I0_ssim = SunEllipticFactor.*super_resample(I0_orig(:,2),I0_orig(:,1),zwvl,ssim_zfwhm);
    Ig_ssim = I0_ssim.*super_resample(t_Ts+t_ts,wvl,zwvl,ssim_zfwhm);
    Ig_NIS = I0.*super_resample(t_Ts+t_ts,wvl,neon_wvl,neon_fwhm);
    Ig_ssim_resampled = super_resample(I0_ssim(17:end).*super_resample(t_Ts+t_ts,wvl,zwvl(17:end),ssim_zfwhm(17:end)),zwvl(17:end),neon_wvl,neon_fwhm);
    plt.figure;
    # plot(I0_orig(:,1),I0_orig(:,2),'k')
    plt.plot(neon_wvl,Ig_NIS,'Color',[0,0,0])
    plt.plot(zwvl,Ig_ssim,'Color',[33,102,172]./255)
    plt.plot(neon_wvl,Ig_ssim_resampled,'Color',[178,24,43]./255)



    plt.legend('NIS Sampling','SSIM Sampling','SSIM Resampled to NIS')
    plt.ylabel('Spectral Irradiance [W m^{-2} nm^{-1}]');
    plt.xlabel('Wavelength [nm]')
    plt.title('Modelled Downwelling Irradiance')
    plt.set(gca,'FontSize',18)
    plt.axis([720,800,0.7,1.3]);grid on;box on
    plt.savefig('IrradianceDiffRespFunc_O2.svg','-dsvg')
    plt.axis([350,1025,0,1.8]);
    plt.savefig('IrradianceDiffRespFunc_Full.svg','-dsvg')

    plt.figure;
    plt.plot(neon_wvl,Ig_NIS-Ig_ssim_resampled,'k')
    plt.ylabel('Spectral Irradiance [W m^{-2} nm^{-1}]');
    plt.xlabel('Wavelength [nm]')
    plt.title('Difference Between NIS and SSIM Sampling')
    plt.axis([720,800,-0.15,0.1]);grid on;box on
    plt.set(gca,'FontSize',18)
    plt.print('DiffRespFunc_O2.svg','-dsvg')

    # Load Ground Reflectances
    sio.loadmat('/Users/wrightad/Documents/Data/NEON/GroundData/ASD_Refl/20150617_ASD_REFL.mat')
    # Average Ground Reflectance Spectra
    tarp3_asd17 = super_resample(mean(tarp03,1)',asd_wvl',neon_wvl,neon_fwhm)./100;
    tarp48_asd17 = super_resample(mean(tarp48,1)',asd_wvl',neon_wvl,neon_fwhm)./100;
    veg_asd17 = super_resample(mean(veg,1)',asd_wvl',neon_wvl,neon_fwhm)./100;
    ewroad_asd17 = super_resample(mean(ewroad,1)',asd_wvl',neon_wvl,neon_fwhm)./100;
    nsroad_asd17 = super_resample(mean(nsroad,1)',asd_wvl',neon_wvl,neon_fwhm)./100;

    sio.loadmat('/Users/wrightad/Documents/Data/NEON/GroundData/ASD_Refl/20150616_ASD_REFL.mat')
    # Average Ground Reflectance Spectra
    tarp3_asd16 = super_resample(mean(tarp03,1)',asd_wvl',neon_wvl,neon_fwhm)./100;
    tarp48_asd16 = super_resample(mean(tarp48,1)',asd_wvl',neon_wvl,neon_fwhm)./100;
    veg_asd16 = super_resample(mean(veg,1)',asd_wvl',neon_wvl,neon_fwhm)./100;
    ewroad_asd16 = super_resample(mean(ewroad,1)',asd_wvl',neon_wvl,neon_fwhm)./100;
    nsroad_asd16 = super_resample(mean(nsroad,1)',asd_wvl',neon_wvl,neon_fwhm)./100;


    sio.loadmat('/Users/wrightad/Documents/Data/NEON/GroundData/ASD_Refl/20150608_ASD_REFL.mat')
    # Average Ground Reflectance Spectra
    tarp3_asd08 = super_resample(mean(tarp03,1)',asd_wvl,neon_wvl,neon_fwhm)./100;
    tarp48_asd08 = super_resample(mean(tarp48,1)',asd_wvl,neon_wvl,neon_fwhm)./100;
    veg_asd08 = super_resample(mean(veg,1)',asd_wvl,neon_wvl,neon_fwhm)./100;
    ewroad_asd08 = super_resample(mean(ewroad,1)',asd_wvl,neon_wvl,neon_fwhm)./100;
    nsroad_asd08 = super_resample(mean(nsroad,1)',asd_wvl,neon_wvl,neon_fwhm)./100;

    plt.figure;plt.subplot(2,2,1);plt.plot(neon_wvl,tarp3_asd17,'k--');hold on;
    plt.plot(neon_wvl,tarp3_asd16,'k:');
    plt.plot(neon_wvl,tarp3_asd08,'k');
    plt.subplot(2,2,2);plt.plot(neon_wvl,tarp48_asd17,'b--');hold on
    plt.plot(neon_wvl,tarp48_asd16,'b:');
    plt.plot(neon_wvl,tarp48_asd08,'b');
    plt.subplot(2,2,3);plt.plot(neon_wvl,veg_asd17,'Color',[34,139,34]./255,'LineStyle','--');hold on
    plt.plot(neon_wvl,veg_asd16,'Color',[34,139,34]./255,'LineStyle',':');
    plt.plot(neon_wvl,veg_asd08,'Color',[34,139,34]./255);
    plt.subplot(2,2,4);plot(neon_wvl,ewroad_asd17,'Color',[156,102,31]./255,'LineStyle','--');hold on
    plt.plot(neon_wvl,ewroad_asd16,'Color',[156,102,31]./255,'LineStyle',':');
    plt.plot(neon_wvl,ewroad_asd08,'Color',[156,102,31]./255);
    plt.legend('Jun 17','Jun 16','Jun 8')

    # plot(neon_wvl,nsroad_asd17,'Color',[156,102,31]./255);
    # plot(neon_wvl,nsroad_asd16,'Color',[156,102,31]./255);
    # plot(neon_wvl,nsroad_asd08,'Color',[156,102,31]./255);

    sio.loadmat(asdfile)
    # asd_wvl = asd_wvl';
    if date == '20150608':
        asd_wvl = asd_wvl.T

    # Modelled vs Measured Irradiances
    obs_time = np.reshape(nistime(tarp3_coord(2,1):tarp3_coord(2,2),tarp3_coord(1,1):tarp3_coord(1,2)),1,[]);
    [~,ssirtime_ind] = np.min(np.abs(ssirtime-obs_time(1)));
    Ifdn_obs_trim = zspect(ssirtime_ind,:);

    # Extract Pixel/Observation Radiance Spectrum
    rad0 = super_resample(r0,wvl_7sc,neon_wvl,neon_fwhm);
    #
    # Ifup_TOA = flx_fup + (mu.*I0.*(Tsts_Tt).*veg_asd08)./(1-sph_alb.*veg_asd08);
    # Ifdn_TOA_B = mu.*I0.*(ts(:,2)+Ts(:,2))+Ifup_TOA.*sph_alb;
    #
    #     # Observation Retrieval

    #

    #     Ifup_obs_trim = super_resample(nspect(ssirtime_ind,:)',nwvl,neon_wvl,neon_fwhm);
    #

    plt.figure;hold on;
    plt.plot(neon_wvl0,super_resample(flx_fdn/mu,flx_wvl,neon_wvl0,neon_fwhm0)*10000,'r');
    plt.plot(zwvl,Ifdn_obs_trim,'k');

    model_Ifdn = super_resample(flx_fdn/mu,flx_wvl,neon_wvl0,neon_fwhm0)*10000;
    obs_Ifdn = Ifdn_obs_trim;

    sio.save(date.format('%s')+'_downirradcomp.mat'),{'model_Ifdn':model_Ifdn,'obs_Ifdn':obs_Ifdn,'neon_wvl0':neon_wvl0,'zwvl':zwvl)

    ## Section 2: Processing
    # Extract Ground Spectra from NIS Scene
    tarp3_spc = np.reshape(data(tarp3_coord(2,1):tarp3_coord(2,2),tarp3_coord(1,1):tarp3_coord(1,2),:),[],426);
    tarp48_spc = np.reshape(data(tarp48_coord(2,1):tarp48_coord(2,2),tarp48_coord(1,1):tarp48_coord(1,2),:),[],426);
    veg_spc = np.reshape(data(veg_coord(2,1):veg_coord(2,2),veg_coord(1,1):veg_coord(1,2),:),[],426);
    EWroad_spc = np.reshape(data(EWroad_coord(2,1):EWroad_coord(2,2),EWroad_coord(1,1):EWroad_coord(1,2),:),[],426);
    NSroad_spc = np.reshape(data(NSroad_coord(2,1):NSroad_coord(2,2),NSroad_coord(1,1):NSroad_coord(1,2),:),[],426);

    # Set Color Tables
    enhanced = [178,24,43]./255;
    light_enhanced = [253,219,199]./255;
    intermediate = [211,84,0]./255;
    light_intermediate = [235,152,78]./255;
    standard = [33,102,172]./255;
    light_standard = [146,197,222]./255;

    # Full Resolution Calculation
    Tsts_Tt = super_resample((t_Ts + t_ts).*(t_T+t_t),wvl,neon_wvl,neon_fwhm);
    Tt2 = super_resample((t_T_1km + t_t_1km).^2,wvl,neon_wvl,neon_fwhm);
    t2Tso2 = super_resample(t_t_1km.^2 + 2 * t_t_1km .* sqrt(t_Tso_1km) + t_Tso_1km,wvl,neon_wvl,neon_fwhm);
    tTt = super_resample((t_t_1km + t_T_1km) .* t_t_1km,wvl,neon_wvl,neon_fwhm);
    tTT = super_resample((t_t_1km + t_T_1km) .* t_T_1km,wvl,neon_wvl,neon_fwhm);
    tTso2 = super_resample((t_t_1km + sqrt(t_Tso_1km)).^2,wvl,neon_wvl,neon_fwhm);

    tTsot = super_resample((t_t_1km + sqrt(t_Tso_1km)) .* t_t_1km,wvl,neon_wvl,neon_fwhm);
    tTsoTso = super_resample((t_t_1km + sqrt(t_Tso_1km)) .* sqrt(t_Tso_1km),wvl,neon_wvl,neon_fwhm);

    # Derived Full Resolution Calculation
    Tsts_Tt_derived = super_resample((t_Ts_derived + super_resample(t_ts,wvl,flx_wvl,2.5)).*(t_T_derived_1km+super_resample(t_t,wvl,flx_wvl,2.5)),flx_wvl,neon_wvl,neon_fwhm);
    # Tt2 = super_resample((t_T_1km + t_t_1km).^2,wvl,neon_wvl,neon_fwhm);
    # t2Tso2 = super_resample(t_t_1km.^2 + 2 * t_t_1km .* sqrt(t_Tso_1km) + t_Tso_1km,wvl,neon_wvl,neon_fwhm);
    # tTt = super_resample((t_t_1km + t_T_1km) .* t_t_1km,wvl,neon_wvl,neon_fwhm);
    # tTT = super_resample((t_t_1km + t_T_1km) .* t_T_1km,wvl,neon_wvl,neon_fwhm);
    # tTso2 = super_resample((t_t_1km + sqrt(t_Tso_1km)).^2,wvl,neon_wvl,neon_fwhm);
    #
    # tTsot = super_resample((t_t_1km + sqrt(t_Tso_1km)) .* t_t_1km,wvl,neon_wvl,neon_fwhm);
    # tTsoTso = super_resample((t_t_1km + sqrt(t_Tso_1km)) .* sqrt(t_Tso_1km),wvl,neon_wvl,neon_fwhm);

    # SSIM Full Resolutions Calculations
    # T_1km_SSIM = super_resample(t_T_1km,wvl,ssim_wvl,ssim_zfwhm);
    # t_1km_SSIM = super_resample(t_t_1km,wvl,ssim_wvl,ssim_zfwhm);
    # Ts_1km_SSIM = super_resample(t_Ts_1km,wvl,ssim_wvl,ssim_zfwhm);
    # ts_1km_SSIM = super_resample(t_ts_1km,wvl,ssim_wvl,ssim_zfwhm);
    # Tso_1km_SSIM = super_resample(t_Tso_1km,wvl,ssim_wvl,ssim_zfwhm);
    # sph_alb_1km_SSIM = super_resample(t_sph_alb_1km,wvl,ssim_wvl,ssim_zfwhm);
    #
    # Tsts_Tt_SSIM = super_resample((t_Ts + t_ts).*(t_T+t_t),wvl,ssim_wvl,ssim_zfwhm);
    # Tt2_SSIM = super_resample((t_T_1km + t_t_1km).^2,wvl,ssim_wvl,ssim_zfwhm);
    # t2Tso2_SSIM = super_resample(t_t_1km.^2 + 2 * t_t_1km .* sqrt(t_Tso_1km) + t_Tso_1km,wvl,ssim_wvl,ssim_zfwhm);
    # tTt_SSIM = super_resample((t_t_1km + t_T_1km) .* t_t_1km,wvl,ssim_wvl,ssim_zfwhm);
    # tTT_SSIM = super_resample((t_t_1km + t_T_1km) .* t_T_1km,wvl,ssim_wvl,ssim_zfwhm);
    # tTso2_SSIM = super_resample((t_t_1km + sqrt(t_Tso_1km)).^2,wvl,ssim_wvl,ssim_zfwhm);
    #
    # tTsot_SSIM = super_resample((t_t_1km + sqrt(t_Tso_1km)) .* t_t_1km,wvl,ssim_wvl,ssim_zfwhm);
    # tTsoTso_SSIM = super_resample((t_t_1km + sqrt(t_Tso_1km)) .* sqrt(t_Tso_1km),wvl,ssim_wvl,ssim_zfwhm);

    ## Retrieve Surface Reflectance
    ##### 3# Tarp #############################################################
    # Match NIS & SSIR Time
    obs_time = np.reshape(nistime(tarp3_coord(2,1):tarp3_coord(2,2),tarp3_coord(1,1):tarp3_coord(1,2)),1,[]);
    n = np.size(obs_time);

    for j = 1:n:
        # Extract Pixel/Observation Radiance Spectrum
        Rob = tarp3_spc(j,1:152).T;
        rad0 = super_resample(r0,wvl_7sc,neon_wvl,neon_fwhm);
        # TOA Retrieval
        tarp3_ret_TOA[j,:] = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*((Ts + ts).*(T + t)));
        tarp3_ret_TOA_B[j,:] = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*Tsts_Tt);

        # Observation Retrieval
        [~,ssirtime_ind] = np.min(np.abs(ssirtime-obs_time[j]));

        Ifdn_obs_trim = super_resample(zspect(ssirtime_ind,:)',zwvl,neon_wvl,neon_fwhm);
        Ifup_obs_trim = super_resample(nspect(ssirtime_ind,:)',nwvl,neon_wvl,neon_fwhm);

        Ifup0 = super_resample(rho_a_I,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;
        rad0 = super_resample(rho_a_R,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;

        tarp3_ret_noadj(j,:) = pi .* (Rob - rad0) ./ (sph_alb_1km .* pi .* ...
            (Rob - rad0) + Ifdn_obs_trim .* (t_1km.^2 + 2 * t_1km .* sqrt(Tso_1km) + Tso_1km));

        rho_bar_ret = (Ifup_obs_trim - Ifup0) ./ ((Ifdn_obs_trim .* ...
            ((t_1km + sqrt(Tso_1km)).^2) + (Ifup_obs_trim - Ifup0) .* sph_alb_1km));
        tarp3_ret(j,:) = (Rob - rad0 - (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* t_1km) ./ ...
            pi .* (rho_bar_ret ./ (1 - rho_bar_ret .* sph_alb_1km))) .* ...
            (pi ./ (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* sqrt(Tso_1km))) .* (1 - rho_bar_ret .* sph_alb_1km);

        flight_alb(j,:) = Rob*pi./Ifdn_obs_trim;

    plt.figure;hold on
    # plot(neon_wvl,tarp3_ret_TOA,'Color',light_standard)
    # plot(neon_wvl,tarp3_ret_noadj,'Color',light_intermediate) # Intermediate
    # plot(neon_wvl,tarp3_ret,'Color',light_enhanced)
    plt.plot(asd_wvl,mean(tarp03,1)'./100,'k','LineWidth',1.5)
    # plot(neon_wvl,mean(tarp3_ret_TOA,1),'Color',standard,'LineWidth',1.5)
    plt.plot(neon_wvl,mean(flight_alb,1),'Color',[153,163,164]./255,'LineWidth',1.5)
    plt.plot(neon_wvl,mean(tarp3_ret_TOA_B,1),'-','Color',standard,'LineWidth',1.5)
    plt.plot(neon_wvl,mean(tarp3_ret_noadj,1),'Color',intermediate,'LineWidth',1.5) # Intermediate
    plt.plot(neon_wvl,mean(tarp3_ret,1),'Color',enhanced,'LineWidth',1.5)
    plt.ylabel('Reflectance','FontSize',12);xlabel('Wavelength [nm]','FontSize',12);
    plt.axis([350,1025,0,0.1])
    plt.set(gca,'FontSize',12,'XColor','k','YColor','k');box on;
    fig = plt.gcf;
    plt.set(fig,'PaperUnits','inches','PaperPosition',[0 0 4 4]);
    plt.savefig(fig,'-depsc','-r300','TestTarp3Retrieval.eps')

    wv_bands = xor(neon_wvl < 928, neon_wvl > 989);
    wv_bands = xor(wv_bands,neon_wvl > 1095);

    tarp03_r = super_resample(mean(tarp03,1),asd_wvl,neon_wvl,neon_fwhm)'./100;

    plt.figure;hold on;plot(neon_wvl,(mean(tarp3_ret,1)-tarp03_r)*100,'Color',enhanced)
    plt.plot(neon_wvl,(mean(tarp3_ret_noadj,1)-tarp03_r)*100,'Color',intermediate)
    plt.plot(neon_wvl,(mean(tarp3_ret_TOA_B,1)-tarp03_r)*100,'Color',standard);grid on;
    RMSE_tarp3 = [sqrt(sum(((mean(tarp3_ret,1)-tarp03_r)*100).^2)./numel(tarp03_r));
        sqrt(sum(((mean(tarp3_ret_noadj,1)-tarp03_r)*100).^2)./numel(tarp03_r));
        sqrt(sum(((mean(tarp3_ret_TOA_B,1)-tarp03_r)*100).^2)./numel(tarp03_r))];
    RMSE_tarp3_no_wv = [sqrt(sum(((mean(tarp3_ret(:,wv_bands),1)-tarp03_r(wv_bands))*100).^2)./numel(tarp03_r(wv_bands)));
        sqrt(sum(((mean(tarp3_ret_noadj(:,wv_bands),1)-tarp03_r(wv_bands))*100).^2)./numel(tarp03_r(wv_bands)));
        sqrt(sum(((mean(tarp3_ret_TOA_B(:,wv_bands),1)-tarp03_r(wv_bands))*100).^2)./numel(tarp03_r(wv_bands)))];
    plt.title('Tarp 3#')

    ##### 48# Tarp ############################################################
    # Match NIS & SSIR Time
    obs_time = nistime(tarp48_coord(2,1):tarp48_coord(2,2),tarp48_coord(1,1):tarp48_coord(1,2));
    n = numel(obs_time);
    for j = 1:n
        # Extract Pixel/Observation Radiance Spectrum
        Rob = tarp48_spc(j,1:152)';
        rad0 = super_resample(r0,wvl_7sc,neon_wvl,neon_fwhm);
        # TOA Retrieval
        tarp48_ret_TOA(j,:) = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*((Ts + ts).*(T + t)));
        tarp48_ret_TOA_B(j,:) = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*(Tsts_Tt));

        # Observation Retrieval
        [~,ssirtime_ind] = min(abs(ssirtime-obs_time(j)));

        Ifdn_obs_trim = super_resample(zspect(ssirtime_ind,:)',zwvl,neon_wvl,neon_fwhm);
        Ifup_obs_trim = super_resample(nspect(ssirtime_ind,:)',nwvl,neon_wvl,neon_fwhm);

        Ifup0 = super_resample(rho_a_I,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;
        rad0 = super_resample(rho_a_R,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;

        tarp48_ret_noadj(j,:) = pi .* (Rob - rad0) ./ (sph_alb_1km .* pi .* ...
            (Rob - rad0) + Ifdn_obs_trim .* (t_1km.^2 + 2 * t_1km .* sqrt(Tso_1km) + Tso_1km));

        rho_bar_ret = (Ifup_obs_trim - Ifup0) ./ ((Ifdn_obs_trim .* ...
            ((t_1km + sqrt(Tso_1km)).^2) + (Ifup_obs_trim - Ifup0) .* sph_alb_1km));
        tarp48_ret(j,:) = (Rob - rad0 - (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* t_1km) ./ ...
            pi .* (rho_bar_ret ./ (1 - rho_bar_ret .* sph_alb_1km))) .* ...
            (pi ./ (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* sqrt(Tso_1km))) .* (1 - rho_bar_ret .* sph_alb_1km);

        # Mu Correction
        tarp48_ret_noadj_mu(j,:) = pi .* (Rob - rad0) ./ (sph_alb_1km .* pi .* ...
            (Rob - rad0) + Ifdn_obs_trim .* (t_1km.^2 + (1 + mu) * t_1km .* sqrt(Tso_1km) + mu.*Tso_1km));

        rho_bar_ret_mu = (Ifup_obs_trim - Ifup0) ./ ((Ifdn_obs_trim .* ...
            (t_1km.^2 + (1 + mu) * t_1km .* sqrt(Tso_1km) + mu.*Tso_1km) + (Ifup_obs_trim - Ifup0) .* sph_alb_1km));
        tarp48_ret_mu(j,:) = (Rob - rad0 - (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* t_1km) ./ ...
            pi .* (rho_bar_ret_mu ./ (1 - rho_bar_ret_mu .* sph_alb_1km))) .* ...
            (pi ./ (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* mu .* sqrt(Tso_1km))) .* (1 - rho_bar_ret_mu .* sph_alb_1km);

        flight_alb(j,:) = Rob*pi./Ifdn_obs_trim;
    end
    figure;hold on
    # plot(neon_wvl,tarp48_ret_TOA,'Color',light_standard)
    # plot(neon_wvl,tarp48_ret_noadj,'Color',light_intermediate) # Intermediate
    # plot(neon_wvl,tarp48_ret,'Color',light_enhanced)
    plot(asd_wvl,mean(tarp48,1)'./100,'k','LineWidth',1.5)
    # plot(neon_wvl,mean(tarp48_ret_TOA,1),'Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(flight_alb,1),'Color',[153,163,164]./255,'LineWidth',1.5)
    plot(neon_wvl,mean(tarp48_ret_TOA_B,1),'-','Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(tarp48_ret_noadj,1),'Color',intermediate,'LineWidth',1.5) # Intermediate
    plot(neon_wvl,mean(tarp48_ret,1),'Color',enhanced,'LineWidth',1.5)
    ylabel('Reflectance','FontSize',12);xlabel('Wavelength [nm]','FontSize',12);
    axis([350,1025,0,1])
    set(gca,'FontSize',12,'XColor',[0,0,205]./255,'YColor',[0,0,205]./255);box on;
    fig = gcf;
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 4 4]);
    print(gcf,'-depsc','-r300','TestTarp48Retrieval.eps')
    # print(gcf,'-dpng','-r300','TestTarp48Retrieval.png')

    figure;hold on
    plot(asd_wvl,mean(tarp48,1)'./100,'k','LineWidth',1.5)
    plot(neon_wvl,mean(tarp48_ret_TOA_B,1),'-','Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(tarp48_ret_noadj_mu,1),'Color',intermediate,'LineWidth',1.5) # Intermediate
    plot(neon_wvl,mean(tarp48_ret_mu,1),'Color',enhanced,'LineWidth',1.5)
    ylabel('Reflectance','FontSize',12);xlabel('Wavelength [nm]','FontSize',12);
    axis([350,1025,0,1])
    set(gca,'FontSize',12,'XColor',[0,0,205]./255,'YColor',[0,0,205]./255);box on;
    fig = gcf;
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 4 4]);
    print(gcf,'-depsc','-r300','TestTarp48Retrieval_MU.eps')



    tarp48_r = super_resample(mean(tarp48,1),asd_wvl,neon_wvl,neon_fwhm)'./100;

    figure;hold on;plot(neon_wvl,(mean(tarp48_ret,1)-tarp48_r)*100,'Color',enhanced)
    plot(neon_wvl,(mean(tarp48_ret_noadj,1)-tarp48_r)*100,'Color',intermediate)
    plot(neon_wvl,(mean(tarp48_ret_TOA_B,1)-tarp48_r)*100,'Color',standard);grid on;
    RMSE_tarp48 = [sqrt(sum(((mean(tarp48_ret,1)-tarp48_r)*100).^2)./numel(tarp48_r));
        sqrt(sum(((mean(tarp48_ret_noadj,1)-tarp48_r)*100).^2)./numel(tarp48_r));
        sqrt(sum(((mean(tarp48_ret_TOA_B,1)-tarp48_r)*100).^2)./numel(tarp48_r))];
    RMSE_tarp48_no_wv = [sqrt(sum(((mean(tarp48_ret(:,wv_bands),1)-tarp48_r(wv_bands))*100).^2)./numel(tarp48_r(wv_bands)));
        sqrt(sum(((mean(tarp48_ret_noadj(:,wv_bands),1)-tarp48_r(wv_bands))*100).^2)./numel(tarp48_r(wv_bands)))
        sqrt(sum(((mean(tarp48_ret_TOA_B(:,wv_bands),1)-tarp48_r(wv_bands))*100).^2)./numel(tarp48_r(wv_bands)))];
    title('Tarp 48#')

    figure;hold on;
    plot(neon_wvl,mean(tarp48_ret,1),'Color',enhanced)
    ylabel('Reflectance');
    xlabel('Wavelength [nm]')
    title('Retrieved Reflectance in O_2 A-Band')
    axis([720,800,0.4,0.5]);grid on;box on
    set(gca,'FontSize',18)
    print('O2_BandReflectance.svg','-dsvg')


    # save('Emp_Fit_SSIM_Jun16_Tarp48.mat','neon_wvl','tarp48_ret','tarp48_ret_TOA')

    # Test "Empirical Smoothing"
    # "Best" fit NIS offset is 1.28
    # dlambda_NIS = -1.3:0.01:-1.2;#1.28;#1.2:0.01:1.4;
    # neon_wvl2 = repmat(neon_wvl,1,length(dlambda_NIS)) + repmat(dlambda_NIS,152,1);
    # dlambda_SSIM_zen = -3:0.1:3;
    # zwvl2 = repmat(zwvl,1,length(dlambda_SSIM_zen)) + repmat(dlambda_SSIM_zen,439,1);
    # for x = 1:length(dlambda_NIS)
    #     rad0 = super_resample(r0,flx_wvl,neon_wvl2(:,x),neon_fwhm);
    #     T_wvl = super_resample(t_T_1km,wvl,neon_wvl2(:,x),neon_fwhm);
    #     t_wvl = super_resample(t_t_1km,wvl,neon_wvl2(:,x),neon_fwhm);
    #     Ts_wvl = super_resample(t_Ts_1km,wvl,neon_wvl2(:,x),neon_fwhm);
    #     ts_wvl = super_resample(t_ts_1km,wvl,neon_wvl2(:,x),neon_fwhm);
    #     Tso_wvl = super_resample(t_Tso_1km,wvl,neon_wvl2(:,x),neon_fwhm);
    #     sph_alb_layer = super_resample(t_sph_alb_1km,wvl,neon_wvl2(:,x),neon_fwhm);
    #     Tsts_Tt = super_resample((t_Ts + t_ts).*(t_T+t_t),wvl,neon_wvl2(:,x),neon_fwhm);
    #
    #     tarp48_ret_TOA_WVLS(:,x) = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*(Tsts_Tt));
    #
    #     for y = 1:length(dlambda_SSIM_zen)
    #         Ifdn_obs_trim = super_resample(zspect(ssirtime_ind,:)',zwvl2(:,y),neon_wvl2(:,x),neon_fwhm);
    #         Ifup_obs_trim = super_resample(nspect(ssirtime_ind,:)',nwvl,neon_wvl2(:,x),neon_fwhm);
    #         Ifup0 = super_resample(rho_a_I,flx_wvl,neon_wvl2(:,x),neon_fwhm).*Ifdn_obs_trim;
    #
    #         tarp48_ret_noadj_WVLS(:,x,y) = pi .* (Rob - rad0) ./ (sph_alb_layer .* pi .* ...
    #             (Rob - rad0) + Ifdn_obs_trim .* (T_wvl + t_wvl).^2);
    #         tarp48_ret_noadj_WVLS2(:,x,y) = pi .* (Rob - rad0) ./ (sph_alb_layer .* pi .* ...
    #             (Rob - rad0) + Ifdn_obs_trim .* (t_wvl.^2 + 2 * t_wvl .* sqrt(Tso_wvl) + Tso_wvl));
    #     end
    # end



    # dlambda_SSIM_zen = -5:0.1:3;
    # zwvl2 = repmat(zwvl,1,length(dlambda_SSIM_zen)) + repmat(dlambda_SSIM_zen,439,1);
    # for x = 1:length(dlambda_SSIM_zen)
    #     Ifdn_obs_trim = super_resample(zspect(ssirtime_ind,:)',zwvl2(:,y),neon_wvl,neon_fwhm);
    #     tarp48_ret_TEST(:,x) = pi.*(Rob-rad0)./(sph_alb_layer.*pi.*(Rob-rad0) + Ifdn_obs_trim.*(Tt2_d(1:152)));
    # end

    # dlambda_SSIM_nad = -3:0.1:3;
    # # dfwhm = -3:0.1:3;
    # neon_wvl2 = repmat(neon_wvl,1,length(dlambda_NIS)) + repmat(dlambda_NIS,152,1);
    # nwvl2 = repmat(nwvl,1,length(dlambda_SSIM_nad)) + repmat(dlambda_SSIM_nad,256,1);
    # zwvl2 = repmat(zwvl,1,length(dlambda_SSIM_zen)) + repmat(dlambda_SSIM_zen,439,1);
    # neon_fwhm2 = repmat(neon_fwhm,1,length(dfwhm)) + repmat(dfwhm,152,1);
    # Rob = mean(tarp48_spc(:,1:152),1)';
    # for x = 1:length(dlambda_SSIM_zen)
    # #     rad_resample = super_resample(data_rad(:,2)./100,data_rad(:,1),neon_wvl2(:,x),neon_fwhm);
    # #     sph_alb_resample(:,1) = super_resample(sph_alb_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm);
    # #     sph_alb_resample(:,2) = super_resample(sph_alb_orig(:,2),wvl,neon_wvl2(:,x),neon_fwhm);
    # #     Ts_resample = super_resample(Ts_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm);
    # #     ts_resample = super_resample(ts_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm);
    # #     T_resample = super_resample(T_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm);
    # #     t_resample = super_resample(t_orig(:,1),wvl,neon_wvl2(:,x),neon_fwhm);
    # #     sph_alb_layer_resample = (sph_alb_resample(:,1)-sph_alb_resample(:,2).*(t_resample(:,1)+T_resample(:,1)))./(1-sph_alb_resample(:,2).*(t_resample(:,1)+T_resample(:,1)));
    #     for y = 1:length(dlambda_SSIM_nad)
    #         Ifdn_obs_trim_NIS_SSIM = super_resample(Ifdn_obs',zwvl2(:,x),neon_wvl,neon_fwhm);
    #         Ifup_obs_trim_NIS_SSIM = super_resample(Ifup_obs',nwvl2(:,y),neon_wvl,neon_fwhm);
    #         Ifup0_resample = rho_a_I.*Ifdn_obs_trim_NIS_SSIM;
    #         rad0_resample = rho_a_R.*Ifdn_obs_trim_NIS_SSIM;
    #         rho_bar_ret_NIS_SSIM = (Ifup_obs_trim_NIS_SSIM-Ifup0_resample)./((Ifdn_obs_trim_NIS_SSIM.*(T(:,1)+t(:,1)).^2)+(Ifup_obs_trim_NIS_SSIM-Ifup0_resample).*sph_alb_layer);
    #
    #         tarp48_ret_SSIM_TEST(:,x,y) = (Rob - rad0_resample - (Ifdn_obs_trim_NIS_SSIM.*(T(:,1)+t(:,1)).*t(:,1))./pi.*(rho_bar_ret_NIS_SSIM./(1-rho_bar_ret_NIS_SSIM.*sph_alb_layer))).*(pi./(Ifdn_obs_trim_NIS_SSIM.*(T(:,1)+t(:,1)).*T(:,1))).*(1-rho_bar_ret_NIS_SSIM.*sph_alb_layer);
    #
    #     end
    # end


    ##### Vegetation ##########################################################
    # Match NIS & SSIR Time
    obs_time = nistime(veg_coord(2,1):veg_coord(2,2),veg_coord(1,1):veg_coord(1,2));
    n = numel(obs_time);
    for j = 1:n
        # Extract Pixel/Observation Radiance Spectrum
        Rob = veg_spc(j,1:152)';
        rad0 = super_resample(r0,wvl_7sc,neon_wvl,neon_fwhm);
        # TOA Retrieval
        veg_ret_TOA(j,:) = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*((Ts + ts).*(T + t)));
        veg_ret_TOA_B(j,:) = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*(Tsts_Tt));

        # Observation Retrieval
        [~,ssirtime_ind] = min(abs(ssirtime-obs_time(j)));

        Ifdn_obs_trim = super_resample(zspect(ssirtime_ind,:)',zwvl,neon_wvl,neon_fwhm);
        Ifup_obs_trim = super_resample(nspect(ssirtime_ind,:)',nwvl,neon_wvl,neon_fwhm);

        Ifup0 = super_resample(rho_a_I,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;
        rad0 = super_resample(rho_a_R,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;

        veg_ret_noadj(j,:) = pi .* (Rob - rad0) ./ (sph_alb_1km .* pi .* ...
            (Rob - rad0) + Ifdn_obs_trim .* (t_1km.^2 + 2 * t_1km .* sqrt(Tso_1km) + Tso_1km));

        rho_bar_ret = (Ifup_obs_trim - Ifup0) ./ ((Ifdn_obs_trim .* ...
            ((t_1km + sqrt(Tso_1km)).^2) + (Ifup_obs_trim - Ifup0) .* sph_alb_1km));
        veg_ret(j,:) = (Rob - rad0 - (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* t_1km) ./ ...
            pi .* (rho_bar_ret ./ (1 - rho_bar_ret .* sph_alb_1km))) .* ...
            (pi ./ (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* sqrt(Tso_1km))) .* (1 - rho_bar_ret .* sph_alb_1km);

        # TestTarp3Retrieval

        # Mu Correction
        veg_ret_noadj_mu(j,:) = pi .* (Rob - rad0) ./ (sph_alb_1km .* pi .* ...
            (Rob - rad0) + Ifdn_obs_trim .* (t_1km.^2 + (1 + mu) * t_1km .* sqrt(Tso_1km) + mu.*Tso_1km));

        rho_bar_ret_mu = (Ifup_obs_trim - Ifup0) ./ ((Ifdn_obs_trim .* ...
            (t_1km.^2 + (1 + mu) * t_1km .* sqrt(Tso_1km) + mu.*Tso_1km) + (Ifup_obs_trim - Ifup0) .* sph_alb_1km));
        veg_ret_mu(j,:) = (Rob - rad0 - (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* t_1km) ./ ...
            pi .* (rho_bar_ret_mu ./ (1 - rho_bar_ret_mu .* sph_alb_1km))) .* ...
            (pi ./ (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* mu .* sqrt(Tso_1km))) .* (1 - rho_bar_ret_mu .* sph_alb_1km);

        flight_alb(j,:) = Rob*pi./Ifdn_obs_trim;
    end
    figure;hold on
    # plot(neon_wvl,veg_ret_TOA,'Color',light_standard)
    # plot(neon_wvl,veg_ret_noadj,'Color',light_intermediate) # Intermediate
    # plot(neon_wvl,veg_ret,'Color',light_enhanced)
    plot(asd_wvl,mean(veg,1)'./100,'k','LineWidth',1.5)
    # plot(neon_wvl,mean(veg_ret_TOA,1),'Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(flight_alb,1),'Color',[153,163,164]./255,'LineWidth',1.5)
    plot(neon_wvl,mean(veg_ret_TOA_B,1),'-','Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(veg_ret_noadj,1),'Color',intermediate,'LineWidth',1.5) # Intermediate
    plot(neon_wvl,mean(veg_ret,1),'Color',enhanced,'LineWidth',1.5)
    ylabel('Reflectance','FontSize',12);xlabel('Wavelength [nm]','FontSize',12);
    axis([350,1025,0,0.4])
    set(gca,'FontSize',12,'XColor',[34,139,34]./255,'YColor',[34,139,34]./255);box on
    # title({'June 8th Vegetated Surface Retrieval','with Solar Zenith Angle Correction'},'Color','k')
    fig = gcf;
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 4 4]);
    print(gcf,'-depsc','-r300','TestVegRetrieval.eps')

    figure;hold on
    plot(asd_wvl,mean(veg,1)'./100,'k','LineWidth',1.5)
    plot(neon_wvl,mean(veg_ret_TOA_B,1),'-','Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(veg_ret_noadj_mu,1),'Color',intermediate,'LineWidth',1.5) # Intermediate
    plot(neon_wvl,mean(veg_ret_mu,1),'Color',enhanced,'LineWidth',1.5)
    ylabel('Reflectance','FontSize',12);xlabel('Wavelength [nm]','FontSize',12);
    axis([350,1025,0,0.4])
    set(gca,'FontSize',12,'XColor',[34,139,34]./255,'YColor',[34,139,34]./255);box on
    title({'June 8th Vegetated Surface Retrieval','with Solar Zenith Angle Correction'},'Color','k')
    fig = gcf;
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
    print(gcf,'-depsc','-r300','TestVegRetrieval_MU.eps')

    veg_r = super_resample(mean(veg,1),asd_wvl,neon_wvl,neon_fwhm)'./100;

    disp('NoAdj Mu:');disp(sqrt(sum(((mean(veg_ret_noadj_mu,1)-veg_r)*100).^2)./numel(veg_r)))
    disp('Ret Mu:');disp(sqrt(sum(((mean(veg_ret_mu,1)-veg_r)*100).^2)./numel(veg_r)))

    figure;hold on;plot(neon_wvl,(mean(veg_ret,1)-veg_r)*100,'Color',enhanced)
    plot(neon_wvl,(mean(veg_ret_noadj,1)-veg_r)*100,'Color',intermediate)
    plot(neon_wvl,(mean(veg_ret_TOA_B,1)-veg_r)*100,'Color',standard);grid on;
    RMSE_veg = [sqrt(sum(((mean(veg_ret,1)-veg_r)*100).^2)./numel(veg_r));
        sqrt(sum(((mean(veg_ret_noadj,1)-veg_r)*100).^2)./numel(veg_r));
        sqrt(sum(((mean(veg_ret_TOA_B,1)-veg_r)*100).^2)./numel(veg_r))];
    RMSE_veg_no_wv = [sqrt(sum(((mean(veg_ret(:,wv_bands),1)-veg_r(wv_bands))*100).^2)./numel(veg_r(wv_bands)));
        sqrt(sum(((mean(veg_ret_noadj(:,wv_bands),1)-veg_r(wv_bands))*100).^2)./numel(veg_r(wv_bands)));
        sqrt(sum(((mean(veg_ret_TOA_B(:,wv_bands),1)-veg_r(wv_bands))*100).^2)./numel(veg_r(wv_bands)))];
    title('Vegetation')

    ##### EW Road #############################################################
    # Match NIS & SSIR Time
    obs_time = nistime(EWroad_coord(2,1):EWroad_coord(2,2),EWroad_coord(1,1):EWroad_coord(1,2));
    n = numel(obs_time);
    for j = 1:n
        # Extract Pixel/Observation Radiance Spectrum
        Rob = EWroad_spc(j,1:152)';
        rad0 = super_resample(r0,wvl_7sc,neon_wvl,neon_fwhm);
        # TOA Retrieval
        EWroad_ret_TOA(j,:) = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*((Ts + ts).*(T + t)));
        EWroad_ret_TOA_B(j,:) = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*(Tsts_Tt));

        # Observation Retrieval
        [~,ssirtime_ind] = min(abs(ssirtime-obs_time(j)));

        Ifdn_obs_trim = super_resample(zspect(ssirtime_ind,:)',zwvl,neon_wvl,neon_fwhm);
        Ifup_obs_trim = super_resample(nspect(ssirtime_ind,:)',nwvl,neon_wvl,neon_fwhm);

        Ifup0 = super_resample(rho_a_I,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;
        rad0 = super_resample(rho_a_R,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;

        EWroad_ret_noadj(j,:) = pi .* (Rob - rad0) ./ (sph_alb_1km .* pi .* ...
            (Rob - rad0) + Ifdn_obs_trim .* (t_1km.^2 + 2 * t_1km .* sqrt(Tso_1km) + Tso_1km));

        rho_bar_ret = (Ifup_obs_trim - Ifup0) ./ ((Ifdn_obs_trim .* ...
            ((t_1km + sqrt(Tso_1km)).^2) + (Ifup_obs_trim - Ifup0) .* sph_alb_1km));
        EWroad_ret(j,:) = (Rob - rad0 - (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* t_1km) ./ ...
            pi .* (rho_bar_ret ./ (1 - rho_bar_ret .* sph_alb_1km))) .* ...
            (pi ./ (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* sqrt(Tso_1km))) .* (1 - rho_bar_ret .* sph_alb_1km);

        flight_alb(j,:) = Rob*pi./Ifdn_obs_trim;
    end
    figure;hold on
    # plot(neon_wvl,EWroad_ret_TOA,'Color',light_standard)
    # plot(neon_wvl,EWroad_ret_noadj,'Color',light_intermediate) # Intermediate
    # plot(neon_wvl,EWroad_ret,'Color',light_enhanced)
    plot(asd_wvl,mean(ewroad,1)'./100,'k','LineWidth',1.5)
    # plot(neon_wvl,mean(EWroad_ret_TOA,1),'Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(flight_alb,1),'Color',[153,163,164]./255,'LineWidth',1.5)
    plot(neon_wvl,mean(EWroad_ret_TOA_B,1),'-','Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(EWroad_ret_noadj,1),'Color',intermediate,'LineWidth',1.5) # Intermediate
    plot(neon_wvl,mean(EWroad_ret,1),'Color',enhanced,'LineWidth',1.5)
    ylabel('Reflectance','FontSize',12);xlabel('Wavelength [nm]','FontSize',12);
    axis([350,1025,0,0.60])
    set(gca,'FontSize',12,'XColor',[156,102,31]./255,'YColor',[156,102,31]./255);box on
    fig = gcf;
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 4 4]);
    print(gcf,'-depsc','-r300','TestEWroadRetrieval.eps')

    ewroad_r = super_resample(mean(ewroad,1),asd_wvl,neon_wvl,neon_fwhm)'./100;

    figure;hold on;plot(neon_wvl,(mean(EWroad_ret,1)-ewroad_r)*100,'Color',enhanced)
    plot(neon_wvl,(mean(EWroad_ret_noadj,1)-ewroad_r)*100,'Color',intermediate)
    plot(neon_wvl,(mean(EWroad_ret_TOA_B,1)-ewroad_r)*100,'Color',standard);grid on
    RMSE_EWroad = [sqrt(sum(((mean(EWroad_ret,1)-ewroad_r)*100).^2)./numel(ewroad_r));
        sqrt(sum(((mean(EWroad_ret_noadj,1)-ewroad_r)*100).^2)./numel(ewroad_r));
        sqrt(sum(((mean(EWroad_ret_TOA_B,1)-ewroad_r)*100).^2)./numel(ewroad_r))];
    RMSE_EWroad_no_wv = [sqrt(sum(((mean(EWroad_ret(:,wv_bands),1)-ewroad_r(wv_bands))*100).^2)./numel(ewroad_r(wv_bands)));
        sqrt(sum(((mean(EWroad_ret_noadj(:,wv_bands),1)-ewroad_r(wv_bands))*100).^2)./numel(ewroad_r(wv_bands)));
        sqrt(sum(((mean(EWroad_ret_TOA_B(:,wv_bands),1)-ewroad_r(wv_bands))*100).^2)./numel(ewroad_r(wv_bands)))];
    title('EW Road')

    ##### NS Road #############################################################
    # Match NIS & SSIR Time
    obs_time = nistime(NSroad_coord(2,1):NSroad_coord(2,2),NSroad_coord(1,1):NSroad_coord(1,2));
    n = numel(obs_time);
    for j = 1:n
        # Extract Pixel/Observation Radiance Spectrum
        Rob = NSroad_spc(j,1:152)';
        rad0 = super_resample(r0,wvl_7sc,neon_wvl,neon_fwhm);
        # TOA Retrieval
        NSroad_ret_TOA(j,:) = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*((Ts + ts).*(T + t)));
        NSroad_ret_TOA_B(j,:) = pi.*(Rob-rad0)./(sph_alb.*pi.*(Rob-rad0)+mu.*I0.*(Tsts_Tt));

        # Observation Retrieval
        [~,ssirtime_ind] = min(abs(ssirtime-obs_time(j)));

        Ifdn_obs_trim = super_resample(zspect(ssirtime_ind,:)',zwvl,neon_wvl,neon_fwhm);
        Ifup_obs_trim = super_resample(nspect(ssirtime_ind,:)',nwvl,neon_wvl,neon_fwhm);

        Ifup0 = super_resample(rho_a_I,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;
        rad0 = super_resample(rho_a_R,flx_wvl,neon_wvl,neon_fwhm).*Ifdn_obs_trim;

        NSroad_ret_noadj(j,:) = pi .* (Rob - rad0) ./ (sph_alb_1km .* pi .* ...
            (Rob - rad0) + Ifdn_obs_trim .* (t_1km.^2 + 2 * t_1km .* sqrt(Tso_1km) + Tso_1km));

        rho_bar_ret = (Ifup_obs_trim - Ifup0) ./ ((Ifdn_obs_trim .* ...
            ((t_1km + sqrt(Tso_1km)).^2) + (Ifup_obs_trim - Ifup0) .* sph_alb_1km));
        NSroad_ret(j,:) = (Rob - rad0 - (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* t_1km) ./ ...
            pi .* (rho_bar_ret ./ (1 - rho_bar_ret .* sph_alb_1km))) .* ...
            (pi ./ (Ifdn_obs_trim .* (t_1km + sqrt(Tso_1km)) .* sqrt(Tso_1km))) .* (1 - rho_bar_ret .* sph_alb_1km);

        flight_alb(j,:) = Rob*pi./Ifdn_obs_trim;
    end
    figure;hold on
    # plot(neon_wvl,NSroad_ret_TOA,'Color',light_standard)
    # plot(neon_wvl,NSroad_ret_noadj,'Color',light_intermediate) # Intermediate
    # plot(neon_wvl,NSroad_ret,'Color',light_enhanced)
    plot(asd_wvl,mean(nsroad,1)'./100,'k','LineWidth',1.5)
    # plot(neon_wvl,mean(NSroad_ret_TOA,1),'Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(flight_alb,1),'Color',[153,163,164]./255,'LineWidth',1.5)
    plot(neon_wvl,mean(NSroad_ret_TOA_B,1),'-','Color',standard,'LineWidth',1.5)
    plot(neon_wvl,mean(NSroad_ret_noadj,1),'Color',intermediate,'LineWidth',1.5) # Intermediate
    plot(neon_wvl,mean(NSroad_ret,1),'Color',enhanced,'LineWidth',1.5)
    ylabel('Reflectance','FontSize',12);xlabel('Wavelength [nm]','FontSize',12);
    axis([350,1025,0,0.6])
    set(gca,'FontSize',12,'XColor',[156,102,31]./255,'YColor',[156,102,31]./255);box on
    fig = gcf;
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 4 4]);
    print(gcf,'-depsc','-r300','TestNSroadRetrieval.eps')

    nsroad_r = super_resample(mean(nsroad,1),asd_wvl,neon_wvl,neon_fwhm)'./100;

    figure;hold on;plot(neon_wvl,(mean(NSroad_ret,1)-nsroad_r)*100,'Color',enhanced)
    plot(neon_wvl,(mean(NSroad_ret_noadj,1)-nsroad_r)*100,'Color',intermediate)
    plot(neon_wvl,(mean(NSroad_ret_TOA_B,1)-nsroad_r)*100,'Color',standard);grid on;
    RMSE_NSroad = [sqrt(sum(((mean(NSroad_ret,1)-nsroad_r)*100).^2)./numel(nsroad_r));
        sqrt(sum(((mean(NSroad_ret_noadj,1)-nsroad_r)*100).^2)./numel(nsroad_r));
        sqrt(sum(((mean(NSroad_ret_TOA_B,1)-nsroad_r)*100).^2)./numel(nsroad_r))];
    RMSE_NSroad_no_wv = [sqrt(sum(((mean(NSroad_ret(:,wv_bands),1)-nsroad_r(wv_bands))*100).^2)./numel(nsroad_r(wv_bands)));
        sqrt(sum(((mean(NSroad_ret_noadj(:,wv_bands),1)-nsroad_r(wv_bands))*100).^2)./numel(nsroad_r(wv_bands)));
        sqrt(sum(((mean(NSroad_ret_TOA_B(:,wv_bands),1)-nsroad_r(wv_bands))*100).^2)./numel(nsroad_r(wv_bands)))];
    title('NS Road')

    # Write file of RMS Error
    fid = fopen('RMSE_Output_No_WV.txt','w');
    fprintf(fid,'\t #s \t #s \t #s \n #s \t #4.2f \t #4.2f \t #4.2f \n #s \t #4.2f \t #4.2f \t #4.2f \n #s \t #4.2f \t #4.2f \t #4.2f \n #s \t #4.2f \t #4.2f \t #4.2f \n #s \t #4.2f \t #4.2f \t #4.2f',...
        'Standard','Intermediate','Enhanced',...
        'Tarp 3#',RMSE_tarp3_no_wv(3),RMSE_tarp3_no_wv(2),RMSE_tarp3_no_wv(1),...
        'Tarp 48#',RMSE_tarp48_no_wv(3),RMSE_tarp48_no_wv(2),RMSE_tarp48_no_wv(1),...
        'Vegetation',RMSE_veg_no_wv(3),RMSE_veg_no_wv(2),RMSE_veg_no_wv(1),...
        'EW Road',RMSE_EWroad_no_wv(3),RMSE_EWroad_no_wv(2),RMSE_EWroad_no_wv(1),...
        'NS Road',RMSE_NSroad_no_wv(3),RMSE_NSroad_no_wv(2),RMSE_NSroad_no_wv(1));
    fclose(fid);

    fid = fopen('RMSE_Output.txt','w');
    fprintf(fid,'\t #s \t #s \t #s \n #s \t #4.2f \t #4.2f \t #4.2f \n #s \t #4.2f \t #4.2f \t #4.2f \n #s \t #4.2f \t #4.2f \t #4.2f \n #s \t #4.2f \t #4.2f \t #4.2f \n #s \t #4.2f \t #4.2f \t #4.2f',...
        'Standard','Intermediate','Enhanced',...
        'Tarp 3#',RMSE_tarp3(3),RMSE_tarp3(2),RMSE_tarp3(1),...
        'Tarp 48#',RMSE_tarp48(3),RMSE_tarp48(2),RMSE_tarp48(1),...
        'Vegetation',RMSE_veg(3),RMSE_veg(2),RMSE_veg(1),...
        'EW Road',RMSE_EWroad(3),RMSE_EWroad(2),RMSE_EWroad(1),...
        'NS Road',RMSE_NSroad(3),RMSE_NSroad(2),RMSE_NSroad(1));
    fclose(fid);

    # Names = {'Tarp 3#';'Tarp 48#';'Vegetation';'EW Road';'NS Road'};
    # Standard_RMSE = [RMSE_tarp3_no_wv(2);RMSE_tarp48_no_wv(2);RMSE_veg_no_wv(2);RMSE_EWroad_no_wv(2);RMSE_NSroad_no_wv(2)];
    # Enhanced_RMSE = [RMSE_tarp3_no_wv(1);RMSE_tarp48_no_wv(1);RMSE_veg_no_wv(1);RMSE_EWroad_no_wv(1);RMSE_NSroad_no_wv(1)];
    # T = table(Standard_RMSE,Enhanced_RMSE,'RowNames',Names);
    # writetable(T,'RMSE_Output_No_WV.txt')

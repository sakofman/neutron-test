"""
Module Name: exp_data_module

Description
-----------
Provides ExpData class for storing and manipulating the data from the FITS datacube.

Date
----
06/12/2024

Corrections
-----------
    (future) Handle FITS processing so imageJ is not needed. Should support manual selection for regions of image.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from skimage import exposure

class ExpData:
    # """
    # A class to represent the experimentally obtained data from sample under invesigation
    # """
    # def __init__(self, datacube_filepath):
    #     """
    #     Constructor for an ExpData instance.
        
    #     Arguments
    #     ---------
    #     datacube_filepath : str
    #         Path to FITS datacube.
    #     """
        # set attr sing fio again as in beam class

    # 
    def ReadExpTransmission(FileName) :
        # read experimental data
        Measured = pd.read_table(FileName)
        # create numpy arrays
        T = Measured[Measured.columns[0]].to_numpy()*1e6  # convert to us
        Tr = Measured[Measured.columns[1]].to_numpy()
        return T, Tr
    
    # TODO: pasted from prev code. adjust/rewrite to fit with new code
    # -------------------------------------------------------------------------------
    # Calculate array of Tbn values from the experimental TOF array.
    # The left boundary is saved in Pixelman code, Tbin starts from TOF value to the right (adding time)
    # Input: Texp - TOF array from experimental data file
    # Returns: Tbin Array
    # -------------------------------------------------------------------------------
    def CalcTbinArray(Texp):
        # Calculate Tbin array automatically
        # need to take care of the readout gaps
        TexpDif1 = np.diff(Texp)
        TexpDif2 = np.diff(TexpDif1)

        # check that there are not too many jumps in Tbin values (due to gaps), otherwise something is not right here
        if len(np.where(TexpDif2>TexpDif1[1:]*0.05)) > 100 : # compared to *0.05 - Spectra.txt file has limited number of digits and often Tbins vary slightly, e.g. 10.24 us binning at 10 ms has 10000.2, then 10001.3, etc
            raise Exception("Cannot calculate Tbn array, something is wrong with TOF array!!!!\n Too many large jumps in TOF array\n Exit here")
        print('Found gaps in TOF array at indices',np.where(TexpDif2>TexpDif1[1:]*0.05))

        Tbin = np.zeros(len(Texp))
        # first cell is set to second Diff cell
        Tbin[0] = TexpDif1[0]
        # rest of Tbins are the differences from the cell to the previos one
        Tbin[1:] = TexpDif1
        # now correct those large jumps in TexpDiff1 (e.g. gaps for readout) - set them to proper values
        for i in np.where(TexpDif2>TexpDif1[1:]*0.05):
            Tbin[i+2] = TexpDif1[i+2]

        #plt.plot(Texp, Tbin, label='Fitted', color='m')
        #plt.xlim([0,100])
        #plt.ylim([0,0.1])

        return Tbin

    # -------------------------------------------------------------------------------
    # Calculate theoretical transmission for Experimental TOF values
    # -----INPUT
    # Par - paremeters of the model
    # Texp - time array of experimental data
    # TrExp - measured transmission values
    # TrConvld - Theoretical transmission convolved with beam profile
    # NsubCells - how many subcells to use in averagin of theoretical transm for TOF bin in experiment
    # Plot - optional, whether to plot the graph or not
    # -----RETURNS
    # square of differences between Exp and Theoretical values
    # TODO: pasted from prev code. adjust/rewrite to fit with new code
    # -------------------------------------------------------------------------------
    def CalcTransmForExpPoints(Par, Texp, TrExp, TrConvld,  Tmax, NsubCells, Plot=False):
        #Tbin = Par['Time bin']
        Tbin = CalcTbinArray(Texp)
        dT = Par['Trigger delay']

        # select only exp. points which are above the Emin(Tmax) value and below Emax
        # T max was reduced by convolution from corresponding Emin value!!!
        Tmin   = Get_TOF_FromE(Par['Maximum E'], Par['Flight path'])
        #Tmax   = Get_TOF_FromE(Par['Minimum E'], Par['Flight path'])
        Tmin1 = max(Tmin,(Tmin - dT))
        Tmax1 = min(Tmax,(Tmax - Tbin[-1] - dT)) # use last Tbin here, the longest TOF to be used in that estimate
        #print('Emin',Get_E_FromTOF(Tmax1, Par['Flight path']), 'Emax=', Get_E_FromTOF(Tmin1, Par['Flight path']))
        Texp1  = Texp  [np.where(((Texp > Tmin1) & (Texp < Tmax1)))]
        TrExp1 = TrExp [np.where(((Texp > Tmin1) & (Texp < Tmax1)))]

        # calculate the values of Theoretical transmission averaged over the time bin used in experiment
        TrTheor = np.zeros(len(Texp1))
        # set TtempArr to first Tbin value, then will be chekcing if it needs to be recalculated
        TbinCurrent = Tbin[0]
        TtempArr = np.arange(NsubCells) * Tbin[0]/ NsubCells
        for i in range(len(Texp1)):
            #check if Tbin has changed, need to recalculate TtempArr
            if TbinCurrent != Tbin[i] :
                TbinCurrent = Tbin[i]
                TtempArr = np.arange(NsubCells) * Tbin[i] / NsubCells
            #Calculate transmission in all points in TtempArr shifted by dT and Texp[i]
            TrTempArr = TrConvld(TtempArr + Texp1[i] + dT)
            TrTheor[i] = np.average(TrTempArr)

        ChiSq = np.sum((TrTheor - TrExp1) ** 2)
        #ChiSq = -np.sum( (TrTheor * TrExp1) )
        print('Chi square=', ChiSq,'\n**********************************\n')

        if Plot:
            plt.figure(figsize=(14,3.5))
            plt.subplot(121)
            plt.scatter(Texp1, TrExp1, label='Measured', color='g', marker='.')
            plt.plot(Texp1, TrTheor, label='Fitted', color='m')
            plt.legend()
            plt.grid(True)
            plt.xlabel("Time (us)")
            plt.xscale('log')
            plt.ylabel("Transmission")
            # plt.xlim([50,200])
            #--------------------------------
            plt.subplot(122)
            Earr = Get_E_FromTOF(Texp1+dT, Par['Flight path'])
            plt.scatter(Earr, TrExp1, label='Measured', color='g', marker='.')
            plt.plot(Earr, TrTheor, label='Fitted', color='m')
            plt.legend()
            plt.grid(True)
            plt.xlabel("Energy (eV)")
            plt.xscale('log')
            plt.ylabel("Transmission")

        return ChiSq


    # def store_spectral_data(self):
    #     """ 
    #     Stores the TOF for the spectral data in a csv file.

    #     Returns
    #     -------
    #     str
    #         Path to created csv file.
    #     """
    #     # TODO:

    # def extract_spectral_data(self):
    #     """
    #     Prompts user to select a section of the image (using matplotlib) and stores spectral data for region in a csv file.

    #     Returns
    #     -------
    #     str
    #         Path to created csv file.
    #     """
    #     # TODO:

    # # FITS processing code from MCP Imaging Automation Script
    # # TODO: adjust to fit our goals, delete prints and plots
    # def process_file(imageFile):
    #     # image_data
    #     hdu_list = fits.open(imageFile)
    #     hdu_list.info()
    #     image_data = hdu_list[0].data
    #     print(type(image_data))
    #     print(image_data.shape)
    #     hdu_list.close()

    #     # some vals
    #     print('Min:', np.min(image_data))
    #     print('Max:', np.max(image_data))
    #     print('Mean:', np.mean(image_data))

    #     # equalize img data
    #     image_data_equalized = exposure.equalize_hist(image_data)
    #     image_data_gamma_adjusted = exposure.adjust_gamma(image_data_equalized) #idk if this helped it

    #     plt.figure(figsize=(15, 5))

    #     plt.subplot(1, 3, 1)
    #     plt.imshow(image_data, cmap='gray')
    #     plt.title('Original Image')
    #     plt.colorbar()

    #     plt.subplot(1, 3, 2)
    #     plt.imshow(image_data_equalized, cmap='gray')
    #     plt.title('Histogram Equalized Image')
    #     plt.colorbar()

    #     plt.subplot(1, 3, 3)
    #     plt.imshow(image_data_gamma_adjusted, cmap='gray')
    #     plt.title('Gamma Adjusted Image')
    #     plt.colorbar()

    #     plt.show()
    #     image_data_log = np.log1p(image_data)

    #     plt.imshow(image_data_log, cmap='gray')
    #     plt.colorbar()
    #     plt.show()

    #     histogram = plt.hist(image_data_equalized.flatten(), bins=800)
    #     return
    


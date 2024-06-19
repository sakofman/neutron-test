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
    Fix TODOs
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from skimage import exposure

class ExpData:
    """
    A class to represent the experimentally obtained data from sample under invesigation
    """
    def __init__(self, datacube_filepath):
        """
        Constructor for an ExpData instance.
        
        Arguments
        ---------
        datacube_filepath : str
            Path to FITS datacube.
        """
        # set attr sing fio again as in beam class

    def store_spectral_data(self):
        """ 
        Stores the TOF for the spectral data in a csv file.

        Returns
        -------
        str
            Path to created csv file.
        """
        # TODO:

    def extract_spectral_data(self):
        """
        Prompts user to select a section of the image (using matplotlib) and stores spectral data for region in a csv file.

        Returns
        -------
        str
            Path to created csv file.
        """
        # TODO:

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
    
    # TODO:
    def file_to_arr(filename):
        """ 
        Read measured transmission data from a TXT or CSV file into pandas frame
        then convert it to numpy arrays
        Returns Time and Tr arrays for the experimental data
        """
        # read experimental data
        Measured = pd.read_table(filename)
        # create numpy arrays
        T = Measured[Measured.columns[0]].to_numpy()*1e6  # convert to us
        Tr = Measured[Measured.columns[1]].to_numpy()
        return T, Tr


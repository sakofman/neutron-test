"""
Module Name: process_module

Description
-----------
- Load locally stored data such as constants, cross section data, and isotope data.
- Store list of isotopes to be fitted to a particular sample (should be a user input -- perhaps checkbox)
- Import locally stored isotope data into a pandas dataframe (filtered to only include the isotopes that are needed)

Date
----
06/18/2024

Corrections
-----------
"""
import pandas as pd
import numpy as np
import scipy.interpolate as intrp
import matplotlib.pyplot as plt

import expdata_module as exp

# Paths to import data files from the fitting_data folder
# INPUTS
ISOTOPES_TO_FIT_PATH = 'data/fitting_data/isotopes_to_fit.txt'
ISOTOPE_DATA_PATH = 'data/fitting_data/isotope_data.txt'
CONSTANTS_PATH = 'data/fitting_data/constants.txt'
CROSS_SECTION_PATH = 'data/cross_section_data/'

# OUTPUTS
ISOTOPE_OUTPUT_PATH = 'data/fitting_data/isotope_ref.txt'

def load_constant_data():
    """
    Loads the constants specified in constants.txt into a dictionary
    """
    dict = {}
    with open(CONSTANTS_PATH, mode='r') as file:
        header = file.readline().split(',')
        vals = file.readline().split(',')
        for i in range(len(header)):
            dict[header[i].strip()] = vals[i].strip()
    return dict

# CONSTANTS
# proton_pulse_gap, flight_path, trigger_delay, time_bin, min_E, max_E, crs_skip_pts
# contained in the CONST_DICT
CONST_DICT = load_constant_data()

# -------------------------------------------------------------------------------
def create_isotope_dataframe():
    """ 
    Returns pandas dataframe populated with isotope data,
    for isotopes to be fitted (from the parameter file). 
    """

    # Store the names of isotopes to be fitted from isotopes.txt in a list
    with open(ISOTOPES_TO_FIT_PATH, mode='r') as file:
        isotope_fitting_list = [line.strip() for line in file.readlines()]

    print("Isotopes to fit list:", isotope_fitting_list)

    # Crop isotope_data.txt file to include only the isotopes in the fitting ilst
    # and save as isotope_reference.txt
    with open(ISOTOPE_DATA_PATH, mode='r') as infile, open(ISOTOPE_OUTPUT_PATH, mode='w') as outfile:
        header = infile.readline()  # Read header
        outfile.write(header)  # Write the header to the output file
        for line in infile:
            file_name = line.split(',')[0]
            if file_name in isotope_fitting_list:
                outfile.write(line)

    # Create a dataframe from the new csv file (isotope_ref.txt) containing only the isotopes to be fitted
    df = pd.read_csv(ISOTOPE_OUTPUT_PATH, delimiter=',')
    return df, isotope_fitting_list

# Taken from transmission module, slightly modified
# Largely Anton's code, but with renamed variables for clarity    
def load_cross_section(isotope_name):
    CSdata = pd.read_table(CROSS_SECTION_PATH + isotope_name)
    # find the index of first row, where E exceeds Emin
    indMin = next(x for x, val in enumerate(CSdata[CSdata.columns[0]]) if val > float(CONST_DICT['min_E'])*1e-6)
    # find the index of first row, where E exceeds Emax
    indMax = next(x for x, val in enumerate(CSdata[CSdata.columns[0]]) if val > float(CONST_DICT['max_E'])*1e-6)
    # select only needed part of table
    CSdata = CSdata.iloc[indMin-1:indMax+1]

    # drop duplicate rows if such exist
    CSdata = CSdata.drop_duplicates(subset=CSdata.columns[0])

    # create numpy arrays
    Earr = CSdata[CSdata.columns[0]].to_numpy()*1e6
    CSarr = CSdata[CSdata.columns[1]].to_numpy()

    #interpolate the cross section to allow for any E value
    CSinterpolated = intrp.interp1d(Earr, CSarr, kind='linear')
    return CSinterpolated

def create_crs_dict(isotope_df):
    # Load interpolated cross section data, between min_E and max_E (specified in the constants)
    # Store in a dictionary with key as the cross-section name, and value as the interpolated function
    interp_crs = {}
    for isotope in isotope_df['file_name']:
        print(isotope)
        interp_crs[isotope] = load_cross_section(isotope)
    return interp_crs


# if __name__ == "__main__":
#     df = create_isotope_dataframe()
#     print("Isotope ref data:")
#     print(df.to_string())
#     print(CONST_DICT)
#     returnthing = load_cross_section(df.iloc[1]['file_name'])
    
#     # Generate a range of energy values to evaluate the interpolated function
#     E_values = np.linspace(0.001, 500, 500)

#     # Compute the cross-section values at these energy points
#     cross_section_values = returnthing(E_values)

#     # Plot the energy values against the cross-section values
#     plt.figure(figsize=(10, 6))
#     plt.plot(E_values, cross_section_values, label='Interpolated Cross Section')
#     plt.xlabel('Energy (eV)')
#     plt.ylabel('Cross Section')
#     plt.title('Interpolated Cross Section vs. Energy')
#     plt.xscale('log')
#     plt.yscale('log')
#     plt.grid(True)
#     plt.legend()
#     plt.show()

def E_from_TOF(TOF, L) :
    #Lamda = (TOF * 3.956 / L / 1000)
    Lamda = TOF * 3.956e-3 / L 
    #E = (6.626**2) / 2 / 1.675 / (Lamda**2) /1.6 / 100
    E = 0.0819102164 / (Lamda**2) 
    return E

def TOF_from_E(E, L) :
    TOF = 72.3457051 * L / np.sqrt(E)
    return TOF

#TODO: make a create_E_mesh() func
# what is this range and steps?
# creates array of E values with proper steps
def Get_E_Array_Properly_Spaced(Emin, Emax) :

    Erange = np.array( [0.0001, 0.001, 0.01, 0.1,  1,    10,   100,  1e3, 1e4, 1e5, 1e6] )
    Esteps = np.array( [5e-5,   5e-4,  5e-3, 5e-3, 1e-2, 1e-2, 1e-2, 0.1, 1,   10,  100] )

    if Emin < Erange[0]:
        raise Exception("Emin is too small!!!!\n Input value was {0} lower than Emin = {1}\n Exit here".format(Emin, Erange[0]))

    if Emax > 20e6:  # 20 MeV is max value in most cross section data
        raise Exception("Emax is too large!!!!\n Input value was {0} larger than Emax = {1}\n Exit here".format(Emax, Erange[-1]))

    # TODO: wrap in a function for clarity. Calculates the index of the maximum energy point in the Erange array. 
    MaxInd = len(Erange) - 1
    ind = MaxInd
    for i in range(0, MaxInd) :
        if Emin >= Erange[i] and Emin < Erange[i+1] :
            ind = i

    # if the max index is at the end if the erange array
    if ind == MaxInd:
        steps = int((Emax-Emin)/Esteps[ind]) #Esteps[ind] = Esteps[maxind] = 100 in this case always.
        Earr = np.linspace(Emin, Emax, num=steps)
        #Earr = np.arange(start=Emin, stop=Emax, step=Esteps[ind])
    else :
        # otherwise use the provided value of Emax to put an upper bound on the energy (unless the Erange[ind+1] is smaller)
        right = min(Emax,Erange[ind+1])
        left = Emin
        steps = int((right-left)/Esteps[ind])
        Earr = np.linspace(left, right, num=steps)
        #Earr = np.arange(start=Emin, stop=min(Emax,Erange[ind+1]), step=Esteps[ind])
        for i in range(ind+1, MaxInd) :
            right = min(Emax, Erange[i + 1])
            left = Erange[i]
            if right > left :
                steps = int((right - left) / Esteps[i])
                E1 = np.linspace(left, right, num=steps)
                #E1 = np.arange(start=Erange[i], stop=min(Emax,Erange[i+1]), step=Esteps[i])
                Earr = np.append(Earr, E1)
        if Emax > Erange[MaxInd] :
            right = Emax
            left = Erange[MaxInd]
            steps = int((right - left) / Esteps[MaxInd])
            E1 = np.linspace(left, right, num=steps)
            #E1 = np.arange(start=Erange[MaxInd], stop=Emax, step=Esteps[MaxInd])
            Earr = np.append(Earr, E1)

    Earr = np.unique(Earr)  #remve duplicates, if they exist
    #plt.plot(Earr[1:], np.diff(Earr), color='g',marker='.')
    #plt.grid(True)
    #plt.xscale('log')
    #plt.yscale('log')
    ##plt.xlim([0.08, 2])
    #plt.show()
    return Earr

# TODO: make a create_T_mesh() func
# creates array of T values with proper steps for the rane in between Emin and Emax (eV)
def Get_T_Array_Properly_Spaced(Emin, Emax, L) :

    Tmin = TOF_from_E(Emax, L)
    Tmax = TOF_from_E(Emin, L)
    
    Erange = np.array( [0.0001, 0.001, 0.01, 0.1,  1,    10,   100,  1e3, 1e4, 1e5, 1e6] )
    # delme section----------------
    #Esteps = np.array( [5e-5,   5e-4,  5e-3, 5e-3, 1e-2, 1e-2, 1e-2, 0.1, 1,   10,  100] )
    #Trange = Get_TOF_FromE(Erange[::-1], L)
    #Tsteps = Trange - Get_TOF_FromE(Erange[::-1]+Esteps[::-1], L) 
    #print('Tsteps', Tsteps)
    #print('Trange', Trange)
    # end delme section----------------

    # set Tseps here. Above code calculates these from Esteps
    Tsteps = np.array( [5e-3, 0.01, 0.05, 0.002, 0.005, 0.05, 0.5, 10, 100, 500, 1000] )
    # convert these values relative to L=14.6 m, which was used for calibration
    Tsteps = Tsteps * L / 14.6

    MaxInd = len(Trange) - 1
    ind = MaxInd
    for i in range(0, MaxInd) :
        if Tmin >= Trange[i] and Tmin < Trange[i+1] :
            ind = i

    if ind == MaxInd:
        steps = int((Tmax-Tmin)/Tsteps[ind])
        Tarr = np.linspace(Tmin, Tmax, num=steps)
        #Tarr = np.arange(start=Tmin, stop=Tmax, step=Tsteps[ind])
    else :
        right = min(Tmax,Trange[ind+1])
        left = Tmin
        steps = int((right-left)/Tsteps[ind])
        Tarr = np.linspace(left, right, num=steps)
        #Tarr = np.arange(start=Tmin, stop=min(Tmax,Trange[ind+1]), step=Tsteps[ind])
        for i in range(ind+1, MaxInd) :
            right = min(Tmax, Trange[i + 1])
            left = Trange[i]
            if right > left :
                steps = int((right - left) / Tsteps[i])
                T1 = np.linspace(left, right, num=steps)
                #T1 = np.arange(start=Trange[i], stop=min(Tmax,Trange[i+1]), step=Tsteps[i])
                Tarr = np.append(Tarr, T1)
        if Tmax > Trange[MaxInd] :
            right = Tmax
            left = Trange[MaxInd]
            steps = int((right - left) / Tsteps[MaxInd])
            T1 = np.linspace(left, right, num=steps)
            #T1 = np.arange(start=Trange[MaxInd], stop=Tmax, step=Tsteps[MaxInd])
            Tarr = np.append(Tarr, T1)

    Tarr = np.unique(Tarr)  #remve duplicates, if they exist
    ##plt.plot(Tarr[1:], np.diff(Tarr), color='g', marker='.')
    #plt.plot(Tarr, Tarr, color='g', marker='.')
    #plt.grid(True)
    #plt.xscale('log')
    ##plt.yscale('log')
    ##plt.xlim([70, 70.1])

    print('Created time array with',len(Tarr),'Time cells')
    return Tarr

def ideal_transmission_solid(isotope_df, E_mesh):
    # ExpTerm = Rho * d * w * A * S2 / (w * m * A) / 1.6605389e4
    ExpTerm = np.zeros(len(E_mesh))  # create array with 0 values for Exponential Term

    # # from transmission_module #TODO: did not edit yet
    # numAtoms = 0 #number of atoms
    # effThickness = self.Parameters.thickness #effective thickness = thickness * fraction
    # density = self.isotopeList[0].Density
    
    # denominator = 0
    # isotopeSum = 0
    
    # for isotope in self.isotopeList:
    #     denominator += isotope.AtomicMass * isotope.Abundance
        
    #     IdealCrossSection(isotope).attenuation(E) * isotope.Abundance
    
    # numAtoms = density / denominator
                        
    # #not exponentiated bc we don't know the effective thickness. should be determined from fitting
    # return - numAtoms * self.Parameters.thickness * isotopeSum


#  ----------------------------------------------------------------------------------
# This function calculates Theoretical transmission on a propoerly spaced Earr mesh
# S - dictionary of cross sections (key=IsotopeName as in data file name). Cross sect. is interpolated function!
# Par - parameters read from the Parameter file (including path length, etc. and the elemental composition table)
# Earr - existing array of E values on which Tr will be calculated
# Returns the calculated transmission for the entire set of elements
# ----------------------------------------------------------------------------------
def CalcIdealTransmission(S, Par, Earr, Plot=False):
    # ExpTerm = Rho * d * w * A * S2 / (w * m * A) / 1.6605389e4
    # calculate exponential term for transmission, group by group
    Elem = Par['Elmnts']
    ExpTerm = np.zeros(len(Earr))  # create array with 0 values for Exponential Term

    for GrpN in range(0, max(Elem['GrupN']) + 1):  # iterate over group numbers here
        GrpSubset = Elem.loc[Elem['GrupN'] == GrpN]  # select subset of rows with GrpN
        if len(GrpSubset) > 0:
            A = np.zeros(len(Earr))  # empty array for that group
            Mass = 0
            for i in range(len(GrpSubset)):  # iterate within single group here
                row = GrpSubset.iloc[i]
                A += row['Abundance'] * row['AtomicFraction'] * S[row['Isotope Name']](Earr) * row['Thickness (um)'] * row['Density (g/cm3)'] / 1.6605389e4
                Mass += row['AtomicFraction'] * row['AtomicMass']
                print(row['Isotope Name'])
            print('Mass=',Mass)
            ExpTerm += A / Mass
            print('-----------------')

    Tr1 = np.exp(-ExpTerm)

    if Plot :
        plt.figure(figsize=(6,5))
        plt.subplot(211)
        plt.plot(Earr, Tr1, label='Ideal theory', color='b')
        plt.xlabel("Energy (eV)")
        plt.xscale('log')
        plt.ylabel("Theor. transm.")
        plt.grid(True)
        #-------------
        plt.subplot(212)
        plt.plot(Get_TOF_FromE(Earr, Par['Flight path']), Tr1, label='Ideal theory', color='g')
        plt.xlabel("TOF (us)")
        plt.xscale('log')
        plt.ylabel("Theor. transm.")
        plt.grid(True)
        # plt.yscale('log')
        # plt.xlim([500,1200])
        plt.subplots_adjust(top=1, bottom=0.0, left=0.0, right=1, hspace=0.3, wspace=0.1)

        plt.show()

    return Tr1

# -------------------------------------------------------------------------------
# Convolve idealtransmission with the neutron pulse time profile
# EarrIdeal - array of energy values for which Tr needs to be convolved
# TrIdeal   - function interpolated for the ideal theoretical transmission
# Returns the function interpolated for convlolved transmission on interval [Emin,Emax]
# -------------------------------------------------------------------------------
def ConvolveTrWithPulseProfile(Par, EarrIdeal, TrIdeal, BEAM_ARR_DIM, Plot=False, FileOutput=False):
    # create the time array corresponding to energy array of theoretical cross ection
    Tarr1 = Get_TOF_FromE(EarrIdeal, Par['Flight path'])

    # interpolate transmission over new Tarr
    # Tr = interp1d(Tarr1, TrIdeal, kind='linear', fill_value='extrapolate')
    Tr1 = interp1d(Tarr1, TrIdeal, kind='linear')

    # create working arrays to be used in convolution, only needed here
    BeamProfile_T = np.zeros(BEAM_ARR_DIM)  # just make an empty array to be filled up
    BeamProfile_Amp = np.zeros(BEAM_ARR_DIM)  # just make an empty array to be filled up
    TrTemp = np.zeros(BEAM_ARR_DIM)  # just make an empty array to be filled up

    # select only fraction of Earr between Emin and Emax, otherwise interpolation does not work for BeamProfile width added on the right side
    #Earr2 = EarrIdeal[np.where((EarrIdeal >= Par['Minimum E']) & (EarrIdeal <= Par['Maximum E']))]
    # leave only those elements which have enough data for convolution  - cut the left side of the array
    # subtract the time length equal to beam profile width at min Energy
    Tmax = Get_TOF_FromE(min(EarrIdeal), Par['Flight path']) - bm.BeamProfileWidth(min(EarrIdeal))*1.05  # 1.05 is the safety margin
    Emin = Get_E_FromTOF(Tmax, Par['Flight path'])
    Earr2 = EarrIdeal[np.where((EarrIdeal >= Emin))]
    Earr2 = Earr2[0::Par['CROSS_SECT_SKIP_POINTS']]  # take each N-th element in the Energy properly spaced array to speed up calculation
    #include last point Emax in Earr2
    if Earr2[-1] < EarrIdeal[-1]: Earr2 = np.append(Earr2, max(EarrIdeal))

    Tr2 = np.zeros(len(Earr2))  # just make an empty array to be filled up
    Tarr2 = Get_TOF_FromE(Earr2, Par['Flight path']) # pay attention here that Tarr2 is in reverse order, descending

    print('Will calculate', len(Earr2), 'points in time between', Earr2[0], 'eV and', Earr2[-1], 'eV')

    for i in tqdm_notebook(range(len(Earr2)), total=len(Earr2), desc="Convolving with beam profile", bar_format="{l_bar}{bar} [ time left: {remaining} ]") :
#        for i in range(len(Earr2)):
        E = Earr2[i]
        bm.BeamProfileArrayCalculated(Par['Proton pulse gap'], E, BeamProfile_T, BeamProfile_Amp, BEAM_ARR_DIM)
        TrTemp = Tr1(BeamProfile_T + Tarr2[i])  # these are arrays, not scalars, times as BmProfile shifted by Tarr4[i]
        Tr2[i] = np.sum(np.multiply(BeamProfile_Amp, TrTemp))
        len(Earr2)

    # interpolate convolved transmission
    # TrConvolved = interp1d(Tarr2, Tr2, kind='linear', fill_value='extrapolate')
    TrConvolved = interp1d(Tarr2, Tr2, kind='linear')

    #print('DONE with calucluations!\n-----------')

    if Plot :
        plt.figure(figsize=(6,5))
        plt.subplot(211)
        plt.plot(Get_E_FromTOF(Tarr2, Par['Flight path']), Tr2, label='Exper. predicted', color='b')
        plt.xlabel("Energy (eV)")
        plt.xscale('log')
        plt.ylabel("Theor. transm.")
        plt.grid(True)
        #-------------
        plt.subplot(212)
        plt.plot(Tarr2, Tr2, label='Exper. predicted', color='g')
        plt.xlabel("TOF (us)")
        plt.xscale('log')
        plt.ylabel("Theor. transm.")
        plt.grid(True)
        # plt.yscale('log')
        # plt.xlim([500,1200])
        plt.subplots_adjust(top=1, bottom=0.0, left=0.0, right=1, hspace=0.3, wspace=0.1)

        plt.show()

    if FileOutput :
        FileArr = np.array([EarrIdeal, TrIdeal])
        FileArr = FileArr.T
        np.savetxt('TrIdeal.csv', FileArr, fmt='%f', delimiter=',', header='E (eV),Tr Ideal')
        FileArr = np.array([Get_E_FromTOF(Tarr2, Par['Flight path']), Tr2])
        FileArr = FileArr.T
        np.savetxt('TrBeamline.csv', FileArr, fmt='%f', delimiter=',', header='E (eV),Tr beamline')

    return TrConvolved, Tarr2[0]


# TODO: pasted from old code, adjust/rewrite to fit into new code
# -------------------------------------------------------------------------------
# Fitting function for L and dT
# FitPar - list of two parameters L and dT in FitPar[0] and FitPar[1]
# Texp, TrExp_T - experimental data cross section
# EarrIdeal, TrIdeal_E - energy array and ideal theoretical cross section function
# Par - parameters of the model
# -------------------------------------------------------------------------------
def FuncToMinimize_PathCalibration(FitPar, Texp, TrExp_T, EarrIdeal, TrIdeal_E, Par):
    Par['Flight path'] = FitPar[0]
    Par['Trigger delay'] = FitPar[1]
    print('L=', FitPar[0], 'dT=', FitPar[1])
    # ----------------------------------------------------------------------------------
    # Convolve ideal Theoretical transmission with beam pulse profile
    # this depends on L and dT parameters, needs to be within minimizer
    TrConv_T, Tmax = ConvolveTrWithPulseProfile(Par, EarrIdeal, TrIdeal_E, BEAM_ARR_DIM)
    # Calculate the difference from experimental points
    ChiSq = CalcTransmForExpPoints(Par, Texp, TrExp_T, TrConv_T, Tmax, NsubCells)

    return ChiSq



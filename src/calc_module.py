



#Largely Anton's code, but with renamed variables for clarity    
def load_cross_section(self):
    CSdata = pd.read_table('CrossSections_BeamProfiles\\'+self.isotopeName+'.txt')
    
    # find the index of first row, where E exceeds Emin
    indMin = next(x for x, val in enumerate(CSdata[CSdata.columns[0]]) if val > self.Emin*1e-6)
    # find the index of first row, where E exceeds Emax
    indMax = next(x for x, val in enumerate(CSdata[CSdata.columns[0]]) if val > self.Emax*1e-6)
    # select only needed part of table
    CSdata = CSdata.iloc[indMin-1:indMax+1]

    # drop duplicate rows if such exist
    CSdata = CSdata.drop_duplicates(subset=CSdata.columns[0])

    # create numpy arrays
    Earr = CSdata[CSdata.columns[0]].to_numpy()*1e6
    CSarr = CSdata[CSdata.columns[1]].to_numpy()

    #interpolate the cross section to allow for any E value
    CSinterpolated =  interp1d(Earr, CSarr, kind='linear')

    return CSinterpolated


def calc_ideal_tr():















# 
# 
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



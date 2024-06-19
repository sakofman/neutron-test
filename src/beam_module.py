"""
Module Name: beam_module

Description
-----------
Provides Beam class for representing and manipulating beam parameters.

Date
----
06/18/2024

Corrections
-----------
    - Profile Array Calculation (?)
    - Implement function to calculate C by integrating non-normalized function from 0 to t_max (?)
    - (later) Attribute setting method for Beam constructor to interface with GUI 
"""
import numpy as np
import scipy.integrate as integrate

class Beam:
    """
    A class to represent a beam specified by the parameters of the Cole-Windsor function.
    """
    def __init__(self):
        """
        Constructor for a beam instance.
        """
        # TODO: Parser code for attribute setting (decide what would best interface with a GUI module, for now keeping as code constants)
        self.t_0 = lambda E: 2.27e-2 + 2.03 * pow(E, -0.46)
        self.gamma_1 = lambda E: 2.95e-2 + 0.905 * pow(E, 0.343)
        self.gamma_2 = lambda E: 6.78e-2 + 9.77e-2 * pow(E, 0.447)
        self.sigma_1 = lambda E: 6.78e-3 + 0.658 * pow(E, -0.468)
        self.sigma_2 = lambda E:3.15e-2 + 1.71 * pow(E, -0.476)
        self.R = lambda E: 0.404 - 0.29 * np.exp(-2.78e-4 * E)
        self.C = 1 # Without normalization

        # TODO: Check if integration calc for C is ok -- for now keeping as a provided constant
        # self.C = self.normalization_factor
        self.C = lambda E: 8.03e18 * np.exp((-8.21*pow(E, 0.0542)))

    # def normalization_factor(self, E): 
    #     """
    #     Returns C (the normalization factor) computed by integrating the beam function from 0 to t_max, keeping energy as a variable.
    #     """
    #     return integrate.quad(self.pulse_shape(E), 0, self.profile_width(E))
    
    def profile_width(self, E) :
        # t_max is from the calibration I (Anton) did for NOBORU function, see BeamProfileWidth_Vs_Energy_June2020.xlsx
        E_1 = np.log(E)
        t_max = ((((5.113247E-06 * E_1 - 8.185929E-05) * E_1 - 3.003155E-04) * E_1 + 1.816361E-02) * E_1 - 3.318026E-01) * E_1 + 2.308175E+00
        t_max = np.exp(t_max)
        return t_max
    
    def pulse_shape(self, t, E):
        """
        Returns the Cole-Windsor function value at time t and energy E with its parameters 
        set according to the attributes of the instance.

        Arguments
        ---------
        t : float
            Time at which to evaluate the pulse shape.
        E : float
            Energy at which to evaluate the pulse shape.

        Returns
        -------
        float
            The value of the Cole-Windsor function at the given time t and energy E.
        """
        delta_t = t - self.t_0(E)
        if (t < self.t_0(E)):
            F_1 = np.exp(-0.5 * ((delta_t / self.sigma_1(E)) ** 2))
            F_2 = F_1
        else:
            to_compare_1 = self.t_0(E) + self.gamma_1(E) * (self.sigma_2(E) ** 2)
            to_compare_2 = self.t_0(E) + self.gamma_2(E) * (self.sigma_2(E) ** 2)
            if (t < to_compare_1):
                F_1 = np.exp(-0.5 * ((delta_t / self.sigma_2(E)) ** 2))
            else:
                F_1 = np.exp((((self.gamma_1(E) * self.sigma_2(E))**2) / 2) - (self.gamma_1(E) * delta_t))
            if (t < to_compare_2):
                F_2 = np.exp(-0.5 * ((delta_t / self.sigma_2(E)) ** 2))
            else:
                F_2 = np.exp((((self.gamma_2(E) * self.sigma_2(E))**2) / 2) - (self.gamma_2(E) * delta_t)) # TODO: check if order of pow/mult matter here for our range of #s
                                                                                                  # or just raise to powers first as in the prev code
        F = self.C(E) * ( (1 - self.R) * F_1 + self.R * F_2)
        return F


# TODO: ProfileArray Calculation

# #---------------------------------------------------------------------------------
# # Function approximating pulse shape (See Hasemi publication for NOBORU pulse
# # Time is the time of the pulse relative to t0 (us) and E is energy of the neutron in eV
# #---------------------------------------------------------------------------------
# def BeamPulseShape(Time, E): 

#     C = 8.03e18 * np.exp((-8.21*pow(E, 0.0542)))
#     to = 2.27e-2 + 2.03 * pow(E, -0.46)
#     gamma1 = 2.95e-2 + 0.905 * pow(E, 0.343)
#     gamma2 = 6.78e-2 + 9.77e-2 * pow(E, 0.447)
#     sigma1 = 6.78e-3 + 0.658 * pow(E, -0.468)
#     sigma2 = 3.15e-2 + 1.71 * pow(E, -0.476);
#     R = 0.404 - 0.29 * np.exp(-2.78e-4 * E)

#     if(Time < to):
#         F1 = np.exp((-0.5*((Time-to)**2)/(sigma1**2)))
#         F2 = F1
#     else :
#         if Time < to+gamma1*(sigma2**2): 
#             F1 = np.exp((-0.5*((Time-to)**2)/(sigma2**2)))
#         else:
#             F1 = np.exp((0.5*(gamma1**2)*(sigma2**2) - gamma1*(Time-to)))
            
#         if Time < to+gamma2*(sigma2**2):
#             F2 = np.exp((-0.5*((Time-to)**2)/(sigma2**2)))
#         else:
#             F2 = np.exp((0.5*(gamma2**2)*(sigma2**2) - gamma2*(Time-to)))

#     F = C * ( (1-R)*F1 + R*F2)

#     return F

# # ---------------------------------------------------------------------------------------
# # Calculates Beam profile width in us
# # uses pre-calibrated for NOBORU values, which I approximated by 5th degree polynomial
# # see file BeamProfileWidth_Vs_Energy_June2020.xlsx
# # returns the beam width in us
# # Input parameter is E in eV
# # ---------------------------------------------------------------------------------
# def BeamProfileWidth(E) :
#     # Tmax is from the calibration I did for NOBORU function, see BeamProfileWidth_Vs_Energy_June2020.xlsx
#     E1 = np.log(E)
#     Tmax = ((((5.113247E-06 * E1 - 8.185929E-05) * E1 - 3.003155E-04) * E1 + 1.816361E-02) * E1 - 3.318026E-01) * E1 + 2.308175E+00
#     Tmax = np.exp(Tmax)
#     return Tmax

# #---------------------------------------------------------------------------------
# # Calculates the pulse shape from the moderator. Can be two pulses with gap between them
# # ProtonPulseGap - us, E in eV, 
# # T - (us) array to be filled in of time values for which beam profile to be calculated, relative to T0
# # BmProf - array to be filled in, iBmProfDim - size of the arrays
# # Time is the time of the pulse relative to t0 (us) and E is energy of the neutron in eV
# #---------------------------------------------------------------------------------
# def BeamProfileArrayCalculated(ProtonPulseGap, E, T, BmProf, iBmProfDim) :
    
#     # Tmax is the width of beam pulse
#     Tmax = BeamProfileWidth (E)

#     Tstep = Tmax / iBmProfDim   # time step of the array to be calculated
    
#     BmProf_Pulse1 = np.zeros(iBmProfDim)

#     # fill in T array of times
#     for i in range(iBmProfDim):
#         T[i] = Tstep * i
    
#     #calculate the single pulse shape here
#     for i in range(iBmProfDim):
#         BmProf_Pulse1[i]  = BeamPulseShape( T[i], E)

#     #interpolate the Shape of one pulse for the second pulse
#     BmProf_Pulse2 =  interp1d(T, BmProf_Pulse1, kind='cubic')
    
#     #add the second pulse to the first one after ProtonPulseGap delay
#     for i in range(iBmProfDim):
#         if T[i] > ProtonPulseGap:
#             BmProf[i] = BmProf_Pulse1[i] + BmProf_Pulse2(T[i]-ProtonPulseGap)
#         else:
#             BmProf[i] = BmProf_Pulse1[i]
            
#     #normalize the BeamPulse shape to unity
#     sum = np.sum(BmProf)
#     BmProf /= sum

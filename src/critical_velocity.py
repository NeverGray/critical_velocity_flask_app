# Project Name: Beyer-Stacey=Brenn CV Calculator
# Description: Calculate Critical Velocity
# Copyright (c) 2025 Justin Edenbaum, Never Gray
#
# This file is licensed under the MIT License.
# You may obtain a copy of the license at https://opensource.org/licenses/MIT

from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Constants and empirical parameters
g = 9.81 #Gravitaional acceleration (m/s^2)
f_h = 0.176 #Empirical factor relating L_n to the length of that part of the fire that influences backlayering propensity (-)
f_a = 0.95 #Empirical constant used in fire length function (-)
f_o = 0.88 #Empirical constant used in fire length function (-
gamma_L = 2.2 #Empirical constant of fire length function (-)
f_g = 0.58 #Empirical constant to bias the 'plume blockage' function (-)
a_w = 1.82 #Empirical constant to offset the 'plume blockage' function (-)
m = 0.35 #Spreading rate of rising plume relevant for interaction of ambient air with plume (-)


@dataclass
class Fire:
    """A class to store fire parameters"""
    hrr: float #Toal fire heat release rate (W)
    intensity: float #Fire intensity (W/m^2)
    width: float #Maximum width of the fire (m)
    epsilon: float #Heat reduction due to imperfect combustion and thermal radiation (-)
    eta: float #Reduction of the convective heat release rate due to water mist systems or sprinklers (-)

    def length(self):
        """Correlation to calculate the total length of the fire source (m), Equation (34)"""
        return self.hrr / (self.intensity * self.width)


@dataclass
class Tunnel:
    """A class to store tunnel parameters"""
    name: str
    height: float #Tunnel height for simplified equation and design purpose or tunnel height from the base of the fire to highest point at the ceiling (m)
    area: float #Tunnel cross-sectional area at the fire site (m^2)
    hydraulic_diameter: float #Tunnel hydraulic diameter at the fire location (m)


def cp_air_low(T):
    """Polynomial for Cp of air from 273.15 K to 1000 K"""
    return 1.0484001086e3 \
           - 3.5085596831e-1 * T \
           + 8.1471072329e-4 * T ** 2.0 \
           - 3.7124712609e-7 * T ** 3.0
           
def cp_air_high(T):
    """Polynomial for Cp of air from 1000 K to 3000 K"""
    return 8.4977500000e2 \
           + 4.2107968837e-1 * T \
           - 1.4924271734e-4 * T ** 2.0 \
           + 1.9395639273e-8 * T ** 3.0

def air_specific_heat(T1, T2):
    """Function to determine the specific heat capacity of air between two 
    temperatures used to derive the empirical parameters"""
    if (T1+T2)/2 <= 1000:
        return cp_air_low((T1+T2)/2)
    elif ((T1+T2)/2 > 1000 and (T1+T2)/2 <= 3000):
        return cp_air_high((T1+T2)/2)
    else:
        print("WARNING!  Temperature is too high.  Temperature is clipped to 3000 K")
        return cp_air_high(3000.0)
        

"""This thermodynamicly more appropriate function for calculating the specific heat capacity is not used.  
    Using this function would require adjustment of the empirical parameters as these parameters were derived 
    based on fitting the equation to real data using the more simplistic correlation above"""
#def air_specific_heat(T1, T2):
#    """Function to determine the specific heat capacity of air between two 
#    temperatures"""
#    if (T1 <= 1000 and T2 <= 1000):
#        return (cp_air_low(T2)*(T2-273.15) - cp_air_low(T1)*(T1-273.15))/(T2-T1)
#    elif (T2 > 1000 and T2 <= 3000 and T1<=1000):
#        return (cp_air_high(T2)*(T2-273.15) - cp_air_low(T1)*(T1-273.15))/(T2-T1)
#    else:
#        print("WARNING!  Temperature is too high.  Temperature is clipped to 3000 K")
#       return (cp_air_high(3000.0)*(3000.0-273.15) - cp_air_low(T1)*(T1-273.15))/(3000-T1)
        

def calculate_delta_t(fire: Fire, tunnel: Tunnel, critical_velocity, ambient_density, specific_heat):
    """Function to calculate delta T_p (relevant for buoyancy force) based on initial/previous critical velocity value.
    Based on Equations (32) and (36)
    """
    effective_fire_length = (f_a + f_o) \
                          / (f_a + np.exp(-(fire.length()) / (gamma_L * tunnel.height)))

    return min(fire.hrr,
               fire.intensity * fire.width * tunnel.height * f_h * effective_fire_length) \
           * (1 - fire.epsilon)*(1-fire.eta) \
           / (min(tunnel.area, tunnel.height * (fire.width + m * tunnel.height)) \
           * critical_velocity * ambient_density * specific_heat)


def calculate_delta_t_simple(fire, tunnel, critical_velocity):
    """Simplified function to calculate delta T_p (relevant for buoyancy force) based on initial/previous critical velocity value.
    Based on Equation (38)
    """
    return (2.92e-4 * critical_velocity * (7.14 + tunnel.height)*(0.95 + np.exp(-8.0e-8 * fire.hrr / tunnel.height))) ** -1


def calculate_critical_velocity(fire, tunnel, delta_T, ambient_T):
    """Function to calculate critical velocity. 
    Based on equation (30) and (35)
    """
    K_F = f_g * a_w * (1.0 - min(1.0, (fire.width * tunnel.height)/tunnel.area)) + 1.0 - f_g

    return K_F * np.sqrt((g * tunnel.height ** 3.0 / tunnel.hydraulic_diameter ** 2.0) * delta_T / (ambient_T + delta_T) )


def calculate_critical_velocity_simple(fire, tunnel, delta_t):
    """Simplified function to calculate critical velocity.
    Based on equation (37)
    """
    return (1.476 - 2.64 * tunnel.height / tunnel.area) \
         * ((g * tunnel.height ** 3.0 / tunnel.hydraulic_diameter ** 2.0) * delta_t / (294.0 + delta_t)) ** 0.5


def iterate_critical_velocity(fire, tunnel, critical_velocity, ambient_T, ambient_density, epsilon, eta, K_g, tol):
    """Iterative critical velocity calculation.
    """
    fire.epsilon = epsilon
    fire.eta = eta
    relaxation = 0.7 #Relaxation factor to reduce number of iterations
    specific_heat = 1007.0
    while True:
        fire_dt = fire.hrr * (1.0 - fire.epsilon)*(1.0 - fire.eta) / (specific_heat * ambient_density * tunnel.area * critical_velocity)
        specific_heat = air_specific_heat(ambient_T, ambient_T + fire_dt)
        old_critical_velocity = critical_velocity
        dt = calculate_delta_t(fire, tunnel, critical_velocity, ambient_density, specific_heat)
        critical_velocity = relaxation*calculate_critical_velocity(fire, tunnel, dt, ambient_T)+ (1.0-relaxation)*old_critical_velocity
        if abs(critical_velocity - old_critical_velocity) < tol:
            return critical_velocity*K_g, dt #Grade factor applied after the iteration according to Equation (28)


def iterate_critical_velocity_simple(fire, tunnel, critical_velocity, tol):
    """Iterative critical velocity solver based on simplified equations.
    Based on suggested parameters (see at the end of the script) and constant specific heat capacity of 1007 J/(kg K)
    """
    relaxation = 0.7 #Relaxation factor to reduce number of iterations
    while True:
        old_critical_velocity = critical_velocity
        dt = calculate_delta_t_simple(fire, tunnel, critical_velocity)
        critical_velocity = relaxation*calculate_critical_velocity_simple(fire, tunnel, dt)+ (1.0-relaxation)*old_critical_velocity
        if abs(critical_velocity - old_critical_velocity) < tol:
            return critical_velocity, dt


def oxygen_depletion(fire_hrr, tunnel, ambient_density):
    """Correlation for minimum required oxygen demand (minimum upstream air velocity) based on tunnel parameters and HRR
    """
    air_demand = 14.16 # Based on C19H30 (kg_air/kg_fuel)
    heat_of_combustion = 42.68 #Based on Memorial Tunnel test (MJ/kg)
    combustion_eff = 0.9 #Combustion efficiency due to imperfect combustion (-)
    vel_depletion = air_demand*fire_hrr*1.0e-6/(ambient_density*tunnel.area*heat_of_combustion*combustion_eff)
    return vel_depletion


def plot_critical_velocity(fig, ax, tunnel, min_hrr, max_hrr, ambient_temp, ambient_density, fire_intensity, fire_width, epsilon, eta, K_g, tol, for_web=True):
    """Function to create a plot of critical velocity against HRR and return the figure and axis
    """
    fire_hrrs = np.linspace(min_hrr, max_hrr, 1491)
    critical_velocity = np.empty_like(fire_hrrs)
    critical_velocity[-1] = 1.0
    delta_t = np.empty_like(fire_hrrs)
    marker=0
    for i, fire_hrr in enumerate(fire_hrrs):
        fire = Fire(fire_hrr, fire_intensity, fire_width, epsilon, eta)
        critical_velocity[i], delta_t[i] = iterate_critical_velocity(fire, tunnel, critical_velocity[i-1], ambient_temp, ambient_density, epsilon, eta, K_g, tol)
        
        #HRR cut off if minimum oxygen requirement not met
        if (critical_velocity[i] <= oxygen_depletion(fire_hrr, tunnel, ambient_density) and marker == 0):
            vel_depletion = oxygen_depletion(fire_hrr, tunnel, ambient_density)
            critical_velocity[i] = vel_depletion
            hrr_depletion =  fire_hrr
            fire_hrrs[i] = hrr_depletion
            marker =1
            print("!!!!HRR cut off due to oxygen depletion!!!!")
        elif marker != 0:
            critical_velocity[i] = vel_depletion
            fire_hrrs[i] = hrr_depletion
            
    ax.plot(fire_hrrs*1e-6, critical_velocity*K_g, label=tunnel.name)
    
# Save CSV file with HRR (MW) and critical velocity (m/s) for each listed tunnel
    if not for_web:
        DF = pd.DataFrame(critical_velocity, fire_hrrs*1.0e-6)
        DF.to_csv(f"{tunnel.name} .csv")
    return fig, ax



def plot_critical_velocity_simple(fig, ax, tunnel, min_hrr, max_hrr, tol, to_csv=False):
    """Function to create a plot of critical velocity  against HRR based on simplified equations and return the figure and axis
    """
    fire_hrrs = np.linspace(min_hrr, max_hrr, 1491)
    critical_velocity = np.empty_like(fire_hrrs)
    critical_velocity[-1] = 1
    delta_t = np.empty_like(fire_hrrs)
    marker=0
    for i, fire_hrr in enumerate(fire_hrrs):
        fire = Fire(fire_hrr, 2.25e6, 2.5, 0.2, 0.0)
        critical_velocity[i], delta_t[i] = iterate_critical_velocity_simple(fire, tunnel, critical_velocity[i-1], tol)
        
        #HRR cut off if minimum oxygen requirement not met
        if (critical_velocity[i] <= oxygen_depletion(fire_hrr, tunnel, 1.2) and marker == 0):
            vel_depletion = oxygen_depletion(fire_hrr, tunnel, 1.2)
            critical_velocity[i] = vel_depletion
            hrr_depletion =  fire_hrr
            fire_hrrs[i] = hrr_depletion
            marker =1
            print("!!!!HRR cut off due to oxygen depletion (simplified)!!!!")
        elif marker != 0:
            critical_velocity[i] = vel_depletion
            fire_hrrs[i] = hrr_depletion
            
    ax.plot(fire_hrrs*1.0e-6, critical_velocity, label=f"{tunnel.name} simplified", linestyle="--")
     
    # Save CSV file with HRR (MW) and critical velocity (m/s) based on smplified equations for each listed tunnel 
    DF = pd.DataFrame(critical_velocity, fire_hrrs*1.0e-6)
    DF.to_csv(f"{tunnel.name} simplified.csv")
    return fig, ax

def fire_response(tunnels, min_hrr, max_hrr, tol, ambient_temp, ref_pressure, fire_intensity, fire_width, epsilon, eta, K_g, for_web=False):
    """Function to summarise critical velocity calculation of a list of tunnels
    """

    ambient_density = ref_pressure/(287.0 * ambient_temp)  #Ambient density calculation based on reference pressure (kg/m^3)
    

    fig, ax = plt.subplots()

    for tunnel in tunnels:
        fig, ax = plot_critical_velocity(fig, ax, tunnel, min_hrr, max_hrr, ambient_temp, ambient_density, fire_intensity, fire_width, epsilon, eta, K_g, tol)
        if not for_web:
            if (tunnel.height < 9.0 and tunnel.area > tunnel.height*(2.5+0.35*tunnel.height)):
                fig, ax = plot_critical_velocity_simple(fig, ax, tunnel, min_hrr if min_hrr>10.0e6 else 10.0e6, max_hrr, tol)

    ax.set(title="Critical Velocity",
           xlabel="Total Fire Heat Release Rate (MW)",
           xticks=np.arange(0.0, max_hrr*1.0e-6+1, 10.0),
           xlim=[0.0, max_hrr*1.0e-6],
           ylabel="Critical Velocity (m/s)",
           yticks=np.arange(0.0, 4.1, 0.5),
           ylim=[0.0, 4.0])
    ax.grid(True)

    ax.legend()
    fig.tight_layout()
    fig.set_size_inches(8.0, 4.5)
    if not for_web:
        fig.savefig("Fire_response")
        plt.show()
    x=2.0*10**-7
    print(1+x-1)
    if for_web:
        return fig

def plot_values():
# Input: Tunnel name, tunnel height (Ln or H) (m), tunnel area (m^2), hydraulic diameter (m)
    Tunnel_A = Tunnel("Tunnel-A", 6.0, 30.0, 5.45)
    Tunnel_B = Tunnel("Tunnel-B", 6.0, 30.0, 5.45)
    #Table 8 Typical TBM or Arched tunnel profiles
    TBM1 = Tunnel("TBM-1", 6.80, 85.48, 9.06)
    TBM2 = Tunnel("TBM-2", 8.07, 85.60, 9.28)
    TBM3 = Tunnel("TBM-3", 9.82, 106.0, 10.7)
    TBM4 = Tunnel("TBM-4", 5.00, 40.00, 5.99)
    TBM5 = Tunnel("TBM-5", 5.89, 40.00, 6.44)
    TBM6 = Tunnel("TBM-6", 7.86, 59.95, 8.13)
    TBMs = [TBM1, TBM2, TBM3, TBM4, TBM5, TBM6]
# Input: List of tunnel names, minimum HRR (W), maximum HRR (W), tolerance of iterative calculations, ambient temperature (K), 
#        reference pressure (Pa), fire intensity (MW/m^2), fire width (m), epsilon (-), eta (-), grade factor (-)
#Suggested parameters: [Tunnel_A], 1e6, 150e6, 1e-6, 294, 101325.0, 2.25e6, 2.5, 0.2, 0.0, 1.0
    fire_response(TBMs, 1e6, 200e6, 1e-5, 294, 101325.0, 2.25e6, 2.5, 0.2, 0.0, 1.0)
    
if __name__ == "__main__":
    #Fire Parameters
    hrr = 50.0e6 #Toal fire heat release rate (W)
    intensity= 2.25e6 #Fire intensity (W/m^2)
    width = 2.5 #Maximum width of the fire (m)
    epsilon = 0.2 #Heat reduction due to imperfect combustion and thermal radiation (-)
    eta = 0.0 #Reduction of the convective heat release rate due to water mist systems or sprinklers (-)
    fire = Fire(hrr, intensity, width, epsilon, eta)
    #Tunnel Paraemters
    name = "TBM-1"
    height = 6.8 #Tunnel height for simplified equation and design purpose or tunnel height from the base of the fire to highest point at the ceiling (m)
    area = 85.48 #Tunnel cross-sectional area at the fire site (m^2)
    hydraulic_diameter = 9.06 #Tunnel hydraulic diameter at the fire location (m)
    tunnel = Tunnel(name, height, area, hydraulic_diameter)
    critical_velocity = 1.0 #Initial guess
    ambient_T = 294.0 #Kelvin
    ambient_density = 1.2 #kg/m^3
    epsilon = 0.2 
    eta = 0.0
    K_g = 1.0
    tol = 6.8e-6
    critical_velocity, dt = iterate_critical_velocity(fire, tunnel, critical_velocity, ambient_T, ambient_density, epsilon, eta, K_g, tol)
    print(f"Critical Velocity: {critical_velocity:.3f} m/s, Delta T: {dt:.1f} K")
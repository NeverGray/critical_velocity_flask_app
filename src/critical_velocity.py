# Project Name: Beyer-Stacey-Brenn Critical Velocity (BSB-CV) Calculator
# Description: Calculate and plot values for Critical Velocity
# Copyright (c) 2025 Justin Edenbaum, Never Gray
#
# This file is licensed under the MIT License.
# License is available at https://opensource.org/licenses/MIT

from dataclasses import dataclass
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.lines as mlines
import copy

VERSION_NUMBER = "1.1.0"
# Constants and empirical parameters
g = 9.81 #Gravitaional acceleration (m/s^2)
f_h = 0.176 #Empirical factor relating L_n to the length of that part of the fire that influences backlayering propensity (-)
f_a = 0.95 #Empirical constant used in fire length function (-)
f_o = 0.88 #Empirical constant used in fire length function (-)
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

def round_up_nice(val):
    """ Function to round up to a "nice" value: 1, 1.5, 2, 2.5, 5, 7.5, 10 × 10^k"""
    if val <= 0:
        return 1.0
    exp = np.floor(np.log10(val))
    frac = val / (10**exp)
    for step in [1, 1.5, 2, 2.5, 5, 7.5, 10]:
        if frac <= step:
            return step * (10**exp)
    return 10 * (10**exp)

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

def calculate_critical_velocity(fire, tunnel, delta_T, ambient_T):
    """Function to calculate critical velocity. 
    Based on equation (30) and (35)
    """
    K_F = f_g * a_w * (1.0 - min(1.0, (fire.width * tunnel.height)/tunnel.area)) + 1.0 - f_g

    return K_F * np.sqrt((g * tunnel.height ** 3.0 / tunnel.hydraulic_diameter ** 2.0) * delta_T / (ambient_T + delta_T) )

def iterate_critical_velocity(fire, tunnel, critical_velocity, ambient_T, ambient_density, K_g=1.0, tol=6.8e-6):
    """Iterative critical velocity calculation.
    """
    relaxation = 0.7 #Relaxation factor to reduce number of iterations
    specific_heat = 1007.0 #Starting value for specific heat capacity of air (J/kg/K)
    iter = 0
    while True:
        fire_dt = fire.hrr * (1.0 - fire.epsilon)*(1.0 - fire.eta) / (specific_heat * ambient_density * tunnel.area * critical_velocity)
        
        specific_heat = air_specific_heat(ambient_T, ambient_T + fire_dt)
        old_critical_velocity = critical_velocity
        dt = calculate_delta_t(fire, tunnel, critical_velocity, ambient_density, specific_heat)
        critical_velocity = relaxation*calculate_critical_velocity(fire, tunnel, dt, ambient_T)+ (1.0-relaxation)*old_critical_velocity
        iter += 1
        if abs(critical_velocity - old_critical_velocity) < tol or iter > 100:
            return critical_velocity*K_g, dt, iter #Grade factor applied after the iteration according to Equation (28)

def oxygen_depletion(fire_hrr, tunnel, ambient_density):
    """Correlation for minimum required oxygen demand (minimum upstream air velocity) based on tunnel parameters and HRR
    """
    air_demand = 14.16 # Based on C19H30 (kg_air/kg_fuel)
    heat_of_combustion = 42.68 #Based on Memorial Tunnel test (MJ/kg)
    combustion_eff = 0.9 #Combustion efficiency due to imperfect combustion (-)
    vel_depletion = air_demand*fire_hrr*1.0e-6/(ambient_density*tunnel.area*heat_of_combustion*combustion_eff)
    return vel_depletion

def plot_critical_velocity(tunnel, fire, ambient_temp, ambient_pressure, for_web=True):
    """Generate a plot of fire heat release rate vs critical velocity - optimized version"""
    # Calculate the density for the ambient pressure
    ambient_density = ambient_pressure/(287.0 * ambient_temp)  #Ambient density calculation based on ambient pressure (kg/m^3)
    
    # Get data and calculate maximums in one step
    fire_hrrs, critical_velocities, sufficient_oxygen, converging_msg = hrrs_vs_critical_velocities(tunnel, fire, ambient_temp, ambient_density, for_web=for_web)
    max_critical_velocity = max(critical_velocities)
    max_fire_hrr = max(fire_hrrs)

    # Determine oxygen status and critical velocity
    iter = 0
    if sufficient_oxygen:
        oxygen_depletion_msg = "" # No message if sufficient oxygen
        critical_velocity, dt, iter = iterate_critical_velocity(fire, tunnel, max_critical_velocity, ambient_temp, ambient_density)
    else:
        if fire.hrr < max_fire_hrr:
            oxygen_depletion_msg = f"There is insufficient oxygen for fires greater than {max_fire_hrr *1e-6:.1f} MW"
            critical_velocity, dt, iter = iterate_critical_velocity(fire, tunnel, max_critical_velocity, ambient_temp, ambient_density)   
        else:
            oxygen_depletion_msg = f"HRR cut off due to oxygen depletion at {max_fire_hrr*1e-6:.1f} MW"
            fire.hrr = max_fire_hrr 
            critical_velocity = max_critical_velocity

    if iter > 100:
        converging_msg = "WARNING! Critical velocity calculation did not converge. Please check input parameters."
    
    # Create plot with streamlined formatting
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    ax.plot(fire_hrrs*1e-6, critical_velocities)
    
    # Create parameter legend (same as before but more compact)
    param_lines = [
        mlines.Line2D([], [], color='none', label=f"{'Fire intensity:':<16} {fire.intensity/1e6:>7.2f} MW/m²"),
        mlines.Line2D([], [], color='none', label=f"{'Fire width:':<16} {fire.width:>7.2f} m"),
        mlines.Line2D([], [], color='none', label=f"{'Heat reduction:':<16} {fire.epsilon:>7.2f}"),
        mlines.Line2D([], [], color='none', label=f"{'Tunnel Area:':<16} {tunnel.area:>7.2f} m²"),
        mlines.Line2D([], [], color='none', label=f"{'Tunnel Height:':<16} {tunnel.height:>7.2f} m"),
        mlines.Line2D([], [], color='none', label=f"{'Hydraulic Diam.:':<16} {tunnel.hydraulic_diameter:>7.2f} m"),
        mlines.Line2D([], [], color='none', label=""),
        mlines.Line2D([], [], color='none', label="Never Gray CV Calculator"),
        mlines.Line2D([], [], color='none', label=f"{'MIT License, Version'} {VERSION_NUMBER}"),
    ]
    ax.legend(param_lines, [l.get_label() for l in param_lines], loc='lower right', 
              prop={'family': 'monospace', 'size': 9})

    # Apply common formatting
    format_cv_plot_axes(ax, max_fire_hrr * 1e-6, max_critical_velocity, "Critical Velocity")
    fig.tight_layout()

    if for_web:
        return fig, fire.hrr, critical_velocity, oxygen_depletion_msg, converging_msg
    else:
        # Save files efficiently
        pd.DataFrame(critical_velocities, fire_hrrs*1.0e-6).to_csv(f"{tunnel.name}_critical_velocity.csv")
        fig.savefig(f"{tunnel.name}_critical_velocity.png")
        plt.close(fig)
        return fire.hrr, critical_velocity, oxygen_depletion_msg, converging_msg

def hrrs_vs_critical_velocities(tunnel, fire, ambient_temp, ambient_density, min_hrr_input=0.001e6, max_hrr_input=150e6, for_web=True):
    """Function to create a table of critical velocities against HRR. The function checks if there is sufficient oxygen for the range of HRR.
    """
    # Initial the range
    min_hrr = min(min_hrr_input, max_hrr_input, fire.hrr) #HRR close to zero create high temperatures and lead to a specific heat capacity cut off
    max_hrr = min(3*fire.hrr, max(max_hrr_input, fire.hrr)) #Maximum total fire heat release rate for plot (W)
    acceptable_resolution = 1000 #Minimum number of plot points for acceptable resolution
    insufficient_resolution = True #Initial resolution set to True to start calculation
    converging_msg = ""
    while insufficient_resolution:
        fire_hrrs = np.linspace(min_hrr, max_hrr, 2001)
        critical_velocities = np.empty_like(fire_hrrs)
        critical_velocities[-1] = 2.0 #First Guess
        delta_t = np.empty_like(fire_hrrs)
        sufficient_oxygen = True
        fire_copy = copy.deepcopy(fire) #Creates a copy so the originaly value of fire.hrr is retained.
        for i, fire_hrr in enumerate(fire_hrrs):
            fire_copy.hrr = fire_hrr
            if sufficient_oxygen:
                critical_velocities[i], delta_t[i], iter = iterate_critical_velocity(fire_copy, tunnel, critical_velocities[i-1], ambient_temp, ambient_density)
                vel_depletion = oxygen_depletion(fire_hrr, tunnel, ambient_density)
                #HRR cut off if minimum oxygen requirement not met
                if (critical_velocities[i] <= vel_depletion):
                    print("!!!!HRR cut off due to oxygen depletion!!!!")
                    sufficient_oxygen = False
                    # Slice arrays to keep only valid data less than oxygen depletion
                    fire_hrrs = fire_hrrs[:i+1]
                    critical_velocities = critical_velocities[:i+1]
                    delta_t = delta_t[:i+1]
                    max_hrr = fire_hrrs[i]
                    break
            if iter > 100:
                converging_msg = "WARNING! Critical velocity calculation did not converge. Please check input parameters."
        if len(fire_hrrs) > acceptable_resolution: insufficient_resolution = False
    return fire_hrrs, critical_velocities, sufficient_oxygen, converging_msg

def format_cv_plot_axes(ax, max_hrr_mw, max_cv, title):
    """Common formatting for critical velocity plots"""
    x_top = round_up_nice(max_hrr_mw)
    y_top = round_up_nice(1.05 * max_cv)
    
    ax.xaxis.set_major_locator(MaxNLocator(nbins=10, steps=[1, 1.5, 2, 2.5, 5, 7.5, 10]))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=10, steps=[1, 2, 2.5, 5, 7.5, 10]))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    ax.set(title=title, xlabel="Total Fire Heat Release Rate (MW)",
           xlim=[0.0, x_top], ylabel="Critical Velocity (m/s)", ylim=[0.0, y_top])
    ax.grid(True)
    return ax

def create_comparison_plot(tunnels, fire, ambient_temp, ambient_density):
    """Create comparison plot by combining tunnel data - single call per tunnel"""
    fig, ax = plt.subplots(figsize=(10.0, 6.0))
    
    max_hrr = 0
    max_cv = 0
    
    # Single call per tunnel - get data, plot, and track maximums
    for tunnel in tunnels:
        fire_hrrs, critical_velocities, _, _ = hrrs_vs_critical_velocities(
            tunnel, fire, ambient_temp, ambient_density, for_web=True)
        ax.plot(fire_hrrs*1e-6, critical_velocities, label=tunnel.name, linewidth=2)
        
        # Track maximums during the same loop
        max_hrr = max(max_hrr, max(fire_hrrs) * 1e-6)
        max_cv = max(max_cv, max(critical_velocities))
    
    # Apply common formatting
    format_cv_plot_axes(ax, max_hrr, max_cv, "Critical Velocity Comparison")
    ax.legend(loc='best')
    fig.tight_layout()
    
    fig.savefig("All_Tunnels_critical_velocity_comparison.png", dpi=300)
    plt.close(fig)
    return "All_Tunnels_critical_velocity_comparison.png"

if __name__ == "__main__":
    #Table 8 Typical TBM or Arched tunnel profiles
    TBM1 = Tunnel("TBM-1", 6.80, 85.38, 9.06)
    TBM2 = Tunnel("TBM-2", 8.07, 85.60, 9.28)
    TBM3 = Tunnel("TBM-3", 9.82, 106.0, 10.7)
    TBM4 = Tunnel("TBM-4", 5.00, 40.00, 5.99)
    TBM5 = Tunnel("TBM-5", 5.89, 40.00, 6.44)
    TBM6 = Tunnel("TBM-6", 7.86, 59.95, 8.13)
    Cut_off = Tunnel("Cut-off", 3, 9, 3)
    Tunnels = [TBM1, TBM2, TBM3, TBM4, TBM5, TBM6, Cut_off]
    
    #Fire parameters
    hrr = 50.0e6 #Total fire heat release rate (W)
    intensity = 2.25e6 #Fire intensity (W/m^2)
    width = 2.5 #Maximum width of the fire (m)
    epsilon = 0.2 #Heat reduction due to imperfect combustion and thermal radiation (-)
    eta = 0.0 #Reduction of the convective heat release rate due to water mist systems or sprinklers (-)
    fire = Fire(hrr, intensity, width, epsilon, eta)

    #Ambient and initial parameters
    ambient_temp = 294.0 #Ambient temperature in Kelvin (K)
    ambient_pressure = 101325 #Ambient pressure (Pa)
    critical_velocity = 2.0 #Initial guess of critical velocity to start iteration (m/s)
    
    # Generate individual plots for each tunnel
    for tunnel in Tunnels:
        hrr, critical_velocity, oxygen_depletion_msg, converging_msg = plot_critical_velocity(tunnel, fire, ambient_temp, ambient_pressure, for_web=False)
        print(f"Tunnel : {tunnel.name}, Critical Velocity: {critical_velocity:.3f} m/s, HRR: {hrr*1e-6:.1f} MW")
    
    # Generate comparison plot with all tunnels
    print("\nGenerating comparison plot for all tunnels...")
    comparison_plot_file = create_comparison_plot(Tunnels, fire, ambient_temp, ambient_pressure)
    print(f"Comparison plot saved as: {comparison_plot_file}")

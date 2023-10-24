#!/bin/python
import math
###
# Purpose: This program is designed to calculate the effects of ventilating a building and compares the energy loss to that of a dehumidifier
###
###
# Specify Calculation Variables
###
# Ventilation variables
duration_of_ventilation_in_minutes = 10
number_of_ventilations_per_day = 3
# Building variables
building_area_in_square_meters = 84 # float
building_height_in_meters = 2.5 # float
building_volume_in_cubic_meters = building_area_in_square_meters*building_height_in_meters# float
# Temperature variables
#outside_temperature_in_degrees = 11.2 # TODO USE AS PREDICTION INPUT VALUE
inside_temperature_at_114_centimeters_room_elevation_in_degrees = 22.4
# measured_temperature_after_duration_at_114_centimeters_room_elevation_in_degrees = 20.4
measured_temperature_after_duration_at_114_centimeters_room_elevation_in_degrees = 18.4 # last measured difference was 4 degrees
specific_heat_capacity_of_air_in_joules_kilogram=1005
# Humidity variables
inside_humidity_at_114_centimeters_elevation_in_percentage_relative_humidity = 68
#outside_humidity_in_percentage_relative_humidity = 92 # WARNING THIS SENSOR IS NOT WORKING CORRECTLY
measured_humidity_after_duration_in_percentage_relative_humidity = 54
# Dehumidifier variables
# Floor area: ~80 m^2
# dehumidifier_wh_per_litre_removed=180.0 # https://www.amazon.de/dp/B096SQ1Q5K 30m^2
dehumidifier_wh_per_litre_removed=297.0 # https://www.amazon.de/dp/B096SQ1Q5K 80m^2
# dehumidifier_wh_per_litre_removed=324.0 # https://www.amazon.de/dp/B082FWPWCR 70m^2
# Heater Efficiency
# Efficiency range is est. 85-95%:
# hydronic_radiant_efficiency = 0.85 
hydronic_radiant_efficiency = 0.95 

###
# Define functions
###
# Returns the energy lost due to loss of heat energy in joules
def calc_thermal_energy_loss_in_joules(specific_heat_cap_joule_kilogram: float=1005.0, volume_in_cubic_meters: float=0, temperature_difference: float=0):
    energy_loss_in_joules = specific_heat_cap_joule_kilogram * volume_in_cubic_meters * temperature_difference
    return energy_loss_in_joules

# Converts joules to watts
def calc_joules_to_watt_hours(joules):
    return joules/3600

def calc_fast_absolute_humidity_in_gram_cubic_meter(rh,temperature):
    """
    Calculate the fast approximation of absolute humidity using Antoine equation.
    
    Parameters:
        rh (float): Relative Humidity in percentage
        temperature (float): Temperature in Celsius
    
    Returns:
        float: Absolute humidity in g/m^3
    """
    # Antoine equation constants for water
    A_water = 8.07131
    B_water = 1730.63
    C_water = 233.426
    # Calculate vapor pressure using Antoine equation (in mmHg)
    P_mmHg = 10**(A_water - (B_water / (temperature + C_water)))
    # Convert vapor pressure to Pascals
    P_Pa = P_mmHg * 133.322
    # Calculate the partial pressure of water vapor (in Pa)
    partial_pressure = rh / 100 * P_Pa
    # Convert temperature to Kelvin
    temperature_K = temperature + 273.15
    # Universal Gas constant (R) in J/(mol*K)
    R = 8.314
    # Molar mass of water in kg/mol
    M_water = 18.01528 / 1000
    # Calculate absolute humidity in g/m^3
    absolute_humidity = (partial_pressure * M_water) / (R * temperature_K) * 1e3
    return absolute_humidity

def calc_accurate_absolute_humidity_in_gram_cubic_meter(rh,temperature):
    """
    Calculate the accurate absolute humidity using Wagner and Pruss equation.
    
    Parameters:
        rh (float): Relative Humidity in percentage
        temperature (float): Temperature in Celsius
    
    Returns:
        float: Absolute humidity in g/m^3
    """
    # Specific gas constant for water vapor in J/(kg*K)
    Rw = 461.5
    # Critical pressure for water in Pa
    Pc = 22.064e6
    # Critical temperature for water in K
    Tc = 647.096
    # Empirical constants for Wagner and Pruss equation
    a1 = -7.85951783
    a2 = 1.84408259
    a3 = -11.7866497
    a4 = 22.6807411
    a5 = -15.9618719
    a6 = 1.80122502
    # Temperature in Kelvin
    T = temperature + 273.15
    # tau for Wagner and Pruss equation
    tau = 1 - T / Tc
    # Saturation vapor pressure using Wagner and Pruss equation in Pa
    Ps = Pc * math.exp(Tc / T * (a1 * tau + a2 * tau**1.5 + a3 * tau**3 + a4 * tau**3.5 + a5 * tau**4 + a6 * tau**7.5))
    # Calculate absolute humidity in g/m^3
    AH = rh * Ps / (Rw * T) * 10
    return AH

def rh_to_ah_in_liters(rh, temperature, volume):
    """
    Convert relative humidity to absolute humidity in liters for a given volume.
    
    Parameters:
        rh (float): Relative Humidity in percentage
        temperature (float): Temperature in Celsius
        volume (float): Volume in cubic meters
    
    Returns:
        float: Absolute humidity in liters
    """
    # Calculate accurate absolute humidity in g/m^3
    ah_gram_per_cubic_meter = calc_accurate_absolute_humidity_in_gram_cubic_meter(rh, temperature)
    # Calculate the total water content in the given volume in grams
    total_water_content_grams = ah_gram_per_cubic_meter * volume
    # Convert the total water content to liters (since 1 liter of water is 1000 grams)
    total_water_content_liters = total_water_content_grams / 1000
    return total_water_content_liters

#  Calculates various effects of velntilating a building 
def simulate_and_compare(duration: float=0.0, 
                volume: float=0.0, 
                inside_temp: float=0.0,
                measured_temp: float=0.0,
                inside_hum_rh: float=0,
                measured_hum_rh: float=0,
                ):
    heat_loss_wh=0
    hum_reduction=0
    # Volume
    print(f"-- ENERGY LOSS CALCULATION --\n\n\tRoom volume:\t\t\t{volume} m^3")
    # Calculate Energy loss due to temperature difference for building volume
    if measured_temp<inside_temp: # If temperature fell during ventilation
        temperature_difference=(inside_temp-measured_temp)*number_of_ventilations_per_day
        heat_loss_in_jules=calc_thermal_energy_loss_in_joules(specific_heat_capacity_of_air_in_joules_kilogram,
                                                                building_volume_in_cubic_meters,
                                                                temperature_difference)
        heat_loss_wh=calc_joules_to_watt_hours(heat_loss_in_jules)
        print(f"\n\tHeat loss({temperature_difference:.1f} C):\t\t{heat_loss_in_jules:.2f} Joules\n\t\t\t\t\t\t\t{heat_loss_wh:.2f} Wh")
        # Compensate for hydronic radiant floor heating energy conversion losses
        heat_loss_wh = heat_loss_wh / hydronic_radiant_efficiency
        print(f"\n\tHeat loss({temperature_difference:.1f} C) after energy conversion compensation ({hydronic_radiant_efficiency} efficiency):\t\t{heat_loss_wh:.2f} Wh")
    else:
        print(f"\n\tHeat loss:\t\t\t\tNONE") 
    # Calculate humidity difference
    if measured_hum_rh<inside_hum_rh: # If humidity fell
        humidity_difference=inside_hum_rh-measured_hum_rh
        print(f"\n\tHumidity reduction:\t\t{humidity_difference}%\trH")
        hum_reduction=humidity_difference
    else:
        print(f"\n\tHumidity reduction:\t\tNONE")
    # Compare energy loss to energy consumption of dehumidifier
    if hum_reduction>0 and heat_loss_wh>0:
        # Project energy loss
        litres_of_water_removed=rh_to_ah_in_liters(hum_reduction,inside_temperature_at_114_centimeters_room_elevation_in_degrees,volume)
        print(f"\n-- COMPARISON --\n\nComparing airing to dehumidifier for {hum_reduction} % rH in {volume} m^3 reduction (litres removed: {litres_of_water_removed:.3f} L):")
        projected_dehumidifier_energy_loss_wh=dehumidifier_wh_per_litre_removed*litres_of_water_removed
        print(f"\nProjected dehumidifier consumption:\t{projected_dehumidifier_energy_loss_wh:.3f} Wh")
        # Compare efficiency
        energy_loss_delta_ventilation_to_dehum=heat_loss_wh-projected_dehumidifier_energy_loss_wh
        print("\n-- VERDICT --\n")
        if energy_loss_delta_ventilation_to_dehum>0:
            print(f"Dehumidifier consumes {abs(energy_loss_delta_ventilation_to_dehum):.3f} Wh less")
        else:
            print(f"Dehumidifier consumes {abs(energy_loss_delta_ventilation_to_dehum):.3f} Wh more")
###
# Run calculation
###
simulate_and_compare(duration_of_ventilation_in_minutes,
            building_volume_in_cubic_meters,
            inside_temperature_at_114_centimeters_room_elevation_in_degrees,
            measured_temperature_after_duration_at_114_centimeters_room_elevation_in_degrees,
            inside_humidity_at_114_centimeters_elevation_in_percentage_relative_humidity,
            measured_humidity_after_duration_in_percentage_relative_humidity)

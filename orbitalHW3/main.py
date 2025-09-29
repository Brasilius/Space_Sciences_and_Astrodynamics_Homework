import math
import numpy as np

def main():
    print("Hello from orbitalhw3!")


    print("Problem 1:  spacecraft is in a circular orbit of Mars at an altitude of 200 km. Calculate its speed and its period. ")
    altitude = 200000 #noted as M
    mars_radius = 3.4e6 #noted as M
    mars_gravity = 6.67e-11 #noted as m/s^2
    mars_mass = 6.417e23 #noted as kg
    orbital_speed = math.sqrt((mars_gravity*mars_mass)/(mars_radius+altitude))
    orbital_period = (2*math.pi*math.sqrt((mars_radius+altitude)**3/(mars_gravity*mars_mass)))
    print("The orbital speed is: ", orbital_speed, "m/s")
    print("The orbital period is: ", orbital_period, "s")
    #conclude problem 1
    #start problem 2
    print("Problem 2: A spacecraft is in a 400-km-by-600-km low Earth orbit. How long (in minutes) does it take to coast from the perigee to the apogee? ")
    earth_radius = 6.371e6 #noted as M
    earth_gravity = 6.67e-11 #noted as m/s^
    earth_mass = 5.972e24 #noted as kg
    perigee = 400000 + earth_radius #noted as m
    apogee = 600000 + earth_radius #noted as m
    semi_major_axis = (perigee + apogee)/2 #noted as
    orbital_period_elliptical = (2*math.pi*math.sqrt((semi_major_axis)**3/(earth_gravity*earth_mass)))
    time_to_coast = orbital_period_elliptical/2
    time_to_coast_minutes = time_to_coast/60
    print("The time to coast from perigee to apogee is: ", time_to_coast_minutes, "minutes")
    #conclude problem 2
    #start problem 3
    print("he altitude of a satellite in an elliptical orbit around the Earth is 2000 km at the apogee and 500 km at the perigee. Determine the eccentricity of the orbit, the orbital speeds at perigee and apogee and the period of the orbit")
    perigee_3 = 500000 + earth_radius #noted as m
    apogee_3 = 2000000 + earth_radius #noted as m
    semi_major_axis_3 = (perigee_3 + apogee_3)/2 #noted as m
    eccentricity = (apogee_3 - perigee_3)/(apogee_3 + perigee_3)
    orbital_speed_perigee = math.sqrt(earth_gravity*earth_mass*((2/perigee_3)-(1/semi_major_axis_3)))
    orbital_speed_apogee = math.sqrt(earth_gravity*earth_mass*((2/apogee_3)-(1/semi_major_axis_3)))
    orbital_period_3 = (2*math.pi*math.sqrt((semi_major_axis_3)**3/(earth_gravity*earth_mass)))
    print("The eccentricity of the orbit is: ", eccentricity)
    print("The orbital speed at perigee is: ", orbital_speed_perigee, "m/s")
    print("The orbital speed at apogee is: ", orbital_speed_apogee, "m/s")
    print("The orbital period is: ", orbital_period_3, "s")
    #conclude problem 3





if __name__ == "__main__":
    main()

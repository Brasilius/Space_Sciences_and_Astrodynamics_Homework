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

    print("Problem 2: A spacecraft is in a 400-km-by-600-km low Earth orbit. How long (in minutes) does it take to coast from the perigee to the apogee? ")
    
if __name__ == "__main__":
    main()

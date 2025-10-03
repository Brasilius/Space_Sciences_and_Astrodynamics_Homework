import math


def main() -> None:
    # Problem 1
    print("Problem 1:  spacecraft is in a circular orbit of Mars at an altitude of 200 km. Calculate its speed and its period. ")
    altitude1 = 200000  # noted as m
    mu_mars = 4.282837e13
    mars_radius = 3396000  # noted as m
    orbital_speed1 = math.sqrt(mu_mars / (mars_radius + altitude1))
    orbital_period1 = 2 * math.pi * math.sqrt((mars_radius + altitude1) ** 3 / mu_mars)
    print("The orbital speed is: ", orbital_speed1, "m/s")
    print("The orbital period is: ", orbital_period1, "s")

    # Problem 2
    print("Problem 2: A spacecraft is in a 400-km-by-600-km low Earth orbit. How long (in minutes) does it take to coast from the perigee to the apogee? ")
    mu_earth = 3.986004418e14
    earth_radius = 6378137  # noted as m
    perigee2 = 400000 + earth_radius  # noted as m
    apogee2 = 600000 + earth_radius  # noted as m
    semi_major2 = (perigee2 + apogee2) / 2
    orbital_period2 = 2 * math.pi * math.sqrt(semi_major2 ** 3 / mu_earth)
    time_to_coast2 = orbital_period2 / 2
    print("The time to coast from perigee to apogee is: ", time_to_coast2 / 60, "minutes")

    # Problem 3
    print("Problem 3: The altitude of a satellite in an elliptical orbit around the Earth is 2000 km at the apogee and 500 km at the perigee. Determine the eccentricity of the orbit, the orbital speeds at perigee and apogee and the period of the orbit")
    perigee3 = 500000 + earth_radius  # noted as m
    apogee3 = 2000000 + earth_radius  # noted as m
    semi_major3 = (perigee3 + apogee3) / 2
    eccentricity3 = (apogee3 - perigee3) / (apogee3 + perigee3)
    orbital_speed_perigee3 = math.sqrt(mu_earth * (2 / perigee3 - 1 / semi_major3))
    orbital_speed_apogee3 = math.sqrt(mu_earth * (2 / apogee3 - 1 / semi_major3))
    orbital_period3 = 2 * math.pi * math.sqrt(semi_major3 ** 3 / mu_earth)
    print("The eccentricity of the orbit is: ", eccentricity3)
    print("The orbital speed at perigee is: ", orbital_speed_perigee3, "m/s")
    print("The orbital speed at apogee is: ", orbital_speed_apogee3, "m/s")
    print("The orbital period is: ", orbital_period3, "s")

    # Problem 4
    print("Problem 4: An earth satellite has a speed of 7.5 km/s and a flight path angle of 10 degrees when its radius is 8000 km")
    v4 = 7500  # m/s
    phi4 = math.radians(10)  # noted as degrees converted to radians
    r4 = 8000000  # noted as m
    eps4 = v4 ** 2 / 2 - mu_earth / r4
    a4 = -mu_earth / (2 * eps4)
    h4 = r4 * v4 * math.cos(phi4)
    e4 = math.sqrt(1 + 2 * eps4 * h4 ** 2 / mu_earth ** 2)
    cos_nu4 = (h4 ** 2 / (mu_earth * r4) - 1) / e4
    cos_nu4 = max(-1.0, min(1.0, cos_nu4))
    nu4 = math.degrees(math.acos(cos_nu4))
    perigee4 = a4 * (1 - e4)
    apogee4 = a4 * (1 + e4)
    print("The semi-major axis is: ", a4 / 1000, "km")
    print("The eccentricity is: ", e4)
    print("The true anomaly is: ", nu4, "degrees")
    print("The perigee is: ", perigee4 / 1000, "km")
    print("The apogee is: ", apogee4 / 1000, "km")

    # Problem 5
    print("Problem 5: For an Earth satellite, the specific angular momentum is 70,000 km^2/s and the specific energy is −10 km^2/s^2. Calculate the apogee and perigee altitudes.")
    h5 = 70000  # km^2/s
    eps5 = -10  # km^2/s^2
    mu_earth_km = mu_earth / 1e9
    semi_major5 = -mu_earth_km / (2 * eps5)
    e5 = math.sqrt(1 + 2 * eps5 * h5 ** 2 / mu_earth_km ** 2)
    perigee5 = semi_major5 * (1 - e5)
    apogee5 = semi_major5 * (1 + e5)
    perigee_alt5 = perigee5 - earth_radius / 1000
    apogee_alt5 = apogee5 - earth_radius / 1000
    print("The perigee altitude is: ", perigee_alt5, "km")
    print("The apogee altitude is: ", apogee_alt5, "km")

    # Problem 6
    print("Problem 6: A hyperbolic Earth departure trajectory has a perigee altitude of 250 km and a perigee speed of 11 km/s.")
    r_p6 = earth_radius / 1000 + 250  # km
    v_p6 = 11  # km/s
    eps6 = v_p6 ** 2 / 2 - mu_earth_km / r_p6
    v_inf6 = math.sqrt(2 * eps6)
    h6 = r_p6 * v_p6
    e6 = math.sqrt(1 + 2 * eps6 * h6 ** 2 / mu_earth_km ** 2)
    a6 = -mu_earth_km / (2 * eps6)
    nu6 = math.radians(100)
    r6 = h6 ** 2 / (mu_earth_km * (1 + e6 * math.cos(nu6)))
    v_r6 = mu_earth_km / h6 * e6 * math.sin(nu6)
    v_t6 = h6 / r6
    print("The hyperbolic excess speed is: ", v_inf6, "km/s")
    print("The radius at true anomaly 100 degrees is: ", r6, "km")
    print("The radial velocity component is: ", v_r6, "km/s")
    print("The transverse velocity component is: ", v_t6, "km/s")

    # Problem 7
    print("Problem 7: A meteoroid is first observed approaching Earth when it is 402,000 km from the center of the Earth with a true anomaly of 150 degrees. If the speed of the meteoroid at that time is 2.23 km/s, calculate the eccentricity of the trajectory, the altitude at closest approach and the speed at the closest approach.")
    r_obs7 = 402000  # km
    v_obs7 = 2.23  # km/s
    nu_obs7 = math.radians(150)
    eps7 = v_obs7 ** 2 / 2 - mu_earth_km / r_obs7
    a7 = -mu_earth_km / (2 * eps7)
    A7 = a7
    B7 = r_obs7 * math.cos(nu_obs7)
    C7 = -(a7 - r_obs7)
    disc7 = B7 ** 2 - 4 * A7 * C7
    e7_1 = (-B7 + math.sqrt(disc7)) / (2 * A7)
    e7_2 = (-B7 - math.sqrt(disc7)) / (2 * A7)
    e7 = e7_1 if e7_1 >= 1 else e7_2
    h7 = math.sqrt(r_obs7 * mu_earth_km * (1 + e7 * math.cos(nu_obs7)))
    perigee7 = a7 * (1 - e7)
    perigee_alt7 = perigee7 - earth_radius / 1000
    v_per7 = math.sqrt(mu_earth_km * (2 / perigee7 - 1 / a7))
    print("The eccentricity is: ", e7)
    print("The altitude at closest approach is: ", perigee_alt7, "km")
    print("The speed at closest approach is: ", v_per7, "km/s")


if __name__ == "__main__":
    main()

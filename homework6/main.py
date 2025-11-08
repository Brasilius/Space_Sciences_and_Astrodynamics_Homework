import numpy as np
import math

def main():
    print("Hello from homework6!")
    # start problem 1
    r1i = 2500
    r1j = 16000
    r1z = 4000
    v1i = -3
    v1j = -1
    v1z = 5

    # vectors
    r_vec = np.array([r1i, r1j, r1z], dtype=float)
    v_vec = np.array([v1i, v1j, v1z], dtype=float)

    # magnitudes
    r1_magnitude = np.linalg.norm(r_vec)
    v1_magnitude = np.linalg.norm(v_vec)

    # radial velocity
    rv_dot = np.dot(r_vec, v_vec)
    radial_velocity = rv_dot / r1_magnitude

    print("problem 1 solutions: ")

    # specific angular momentum
    h = np.cross(r_vec, v_vec)
    h_magnitude = np.linalg.norm(h)
    print("specific angular momentum vector h:", h, "km^2/s")
    print("specific angular momentum magnitude |h|:", h_magnitude, "km^2/s")

    # inclination i
    hz = h[2]
    i = math.acos(hz / h_magnitude)
    i_deg = math.degrees(i)
    print("inclination i:", i_deg, "deg")

    # node vector N = k × h
    k_vec = np.array([0.0, 0.0, 1.0])
    N = np.cross(k_vec, h)
    N_mag = np.linalg.norm(N)
    print("node vector N:", N)
    print("|N|:", N_mag)

    # RAAN Ω
    if N_mag != 0:
        Omega = math.degrees(math.acos(N[0] / N_mag))
        if N[1] < 0:      # quadrant check
            Omega = 360.0 - Omega
    else:
        Omega = 0.0
    print("RAAN Ω:", Omega, "deg")

    # gravitational parameter μ
    mu = 398600.0

    # eccentricity vector using:
    # e = (1/μ)[(v^2 - μ/r) r - (r·v) v]
    term1 = (v1_magnitude**2 - mu / r1_magnitude) * r_vec
    term2 = rv_dot * v_vec
    e_vec = (term1 - term2) / mu
    emag = np.linalg.norm(e_vec)
    print("eccentricity vector e:", e_vec)
    print("eccentricity magnitude e:", emag)

    # argument of perigee ω (angle from N to e)
    if emag != 0 and N_mag != 0:
        cos_omega = np.dot(N, e_vec) / (N_mag * emag)
        # numerical safety
        cos_omega = max(min(cos_omega, 1.0), -1.0)
        omega = math.degrees(math.acos(cos_omega))
        if e_vec[2] < 0:   # quadrant check using e_z
            omega = 360.0 - omega
    else:
        omega = 0.0
    print("argument of perigee ω:", omega, "deg")

    # true anomaly θ (angle from e to r)
    if emag != 0:
        cos_theta = np.dot(e_vec, r_vec) / (emag * r1_magnitude)
        # numerical safety
        cos_theta = max(min(cos_theta, 1.0), -1.0)
        theta = math.degrees(math.acos(cos_theta))
        if radial_velocity < 0:   # quadrant check using r·v
            theta = 360.0 - theta
    else:
        theta = 0.0
    print("true anomaly θ:", theta, "deg")

    # semi-major axis a from vis-viva: 1/a = 2/r - v^2/μ
    a = 1.0 / (2.0 / r1_magnitude - v1_magnitude**2 / mu)
    print("semi-major axis a:", a, "km")


    print("Problem 2 solutions:")

    # Given state
    r_vec = np.array([0.0, 0.0, -13000.0])
    v_vec = np.array([4.0, 5.0, 6.0])

    mu = 398600.0

    # Magnitudes
    r_mag = np.linalg.norm(r_vec)
    v_mag = np.linalg.norm(v_vec)

    # Radial velocity
    rv_dot = np.dot(r_vec, v_vec)
    radial_velocity = rv_dot / r_mag

    print("r magnitude:", r_mag, "km")
    print("v magnitude:", v_mag, "km/s")
    print("radial velocity:", radial_velocity, "km/s")

    # Specific angular momentum
    h = np.cross(r_vec, v_vec)
    h_mag = np.linalg.norm(h)
    print("h vector:", h, "km^2/s")
    print("|h|:", h_mag, "km^2/s")

    # Inclination i
    hz = h[2]
    i = math.acos(hz / h_mag)
    i_deg = math.degrees(i)
    print("inclination i:", i_deg, "deg")

    # Node vector N = k × h
    k_vec = np.array([0.0, 0.0, 1.0])
    N = np.cross(k_vec, h)
    N_mag = np.linalg.norm(N)
    print("node vector N:", N)
    print("|N|:", N_mag)

    # RAAN Ω
    if N_mag != 0:
        Omega = math.degrees(math.acos(N[0] / N_mag))
        if N[1] < 0:
            Omega = 360.0 - Omega
    else:
        Omega = 0.0
    print("RAAN Ω:", Omega, "deg")

    # Eccentricity vector:
    # e = (1/μ)[(v^2 - μ/r) r - (r·v) v]
    term1 = (v_mag**2 - mu / r_mag) * r_vec
    term2 = rv_dot * v_vec
    e_vec = (term1 - term2) / mu
    e_mag = np.linalg.norm(e_vec)
    print("eccentricity vector e:", e_vec)
    print("eccentricity magnitude e:", e_mag)

    # Argument of perigee ω (angle from N to e)
    if e_mag != 0 and N_mag != 0:
        cos_omega = np.dot(N, e_vec) / (N_mag * e_mag)
        cos_omega = max(min(cos_omega, 1.0), -1.0)  # clamp
        omega = math.degrees(math.acos(cos_omega))
        if e_vec[2] < 0:
            omega = 360.0 - omega
    else:
        omega = 0.0
    print("argument of perigee ω:", omega, "deg")

    # True anomaly θ (angle from e to r)
    if e_mag != 0:
        cos_theta = np.dot(e_vec, r_vec) / (e_mag * r_mag)
        cos_theta = max(min(cos_theta, 1.0), -1.0)
        theta = math.degrees(math.acos(cos_theta))
        if rv_dot < 0:
            theta = 360.0 - theta
    else:
        theta = 0.0
    print("true anomaly θ:", theta, "deg")

    # Semi-major axis a (will be negative for hyperbolic)
    a = 1.0 / (2.0 / r_mag - v_mag**2 / mu)
    print("semi-major axis a:", a, "km")

    # Now hyperbolic eccentric anomaly F and hyperbolic mean anomaly Mh
    # (since e > 1 we use hyperbolic formulas)
    nu = math.radians(theta)

    if e_mag > 1.0:
        # tanh(F/2) = sqrt((e - 1)/(e + 1)) * tan(nu/2)
        A = math.sqrt((e_mag - 1.0) / (e_mag + 1.0))
        tanhF2 = A * math.tan(nu / 2.0)
        F = 2.0 * math.atanh(tanhF2)

        # Mh = e*sinh(F) - F
        Mh = e_mag * math.sinh(F) - F

        print("hyperbolic eccentric anomaly F:", F, "rad")
        print("hyperbolic mean anomaly Mh:", Mh, "rad")
    else:
        # (Not used here, but for elliptical orbits:)
        # E = 2*atan( sqrt((1-e)/(1+e)) * tan(nu/2) )
        # Mh = E - e*sin(E)
        print("Orbit is not hyperbolic; use elliptical formulas here.")
    # Given vectors
    r = np.array([-6000.0, -1000.0, -5000.0])  # km
    e = np.array([0.4, 0.5, 0.6])

    # Magnitudes
    r_mag = np.linalg.norm(r)
    e_mag = np.linalg.norm(e)

    # Compute true anomaly from e·r = e r cos(theta)
    cos_theta = np.dot(e, r) / (e_mag * r_mag)

    # Numerical safety: clamp to [-1, 1]
    cos_theta = max(min(cos_theta, 1.0), -1.0)

    # Principal value in [0, π]
    theta = np.arccos(cos_theta)

    # Satellite is approaching perigee ⇒ use angle in (π, 2π)
    theta = 2.0 * np.pi - theta

    # Convert to degrees
    theta_deg = np.degrees(theta)
    print("True anomaly θ:", theta_deg, "deg")
    # Given vectors
    r = np.array([-6600.0, -1300.0, -5200.0])  # km
    e = np.array([-0.4, -0.5, -0.6])

    # The orbit's angular momentum vector is perpendicular to both r and e
    # h = r × e  (direction of orbital angular momentum)
    h = np.cross(r, e)

    # Inclination: angle between h and z-axis (k-hat)
    hz = h[2]
    h_mag = np.linalg.norm(h)

    i = np.arccos(hz / h_mag)
    i_deg = np.degrees(i)
    print("Inclination i:", i_deg, "deg")


if __name__ == "__main__":
    main()

import math
from typing import Tuple

def earth_escape_delta_v() -> Tuple[float, float, float]:
    mu_sun = 1.32712440018e11
    r = 1.495978707e8
    v_circ = math.sqrt(mu_sun / r)
    v_escape = math.sqrt(2 * mu_sun / r)
    delta_v = v_escape - v_circ
    return v_circ, v_escape, delta_v

def hyperbolic_eccentricity(r_p: float, v_inf: float, mu: float) -> float:
    return 1.0 + (r_p * v_inf * v_inf) / mu

def true_anomaly_from_state(r: float, v: float, angle_rad: float, mu: float) -> float:
    v_r = v * math.cos(angle_rad)
    v_theta = v * math.sin(angle_rad)
    h = r * v_theta
    epsilon = v * v / 2.0 - mu / r
    e = math.sqrt(1.0 + 2.0 * epsilon * h * h / (mu * mu))
    p = h * h / mu
    cos_nu = (p / r - 1.0) / e
    cos_nu = max(min(cos_nu, 1.0), -1.0)
    nu = math.acos(cos_nu)
    sin_nu = (v_r * h) / (e * mu)
    if sin_nu < 0:
        nu = -nu
    return math.degrees(nu)

def propagate_true_anomaly(
    r_vec: Tuple[float, float, float],
    v_vec: Tuple[float, float, float],
    delta_nu_deg: float,
    mu: float = 398600.0,
) -> Tuple[Tuple[float, float, float], float]:
    import numpy as np
    r0 = np.array(r_vec, dtype=float)
    v0 = np.array(v_vec, dtype=float)
    r = np.linalg.norm(r0)
    v = np.linalg.norm(v0)
    h_vec = np.cross(r0, v0)
    h = np.linalg.norm(h_vec)
    k_vec = np.array([0.0, 0.0, 1.0])
    n_vec = np.cross(k_vec, h_vec)
    n = np.linalg.norm(n_vec)
    e_vec = (np.cross(v0, h_vec) / mu) - (r0 / r)
    e = np.linalg.norm(e_vec)
    energy = v * v / 2.0 - mu / r
    if abs(energy) > 1e-10:
        a = -mu / (2.0 * energy)
    else:
        a = float('inf')
    p = h * h / mu
    inc = math.degrees(math.acos(h_vec[2] / h))
    if n != 0:
        RAAN = math.degrees(math.acos(n_vec[0] / n))
        if n_vec[1] < 0.0:
            RAAN = 360.0 - RAAN
    else:
        RAAN = 0.0
    if n != 0 and e != 0:
        argp = math.degrees(math.acos(np.dot(n_vec, e_vec) / (n * e)))
        if e_vec[2] < 0.0:
            argp = 360.0 - argp
    else:
        argp = math.degrees(math.atan2(e_vec[1], e_vec[0]))
        if argp < 0.0:
            argp += 360.0
    if e != 0:
        nu0 = math.degrees(math.acos(np.dot(e_vec, r0) / (e * r)))
        if np.dot(r0, v0) < 0.0:
            nu0 = 360.0 - nu0
    else:
        nu0 = math.degrees(math.atan2(r0[2] / math.sin(math.radians(inc)), r0[0]))
    nu1 = math.radians(nu0 + delta_nu_deg)
    r1 = p / (1.0 + e * math.cos(nu1))
    r_pf = np.array([r1 * math.cos(nu1), r1 * math.sin(nu1), 0.0])
    O = math.radians(RAAN)
    i = math.radians(inc)
    w = math.radians(argp)
    R3_O = np.array([[math.cos(O), -math.sin(O), 0.0],
                     [math.sin(O),  math.cos(O), 0.0],
                     [0.0,          0.0,         1.0]])
    R1_i = np.array([[1.0, 0.0, 0.0],
                     [0.0, math.cos(i), -math.sin(i)],
                     [0.0, math.sin(i),  math.cos(i)]])
    R3_w = np.array([[math.cos(w), -math.sin(w), 0.0],
                     [math.sin(w),  math.cos(w), 0.0],
                     [0.0,          0.0,         1.0]])
    Q = R3_O @ R1_i @ R3_w
    r_new = Q @ r_pf
    return tuple(r_new), nu0

def ballistic_descent(
    altitude1: float,
    velocity1: float,
    gamma1_deg: float,
    altitude2: float,
    mu: float = 398600.0,
    R_body: float = 6378.0,
) -> Tuple[float, float]:
    r1 = R_body + altitude1
    r2 = R_body + altitude2
    gamma1 = math.radians(gamma1_deg)
    v_r1 = velocity1 * math.sin(gamma1)
    v_theta1 = velocity1 * math.cos(gamma1)
    h = r1 * v_theta1
    energy1 = velocity1 * velocity1 / 2.0 - mu / r1
    v2 = math.sqrt(2.0 * (energy1 + mu / r2))
    cos_gamma2 = h / (r2 * v2)
    cos_gamma2 = max(min(cos_gamma2, 1.0), -1.0)
    gamma2 = math.acos(cos_gamma2)
    gamma2_deg = -math.degrees(gamma2)
    return v2, gamma2_deg

def classify_conic_from_energy(E: float) -> str:
    if abs(E) < 1e-10:
        return 'parabolic'
    elif E < 0.0:
        return 'elliptical'
    else:
        return 'hyperbolic'

def classify_orbits() -> dict:
    result = {}
    mu = 1.0
    r_a = 3.0
    v_a = 1.5
    E_a = v_a * v_a / 2.0 - mu / r_a
    type_a = classify_conic_from_energy(E_a)
    result['a'] = type_a
    r_p = 1.5
    p_b = 3.0
    e_b = p_b / r_p - 1.0
    if abs(e_b - 1.0) < 1e-10:
        type_b = 'parabolic'
    elif e_b < 1.0:
        type_b = 'elliptical'
    else:
        type_b = 'hyperbolic'
    result['b'] = type_b
    E_c = -1.0 / 3.0
    p_c = 1.5
    a_c = -mu / (2.0 * E_c)
    e_sq_c = max(0.0, 1.0 - p_c / a_c)
    e_c = math.sqrt(e_sq_c)
    if e_c < 1e-10:
        type_c = 'circular'
    else:
        type_c = 'elliptical'
    result['c'] = type_c
    import numpy as np
    r_d = np.array([0.0, 1.0, 0.2])
    v_d = np.array([0.9, 0.0, 0.123])
    r_norm_d = np.linalg.norm(r_d)
    v_norm_d = np.linalg.norm(v_d)
    E_d = v_norm_d * v_norm_d / 2.0 - mu / r_norm_d
    e_d = math.sqrt(1.0 + 2.0 * E_d * (np.linalg.norm(np.cross(r_d, v_d)) ** 2) / (mu * mu))
    if E_d < 0.0:
        if abs(e_d) < 1e-10:
            type_d = 'circular'
        else:
            type_d = 'elliptical'
    elif abs(E_d) < 1e-10:
        type_d = 'parabolic'
    else:
        type_d = 'hyperbolic'
    result['d'] = type_d
    r_e = np.array([0.0, 0.0, 1.01])
    v_e = np.array([1.0, 0.0, 1.4])
    r_norm_e = np.linalg.norm(r_e)
    v_norm_e = np.linalg.norm(v_e)
    E_e = v_norm_e * v_norm_e / 2.0 - mu / r_norm_e
    e_e = math.sqrt(1.0 + 2.0 * E_e * (np.linalg.norm(np.cross(r_e, v_e)) ** 2) / (mu * mu))
    if E_e < 0.0:
        if abs(e_e) < 1e-10:
            type_e = 'circular'
        else:
            type_e = 'elliptical'
    elif abs(E_e) < 1e-10:
        type_e = 'parabolic'
    else:
        type_e = 'hyperbolic'
    result['e'] = type_e
    return result

def main() -> None:
    v_circ, v_esc, delta_v = earth_escape_delta_v()
    print('Problem 1:')
    print(f'  Earth orbital speed at 1 AU = {v_circ:.2f} km/s')
    print(f'  Parabolic escape speed at 1 AU = {v_esc:.2f} km/s')
    print(f'  Velocity increase required = {delta_v:.2f} km/s\n')
    e_demo = hyperbolic_eccentricity(10000.0, 3.0, 398600.0)
    print('Problem 2:')
    print(f'  Example eccentricity for r_p=10,000 km, v_inf=3 km/s: e = {e_demo:.3f}\n')
    mu_earth = 398600.0
    r3 = 10000.0
    v3 = 10.0
    nu3 = true_anomaly_from_state(r3, v3, math.radians(120.0), mu_earth)
    print('Problem 3:')
    print(f'  True anomaly = {nu3:.2f}°')
    print('  (Values >360° or negative can be interpreted modulo 360°)\n')
    r0 = (7000.0, 0.0, 0.0)
    v0 = (7.0, 7.0, 0.0)
    r_after, nu0 = propagate_true_anomaly(r0, v0, 90.0)
    print('Problem 4:')
    print(f'  Initial true anomaly = {nu0:.2f}°')
    print('  Position vector after +90°:')
    print(f'    x = {r_after[0]:.1f} km')
    print(f'    y = {r_after[1]:.1f} km')
    print(f'    z = {r_after[2]:.1f} km\n')
    v2, gamma2 = ballistic_descent(100.0, 7.62, -60.0, 200.0)
    print('Problem 5:')
    print(f'  Speed at 200 km = {v2:.2f} km/s')
    print(f'  Flight-path angle at 200 km = {gamma2:.2f}°\n')
    types = classify_orbits()
    print('Problem 6:')
    for part, typ in types.items():
        print(f'  Part {part}: {typ} orbit')

if __name__ == '__main__':
    main()

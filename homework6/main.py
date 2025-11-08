import numpy as np
import math

def main():
    print("Hello from homework6!")
    #start problem 1
    r1i = 2500
    r1j = 16000
    r1z = 4000
    v1i = -3
    v1j = -1
    v1z = 5
    #magnitude of r1
    r1_magnitude = math.sqrt(r1i**2 + r1j**2 + r1z**2)
    #magnitude of v1
    v1_magnitude = math.sqrt(v1i**2 + v1j**2 + v1z**2)
    radial_velocity = (r1i*v1i + r1j*v1j + r1z*v1z)/r1_magnitude
    print("problem 1 solutions: ")
    h = np.cross([r1i, r1j, r1z], [v1i, v1j, v1z])
    print("specific angular momentum vector: ", h)
    h_magnitude = np.linalg.norm(h)
    print("specific angular momentum magnitude: ", h_magnitude)
    hz = h[2]
    i = np.arccos(hz/h_magnitude)
    i2 = np.rad2deg(i)
    print("angle of orbit: ", i2)
    N = np.cross([0,0,1],[h[0],h[1],h[2]])
    print(N)
    N_mag = np.linalg.norm(N)
    print("Nmag:", N_mag)
    omega = np.arccos(N[0]/N_mag)
    omega = np.rad2deg(omega)
    print("omega is: " , omega , "deg")
    mew = 398600
    ei = ((1/mew) * (v1_magnitude**2 - mew/r1_magnitude) * v1i)
    ej = ((1/mew) * (v1_magnitude**2 - mew/r1_magnitude) * v1j)
    ez = ((1/mew) * (v1_magnitude**2 - mew/r1_magnitude) * v1z)                  
    print(ei, ej, ez)
    emag = math.sqrt(ei**2 + ej**2 + ez**2)
    print(emag)
    
if __name__ == "__main__":
    main()

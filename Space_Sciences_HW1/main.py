import numpy as np


def main():
    print("Hello from space-sciences-hw1!")
    # Define vectors
    A = np.array([8, 9, 12])
    B = np.array([9, 6, 1])
    C = np.array([3, 5, 10])

    # Normalize to the plane
    n = np.cross(A, B)

    # Projection of C onto plane
    C_proj = C - np.dot(C, n) / np.dot(n, n) * n

    # Scalar projection magnitude
    scalar_proj = np.linalg.norm(C_proj)

    print("Projection vector of C onto plane:", C_proj)
    print("Scalar projection CAB =", scalar_proj)

    # starting problem 2 here with defining the variable t as equal to 3.
    t = 3
    #defining functions for the components
    x = np.sin(3*t)
    y = np.cos(t)
    z = np.sin(2*t)
    r = np.array([x, y, z])
    #calculating derivatives
    vx = 3*np.cos(3*t)
    vy = -np.sin(t)
    vz = 2*np.cos(2*t)
    v = np.array([vx, vy, vz])
    # Calculating Speed
    speed = np.linalg.norm(v)
    # calculating unit tangent
    ut = v / speed
    #Calculating Acceleration
    ax = -9*np.sin(3*t)
    ay = -np.cos(t)
    az = -4*np.sin(2*t)
    a = np.array([ax, ay, az])
    #calculating unit binomial and the normal vectors
    ub = np.cross(v, a)
    ub /= np.linalg.norm(ub)
    un = np.cross(ub, ut)
    #calculating the angles
    theta = np.degrees(np.arccos(v / speed))
    #doing the same with acceleration
    a_mag = np.linalg.norm(a)
    phi = np.degrees(np.arccos(a / a_mag))
    #calculating tangential and normal accelerations
    at = np.dot(a, ut)
    an = np.sqrt(a_mag**2 - at**2)
    #calculating the radius of curvature
    R = speed**2 / an
    #calculating the center of curvature
    center = r + R*un
    #printing question 2 results
    print("Velocity v:", v)
    print("Speed:", speed)
    print("Unit tangent ut:", ut)
    print("Angles (theta_x, theta_y, theta_z):", theta)
    print("Acceleration a:", a)
    print("Unit binormal ub:", ub)
    print("Unit normal un:", un)
    print("Angles (phi_x, phi_y, phi_z):", phi)
    print("Tangential accel at:", at)
    print("Normal accel an:", an)
    print("Radius of curvature:", R)
    print("Center of curvature:", center)
    #Moving on to question 3: starting by declaring my constants for graviation etc:
    G = 6.67430e-11  # grav const (im surprised e worked LOL)
    m1, m2 = 80, 50  #le mass
    r = 0.5
    #gravi equation
    F = G * m1 * m2 / r**2
    #printing the result
    print("Gravitational force between both people: ", float(F))
    #Moving on to question 4: starting with declaring my constants again this time for gravitational values for various bodies
    mars_g = 3.71/9.81 #note im doing the normalization for the weights inside the declarations directly since im lazy
    moon_g = 1.62/9.81
    jupiter_g = 24.79/9.81


    print("Weight on the moon with respect to W: " + str(moon_g)+"W") #I LOVE CASTING TO STRING URAAAH
    print("Weight on mars with respect to W: " + str(mars_g)+"W")
    print("Weight on Jupiter with respect to W: " + str(jupiter_g)+"W")



if __name__ == "__main__":
    main()

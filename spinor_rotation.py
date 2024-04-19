import sympy as sp
from qiskit.visualization import plot_bloch_vector
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Defining Identity matrix & Pauli matrices
I = sp.eye(2)
sigma_x = sp.Matrix([[0, 1], [1, 0]])
sigma_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
sigma_z = sp.Matrix([[1, 0], [0, -1]])

# Function to convert from sympy data type to native float type
def sympy_to_normal(sympy_matrix):
    return [1, float(sp.N(sympy_matrix[1])), float(sp.N(sympy_matrix[2]))]

# Function to rotate spinor by a given angle around a specified axis
def rotate_spinor(spinor, angle, axis):
    print("#### ROTATING SPINOR ONCE####\n")
    axis_unit_vector = [sp.sin(axis[1])*sp.cos(axis[2]),
                        sp.sin(axis[1])*sp.sin(axis[2]),
                        sp.cos(axis[1])
                        ]
    print("Axis unit vector: ", axis_unit_vector)

    spinor_matrix = sp.Matrix([[sp.cos(spinor[1]/2)],
                               [sp.sin(spinor[1]/2)*sp.exp(sp.I*spinor[2])]
                               ])
    print("Spinor matrix: ", spinor_matrix)
    # dot product of Pauli matrices and axis unit vector
    pauli_dot_axis = sigma_x*axis_unit_vector[0] + sigma_y*axis_unit_vector[1] + sigma_z*axis_unit_vector[2]

    rotation_matrix = sp.cos(angle/2)*I - sp.sin(angle/2)*pauli_dot_axis*sp.I
    print("Rotation matrix: ", rotation_matrix)

    # rotated spinor is the dot product of rotation matrix and spinor matrix
    rotated_spinor = rotation_matrix*spinor_matrix
    print("Rotated spinor: ", sp.N(rotated_spinor))

    # extract the rotated spinor's θ and φ
    if sp.re(rotated_spinor[0,0])== 0:
        if sp.im(rotated_spinor[0,0]) >= 0:
            A = sp.pi/2
        else: # if sine component is negative
            A = -sp.pi/2
        A = sp.pi/2
    else:
        A = sp.atan(sp.im(rotated_spinor[0,0])/sp.re(rotated_spinor[0,0]))

    if sp.re(rotated_spinor[1,0] )== 0:
        if sp.im(rotated_spinor[1,0]) >= 0:
            B = sp.pi/2
        else: # if sine component is negative
            B = -sp.pi/2
    else:
        B = sp.atan(sp.im(rotated_spinor[1,0])/sp.re(rotated_spinor[1,0]))

    spinor_theta = 2*sp.acos(sp.sqrt((sp.re(rotated_spinor[0,0])**2 + (sp.im(rotated_spinor[0,0])**2))))
    spinor_phi = B - A
    print("Rotated spinor theta: ", sp.N(spinor_theta))
    print("Rotated spinor phi: ", spinor_phi)

    return [1, spinor_theta, spinor_phi]

# Function to find single equivalent rotation from two rotations
def find_single_rotation(spinor, angle_1, axis_1, angle_2, axis_2):
    print("#### FINDING SINGLE EQUIVALENT ROTATION ####\n")
    # finding unit vector defining the two axes of rotation
    axis_1_unit_vector = sp.Matrix([sp.sin(axis_1[1]*sp.cos(axis_1[2])),
                                    sp.sin(axis_1[1]*sp.sin(axis_1[2])),
                                    sp.cos(axis_1[1])
                                    ])
    axis_2_unit_vector = sp.Matrix([sp.sin(axis_2[1]*sp.cos(axis_2[2])),
                                    sp.sin(axis_2[1]*sp.sin(axis_2[2])),
                                    sp.cos(axis_2[1])
                                    ])
    print("Axis 1 unit vector: ", axis_1_unit_vector)
    print("Axis 2 unit vector: ", axis_2_unit_vector)

    # finding the real component of the rotation matrix (cosine term)
    cosine_term = sp.cos(angle_1/2)*sp.cos(angle_2/2) - sp.sin(angle_1/2)*sp.sin(angle_2/2)*axis_1_unit_vector.dot(axis_2_unit_vector)
    print("Cosine term: ", cosine_term)

    # finding the imaginary component of the rotation matrix (sine term)
    sine_term = sp.sin(angle_1/2)*sp.cos(angle_2/2)*axis_1_unit_vector + sp.cos(angle_1/2)*sp.sin(angle_2/2)*axis_2_unit_vector - sp.sin(angle_1/2)*sp.sin(angle_2/2)*axis_1_unit_vector.cross(axis_2_unit_vector)
    print("Sine term: ", sine_term)

    single_angle_of_rotation = 2*sp.acos(cosine_term)
    print("Single angle of rotation: ", single_angle_of_rotation)

    if sine_term != 0:
        axis_of_rotation = sine_term/sp.sin(single_angle_of_rotation/2)
    else: # case when angle of rotation is 0, making sine term 0
        axis_of_rotation = [0,0,0]
    print("Axis of rotation: ", axis_of_rotation)
    # finding cartesian coordinates of the axis of rotation
    axis_theta = sp.acos(axis_of_rotation[2])

    if axis_of_rotation[0] == 0:
        if axis_of_rotation[2] != 0:
            axis_phi = 0
        else:
            axis_phi = sp.pi/2
    else:
        axis_phi = sp.atan(axis_of_rotation[1]/axis_of_rotation[0])

    print("Axis theta: ", axis_theta)
    print("Axis phi: ", axis_phi)

    spinor_after_one_rotation = rotate_spinor(spinor, single_angle_of_rotation, [1, axis_theta, axis_phi])
    print("Spinor after single rotation: ", spinor_after_one_rotation)

    return [1, axis_theta, axis_phi], spinor_after_one_rotation


# Function to plot Bloch spheres
def plot_bloch_spheres(sphere_list, title_list):
    fig = plt.figure(figsize = [15,12])

    positions = [
        [0,0.5],
        [0.333, 0.5],
        [0.667, 0.5],
        [0, 0],
        [0.333, 0]
    ]

    for m in range(len(sphere_list)):
        # Add the bloch sphere
        ax = fig.add_axes([positions[m][0], positions[m][1],0.333,0.5], axes_class = Axes3D)
        plot_bloch_vector(sphere_list[m], ax = ax, coord_type='spherical')
        ax.set_title(title_list[m], pad=25)

    plt.show()

# Input initial spinor and rotation axes
# spinor_theta = sp.rad(int(input("Enter the initial spinor's theta: ")))
# spinor_phi = sp.rad(int(input("Enter the initial spinor's phi: ")))
# spinor = [1, spinor_theta, spinor_phi]
#
# angle1 = sp.rad(int(input("Enter the 1st angle of rotation: ")))
# axis1_theta = sp.rad(int(input("Enter the 1st axis of rotation's theta: ")))
# axis1_phi = sp.rad(int(input("Enter the 1st axis of rotation's phi: ")))
# axis1 = [1, axis1_theta, axis1_phi]
#
# angle2 = sp.rad(int(input("Enter the 2nd angle of rotation: ")))
# axis2_theta = sp.rad(int(input("Enter the 2nd axis of rotation's theta: ")))
# axis2_phi = sp.rad(int(input("Enter the 2nd axis of rotation's phi: ")))
# axis2 = [1, axis2_theta, axis2_phi]

spinor = [1,sp.pi/4,0]
axis1 = [1,0,0]
axis2 = [1,0,0]
angle1 = sp.pi/4
angle2 = sp.pi/4

# Rotating spinors
rotated_spinor_1 = rotate_spinor(spinor, angle1, axis1)
rotated_spinor_2 = rotate_spinor(rotated_spinor_1, angle2, axis2)
# rotated_spinor_2 = rotate_spinor(spinor, sp.pi, axis2)

single_rot_axis, single_rot_spinor = find_single_rotation(spinor, sp.pi/4, axis1, sp.pi/4, axis2)

# Putting all spinors in a list
spheres = [sympy_to_normal(spinor),
           sympy_to_normal(rotated_spinor_1),
           sympy_to_normal(rotated_spinor_2),
           sympy_to_normal(spinor),
           sympy_to_normal(single_rot_spinor)
           ]

print(spheres)

# Defining titles for the plots
titles = ["Initial Spinor",
          "Spinor after 1st Rotation",
          "Spinor after 2nd Rotation",
          "Initial Spinor",
          "After Single Rotation"
          ]

# Plotting Bloch spheres
plot_bloch_spheres(spheres, titles)

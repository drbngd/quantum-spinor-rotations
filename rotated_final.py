import sympy as sp
import numpy as np
from scipy.linalg import expm
from qiskit.visualization import plot_bloch_vector
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define constants and initial matrices
I = sp.eye(2)
q_0 = sp.Matrix([[1], [0]])
q_1 = sp.Matrix([[0], [1]])
sigma_x = sp.Matrix([[0, 1], [1, 0]])
sigma_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
sigma_z = sp.Matrix([[1, 0], [0, -1]])

# Setup initial spinor and rotation axes
spinor = [1, 0, 0]
axis1 = [1, sp.pi / 2, 0]
axis2 = [1, sp.pi / 2, 0]

def rotate_spinor(spinor, angle, axis):
    """Rotate spinor by a given angle around a specified axis."""
    n_vector = [sp.sin(axis[1]) * sp.cos(axis[2]), sp.sin(axis[1]) * sp.sin(axis[2]), sp.cos(axis[1])]
    spinor_matrix = sp.Matrix([[sp.cos(spinor[1] / 2)], [sp.sin(spinor[1] / 2) * sp.exp(sp.I * spinor[2])]])

    sigma_dot_axis = sigma_x * n_vector[0] + sigma_y * n_vector[1] + sigma_z * n_vector[2]
    rot_matrix = sp.cos(angle / 2) * I - sp.I * sp.sin(angle / 2) * sigma_dot_axis
    rotated_spinor = rot_matrix * spinor_matrix

    # Calculate new angles theta and phi from the rotated spinor
    theta = 2 * sp.acos(sp.sqrt(sp.re(rotated_spinor[0, 0]) ** 2 + sp.im(rotated_spinor[0, 0]) ** 2))
    phi = sp.atan2(sp.im(rotated_spinor[1, 0]), sp.re(rotated_spinor[1, 0])) - sp.atan2(sp.im(rotated_spinor[0, 0]), sp.re(rotated_spinor[0, 0]))

    return [1, theta, phi]

# Rotate spinors
spinor_r1 = rotate_spinor(spinor, sp.pi / 2, axis1)
spinor_r2 = rotate_spinor(spinor_r1, sp.pi / 2, axis2)

def sympy_to_numpy(sympy_matrix):
    """Convert a sympy matrix to a numpy array."""
    return [1, float(sp.N(sympy_matrix[1])), float(sp.N(sympy_matrix[2]))]

def find_single_rotation(spinor, angle_1, axis_1, angle_2, axis_2):
    """Find a single equivalent rotation from two rotations."""
    n1_vect = sp.Matrix([sp.sin(axis_1[1]) * sp.cos(axis_1[2]), sp.sin(axis_1[1]) * sp.sin(axis_1[2]), sp.cos(axis_1[1])])
    n2_vect = sp.Matrix([sp.sin(axis_2[1]) * sp.cos(axis_2[2]), sp.sin(axis_2[1]) * sp.sin(axis_2[2]), sp.cos(axis_2[1])])

    cosine_term = sp.cos(angle_1 / 2) * sp.cos(angle_2 / 2) - sp.sin(angle_1 / 2) * sp.sin(angle_2 / 2) * n1_vect.dot(n2_vect)
    single_angle_of_rotation = 2 * sp.acos(cosine_term)

    sine_term = sp.sin(angle_1 / 2) * sp.cos(angle_2 / 2) * n1_vect + sp.cos(angle_1 / 2) * sp.sin(angle_2 / 2) * n2_vect - sp.sin(angle_1 / 2) * sp.sin(angle_2 / 2) * n1_vect.cross(n2_vect)
    axis_of_rotation = sine_term / sp.sin(single_angle_of_rotation / 2) if sine_term != 0 else [0, 0, 0]

    return [1, sp.acos(axis_of_rotation[2]), sp.atan(axis_of_rotation[1] / axis_of_rotation[0])], rotate_spinor(spinor, single_angle_of_rotation, [1, sp.acos(axis_of_rotation[2]), sp.atan(axis_of_rotation[1] / axis_of_rotation[0])])

# Perform and display the results of single rotations
single_rot_axis, single_rot_spinor = find_single_rotation(spinor, sp.pi / 2, axis1, sp.pi / 2, axis2)
states = [sympy_to_numpy(spinor), sympy_to_numpy(spinor_r1), sympy_to_numpy(spinor_r2), sympy_to_numpy(spinor), sympy_to_numpy(single_rot_spinor)]

# Plotting results
fig = plt.figure(figsize=[12, 8])
positions = [[0, 0.5], [0.333, 0.5], [0.667, 0.5], [0, 0], [0.333, 0]]
for m, state in enumerate(states):
    ax = fig.add_axes([positions[m][0], positions[m][1], 0.333, 0.5], projection='3d')
    plot_bloch_vector(state, ax=ax, coord_type='spherical')
plt.show()

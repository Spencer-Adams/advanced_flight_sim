import numpy as np # type: ignore
import os # type: ignore
import math # type: ignore
import sympy as sp # type: ignore
from sympy import symbols, cos, sin, I, re, im, sqrt # type: ignore
import time # type: ignore
import threading # type: ignore
import matplotlib.pyplot as plt # type: ignore
import matplotlib.image as mpimg # type: ignore
import matplotlib.animation as animation  # type: ignore
from matplotlib.animation import PillowWriter  # type: ignore
from matplotlib.offsetbox import OffsetImage, AnnotationBbox # type: ignore
from io import BytesIO
from tqdm import tqdm # type: ignore
# Helper functions below

def euler_to_quat(euler):
    # Eq. 1.6.2 in the book. 
    phi = euler[0]
    theta = euler[1]
    psi = euler[2]
    hphi = 0.5*phi
    htheta = 0.5*theta
    hpsi = 0.5*psi
    cos_phi = np.cos(hphi)
    cos_theta = np.cos(htheta)
    cos_psi = np.cos(hpsi)
    sin_phi = np.sin(hphi)
    sin_theta = np.sin(htheta)
    sin_psi = np.sin(hpsi)
    cos_phi_cos_theta = cos_phi*cos_theta
    sin_phi_sin_theta = sin_phi*sin_theta
    sin_phi_cos_theta = sin_phi*cos_theta
    cos_phi_sin_theta = cos_phi*sin_theta
    quat = np.array([0.0,0.0,0.0,0.0])
    quat[0] = cos_phi_cos_theta*cos_psi + sin_phi_sin_theta*sin_psi
    quat[1] = sin_phi_cos_theta*cos_psi - cos_phi_sin_theta*sin_psi
    quat[2] = cos_phi_sin_theta*cos_psi + sin_phi_cos_theta*sin_psi
    quat[3] = cos_phi_cos_theta*sin_psi - sin_phi_sin_theta*cos_psi
    return quat

def quat_to_euler(quat): 
    e0 = quat[0]
    ex = quat[1]
    ey = quat[2]
    ez = quat[3]
    euler = np.array([0.0,0.0,0.0])
    cos_pi_4 = ex/np.cos(np.pi*0.25)
    arcsin_quat_cos_pi_4 = np.arcsin(cos_pi_4)
    if (e0*ey-ex*ez == 0.5):
        euler[0] = 2*arcsin_quat_cos_pi_4
        euler[1] = 0.5*np.pi
        euler[2] = 0.0
    elif (e0*ey-ex*ez == 0.5):
        euler[0] = 2*arcsin_quat_cos_pi_4
        euler[1] = -0.5*np.pi
        euler[2] = 0.0
    else: 
        euler[0] = np.arctan2(2*(e0*ex + ey*ez), (e0**2 + ez**2 -ex**2 -ey**2))
        euler[1] = np.arcsin(2*(e0*ey-ex*ez))
        euler[2] = np.arctan2(2*(e0*ez + ex*ey), (e0**2 + ex**2 -ey**2 -ez**2))
    return euler 

def quat_base_to_dependent(vec, quat):
    # Eq. (1.5.4 in the book which uses the quaternion product twice to get Eq. 1.5.7) Can be used to go from earth-fixed to body-fixed coordinates.
    vx = vec[0]
    vy = vec[1]
    vz = vec[2]
    e0 = quat[0]
    ex = quat[1]
    ey = quat[2]
    ez = quat[3]
    t0 = -vx*ex - vy* ey - vz*ez 
    tx =  vx*e0 + vy* ez - vz* ey 
    ty =  -vx* ez + vy* e0 + vz* ex 
    tz = vx* ey - vy* ex + vz* e0
    v2 = np.array([0.0,0.0,0.0])
    v2[0] = e0*tx - ex*t0 - ey*tz + ez*ty 
    v2[1] = e0*ty + ex*tz - ey*t0 - ez*tx 
    v2[2] = e0*tz - ex*ty + ey*tx - ez*t0
    return v2 

def quat_dependent_to_base(v2, quat): # can be used to go from body-fixed to earth-fixed coordinates
    vx = v2[0]
    vy = v2[1]
    vz = v2[2]
    e0 = quat[0]
    ex = quat[1]
    ey = quat[2]
    ez = quat[3]
    t0 = vx*ex + vy*ey + vz*ez 
    tx =  vx*e0 - vy*ez + vz*ey 
    ty =  vx*ez + vy*e0 - vz*ex 
    tz = -vx*ey + vy*ex + vz*e0
    v1 = np.array([0.0,0.0,0.0])
    v1[0] = e0*tx + ex*t0 + ey*tz - ez*ty 
    v1[1] = e0*ty - ex*tz + ey*t0 + ez*tx 
    v1[2] = e0*tz + ex*ty - ey*tx + ez*t0
    return v1

def parse_dictionary_or_return_default(dictionary, keys, default):
    """Safely get nested dictionary values. It makes it so that if a json key is not found, it returns a default value instead of throwing an error."""
    for key in keys:
        dictionary = dictionary.get(key, {})
        if not isinstance(dictionary, dict) and key != keys[-1]:
            return default
    if dictionary == {} and default is not None:
        print("The dictionary ", keys, " is empty. Returning default value of ", default, ".")
        return default
    elif dictionary != {}:
        return dictionary 
    else:
        raise ValueError(f"Key {keys} not found in dictionary and no default value provided.")
    
def plot_xy_array(xy_array: np.array, color: str = 'black', linewidth: float = 0.5):
    plt.plot(xy_array[:, 0], xy_array[:, 1], color=color, linewidth=linewidth)

def compute_arc_length(points: np.ndarray) -> float:
    """
    Compute the arc length of a curve defined by a sequence of points.
    
    Args:
        points: np.ndarray of shape (N, D), where D is typically 2 or 3

    Returns:
        Total arc length as a float.
    """
    # Compute differences between consecutive points
    deltas = np.diff(points, axis=0)
    # Compute Euclidean distances between consecutive points
    segment_lengths = np.linalg.norm(deltas, axis=1)
    # Sum to get total arc length
    return np.sum(segment_lengths)

def vector_magnitude(vector: np.array): # Used to help A6 in the project. The streamlines need to be integrated forward using the unit velocity vector
    """
    Calculate the magnitude of a vector.

    Parameters:
    vector (list): A list representing the vector.

    Returns:
    float: The magnitude of the vector.
    """
    magnitude = math.sqrt(sum([element**2 for element in vector]))
    return magnitude

def unit_vector(vector: np.array): # Used to help A6 in the project. The streamlines need to be integrated forward using the unit velocity vector
    """Calculate the unit vector of a given vector."""
    # make sure that vector is a list of two or three float elements
    if len(vector) < 2 or not all(isinstance(coord, float) for coord in vector):
        raise ValueError("The vector must be a list of two or three float elements.")
    # calculate the magnitude of the vector
    magnitude = vector_magnitude(vector)
    # calculate the unit vector
    unit_vector = [element/magnitude for element in vector]
    return unit_vector

def calc_tangent(main: np.ndarray, prev: np.ndarray, next_: np.ndarray) -> np.ndarray:
    """
    Compute an approximate tangent vector at `main` based on the adjacent points.
    
    Args:
        main: np.ndarray of shape (2,), the main point [x, y]
        prev: np.ndarray of shape (2,), the previous point [x, y]
        next_: np.ndarray of shape (2,), the next point [x, y]

    Returns:
        Tangent vector (not yet normalized)
    """
    # Tangent is approximated by the vector from prev to next
    tangent = next_ - prev
    return tangent

def calc_normal(tangent: np.ndarray) -> np.ndarray:
    """
    Compute a normal vector (90 degrees CCW) from a 2D tangent vector.

    Args:
        tangent: np.ndarray of shape (2,), the tangent vector [dx, dy]

    Returns:
        Normal vector (not yet normalized)
    """
    dx, dy = tangent
    normal = np.array([-dy, dx])  # Rotate CCW
    return normal

def rk4(start: np.array, direction: int, step_size: float,  move_off_direction_func: callable, function: callable):
    """This function performs a Runge-Kutta 4th order integration

    Args:
        - start (list): A list of two float values representing x and y coordinates.
        - direction (int): An integer value representing the direction of integration.
        - step_size (float): The step size for the integration.
        - function (callable): A function that calculates the derivative at a given point.

    Returns:
        list: A numpy array of float values representing the new coordinates after integration.
    """
    # make sure that start is a list of two or three float elements
    if len(start) < 2 or not all(isinstance(coord, float) for coord in start):
        raise ValueError("The start must be a list of two float values representing x and y coordinates.")
    # make sure that direction is an integer
    if not isinstance(direction, int):
        raise TypeError("The direction must be an integer. Please provide an integer value.")
    # set the step size
    # print(direction)
    point = start
    h = direction*step_size
    # if all of the function values in the np.array are really close to zero, step the point off in the move off direction. 
    while np.all(np.abs(function(point)) < 1e-12):
        for i in range(len(point)):
            print("The function values are too small. Stepping the point off in the move off direction.")
            print("point: ", point)
            print("velocity: ", function(point))
            move_off_direction = move_off_direction_func(point[0])[i]
            point[i] = point[i] + move_off_direction[i]*1e-6
            break
    # set the initial values of k1, k2, k3, and k4
    k1 = np.array(function(point))
    k2 = np.array(function(point + 0.5 * h * k1))
    k3 = np.array(function(point + 0.5 * h * k2))
    k4 = np.array(function(point + h * k3))
    point_new = point + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    point_new = np.array(point_new)
    return point_new

def bisection(f, a, b, tol=1e-8, max_iter=100):
    """
    Perform the bisection method to find a root of the function `f` in the interval [a, b].

    Args:
        f: A callable that accepts a float and returns a float.
        a: Lower bound of the interval.
        b: Upper bound of the interval.
        tol: Desired tolerance for the root (default 1e-8).
        max_iter: Maximum number of iterations to prevent infinite loops.

    Returns:
        Approximate root of f in [a, b].

    Raises:
        ValueError: If f(a) and f(b) do not have opposite signs.
    """
    fa = f(a)
    fb = f(b)
    
    if fa * fb > 0:
        raise ValueError("Function must have opposite signs at the endpoints a and b.")

    for i in range(max_iter):
        c = 0.5 * (a + b)
        fc = f(c)

        if abs(fc) < tol or (b - a) < tol:
            return c

        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc

    raise RuntimeError("Bisection method did not converge within the maximum number of iterations.")

def sort_points_counterclockwise(points: np.ndarray) -> np.ndarray:
    """
    Sort an array of 2D points in counterclockwise order around their centroid.

    Parameters:
    - points (np.ndarray): shape (n, 2), each row is an (x, y) point.

    Returns:
    - sorted_points (np.ndarray): shape (n, 2), points in CCW order.
    """
    points = np.asarray(points)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("Input must be an (n, 2) array of xy points.")
    # Step 1: Compute centroid
    centroid = np.mean(points, axis=0)
    # Step 2: Compute angle to each point
    angles = np.arctan2(points[:,1] - centroid[1], points[:,0] - centroid[0])
    # Step 3: Sort points by angle
    sorted_indices = np.argsort(angles)
    sorted_points = points[sorted_indices]
    return sorted_points


def list_to_range(three_element_list: list):
    """First element is the start, second element is the end, and the third element is the step size."""
    if len(three_element_list) != 3:
        raise ValueError("The list must contain three elements. The first element is the start, the second element is the end, and the third element is the step size.")
    
    start = three_element_list[0]
    end = three_element_list[1]
    step_size = three_element_list[2]
    
    if start != end:
        # Calculate the number of points to include in the range
        num_points = int(round((end - start) / step_size)) + 1
        values = np.linspace(start, end, num_points)
    else:
        values = np.array([start])
    
    return values

def xy_to_r_theta(x: float, y: float):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta

def r_theta_to_xy(r: float, theta: float):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def central_difference(f_plus: float, f_minus: float, step_size: float):
    derivative = (f_plus - f_minus)/(2*step_size)
    return derivative

def polar_vector(theta, cartesian_vector):
    """This function converts the cartesian velocity to polar velocity (can go from z, zeta, or chi plane to polar velocity)"""
    r = cartesian_vector[0]*np.cos(theta) + cartesian_vector[1]*np.sin(theta)
    # print("\nradial velocity", velocity_r)
    theta = cartesian_vector[1]*np.cos(theta) - cartesian_vector[0]*np.sin(theta)
    # print("theta velocity", velocity_theta)
    polar_vector = np.array([r, theta])
    return polar_vector

# Function to remove the middle element
def remove_middle_element(arr):
    mid_idx = len(arr) // 2
    return np.delete(arr, mid_idx, axis=0)

def numerical_derivative(func, x, r_values=np.array([0.0]), theta_values=np.array([0.0]), h=1e-6, D=0.0):
    """ Compute first derivative using central difference. """
    f_plus = func(x + h, r_values, theta_values, D)
    f_minus = func(x - h, r_values, theta_values, D)
    # print(f"f_plus: {f_plus}, f_minus: {f_minus}")  # Debug print
    return (f_plus - f_minus) / (2 * h)

def numerical_second_derivative(func, x, r_values=np.array([0.0]), theta_values=np.array([0.0]), h=1e-6, D=0.0):
    """ Compute second derivative using central difference. """
    f_plus = func(x + h, r_values, theta_values, D)
    f = func(x, r_values, theta_values, D)
    f_minus = func(x - h, r_values, theta_values, D)
    # print(f"f_plus: {f_plus}, f: {f}, f_minus: {f_minus}")  # Debug print
    return (f_plus - 2 * f + f_minus) / (h**2)

def newtons_method(func, x0, r_values=np.array([0.0]), theta_values=np.array([0.0]), tol=1e-10, max_iter=1000, D=0.0):
    """ Newton's method to find extrema of a function. """
    x = x0
    
    for i in range(max_iter):
        f_prime = numerical_derivative(func, x, r_values, theta_values)
        f_double_prime = numerical_second_derivative(func, x, r_values, theta_values)
        
        if abs(f_double_prime) < 1e-12:
            raise ValueError("Second derivative is too small — possible inflection point or flat region.")
        
        x_new = x - f_prime / f_double_prime
        epsilon = abs(x_new - x)
        # stop if epsilon is less than tolerance or if the previous x is the same as the new x out to 14 decimal places
        if epsilon < tol:# or np.isclose(x, x_new, atol=1e-16):
            if epsilon < tol:
                print(f"Converged in {i+1} iterations.")
                # print epsilon at convergence
                print(f"epsilon = {epsilon:.16f}")
            elif np.isclose(x, x_new, atol=1e-16):
                print(f"Converged in {i+1} iterations because the previous iteration was the same as this one to 16 digits.")
                print(f"epsilon = {epsilon:.16f}")
            return x_new, func(x_new, r_values, theta_values, D)  # Return the extremum and its value
        
        # print iteration, and epsilon compared to tol every 10 steps
        if i % 5 == 0:
            print(f"Iteration {i+1}: x = {x_new:.16f}, epsilon = {epsilon:.16f}")
            # print(f"Iteration {i+1}: x = {x_new:.6f}")
        
        x = x_new

    raise ValueError("Newton's method did not converge.")

def polyfit(func, r_values, theta_values, gamma_vals, order_of_polynomial, is_plot, xlabel, ylabel, plot_title, D, is_get_poly_coeffs_and_export_function_to_10000_points = False):
    # time_1 = time.time()
    # print("length of gamma_vals", len(gamma_vals))
    
    # Evaluate function at all gamma values
    appellian_vals = np.array([[gamma, func(gamma, r_values, theta_values)] for gamma in gamma_vals])

    # Fit a polynomial of specified order
    coeffs = np.polyfit(appellian_vals[:, 0], appellian_vals[:, 1], order_of_polynomial)

    # Warn if leading coefficient is negative
    if coeffs[0] < 0:
        print("Warning: The highest order coefficient is negative. This may indicate that the polynomial is not a good fit for the data.")
        print("This warning is at D =", D)

    # Derivative and roots (extrema)
    derivative_coeffs = np.polyder(coeffs)
    extrema = np.roots(derivative_coeffs)
    # print(f"Extrema found: {extrema}")
    # Separate real and complex extrema (with tolerance)
    real_extrema = [x.real for x in extrema if np.isclose(x.imag, 0)]
    complex_extrema = [x for x in extrema if not np.isclose(x.imag, 0)]

    # Determine which extrema to evaluate
    if len(real_extrema) == 1 and len(complex_extrema) == 2:
        # print("Case: One real root and a complex conjugate pair — selecting the real root as the minimum.")
        selected_extrema = [real_extrema[0]]
    elif len(real_extrema) == 3 and len(set(np.round(real_extrema, 12))) == 1:
        # print("Case: Triple real root — all roots identical. Choosing one.")
        selected_extrema = [real_extrema[0]]
    elif len(real_extrema) >= 1:
        # print("Case: All real distinct roots — evaluating all to find the minimum.")
        selected_extrema = real_extrema
    else:
        print("\n\n\n")
        print("coeffs:", coeffs)
        print("derivative_coeffs:", derivative_coeffs)
        print("extrema!", extrema)
        print("\n\n\n")
        raise ValueError("Unexpected root configuration in polynomial derivative.")
    # print(f"Selected extrema: {selected_extrema}")
    # Evaluate the true function and polynomial at extrema
    # print("length of selected extrema", len(selected_extrema))
    func_vals = np.array([func(x, r_values, theta_values) for x in selected_extrema])
    poly_vals = np.array([np.polyval(coeffs, x) for x in selected_extrema])

    # for i, (f_val, p_val, x_val) in enumerate(zip(func_vals, poly_vals, selected_extrema)):
    #     if abs(f_val - p_val) > 1e-12:
    #         print("\nWarning: The function and polynomial values differ at extremum.")
    #         print(f"This warning is at D = {D}")
    #         print(f"extremum = {x_val}")
    #         print(f"func_val = {f_val}, poly_val = {p_val}")
    #         print(f"absolute difference = {abs(f_val - p_val)}\n")

    # Choose minimum among selected extrema
    min_index = np.argmin(func_vals)
    extremum_min = selected_extrema[min_index]
    value_at_min = func_vals[min_index]

    # if is_get_poly_coeffs_and_export_function_to_10000_points:
        # gamma_fine, poly_vals_fine = get_poly_coeffs_and_export_function_to_10000_points(coeffs, gamma_vals, D)

    # Optional plot
    # if is_plot:
    #     plt.plot(appellian_vals[:, 0], appellian_vals[:, 1], marker='o', linestyle='-', color='black', markersize=3)
    #     plt.xlabel(xlabel)
    #     plt.ylabel(ylabel)
    #     plt.title(plot_title)

    #     # Add extrema to plot
    #     y_range = plt.ylim()[1] - plt.ylim()[0]
    #     y_offset = 0.05 * y_range
    #     x_text = 0.5 * (plt.xlim()[0] + plt.xlim()[1])
    #     y_text = plt.ylim()[1] - y_offset

    #     for i, x_val in enumerate(selected_extrema):
    #         plt.text(
    #             x_text,
    #             y_text - i * y_offset,
    #             f"({x_val:.8f}, {func_vals[i]:.2f})",
    #             fontsize=8,
    #             ha='center',
    #             va='top'
    #         )
    #     plt.show()

    # time_2 = time.time()
    # print("Time taken to find extremum using polyfit:", time_2 - time_1)

    return extremum_min, value_at_min

def get_poly_coeffs_and_export_function_to_10000_points(poly_coeffs, gamma_vals, D):
    """
    This function takes the polynomial coefficients (assumed 4th order: [B1, B2, B3, B4, B5])
    and plots the polynomial function over 10,000 points in the range of gamma_vals.
    Calculates values as: val = B1*x^4 + B2*x^3 + B3*x^2 + B4*x + B5
    """
    # Create 10,000 points in the range of gamma_vals
    gamma_max = 25#np.max(gamma_vals) 
    gamma_min = -20#np.min(gamma_vals)
    gamma_fine = np.linspace(gamma_min, gamma_max, 10000)
    # Unpack coefficients (assume length 5)
    B1, B2, B3, B4, B5 = poly_coeffs
    # Evaluate the polynomial manually
    poly_vals = B1 * gamma_fine**4 + B2 * gamma_fine**3 + B3 * gamma_fine**2 + B4 * gamma_fine + B5
    gamma_fine /= (4*np.pi*10.0)
    gamma_fine /= 0.08715574
    gamma_vals /= (4*np.pi*10.0) 
    gamma_vals /= 0.08715574
    poly_vals /= 10000.0
    # Make an array of shape (10000, 2) where the first column is gamma_fine and the second column is poly_vals
    poly_array = np.column_stack((gamma_fine, poly_vals))
    # Save the array to a text file
    np.savetxt(f"polynomial_fit_10000_points_D{D}_gammahat_divided_by_gamma_kutta_vals_{gamma_vals}.txt", poly_array, header="gamma_hat, s_hat", delimiter=",")
    # Plot the polynomial
    plt.plot(gamma_fine, poly_vals, color='blue', linewidth=1)
    plt.xlabel("gamma")
    plt.ylabel("Polynomial Value")
    plt.title(f"Polynomial Fit over 10,000 Points at D = {D}")
    plt.show()
    return gamma_fine, poly_vals

def combine_two_plots(plot1_filename, plot2_filename, output_filename):
    """
    Combine two plots into one figure with two subplots stacked vertically.
    
    Parameters:
    - plot1_filename: str, path to the first plot image file.
    - plot2_filename: str, path to the second plot image file.
    - output_filename: str, path to save the combined figure.
    """
    # Create a figure with two subplots
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 12))
    
    # Load and display the first plot
    img1 = mpimg.imread(plot1_filename)
    axes[0].imshow(img1)
    axes[0].axis('off')  # Turn off axes for the image
    # axes[0].set_title("Plot 1", fontsize=16)

    # Load and display the second plot
    img2 = mpimg.imread(plot2_filename)
    axes[1].imshow(img2)
    axes[1].axis('off')  # Turn off axes for the image
    # axes[1].set_title("Plot 2", fontsize=16)

    # Adjust layout and save the combined figure
    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close(fig)

def combine_plots(plot_names_array, output_filename):
    num_plots = len(plot_names_array)
    num_rows = math.ceil(num_plots / 2)

    fig, axes = plt.subplots(
        nrows=num_rows,
        ncols=2,
        figsize=(6.5, 7), # size of the figure in inches 
        gridspec_kw={'wspace': 0.05, 'hspace': -0.45} # wspace adjusts the space between columns, hspace adjusts the space between rows
    )

    axes = axes.flatten()  # flatten to make indexing simpler

    for i, plot_name in enumerate(plot_names_array):
        img = mpimg.imread(plot_name)
        axes[i].imshow(img)
        axes[i].axis('off')

    # Turn off any unused axes if number of plots is odd
    for j in range(len(plot_names_array), len(axes)):
        axes[j].axis('off')

    # plt.tight_layout(pad=0.5, w_pad=0.1, h_pad=0.1)
    plt.savefig(output_filename, bbox_inches='tight')
    plt.close(fig)

def create_animation_from_all_figs_in_folder(folder, output_filename):
    """
    Create an animation from all .png figures in a specified folder.

    Parameters:
    - folder: str, path to the folder containing the .png figures.
    - output_filename: str, path to save the animation.
    """

    # Create a figure and axis for the animation
    fig, ax = plt.subplots(figsize=(10, 8))

    # List to store sorted image filenames
    image_filenames = []

    # Find all .png files in the folder that end with a number
    for filename in os.listdir(folder):
        if filename.endswith(".png") and filename[-5].isdigit():
            image_filenames.append(os.path.join(folder, filename))

    # Sort the filenames based on the numeric suffix
    image_filenames.sort(key=lambda x: int(x.split("__")[-1].split(".")[0]))

    # Load .png images into memory for animation
    png_images = [plt.imread(img_file) for img_file in image_filenames]

    # Function to update the frame in the animation
    def update_frame(frame_idx):
        ax.clear()
        ax.imshow(png_images[frame_idx])
        ax.axis("off")  # Turn off axes for clean animation

    # Set interval to 125ms for 1/8 second per frame
    interval = 125
    fps = 8  # 1000 / 125 ms = 8 fps

    # Create an animation from the images
    ani = animation.FuncAnimation(fig, update_frame, frames=len(png_images), interval=interval)

    # Save the animation as a GIF
    ani.save(output_filename, writer="pillow", fps=fps)

    print("Animation created successfully!")

    # Clean up the images
    for filename in image_filenames:
        os.remove(filename)

def riemann_integration(r_values, theta_values, function: callable):
    r0, r1 = r_values
    theta0, theta1 = theta_values
    dtheta = theta1 - theta0
    if dtheta <= 0:
        dtheta += 2 * np.pi
    if r0 == r1:  # 1D
        f0 = function([r0, theta0])
        return f0 * r0 * dtheta
    else: 
        dr = r1 - r0
        f = function([r0, theta0])
        return f * r0 * dr * dtheta
    
def trapezoidal_integration(r_values, theta_values, function: callable):
    r0, r1 = r_values
    theta0, theta1 = theta_values
    dtheta = theta1 - theta0
    if dtheta <= 0:
        dtheta += 2 * np.pi
    if r0 == r1:  # 1D
        f0 = function([r0, theta0])
        f1 = function([r0, theta1])
        return 0.5 * (f0 + f1) * r0 * dtheta
    else:  # 2D with midpoints
        dr = r1 - r0
        r_mid = 0.5 * (r0 + r1)
        theta_mid = 0.5 * (theta0 + theta1)
        # Corner evaluations
        f00 = function([r0, theta0]) * r0
        f01 = function([r0, theta1]) * r0
        f10 = function([r1, theta0]) * r1
        f11 = function([r1, theta1]) * r1
        fmm = function([r_mid, theta_mid]) * r_mid
        # Average with midpoint
        avg = (f00 + f01 + f10 + f11 + 4 * fmm) / 8
        return avg * dr * dtheta
        
def simpson_13_integration(r_values, theta_values, function: callable):
    """Simpson's 1/3 rule for 1D or single 2D cell in polar coords."""
    r0, r1 = r_values
    theta0, theta1 = theta_values
    # Ensure positive theta difference
    dtheta = theta1 - theta0
    if dtheta <= 0:
        dtheta += 2 * np.pi
    if r0 == r1:
        # 1D case: integrate over theta at fixed r
        theta_mid = 0.5 * (theta0 + theta1)
        f0 = function([r0, theta0])
        f1 = function([r0, theta_mid])
        f2 = function([r0, theta1])
        return (f0 + 4*f1 + f2) * r0 * dtheta / 6
    else:
        # 2D case
        dr = r1 - r0
        r_mid = 0.5 * (r0 + r1)
        theta_mid = 0.5 * (theta0 + theta1)
        f00 = function([r0, theta0]) * r0
        f01 = function([r0, theta1]) * r0
        f10 = function([r1, theta0]) * r1
        f11 = function([r1, theta1]) * r1
        fmm = function([r_mid, theta_mid]) * r_mid
        f0m = function([r0, theta_mid]) * r0
        f1m = function([r1, theta_mid]) * r1
        fm0 = function([r_mid, theta0]) * r_mid
        fm1 = function([r_mid, theta1]) * r_mid
        return (f00 + f01 + f10 + f11 + 4*fmm + 2*(f0m + f1m + fm0 + fm1)) * dr * dtheta / 36

def example_polar_function(rtuple):
    """Example polar function for testing."""
    # Example: r * cos(theta) + D
    r, theta = rtuple
    return r * np.cos(theta)

def analytic_integration_in_theta_of_example_polar_function(r_values, theta_values):
    """Analytic integration of the example polar function."""
    r0 = r_values[0]
    theta0, theta1 = theta_values[0], theta_values[-1]
    if r0 == 0:
        raise ValueError("r_values must not contain zero to avoid division by zero in the integration.")
    # Integral of r^2 * cos(theta) dtheta from theta0 to theta1
    integral = r0**2*(np.sin(theta1) - np.sin(theta0)) 
    return integral 

def analytic_integration_in_r_and_theta_example_function(r_values, theta_values):
    r0, r1 = r_values[0], r_values[-1]
    theta0, theta1 = theta_values[0], theta_values[-1]
    return ((r1**3/3)-(r0**3/3))*(np.sin(theta1)-np.sin(theta0))

def single_romberg_integration(integral_low, integral_high):
    """Takes the two fidelities of integration to make a better estimate"""
    return 4/3 * integral_high - 1/3 * integral_low

if __name__ == "__main__":
    # combine_two_plots("figures/zeta_-0.25_D_sweep_alpha_5.png", "figures/zeta_-0.25_D_0_alpha_5.png", "combined_plot.svg")
    theta0 = np.pi/4 
    theta1 = 3*np.pi/8
    r0 = 1.6
    r1 = 1.8
    # split up the interval between theta0 and theta1 into N intervals, and add the trapezoidal integration of the example polar function over that interval
    N = 111
    theta_values = np.linspace(theta0, theta1, N)
    r_values = np.array([r0, r1])
    total_riemann_integral = 0.0
    total_trapezoidal_integral = 0.0
    total_simpson_integral = 0.0
    total_romberg_integral = 0.0
    if r0 == r1:
        print("Calculating trapezoidal integrals for each subinterval...")
        analytic_integral = analytic_integration_in_theta_of_example_polar_function(r_values=np.array([r0, r1]),theta_values=np.array([theta0, theta1]))
        for i in tqdm(range(len(theta_values) - 1), desc = "Calculating thetas"):
            theta0 = theta_values[i]
            theta1 = theta_values[i + 1]
            thetamid = (theta0+theta1)/2
            refined_trapezoidal = 0.0
            thetas_unrefined = np.array([theta0, theta1])
            thetas_refined = np.array([theta0,thetamid,theta1])
            riemann_integral = riemann_integration(r_values, theta_values=thetas_unrefined, function = example_polar_function)
            trapezoidal_integral = trapezoidal_integration(r_values=r_values, theta_values=thetas_unrefined, function=example_polar_function)
            for j in range(len(thetas_refined)-1):
                refined_trapezoidal += trapezoidal_integration(r_values=r_values,theta_values = np.array([thetas_refined[j], thetas_refined[j+1]]), function=example_polar_function)
            romberg_integral = single_romberg_integration(trapezoidal_integral, refined_trapezoidal)
            simpson_integral = simpson_13_integration(r_values = r_values, theta_values=np.array([theta0, theta1]), function = example_polar_function)
            total_riemann_integral += riemann_integral
            total_trapezoidal_integral += trapezoidal_integral
            total_simpson_integral += simpson_integral 
            total_romberg_integral += romberg_integral
        percent_error_riemann = abs(total_riemann_integral-analytic_integral)/(analytic_integral)*100
        percent_error_trapezoidal = abs(total_trapezoidal_integral-analytic_integral)/(analytic_integral)*100
        percent_error_simpson13 = abs(total_simpson_integral-analytic_integral)/(analytic_integral)*100
        percent_error_romberg = abs(total_romberg_integral-analytic_integral)/(analytic_integral)*100
        print("Analytic integration result:   ", analytic_integral)
        print("Riemann integration result:    ", total_riemann_integral)
        print("Riemann Percent Difference:    ", percent_error_riemann)
        print("Trapezoidal integration result:", total_trapezoidal_integral)
        print("Trapezoidal Percent Difference:", percent_error_trapezoidal)
        print("Simpson 1/3 integration result:", total_simpson_integral)
        print("Simpson 1/3 Percent Difference:", percent_error_simpson13)
        print("Romberg Integration Result:    ", total_romberg_integral)
        print("Romberg Percent Difference:    ", percent_error_romberg)
    else:
        analytic_integral = analytic_integration_in_r_and_theta_example_function(r_values=np.array([r0, r1]),theta_values=np.array([theta0, theta1]))
        r_values_grid = np.linspace(r0, r1, N)
        total_riemann_integral = 0.0
        total_trapezoidal_integral = 0.0
        total_simpson_integral = 0.0
        total_romberg_integral = 0.0
        for i in tqdm(range(len(r_values_grid) - 1), desc="r values"):
            r0_cell = r_values_grid[i]
            r1_cell = r_values_grid[i + 1]
            for j in range(len(theta_values) - 1):
                theta0_cell = theta_values[j]
                theta1_cell = theta_values[j + 1]
                thetamid = 0.5 * (theta0_cell + theta1_cell)
                rmid = 0.5 * (r0_cell + r1_cell)
                # Unrefined and refined theta intervals for Romberg
                thetas_unrefined = np.array([theta0_cell, theta1_cell])
                thetas_refined = np.array([theta0_cell, thetamid, theta1_cell])
                r_values_unrefined = np.array([r0_cell, r1_cell])
                r_values_refined = np.array([r0_cell, rmid, r1_cell])
                refined_trapezoidal = 0.0
                for k in range(len(r_values_refined) - 1):
                    for M in range(len(thetas_refined)-1):
                        refined_trapezoidal += trapezoidal_integration(r_values=np.array([r_values_refined[k], r_values_refined[k+1]]),theta_values=np.array([thetas_refined[M], thetas_refined[M + 1]]),function=example_polar_function)
                # Compute integrations on this 2D cell
                riemann = riemann_integration(r_values=np.array([r0_cell, r1_cell]),theta_values=thetas_unrefined,function=example_polar_function)
                trapezoidal = trapezoidal_integration(r_values=np.array([r0_cell, r1_cell]),theta_values=thetas_unrefined,function=example_polar_function)
                simpson = simpson_13_integration(r_values=np.array([r0_cell, r1_cell]),theta_values=thetas_unrefined,function=example_polar_function)
                romberg = single_romberg_integration(trapezoidal, refined_trapezoidal)
                total_riemann_integral += riemann
                total_trapezoidal_integral += trapezoidal
                total_simpson_integral += simpson
                total_romberg_integral += romberg
        # Compute percent errors
        percent_error_riemann = abs(total_riemann_integral - analytic_integral) / abs(analytic_integral) * 100
        percent_error_trapezoidal = abs(total_trapezoidal_integral - analytic_integral) / abs(analytic_integral) * 100
        percent_error_simpson13 = abs(total_simpson_integral - analytic_integral) / abs(analytic_integral) * 100
        percent_error_romberg = abs(total_romberg_integral - analytic_integral) / abs(analytic_integral) * 100
        # Output results
        print("Analytic integration result:   ", analytic_integral)
        print("Riemann integration result:    ", total_riemann_integral)
        print("Riemann Percent Difference:    ", percent_error_riemann)
        print("Trapezoidal integration result:", total_trapezoidal_integral)
        print("Trapezoidal Percent Difference:", percent_error_trapezoidal)
        print("Simpson 1/3 integration result:", total_simpson_integral)
        print("Simpson 1/3 Percent Difference:", percent_error_simpson13)
        print("Romberg Integration Result:    ", total_romberg_integral)
        print("Romberg Percent Difference:    ", percent_error_romberg)


        




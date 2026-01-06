import numpy as np
import sys
import pygame 
from connection_m import connection
import matplotlib.pyplot as plt 
import json 
import helper as hlp # type: ignore
import time # type: ignore
from controller import get_errors, get_commanded
import matplotlib as mpl

class view_plane:
    """"""
    def __init__(self, json_file):
        self.viewplane_json_file = json_file
        self.parse_viewplane_json()

    def parse_viewplane_json(self):
            """This function reads the json file and stores the cylinder specific data in a dictionary"""
            with open(self.viewplane_json_file, 'r') as json_handle:
                input = json.load(json_handle)
                # appellian stuff
                self.dx = 0.01
                self.dy = 0.01
                self.distance_observe_to_viewplane = hlp.parse_dictionary_or_return_default(input, ["camera", "view_plane", "distance[ft]"], 1.0)
                self.observation_angle = np.deg2rad(hlp.parse_dictionary_or_return_default(input, ["camera", "view_plane", "angle[deg]"], 45.0))
                self.viewplane_RA = hlp.parse_dictionary_or_return_default(input, ["camera", "view_plane", "aspect_ratio"], 2.0)
                self.camera_location_xyz = np.array(hlp.parse_dictionary_or_return_default(input, ["camera", "location_xyz_from_vehicle[ft]"], [-200.0,0.0,0.0]))
                self.original_camera_location_xyz = np.array(self.camera_location_xyz)
                self.camera_orientation_phi_theta_psi = np.array(np.deg2rad(hlp.parse_dictionary_or_return_default(input, ["camera", "orientation_phi_theta_psi[deg]"], [0.0,0.0,0.0])))
                self.camera_quaternion = hlp.euler_to_quat(self.camera_orientation_phi_theta_psi)
                # ground grid
                self.ground_altitude = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "altitude[ft]"], 0.0)
                self.ground_grid_number = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "grid_number"], 0.0)
                self.ground_grid_scale = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "grid_scale[ft]"], 0.0)
                self.ground_grid_color = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "color"], 0.0)
                # vehicle stuff 
                self.vehicle_file = hlp.parse_dictionary_or_return_default(input, ["scene", "vehicle", "vtk_file"], "F16_coarse.vtk")
                self.vehicle_location_xyz = np.array(hlp.parse_dictionary_or_return_default(input, ["scene", "vehicle", "location_xyz[ft]"], [0.0,0.0,-2000.0]))
                self.camera_location_xyz += self.vehicle_location_xyz
                print("Vehicle Location\n", self.vehicle_location_xyz)
                print("Camera Location\n", self.camera_location_xyz)
                self.ground_line = []#ax.plot([], [], color=viewplane_object.ground_grid_color)
                self.vehicle_line = []#ax.plot([], [], color='black')  # or any color
                # controller stuff 
                self.vehicle_file = hlp.parse_dictionary_or_return_default(input, ["scene", "vehicle", "vtk_file"], "F16_coarse.vtk")
                self.is_use_controler = hlp.parse_dictionary_or_return_default(input, ["controller", "is_use_controller"], True)
                self.rollRateControl = hlp.parse_dictionary_or_return_default(input, ["controller", "rollRateControl"], [0.0,0.0,0.0,0.0])
                self.bankAngleControl = hlp.parse_dictionary_or_return_default(input, ["controller", "bankAngleControl"], [0.0,0.0,0.0,0.0])
                self.pitchRateControl = hlp.parse_dictionary_or_return_default(input, ["controller", "pitchRateControl"], [0.0,0.0,0.0,0.0])
                self.elevationAngleControl = hlp.parse_dictionary_or_return_default(input, ["controller", "elevationAngleControl"], [0.0,0.0,0.0,0.0])
                self.yawRateControl = hlp.parse_dictionary_or_return_default(input, ["controller", "yawRateControl"], [0.0,0.0,0.0,0.0])
                self.velocityControl = hlp.parse_dictionary_or_return_default(input, ["controller", "velocityControl"], [0.0,0.0,0.0,0.0])

    def parse_vtk(self):
        filename = self.vehicle_file
        with open(filename, 'r') as f:
            lines = f.readlines()
        points_idx = None
        for i, line in enumerate(lines):
            if line.startswith("POINTS"):
                points_idx = i
                break
        if points_idx is None:
            raise RuntimeError("VTK file missing POINTS section")
        # points
        parts = lines[points_idx].split()
        n_points = int(parts[1])
        vehicle_points = np.zeros((n_points, 3))
        for i in range(n_points):
            x, y, z = map(float, lines[points_idx + 1 + i].split())
            vehicle_points[i] = [x, y, z]
        # lines 
        lines_idx = None
        for i, line in enumerate(lines):
            if line.startswith("LINES"):
                lines_idx = i
                break
        if lines_idx is None:
            raise RuntimeError("VTK file missing LINES section")
        parts = lines[lines_idx].split()
        n_lines = int(parts[1])
        vehicle_lines = np.zeros((n_lines, 2), dtype=int)
        for k in range(n_lines):
            # 2 i j
            row = lines[lines_idx + 1 + k].split()
            if int(row[0]) != 2:
                raise RuntimeError("This parser only supports 2-point LINES entries (2 i j)")
            i1 = int(row[1])
            i2 = int(row[2])
            vehicle_lines[k] = [i1, i2]
        vehicle_points += self.vehicle_location_xyz
        self.vehicle_points = vehicle_points
        # self.vehicle_points += 
        self.vehicle_lines  = vehicle_lines
        self.vehicle_num_points = n_points
        self.vehicle_num_lines  = n_lines
        # 2D projected version — same style as ground grid
        self.vehicle_lines2D = np.full((self.vehicle_num_lines * 3, 2), None, dtype=object)

    def camera_set_state(self, camera_location, quat):
        """"""
        self.camera_location_xyz = camera_location
        self.camera_quaternion = quat 
        # Calculate the coordinates of the viewplane corners 
        self.calc_coordinates_of_viewplane_corners()
        # convert those to a 2D thing
        self.convert_3d_corner_coordinates_into_xy_plane()
        # calc corners in earth-fixed coordinates
        self.corners_in_earth_fixed()
        # calc PO
        self.calc_PO()
        # calc P1_P2
        self.calc_P1_P2()
        # calc_normal_vec_viewplane
        self.calc_normal_vec_viewplane()

    def calc_width_viewplane(self):
        """Calculates the width of the viewplane based on self.distance_observe_to_viewplane, self.observation_angle"""
        self.width_viewplane = 2*self.distance_observe_to_viewplane*np.tan(self.observation_angle) # Eq. 11.2.1 in flight sim book

    def calc_height_viewplane(self):
        """Calculates the height of the viewplane based on self.viewplane_RA and self.width_viewplane"""
        self.height_viewplane = self.width_viewplane/self.viewplane_RA

    def calc_coordinates_of_viewplane_corners(self):
        """Calculates the corners of the viewplane corners based on the viewplane distance away from the object, the viewplane height, and the viewplan width"""
        # create 3 4-element zero arrays for the x, y, and z coordinates of the 4 different corners of the viewplane going from top left, bottom left, bottom right, top right. (This is arbitrary, but should be kept track of)
        self.x_camera_viewplane_corners = np.array([self.distance_observe_to_viewplane, self.distance_observe_to_viewplane, self.distance_observe_to_viewplane, self.distance_observe_to_viewplane])
        self.y_camera_viewplane_corners = 0.5*np.array([-self.width_viewplane, -self.width_viewplane, self.width_viewplane, self.width_viewplane])
        self.z_camera_viewplane_corners = 0.5*np.array([-self.height_viewplane, self.height_viewplane, self.height_viewplane, -self.height_viewplane])
        # these coordinates are from Eq. 11.2.3 in the book
    
    def convert_3d_corner_coordinates_into_xy_plane(self):
        """This function converts the 3D coordinates such that they appear in an xy plane for visualization"""
        self.x_2D_corners = self.y_camera_viewplane_corners
        self.y_2D_corners = -self.z_camera_viewplane_corners
    
    def plot_viewplane_in_2D(self, lambda_array, points, lines, num_lines, projected_xy_array, lines2D, line_type):
        """plots the corners in 2D"""
        for i in range(num_lines):
            first_point = lines[i,0]
            second_point = lines[i,1]
            if lambda_array[first_point] > 0 and lambda_array[second_point] > 0:
                lines2D[3*i, :] = projected_xy_array[first_point,:]
                lines2D[3*i+1,:] = projected_xy_array[second_point,:]
            elif lambda_array[first_point] < 0 and lambda_array[second_point] > 0:
                lca = points[second_point]  - points[first_point] 
                if not np.isclose(np.dot(lca, self.norm_viewplane), 0.0, atol=1e-8): # prevents division by zero
                    lamb =  np.dot(self.P0 - points[first_point], self.norm_viewplane) / np.dot(lca, self.norm_viewplane) # this is lambda from eq. 11.3.5. The numerator has to be recomputed for the new Lca
                    pb = points[first_point] + lamb * lca                         # 3D intersection point on plane
                    lca_pb = pb - self.camera_location_xyz   # 3D vector from camera to pb
                    lambda_pb = self.lambda_numerator/np.dot(lca_pb, self.norm_viewplane)
                    Rotation_vec = self.rotation_matrix_earth_to_body(lambda_pb, lca_pb) # 1.0 because we want the full vector from camera to pb
                    projected_point = np.array([Rotation_vec[1], -Rotation_vec[2]])
                    # set the visible end to the intersection point, the other end to the visible point's projection
                    lines2D[3*i, :] = projected_point
                    lines2D[3*i+1, :] = projected_xy_array[second_point, :]
            elif lambda_array[first_point] > 0 and lambda_array[second_point] < 0:
                lca = points[second_point]  - points[first_point] 
                if not np.isclose(np.dot(lca, self.norm_viewplane), 0.0, atol=1e-8): # prevents division by zero
                    lamb =  np.dot(self.P0 - points[first_point], self.norm_viewplane) / np.dot(lca, self.norm_viewplane) # this is lambda from eq. 11.3.5. The numerator has to be recomputed for the new Lca
                    pb = points[first_point] + lamb * lca # 3D intersection point on plane
                    lca_pb = pb - self.camera_location_xyz   # 3D vector from camera to pb
                    lambda_pb = self.lambda_numerator/np.dot(lca_pb, self.norm_viewplane)
                    Rotation_vec = self.rotation_matrix_earth_to_body(lambda_pb, lca_pb) # 1.0 because we want the full vector from camera to pb
                    projected_point = np.array([Rotation_vec[1], -Rotation_vec[2]])
                    # set the visible end to the visible point's projection, the other end to the intersection point
                    lines2D[3*i, :] = projected_xy_array[first_point, :]
                    lines2D[3*i+1, :] = projected_point
        line_type.set_data(lines2D[:,0], lines2D[:,1])
        # if overwrite_data:
        #     self.ground_line.set_data(lines2D[:,0], lines2D[:,1])
        # else:
        #     self.vehicle_line.set_data(lines2D[:,0], lines2D[:,1])
            

    def rotation_matrix_for_body_fixed_to_earth_fixed(self, x_cp, y_cp, z_cp):
        """"""
        vec_camera = np.array([x_cp, y_cp, z_cp])
        return hlp.quat_dependent_to_base(vec_camera, self.camera_quaternion) + self.camera_location_xyz
    
    def corners_in_earth_fixed(self): 
        """"""
        first = self.rotation_matrix_for_body_fixed_to_earth_fixed(self.x_camera_viewplane_corners[0], self.y_camera_viewplane_corners[0], self.z_camera_viewplane_corners[0])
        second = self.rotation_matrix_for_body_fixed_to_earth_fixed(self.x_camera_viewplane_corners[1], self.y_camera_viewplane_corners[1], self.z_camera_viewplane_corners[1])
        third = self.rotation_matrix_for_body_fixed_to_earth_fixed(self.x_camera_viewplane_corners[2], self.y_camera_viewplane_corners[2], self.z_camera_viewplane_corners[2])
        fourth = self.rotation_matrix_for_body_fixed_to_earth_fixed(self.x_camera_viewplane_corners[3], self.y_camera_viewplane_corners[3], self.z_camera_viewplane_corners[3])
        self.x_camera_viewplane_corners_earth_fixed = np.array([first[0], second[0], third[0], fourth[0]])
        self.y_camera_viewplane_corners_earth_fixed = np.array([first[1], second[1], third[1], fourth[1]])
        self.z_camera_viewplane_corners_earth_fixed = np.array([first[2], second[2], third[2], fourth[2]])

    def rotation_matrix_earth_to_body(self, lambda_scalar, lca_element_vector):
        """"""
        return hlp.quat_base_to_dependent(lambda_scalar*lca_element_vector,self.camera_quaternion)

    def calc_PO(self):
        """"""
        x_avg_location = np.average(self.x_camera_viewplane_corners_earth_fixed)
        y_avg_location = np.average(self.y_camera_viewplane_corners_earth_fixed)
        z_avg_location = np.average(self.z_camera_viewplane_corners_earth_fixed)
        self.P0 = np.array([x_avg_location, y_avg_location, z_avg_location])

    def calc_P1_P2(self):
        """"""
        self.P1_array = np.array([self.x_camera_viewplane_corners_earth_fixed[3],self.y_camera_viewplane_corners_earth_fixed[3],self.z_camera_viewplane_corners_earth_fixed[3]])
        self.P2_array = np.array([self.x_camera_viewplane_corners_earth_fixed[0],self.y_camera_viewplane_corners_earth_fixed[0],self.z_camera_viewplane_corners_earth_fixed[0]])

    def calc_normal_vec_viewplane(self):
        """"""
        self.norm_viewplane = np.cross(self.P1_array-self.P0, self.P2_array-self.P0)
        self.P0_minus_Pc = self.calc_length_a_minus_c(self.P0)

    def calc_length_a_minus_c(self, point_a):
        return point_a - self.camera_location_xyz
    
    def calc_l_ca_array(self, array_of_points):
        """"""
        l_ca_array_points = np.copy(array_of_points)
        for i in range(len(array_of_points)):
            l_ca_array_points[i] = self.calc_length_a_minus_c(array_of_points[i])
        return l_ca_array_points
    
    def calc_lambda_numerator(self):
        """"""
        self.lambda_numerator = np.dot(self.P0_minus_Pc, self.norm_viewplane)

    def calc_lambda_denominator(self, lca):
        """"""
        return np.dot(lca, self.norm_viewplane)
    
    def calc_lambda_array(self, array_of_lca):
        """"""
        lambda_array = np.zeros(len(array_of_lca))
        for i in range(len(array_of_lca)):
            lambda_array[i] = self.lambda_numerator/(self.calc_lambda_denominator(array_of_lca[i]))
        return lambda_array
    
    def calc_xy_projection_onto_viewplane(self, array_of_lca, array_lambda):
        """"""
        x_y_array = np.zeros((len(array_of_lca), 2))
        for i in range(len(array_of_lca)):
            Rotation_vec = self.rotation_matrix_earth_to_body(array_lambda[i], array_of_lca[i])
            x_y_array[i] = [Rotation_vec[1], -Rotation_vec[2]]
        return x_y_array

    def calc_ground_grid(self):
        scale  = self.ground_grid_scale  # spacing between lines
        N  = self.ground_grid_number # number of positive/negative steps
        Z  = -self.ground_altitude
        # calculate grid number based on ground altitude
        # if abs(Z-self.camera_location_xyz[2]) < 100.0:
        #     scale = 150.0
        # elif 100.0 <= abs(Z-self.camera_location_xyz[2]) < 2000.0:
        #     scale = 0.10101*abs(Z-self.camera_location_xyz[2]) + 150.0
        # else: 
        #     scale = 0.10101*2000.0 + 150.0
        # Aircraft current x,y
        cam_x, cam_y = self.camera_location_xyz[:2]
        # Determine which grid cell aircraft sits in
        cx = np.floor(cam_x / scale)
        cy = np.floor(cam_y / scale)
        # The new "center cell" of the grid
        center_x = cx * scale
        center_y = cy * scale
        # Range of line coordinates
        offsets = (np.arange(-N, N+1) * scale)
        # Total number of lines in each direction
        ground_num_lines = 2 * N + 1
        self.gound_n_lines = 2 * ground_num_lines
        # Preallocate (same shape as original code)
        self.ground_points = np.zeros((4 * ground_num_lines, 3))
        self.ground_lines  = np.zeros((2 * ground_num_lines, 2), dtype=int)
        # X-parallel lines (horizontal) — vary in y
        for i in range(ground_num_lines):
            y = center_y + offsets[i]
            x_left  = center_x - N * scale
            x_right = center_x + N * scale
            self.ground_points[2*i, :]   = [x_left,  y, Z]
            self.ground_points[2*i+1, :] = [x_right, y, Z]
            self.ground_lines[i, :]      = [2*i, 2*i+1]
        # Y-parallel lines (vertical) — vary in x
        base = 2 * ground_num_lines
        for i in range(ground_num_lines):
            x = center_x + offsets[i]
            y_bottom = center_y - N * scale
            y_top    = center_y + N * scale
            self.ground_points[base + 2*i,   :] = [x, y_bottom, Z]
            self.ground_points[base + 2*i+1, :] = [x, y_top,    Z]
            self.ground_lines[ground_num_lines + i, :] = [base + 2*i, base + 2*i + 1]
        # Final bookkeeping
        self.ground_num_points = len(self.ground_points)
        self.ground_num_lines  = len(self.ground_lines)
        self.lines2D = np.full((self.ground_num_lines * 3, 2), None, dtype=object)

class HUD:
    def __init__(self,ax,camera):
        color = "black"
        box_background_color = 'lightgrey'
        dx = camera.dx
        dy = camera.dy
        self.altitude_box_offset = np.array([-0.2*dx, 0.0])
        self.velocity_box_offset = np.array([-0.61*dx, 0.0])
        self.heading_box_offset  = np.array([0.0, 0.27*dy])
        self.pilot_controls_box_offset = np.array([-0.4*dx, 0.0*dy])
        self.final_controls_box_offset = np.array([-0.4*dx, -0.1*dy])
        # Altitude
        self.altitude_minor = TickerTape(ax,color,'vertical', 0.35*dy,10,100,0.40*dx, -0.02*dx, True, perp_offset=0.2*dx, int_divid_by=500)
        # apply altitude offset
        a_ox, a_oy = self.altitude_box_offset
        alt_x = np.array([0.4*dx, 0.42*dx, 0.5*dx, 0.5*dx, 0.42*dx, 0.4*dx]) + a_ox
        alt_y = np.array([0.0, 0.03*dy,0.03*dy, -0.03*dy, -0.03*dy,0.0]) + a_oy
        ax.fill(alt_x.tolist(), alt_y.tolist(), facecolor=box_background_color,edgecolor=color,linewidth=1,zorder=100)
        self.altitude_box = ax.text(0.415*dx + a_ox, -0.02*dy + a_oy, str("{:0.0f}".format(0.0)), color=color,zorder=101)
        # Velocity
        self.velocity_minor = TickerTape(ax,color,'vertical', 0.35*dy,5,25,0.30*dx, -0.02*dx, True, label_dx=-0.02*dx, perp_offset=-0.21*dx, int_divid_by=50)
        # apply velocity offset
        v_ox, v_oy = self.velocity_box_offset
        vel_x = np.array([0.3*dx, 0.3*dx, 0.38*dx,0.4*dx, 0.38*dx]) + v_ox
        vel_y = np.array([0.03*dy, -0.03*dy, -0.03*dy,0.0, 0.03*dy]) + v_oy
        ax.fill(vel_x.tolist(), vel_y.tolist(), facecolor=box_background_color,edgecolor=color,linewidth=1,zorder=100)
        self.velocity_box = ax.text(0.315*dx + v_ox, -0.02*dy + v_oy, str("{:0.0f}".format(0.0)), color=color,zorder=101)
        # Heading
        self.heading_major = TickerTape(ax,color,'horizontal',0.2*dx,2,2.5,-0.48*dy, 0.05*dy,True, perp_offset=-0.2*dy, int_divid_by = 5)
        # apply heading offset
        h_ox, h_oy = self.heading_box_offset
        head_x = np.array([0.0, -0.02*dx, -0.02*dx, 0.02*dx, 0.02*dx, 0.0]) + h_ox # these are the coord
        head_y = np.array([-0.48*dy,-0.5*dy, -0.57*dy, -0.57*dy, -0.5*dy,-0.48*dy]) + h_oy
        ax.fill(head_x.tolist(), head_y.tolist(), facecolor=box_background_color,edgecolor=color,linewidth=1,zorder=100)
        self.heading_box = ax.text(-0.03*dx + h_ox+0.08, -0.5*dy + h_oy -0.05, str("{:0.0f}".format(0.0)),color=color,zorder=101)    
        # Pilot Controls Position
        p_ox, p_oy = self.pilot_controls_box_offset
        # pilot_x = np.array([-0.02*dx, -0.02*dx, 0.1*dx, 0.1*dx]) + p_ox
        # pilot_y = np.array([0.39*dy, 0.45*dy, 0.45*dy, 0.39*dy]) + p_oy
        # ax.fill(pilot_x.tolist(), pilot_y.tolist(), facecolor=box_background_color,edgecolor=color,linewidth=1,zorder=100)
        self.pilot_controls_box = ax.text(-0.015*dx + p_ox, 0.405*dy + p_oy, str("{:0.4f}".format(0.0)), color=color,zorder=101)
        # Final Controls Position
        f_ox, f_oy = self.final_controls_box_offset
        # final_x = np.array([-0.02*dx, -0.02*dx, 0.1*dx, 0.1*dx]) + f_ox
        # final_y = np.array([0.39*dy, 0.45*dy, 0.45*dy, 0.39*dy]) + f_oy
        # ax.fill(final_x.tolist(), final_y.tolist(), facecolor=box_background_color,edgecolor=color,linewidth=1,zorder=100)
        self.final_controls_box = ax.text(-0.015*dx + f_ox, 0.405*dy + f_oy, str("{:0.4f}".format(0.0)), color=color,zorder=101)

    def draw(self, camera, euler, u_velocity, pilot_controls, final_controls):
        dx = camera.dx
        dy = camera.dy
        # Altitude Ticker
        self.altitude_minor.update(-camera.camera_location_xyz[2])
        self.altitude_box.set_text(str("{:0.0f}".format(-camera.camera_location_xyz[2])))
        # Velocity Ticker 
        self.velocity_minor.update(u_velocity)
        self.velocity_box.set_text(str("{:0.0f}".format(u_velocity)))
        # Heading Ticker
        # self.heading_minor.update(euler[2]*180.0/np.pi)
        self.heading_major.update(euler[2]*180.0/np.pi,True)
        self.heading_box.set_text(str("{:0.0f}".format(euler[2]*180.0/np.pi)))
        # Pilot Controls Position using the first 3 elements of pilot_controls
        pd = np.rad2deg(pilot_controls[:3])
        pd_str = ", ".join(f"{x:0.4f}" for x in pd)
        self.pilot_controls_box.set_text("da de dr pilot = " + pd_str)
        # Final Controls Position (first 3 elements)
        fd = np.rad2deg(final_controls[:3])
        fd_str = ", ".join(f"{x:0.4f}" for x in fd)
        self.final_controls_box.set_text("da de dr final = " + fd_str)

class TickerTape:
    def __init__(self, ax, color,orientation,display_length,tick_spacing_display,tick_value_increment,x_pos, y_pos,major=False,label_dx=0.0,label_dy=0.0, perp_offset=0.0, int_divid_by = 1000):
        """
        ax: Matplotlib axes where everything is drawn
        color: tick and text color
        orientation: 'vertical' or 'horizontal'
        display_length: physical length of scrolling window (in axis units)
        tick_spacing_display: distance between tick marks on screen
        tick_value_increment: real-world units represented per tick spacing
        x_pos, y_pos: anchor point where the tape is drawn
        major: if True, draw longer ticks + labels
        label_dx, label_dy: label offset from tick marks
        """
        self.int_divid_by = int_divid_by
        self.ax = ax
        self.color = color
        self.orientation = orientation
        self.display_length = display_length
        # How far apart tick marks are (in screen units)
        self.tick_spacing_display = tick_spacing_display
        # Numerical increment between tick values
        self.tick_value_increment = tick_value_increment
        # Position on screen where the tape is anchored
        self.x_pos = x_pos
        self.y_pos = y_pos
        # Major tick control
        self.major = major
        self.label_dx = label_dx
        self.label_dy = label_dy
        # Perpendicular offset (axis data units): for horizontal tapes this
        # shifts the tape in Y; for vertical tapes it shifts the tape in X.
        # Can be changed after construction (e.g. `hud.heading_minor.perp_offset = 0.02*dx`).
        self.perp_offset = perp_offset
        # A list of matplotlib line and text objects
        self.tick_lines = []
        self.tick_labels = []
        # The most recent real-world value (altitude or heading)
        self.current_value = 0.0
        self.init_graphics()

    def init_graphics(self):
        """Create the initial set of ticks visible in the display region."""
        # Normalize `tick_spacing_display` into an effective spacing in axis
        # units. If the caller passed a large integer (e.g. 4 or 10) while
        # `display_length` is a small axis fraction (e.g. 0.2*dx), assume the
        # integer was intended as a number of ticks and convert it to spacing.
        if self.tick_spacing_display > self.display_length:
            # treat tick_spacing_display as a tick count
            tick_count = int(self.tick_spacing_display) if self.tick_spacing_display > 0 else 1
            self.spacing_display = max(self.display_length / tick_count, 1e-6)
        else:
            self.spacing_display = self.tick_spacing_display
        # Compute a starting anchor centered in the axes. The exact anchor is
        # recomputed each update so dynamic axis changes are handled; this
        # initial computation avoids a None state before the first update.
        if self.orientation == 'horizontal':
            xlim = self.ax.get_xlim()
            # center of the axes in data units
            self.x_pos = 0.5 * (xlim[0] + xlim[1])
            ylim = self.ax.get_ylim()
            self.y_pos = 0.5 * (ylim[0] + ylim[1]) + self.perp_offset
        else:
            ylim = self.ax.get_ylim()
            self.y_pos = 0.5 * (ylim[0] + ylim[1])
            xlim = self.ax.get_xlim()
            self.x_pos = 0.5 * (xlim[0] + xlim[1]) + self.perp_offset
        # draw ticks dynamically in update() so they can scroll and
        # appear/disappear at the display edges. Start with empty
        self.tick_lines = []
        self.tick_labels = []
        # half span of visible area
        self._half_span = 0.5 * self.display_length

    def update(self, new_value, update_labels=False):
        """
        Shift the ticker tape so that the tick corresponding to new_value
        appears centered in the box region.
        """
        self.current_value = new_value
        # Recompute center anchors each update so the tape stays centered when
        # the axes change size or limits and so `perp_offset` can be changed
        # at runtime and immediately observed.
        if self.orientation == 'horizontal':
            xlim = self.ax.get_xlim()
            self.x_pos = 0.5 * (xlim[0] + xlim[1])
            ylim = self.ax.get_ylim()
            self.y_pos = 0.5 * (ylim[0] + ylim[1]) + self.perp_offset
        else:
            ylim = self.ax.get_ylim()
            self.y_pos = 0.5 * (ylim[0] + ylim[1])
            xlim = self.ax.get_xlim()
            self.x_pos = 0.5 * (xlim[0] + xlim[1]) + self.perp_offset
        # Compute how far the tape should shift in display units so that the
        # tick corresponding to `new_value` is centered. This mirrors the
        # original behavior where `shift_ticks = new_value / tick_value_increment`.
        shift_ticks = new_value / float(self.tick_value_increment) if self.tick_value_increment != 0 else 0.0
        shift_display = shift_ticks * self.spacing_display
        
        # REMOVE ALL PREVIOUSLY VISIBLE LINES AND LABELS
        # This includes ticks that have scrolled out of the visible 
        # window as well as those that will be repositioned
        for ln in self.tick_lines:
            try:
                ln.remove() # remove is a matplotlib method to remove the line from the axes
            except Exception: # this happens if the line was already removed
                pass
        for lbl in self.tick_labels:
            try:
                lbl.remove()
            except Exception: # this happens if the label was already removed
                pass
        self.tick_lines = []
        self.tick_labels = []
        
        # COMPUTE RANGE OF TICKS CURRENTLY IN VISIBLE WINDOW
        # Determines which tick indices (m) fall within the display region
        spacing = self.spacing_display
        half = self._half_span
        if spacing <= 0:
            return
        m_min = int(np.floor((shift_display - half) / spacing)) # 
        m_max = int(np.ceil((shift_display + half) / spacing))
        
        # CREATE NEW LINES AND LABELS FOR ALL VISIBLE TICKS
        # This includes ticks that were previously visible (repositioned)
        # and ticks that are newly entering the visible window
        for m in range(m_min, m_max + 1):
            pos_offset = m * spacing - shift_display
            if self.orientation == 'vertical':
                x0 = self.x_pos
                x1 = self.x_pos + (0.01 if self.major else 0.005)
                y = self.y_pos + pos_offset
                line = self.ax.plot([x0, x1], [y, y], color=self.color, linewidth=1)[0]
                self.tick_lines.append(line)
                if self.major:
                    label_value = m * self.tick_value_increment
                    if label_value % self.int_divid_by == 0:
                        # align label based on label_dx sign: positive => left, negative => right
                        ha = 'left' if self.label_dx >= 0 else 'right'
                        label = self.ax.text(x1 + self.label_dx, y + self.label_dy, str(int(round(label_value))), color=self.color, ha=ha, va='center')
                        self.tick_labels.append(label)
            else:
                y0 = self.y_pos
                y1 = self.y_pos - (0.01 if self.major else 0.005)
                x = self.x_pos + pos_offset
                line = self.ax.plot([x, x], [y0, y1], color=self.color, linewidth=1)[0]
                self.tick_lines.append(line)
                if self.major:
                    label_value = m * self.tick_value_increment
                    if label_value % self.int_divid_by == 0:
                        label = self.ax.text(x + self.label_dx, y1 + self.label_dy, str(int(round(label_value))), color=self.color, ha='center', va='top')
                        self.tick_labels.append(label)

if __name__ == "__main__":
    mpl.rcParams['font.weight'] = 'bold'
    mpl.rcParams['axes.labelweight'] = 'bold'
    np.set_printoptions(formatter={'float': lambda x: f"{x:.12g}"})
    print_stuff = False
    plot_stuff = True
    # make viewplane object (which sets self.distance_observe_to_viewplane, self.observation_angle, and self.viewplane_RA)
    viewplane_object = view_plane("graphics.json")
    # Bank Angle settings
    desiredBankAngle = viewplane_object.bankAngleControl[0]
    integral_bank_angle = 0.0
    error_bank_angle = 0.5
    P_bankAngle = viewplane_object.bankAngleControl[1]
    I_bankAngle = viewplane_object.bankAngleControl[2]
    D_bankAngle = viewplane_object.bankAngleControl[3]
    # Roll Rate settings
    desiredRollRate = viewplane_object.rollRateControl[0]
    integral_roll_rate = 0.0
    error_roll_rate = 0.5
    P_rollRate = viewplane_object.rollRateControl[1]
    I_rollRate = viewplane_object.rollRateControl[2]
    D_rollRate = viewplane_object.rollRateControl[3]
    # Elevation Angle settings
    desiredElevationAngle = viewplane_object.elevationAngleControl[0]
    integral_elevation_angle = 0.0
    error_elevation_angle = 0.5
    P_elevationAngle = viewplane_object.elevationAngleControl[1]
    I_elevationAngle = viewplane_object.elevationAngleControl[2]
    D_elevationAngle = viewplane_object.elevationAngleControl[3]
    # Pitch Rate settings
    desiredPitchRate = viewplane_object.pitchRateControl[0]
    integral_pitch_rate = 0.0
    error_pitch_rate = 0.5
    P_pitchRate = viewplane_object.pitchRateControl[1]
    I_pitchRate = viewplane_object.pitchRateControl[2]
    D_pitchRate = viewplane_object.pitchRateControl[3]
    # Yaw Rate settings
    desiredYawRate = viewplane_object.yawRateControl[0]
    integral_yaw_rate = 0.0
    error_yaw_rate = 0.5
    P_yawRate = viewplane_object.yawRateControl[1]
    I_yawRate = viewplane_object.yawRateControl[2]
    D_yawRate = viewplane_object.yawRateControl[3]
    # Velocity settings
    desiredVelocity = viewplane_object.velocityControl[0]
    integral_velocity = 0.0
    error_velocity = 0.5
    P_velocity = viewplane_object.velocityControl[1]
    I_velocity = viewplane_object.velocityControl[2]
    D_velocity = viewplane_object.velocityControl[3]
    with open(viewplane_object.viewplane_json_file, 'r') as json_handle:
        file_loc = json.load(json_handle)
    states_connection = connection(file_loc["connections"]["receive_states"])
    controls_connection = connection(file_loc["connections"]["send_states"])
    pygame.init()
    pygame.joystick.init()
    joy = pygame.joystick.Joystick(0) #####
    joy.init() #####
    ## axis 0: Rudder
    ## axis 1: Throttle
    ## axis 2: Aileron 
    ## axis 3: elevator
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit() 
    states = [0]*14
    frame = 0
    fps = 0.0
    # parse vehicle vtk
    viewplane_object.parse_vtk()
    # calculate the width of the viewplane using the distance 
    viewplane_object.calc_width_viewplane()
    # calculate height of the view plane based on the width and aspect ratio of the viewplane
    viewplane_object.calc_height_viewplane()
    viewplane_object.camera_set_state(viewplane_object.camera_location_xyz, viewplane_object.camera_quaternion)
    viewplane_object.calc_ground_grid()
    l_ca_array = viewplane_object.calc_l_ca_array(viewplane_object.ground_points)
    l_ca_vehicle_array = viewplane_object.calc_l_ca_array(viewplane_object.vehicle_points)
    viewplane_object.calc_lambda_numerator()
    lambda_array = viewplane_object.calc_lambda_array(l_ca_array)
    lambda_vehicle_array = viewplane_object.calc_lambda_array(l_ca_vehicle_array)
    ground_xy_projected_on_viewplane = viewplane_object.calc_xy_projection_onto_viewplane(l_ca_array, lambda_array)
    vehicle_xy_projected_on_viewplane = viewplane_object.calc_xy_projection_onto_viewplane(l_ca_vehicle_array, lambda_vehicle_array)

    if print_stuff:
        print("")
        print("Distance observer to viewplane = ", viewplane_object.distance_observe_to_viewplane)
        print("Observation angle = ", np.rad2deg(viewplane_object.observation_angle))
        print("Viewplane aspect ratio = ", viewplane_object.viewplane_RA)
        print("Viewplane width = ", viewplane_object.width_viewplane)
        print("Viewplane height = ", viewplane_object.height_viewplane)
        print("x corners 3d = ", viewplane_object.x_camera_viewplane_corners)
        print("y corners 3d = ", viewplane_object.y_camera_viewplane_corners)
        print("z corners 3d = ", viewplane_object.z_camera_viewplane_corners)
        print("2D x corners = ", viewplane_object.x_2D_corners)
        print("2D y corners = ", viewplane_object.y_2D_corners)
        print("x corners earth fixed = ", viewplane_object.x_camera_viewplane_corners_earth_fixed)
        print("y corners earth fixed = ", viewplane_object.y_camera_viewplane_corners_earth_fixed)
        print("z corners earth fixed = ", viewplane_object.z_camera_viewplane_corners_earth_fixed)
        print("P0 = ", viewplane_object.P0)
        print("P1 = ", viewplane_object.P1_array)
        print("P2 = ", viewplane_object.P2_array)
        print("Normal Vec = ", viewplane_object.norm_viewplane)
        print("l_ca_array = \n", l_ca_array)
        print("lambda_array = \n", lambda_array)
        print("ground_xy_projected_on_viewplane = \n", ground_xy_projected_on_viewplane)
        print("")
    if plot_stuff:
        fig = plt.figure(figsize=(viewplane_object.viewplane_RA*5.0, 5.0))
        ax = fig.add_subplot(111)
        # Store the Axes object
        viewplane_object.ax = ax
        plt.subplots_adjust(top=1, bottom=0, left=0, right=1)
        plt.axis('off')
        # Create separate line handles
        viewplane_object.ground_line, = ax.plot([], [], color=viewplane_object.ground_grid_color)
        viewplane_object.vehicle_line, = ax.plot([], [], color='black')  # or any color
        ax.set_xlim(viewplane_object.x_2D_corners[0], viewplane_object.x_2D_corners[2])
        viewplane_object.dx = viewplane_object.x_2D_corners[2]-viewplane_object.x_2D_corners[0]
        ax.set_ylim(viewplane_object.y_2D_corners[1], viewplane_object.y_2D_corners[0])
        viewplane_object.dy = viewplane_object.y_2D_corners[0]-viewplane_object.y_2D_corners[1]
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect('equal')
        # Draw ground
        plt.show(block = False)
        viewplane_object.plot_viewplane_in_2D(lambda_array,viewplane_object.ground_points,viewplane_object.ground_lines,viewplane_object.ground_num_lines,ground_xy_projected_on_viewplane,viewplane_object.lines2D,line_type = viewplane_object.ground_line)
        # Draw vehicle
        viewplane_object.plot_viewplane_in_2D(lambda_vehicle_array,viewplane_object.vehicle_points,viewplane_object.vehicle_lines,viewplane_object.vehicle_num_lines,vehicle_xy_projected_on_viewplane,viewplane_object.vehicle_lines2D,line_type = viewplane_object.vehicle_line)
        # ax_HUD = fig.add_subplot(111)
        hud = HUD(ax, viewplane_object)
        # hud.ax = ax_HUD
        hud.draw(viewplane_object, [0.0,0.0,0.0], 0.0, [0.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0])
        # plt.show()
        dt = 0.001
        max_da_rad = 0.375246
        max_de_rad = 0.436332
        max_dr_rad = 0.523599
        max_tau = 1.0
        Controls = np.zeros(4)
        # At top of loop (or once before loop)
        plot_viewplane = viewplane_object.plot_viewplane_in_2D
        calc_l_ca_array = viewplane_object.calc_l_ca_array
        calc_lambda_array = viewplane_object.calc_lambda_array
        calc_xy_projection = viewplane_object.calc_xy_projection_onto_viewplane
        send_states = states_connection.recv
        send_controls = controls_connection.send
        calc_lambda_numerator = viewplane_object.calc_lambda_numerator
        calc_ground_grid = viewplane_object.calc_ground_grid
        while(frame<1000.0):
            time_start = time.time()
            states = states_connection.recv()
            # print("length_states = ", len(states))
            phiThetaPsi = hlp.quat_to_euler(states[10:14])
            plot_viewplane(lambda_array, viewplane_object.ground_points, viewplane_object.ground_lines, viewplane_object.ground_num_lines, ground_xy_projected_on_viewplane, viewplane_object.lines2D, viewplane_object.ground_line)
            # viewplane_object.plot_viewplane_in_2D(lambda_vehicle_array, viewplane_object.vehicle_points, viewplane_object.vehicle_lines, viewplane_object.vehicle_num_lines, vehicle_xy_projected_on_viewplane, viewplane_object.vehicle_lines2D)
            fig.canvas.draw()
            fig.canvas.flush_events() ####
            Controls[0] = np.round(-joy.get_axis(2)**3*max_da_rad,1)
            Controls[1] = np.round(-joy.get_axis(3)**3*max_de_rad,1)
            Controls[2] = np.round(joy.get_axis(0)**3*max_dr_rad,1)
            Controls[3] = np.round(-0.5*joy.get_axis(1)+0.5,1)
            dapilot = Controls[0]
            depilot = Controls[1]
            drpilot = Controls[2]
            taupilot = 0.0#Controls[3]
            alpha = np.arctan2(states[3], states[1])
            desiredElevationAngle = alpha
            # --- OUTER LOOP: bank and elevation angle PID - produces roll_rate_command and pitch_rate_command (setpoints for inner loop)
            integral_bank_angle, error_bank_angle, error_bank_P, error_bank_I, error_bank_D = get_errors(dt, desiredBankAngle, phiThetaPsi[0], P_bankAngle, I_bankAngle, D_bankAngle, integral_bank_angle, error_bank_angle)
            integral_elevation_angle, error_elevation_angle, error_elevation_P, error_elevation_I, error_elevation_D = get_errors(dt, desiredElevationAngle,phiThetaPsi[1],P_elevationAngle,I_elevationAngle,D_elevationAngle,integral_elevation_angle,error_elevation_angle)
            roll_rate_command = get_commanded(error_bank_P, error_bank_I, error_bank_D, dapilot, desiredRollRate, False)
            pitch_rate_command = get_commanded(error_elevation_P, error_elevation_I, error_elevation_D, depilot, desiredPitchRate, False)
            # Inner loop: Includes yaw and throttle (no inner loop stuff for this)
            integral_roll_rate, error_roll_rate, error_roll_P, error_roll_I, error_roll_D = get_errors(dt, roll_rate_command, states[4], P_rollRate, I_rollRate, D_rollRate, integral_roll_rate, error_roll_rate)
            integral_pitch_rate, error_pitch_rate, error_pitch_P, error_pitch_I, error_pitch_D = get_errors(dt, pitch_rate_command, states[5], P_pitchRate, I_pitchRate, D_pitchRate, integral_pitch_rate, error_pitch_rate)
            integral_yaw_rate, error_yaw_rate, error_yaw_P, error_yaw_I, error_yaw_D = get_errors(dt, desiredYawRate, states[6], P_yawRate, I_yawRate, D_yawRate, integral_yaw_rate, error_yaw_rate)
            integral_velocity, error_velocity, error_vel_P, error_vel_I, error_vel_D = get_errors(dt, desiredVelocity, states[1], P_velocity, I_velocity, D_velocity, integral_velocity, error_velocity)
            # final aileron command: when getting_final_command=True -> pilot overrides final actuator, otherwise use PID sum.
            Controls[0] = np.clip(get_commanded(error_roll_P, error_roll_I, error_roll_D, dapilot, roll_rate_command, True), -max_da_rad, max_da_rad) 
            Controls[1] = np.clip(get_commanded(error_pitch_P, error_pitch_I, error_pitch_D, depilot, pitch_rate_command, True), -max_de_rad, max_de_rad)
            Controls[2] = np.clip(get_commanded(error_yaw_P, error_yaw_I, error_yaw_D, drpilot, desiredYawRate, True), -max_dr_rad, max_dr_rad)
            Controls[3] = np.clip(get_commanded(error_vel_P, error_vel_I, error_vel_D, taupilot, desiredVelocity, True), 0.0, max_tau)
            Controls_sent = controls_connection.send(Controls)
            # print(states)
            hud.draw(viewplane_object, phiThetaPsi, states[1], [dapilot, depilot, drpilot, taupilot], Controls)
            # viewplane_object.vehicle_location_xyz[:] = states[7:10]
            viewplane_object.camera_location_xyz[:] = states[7:10] #+ viewplane_object.original_camera_location_xyz
            viewplane_object.camera_quaternion[:] = states[10:14]
            viewplane_object.camera_set_state(viewplane_object.camera_location_xyz, viewplane_object.camera_quaternion)
            calc_ground_grid()
            l_ca_array = calc_l_ca_array(viewplane_object.ground_points)
            # l_ca_vehicle_array = viewplane_object.calc_l_ca_array(viewplane_object.vehicle_points)
            calc_lambda_numerator()
            lambda_array = calc_lambda_array(l_ca_array)
            # lambda_vehicle_array = viewplane_object.calc_lambda_array(l_ca_vehicle_array)
            ground_xy_projected_on_viewplane = calc_xy_projection(l_ca_array, lambda_array)
            time_end = time.time()
            dt = time_end-time_start
            fps = 1/(dt)
            print("      update hz = ", np.round(fps,3), " alpha  = ", np.round(np.rad2deg(alpha),3))
            # print("        aileron = ", np.round(Controls[0]*180/np.pi,2), " elevator = ", np.round(Controls[1]*180/np.pi,2), " rudder = ", np.round(Controls[2]*180/np.pi,2), " throttle = ", np.round(Controls[3],2), "update hz = ", np.round(fps,3), " velocity: ", np.round(states[1],3), " desired V", np.round(desiredVelocity,1))
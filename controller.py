import numpy as np

def get_errors(dt, desired_value, actual_value, gain_P, gain_I, gain_D, integral, error): 
    # calculate the error
    newError = desired_value - actual_value
    # error in P is gain_P multiplied by that error
    error_P = gain_P*newError
    # integral using trapezoidal rule.
    integral = integral + 0.5 * (newError + error) * dt
    # protect from integral windup by setting bounds of what integral can be
    if integral > 0.5:
        integral = 0.5
    elif integral <-0.5: 
        integral = -0.5
    error_I = gain_I * integral # integral error
    error_D = gain_D * (newError-error)/dt # derivative error
    error = newError
    return integral, error, error_P, error_I, error_D

def get_commanded(error_P, error_I, error_D, pilot_yoke, command, getting_final_command):
    """""" 
    if getting_final_command:
        if abs(pilot_yoke) <= 0.01:
            controller_output = error_P + error_I + error_D
        else:
            controller_output = pilot_yoke
    else:
        if abs(pilot_yoke) <= 0.01:
            controller_output = error_P + error_I + error_D
        else:
            controller_output = command
    return controller_output

def PID_controller(states, dt, desired_value, actual_value, gain_P, gain_I, gain_D, integral, error,):
    """"""
    
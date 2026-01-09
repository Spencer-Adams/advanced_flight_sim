import pygame
import sys
import time
import numpy as np
from connection_m import connection

pygame.init()
pygame.joystick.init()

print("Looking for controllers...")

num_joy = pygame.joystick.get_count()
if num_joy == 0:
    print("No controllers found!")
    print("No controllers found!")
    print("No controllers found!")
    print("No controllers found!")
    print("No controllers found!")
    print("No controllers found!")
    sys.exit()

# Initialize the first joystick
joy = pygame.joystick.Joystick(0)
joy.init()
Controls = np.array([0.0,0.0,0.0,0.0])
print(f"Controller found: {joy.get_name()}")
print(f"  Axes: {joy.get_numaxes()}")
print(f"  Buttons: {joy.get_numbuttons()}")
print(f"  Hats: {joy.get_numhats()}")
print("\nMove sticks / press buttons to see output.\n")

# Main loop
clock = pygame.time.Clock()
# axis 0: Rudder
# axis 1: Throttle
# axis 2: Aileron 
# axis 3: elevator 
for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()

Controls = np.array([
    joy.get_axis(2),  # Aileron
    joy.get_axis(3),  # Elevator
    joy.get_axis(0),  # Rudder
    joy.get_axis(1),  # Throttle
])


while True:
    # Process pygame events
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()

    # Read axes
    axis_values = [joy.get_axis(i) for i in range(joy.get_numaxes())]
    # Read buttons
    button_values = [joy.get_button(i) for i in range(joy.get_numbuttons())]
    # Read hats (D-pad)
    hat_values = [joy.get_hat(i) for i in range(joy.get_numhats())]

    # Print current state
    print("\rAxes: " + str(axis_values) + 
          " | Buttons: " + str(button_values) +
          " | Hats: " + str(hat_values),
          end="")

    # Slow down the loop slightly
    clock.tick(30)
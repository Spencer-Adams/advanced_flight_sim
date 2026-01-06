import os
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import re
from mpl_toolkits.mplot3d import Axes3D # type: ignore
import pandas as pd

plt.rcParams["font.family"] = "Serif"
plt.rcParams["font.size"] = 17.0
plt.rcParams["axes.labelsize"] = 17.0
plt.rcParams['lines.linewidth'] = 1.0 # 1.0
plt.rcParams["xtick.minor.visible"] = True 
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.direction"] = plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.bottom"] = plt.rcParams["xtick.top"]= True 
plt.rcParams["ytick.left"] = plt.rcParams["ytick.right"] =True
plt.rcParams["xtick.major.width"] = plt.rcParams["ytick.major.width"] = 0.75
plt.rcParams["xtick.minor.width"] = plt.rcParams["ytick.minor.width"] = 0.75
plt.rcParams["xtick.major.size"] = plt.rcParams["ytick.major.size"] = 5.0
plt.rcParams["xtick.minor.size"] = plt.rcParams["ytick.minor.size"] = 2.5
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams['figure.dpi'] = 300.0
## change legend parameters
plt.rcParams["legend.fontsize"] = 17.0
plt.rcParams["legend.frameon"] = True
subdict = {"figsize" : (3.25,3.5),"constrained_layout" : True,"sharex" : True}
all_linestyles = [
    "solid",            # solid
    (0, (5, 5)),        # dashed
    "dashdot",          # dash-dot (standard Matplotlib style)
    (0, (3, 5, 1, 5)),  # dash-dot with spacing
    (0, (1, 5)),        # dotted (sparse)
    (0, (1, 3)),        # dotted (denser)
    (0, (1, 1))         # very fine dotted
]
color = ["black"]


##### First Plot Soccer ball height vs time for different initial velocities #####
file_location = "figures/"
x_is_log_scale = False
y_is_log_scale = False
labels = ["$V_0$ = 25.0 [ft/s]", "$V_0$ = 50.0 [ft/s]", "$V_0$ = 75.0 [ft/s]", "$V_0$ = 100.0 [ft/s]"] # labels for the different Velocities
is_legend = True
is_legend_below_fig = True
bbox_to_anchor = (0.5, -0.38)
x_axis_title = "time (s)"
y_axis_title = "Altitude (ft)"
labelpad = 1
x_limit = (0.0, 4.0)
x_ticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
x_tick_labels = ["0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "3.5", "4.0"]
is_move_x_tick_label_right = False
x_shift = 0.05
y_shift = -0.28
is_move_first_x_tick_label_right = True
y_limit = (0.0, 50.0)
y_ticks = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]
y_tick_labels = ["0", "10","20",  "30",  "40",  "50"]
is_show_plot = False
is_save_plot = True
plot_name = "Soccer_ball_Height_vs_Time_for_Different_Initial_Velocities"
first_column_for_plotting = 0
second_column_for_plotting = 9
number_of_rows_to_skip = 2
y_rotation = False
plot_title = "Ball Height vs Time for Different Initial Velocities"

def convert_excels_to_csv(folder):
    for filename in os.listdir(folder):
        if filename.lower().endswith(('.xls', '.xlsx')):
            excel_path = os.path.join(folder, filename)
            csv_path = os.path.splitext(excel_path)[0] + ".csv"
            try:
                df = pd.read_excel(excel_path)
                df.to_csv(csv_path, index=False)
                print(f"Converted {filename} to {os.path.basename(csv_path)}")
            except Exception as e:
                print(f"Failed to convert {filename}: {e}")

def read_file(file):
    """reads in a file and returns its content as a numpy array"""
    ext = os.path.splitext(file)[1].lower()

    # Convert Excel to CSV first
    if ext in ['.xls', '.xlsx']:
        convert_excels_to_csv(os.path.dirname(file))
        file = os.path.splitext(file)[0] + ".csv"

    # Decide delimiter based on extension
    if ext in [".csv"]:
        return np.loadtxt(file, delimiter=",", skiprows=number_of_rows_to_skip)
    else:
        # For .txt and others → use any whitespace as delimiter
        return np.loadtxt(file, skiprows=number_of_rows_to_skip)


def plot_files(file_location, first_column_for_plotting: int = 0, second_column_for_plotting: int = 1, is_legend: bool = True, is_legend_below_fig: bool = False, bbox_to_anchor: tuple = (0.5, -0.5), is_move_x_tick_label_right = False, y_rotation = True, plot_title: str = "", labelpad: int = 1):
    files = [f for f in os.listdir(file_location) if f.endswith(".csv") or f.endswith(".txt") or f.endswith(".xlsx")]
    files_sorted = sorted(files, key=lambda x: int(x.split('_')[0]))
    local_labels = labels.copy()
    n_lines = len(all_linestyles)
    for i, filename in enumerate(files_sorted):
        file_path = os.path.join(file_location, filename)
        data = read_file(file_path)
        if y_is_log_scale:
            y_data = data[:, second_column_for_plotting]
            # Use np.isclose to 0 with high precision (16 digits)
            y_data = np.where(np.isclose(y_data, 0.0, atol=1e-16), 2e-16, y_data)
            data = data.copy()
            data[:, second_column_for_plotting] = y_data
        if len(files_sorted) == len(local_labels):
            linestyle = all_linestyles[i % n_lines]
            if i < n_lines:
                plt.plot(
                    data[:, first_column_for_plotting],
                    -data[:, second_column_for_plotting],
                    label=local_labels[i],
                    color="black",
                    linestyle=linestyle,
                )
            else:
                plt.plot(
                    data[:, first_column_for_plotting],
                    -data[:, second_column_for_plotting],
                    label=local_labels[i],
                    color="gray",
                    linestyle=linestyle,
                )
        else:
            raise ValueError("Mismatch between number of files and labels")
    if y_rotation is True:
        plt.ylabel(y_axis_title, rotation=0)
    else:
        plt.ylabel(y_axis_title)
    plt.ylabel(y_axis_title, labelpad=labelpad)
    plt.xlabel(x_axis_title, labelpad=labelpad)
    plt.xlim(x_limit)
    plt.ylim(y_limit)
    if x_is_log_scale:
        plt.xscale("log")
        plt.xticks(x_ticks, labels=x_tick_labels)
    else:
        plt.xticks(x_ticks, labels=x_tick_labels)
    if y_is_log_scale:
        plt.yscale("log")
        plt.yticks(y_ticks, labels=y_tick_labels)
    else:
        plt.yticks(y_ticks, labels=y_tick_labels)
    if is_legend:
        handles, legend_labels = plt.gca().get_legend_handles_labels()
        if len(handles) == 7:
            # Reorder: group by rows instead of columns for 3 columns
            reordered = [
                handles[0], handles[3], handles[6],
                handles[1], handles[4],
                handles[2], handles[5]
            ]
            reordered_labels = [
                legend_labels[0], legend_labels[3], legend_labels[6],
                legend_labels[1], legend_labels[4],
                legend_labels[2], legend_labels[5]
            ]
            if is_legend_below_fig:

                plt.legend(
                    reordered, reordered_labels,
                    loc="lower center",
                    bbox_to_anchor= bbox_to_anchor,  # Move legend further down
                    ncol=3
                )
            else:
                plt.legend(reordered, reordered_labels, ncol=3)
        else:
            # Default legend order for any other number of curves
            if is_legend_below_fig:
                plt.legend(
                    loc="lower center",
                    bbox_to_anchor=bbox_to_anchor,  # Move legend further down
                    ncol=2
                )
            else:
                plt.legend(ncol=1)

    if is_move_x_tick_label_right:
        ax = plt.gca()
        ticklabels = ax.get_xticklabels()
        xticks = ax.get_xticks()
        if ticklabels and len(xticks) > 0:
            # Remove the first label
            new_labels = ["" if i == 0 else label.get_text() for i, label in enumerate(ticklabels)]
            ax.set_xticklabels(new_labels)
            # Add custom text at a shifted position
            # You can adjust the y value (vertical position) as needed
            # Use the same y as the original tick label
            x, y = ticklabels[0].get_position()
            ax.text(xticks[0] + x_shift * (xticks[1] - xticks[0]), y+y_shift, x_tick_labels[0], ha='center', va='center')
            # move yaxis title down a bit
            ax.yaxis.label.set_position((ax.yaxis.label.get_position()[0], ax.yaxis.label.get_position()[1] + 0.05))
    # plt.title(plot_title)
    # move the plot title up a bit
    plt.title(plot_title, pad=15)
    if is_save_plot:
        plt.savefig(os.path.join(file_location, f"{plot_name}.svg"), bbox_inches="tight")
    if is_show_plot:
        plt.show()

# run the plotting function
# plot_files(file_location, first_column_for_plotting=first_column_for_plotting, second_column_for_plotting=second_column_for_plotting, is_legend=is_legend, is_legend_below_fig=is_legend_below_fig, bbox_to_anchor=bbox_to_anchor, is_move_x_tick_label_right=is_move_x_tick_label_right, y_rotation=y_rotation, plot_title=plot_title, labelpad=labelpad)


# 3D version of plot_files: plot_3d_files
def plot_3d_files(
    file_location,
    x_column_for_plotting: int = 0,
    y_column_for_plotting: int = 1,
    z_column_for_plotting: int = 2,
    is_legend: bool = True,
    is_legend_below_fig: bool = False,
    bbox_to_anchor: tuple = (0.5, -0.5),
    plot_title: str = "3D Trajectory",
    labelpad: int = 1,
    skiprows: int = 2,
    delimiter=None,
    save_name: str = "3D_trajectory.svg",
    show: bool = True
):
    print("PLOTTING 3D FILES...")
    files = [f for f in os.listdir(file_location) if f.endswith(".csv") or f.endswith(".txt") or f.endswith(".xlsx")]
    files_sorted = sorted(files, key=lambda x: int(x.split('_')[0]))
    local_labels = labels.copy() if 'labels' in globals() else [f"File {i+1}" for i in range(len(files_sorted))]
    n_lines = len(all_linestyles) if 'all_linestyles' in globals() else 1
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i, filename in enumerate(files_sorted):
        file_path = os.path.join(file_location, filename)
        data = read_file(file_path)
        linestyle = all_linestyles[i % n_lines] if n_lines > 1 else 'solid'
        color_val = "black" if i < n_lines else "gray"
        label_val = local_labels[i] if len(files_sorted) == len(local_labels) else f"File {i+1}"
        ax.plot(
            data[:, x_column_for_plotting],
            data[:, y_column_for_plotting],
            -data[:, z_column_for_plotting],
            label=label_val,
            color=color_val,
            linestyle=linestyle
        )
    ax.set_xlabel('X', labelpad=labelpad)
    ax.set_ylabel('Y', labelpad=labelpad)
    ax.set_zlabel('Z', labelpad=labelpad)
    plt.title(plot_title, pad=15)
    if is_legend:
        if is_legend_below_fig:
            ax.legend(loc="lower center", bbox_to_anchor=bbox_to_anchor, ncol=1)
        else:
            ax.legend(ncol=1)
    if save_name:
        plt.savefig(os.path.join(file_location, save_name), bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)

# Example usage:
plot_3d_files(file_location, x_column_for_plotting=7, y_column_for_plotting=8, z_column_for_plotting=9)


# # 3D trajectory plotting function
# def plot_3d_trajectory(file_path, x_col, y_col, z_col, skiprows=1, delimiter=None, label=None, save_path=None, show=True):
#     """
#     Plots a 3D trajectory from a data file, allowing user to specify which columns are x, y, z.
#     Args:
#         file_path (str): Path to the data file.
#         x_col (int): Index of the x column (0-based).
#         y_col (int): Index of the y column (0-based).
#         z_col (int): Index of the z column (0-based).
#         skiprows (int): Number of rows to skip (for header). Default is 1.
#         delimiter: Delimiter for np.loadtxt. Default None (any whitespace).
#         label (str): Optional label for the trajectory.
#         save_path (str): If provided, saves the plot to this path.
#         show (bool): Whether to display the plot window.
#     """
#     data = np.loadtxt(file_path, skiprows=skiprows, delimiter=delimiter)
#     x = data[:, x_col]
#     y = data[:, y_col]
#     z = data[:, z_col]

#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot(x, y, z, label=label if label else 'Trajectory')
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     if label:
#         ax.legend()
#     if save_path:
#         plt.savefig(save_path, bbox_inches='tight')
#     if show:
#         plt.show()
#     plt.close(fig)

# Example usage:
# plot_3d_trajectory('yourfile.txt', x_col=7, y_col=8, z_col=9, skiprows=1)
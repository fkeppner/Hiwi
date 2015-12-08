__author__ = 'felix'
######
# Script to generate a .proc file for Marc Mentat. The .proc file will generate a lever structure.
# The Python script takes several input values and calculates the position of all necessary points
######


import math
import os
import sys

### Initial Values
theta = math.radians(20.)   # radians
t = 1.5                     # mm
L_ges = 20.                 # mm
h = 5.                      # mm
beta = math.radians(34.)    # radians
rep_x = 3                   # repetitions in x
rep_y = 5                   # repetitions in y
g = 4.                      # mm
output_file_name = "hebelstruktur.proc"



#### Obligatory things
out_file = os.getcwd() + os.sep + output_file_name
pt_id = 1  # Points are numbered continuous.
curve_id = 1  # Curves are also numbered continuous.


class Point:
    """Stores coordinates of points in three variables, x, y and z."""

    def __init__(self, x, y, z=0.0):
        global pt_id
        try:
            self.x = float(x)
            self.y = float(y)
            self.z = float(z)
            self.id = pt_id

            pt_id += 1
        except:
            sys.exit("Bad Coordinates.")


def loop_counterclock(points, conn):
    for p in range(len(points)):
        if p == len(points) - 1:
            conn.append((points[p], points[0]))
        else:
            conn.append((points[p], points[p + 1]))
    return conn


def copy_curves(delta_x, delta_y, rep, curves):
    """
    :param delta_x: displacement in x
    :param delta_y: displacement in y
    :param rep: repetitions
    :param curves: expects a string of curve ids or "all_existing"
    :return: multiple line string for the .proc file
    """
    string = """*set_duplicate_translation x %f
*set_duplicate_translation y %f
*set_duplicate_repetitions %d
*duplicate_curves
%s
#
""" % (delta_x, delta_y, rep, curves)
    return string


## Calculated Values
t_1 = t * 0.5
t_2 = t
L_y_exp = h + 2 * math.sin(theta) * math.sqrt((t_1 ** 2) / (1 - math.sin(beta) ** 2))
h_dash = h - 2 * math.sin(theta) * math.sqrt((t_1 ** 2) / (1 - math.sin(beta) ** 2))
### TODO: Gibt es einen Unterschied zwischen h und L_y_exp in der Zeichnung?
l = L_ges - math.sin(theta) * math.sqrt((t_1 ** 2) / (1 - math.sin(theta) ** 2))
a = l * math.sin(theta)
b = l * math.sqrt(1 - math.sin(theta) ** 2)
g_right = g  # Falls gewuenscht koennen unterschiedliche Werte fuer das rechte und linke g eingefuehrt werden.
g_left = g
x_s_left = g_left + t_2 * 0.5 + t_1
x_s_right = g_right + t_2 * 0.5 + t_1
y_s = h * 0.5 + a - t_1 - h_dash * 0.5

### Checking if g is bigger than possible
g_max = math.fabs(b - t_1 - t_2)
if g_left > g_max or g_right > g_max:
    sys.exit("The value of g is bigger than the maximum.")


# Calculating the coordinates of the unit cell nodes
unit_cell = []

unit_cell.append(Point(0, h * 0.5 + a - t_1))  # P1
unit_cell.append(Point(0, h * 0.5 + a))  # P2
unit_cell.append(Point(0, -1 * (h * 0.5 + a - t_1)))  # P3
unit_cell.append(Point(0, -1 * (h * 0.5 + a)))  # P4

unit_cell.append(Point(- b + t_1, h_dash * 0.5))  # P5
unit_cell.append(Point(- b, h * 0.5))  # P6
unit_cell.append(Point(- b + t_1, -1 * (h_dash * 0.5)))  # P7
unit_cell.append(Point(- b, -1 * (h * 0.5)))  # P8

unit_cell.append(Point(b - t_1, h_dash * 0.5))  # P9
unit_cell.append(Point(b, h * 0.5))  # P10
unit_cell.append(Point(b - t_1, - h_dash * 0.5))  # P11
unit_cell.append(Point(b, - h * 0.5))  # P12

# Creating the fillett (Steg)
g_p13 = - b + t_1 + g_left
g_p15 = g_p13 + t_2
archimedes = math.sin(theta) * math.sqrt(1 / (1 - math.sin(theta) ** 2))

unit_cell.append(Point(g_p13, 0))  # P13
unit_cell.append(Point(g_p15, 0))  # P14 = P15 in drawing

unit_cell.append(Point(b - t_1 - g_right, 0))  # P15 = P19 in drawing
unit_cell.append(Point(b - t_1 - t_2 - g_right, 0))  # P16 = P21 in drawing

unit_cell.append(Point(g_p13, h_dash * 0.5 + math.tan(theta) * g_left))  # P17 = P25 in drawing
unit_cell.append(Point(g_p15, h_dash * 0.5 + math.tan(theta) * (g_left + t_2)))  # P18 = P26 in drawing

unit_cell.append(Point(1 * unit_cell[16].x, -1 * unit_cell[16].y))  # P19 = P27 in drawing
unit_cell.append(Point(1 * unit_cell[17].x, -1 * unit_cell[17].y))  # P20 = P28 in drawing
unit_cell.append(Point(-1 * unit_cell[16].x, 1 * unit_cell[16].y))  # P21 = P29 in drawing
unit_cell.append(Point(-1 * unit_cell[17].x, 1 * unit_cell[17].y))  # P22 = P30 in drawing
unit_cell.append(Point(-1 * unit_cell[16].x, -1 * unit_cell[16].y))  # P23 = P31 in drawing
unit_cell.append(Point(-1 * unit_cell[17].x, -1 * unit_cell[17].y))  # P24 = P32 in drawing

###  Manually setting the connection for every point in counterclockwise mode
connections = []
# inner circle
inner_circle = [1, 18, 14, 20, 3, 24, 16, 22]
connections = loop_counterclock(inner_circle, connections)

# outer circle
outer_circle = [2, 6, 8, 4, 12, 10]
connections = loop_counterclock(outer_circle, connections)

# left small circe
left_small_circle = [17, 5, 7, 19, 13]
connections = loop_counterclock(left_small_circle, connections)

# right small circe
right_small_circle = [9, 21, 15, 23, 11]
connections = loop_counterclock(right_small_circle, connections)


def write_proc(path, points):
    file = open(path, "w")
    file.write("*set_model_length_unit meter\n")
    file.write("*add_points\n")
    for p in points:
        file.write("%.3f %.3f %.3f\n" % (p.x, p.y, p.z))
    file.write("*set_curve_type line\n")
    file.write("*add_curves\n")
    for c in connections:
        file.write("%d %d\n" % (c[0], c[1]))
    file.write("*fill_view\n")
    # Copy unit cell to intermediate layer
    file.write(copy_curves(b, (h + a), 1, "all_existing"))
    # Copy unit cells to all relative positions
    main_unit_cell = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24"
    file.write(copy_curves(2 * b, 0, rep_x - 1, main_unit_cell))
    intermed_unit_cell = "25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48"
    file.write(copy_curves(2 * b, 0, rep_x - 2,
                           intermed_unit_cell))
    main_layer = main_unit_cell
    intermed_layer = intermed_unit_cell
    counter = 48
    for k in range(1, rep_x):
        for i in range(1, 25):
            main_layer += " %s" % (i + counter)
        counter += 24
    for k in range(1, rep_x - 1):
        for i in range(1, 25):
            intermed_layer += " %s" % (i + counter)
        counter += 24
    file.write(copy_curves(0, 2 * (h + a), rep_y - 1,
                           main_layer))
    file.write(copy_curves(0, 2 * (h + a), rep_y - 2,
                           intermed_layer))

    file.write("*fill_view\n")

    file.write("*set_curve_div_type_variable_l1l2\n")
    file.write("*apply_curve_divisions\n")
    file.write("all_existing\n")
    file.write("*dt_planar_trimesh\n")
    file.write("all_existing\n")

    file.write("*fill_view\n")


write_proc(out_file, unit_cell)

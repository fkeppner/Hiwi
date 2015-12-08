__author__ = 'felix'
######
# Script to generate a .proc file for Marc Mentat. The .proc file will generate a lever structure.
# The Python script takes several input values and calculates the position of all necessary points
######

import math
import os
import sys


### Initial Values
theta = math.radians(20.)  # radians
t = 1.5  # mm
L_ges = 20.  # mm
h = 5.  # mm
beta = math.radians(34.)  # radians
rep_x = 5  # repetitions in x
rep_y = 3  # repetitions in y
g = 4.  # mm
output_file_name = "hebelstruktur.proc"

## Calculated Values
t_1 = t * 0.5
t_2 = t
L_y_exp = h + 2 * math.sin(theta) * math.sqrt((t_1 ** 2) / (1 - math.sin(beta) ** 2))
h_dash = h - 2 * math.sin(theta) * math.sqrt((t_1 ** 2) / (1 - math.sin(beta) ** 2))
### TODO: Gibt es einen Unterschied zwischen h und L_y_exp in der Zeichnung?
l = L_ges - math.sin(theta) * math.sqrt((t_1 ** 2) / (1 - math.sin(theta) ** 2))
a = l * math.sin(theta)
b = l * math.sqrt(1 - math.sin(theta) ** 2)
g_right = g
g_left = g
x_s_left = g_left + t_2 * 0.5 + t_1
x_s_right = g_right + t_2 * 0.5 + t_1
y_s = h * 0.5 + a - t_1 - h_dash * 0.5

#### Obligatory things
out_file = os.getcwd() + os.sep + output_file_name
pt_id = 1  # Points are numbered continuous.
curve_id = 1  # Curves are also numbered continuous.


def loop_counterclock(points, translation_dict):
    conn = list()
    for p in range(len(points)):
        if p == len(points) - 1:
            conn.append((translation_dict[points[p]], translation_dict[points[0]]))
        else:
            conn.append((translation_dict[points[p]], translation_dict[points[p + 1]]))
    return conn


def loop_out(points):
    conn_clock = list()
    for p in range(len(points) - 1, -1, -1):
        if p == 0:
            conn_clock.append((points[p], points[-1]))
        else:
            conn_clock.append((points[p], points[p - 1]))
    return conn_clock


def loop_gaps(points):
    # Lists inside lists
    conns = list()
    for gap in points:
        for p in range(len(gap)):
            if p == 0:
                conns.append((gap[p], gap[-1]))
            else:
                conns.append((gap[p], gap[p - 1]))
    return conns


### Checking if g is bigger than possible
g_max = math.fabs(b - t_1 - t_2)
if g_left > g_max or g_right > g_max:
    sys.exit("The value of g is bigger than the maximum.")


# Calculating the coordinates of all points
py_id = 1  # Point id inside this program
coord_dict = dict()  # Dictionary in which all points and their coordinates are stored
translation = dict()  # Dictionary which contains the mapping between py_id und marc_id
connections = list()


def calc_cell(coordinates, internal_id, base_x, base_y):
    coordinates.update({internal_id: [base_x, h * 0.5 + a - t_1 + base_y, 0]})  # P1
    internal_id += 1
    coordinates.update({internal_id: [base_x, h * 0.5 + a + base_y, 0]})  # P2
    internal_id += 1
    coordinates.update({internal_id: [base_x, -1 * (h * 0.5 + a - t_1) + base_y, 0]})  # P3
    internal_id += 1
    coordinates.update({internal_id: [base_x, -1 * (h * 0.5 + a) + base_y, 0]})  # P4
    internal_id += 1
    coordinates.update({internal_id: [- b + t_1 + base_x, h_dash * 0.5 + base_y, 0]})  # P5
    internal_id += 1
    coordinates.update({internal_id: [- b + base_x, h * 0.5 + base_y, 0]})  # P6
    internal_id += 1
    coordinates.update({internal_id: [- b + t_1 + base_x, -1 * (h_dash * 0.5) + base_y, 0]})  # P7
    internal_id += 1
    coordinates.update({internal_id: [- b + base_x, -1 * (h * 0.5) + base_y, 0]})  # P8
    internal_id += 1
    coordinates.update({internal_id: [b - t_1 + base_x, h_dash * 0.5 + base_y, 0]})  # P9
    internal_id += 1
    coordinates.update({internal_id: [b + base_x, h * 0.5 + base_y, 0]})  # P10
    internal_id += 1
    coordinates.update({internal_id: [b - t_1 + base_x, - h_dash * 0.5 + base_y, 0]})  # P11
    internal_id += 1
    coordinates.update({internal_id: [b + base_x, - h * 0.5 + base_y, 0]})  # P12
    internal_id += 1

    # Creating the fillett (Steg)
    g_p13 = - b + t_1 + g_left
    g_p15 = g_p13 + t_2

    coordinates.update({internal_id: [g_p13 + base_x, base_y, 0]})  # P13
    internal_id += 1
    coordinates.update({internal_id: [g_p15 + base_x, base_y, 0]})  # P14 = P15 in drawing
    internal_id += 1

    coordinates.update({internal_id: [b - t_1 - g_right + base_x, base_y, 0]})  # P15 = P19 in drawing
    internal_id += 1
    coordinates.update({internal_id: [b - t_1 - t_2 - g_right + base_x, base_y, 0]})  # P16 = P21 in drawing
    internal_id += 1

    p16_y = h_dash * 0.5 + math.tan(theta) * g_left
    p17_y = h_dash * 0.5 + math.tan(theta) * (g_left + t_2)

    coordinates.update({internal_id: [g_p13 + base_x, p16_y + base_y, 0]})  # P17 = P25 in drawing
    internal_id += 1
    coordinates.update({internal_id: [g_p15 + base_x, p17_y + base_y, 0]})  # P18 = P26 in drawing
    internal_id += 1

    coordinates.update({internal_id: [1 * g_p13 + base_x, -1 * p16_y + base_y, 0]})  # P19 = P27 in drawing
    internal_id += 1
    coordinates.update({internal_id: [1 * g_p15 + base_x, -1 * p17_y + base_y, 0]})  # P20 = P28 in drawing
    internal_id += 1
    coordinates.update({internal_id: [-1 * g_p13 + base_x, 1 * p16_y + base_y, 0]})  # P21 = P29 in drawing
    internal_id += 1
    coordinates.update({internal_id: [-1 * g_p15 + base_x, 1 * p17_y + base_y, 0]})  # P22 = P30 in drawing
    internal_id += 1
    coordinates.update({internal_id: [-1 * g_p13 + base_x, -1 * p16_y + base_y, 0]})  # P23 = P31 in drawing
    internal_id += 1
    coordinates.update({internal_id: [-1 * g_p15 + base_x, -1 * p17_y + base_y, 0]})  # P24 = P32 in drawing
    internal_id += 1

    return coordinates, internal_id


def trim_floats(data):
    # Cancels the floats to 3 digits after comma
    for m in data:
        data[m] = [round(data[m][0], 3), round(data[m][1], 3), round(data[m][2], 3)]
    return data


def find_twin_points(data):
    unique = list()
    result = dict()
    for key in data:
        if data[key] in unique:
            result.update({key: unique.index(data[key]) + 1})
            unique.append([0, 0, 0])
        else:
            unique.append(data[key])
            result.update({key: key})
    return result


def find_twin_curves(conn):
    unique = list()
    counter = 0
    doubles = 0
    for pair in conn:
        # ordered_pair = pair
        ordered_pair = (min(pair), max(pair))
        if ordered_pair not in unique:
            unique.append(ordered_pair)
            counter += 1
        else:
            doubles += 1
    # print(counter, doubles)
    return unique


def create_rim(internal_id):
    points = dict()
    ring = list()
    gaps = list()
    for rep in range(rep_x):
        # Lower Rim
        base_x = 2 * b * rep
        base_y = 0
        points.update({internal_id: [- b + base_x, -1 * (h * 0.5) + base_y - t_1, 0]})  # P8
        ring.append(internal_id)
        internal_id += 1
        points.update({internal_id: [base_x, -1 * (h * 0.5 + a) + base_y - t_1, 0]})  # P4
        ring.append(internal_id)
        internal_id += 1
        points.update({internal_id: [b + base_x, - h * 0.5 + base_y - t_1, 0]})  # P12
        ring.append(internal_id)
        internal_id += 1

    for rep in range(rep_y):
        # Right Side,
        base_x = 2 * b * (rep_x - 1)
        base_y = 2 * (h + a) * rep
        points.update({internal_id: [b + base_x + t_1, - h * 0.5 + base_y - t_1 * math.sin(beta) * math.sqrt(
            1 / (1 - math.sin(beta) ** 2)), 0]})  # P12
        ring.append(internal_id)
        internal_id += 1
        points.update({internal_id: [b + base_x + t_1,
                                     h * 0.5 + base_y + t_1 * math.sin(beta) * math.sqrt(1 / (1 - math.sin(beta) ** 2)),
                                     0]})  # P10
        ring.append(internal_id)
        internal_id += 1
        # Lever on the intermediate layer
        if rep != rep_y - 1:
            p17_y = h_dash * 0.5 + math.tan(theta) * (g_left + t_2)
            g_p15 = - b + t_1 + g_left + t_2
            base_x = 2 * b * (rep_x - 2) + b
            base_y = 2 * (h + a) * rep + h + a
            points.update({internal_id: [-1 * g_p15 + base_x + 2 * t_2 + 2 * g + 2 * t_1, -1 * p17_y + base_y,
                                         0]})  # P24 = P32 in drawing
            ring.append(internal_id)
            internal_id += 1
            points.update({internal_id: [-1 * g_p15 + base_x + 2 * t_2 + 2 * g + 2 * t_1, 1 * p17_y + base_y,
                                         0]})  # P22 = P30 in drawing
            ring.append(internal_id)
            internal_id += 1

    for rep in range(rep_x - 1, -1, -1):
        # Upper Rim
        base_x = 2 * b * rep
        base_y = 2 * (h + a) * (rep_y - 1)
        points.update({internal_id: [b + base_x, h * 0.5 + base_y + t_1, 0]})  # P10
        ring.append(internal_id)
        internal_id += 1
        points.update({internal_id: [base_x, h * 0.5 + a + base_y + t_1, 0]})  # P2
        ring.append(internal_id)
        internal_id += 1
        points.update({internal_id: [- b + base_x, h * 0.5 + base_y + t_1, 0]})  # P6
        ring.append(internal_id)
        internal_id += 1

    for rep in range(rep_y - 1, -1, -1):
        # Left Side
        base_x = 0
        base_y = 2 * (h + a) * rep
        points.update({internal_id: [- b + base_x - t_1,
                                     h * 0.5 + base_y + t_1 * math.sin(beta) * math.sqrt(1 / (1 - math.sin(beta) ** 2)),
                                     0]})  # P6
        ring.append(internal_id)
        internal_id += 1
        points.update({internal_id: [- b + base_x - t_1, -1 * (h * 0.5) + base_y - t_1 * math.sin(beta) * math.sqrt(
            1 / (1 - math.sin(beta) ** 2)), 0]})  # P8
        ring.append(internal_id)
        internal_id += 1
        # Lever on the intermediate layer
        if rep != 0:
            p17_y = h_dash * 0.5 + math.tan(theta) * (g_left + t_2)
            g_p15 = - b + t_1 + g_left + t_2
            base_x = 2 * b * 0 + b
            base_y = 2 * (h + a) * (rep - 1) + h + a
            points.update(
                {internal_id: [g_p15 + base_x - 2 * t_2 - 2 * g - 2 * t_1, p17_y + base_y, 0]})  # P18 = P26 in drawing
            ring.append(internal_id)
            internal_id += 1
            points.update({internal_id: [1 * g_p15 + base_x - 2 * t_2 - 2 * g - 2 * t_1, -1 * p17_y + base_y,
                                         0]})  # P20 = P28 in drawing
            ring.append(internal_id)
            internal_id += 1

    for rep in range(rep_y - 1):
        # Gap Holes - left
        base_x = 2 * b * 0 + b
        base_y = 2 * (h + a) * rep + h + a
        left = list()
        right = list()
        g_p13 = - b + t_1 + g_left
        p16_y = h_dash * 0.5 + math.tan(theta) * g_left

        points.update({internal_id: [g_p13 + base_x - 2 * g - t, p16_y + base_y, 0]})  # P17 = P25 in drawing
        left.append(internal_id)
        internal_id += 1
        points.update({internal_id: [1 * g_p13 + base_x - 2 * g - t, -1 * p16_y + base_y, 0]})  # P19 = P27 in drawing
        left.append(internal_id)
        internal_id += 1
        points.update({internal_id: [- b + t_1 + base_x - t, -1 * (h_dash * 0.5) + base_y, 0]})  # P7
        left.append(internal_id)
        internal_id += 1
        points.update({internal_id: [- b + t_1 + base_x - t, h_dash * 0.5 + base_y, 0]})  # P5
        left.append(internal_id)
        internal_id += 1
        gaps.append(left)

        # Gap Holes - right
        base_x = 2 * b * (rep_x - 2) + b
        base_y = 2 * (h + a) * rep + h + a

        points.update({internal_id: [-1 * g_p13 + base_x + 2 * g + t, 1 * p16_y + base_y, 0]})  # P21 = P29 in drawing
        right.append(internal_id)
        internal_id += 1
        points.update({internal_id: [-1 * g_p13 + base_x + 2 * g + t, -1 * p16_y + base_y, 0]})  # P23 = P31 in drawing
        right.append(internal_id)
        internal_id += 1
        points.update({internal_id: [b - t_1 + base_x + t, - h_dash * 0.5 + base_y, 0]})  # P11
        right.append(internal_id)
        internal_id += 1
        points.update({internal_id: [b - t_1 + base_x + t, h_dash * 0.5 + base_y, 0]})  # P9
        right.append(internal_id)
        internal_id += 1
        gaps.append(right)

    return points, ring, gaps


def write_proc(path, coordinates, conn, rim):
    file = open(path, "w")
    file.write("*set_model_length_unit meter\n")
    file.write("*add_points\n")
    for p in coordinates:
        file.write("%.3f %.3f %.3f\n" % (coordinates[p][0], coordinates[p][1], coordinates[p][2]))
    for p in rim:
        file.write("%.3f %.3f %.3f\n" % (rim[p][0], rim[p][1], rim[p][2]))
    file.write("*fill_view\n")
    file.write("*set_curve_type line\n")
    file.write("*add_curves\n")
    for c in conn:
        # pass
        file.write("%d %d\n" % (c[0], c[1]))
    file.write("*sweep_all\n")
    file.write("*remove_unused_points\n")
    file.write("*set_relative_tol 0.01\n")
    # file.write("*set_curve_div_type_variable_l1l2\n")
    # file.write("*apply_curve_divisions\n")
    # file.write("all_existing\n")
    # file.write("*dt_planar_trimesh\n")
    # file.write("all_existing\n")
    file.write("*fill_view\n")


# The coordinates of every point of every cell is calculated here
cells = 0
for repx in range(rep_x):
    for repy in range(rep_y):
        # Main Layer
        coord_dict, py_id = calc_cell(coord_dict, py_id, 2 * b * repx, 2 * (h + a) * repy)
        cells += 1
for repx in range(rep_x - 1):
    for repy in range(rep_y - 1):
        # Intermediate Layer
        coord_dict, py_id = calc_cell(coord_dict, py_id, 2 * b * repx + b, 2 * (h + a) * repy + h + a)
        cells += 1
trim_floats(coord_dict)
translation.update(find_twin_points(coord_dict))

for k in range(cells):
    #  Manually setting the connection for every cell in counterclockwise mode
    # inner circle
    factor = k * 24
    inner_circle = [1 + factor, 18 + factor, 14 + factor, 20 + factor, 3 + factor, 24 + factor, 16 + factor,
                    22 + factor]
    connections.extend(loop_counterclock(inner_circle, translation))

    # outer circle
    # outer_circle = [2 + factor, 6 + factor, 8 + factor, 4 + factor, 12 + factor, 10 + factor]
    # connections.extend(loop_counterclock(outer_circle, translation))

    # left small circe
    left_small_circle = [17 + factor, 5 + factor, 7 + factor, 19 + factor, 13 + factor]
    connections.extend(loop_counterclock(left_small_circle, translation))

    # right small circe
    right_small_circle = [9 + factor, 21 + factor, 15 + factor, 23 + factor, 11 + factor]
    connections.extend(loop_counterclock(right_small_circle, translation))


# Deleting double curves and adding the outer ring to the structure
conn_single = find_twin_curves(connections)
outer_rim, outer_ring, extra_gaps = create_rim(py_id)
conn_single.extend(loop_out(outer_ring))
conn_single.extend(loop_gaps(extra_gaps))
write_proc(out_file, coord_dict, conn_single, outer_rim)



def plot(data):
    # Plots all points, for debugging
    import matplotlib.pyplot as plt
    x_data = list()
    y_data = list()
    for m in data:
        x_data.append(data[m][0])
        y_data.append(data[m][1])
    plt.plot(x_data, y_data, '.')
    plt.title('Structure')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()


# coord_dict.update(outer_rim)
# plot(coord_dict)

print("Struktur in der Datei {} gespeichert.".format(output_file_name))

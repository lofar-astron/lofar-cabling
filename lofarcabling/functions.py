import sys
import numpy as np
import scipy.sparse
import shapely
from scipy.sparse import csgraph
import matplotlib.pyplot as plt
from numpy.linalg import norm
from scipy.spatial.distance import cdist
import os
from matplotlib.patches import Rectangle
from shapely.geometry import LineString, Polygon
import descartes
import csv
import math

def read_NMS_pts(ptscsv):
    """
    Reads out a set of points from .csv and puts them in a dictionary
    Args:
        ptscsv(str): filename of the .csv
    Rtrns:
        pts(dict): A dictionary with all the relevant information from the .csv
    """
    pts = {}
    dud_lines = ['', 'NAME', 'STATION-P', 'STATION-Q', 'STATION-R',
                 'NAME;STATION-P;STATION-Q;STATION-R']
    gud_chars = ['M', 'D', 'L', 'Q']
    with open(ptscsv) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            if (row[0][0] in gud_chars and row[0] not in dud_lines):
                pts[row[0]] = float(row[1]), float(row[2])

    return pts


def read_layout(layout):
    """
    Reads out a layout from .csv and puts it in a dictionary
    Args:
        layout(str): filename of the .csv
    Rtrns:
        pts(dict): A dictionary with all the relevant information from the .csv
    """
    pts = {}
    with open(layout) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            if row[0][0] == 'D':
                pts[row[0]] = float(row[1]), float(row[2])
            if row[0][0] == 'T':
                pts[row[0]] = row[1], row[2]
            if row[0][0] == 'U':
                pts[row[0]] = int(row[1]), int(row[2])
    return pts


def to_points_csv(points, ptsfile):
    """
    Writes a set of points to a .csv file like:
    ptsfile:
            'M1',19,20
            'M2',13,14
    Args:
        points(dict): dict of points {id: (x, y)}
        ptsfile(str): filename to write points csv to
    """
    pfile = open(ptsfile, 'w')
    for point in points:
        pfile.write(str(point) + ',')
        pfile.write(str(points[point][0]) + ',')
        pfile.write(str(points[point][1]) + '\n')
    pfile.close()


def to_layout_csv(layout, layoutfile):
    """
    Writes a layout to a .csv files like:
    lnsfile:
            'D1',12,34
            'D12','D13'
            'D14','D15'
    Args:
        layout(dict): dict of lines {id: (id first point, id second point)}
        layoutfile(str): filename to write lines csv to
    """
    lfile = open(layoutfile, 'w')
    for line in layout:
        lfile.write(str(line) + ',')
        lfile.write(str(layout[line][0]) + ',')
        lfile.write(str(layout[line][1]) + '\n')
    lfile.close()


def layout_matrix(points, layout, start_point):
    """
    creates a matrix based on points and a layout and a dictionary with
    the distance for every point to start_point.
    Args:
        points(dict): contains the cableorigin, and all antenna positions.
        layout(dict): contains all points and lines.
        start_point(int) : 0, if the cablehouse is set as the first of points.
    Rtrns:
        distance(dict): dict with distances for every point.
        mat: matrix based on points and layout.
    """
    dist_pid = 0
    dist_points = {}
    dist_lines = {}
    distance = {}
    for point in points:
        if point[0] != 'R':
            dist_points[point] = points[point][0], points[point][1], dist_pid
            dist_pid += 1
    for point in layout:
        if point[0] == 'D':
            dist_points[point] = layout[point][0], layout[point][1], dist_pid
            dist_pid += 1
    for line in layout:
        if line[0] == 'T':
            dist_lines[line] = layout[line][0], layout[line][1]
    mat = np.zeros((len(dist_points), len(dist_points)))
    dist_points_array = np.zeros((len(dist_points), 2))
    for pointId in dist_points:
        dist_points_array[dist_points[
                          pointId][2], :] = dist_points[pointId][:2]
    distmat = cdist(dist_points_array, dist_points_array)
    for lineId in dist_lines:
        mat[dist_points[dist_lines[lineId][0]][2],
            dist_points[dist_lines[lineId][1]][2]] = distmat[
            dist_points[dist_lines[lineId][0]][2],
                dist_points[dist_lines[lineId][1]][2]]
    for point in dist_points:
        if point[0] != 'D':
            dist = scipy.sparse.csgraph.dijkstra(
                    mat)[start_point, dist_points[point][2]]
            distance[point] = (dist)
    return (distance, mat)


def trench_tot(points, layout, point):
    """
    Calculates the total meters of digging that has been done for given field.
        Args:
        points(dict): contains the origin of all cables,
                      and all antenna positions.
        layout(dict): contains all points and lines that the path consists of.
        point(int): pointid of the cable origin (supposed to be on 0).
    Returns:
        (float) total trenches in meters.
    """
    processedlines = {}
    pid = 0
    pointsids = {}
    trenchtot = 0
    matrix = layout_matrix(points, layout, point)[1]

    for point in points:
        if point[0] == 'Q' or point[0] == 'M':
            pointsids[point] = points[point][0], points[point][1], pid
            pid += 1
    for l in layout:
        if l[0] == 'D':
            pointsids[l] = layout[l][0], layout[l][1], pid
            pid += 1

    for l in layout:
        if l[0] == 'T':
            for p in processedlines:
                if (layout[l][0] == processedlines[p][0] and
                    layout[l][1] == processedlines[p][1]) or \
                    (layout[l][0] == processedlines[p][1] and
                     layout[l][1] == processedlines[p][0]):
                    trenchtot -= scipy.sparse.csgraph.dijkstra(matrix)[
                                 pointsids[layout[l][0]][2],
                                 pointsids[layout[l][1]][2]]
            trenchtot += scipy.sparse.csgraph.dijkstra(matrix)[
                         pointsids[layout[l][0]][2],
                         pointsids[layout[l][1]][2]]
            processedlines[l] = layout[l][0], layout[l][1]
    return trenchtot


def cable_len(points, layout, point):
    """
    Calculates amount of cables that are long and short.
    Args:
        points(dict): contains the origin of all cables,
                      and all antenna positions.
        layout(dict): contains all points and lines that the path consists of.
        point(int): pointid of the cable origin (supposed to be on 0).
    Returns:
       (int) amount short cables, (int) amount long cables.
    """
    distances = layout_matrix(points, layout, point)[0]
    short = 0
    long = 0

    for d in distances:
        if d[0] != 'D':
            if distances[d] > 75:
                long += 1
            else:
                short += 1
    return short, long


def cost(points, layout):
    """
    Calculates € cost of a LOFAR field.

    Args:
        points(dict): contains the origin of all cables,
                      and all antenna positions.
        layout(dict): contains all points and lines that the path consists of.
    Returns:
        (int) Calculated value the field will be costing based on €.
    """
    c = cable_len(points, layout, 0)
    t = trench_tot(points, layout, 0)

    return(int(c[0]*170 + c[1]*230 + t*20))


def get_rad(points_csv):
    """
    only reads the radians given in a csv as:
                                R0 : {rad} : 0
    Args:
        points_csv(string): file name of a point set.
    Rtrns:
        rad(int): The integer that was read from the given .csv file.
    """
    with open(points_csv) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            if (row[0][0] == 'R'):
                rad = int(row[1])
    return rad


def center_field(points_csv, radians):
    """
    Changes given field to be centered around 'M0' at 0, 0.
    Then rotates the field clockwise.
    Args:
        points_csv(string): file name of a pointset.
        radians: Amount of degrees the field is to be rotated clockwise.
    Rtrns:
        points(dict): The rotated field.
    """
    rad = np.radians(radians)
    points = read_NMS_pts(points_csv)

    for p in points:
        pon = points[p]
        a = rotate_origin(pon, rad)
        points[p] = a
    return points


def center_layout(layout_csv, radians):
    """
    Rotates a layout given amount of degrees (rotates around 0,0).
    Args:
        layout_csv(string): file name of a layout.
        radians: amount of degrees the field is to be rotated.
    Rtrns:
        points(dict): The rotated layout.
    """
    rad = np.radians(radians)
    points = read_layout(layout_csv)
    for p in points:
        if p[0] == 'D':
            pon = points[p]
            a = rotate_origin(pon, rad)
            points[p] = a
    return points


def rotate_origin(xy, radians):
    """
    Rotate a point around the origin (0, 0).

    Args:
        xy: X and Y cordinates like (x, y).
        radians: np.radians(degrees rotation).
    Returns:
        xx: new x.
        yy: new y.
    """
    x, y = xy
    xx = x * math.cos(radians) + y * math.sin(radians)
    yy = -x * math.sin(radians) + y * math.cos(radians)

    return xx, yy


def draw_field(points, lay):
    """
    Plots a set of points and layout
    Args:
        points(dict): dictionary with points that are compatably formatted
        lay(dict): dictionary with a layout that is compatably formatted
    """
    minx = 999
    maxx = -999
    miny = 999
    maxy = -999
    fig, ax = plt.subplots(figsize=(8, 8))
    for point in points:
        if point[0] == 'M' or point[0] == 'L':
            ax.plot(points[point][0], points[point][1], 'k.', markersize=3)
            ax.add_artist(Rectangle(((float(points[point][0]) - 1.5),
                                     (float(points[point][1]) - 1.5)), width=3,
                                    height=3, facecolor='lightgray', zorder=0))
            for polygon in herrings(points[point]):
                ax.add_patch(descartes.PolygonPatch(polygon,  fc='darkorchid'))
        if point[0] == 'Q':
            ax.plot(points[point][0], points[point][1], 'r.')
            ax.add_artist(Rectangle(((float(points[point][0]) - 1.5),
                                     (float(points[point][1]) - 1.5)), width=3,
                                    height=3, facecolor='red', zorder=0))
    for l in lay:
        if l[0] == 'D':
            ax.plot(lay[l][0], lay[l][1], 'y.', markersize=3)
        elif l[0] == 'T':
            if lay[l][0][0] == 'D':
                fstx = lay[lay[l][0]][0]
                fsty = lay[lay[l][0]][1]
            else:
                fstx = points[lay[l][0]][0]
                fsty = points[lay[l][0]][1]
            if lay[l][1][0] == 'D':
                scdx = lay[lay[l][1]][0]
                scdy = lay[lay[l][1]][1]
            else:
                scdx = points[lay[l][1]][0]
                scdy = points[lay[l][1]][1]
            ax.plot((fstx, scdx), (fsty, scdy), color='black',
                    linewidth=0.5, zorder=0.5)
    for point in points:
        if point[0] != 'R':
            if points[point][0] < minx:
                minx = points[point][0]
            if points[point][0] > maxx:
                maxx = points[point][0]
            if points[point][1] < miny:
                miny = points[point][1]
            if points[point][1] > maxy:
                maxy = points[point][1]
    ax.set_xlim(minx-5, maxx+5)
    ax.set_ylim(miny-5, maxy+5)
    ax.set_aspect(1)


def run_layouts(pts):
    """
    tries every possible loadout and saves the loadout with the lowest cost.
    Args:
        Pts(dict): dict containing 'Q1' as cableorigin,
        'M0' - 'M95' as antennas and 'R0' as rotation.
    Rtrns:
        layout(dict): the layout that gave the best result paired with pts.
    """
    location = (file_prefix() + '/share/lofarcabling/layouts' + os.path.sep)
    layouts = [[location + 'layout602.csv'], [location + 'layout603.csv'],
               [location + 'layout604.csv'], [location + 'layout605.csv'],
               [location + 'layout606.csv'], [location + 'layout607.csv'],
               [location + 'layout608.csv'], [location + 'layout609.csv'],
               [location + 'layout610.csv'], [location + 'layout611.csv'],
               [location + 'layout612.csv'], [location + 'layout613.csv'],
               [location + 'layout614.csv'], [location + 'layout615.csv'],
               [location + 'layout616.csv'], [location + 'layout617.csv'],
               [location + 'layout618.csv'], [location + 'layout619.csv'],
               [location + 'layout620.csv'], [location + 'layout621.csv'],
               [location + 'layout622.csv'], [location + 'layout623.csv']]

    bestcost = 999999
    costs = []
    bestlayout = ""

    for l in layouts:
        rad = get_rad(l[0])
        zero = 360 - rad
        lay = center_layout(l[0], zero)

        c = cost(pts, lay)
        costs.append(c)
        if c < bestcost:
            bestcost = c
            bestlayout = l[0]
    return bestcost, bestlayout, costs


def herrings(pt):
    """
    Give 8 herring locations for inserted antenna location
    Args:
        pt: id of the requested point
    Returns:
        herrings(array): array of the herring locations
   """
    herrings = []
    oneherring = np.array([[-0.225, 0.225], [0.225, 0.225],
                          [0.225, -0.225], [-0.225, -0.225]])
    for xoffset in [-1.5+0.225, 0, 1.5-0.225]:
        for yoffset in [-1.5+0.225, 0, 1.5-0.225]:
            if xoffset == 0 and yoffset == 0:
                continue
            offset = np.array([xoffset, yoffset])
            herrings += [Polygon(np.array(pt) + offset + oneherring)]
    return herrings


def herring_intersections(points, layout):
    """
    finds all points where lines intersect with herrings
    Args:
        points(str): file location with .csv of points
        layout(str): file location with .csv of layout
    """
    pts = points
    lay = layout
    intersections = {}
    intersection_id = 0
    for l in lay:
        if l[0] == 'T':
            if lay[l][0][0] == 'D':
                startp = lay[lay[l][0]]
            else:
                startp = pts[lay[l][0]]

            if lay[l][1][0] == 'D':
                endp = lay[lay[l][1]]
            else:
                endp = pts[lay[l][1]]

            line = LineString((startp, endp))
            for point in pts:
                if point[0] == 'M':
                    for her in herrings(pts[point]):
                        if line.intersection(her).length >= 0.025:
                            for mok in lay:
                                if lay[mok][1] == lay[l][0]:
                                    intersections[intersection_id] = \
												point, pts[point], \
                                                line.intersection(her).length, \
                                                lay[l][0], lay[mok][0]
                                    intersection_id += 1
    return intersections


def fix_herrings(points, layout):
    """
    Changes all the lines in the layout that intersect with the
    herrings that belong to the antennas, changes them so
    they no longer intersect.
    Args:
        points(str): file location with .csv of points
        layout(str): file location with .csv of layout
    """
    intersect = herring_intersections(points, layout)
    x, y = 0, 0
    tot = 0
    cur = -1
    for i in intersect:
        if i > tot:
            tot = i

    for l in layout:
        if l[0] == 'D':
            for t in layout:
                if t[0] == 'T':
                    if layout[t][1] == l:
                        intersect = herring_intersections(points, layout)
                        for itrsct in intersect:
                            if itrsct in intersect:
                                if layout[t][1] == intersect[itrsct][3] and layout[t][0] == intersect[itrsct][4]:
                                    x_list, y_list = [], []
                                    k = 0.075
                                    i = 0
                                    if layout[t][0][0] == 'D' and \
                                       layout[t][1][0] == 'D':
                                        f = ((layout[intersect[itrsct][3]][0] -
                                             layout[intersect[itrsct][4]][0]),
                                             (layout[intersect[itrsct][3]][1] -
                                             layout[intersect[itrsct][4]][1]))
                                    elif layout[t][0][0] == 'D' and (
                                       layout[t][1][0] != 'D'):
                                        f = ((points[intersect[itrsct][3]][0] -
                                             layout[intersect[itrsct][4]][0]),
                                             (points[intersect[itrsct][3]][1] -
                                             layout[intersect[itrsct][4]][1]))
                                    elif layout[t][0][0] != 'D' and (
                                       layout[t][1][0] == 'D'):
                                        f = ((layout[intersect[itrsct][3]][0] -
                                             points[intersect[itrsct][4]][0]),
                                             (layout[intersect[itrsct][3]][1] -
                                             points[intersect[itrsct][4]][1]))
                                    while k <= 0.825:
                                        x_list.append(layout[intersect[
                                                      itrsct][3]][0] +
                                                      f[0]*(0-k))
                                        y_list.append(layout[intersect[
                                                      itrsct][3]][1] +
                                                      f[1]*(0-k))
                                        x_list.append(layout[intersect[
                                                      itrsct][3]][0] +
                                                      f[0]*(0+k))
                                        y_list.append(layout[intersect[
                                                      itrsct][3]][1] +
                                                      f[1]*(0+k))
                                        k += 0.075
                                    unsolved, stillhere = True, True
                                    many = 0
                                    while unsolved:
                                        if many > 10:
                                            unsolved = False
                                            cur += 1
                                            many = 0
                                            break
                                        if i <= 21:
                                            templayout = layout
                                            x, y = x_list[i], y_list[i]
                                            i += 1
                                            templayout[l] = (x, y)
                                            checkinct = herring_intersections(
                                                        points, templayout)
                                            if checkinct != {}:
                                                stillhere = False
                                                for c in checkinct:
                                                    if checkinct[c][3] == l:
                                                        stillhere = True
                                                if not stillhere:
                                                    layout = templayout
                                                    unsolved = False
                                                    intersect = checkinct
                                                    cur += 1
                                            else:
                                                layout = templayout
                                                unsolved = False
                                                intersect = checkinct
                                                cur += 1
                                        else:
                                            p = 0
                                            for x in x_list:
                                                x_list[p] = x_list[p] + 0.075
                                                y_list[p] = y_list[p] + 0.075
                                                p += 1
                                            i = 0
                                            many += 1


def clean_up(layout):
    """
    Removes useless points from a layout
    Args:
        layout: dict of points and lines
    Returns:
        the same layout with less points different lines
    """
    firstenc = {}
    secondenc = {}
    poplist = {}
    for d in layout:
        count = 0
        if d[0] == 'D':
            for l in layout:
                if l[0] == 'T':
                    if (layout[l][1] == d or layout[l][0] == d):
                        count += 1
                        if layout[l][1] == d:
                            firstenc = layout[l]
                        elif layout[l][0] == d:
                            secondenc = layout[l]
            if count == 2 and secondenc[1][0] == 'D' and firstenc[0][0] != 'Q':
                line1 = np.array(layout[firstenc[0]]) - \
                    np.array(layout[firstenc[1]])
                line2 = np.array(layout[secondenc[0]]) - \
                    np.array(layout[secondenc[1]])
                costheta = line1.dot(line2) / \
                    (np.linalg.norm(line1) * np.linalg.norm(line2))
                if np.rad2deg(np.arccos(costheta)) < 5:
                    poplist[d] = d
                    for l in layout:
                        if l[0] == 'T':
                            if layout[l] == firstenc:
                                layout[l] = (firstenc[0], secondenc[1])
                            if layout[l] == secondenc:
                                poplist[l] = l
    for pop in poplist:
        layout.pop(pop)
    return layout


def find_layout(points_csv, fixherrings):
    """
    finds the best (currently known) layout for the given set of antennas.

    Args:
        points_csv like:    Q1, x, y
                            M0, x, y
                            M1, x, y
                            ........
                            R1, rad, 0
        where q is the cablehouse/point of origin,
        where m is an antenna,
        where r is degrees this field is rotated.
    Rtrns:
        plot of points with best layout,
        name of best layout, and cost for reference.
    """
    radpts = get_rad(points_csv)
    zerop = 360-radpts
    points = center_field(points_csv, zerop)

    result = run_layouts(points)
    radlay = get_rad(result[1])

    if radlay > radpts:
        raddiff = 360-radlay + radpts
    elif radpts > radlay:
        raddiff = radpts - radlay
    else:
        raddiff = 0

    points = center_field(points_csv, 0)
    lay = center_layout(result[1], raddiff)

    if(fixherrings):
        end = fix_herrings(points, lay)

    lay = clean_up(lay)
    draw_field(points, lay)
    return result, lay


def go(rot, loc, namepoints, namelayout, fixherrings):
    """
    Finds a layout for the field that has a rotation "rot" and startpoint "loc"
    Saves a copy of the field under the given name
    Args:
        rot(int): rotation in degrees
        loc(list): x and y coordinates of the startpoint like : (10,10)
        namepoints(str): desired filelocation for points
        namelayout(str): desired filelocation for layout
        fixherrings(bool): True if you want to fix herrings,
        False if you want to do it manually
    Rtrns:
         plot of points with best layout, name of best layout,
         and cost for reference.
    """
    location = (file_prefix() + os.path.sep)
    field = center_field((location +
                         '/share/lofarcabling/layouts/examplefield.csv'),
                         rot)
    field['Q1'] = loc
    field['R0'] = (rot, 0)

    to_points_csv(field, namepoints)

    layout = find_layout(namepoints, fixherrings)[1]

    to_layout_csv(layout, namelayout)


def file_prefix():
    """
    Find the prefix for /lofar-cabling
    Rtrns:
        the prefix as a string
    """
    filepath = os.path.split(__file__)
    filepath = os.path.split(filepath[0])
    return filepath[0]

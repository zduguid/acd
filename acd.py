# /usr/local/bin/python3
#
# - Implementation of the Approximate Convex Decomposition (ACD) algorithm
# - Algorithm developed Jyh-Ming Lien and Nancy Amato
# - Original paper can be found here: http://masc.cs.gmu.edu/wiki/ACD
#
# Author: Zach Duguid
# Last Updated: 10/29/2018

import re
import sys
import pygeoj
import numpy as np


def ACD(input_polygon, tau, length):
    """
    recursively derive the approximate convex decomposition of the input polygon
    :param input_polygon: the input polygon that is to be decomposed
        :requires: the points in the polygon are given in clockwise order
    :param tau: the non-dimensional tolerance parameter that allows for some level of concavity
    :param length: the length scale of the polygon
    :returns: a list of CW polygons that form the decomposition of the input polygon
    """
    # reformat the input polygon to work with ACD algorithm
    polygon = get_reformatted_polygon(input_polygon, length)
    assert len(polygon) >= 4,         'polygon needs at least 4 vertices (has %d)' % len(polygon)
    assert polygon[0] == polygon[-1], 'polygon needs to be closed'
    tau = min(max(tau, 0), 1)

    # get the notch point with greatest concavity from the input polygon
    notch, point_dict, hull = get_concavity(polygon)
    concavity = point_dict[notch]
    
    # return the input polygon if the concavity tolerance is already satisfied
    if concavity <= (tau*length):
        hull.append(hull[0])
        return [hull]

    # resolve the notch associated with the greatest concavity 
    resolved_polygons = get_resolved_polygons(polygon, notch, point_dict)

    # collect all polygons created during the decomposition
    output_polygons = [] 

    # recursively call the approximate convex decomposition algorithm
    for resolved_polygon in resolved_polygons:
        new_polygons = ACD(resolved_polygon, tau, length)
        output_polygons.extend(new_polygons)

    return output_polygons


def get_concavity(polygon):
    """
    calculates the concavity measure of a polygon
    :param polygon: the input polygon to be measured
    :returns: a float that represents the concavity measure of the polygon
    """
    # define open polygon where first and last points are now distinct  
    polygon_open = polygon[:-1]

    # get the convex hull of the polygon 
    #   + convex hull starts with bottom left point first
    #   + convex hull follows clockwise ordering
    hull = get_convex_hull(polygon_open)
    point_dict = {}

    # initialize bride references
    #   + bridge refers to points on convex hull that bound any non-hull points
    left_bridge_index = -1
    right_bridge_index = 0

    # get the concavity measure of each point
    #   + measured by the minimum distance between the point and its corresponding bridge points
    for point in polygon_open:

        # calculates the concavity of all non-hull points
        if point not in hull:
            point_dict[point] = min(get_distance(point, hull[left_bridge_index]), 
                                    get_distance(point, hull[right_bridge_index]))

        # hull points have zero concavity
        else: 
            point_dict[point] = 0

            # advance the bridge indices appropriately
            left_bridge_index += 1
            right_bridge_index += 1
            if right_bridge_index == len(hull) - 1:
                right_bridge_index = 0

    # determine the maximum concavity point of the polygon
    #   + this provides us with the notch that must be resolved when decomposing
    concavity = max(point_dict.values())
    notch     = max(point_dict, key=point_dict.get)

    return notch, point_dict, hull


def get_resolved_polygons(polygon, notch, point_dict):
    """
    decomposes the input polygon into two sub-polygons that resolve the notch, if possible
    :param polygon: the input polygon to be resolved
    :param notch: the point that is to be resolved during decomposition
    :returns: a list of sub-polygons that make up the input polygon
    """
    # define open polygon where first and last points are now distinct  
    polygon_open = polygon[:-1]

    # constants for computing the score value
    Sc = 0.1
    Sd = 1.0
    best_score = 0
    best_point = None
    best_polygons = None

    # iterate through possible points to find the highest score
    for point in polygon_open:
        if point != notch:

            # score represents how good a particular split would be
            score = (1 + Sc*point_dict[point]) / (Sd + get_distance(point, notch))
 
            # if score is greater than current best score, assess validity of split
            if score > best_score:
                resolved_polygons = get_cut_polygon(polygon_open, notch, point)

                # update the current best diagonal if the notch is valid
                if get_valid_resolve(resolved_polygons, notch):
                    best_polygons = resolved_polygons
                    best_score = score
                    best_point = (point[0], point[1])

    # default behavior returns the closed convex hull of the input polygon 
    if best_polygons == None:
        hull = get_convex_hull(polygon_open)
        hull.append(hull[0])
        return [hull]

    # return the closed polygons after they have been resolved
    closed_polygons = []
    for best_polygon in best_polygons:
        best_polygon.append(best_polygon[0])
        closed_polygons.append(best_polygon)

    return closed_polygons


def get_cut_polygon(polygon_open, notch, cut_point):
    """
    divides the polygon into two subpieces using diagonal connecting notch and cut point
    :param polygon_open: the input polygon to be divided into two polygons
    :param notch: one point on the diagonal that divides the polygon into two pieces
    :param cut_point: one point on the diagonal that divides the polygon into two pieces
    :returns: list of two polygons that represent the input polygon after being cut by the diagonal
    """
    # get the indices of points that constitute the diagonal
    for i in range(len(polygon_open)):
        if polygon_open[i] == cut_point:    index1 = i
        elif polygon_open[i] == notch:      index2 = i

    # initialize the two polygons that arise from making the cut
    polygon1 = []
    polygon2 = []

    # derive the two resulting polygons
    for i in range(len(polygon_open)):
        if i < min(index1,index2):
            polygon1.append(polygon_open[i])

        elif i == min(index1,index2):
            polygon1.append(polygon_open[i])
            polygon2.append(polygon_open[i])

        elif i > min(index1,index2) and i < max(index1,index2):
            polygon2.append(polygon_open[i])

        elif i == max(index1,index2):
            polygon1.append(polygon_open[i])
            polygon2.append(polygon_open[i])

        elif i > max(index1,index2):
            polygon1.append(polygon_open[i])

    return polygon1,polygon2


def get_valid_resolve(polygons, notch):
    """
    determines if the polygons have succesfully resolved the notch point
    :param polygons: the input list of polygons
    :param notch: the notch point that must be resolved to be a valid decomposition
    """
    # return False if any polygon has less than three points
    min_length = 3
    if min([len(polygon) for polygon in polygons]) < min_length:    return False

    # check that the interior angle at the notch is less than 180 degrees
    else:

        # iterate through the polygons to make sure each is valid
        for polygon in polygons:

            # get the indices of the notch for the polygon 
            for i in range(len(polygon)):
                if polygon[i] == notch:         index_notch = i

            # get the indices immediately before and aftern the notch
            index_previous = index_notch - 1
            if indexapp_notch == len(polygon)-1:   index_next = 0
            else:                               index_next = index_notch + 1

            # determine the interior angle at the notch
            interior_angle = get_interior_angle(polygon[index_previous], polygon[index_notch], polygon[index_next])

            # if interior angle is negative, this indicates greater than 180 degree interior angle
            if interior_angle < 0:              return False

        return True


def get_convex_hull(polygon):
    """
    computes the convex hull of the input polygon
    :param polygon: the input polygon to find the convex hull of
    :returns: a list of points that represent the resulting convex hull
        + considers all points that lie on the convex hull to be part of the convex hull
    """
    # initialize convex hull list
    hull_list = []
    hull_complete = False

    # find bottom-left point and add it to the hull
    left_point_index = np.argmin([x for (x,y) in polygon])
    left_point = tuple(polygon[left_point_index])
    for point in polygon:
        if ((point[0] == left_point[0]) and (point[1] < left_point[1])):    left_point = point
    hull_list.append(left_point)

    # keep track of the previous point and current heading
    previous_point = (left_point[0], left_point[1])
    heading = 0

    # iterate until the hull is complete
    while not(hull_complete):

        # constants
        current_adjustment = 360
        current_distance = 0

        # iterate through points to determine next to add to hull
        for next_point in polygon:
            if not(next_point == previous_point):

                # get the required adjustment angle to reach the next point
                new_adjustment = get_adjustment_angle(heading, previous_point, next_point)

                # keep track of best candidate hull point seen so far
                if ((new_adjustment == current_adjustment and 
                     get_distance(previous_point, next_point) < current_distance) or 
                    (new_adjustment < current_adjustment)):
                    current_adjustment = new_adjustment
                    current_point = next_point
                    current_distance = get_distance(previous_point, next_point)

        # assess the hull complete condition 
        if current_point == left_point:     hull_complete = True

        # add new point to the hull if not complete
        else:
            hull_list.append(current_point)
            heading = heading + current_adjustment
            previous_point = current_point

    return hull_list


def get_adjustment_angle(heading, previous_point, next_point):
    """
    calculates the adjustment angle needed to traverse between two points 
    :param heading: the current heading while approaching the turn point
    :param previous_point: the previous point in the traversal
    :param next_point: the next point in the reversal
    :returns: the adjustment angle needed to transition between points
    """
    # constants
    deg_in_circle = 360
    left = 270
    right = 90

    # calculate difference between points for each dimension
    delta_x = next_point[0] - previous_point[0]
    delta_y = next_point[1] - previous_point[1]

    # calculate the relative angle between points
    relative_angle = np.arctan2(delta_y, delta_x) * 180/np.pi

    # measure relative angle clockwise from East
    if relative_angle < 0:      relative_angle *= -1
    else:                       relative_angle = deg_in_circle - relative_angle

    # measure relative angle clockwise from North
    if relative_angle < left:   relative_angle += right
    else:                       relative_angle -= left

    # measure adjustment angle in range [0,360] degrees
    adjustment_angle = relative_angle - heading
    adjustment_angle = adjustment_angle % deg_in_circle

    return adjustment_angle


def get_interior_angle(prev_point, middle_point, next_point):
    """
    calculates the interior angle at a specified point in the polygon
    :param prev_point: the point immediately before the point of interest 
    :param middle_point: the point of interest
    :param next_point: the point immediately after the point of interest
    :returns: the interior angle at the point of interest
    """
    # constants
    deg_in_circle = 360
    left = 270
    right = 90
    reverse = 180

    # get heading from previous point to notch point
    heading = np.arctan2(middle_point[1]-prev_point[1], middle_point[0]-prev_point[0])*180/np.pi

    # force heading to be measured as clockwise from North
    if heading < 0:     heading *= -1
    else:               heading = deg_in_circle - heading

    # force heading to be measured as clockwise from North
    if heading < left:  heading += right
    else:               heading -= left

    # calculate the interior angle
    interior_angle = reverse - get_adjustment_angle(heading, middle_point, next_point)
    return interior_angle


def get_reformatted_polygon(polygon, length):
    """
    reorders the polygon coordinates to match the ordering of the convex hull 
        and performs interpolation between points if necessary 
    :param polygon: the input polygon to be reformatted
    :param length: the reference length associated with the input polygon
    """
    # convert the points in the polygon to tuple
    #   + consider open polygons to avoid duplicate points
    reordered_polygons = [tuple(point) for point in polygon[:-1]]

    # find bottom-left point to reorder 
    left_point_index = np.argmin([x for (x,y) in reordered_polygons])
    for i in range(len(reordered_polygons)):
        if ((reordered_polygons[i][0] == reordered_polygons[left_point_index][0]) and 
            (reordered_polygons[i][1] <  reordered_polygons[left_point_index][1])):    
            left_point_index = i

    # compute the reordered polygon
    reordered_polygons = reordered_polygons[left_point_index:] + reordered_polygons[:left_point_index]
    
    # close the polygon
    reordered_polygons.append(reordered_polygons[0])

    return reordered_polygons


def get_distance(point1, point2):
    """
    calculates the distance between two points
    :param point1: the first reference point
    :param point2: the second reference point
    :returns: the euclidian distance between points
    """
    return ((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)**0.5


def get_decomposition(input_polygon, tau, length):
    """
    decomposes the input polygon into smaller convex sub-polygons
    :param input_polygon: the input polygon that is to be decomposed
        :requires: first and last point in the polygon are the same point
        :requires: the polygon is of length 4 or greater
        :requires: the points in the polygon are given in counterclockwise order
    :param tau: the non-dimensional tolerance parameter that allows for some level of concavity
        :requires: tau is within the range [0,1]
    :param length: the length scale of the polygon
    :returns: a list of CCW polygons that form the decomposition of the input polygon
    """
    input_polygon_CW = input_polygon[::-1]
    output_polygons_CW = ACD(input_polygon_CW, tau, length)
    output_polygons_CCW = [polygon[::-1] for polygon in output_polygons_CW]
    return output_polygons_CCW


if __name__ == '__main__':
    """
    runs the ACD algorithm with user inputs from command line
        + throws errors if user inputs are invalid
        + sets defualt tau parameter if necessary 
        + sets reference length parameter for each polygon
        + prints out helpful statistics for the user
    """
    # parse the user input
    try:    file_name = sys.argv[1]
    except  IndexError:
        print('>> ERROR: file name argument is required')
        quit()

    # load the file
    try:    data = pygeoj.load(file_name)
    except  FileNotFoundError:
        print('>> ERROR: file not found or incorrect file format')
        quit()

    # parse tau parameter and set default value if necessary
    tau_default = 0.05
    try:    tau = float(sys.argv[2])
    except  IndexError:
        tau = tau_default
    except ValueError:
        tau = tau_default
    file_output = file_name[:-8] + '_convex' + file_name[-8:]

    # initialize new geojson file
    new_file = pygeoj.new()
    print('>> running ACD on:       %s' % file_name)
    print('>> output saved to:      %s' % file_output)
    print('>> concavity tolerance:  %s' % str(tau))

    # iterate through the features of the file and apply the algorithm
    for feature in data:
        
        # iterate through the polygons of each feature
        for polygon in feature.geometry.coordinates:

            # extract a reference length for this polygon 
            bounding_box = feature.geometry.bbox
            length = min(abs(bounding_box[2]-bounding_box[0]),
                         abs(bounding_box[3]-bounding_box[1]))

            # add a new feature for each decomposed polygon
            convex_polygons = get_decomposition(polygon, tau, length)

            # add each new polygon to the new geojson file
            for convex_polygon in convex_polygons:
                new_file.add_feature(geometry={'type': 'Polygon', 'coordinates': [convex_polygon]})

    # save the new geojson file
    new_file.save(file_output) 
    print('>> input # of points:    %d' % sum([len(feature.geometry.coordinates[0]) for feature in data]))
    print('>> input # of polygons:  %d' % len(data))
    print('>> output # of polygons: %d' % len(new_file))
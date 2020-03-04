import math
import numpy as np
from scipy import interpolate

__author__ = 'Allison Parrish'
__email__ = 'allison@decontextualize.com'
__version__ = '0.0.3'

class Bezier:
    """A class to represent Bezier curves.

    :param start: [x, y] of curve start
    :param cp1: [x, y] of first control point
    :param cp2: [x, y] of second control point
    :param end: [x, y] of curve end
    """
    def __init__(self, start, cp1, cp2, end):
        self.start = start
        self.cp1 = cp1
        self.cp2 = cp2
        self.end = end

    def point(self, t):
        """Get coordinates of point at the given position on the curve.

        :param t: position on the curve (should be between 0 and 1)
        :returns: [x, y] of coordinates at position
        """
        adjusted = 1 - t
        x = (
            pow(adjusted, 3) * self.start[0] +
            3 * pow(adjusted, 2) * t * self.cp1[0] +
            3 * adjusted * pow(t, 2) * self.cp2[0] +
            pow(t, 3) * self.end[0]
        )
        y = (
            pow(adjusted, 3) * self.start[1] +
            3 * pow(adjusted, 2) * t * self.cp1[1] +
            3 * adjusted * pow(t, 2) * self.cp2[1] +
            pow(t, 3) * self.end[1]
        )
        return [x, y]

    def tangent(self, t):
        """Get vector of tangent at the given position on the curve.

        :param t: position on the curve (should be between 0 and 1)
        :returns: [x, y] vector of tangent at position
        """
        adjusted = 1 - t
        x = (
            3 * self.end[0] * pow(t, 2) -
            3 * self.cp2[0] * pow(t, 2) +
            6 * self.cp2[0] * adjusted * t -
            6 * self.cp1[0] * adjusted * t +
            3 * self.cp1[0] * pow(adjusted, 2) -
            3 * self.start[0] * pow(adjusted, 2)
        )
        y = (
            3 * self.end[1] * pow(t, 2) -
            3 * self.cp2[1] * pow(t, 2) +
            6 * self.cp2[1] * adjusted * t -
            6 * self.cp1[1] * adjusted * t +
            3 * self.cp1[1] * pow(adjusted, 2) -
            3 * self.start[1] * pow(adjusted, 2)
        )
        return [x, y]

    def flat_commands(self, moveto=True):
        """Get flat commands to draw the curve.

        Convenience function to return path commands for use with the flat
        library.

        :param moveto: if True (default), include "moveto" to curve start
        :returns: list of flat commands
        """

        from flat import command
        cmds = []
        if moveto:
            cmds.append(command.moveto(self.start[0], self.start[1]))
        cmds.append(command.curveto(self.cp1[0],
                                    self.cp1[1],
                                    self.cp2[0],
                                    self.cp2[1],
                                    self.end[0],
                                    self.end[1]))
        return cmds

    def offsets(self, distances):
        """Get list of [x,y] points offset from curve normals.

        Calculates and returns a list of points offset from the normals of the
        curve. The second parameter should be a list with offset distances; the
        length of this list will determine the "resolution" at which the curve
        is sampled.

        :param distances: list of distances
        :returns: list of [x,y] coordinates
        """
        segs = len(distances) - 1
        tangents = [self.tangent(i/segs) for i in range(segs+1)]
        pts = [self.point(i/segs) for i in range(segs+1)]
        offset_pts = []
        for i, (x, y) in enumerate(tangents):
            theta = math.atan2(y, x) + (math.pi * 0.5)
            newx_in = pts[i][0] + math.cos(theta) * distances[i]
            newy_in = pts[i][1] + math.sin(theta) * distances[i]
            offset_pts.append([newx_in, newy_in])
        return offset_pts

    def offset_polygon(self, thicknesses):
        """Get a polygon (list of [x, y] points) of points offset from normals.

        Calculates and returns a list of points offset from the normals of the
        curve, forming a polygon "around" the curve. The "thickness" of the
        polygon at evenly sampled points along the curve are given as the
        second parameter.  (The length of this list determines the resolution
        of curve sampling.)

        :param thicknesses: list of "thickness" (distance from curve)
        :returns: list of [x, y] coordinates
        """
        inner = self.offsets([i*-0.5 for i in thicknesses])
        outer = list(reversed(self.offsets([i*0.5 for i in thicknesses])))
        return inner + outer

def beziers_from_catmull(vertices, tightness=0.0):
    """Smooth a list of points into Bezier curves using Catmull-Rom.

    Applies the Catmull-Rom splines algorithm to get a list of Bezier curves
    passing through the given vertices. The tightness parameter controls the
    "tightness" of the curves (1.0 results in straight lines, 0.0 is the
    default).

    :param vertices: list of [x, y] points (polyline) to smooth
    :param tightness: "tightness" of resulting curves (default 0.0)
    :returns: list of Bezier objects
    """
    if len(vertices) < 3:
        raise ValueError("at least three vertices needed")
    beziers = []
    s = 1 - tightness
    for i in range(1, len(vertices) - 2):
        v = vertices[i]
        b = [[]] * 4
        b[0] = [v[0], v[1]]
        b[1] = [
          v[0] + (s * vertices[i + 1][0] - s * vertices[i - 1][0]) / 6,
          v[1] + (s * vertices[i + 1][1] - s * vertices[i - 1][1]) / 6
        ]
        b[2] = [
          vertices[i + 1][0] +
            (s * vertices[i][0] - s * vertices[i + 2][0]) / 6,
          vertices[i + 1][1] + (s * vertices[i][1] - s * vertices[i + 2][1]) / 6
        ]
        b[3] = [vertices[i + 1][0], vertices[i + 1][1]]
        beziers.append(Bezier(b[0], b[1], b[2], b[3]))
    return beziers

def beziers_flat_commands(beziers):
    """Get a list of flat commands for the list of beziers.

    :param beziers: list of Bezier objects
    :returns: list of flat commands
    """
    commands = []
    for bez in beziers:
        commands.extend(bez.flat_commands())
    return commands

def smooth_point_path(vertices, tightness=0.0):
    """List of flat commands for Catmull Bezier path from points.

    Convenience function to return the flat path commands for a Catmull-Rom
    curve going through the specified points.

    :param vertices: list of [x, y] points
    :param tightness: tightness of curve
    :returns: list of flat commands
    """
    return beziers_flat_commands(beziers_from_catmull(vertices, tightness))

def resample(src, needed_len, kind):
    x = np.linspace(0, len(src), len(src))
    y = np.array(src, dtype=np.float32)
    newx = np.linspace(0, x.shape[0], needed_len)
    terp = interpolate.interp1d(x, y, kind=kind)
    newy = terp(newx)
    return newy

def beziers_tangent_offsets(beziers, distances, samples_per=6, interp='linear'):
    """Get a list of [x,y] points offset from bezier curves in list.

    Calculates and returns a list of points offset from the normals of all
    curves in the given list at the given distances. Each curve is sampled the
    specified number of times, evenly spaced along the curve. If the length of
    the distances list is not the product of the number of beziers and the
    number of samples per curve, it will be resampled with scipy's interp1d
    function using the specified interpolation (default "linear"; also try
    "nearest", "quadratic", "cubic"). 

    :param beziers: list of Bezier objects
    :param distances: list of offset distances
    :param samples_per: number of evenly-spaced samples per curve
    :param interp: scipy interpolation kind
    :returns: list of [x, y] points (polyline)
    """
    needed_pts = len(beziers) * samples_per
    if len(distances) != needed_pts:
        distances = resample(distances, needed_pts, interp)
    offset_pts = beziers[0].offsets(distances[:samples_per])
    for i, bez in enumerate(beziers[1:]):
        this_offset = bez.offsets(distances[(i+1)*samples_per:(i+2)*samples_per])
        offset_pts.extend(this_offset)
    return offset_pts

def beziers_tangent_offset_polygon(beziers, thicknesses, samples_per=6,
        interp='linear'):
    """Get a polygon (list of [x, y] points) offset from normals of curve list.

    Returns a polygon (list of [x, y] points) from normals of the curve, offset 
    with symmetrical distances (given as a list in the second parameter).
    The "thicknesses" list will be interpolated to the needed length using
    the specified interpolation method.

    :param beziers: list of Bezier objects
    :param thicknesses: list of offset distances
    :param samples_per: number of evenly-spaced samples per curve
    :param interp: scipy interpolation kind
    :returns: list of [x, y] coordinates
    """
    inner = beziers_tangent_offsets(
            beziers, [i*-0.5 for i in thicknesses], samples_per, interp)
    outer = beziers_tangent_offsets(
            beziers, [i*0.5 for i in thicknesses], samples_per, interp)
    return inner + list(reversed(outer))

def fancy_curve(vertices, thicknesses, tightness=0.0, samples_per=6,
        interp='linear'):
    """Get a polygon surrounding a Catmull-Rom curve produced from vertices.

    Takes a list of [x,y] points and draws a polygon around a smooth curve
    through those points, with the polygon's thickness determined by values
    in the second parameter. (The first item of the list is the thickness of
    the polygon at the start of the curve; the last item of the list is the
    thickness of the polygon at the end of the curve.) The "thicknesses" list
    will be interpolated to the needed length using the specified interpolation
    method.

    :param vertices: list of [x, y] points
    :param thicknesses: list of offset distances
    :param tightness: tightness of Catmull-Rom curve
    :param interp: scipy interpolation kind
    :returns: list of [x, y] points (polygon)
    """
    beziers = beziers_from_catmull(vertices, tightness)
    polygon = beziers_tangent_offset_polygon(beziers, thicknesses, samples_per,
            interp)
    return polygon

def flatten(t):
    "Convenience function to flatten lists of pts to format required for flat"
    from itertools import chain
    return list(chain(*t))

if __name__ == '__main__':
    from random import randrange
    from flat import document, rgba, shape
    width, height = (500, 500)
    pts = [[100,100], [100,100], [100,400], [400,400], [400,100], [400,100]]
    bez = beziers_from_catmull(pts)
    poly = fancy_curve(pts, [1, 10, 50, 25, 20, 15, 10], samples_per=24,
            interp='cubic')
    bg_shape = shape().nostroke().fill(rgba(255, 255, 255, 255))
    pts_shape = shape().stroke(rgba(255, 0, 0, 240)).width(2)
    catmull_shape = shape().stroke(rgba(0, 255, 0, 240)).width(2)
    poly_shape = shape().nostroke().fill(rgba(0, 0, 0, 220)).width(4)
    doc = document(width, height, 'mm')
    page = doc.addpage()
    page.place(bg_shape.rectangle(0, 0, width, height))
    page.place(poly_shape.polygon(flatten(poly)))
    page.place(catmull_shape.path(smooth_point_path(pts)))
    page.place(pts_shape.polyline(flatten(pts)))
    page.image(kind='rgba').png('test.png')


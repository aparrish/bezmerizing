import math
import numpy as np
from copy import copy
from scipy import interpolate
from flat import command

__author__ = 'Allison Parrish'
__email__ = 'allison@decontextualize.com'
__version__ = '0.1.0'

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

    def to_path(self, moveto=True):
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
        return Path(cmds)

    def offsets(self, distances):
        """Get a Polyline of [x,y] points offset from curve normals.

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
        return Polyline(offset_pts)

    def offset_polygon(self, thicknesses):
        """Get a Polyline of points offset from normals on both sides.

        Calculates and returns a list of points offset from the normals of the
        curve, forming a polygon "around" the curve. The "thickness" of the
        polygon at evenly sampled points along the curve are given as the
        second parameter.  (The length of this list determines the resolution
        of curve sampling.)

        :param thicknesses: list of "thickness" (distance from curve)
        :returns: Polyline
        """
        inner = self.offsets([i*-0.5 for i in thicknesses])
        outer = self.offsets([i*0.5 for i in thicknesses]).reverse()
        return inner + outer


class Spline:
    """A sequence of Bezier curves."""

    def __init__(self, beziers):
        self.beziers = beziers

    def to_path(self):
        """Get a list of flat commands for this spline's bezier curves.

        :returns: Path object
        """
        commands = []
        for bez in self.beziers:
            commands.extend(bez.to_path())
        return Path(commands)

    def tangent_offsets(self, distances, samples_per=6, interp='linear'):
        """Get a Polyline of [x,y] points offset from bezier curves in spline.

        Calculates and returns a list of points offset from the normals of all
        curves in the given list at the given distances. Each curve is sampled
        the specified number of times, evenly spaced along the curve. If the
        length of the distances list is not the product of the number of
        beziers and the number of samples per curve, it will be resampled with
        scipy's interp1d function using the specified interpolation (default
        "linear"; also try "nearest", "quadratic", "cubic"). 

        :param distances: list of offset distances
        :param samples_per: number of evenly-spaced samples per curve
        :param interp: scipy interpolation kind
        :returns: list of [x, y] points (polyline)
        """
        needed_pts = len(self.beziers) * samples_per
        if len(distances) != needed_pts:
            distances = resample(distances, needed_pts, interp)
        offset_pts = self.beziers[0].offsets(distances[:samples_per])
        for i, bez in enumerate(self.beziers[1:]):
            this_offset = bez.offsets(
                    distances[(i+1)*samples_per:(i+2)*samples_per])
            offset_pts += this_offset
        return offset_pts

    def tangent_offset_polygon(self, thicknesses, samples_per=6,
        interp='linear'):
        """Get a polygon (list of [x, y] points) offset from normals of spline.

        Calculates polygon with symmetrical distances (given as a list in the
        second parameter).  The "thicknesses" list will be interpolated to the
        needed length using the specified interpolation method.

        :param thicknesses: list of offset distances
        :param samples_per: number of evenly-spaced samples per curve
        :param interp: scipy interpolation kind
        :returns: list of [x, y] coordinates
        """
        inner = self.tangent_offsets(
                [i*-0.5 for i in thicknesses], samples_per, interp)
        outer = self.tangent_offsets(
                [i*0.5 for i in thicknesses], samples_per, interp)
        return inner + outer.reverse()


class Path:
    """Represents a path as a series of Flat commands."""

    def __init__(self, commands):
        self.commands = commands

    @classmethod
    def frompoints(cls, points, close=False):
        "creates a Path from a list of 2d lists, e.g. [[0,1],[3,4],[5,6]...]"
        commands = [command.moveto(*points[0])]
        for pt in points[1:]:
            commands.append(command.lineto(*pt))
        if close:
            commands.append(command.closepath)
        return cls(commands)

    def scale(self, x):
        "scales the path by the same factor along both axes"
        return self.scalexy(x, x)

    def scalexy(self, x, y):
        "scales the path by the given values on the x and y axes"
        return Path(
            [copy(cmd).transform(x, 0, 0, y, 0, 0) for cmd in self.commands])

    def translate(self, x, y):
        "translates (moves) the path by x and y units in the respective axes"
        return Path(
            [copy(cmd).transform(1, 0, 0, 1, x, y) for cmd in self.commands])

    def rotate(self, theta):
        "rotates the path around the origin with the given angle (in radians)"
        return Path(
                [copy(c).transform(np.cos(theta), -np.sin(theta),
                    np.sin(theta), np.cos(theta), 0, 0) for c in self.commands])

    def __add__(self, other):
        return Path(self.commands + other.commands)

    def __iter__(self):
        return iter(self.commands)

    def __repr__(self):
        return "Path([{0}])".format(
                ", ".join([command_repr(cmd) for cmd in self.commands]))


def command_repr(cmd):
    "Returns a reasonable string representation of a Flat command."
    if isinstance(cmd, command.closepath.__class__):
        return "closepath"
    class_name = type(cmd).__name__
    s = ", ".join(["{k}={v:0.4f}".format(k=slot, v=getattr(cmd, slot))
        for slot in cmd.__slots__])
    return class_name + "(" + s + ")"


class Polyline:

    """Represents a polyline (or polygon) as a 2d array of points."""

    def __init__(self, vertices):
        self.vertices = np.array(vertices)
        if self.vertices.shape[1] != 2:
            raise ValueError("2d data only, sorry")

    def scale(self, x):
        return self.scalexy(x, x)

    def scalexy(self, x, y):
        return Polyline(self.vertices * np.array([x, y]))

    def translate(self, x, y):
        return Polyline(self.vertices + np.array([x, y]))

    def rotate(self, theta):
        return Polyline(self.vertices.dot(
                np.array([[np.cos(theta), -np.sin(theta)],
                          [np.sin(theta), np.cos(theta)]])))

    def reverse(self):
        return Polyline(np.flipud(self.vertices))

    def augment(self):
        """Returns this Polyline with its first and last items repeated

        This is helpful as a brute-force method for calculating a spline from
        pre-existing data where you want the spline to run through the first
        and last points, which are otherwise "eaten" by the algorithm.
        """
        poly = Polyline(np.concatenate([
            [self.vertices[0]],
            self.vertices,
            [self.vertices[-1]]]))
        return poly

    def catmull_spline(self, tightness=0.0):
        """Smooth this polyline into a Bezier spline using Catmull-Rom.

        Applies the Catmull-Rom splines algorithm to get a list of Bezier curves
        passing through the given vertices. The tightness parameter controls the
        "tightness" of the curves (1.0 results in straight lines, 0.0 is the
        default).

        :param vertices: list of [x, y] points (polyline) to smooth
        :param tightness: "tightness" of resulting curves (default 0.0)
        :returns: list of Bezier objects
        """
        if len(self.vertices) < 3:
            raise ValueError("at least three vertices needed")
        beziers = []
        s = 1 - tightness
        for i in range(1, len(self.vertices) - 2):
            v = self.vertices[i]
            b = [[]] * 4
            b[0] = [v[0], v[1]]
            b[1] = [
                v[0]+(s*self.vertices[i+1][0]-s*self.vertices[i-1][0]) / 6,
                v[1]+(s*self.vertices[i+1][1]-s*self.vertices[i-1][1]) / 6
            ]
            b[2] = [
                self.vertices[i + 1][0] +
                    (s * self.vertices[i][0] - s * self.vertices[i+2][0]) / 6,
                self.vertices[i + 1][1] + 
                    (s * self.vertices[i][1] - s * self.vertices[i+2][1]) / 6
            ]
            b[3] = [self.vertices[i + 1][0], self.vertices[i + 1][1]]
            beziers.append(Bezier(b[0], b[1], b[2], b[3]))
        return Spline(beziers)

    def smooth_path(self, tightness=0.0):
        """Path object for Catmull Bezier spline from this polyline's points.

        Convenience function to return the flat path commands for a Catmull-Rom
        curve going through the specified points.

        :param tightness: tightness of curve
        :returns: list of flat commands
        """
        return Path(self.catmull_spline(tightness).to_path())

    def fancy_curve(self, thicknesses, tightness=0.0, samples_per=6,
            interp='linear'):
        """Get a polygon surrounding a Catmull-Rom curve produced from vertices.

        Returns the polygon formed from a Catmull-Rom spline through this
        polyline's points, with the polygon's thickness determined by values
        in the firstparameter. (The first item of the list is the thickness of
        the polygon at the start of the curve; the last item of the list is the
        thickness of the polygon at the end of the curve.) The "thicknesses"
        list will be interpolated to the needed length using the specified
        interpolation method.

        :param thicknesses: list of offset distances
        :param tightness: tightness of Catmull-Rom curve
        :param interp: scipy interpolation kind
        :returns: list of [x, y] points (polygon)
        """
        bez_spline = self.catmull_spline(tightness)
        polygon = bez_spline.tangent_offset_polygon(thicknesses, samples_per,
                interp)
        return polygon

    def __iter__(self):
        return iter(self.vertices.flatten())

    def __getitem__(self, i):
        return self.vertices[int(i/2),i%2]

    def __len__(self):
        return self.vertices.shape[0] * self.vertices.shape[1]

    def __add__(self, other):
        return Polyline(np.concatenate([self.vertices, other.vertices]))
    
    def __repr__(self):
        return "Polyline([{0}])".format(
                ', '.join(["[%0.4f, %0.4f]" % tuple(item)
                    for item in self.vertices.tolist()]))


def resample(src, needed_len, kind):
    x = np.linspace(0, len(src), len(src))
    y = np.array(src, dtype=np.float32)
    newx = np.linspace(0, x.shape[0], needed_len)
    terp = interpolate.interp1d(x, y, kind=kind)
    newy = terp(newx)
    return newy


def flatten(t):
    "Convenience function to flatten lists of pts to format required for flat"
    from itertools import chain
    return list(chain(*t))


if __name__ == '__main__':
    from random import randrange
    from flat import document, rgba, shape
    width, height = (500, 500)
    pts = [[100,100], [100,100], [100,400], [400,400], [400,100], [400,100]]
    pts_poly = Polyline(pts)
    bez_spline = pts_poly.catmull_spline()
    poly = pts_poly.fancy_curve([1, 10, 50, 25, 20, 15, 10], samples_per=24,
            interp='cubic')
    bg_shape = shape().nostroke().fill(rgba(255, 255, 255, 255))
    pts_shape = shape().stroke(rgba(255, 0, 0, 240)).width(2)
    catmull_shape = shape().stroke(rgba(0, 255, 0, 240)).width(2)
    poly_shape = shape().stroke(rgba(64, 64, 64, 255)).fill(rgba(0, 0, 0, 220)).width(4)
    doc = document(width, height, 'mm')
    page = doc.addpage()
    page.place(bg_shape.rectangle(0, 0, width, height))
    page.place(poly_shape.polygon(poly))
    page.place(catmull_shape.path(pts_poly.smooth_path()))
    page.place(pts_shape.polyline(pts_poly))
    page.image(kind='rgba').png('test.png')


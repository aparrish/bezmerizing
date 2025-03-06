import numpy as np
from copy import copy

from bezmerizing.offsetting import (calculate_offsets,
                                    calculate_spline_offsets,
                                    trim_cusps)

__author__ = 'Allison Parrish'
__email__ = 'allison@decontextualize.com'
__version__ = '0.2.0'

class BaseBezier:
    """Base class for Bezier curves."""

    def __init__(self):
        raise NotImplemented

    def to_svg_path(self):
        raise NotImplemented

    def point(self, t):
        raise NotImplemented

    def split(self, t):
        raise NotImplemented

    def control_point_distance(self):
        raise NotImplemented

    def deriv(self):
        raise NotImplemented

    def slice(self, t1, t2):
        "returns a new Bezier, sliced between t1 and t2 from this Bezier"
        first, second = self.split(t1)
        # pretty sure this is right lmao
        slice_, leftover = second.split((t2-t1)/(1-t1))
        return slice_

    def kappa(self, t):
        "curvature of the Bezier at time t"
        d = self.deriv().point(t)
        p = self.deriv().deriv().point(t)
        return ((d[0]*p[1]) - (p[0]*d[1])) / pow(d[0]**2 + d[1]**2, 3/2)

    def tangent(self, t):
        """Get vector of tangent at the given position on the curve.

        :param t: position on the curve (should be between 0 and 1)
        :returns: [x, y] vector of tangent at position
        """
        return self.deriv().point(t)

    def tangent_unit(self, t):
        "unit vector from tangent at time t"
        vec = np.array(self.tangent(t))
        hat = vec / np.linalg.norm(vec)
        return hat

    def normal(self, t):
        "normal vector at time t"
        vec = np.array(self.tangent(t))
        return np.array([vec[1], vec[0] * -1])

    def normal_unit(self, t):
        "unit vector from normal at time t"
        vec = self.normal(t)
        return vec / np.linalg.norm(vec)

    def flatness(self):
        """a simple measurement of "flatness"

        The ratio of the curve's start/end line segment to the max distance
        of the control point to that segment."""

        start_end_length = np.linalg.norm(
            np.array([self.start[0], self.start[1]]) -
            np.array([self.end[0], self.end[1]]))
        return self.control_point_distance() / start_end_length

    def offsets(self, distances, trim=True, subdivide=True, threshold=0.1,
                depth=0, max_depth=5):
        """Get a Polyline of [x,y] points offset from curve normals.

        Calculates and returns a list of points offset from the normals of the
        curve. The second parameter should be a list with offset distances; the
        length of this list will determine the "resolution" at which the curve
        is sampled.

        :param distances: list of distances
        :param trim: if True, attempts to trim cusps
        :param subdivide: if True, subdivides insufficiently "flat" segments
        :param threshold: threshold for subdivision (between 0 and 1)
        :param depth: current recursion depth (set to 0)
        :param max_depth: recursion limit
        :returns: Polyline of [x,y] coordinates
        """

        offset_pts, sign_change_idx = calculate_offsets(
                curve=self, distances=distances, subdivide=subdivide,
                threshold=threshold, depth=depth, max_depth=max_depth)

        if trim and len(sign_change_idx) > 1:
            return trim_cusps(Polyline(offset_pts), sign_change_idx)
        else:
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

    def to_polyline(self, samples_per=6):
        "returns a polyline generated from samples of the curve"
        pts = []
        for i in np.linspace(0, 1, samples_per):
            pts.append(self.point(i))
        return Polyline(pts)


class Bezier(BaseBezier):

    """A class to represent (cubic) Bezier curves.

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

    def to_svg_path(self, close=False):
        return " ".join([
            f"M {self.start[0]},{self.start[1]}",
            f"C {self.cp1[0]},{self.cp1[1]}",
            f"{self.cp2[0]},{self.cp2[1]}",
            f"{self.end[0]},{self.end[1]}"]) + (" Z" if close else "")

    def split(self, t):
        "splits the curve at time t, returning the resulting two curves"

        p4 = [lerp(self.start[0], self.cp1[0], t),
              lerp(self.start[1], self.cp1[1], t)]
        p5 = [lerp(self.cp1[0], self.cp2[0], t),
              lerp(self.cp1[1], self.cp2[1], t)]
        p6 = [lerp(self.cp2[0], self.end[0], t),
              lerp(self.cp2[1], self.end[1], t)]
        p7 = [lerp(p4[0], p5[0], t),
              lerp(p4[1], p5[1], t)]
        p8 = [lerp(p5[0], p6[0], t),
              lerp(p5[1], p6[1], t)]
        p9 = [lerp(p7[0], p8[0], t),
              lerp(p7[1], p8[1], t)]

        return (Bezier(start=self.start, cp1=p4, cp2=p7, end=p9),
                Bezier(start=p9, cp1=p8, cp2=p6, end=self.end))

    def control_point_distance(self):
        "max distance of the control points to start/end"
        c1 = np.array(self.cp1)
        c2 = np.array(self.cp2)
        start = np.array(self.start)
        end = np.array(self.end)
        c1_dist_sq = ((np.linalg.norm(c1-start)**2) - 
                      (np.dot(c1-start, end-start)**2) / np.dot(end-start,
                                                               end-start))
        c2_dist_sq = ((np.linalg.norm(c2-start)**2) -
                      (np.dot(c2-start, end-start)**2) / np.dot(end-start,
                                                                end-start))
        return np.sqrt(max(c1_dist_sq, c2_dist_sq))

    def deriv(self):
        "returns derivative of the curve"
        return QuadraticBezier(
                start=[
                    (self.cp1[0] - self.start[0]) * 3,
                    (self.cp1[1] - self.start[1]) * 3
                ],
                cp1=[
                    (self.cp2[0] - self.cp1[0]) * 3,
                    (self.cp2[1] - self.cp1[1]) * 3
                ],
                end=[
                    (self.end[0] - self.cp2[0]) * 3,
                    (self.end[1] - self.cp2[1]) * 3
                ]
            )

    def inflections(self):
        "inflection points of this curve"

        x0 = self.start[0]
        x1 = self.cp1[0]
        x2 = self.cp2[0]
        x3 = self.end[0]
        y0 = self.start[1]
        y1 = self.cp1[1]
        y2 = self.cp2[1]
        y3 = self.end[1]

        ax = -x0 + 3*x1 - 3*x2 + x3
        bx = 3*x0 - 6*x1 + 3*x2
        cx = -3*x0 + 3*x1
        dx = x0

        ay = -y0 + 3*y1 - 3*y2 + y3
        by = 3*y0 - 6*y1 + 3*y2
        cy = -3*y0 + 3*y1
        dy = y0

        cusp_t = -0.5 * ((ay*cx - ax*cy) / (ay*bx - ax*by))
        t1 = cusp_t - np.sqrt(
                pow(cusp_t, 2) - (1/3)*((by*cx - bx*cy) / (ay*bx - ax*by)))
        t2 = cusp_t + np.sqrt(
                pow(cusp_t, 2) - (1/3)*((by*cx - bx*cy) / (ay*bx - ax*by)))

        return [t1, t2]

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

    def __repr__(self):
        formatted = {k: "[%0.4f, %0.4f]" % tuple(getattr(self, k)) for k in
            ['start', 'cp1', 'cp2', 'end']}
        return ("Bezier(start={start}, cp1={cp1}, cp2={cp2}, end={end})"
                .format(**formatted))


class QuadraticBezier(BaseBezier):

    """A class to represent quadratic Bezier curves.

    :param start: [x, y] of curve start
    :param cp1: [x, y] of control point
    :param end: [x, y] of curve end
    """

    def __init__(self, start, cp1, end):
        self.start = start
        self.cp1 = cp1
        self.end = end

    def to_svg_path(self, close=False):
        return " ".join([
            f"M {self.start[0]},{self.start[1]}",
            f"Q {self.cp1[0]},{self.cp1[1]}",
            f"{self.end[0]},{self.end[1]}"]) + (" Z" if close else "")

    def control_point_distance(self):
        "max distance of the control points to start/end"
        c1 = np.array(self.cp1)
        start = np.array(self.start)
        end = np.array(self.end)
        c1_dist_sq = ((np.linalg.norm(c1-start)**2) - 
                      (np.dot(c1-start, end-start)**2) / np.dot(end-start,
                                                               end-start))
        return np.sqrt(c1_dist_sq)

    def split(self, t):
        "splits the curve at time t, returning the resulting two curves"
        p4 = [lerp(self.start[0], self.cp1[0], t),
              lerp(self.start[1], self.cp1[1], t)]
        p5 = [lerp(self.cp1[0], self.end[0], t),
              lerp(self.cp1[1], self.end[1], t)]
        mid = [lerp(p4[0], p5[0], t), lerp(p4[1], p5[1], t)]
        return (QuadraticBezier(self.start, p4, mid),
                QuadraticBezier(mid, p5, self.end))

    def deriv(self):
        "derivative of the curve (a line segment)"
        return LineSegment(
            start=[(self.cp1[0] - self.start[0]) * 2,
                   (self.cp1[1] - self.start[1]) * 2],
            end=[(self.end[0] - self.cp1[0]) * 2,
                 (self.end[1] - self.cp1[1]) * 2])

    def point(self, t):
        """Get coordinates of point at the given position on the curve.

        :param t: position on the curve (should be between 0 and 1)
        :returns: [x, y] of coordinates at position
        """
        adjusted = 1 - t
        x = (
            pow(adjusted, 2) * self.start[0] + 
            2 * adjusted * t * self.cp1[0] +
            pow(t, 2) * self.end[0]
        )
        y = (
            pow(adjusted, 2) * self.start[1] + 
            2 * adjusted * t * self.cp1[1] +
            pow(t, 2) * self.end[1]
        )
        return [x, y]

    def to_cubic_bezier(self):
        "elevate curve order to cubic"
        return Bezier(
                start=self.start,
                cp1=[self.start[0]*(1/3) + self.cp1[0]*(2/3),
                     self.start[1]*(1/3) + self.cp1[1]*(2/3)],
                cp2=[self.cp1[0]*(2/3) + self.end[0]*(1/3),
                     self.cp1[1]*(2/3) + self.end[1]*(1/3)],
                end=self.end)

    def __repr__(self):
        formatted = {k: "[%0.4f, %0.4f]" % tuple(getattr(self, k)) for k in
            ['start', 'cp1', 'end']}
        return ("QuadraticBezier(start={start}, cp1={cp1}, end={end})"
                .format(**formatted))

class LineSegment:
    """A class to represent line segments, primarily as the derivative of
    a quadratic Bezier curve."""

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def point(self, t):
        return [self.start[0] + (self.end[0] - self.start[0]) * t,
                self.start[1] + (self.end[1] - self.start[1]) * t]

    def deriv(self):
        # this is a cheat! I don't want to overcomplicate the Bezier
        # kappa method by checking for constants vs. objects with .point()
        # methods, and I don't need a Point class. so we represent
        # the derivative of a LineSegment as a LineSegment whose start and
        # end points are the same.
        pt = [self.end[0] - self.start[0], self.end[1] - self.start[1]]
        return LineSegment(start=pt, end=pt)

    def to_svg_path(self, close=False):
        return (f"M {self.start[0]},{self.start[1]} " 
                + f"L {self.end[0]},{self.end[1]}"
                + (" Z" if close else ""))

    def to_polyline(self):
        return Polyline([self.start, self.end])

    def intersection(self, other):
        "returns intersection of this line segment with another"
        
        x1 = self.start[0]
        y1 = self.start[1]
        x2 = self.end[0]
        y2 = self.end[1]
        
        x3 = other.start[0]
        y3 = other.start[1]
        x4 = other.end[0]
        y4 = other.end[1]

        d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if d == 0:
            return np.nan
        t = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / d
        u = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1-x3)) / d

        if not((0 <= t <= 1) and (0 <= u <= 1)) or d == 0:
            return np.nan
        else:
            return self.point(t)

    def __repr__(self):
        return "LineSegment(start=[{0}, {1}], end=[{2}, {3}])".format(
                self.start[0], self.start[1], self.end[0], self.end[1])


class Spline:
    """A sequence of Bezier curves."""

    def __init__(self, beziers):
        self.beziers = beziers

    def to_polyline(self, samples_per=6):
        "Flatten spline to a single polyline by sampling each component curve"
        return Polyline(
                flatten([s.to_polyline(samples_per) for s in self.beziers]))

    def to_svg_path(self, close=False):
        return " ".join([s.to_svg_path(close) for s in self.beziers])

    def offsets(self, distances, samples_per=6, trim=True, subdivide=True,
                threshold=0.1, depth=0, max_depth=5):
        """Get a Polyline of [x,y] points offset from bezier curves in spline.

        Calculates and returns a list of points offset from the normals of all
        curves in the given list at the given distances. Each curve is sampled
        the specified number of times, evenly spaced along the curve. If the
        length of the distances list is not the product of the number of
        beziers and the number of samples per curve, it will be resampled
        to the needed length.

        :param distances: list of offset distances
        :param samples_per: number of evenly-spaced samples per curve
        :param trim: if True, attempts to trim cusps
        :param subdivide: if True, subdivides insufficiently "flat" segments
        :param threshold: threshold for subdivision (between 0 and 1)
        :param depth: current recursion depth (set to 0)
        :param max_depth: recursion limit
        :returns: Polyline of [x,y] coordinates
        """

        return calculate_spline_offsets(self, distances, samples_per,
                                        trim, subdivide, threshold, depth,
                                        max_depth)

    def offset_polygon(self, thicknesses, samples_per=6):
        """Get a polygon offset from normals of spline.

        Calculates polygon with symmetrical distances (given as a list in the
        second parameter).  The "thicknesses" list will be interpolated to the
        needed length.

        :param thicknesses: list of offset distances
        :param samples_per: number of evenly-spaced samples per curve
        :returns: list of [x, y] coordinates
        """
        inner = self.offsets(
                [i*-0.5 for i in thicknesses], samples_per)
        outer = self.offsets(
                [i*0.5 for i in thicknesses], samples_per)
        return inner + outer.reverse()


class Polyline:

    """Represents a polyline (or polygon) as a 2d array of points."""

    def __init__(self, vertices):
        self.vertices = np.array(vertices)
        if self.vertices.shape[1] != 2:
            raise ValueError("2d data only, sorry")

    def to_svg_path(self, close=False):
        return (
            f"M {self.vertices[0][0]},{self.vertices[0][1]} " 
            + " ".join([f"L {item[0]},{item[1]}" for item in self.vertices])
            + (" Z" if close else ""))

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

        Applies the Catmull-Rom splines algorithm to get a list of Bezier
        curves passing through the given vertices. The tightness parameter
        controls the "tightness" of the curves (1.0 results in straight
        lines, 0.0 is the default).

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

    def fancy_curve(self, thicknesses, tightness=0.0, samples_per=6):
        """Get the polygon of a Catmull-Rom curve produced from vertices.

        Returns the polygon formed from a Catmull-Rom spline through this
        polyline's points, with the polygon's thickness determined by values
        in the firstparameter. (The first item of the list is the thickness of
        the polygon at the start of the curve; the last item of the list is the
        thickness of the polygon at the end of the curve.) The "thicknesses"
        list will be interpolated to the needed length.

        :param thicknesses: list of offset distances
        :param tightness: tightness of Catmull-Rom curve
        :returns: list of [x, y] points (polygon)
        """
        bez_spline = self.catmull_spline(tightness)
        polygon = bez_spline.offset_polygon(thicknesses, samples_per)
        return polygon

    def resample(self, samples_per=6):
        """'Resamples' a polyline

        Returns a new polyline with each constituent line of the original
        replaced with a new polyline having the specified number of points
        along the original line."""
        pts = []
        for i in range(len(self.vertices) - 1):
            if i == len(self.vertices) - 2:
                incl_end = True
            else:
                incl_end = False
            newl = np.linspace(
                self.vertices[i],
                self.vertices[i+1],
                num=samples_per,
                endpoint=incl_end).tolist()
            pts.extend(newl)
        return Polyline(pts)

    def area_signed(self):
        "calculates the signed area of the polygon"
        x = self.vertices[:,0]
        y = self.vertices[:,1]
        area = 0.5 * (np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
        return area

    def area(self):
        "convenience function that gets the absolute value of the area"
        return np.abs(self.area_signed())

    def is_clockwise(self):
        "true if the polyline's points are in clockwise order"
        return bool(self.area_signed() < 0)

    def is_counterclockwise(self):
        "true if the polyline's points are all widdershins"
        return bool(self.area_signed() > 0)

    def intersections(self, other, first_only=False):
        """Calculates the intersections between this polyline and another.

        The return value is a tuple whose first element is a list, where
        each element of the list whose first item is the index of the
        intersecting segment in this polyline, and whose second item is
        the index of the intersecting segment in the other polyline.
        The second element of the return tuple is a list of [x, y] points
        where the intersections occurred, which corresponds index-wise
        with the list in the first element of the tuple.

        :param other: a PolyLine object
        :param first_only: if True, returns after finding one intersection
        :returns: tuple of lists: (indices, coordinates)
        """

        intersection_idx = []
        intersection_pts = []

        # yes this is O(n^2). but short of an rtree I don't know if there's
        # really any low hanging fruit for optimization that is worthwhile?
        for i, (my_pt1, my_pt2) in enumerate(zip(self.vertices[:-1],
                                                 self.vertices[1:])):
            my_seg = LineSegment(my_pt1, my_pt2)
            for j, (their_pt1, their_pt2) in enumerate(
                    zip(other.vertices[:-1], other.vertices[1:])):
                their_seg = LineSegment(their_pt1, their_pt2)
                intersection = my_seg.intersection(their_seg)
                if intersection is not np.nan:
                    intersection_idx.append([i, j])
                    intersection_pts.append(intersection)
                    if first_only: break

        return (intersection_idx, intersection_pts)

    def __getitem__(self, i):
        return self.vertices[i]

    def __len__(self):
        return len(self.vertices)

    def __add__(self, other):
        return Polyline(np.concatenate([self.vertices, other.vertices]))
    
    def __repr__(self):
        return "Polyline([{0}])".format(
                ', '.join(["[%0.4f, %0.4f]" % tuple(item)
                    for item in self.vertices.tolist()]))


class PolylineList:

    def __init__(self, plines):
        self.plines = plines

    def to_svg_path(self, close=False):
        return " ".join([item.to_svg_path(close) for item in self.plines])

    def scale(self, x):
        return self.scalexy(x, x)

    def scalexy(self, x, y):
        return PolylineList([p.scalexy(x, y) for p in self.plines])

    def translate(self, x, y):
        return PolylineList([p.translate(x, y) for p in self.plines])

    def rotate(self, theta):
        return PolylineList([p.rotate(theta) for p in self.plines])

    def map(self, func):
        return PolylineList([func(p) for p in self.plines])

    def __getitem__(self, i):
        return self.plines[i]


def flatten(t):
    "Convenience function to flatten lists of 2d points"
    from itertools import chain
    return list(chain(*t))


def lerp(a, b, t):
    "good ol' lerp"
    return (a * (1 - t)) + (b * t)

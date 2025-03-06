import numpy as np
import bezmerizing
from scipy.interpolate import CubicSpline

def calculate_offsets(curve, distances, subdivide=True, threshold=0.1,
                      depth=0, max_depth=5):
    "calculate offsets of the given Bezier curve at given distances"

    segs = len(distances) - 1
    pts = len(distances)
    offset_pts = []
    sign_change_idx = []
    current_sign = None

    for i in range(segs+1):

        # calculations for this iteration
        t = i / segs
        next_t = (i + 1) / segs
        px, py = curve.point(t)
        nx, ny = curve.normal_unit(t)
        k = curve.kappa(t)
        curve_sign = np.sign(1 + (k * distances[i]))

        # deal with curve sign. if the sign change is happening on the very
        # last sample, we'll leave it "hanging"; downstream filters will handle
        # this as either a cusp between curves (spline) or a line segment to
        # trim back to the first sign change
        if (curve_sign != current_sign and current_sign is not None):
            sign_change_idx.append(len(offset_pts))
        current_sign = curve_sign

        #print(i, sign_change_idx, len(offset_pts))

        # the "typical" case: simply add a point to the output at the
        # unit normal times the specified distance
        offset_pts.append([px + (nx * distances[i]),
                           py + (ny * distances[i])])

        # attempt to subdivide this segment if we're on a convex portion of
        # the curve
        if (subdivide and depth < max_depth and k * distances[i] > 0
            and i < segs):

            slice_ = curve.slice(t, next_t)
            flatness = slice_.flatness()
            if flatness > threshold:

                # split the curve in two
                first, second = slice_.split(0.5)
                # interpolate some new points for each curve
                new_distances = np.linspace(distances[i], distances[i+1], 4)
                # recurse
                first_off, _ = calculate_offsets(curve=first,
                                                distances=new_distances[:2],
                                                subdivide=subdivide,
                                                depth=depth,
                                                max_depth=max_depth)
                second_off, _ = calculate_offsets(curve=second,
                                                 distances=new_distances[2:],
                                                 subdivide=subdivide,
                                                 depth=depth,
                                                 max_depth=max_depth)

                offset_pts.extend(first_off)
                offset_pts.extend(second_off)

    return (offset_pts, sign_change_idx)


def calculate_spline_offsets(spline, distances, samples_per=6, trim=True,
                             subdivide=True, threshold=0.1, depth=0,
                             max_depth=5):
    "calculates a polyline of offsets from the curves in the given spline"

    # calculate array of distances by resampling the given list of distances
    # to the appropriate length
    needed_pts = len(spline.beziers) * samples_per
    if len(distances) != needed_pts:
        distances = resample(distances, needed_pts)

    all_offset_pts = []
    all_sign_change_idx = []

    # flatten each bezier in the spline; add their coords and sign
    # changes to one big array
    for i, bez in enumerate(spline.beziers):
        this_distances = distances[i*samples_per:(i+1)*samples_per]
        this_pts, this_sign_change_idx = calculate_offsets(
                curve=bez, distances=this_distances,
                threshold=threshold, depth=depth, max_depth=max_depth)
        all_offset_pts.append(this_pts)
        all_sign_change_idx.append(this_sign_change_idx)

    flattened = []

    # if we're not going to trim cusps, just throw everything in there
    if not(trim):
        for item in all_offset_pts:
            flattened.extend(item)

    # otherwise, lotsa trimmin' to do
    else:

        # FIXME: this could be way more robust :( someone who knows
        # math please help, my splines are dying

        i = 0
        while i < len(all_offset_pts):

            this_pts = all_offset_pts[i]
            try:
                next_pts = all_offset_pts[i+1]
            except IndexError:
                next_pts = []
            this_chg_idx = all_sign_change_idx[i]
            try:
                next_chg_idx = all_sign_change_idx[i+1]
            except IndexError:
                next_chg_idx = []

            # big if/elif/else block to handle different possibilities
            # related to cusps: inside curves and across curves

            # started with a "half cusp"; snip and move to next pair
            if (i == 0 and len(this_chg_idx) == 1
                    and this_chg_idx[0] < len(this_pts) / 2):
                flattened.extend(this_pts[this_chg_idx[0]:])
                i += 1
                continue

            # this curve has two sign changes: it's a classic
            # curve-internal cusp
            elif (len(this_chg_idx) == 2
                  and this_chg_idx[1] < len(this_pts) - 1):
                trimmed = trim_cusps(bezmerizing.Polyline(this_pts),
                                     this_chg_idx)
                flattened.extend(list(trimmed.vertices))
                i += 1
                continue

            # cusp across segments when there are two sign changes, but
            # the second falls on the last index of the polyline
            elif (len(this_chg_idx) == 2
                  and this_chg_idx[1] == len(this_pts) - 1
                  and len(next_chg_idx) > 0):
                combined = (bezmerizing.Polyline(this_pts)
                            + bezmerizing.Polyline(next_pts))
                trimmed = trim_cusps(
                        combined,
                        [this_chg_idx[1], next_chg_idx[0]+len(this_pts)])
                flattened.extend(list(trimmed.vertices))
                i += 2
                continue

            # cusp across segments
            elif len(this_chg_idx) == 1 and len(next_chg_idx) == 1:
                combined = (bezmerizing.Polyline(this_pts)
                            + bezmerizing.Polyline(next_pts))
                trimmed = trim_cusps(
                        combined,
                        [this_chg_idx[0], next_chg_idx[0]+len(this_pts)])
                flattened.extend(list(trimmed.vertices))
                i += 2
                continue

            # incomplete cusp in this curve, but two cusps in next
            elif len(this_chg_idx) == 1 and len(next_chg_idx) == 2:
                combined = (bezmerizing.Polyline(this_pts)
                            + bezmerizing.Polyline(next_pts))
                trimmed = trim_cusps(
                        combined,
                        [this_chg_idx[0], next_chg_idx[0]+len(this_pts)])
                flattened.extend(list(trimmed.vertices))
                i += 2

            # "half cusp" at the end of the list of points
            elif (i == len(all_offset_pts) - 1
                  and len(this_chg_idx) == 1
                  and this_chg_idx[0] > len(next_pts) / 2):
                flattened.extend(this_pts[:this_chg_idx[0]])
                break

            # done!
            elif i == len(all_offset_pts) - 1:
                flattened.extend(this_pts)
                break

            else:
                flattened.extend(this_pts)
                i += 1

    return bezmerizing.Polyline(flattened)


def trim_cusps(polyline, sign_change_idx):
    "trim cusps of a curve offset polyline, following Elber and Cohen 1991"
    pline1 = bezmerizing.Polyline(polyline.vertices[:sign_change_idx[0]+1])
    pline2 = bezmerizing.Polyline(polyline.vertices[sign_change_idx[1]:])

    intersection = pline1.reverse().intersections(pline2, first_only=True)

    if len(intersection[0]) > 0:
        idxa, idxb = intersection[0][0]
        return (
            bezmerizing.Polyline(pline1.vertices[:sign_change_idx[0]-idxa-1]) +
            bezmerizing.Polyline(intersection[1]) +
            bezmerizing.Polyline(pline2.vertices[idxb+1:]))
    else:
        # if we couldn't find an intersection, the curve is probably doubling
        # back on itself; we'll just join this with a straight line
        return bezmerizing.Polyline(pline1.vertices[:-1]) + pline2

def resample(src, needed_len):
    x = np.linspace(0, 1, len(src))
    y = np.array(src)
    xnew = np.linspace(0, 1, needed_len)
    return CubicSpline(x, y)(xnew)


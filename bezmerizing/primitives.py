from bezmerizing import Polyline
import numpy as np

def arc(x, y, rx, ry, start, end, n=64):
    """Returns a polyline with points along an arc.

    Note that the n parameter includes the final point in the arc, so if
    you want to use this function to produce a regular polygon, pass in
    the number of sides you want minus one.
    
    :param x: the x coordinate of the arc's center
    :param y: the y coordinate of the arc's center
    :param rx: the radius of the arc on the x axis
    :param ry: the radius of the arc on the y axis
    :param start: where the arc should start (in radians)
    :param end: where the arc should end (in radians)
    :param n: the number of points approximating the arc
    """
    pts = []
    for i, theta in enumerate(np.linspace(start, end, n)):
        px = x + np.cos(theta) * rx
        py = y + np.sin(theta) * ry
        pts.append([px, py])
    return Polyline(pts)

def ellipse(x, y, rx, ry, n=64):
    """Returns a polyline approximating an ellipse.
    
    :param x: the x coordinate of the ellipse's center
    :param y: the y coordinate of the ellipse's center
    :param rx: the radius of the arc on the x axis
    :param ry: the radius of the arc on the y axis
    :param n: the number of points approximating the arc
    """
    return arc(x, y, rx, ry, 0, np.pi*2, n=n)

def circle(x, y, rx, n=64):
    """Returns a polyline approximating a circle.
    
    :param x: the x coordinate of the circle's center
    :param y: the y coordinate of the circle's center
    :param rx: the circle's radius
    :param n: the number of points approximating the arc
    """
    return ellipse(x, y, rx, rx, n)

def rect(x, y, w, h):
    """Returns a rectangular polyline.

    :param x: the x coordinate of the rectangle's upper left corner
    :param y: the y coordinate of the rectangle's upper right corner
    :param w: the rectangle's width
    :param h: the rectangle's height
    """
    return Polyline([[x, y], [x+w, y], [x+w, y+h], [x, y+h]])


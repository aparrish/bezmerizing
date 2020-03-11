from bezmerizing import Polyline
import numpy as np

def arc(x, y, rx, ry, start, end, n=64):
    pts = []
    for i, theta in enumerate(np.linspace(start, end, n+1)[:-1]):
        px = x + np.cos(theta) * rx
        py = y + np.sin(theta) * ry
        pts.append([px, py])
    return Polyline(pts)

def ellipse(x, y, rx, ry, n=64):
    return arc(x, y, rx, ry, 0, np.pi*2, n=n)

def circle(x, y, rx, n=64):
    return ellipse(x, y, rx, rx, n)

def rect(x, y, w, h):
    return Polyline([[x, y], [x+w, y], [x+w, y+h], [x, y+h]])


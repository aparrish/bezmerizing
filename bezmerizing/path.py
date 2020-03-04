from flat import command
from copy import copy
from math import sin, cos

class Path:
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
                [copy(cmd).transform(cos(theta), -sin(theta),
                    sin(theta), cos(theta), 0, 0) for cmd in self.commands])

    def __iter__(self):
        return iter(self.commands)

    def __repr__(self):
        return "Path([{0}])".format(
                ", ".join([command_repr(cmd) for cmd in self.commands]))


def command_repr(cmd):
    if isinstance(cmd, command.closepath.__class__):
        return "closepath"
    class_name = type(cmd).__name__
    s = ", ".join(["{k}={v:0.4f}".format(k=slot, v=getattr(cmd, slot))
        for slot in cmd.__slots__])
    return class_name + "(" + s + ")"



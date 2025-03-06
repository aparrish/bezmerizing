# Bezmerizing

By [Allison Parrish](https://www.decontextualize.com/)

*Bezmerizing* is a tiny quirky library with some potentially helpful classes
and functions for working with Bezier curve functions, like:

* Finding points on the curve
* Calculating tangents
* Generating curves from lists of points with Catmull-Rom
* Producing lists of points offset at a certain distance from curve normals

Additionally, the `fancy_curve()` function generates a polygon that traces
"around" a list of Bezier curves, with adjustable thickness along the curve.

## Installation

From this repository:

    pip install https://github.com/aparrish/bezmerizing/archive/main.zip

## Usage

[See the demo notebook](demo.ipynb).

## Change log

### Version 0.2.0

This version includes some breaking changes from previous versions.

* Added more mathy features for the Bézier curve classes (e.g., derivatives,
  normals, inflections, measure of flatness, curve splitting).
* Most objects now have a `.to_svg_path()` method, which returns the [SVG path
  definition](https://developer.mozilla.org/en-US/docs/Web/SVG/Attribute/d)
  necessary to draw the shape in question.
* The `.offsets()` methods now perform some simple normalization techniques on
  offset curves, including cusp trimming and recursive subdivision.
* The `interp` parameter for Polyline curves has been removed. (Interpolation is
  always performed with scipy's CubicSpline class.)
* Added methods to find the area and winding direction (clockwise or
  counterclockwise) of polylines.
* Added method to find intersections of polylines.
* Adopted `pyproject.toml` for all you people that run up-to-date Python
  versions.
* Removed dependency on [flat](https://github.com/xxyxyz/flat/), along with
  explicit integration features (such as the `Path` class).
* Changed license to Anti-Capitalist Software License v1.4.

## Credits and acknowledgements

I've had [Pomax's grand tome](https://pomax.github.io/bezierinfo/) open in a
tab for weeks, and I used [Simon Cozens'
beziers.py](https://github.com/simoncozens/beziers.py) to sanity check some of
my own implementations. Big thanks to [whoever wrote these lecture
notes](https://www.clear.rice.edu/comp360/lectures/old/SubdivisionTextNew.pdf)
and [this Stack Exchange
thread](https://math.stackexchange.com/questions/1952351/bezier-offset-self-intersections).

Much of the original Bezier math was adapted from
[p5.js](https://github.com/processing/p5.js/blob/0.6.1/src/core/curves.js).

## License

See `LICENSE`.


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

    pip install https://github.com/aparrish/bezmerizing/archive/master.zip

The library requires `flat` (for drawing the curves), and `scipy` with `numpy`
(for interpolating the "thickness" parameter of the polygon curves).

## Usage

[See the demo notebook](demo.ipynb).

## Credits

Most of the Bezier math is copy/pasted from
[p5.js](https://github.com/processing/p5.js/blob/0.6.1/src/core/curves.js).

## License

See `LICENSE`.


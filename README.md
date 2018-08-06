# Shapestats
[![Build Status](https://travis-ci.org/ljwolf/shapestats.svg?branch=master)](https://travis-ci.org/ljwolf/shapestats)
[![DOI](https://zenodo.org/badge/143476127.svg)](https://zenodo.org/badge/latestdoi/143476127)

This is a basic library to compute shape statistics about polygons. 
These include:
- isoperimetric quotients (aka the Polsby-Popper measure)
- convex hull area/perimeter measures
- length/width ratios
- reflex angle measures
- maximum contained circle
- minimum bounding circles using Skyum's algorithm:
![minimum bounding circle computation](https://raw.githubusercontent.com/ljwolf/shapestats/master/_img/minbc.gif)

Each of the measures in `compactness` are defined in their docstrings:
```python
help(shapestats.boundary_amplitude)
```

also, you can use the minium bounding circle & minimum contained circle constructors directly on shapely shapes:

```python
shapestats.maxbc.maximum_bounding_circle(polygon)
shapestats.maxbc.minimum_contained_circle(polygon)
```

# usage

```python
import shapestats
import geopandas
df = geopandas.read_file(geopandas.datasets.get_path("nybb"))
df.geometry.apply(shapestats.ipq)
```

# dependencies
`shapely`
`scipy`
`libpysal`

# License
Copyright 2018 Levi John Wolf

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

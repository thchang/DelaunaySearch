# Delaunay Search Black-Box Optimizer

This package contains a Python 3 implementation of an experimental optimization
search method.
The Delaunay search technique is intended for usage with black-box functions
in situations where the user already has a significant database of
pre-evaluated objective-function values.
To refine the best prediction made by these values, a Delaunay triangulation
(unstructured mesh) is implemented over the dataset.
In each iteration:
 - The candidate cells are the simplices incident at the current minima
   (i.e., the umbrella neighborhood).
 - The objective value at the center of each candidate cell is predicted
   based on the average prediction of each of its d+1 neighboring cells
   (via linear extrapolatoin).
 - The candidate cell with the lowest predicted objective value is subdivided
   by way of a function evaluation at its center.

For objective functions with Lipschitz continuous first derivatives, the
linear fit in each Delaunay cell will converge to the average gradient in that
cell almost surely as more data is sampled.
Therefore, the Delaunay search will converge to a local minima.

This is experimental research code, and is not yet suitable for all use
cases.
In particular, if the global minima is not in the convex hull of the provided
data set, then the Delaunay search can never reach it.

## Contents

 - delsearch.py contains the DelaunaySearch subroutine and an optional main
   program which uses a quadratic function and latin hypercube design for
   testing.

## Prerequisites

The DelaunaySearch function is written for Python 3 and depends on several
other python modules that can be installed via pip3.
 - numpy for numerical linear algebra
 - scipy.spatial for their Delaunay qhull wrapper
 - matplotlib.pyplot for visualizing results
 - pyDOE for their lhs (latin hypercube design) function

### User Interface

To use DelaunaySearch, call
``
[minind, x, f] = DelaunaySearch(data, obj_func, budget)
``
where
 - data is a python list containing the input data set (the first d entries
in each row are the coordinates of the points, and the last entry is the
objective value),
 - obj\_func is the objective function,
 - budget (optional) is the function evaluation budget
 - minind is the index of the minimum objective value observed
 - x is a list of design points evaluated (including those provided in data)
 - f is a list of corresponding objective values

### Test cases

To optimize a 2 dimensional quadratic using a latin hypercube design and
visualize the results in matplotlib:

``
python3 delsearch.py
``

## Author

* **Tyler H. Chang** - *Primary author*


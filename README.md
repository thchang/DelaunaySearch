# Delaunay Search Black-Box Optimizer

This repo contains a Python 3 implementation of an experimental Delaunay
mesh refinedment based black-box optimizer.
The Delaunay search technique is intended for usage with expensive black-box
functions in situations where the user already has a significant database of
pre-evaluated cost values.
To refine the current database toward the true minima, a Delaunay triangulation
(unstructured mesh) is implemented over the dataset.
In each iteration:
 - The objective value at the circumcenter of each element in the mesh is
   estimated based on the average linearly extrapolated value predicted by
   its d+1 neighboring cells.
 - The candidate cell with the lowest predicted objective value is subdivided
   by way of a function evaluation at its circumcenter.

For objective functions with Lipschitz continuous first derivatives, the
linear fit in each Delaunay cell will converge to the average gradient in that
cell almost surely as more data is sampled.
Therefore, the Delaunay search will converge to a local minima.

This is experimental research code, and is not yet suitable for all use
cases.
In particular, if the global minimum is not in the convex hull of the provided
data set, then the Delaunay search can never reach it.

## Contents

Note that all test functions have been shifted so that their global minima
has value f(x) = 0.

 - delsearch.py contains the DelaunaySearch subroutine and an optional main
   program which uses a quadratic function and latin hypercube design for
   testing.
 - Rosenbrock.py contains a 2D test function with a single minima at the
   bottom of a shallow "banana shaped" Rosenbrock valley at (0,0).
 - Griewank.py contains a 2D test function with many local minima, and
   a single global minima at (0,0).
 - Michaelwicz.py contains a 2D test function with a single global minima
   at approximately (2.2, 1.57), amid large flat ridges and sharp narrow
   valleys.

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

To do the same for the Rosenbrock, Griewank, and Michaelwicz functions,
execute ``python3 Rosenbrock.py``, ``python3 Griewank.py``, and
``python3 Michaelwicz.py`` respectively.

## Author

* **Tyler H. Chang** - *Primary author*


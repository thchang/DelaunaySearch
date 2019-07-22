# The DelaunaySearch function optimizes a black-box objective function
# by refining over a Delaunay mesh.
#
# It is assumed that the user has already gathered some initial data by
# way of experimental search, which will be refined by the DelaunaySearch
# algorithm.
#
# In each iteration, DelaunaySearch guesses the objective value at the
# barycenter of each simplex in the umbrella neighborhood of the current
# minima, using the planar fit in its d+1 neighbors.
#
# The candidate simplex with the lowest average predicted value is subdivided
# in the subsequent iteration.
#
# This is experimental code and not yet robust for all use cases.
# In particular, if the true minima is outside the convex hull of the provided
# data set, DelaunaySearch cannot reach it since it only subdivides Delaunay
# cells.
#
# Input: obj_func is an objective function f : R^d --> R.
#
# Input: data is a python list containing at least d+1 data points, where
# the first d entries specify a data point and the last entry is the
# corresponding objective value.
#
# Returns: a list [minind, x_vals, f_vals] where minind is the index of the
# minimum value observed, x_vals are the sampled points (including the provided
# input data), and f_vals are the objective values at the x_vals.
#
# Author: Tyler Chang
# Last Update: July, 2019
#
def DelaunaySearch(data, obj_func, budget=10000):
   import numpy as np
   from scipy.spatial import Delaunay
   # Separate the design points and function values.
   np_x = np.array(data)[:,0:-1]
   np_f = np.array(data)[:,-1]
   # Get the dimensions and size of the initial data set.
   [n, d] = np_x.shape
   # Initialize the weight vector.
   center_weights = np.ones(d+1) / (d+1)
   # Initialize the current minima.
   minind = np.argmin(np_f)
   # Loop until the budget is exhausted.
   for i in range(budget):
      # Get the current mesh.
      mesh = Delaunay(np_x)

      # Get the umbrella neighborhood of the minind.
      neighborhood = []
      j = 0
      for elmt in mesh.simplices:
         if any(elmt == minind):
            neighborhood.append(j)
         j += 1

      # Reinitialize the list of points and predictions for this iteration.
      predictions = []
      # Loop over all elements in the umbrella neighborhood.
      for elmt in neighborhood:
         # Compute the center of the current element.
         center = np.transpose(np_x[mesh.simplices[elmt,:]]).dot(center_weights)
         predictions.append([center, 0.0])
         # Loop over all d+1 neighbors of the current element.
         j = 0
         for neighbor in mesh.neighbors[elmt,:]:
            # If on the convex hull, some facets will have no neighbor.
            if neighbor == -1:
               continue
            else:
               # Extrapolate to the current element's center from each neighbor.
               affine_w = np.zeros(d+1)
               affine_w[0:d] = mesh.transform[neighbor,:d,:d].dot(
		   center-mesh.transform[neighbor,d,:])
               affine_w[d] = 1.0 - sum(affine_w[0:d])
               predictions[-1][1] += np_f[mesh.simplices[neighbor]].dot(
                  affine_w)
               j += 1 # Count the number of true neighbors.
         predictions[-1][1] /= j
      # The next sample should be taken at the minimum predicted value.
      next_x = np.array(predictions)[np.argmin(np.array(predictions)[:,1]),0]
      next_f = obj_func(next_x)

      # Append x and f to the data arrays.
      np_x = np.append(np_x, [next_x], axis=0)
      np_f = np.append(np_f, next_f)

      # Update the current minima.
      minind = np.argmin(np_f)
   
   # The budget has exhausted. Return the collected data.
   return [minind, np_x, np_f]

def main():
   import matplotlib.pyplot as plt
   from scipy.spatial import Delaunay
   from pyDOE import lhs
   # Test DelaunaySearch on a quadratic.
   func = lambda x: sum([i ** 2 for i in x])
   # Generate a latin hypercube design in 10d.
   data = lhs(2, 10)
   data = (data - 0.5) * 2.0
   # Uncomment for a grid design in 3d (8 points).
   #data = [[-1, -1, -1], [-1,-1,1], [-1,1,-1], [-1,1,1],
   #        [1,-1,-1], [1,-1,1],[1,1,-1],[1,1,1]]
   data = data.tolist()
   # Generate function values.
   for i in range(len(data)):
      data[i].append(func(data[i]))
   # 50 iteration budget.
   [ind, x, f] = DelaunaySearch(data, func, budget=100)
   # Print the results.
   print('x = ', x[ind,:], '. f = ', f[ind], '.')
   # Display the final mesh on screen.
   mesh = Delaunay(x)
   plt.triplot(x[:,0], x[:,1], mesh.simplices.copy())
   plt.plot(x[:,0], x[:,1], 'o')
   plt.show()
   # Done exit.
   return

# Run as main if called directly.
if __name__ == "__main__":
   main()


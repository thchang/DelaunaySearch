# 2D Rosenbrock's valley function. Minima is at (0,0)
def Rosenbrock(x):
   return x[0]**2.0 + 100.0*(x[1] - x[0]**2.0)**2.0

# Test DelaunaySearch on the Rosenbrock valley.
def main():
   from delsearch import DelaunaySearch
   import matplotlib.pyplot as plt
   from scipy.spatial import Delaunay
   from pyDOE import lhs
   # Generate a latin hypercube design in 2d.
   data = lhs(2, 20)
   data = (data - 0.5) * 2.0
   data = data.tolist()
   # Populate function values.
   for i in range(len(data)):
      data[i].append(Rosenbrock(data[i]))
   # 100 iteration budget.
   [ind, x, f] = DelaunaySearch(data, Rosenbrock, budget=100)
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

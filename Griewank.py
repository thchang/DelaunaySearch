# 2D Griewank function with many local minima, one global minima at (0,0).
def Griewank(x):
   from math import cos as cos, sqrt as sqrt
   return (x[0]**2+x[1]**2)/4000 - cos(x[0])*cos(x[1])/sqrt(2) + 1.0

# Test DelaunaySearch on the Griewank function.
def main():
   from delsearch import DelaunaySearch
   import matplotlib.pyplot as plt
   from scipy.spatial import Delaunay
   from pyDOE import lhs
   # Generate a latin hypercube design in 2d.
   data = lhs(2, 20)
   data = (data - 0.5) * 20
   data = data.tolist()
   # Populate function values.
   for i in range(len(data)):
      data[i].append(Griewank(data[i]))
   # 100 iteration budget.
   [ind, x, f] = DelaunaySearch(data, Griewank, budget=100)
   # Print the results.
   print('x = ', x[ind,:], '. f = ', f[ind], '.')
   # Display the final data and mesh on screen.
   mesh = Delaunay(x)
   t = []
   for i in range(len(f) - len(data)):
      t.append(i)
   plt.subplot(1,2,1)
   plt.title('Iteration vs. Function Value')
   plt.yscale('log')
   plt.scatter(t,f[len(data):])
   plt.subplot(1,2,2)
   plt.title('Final Mesh')
   plt.triplot(x[:,0], x[:,1], mesh.simplices.copy())
   plt.plot(x[:,0], x[:,1], 'o')
   plt.show()
   # Done exit.
   return

# Run as main if called directly.
if __name__ == "__main__":
   main()

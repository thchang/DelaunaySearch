# 2D Michaelwicz's function with flat ridges and steep valleys.
# Minima is at ~(2.2,1.57).
def Michaelwicz(x):
   from math import sin as sin, pi as pi
   return -(sin(x[0]) * sin(x[0]**2/pi)**20 + sin(x[1]) * sin(2*x[1]**2/pi)**20)

# Test DelaunaySearch on the Michealwicz function.
def main():
   from delsearch import DelaunaySearch
   import matplotlib.pyplot as plt
   from scipy.spatial import Delaunay
   from pyDOE import lhs
   from math import pi
   # Generate a latin hypercube design in 2d.
   data = lhs(2, 20)
   data = (data) * pi
   data = data.tolist()
   # Populate function values.
   for i in range(len(data)):
      data[i].append(Michaelwicz(data[i]))
   # 100 iteration budget.
   [ind, x, f] = DelaunaySearch(data, Michaelwicz, budget=100)
   # Shift up so the minima is at zero.
   f[:] += 1.8013
   # Print the results.
   print('x = ', x[ind,:], '. f = ', f[ind], '.')
   print(f)
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

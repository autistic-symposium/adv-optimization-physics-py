# finite-difference implementation of linear advection with the Lax-Wendroff
# 
# We are solving a_t + u a_x = 0
#
# We run at several resolutions and compute the error.  This uses a
# cell-centered finite-difference grid
#
# M. Zingale (2013-03-12)

import numpy
import pylab
import sys

class ccFDgrid:

    def __init__(self, nx, ng, xmin=0.0, xmax=1.0):

        self.xmin = xmin
        self.xmax = xmax
        self.ng = ng
        self.nx = nx

        # python is zero-based.  Make easy intergers to know where the
        # real data lives
        self.ilo = ng
        self.ihi = ng+nx-1

        # physical coords -- cell-centered
        self.dx = (xmax - xmin)/(nx)
        self.x = xmin + (numpy.arange(nx+2*ng)-ng+0.5)*self.dx

        # storage for the solution
        self.a = numpy.zeros((nx+2*ng), dtype=numpy.float64)

    def scratchArray(self):
        """ return a scratch array dimensioned for our grid """
        return numpy.zeros((self.nx+2*self.ng), dtype=numpy.float64)

    def fillBCs(self):
        """ fill the a single ghostcell with periodic boundary conditions """
        self.a[self.ilo-1] = self.a[self.ihi]
        self.a[self.ihi+1] = self.a[self.ilo]

    def norm(self, e):
        """ return the norm of quantity e which lives on the grid """
        if not len(e) == (2*self.ng + self.nx):
            return None

        return numpy.sqrt(self.dx*numpy.sum(e[self.ilo:self.ihi+1]**2))


# define the speed
u = 1.0

N = [16, 32, 64, 128, 256, 512, 1024]
ng = 1

err = []

for nx in N:

    # create the grid
    g = ccFDgrid(nx, ng)

    # CFL number
    C = 0.8

    # time info
    dt = C*g.dx/u
    t = 0.0

    # 5 periods
    tmax = 5.0*(g.xmax - g.xmin)/u


    # initialize the data
    # tophat
    #g.a[numpy.logical_and(g.x >= 0.333, g.x <= 0.666)] = 1.0

    # gaussian
    g.a[:] = numpy.exp(-(g.x - 0.5)**2/0.1**2)

    ainit = g.a.copy()

    # evolution loop
    anew = g.scratchArray()
    
    while (t < tmax):

        # make sure we end right at tmax
        if (t + dt > tmax):
            dt = tmax - t
            C = dt*u/g.dx

        # fill the boundary conditions
        g.fillBCs()

        # loop over zones
        i = g.ilo
        while (i <= g.ihi):
        
            # Lax-Wendroff
            anew[i] = g.a[i] - 0.5*C*(g.a[i+1] - g.a[i-1]) + \
                0.5*C**2*(g.a[i+1] - 2.0*g.a[i] + g.a[i-1])
                
            i += 1

        # store the updated solution
        g.a[:] = anew[:]

        t += dt


    # compute the error
    err.append(g.norm(g.a - ainit))
    print g.dx, nx, err[-1]

    pylab.plot(g.x[g.ilo:g.ihi+1], g.a[g.ilo:g.ihi+1], label="N = %d" %(nx))

pylab.legend(frameon=False, fontsize="small")
pylab.xlabel("x")
pylab.savefig("lw.png")

pylab.clf()

N = numpy.array(N, dtype=numpy.float64)
err = numpy.array(err)

pylab.scatter(N, err, color="r")
pylab.plot(N, err[len(N)-1]*(N[len(N)-1]/N)**2, color="k", label="$\mathcal{O}(\Delta x^2)$")

print N
print err[0]*(N[0]/N)**2

pylab.legend(frameon=False)

ax = pylab.gca()
ax.set_xscale('log')
ax.set_yscale('log')

pylab.xlabel("x")
pylab.ylabel("L2 norm of absolute error")
pylab.savefig("lw-converge.png")









    

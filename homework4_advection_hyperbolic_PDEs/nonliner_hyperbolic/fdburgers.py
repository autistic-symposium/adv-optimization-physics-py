# finite-difference implementation of the inviscid Burger's equation
# compare conservative vs. non-conservative differencing
# 
# We are solving u_t + u u_x = 0
#
# M. Zingale (2013-04-10)

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
        self.u = numpy.zeros((nx+2*ng), dtype=numpy.float64)
        self.uinit = numpy.zeros((nx+2*ng), dtype=numpy.float64)

    def scratchArray(self):
        """ return a scratch array dimensioned for our grid """
        return numpy.zeros((self.nx+2*self.ng), dtype=numpy.float64)

    def fillBCs(self):
        """ fill the a single ghostcell with Neuamnn BCs """
        self.u[self.ilo-1] = self.u[self.ilo]
        self.u[self.ihi+1] = self.u[self.ihi]

    def norm(self, e):
        """ return the norm of quantity e which lives on the grid """
        if not len(e) == (2*self.ng + self.nx):
            return None

        return numpy.sqrt(self.dx*numpy.sum(e[self.ilo:self.ihi+1]**2))



def evolve(nx, C, tmax, conservative=1, init="shock"):

    ng = 1

    # create the grid
    g = ccFDgrid(nx, ng)

    # initial conditions
    if init == "shock":
        g.u[:] = 1.0
        g.u[g.x < 0.5] = 2.0
    elif init == "rarefaction":
        g.u[:] = 1.0
        g.u[g.x > 0.5] = 2.0


    # fill the boundary conditions
    g.fillBCs()

    t = 0.0

    g.uinit = g.u.copy()

    # evolution loop
    unew = g.scratchArray()
    
    while (t < tmax):

        # timestep
        dt = C*g.dx/max(g.u[g.ilo:g.ihi+1])


        # make sure we end right at tmax
        if (t + dt > tmax):
            dt = tmax - t


        # loop over zones
        i = g.ilo
        while (i <= g.ihi):
        
            if conservative == 1:
                unew[i] = g.u[i] - dt*(0.5*g.u[i]**2 - 0.5*g.u[i-1]**2)/g.dx
            else:
                unew[i] = g.u[i] - dt*g.u[i]*(g.u[i] - g.u[i-1])/g.dx

            i += 1

        # store the updated solution
        g.u[:] = unew[:]

        t += dt

        # fill the boundary conditions
        g.fillBCs()


    return g



#----------------------------------------------------------------------------- 
nx = 128
C = 0.5

pylab.clf()

tmax = 0.3

gc = evolve(nx, C, tmax, conservative=1, init="rarefaction")
gnc = evolve(nx, C, tmax, conservative=0, init="rarefaction")

pylab.plot(gc.x, gc.u, color="r", label="conservative")
pylab.plot(gnc.x, gnc.u, color="b", label="non-conservative")


tmax = 0.1

gc = evolve(nx, C, tmax, conservative=1, init="rarefaction")
gnc = evolve(nx, C, tmax, conservative=0, init="rarefaction")

pylab.plot(gc.x, gc.u, color="r", ls="--")
pylab.plot(gnc.x, gnc.u, color="b", ls="--")




pylab.plot(gnc.x, gc.uinit, color="0.5", label="initial conditions")

pylab.xlim(0,1)

#pylab.legend(loc=3, frameon=False, fontsize="small")

pylab.xlabel("x")
pylab.ylabel("u")

pylab.savefig("burger-rarefaction.png")


#-----------------------------------------------------------------------------
tmax = 0.3

pylab.clf()

gc = evolve(nx, C, tmax, conservative=1, init="shock")
gnc = evolve(nx, C, tmax, conservative=0, init="shock")

pylab.plot(gc.x, gc.u, color="r", label="conservative")
pylab.plot(gnc.x, gnc.u, color="b", label="non-conservative")

# analytic shock position, assuming u_l = 2, u_r = 1, and initial
# discontinuity at 0.5
S = 0.5*(1.0 + 2.0)
pylab.plot([S*tmax + 0.5, S*tmax + 0.5], [1,2], 
           color="0.5", ls=":", label="analytic shock position")


tmax = 0.1

gc = evolve(nx, C, tmax, conservative=1, init="shock")
gnc = evolve(nx, C, tmax, conservative=0, init="shock")

pylab.plot(gc.x, gc.u, color="r", ls="--")
pylab.plot(gnc.x, gnc.u, color="b", ls="--")




pylab.plot(gnc.x, gc.uinit, color="0.5", label="initial conditions")

pylab.xlim(0,1)

#pylab.legend(loc=3, frameon=False, fontsize="small")

pylab.xlabel("x")
pylab.ylabel("u")

pylab.savefig("burger-shock.png")

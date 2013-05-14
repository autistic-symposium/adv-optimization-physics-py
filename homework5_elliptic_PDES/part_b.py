"""
    We solve the poisson equation for phi using Gauss-Seidel relaxation.
    The multigrid for cell-centered data, nx, ny is a power of 2 and 
    we make dx = dy
"""

import numpy
import pylab



def rightSide(Ustar, dx, dy, nx, ny, ng, ilo, ihix, ihiy):
    """ takes the vector and gives the divergent of it,
        for out problem, this is the right side of the poisson equation """
    u = Ustar[0].copy()
    v = Ustar[1].copy()

    fillBC(u, ilo, ihix, ihiy)
    fillBC(v, ilo, ihix, ihiy)

    DUij= numpy.zeros((nx + 2*ng, ny + 2*ng), dtype=numpy.float64)
    for i in range(ilo,ihix+1):
        for j in range(ilo,ihiy+1):  
            DUij[i, j]= numpy.array( (u[i+1,j] - u[i-1,j] )/ (2.0*dx) +  ( v[i,j+1] - v[i,j-1] )/ (2.0*dy) )
    
    fillBC(DUij, ilo, ihix, ihiy)

    return DUij



def fillBC(f, ilo, ihix, ihiy):
    """ in a square with periodic boundary conditions, 
        we want the right part to meet the left and the top to meet the bottom"""
    f[ilo-1,:] =  f[ihix,:]
    f[ihix+1,:] =  f[ilo,:]
    f[:,ilo-1] =  f[:, ihiy]
    f[:, ihiy+1] =  f[:, ilo] 

    return f



def computeResidual(ilo, ihix, ihiy, dx, dy, phi, f):
    """ calculates the residual """
    r = numpy.zeros( (len(phi),len(phi)), dtype=numpy.float64)

    r[ilo:ihix+1,ilo:ihix+1] = \
            f[ilo:ihix+1,ilo:ihiy+1] - \
            (phi[ilo-1:ihix  ,ilo  :ihiy+1] + \
             phi[ilo+1:ihix+2,ilo  :ihiy+1] - \
             2.0*phi[ilo:ihix+1,ilo:ihiy+1])/(dx*dx) - \
            (phi[ilo  :ihix+1,ilo-1:ihiy  ] +\
             phi[ilo  :ihix+1,ilo+1:ihiy+2] -\
             2.0*phi[ilo:ihix+1,ilo:ihiy+1])/(dy*dy)
    
    return r



def error(ilo, ihix, ihiy, dx, r):
    """ calculates the L2 norm of elements in the residual """
    e = numpy.sqrt(dx*dx*numpy.sum((r[ilo:ihix+1,ilo:ihiy+1]**2)))
    
    return e


   
def doPlot_phi(phi, phi_true, n, name):   
    """ plot phi numerical and analytical and for many iterations
        to follow the evolution and convergence of the process"""
  
    pylab.clf()
    pylab.cla()
    
    f = pylab.figure() 
    f.text(.5, .95, r"$\phi = \frac{1}{10} \cos(2\pi y) \cos(2\pi x)$ - Iteration #%d  " %n,  horizontalalignment='center')

    pylab.subplot(221)  
    pylab.title("Analytical " r"$\phi$ in 2D" , size = 8)
    pylab.imshow(phi_true)
    pylab.colorbar()
    pylab.xlim(1,32)
    pylab.ylim(1,32)
    pylab.xlabel("Analytical " r"$\phi$ in 1D:" , size = 8)
    pylab.ylabel("# of cells", size = 8)
    
    pylab.subplot(223)
    pylab.plot(phi_true, color="r")
    pylab.xlabel("$ x$", size = 8)
    pylab.ylabel(r"$ \phi(x)$", size = 8)
    pylab.xlim(1,32)
    pylab.ylim(-0.1,0.1)
    
    pylab.subplot(222)  
    pylab.title("Numerical " r"$\phi $ in 2D", size = 8)
    pylab.imshow(phi)
    pylab.xlim(1,32)
    pylab.ylim(1,32)
    pylab.colorbar()
    pylab.xlabel("Numerical " r"$\phi $ in 1D", size = 8)
    
    pylab.subplot(224)
    pylab.plot(phi, color="r")
    pylab.xlabel("$ x$", size = 8)
    pylab.xlim(1,32)
    pylab.ylim(-0.1,0.1)

    outfile = "plots/item_b_" + name + "_num_vs_ana_%d.png" %(n)
    pylab.savefig(outfile)

    return 0



def smoothRun(Ustar, phi_true, nx, ny, ng, dx, dy,DO_PLOTS):

    ilo = ng
    ihix = ng + nx - 1
    ihiy = ng + ny - 1
    
    """ initialize the solution to zero """
    phi = numpy.zeros((nx + 2*ng, ny + 2*ng), dtype=numpy.float64)
    
    """ initialize the RHS with DUstar """
    DUrhs =  rightSide(Ustar, dx, dy, nx, ny, ng, ilo, ihix, ihiy)

    """ error """
    eps = 0.5* 10**(-7)
    norm_S = error(ilo, ihix, ihiy, dx, DUrhs)
    eps = eps*norm_S

    norm_r = 2*eps
    n = 0
    
    """ start smoothing: Gauss-Seidel - here we assume dx = dy 
        and nx and ny is a power of 2  """
    while (norm_r > eps):
        n += 1
        
        """ odd """
        phi[ilo:ihix:2, ilo:ihiy+1:2 ] = \
                0.25*(-dx*dx*DUrhs[ilo:ihix+1:2, ilo:ihiy+1:2]+ \
                phi[ilo+1:ihix+2:2,ilo:ihiy+1:2] + \
                phi[ilo-1:ihix :2, ilo:ihiy+1:2] + \
                phi[ilo:ihix+1:2,ilo+1:ihiy+2:2] + \
                phi[ilo  :ihix+1:2, ilo-1:ihiy:2])

        phi[ilo+1:ihix+1:2,ilo+1:ihiy+1:2] = \
                 0.25*(-dx*dx*DUrhs[ilo+1:ihix+1:2,ilo+1:ihiy+1:2] +\
                 phi[ilo+2:ihix+2:2,ilo+1:ihiy+1:2] + \
                 phi[ilo  :ihix  :2,ilo+1:ihiy+1:2] + \
                 phi[ilo+1:ihix+1:2,ilo+2:ihiy+2:2] + \
                 phi[ilo+1:ihix+1:2,ilo  :ihiy  :2])

        """ fill the ghost cells """
        fillBC(phi, ilo, ihix, ihiy)

        """ even """
        phi[ilo+1:ihix+1:2,ilo:ihiy+1:2] = \
                 0.25*(-dx*dx*DUrhs[ilo+1:ihix+1:2,ilo:ihiy+1:2] + \
                 phi[ilo+2:ihix+2:2,ilo  :ihiy+1:2] + \
                 phi[ilo  :ihix  :2,ilo  :ihiy+1:2] + \
                 phi[ilo+1:ihix+1:2,ilo+1:ihiy+2:2] + \
                 phi[ilo+1:ihix+1:2,ilo-1:ihiy  :2]) \

        phi[ilo:ihix:2,ilo+1:ihiy+1:2] = \
                0.25*(-dx*dx*DUrhs[ilo:ihix+1:2,ilo+1:ihiy+1:2] +
                 phi[ilo+1:ihix+2:2,ilo+1:ihiy+1:2] + \
                 phi[ilo-1:ihix  :2,ilo+1:ihiy+1:2] + \
                 phi[ilo  :ihix+1:2,ilo+2:ihiy+2:2] + \
                 phi[ilo  :ihix+1:2,ilo  :ihiy  :2])

        """ fill the ghost cells """
        fillBC(phi, ilo, ihix, ihiy)  
        
        """ calculates the residual and norm of residual """
        resid = computeResidual(ilo, ihix, ihiy, dx, dy, phi, DUrhs)
        norm_r = error(ilo, ihix, ihiy, dx, resid)
    
        """ plot some of the iterations """
        if (DO_PLOTS == 1):
            if  (n <= 10) | (n > 10) & (n%50 == 0.0):
                doPlot_phi(phi, phi_true, n, "phi") 
    
    """ ghetto way of unpadding"""
    phi_final = phi[ilo:ihix+1]
    
    print "The total number of iterations for the desired |r| < %.8f is %d." %(norm_r, n)
    
    return phi_final




def doPartB(Ustar, phi_true, nx, ny,  ng, dx, dy, DO_PLOTS):
    
    phi_final = smoothRun(Ustar, phi_true, nx, ny, ng, dx, dy,DO_PLOTS)
    
    return phi_final
    

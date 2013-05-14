"""
2-d approximate projection on a cell-centered grid with periodic BCs 

M. Zingale (2013-04-07)

"""

import sys
import numpy
import pylab
import math

class grid:

    def __init__(self, nx, ny, ng=1, 
                 xmin=0.0, xmax=1.0,
                 ymin=0.0, ymax=1.0):

        self.nx = nx
        self.ny = ny
        self.ng = ng
        
        self.ilo = ng
        self.ihi = ng+nx-1

        self.jlo = ng
        self.jhi = ng+ny-1

        self.dx = (xmax - xmin)/nx
        self.dy = (ymax - ymin)/ny

        self.x = (numpy.arange(nx+2*ng) - ng + 0.5)*self.dx + xmin
        self.y = (numpy.arange(ny+2*ng) - ng + 0.5)*self.dy + ymin

        x2d = numpy.repeat(self.x, 2*ng+ny)
        x2d.shape = (2*ng+nx, 2*ng+ny)
        self.x2d = x2d

        y2d = numpy.repeat(self.y, 2*ng+nx)
        y2d.shape = (2*ng+ny, 2*ng+nx)
        y2d = numpy.transpose(y2d)
        self.y2d = y2d


        self.phi = numpy.zeros((nx + 2*ng, ny + 2*ng), dtype=numpy.float64)


    def scratchArray(self):
        return numpy.zeros((self.nx + 2*self.ng, self.ny + 2*self.ng), dtype=numpy.float64)


    def fillBC(self, v):
        """ fill periodic BCs """

        if not self.ng == 1:
            sys.exit("invalid ng")

        v[self.ilo-1,:] = v[self.ihi,:]
        v[self.ihi+1,:] = v[self.ilo,:]

        v[:,self.jlo-1] = v[:,self.jhi]
        v[:,self.jhi+1] = v[:,self.jlo]

        
    def norm(self, e):
        """ L2 norm """

        #return numpy.sqrt(self.dx*self.dy*numpy.sum((e[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]**2).flat))
        return numpy.max(numpy.abs(e))

        
def residual(gr, f):
    """ compute the residual """

    ib = gr.ilo
    ie = gr.ihi
    jb = gr.jlo
    je = gr.jhi

    r = gr.scratchArray()

    r[ib:ie+1,jb:je+1] = f[ib:ie+1,jb:je+1] - \
        (gr.phi[ib-1:ie,jb:je+1] - 2.0*gr.phi[ib:ie+1,jb:je+1] + gr.phi[ib+1:ie+2,jb:je+1])/gr.dx**2 - \
        (gr.phi[ib:ie+1,jb-1:je] - 2.0*gr.phi[ib:ie+1,jb:je+1] + gr.phi[ib:ie+1,jb+1:je+2])/gr.dy**2 
    
    return r



def smooth(gr, f, eps):

    ib = gr.ilo
    ie = gr.ihi
    jb = gr.jlo
    je = gr.jhi

    fnorm = gr.norm(f)
    print fnorm

    rnorm = 1.e10

    i = 0
    while (rnorm > eps*fnorm and i < 25000):

        gr.fillBC(gr.phi)

        gr.phi[ib:ie+1:2,jb:je+1:2] = \
            0.25*(gr.phi[ib+1:ie+2:2,jb  :je+1:2] + \
                  gr.phi[ib-1:ie  :2,jb  :je+1:2] + \
                  gr.phi[ib  :ie+1:2,jb+1:je+2:2] + \
                  gr.phi[ib  :ie+1:2,jb-1:je  :2] - \
                  gr.dx**2*f[ib:ie+1:2,jb:je+1:2])

        gr.phi[ib+1:ie+1:2,jb+1:je+1:2] = \
            0.25*(gr.phi[ib+2:ie+2:2,jb+1:je+1:2] + \
                  gr.phi[ib  :ie  :2,jb+1:je+1:2] + \
                  gr.phi[ib+1:ie+1:2,jb+2:je+2:2] + \
                  gr.phi[ib+1:ie+1:2,jb  :je  :2] - \
                  gr.dx**2*f[ib+1:ie+1:2,jb+1:je+1:2])

        gr.fillBC(gr.phi)

        gr.phi[ib+1:ie+1:2,jb:je+1:2] = \
            0.25*(gr.phi[ib+2:ie+2:2,jb  :je+1:2] + \
                  gr.phi[ib  :ie  :2,jb  :je+1:2] + \
                  gr.phi[ib+1:ie+1:2,jb+1:je+2:2] + \
                  gr.phi[ib+1:ie+1:2,jb-1:je  :2] - \
                  gr.dx**2*f[ib+1:ie+1:2,jb:je+1:2])

        gr.phi[ib:ie+1:2,jb+1:je+1:2] = \
            0.25*(gr.phi[ib+1:ie+2:2,jb+1:je+1:2] + \
                  gr.phi[ib-1:ie  :2,jb+1:je+1:2] + \
                  gr.phi[ib  :ie+1:2,jb+2:je+2:2] + \
                  gr.phi[ib  :ie+1:2,jb  :je  :2] - \
                  gr.dx**2*f[ib:ie+1:2,jb+1:je+1:2])


        rnorm = gr.norm(residual(gr, f))
        #print rnorm, gr.phi[gr.nx/2,gr.ny/2]

        i += 1


#-----------------------------------------------------------------------------
# an analytic test problem for the smoothing

def true(x,y):
    pi = math.pi
    return numpy.sin(2.0*pi*x)**2*numpy.cos(4.0*pi*y) + \
        numpy.sin(4.0*pi*x)*numpy.cos(2.0*pi*y)**2


# the righthand side                                                            
def frhs(x,y):
    pi = math.pi
    return 8.0*pi**2*numpy.cos(4.0*pi*y)*(numpy.cos(4.0*pi*x) -
                                          numpy.sin(4.0*pi*x)) - \
           16.0*pi**2*(numpy.sin(4.0*pi*x)*numpy.cos(2.0*pi*y)**2 +
                       numpy.sin(2.0*pi*x)**2 * numpy.cos(4.0*pi*y))

#-----------------------------------------------------------------------------
# projection example

def udivfree(gr):
    u = -numpy.sin(math.pi*gr.x2d)**2*numpy.sin(2.0*math.pi*gr.y2d)
    v =  numpy.sin(math.pi*gr.y2d)**2*numpy.sin(2.0*math.pi*gr.x2d)

    return u, v


def phif(gr):
    return 0.1*numpy.cos(2.0*math.pi*gr.y2d)*numpy.cos(2.0*math.pi*gr.x2d)


def gradphi(gr, phi):
    gphi_x = gr.scratchArray()
    gphi_y = gr.scratchArray()

    gphi_x[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1] = \
        0.5*(phi[gr.ilo+1:gr.ihi+2,gr.jlo:gr.jhi+1] - \
             phi[gr.ilo-1:gr.ihi  ,gr.jlo:gr.jhi+1])/gr.dx

    gphi_y[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1] = \
        0.5*(phi[gr.ilo:gr.ihi+1,gr.jlo+1:gr.jhi+2] - \
             phi[gr.ilo:gr.ihi+1,gr.jlo-1:gr.jhi  ])/gr.dy

    return gphi_x, gphi_y


def divU(gr, u, v):
    dU = gr.scratchArray()
    dU[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1] = \
        0.5*(u[gr.ilo+1:gr.ihi+2,gr.jlo:gr.jhi+1] - \
             u[gr.ilo-1:gr.ihi  ,gr.jlo:gr.jhi+1])/gr.dx + \
        0.5*(v[gr.ilo:gr.ihi+1,gr.jlo+1:gr.jhi+2] - \
             v[gr.ilo:gr.ihi+1,gr.jlo-1:gr.jhi  ])/gr.dy

    return dU


#-----------------------------------------------------------------------------

nx = 16
ny = 16

gr = grid(nx, ny)

gr.phi[:] = 0.0


# get the original divergence-free field
ud, vd = udivfree(gr)

udOrig = ud.copy()
vdOrig = vd.copy()


# add the gradient of a scalar
phi = phif(gr)
phiOrig = phi.copy()

gpx, gpy = gradphi(gr, phi)

# pollute the velocity field
ud += gpx
gr.fillBC(ud)

vd += gpy
gr.fillBC(vd)


f = divU(gr, ud, vd)

smooth(gr, f, 1.e-7)

# correct the velocity field to recover the divergence free part
gpx, gpy = gradphi(gr, gr.phi)

unew = ud - gpx
vnew = vd - gpy


#-----------------------------------------------------------------------------
# plots

pylab.subplot(131)

pylab.imshow(numpy.transpose(phiOrig[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]),
             interpolation="nearest", origin="lower")
pylab.colorbar()
pylab.title("original phi")


pylab.subplot(132)

pylab.imshow(numpy.transpose(gr.phi[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]),
             interpolation="nearest", origin="lower")
pylab.colorbar()
pylab.title("new phi")

pylab.subplot(133)

phie = gr.phi - phiOrig

pylab.imshow(numpy.transpose(phie[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]),
             interpolation="nearest", origin="lower")
pylab.colorbar()
pylab.title("phi error")


f = pylab.gcf()
f.set_size_inches(10.0,4.0)

pylab.tight_layout()

pylab.savefig("project-phi.png")


# velocities

pylab.subplot(321)

pylab.imshow(numpy.transpose(udOrig[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]),
             interpolation="nearest", origin="lower")
pylab.colorbar()
pylab.title("original u")


pylab.subplot(322)

pylab.imshow(numpy.transpose(vdOrig[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]),
             interpolation="nearest", origin="lower")
pylab.colorbar()
pylab.title("original v")


pylab.subplot(323)

pylab.imshow(numpy.transpose(ud[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]),
             interpolation="nearest", origin="lower")
pylab.colorbar()
pylab.title("\'polluted\' u")


pylab.subplot(324)

pylab.imshow(numpy.transpose(vd[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]),
             interpolation="nearest", origin="lower")
pylab.colorbar()
pylab.title("\'polluted\' v")


pylab.subplot(325)

pylab.imshow(numpy.transpose(unew[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]),
             interpolation="nearest", origin="lower")
pylab.colorbar()
pylab.title("projected u")


pylab.subplot(326)

pylab.imshow(numpy.transpose(vnew[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]),
             interpolation="nearest", origin="lower")
pylab.colorbar()
pylab.title("projected v")

f = pylab.gcf()
f.set_size_inches(6.0,8.0)

pylab.tight_layout()


pylab.savefig("project-u.png")


# compute the error
eu = unew - udOrig
ev = vnew - vdOrig

print "Nx, Ny, L2 norm of error (u, v): ", nx, ny, \
    numpy.sqrt(gr.dx*gr.dy*numpy.sum((eu[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]**2).flat)), \
    numpy.sqrt(gr.dx*gr.dy*numpy.sum((ev[gr.ilo:gr.ihi+1,gr.jlo:gr.jhi+1]**2).flat))

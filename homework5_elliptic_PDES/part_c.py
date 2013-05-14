"""
    We recover the original divergence-free velocity field via 
        Ud,new = Ustar - Gphi
""" 

import numpy
import pylab
import operator


def do_plots_c(Ud, Unew):
    """ plot Ud,new and Ud with zoom on the bug """ 
    pylab.clf()
    pylab.cla()
    
    f = pylab.figure() 
    f.text(.5, .95, r"$U_{\rm d}$ (left) and $U_{\rm d, new}$ (right) ",  horizontalalignment='center')


    pylab.subplot(221)
    pylab.imshow(Ud[0])
    pylab.ylabel("# of cells", size =8)

    
    pylab.subplot(223)
    pylab.imshow(Ud[1])
    pylab.xlim(1,32)
    pylab.xlabel("# of cells", size =8)
    pylab.ylabel("# of cells", size =8)

    pylab.subplot(222)
    pylab.imshow(Unew[0])
    pylab.ylabel("# of cells", size =8)
    
    pylab.subplot(224) 
    pylab.imshow(Unew[1])
    pylab.xlim(1,32)
    pylab.xlabel("# of cells", size =8)
    pylab.ylabel("# of cells", size =8)

    pylab.savefig("plots/item_c_Udnew.png")


def doPartC(Ustar, phi_num, Ud, nx, ny, xmin, xmax, ymin, ymax, DO_PLOTS):
    """ coordinates of centers """
    dx = (xmax - xmin)/nx
    dy = (ymax - ymin)/ny
    
    """ calcuates the new gradient"""
    Gphi = numpy.gradient(phi_num, dx, dy)


    """ recover Ud, new """
    Unew = map(operator.sub, Ustar,Gphi)
    
    
    if (DO_PLOTS == 1):
        do_plots_c(Ud, Unew)
    
    return 0

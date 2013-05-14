import sys
from util import runparams
import mesh.patch as patch
import numpy
from util import msg

def initData(myPatch):
    """ initialize the Rayleigh-Taylor instability problem """

    msg.bold("initializing the Rayleigh-Taylor instability problem...")

    # make sure that we are passed a valid patch object
    if not isinstance(myPatch, patch.ccData2d):
        print "ERROR: patch invalid in raytay.py"
        print myPatch.__class__
        sys.exit()

    # get the density, momenta, and energy as separate variables
    dens = myPatch.getVarPtr("density")
    xmom = myPatch.getVarPtr("x-momentum")
    ymom = myPatch.getVarPtr("y-momentum")
    ener = myPatch.getVarPtr("energy")

    # get the other input parameters for the problem
    gamma = runparams.getParam("eos.gamma")
    grav = runparams.getParam("compressible.grav")
    dens_up = runparams.getParam("raytay.dens_up")
    dens_down = runparams.getParam("raytay.dens_down")
    A = runparams.getParam("raytay.pert_amplitude_factor")
    p0 = runparams.getParam("raytay.pressure_bottom")
    sigma = runparams.getParam("raytay.sigma")


    # get the grid parameters    
    xmin = runparams.getParam("mesh.xmin")
    xmax = runparams.getParam("mesh.xmax")
    ymin = runparams.getParam("mesh.ymin")
    ymax = runparams.getParam("mesh.ymax")   
    
    yctr = 0.5*(ymin + ymax)


    # initialize the components, remember, that ener here is
    # rho*eint + 0.5*rho*v**2, where eint is the specific
    # internal energy (erg/g)
    xmom[:,:] = 0.0
    ymom[:,:] = 0.0
    dens[:,:] = 0.0
    
    Lx = xmax - xmin

    # set the density and energy to be stratified in the y-direction
    myg = myPatch.grid

    j = myg.jlo
    while j <= myg.jhi:
        if myg.y[j] < yctr :
            dens[:,j] = dens_down
            pres = dens_down*grav*myg.y[j] + p0 
            ener[:,j] = pres/(gamma - 1.0) + 0.5*(xmom[:,j]**2 + ymom[:,j]**2)/dens_down
        else:
            dens[:,j] = dens_up
            pres = p0 + dens_down*grav*yctr + dens_up*grav*(myg.y[j] - yctr)
            ener[:,j] = pres/(gamma - 1.0) + 0.5*(xmom[:,j]**2 + ymom[:,j]**2)/dens_up
        j += 1


    i = myg.ilo
    while i <= myg.ihi:

        j = myg.jlo
        while j <= myg.jhi:

            v_pert = A*numpy.sin(2.0*numpy.pi*myg.x[i]/Lx)*numpy.exp( - (( myg.y[j]-yctr)/sigma)**2 ) 
            
            if myg.y[j] <=  yctr:
                
                # momentum
                ymom[i,j] = dens[i,j]*v_pert
                #xmom[i,j] = 0.0

                # internal energy
                eint = (ener[i,j] - 0.5*(xmom[i,j]**2 - ymom[i,j]**2)/dens[i,j])/dens[i,j]
                
                #pressure
                #pres =  p0 + dens_down*grav*myg.y[j] 
                pres = dens[i,j]*eint*(gamma - 1.0)

                #density
                dens[i,j] = pres/(eint*(gamma - 1.0))
                
                #energy
                ener[i,j] = pres/(gamma - 1.0) + 0.5*(xmom[i,j]**2 + ymom[i,j]**2)/dens[i,j]
            
            else:
                # momentum
                ymom[i,j] = dens[i,j]*v_pert
                xmom[i,j] = 0.0
                
                # internal energy
                eint = (ener[i,j] - 0.5*(xmom[i,j]**2 - ymom[i,j]**2)/dens[i,j])/dens[i,j]
                
                #pressure
                #pres =  p0 + dens_down*grav*yctr + dens_up*grav*(myg.y[j] - yctr)
                pres = dens[i,j]*eint*(gamma - 1.0) 
                
                #density
                dens[i,j] = pres/(eint*(gamma - 1.0))
                
                #energy
                ener[i,j] = pres/(gamma - 1.0) + 0.5*(xmom[i,j]**2 + ymom[i,j]**2)/dens[i,j]
                
                
            j += 1
        i += 1
        
    

    
                             

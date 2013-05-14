

    
    
    



    """ Start Iterations """
    
    i = myg.ilo
    while i <= myg.ihi:

        j = myg.jlo
        while j <= myg.jhi: 

                if myg.y[j] <= yctr : #+ A*numpy.sin(2.0*numpy.pi*myg.x[i]/Lx)*numpy.exp( - (myg.y[j]-yctr)**2/sigma**2 ):
                    print "big"
                    pres = 1.0# p0 + dens_down*grav*yctr + dens_up*grav*(myg.y[i] - yctr)
                    
                    xmom[i,j] = 0.0
                    ymom[i,j] = 1.0# A*numpy.sin(2.0*numpy.pi*myg.x[i]/Lx)*numpy.exp( - (myg.y[j]-yctr)**2/sigma**2 )
  
                    ener[i,j] = 1.0#pres + (0.5*ymom[i,j]**2)/dens_up
                
                else:
                    print "small"
                    pres = 1.0#p0 + dens_down*grav*myg.y[j]
                    
                    xmom[i,j] = 0.0
                    ymom[i,j] = 1.0#A*numpy.sin(2*numpy.pi*myg.x[i]/Lx)*numpy.exp( - (myg.y[j]-yctr)**2/sigma**2 )
           
                    ener[i,j] = 1.0#pres + (0.5*ymom[i,j]**2)/dens_up

                j += 1
                print j
        i += 1
        print i

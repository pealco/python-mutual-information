from numpy import *
from megprocess import *
import numpy
import time

def entropy(counts):
    '''Compute entropy.'''
    ps = counts/float(sum(counts))  # coerce to float and normalize
    ps = ps[nonzero(ps)]            # toss out zeros
    H = -sum(ps * numpy.log2(ps))   # compute entropy
    
    return H

def mi(x, y, bins):
    '''Compute mutual information'''
    counts_xy = histogram2d(x, y, bins=bins, range=[[0, 311335], [0, 311335]])[0]
    counts_x  = histogram(  x,    bins=bins, range=[0, 311335])[0]
    counts_y  = histogram(  y,    bins=bins, range=[0, 311335])[0]
    
    H_xy = entropy(counts_xy)
    H_x  = entropy(counts_x)
    H_y  = entropy(counts_y)
    
    return H_x + H_y - H_xy

def mirror(data):
    n_chan = shape(data)[0]
    for c1 in range(n_chan):
        for c2 in range(n_chan):
            if c1 > c2:
                data[c1,c2] = data[c2,c1]
                data[c1,c2] = data[c2,c1]
    return data

def cross_mi(epochs, interest_region, n_chan, bins):
    mi_ldif = zeros((n_chan, n_chan))
    mi_dlif = zeros((n_chan, n_chan))
    for ch1 in range(n_chan):
        print ch1
        for ch2 in range(n_chan):
            if ch2 >= ch1:
                ldif_a = epochs[7,interest_region,ch1,:].flatten()
                ldif_b = epochs[7,interest_region,ch2,:].flatten()
                dlif_a = epochs[3,interest_region,ch1,:].flatten()
                dlif_b = epochs[3,interest_region,ch2,:].flatten()
                mi_ldif[ch1, ch2] = mi(ldif_a, ldif_b, bins=bins)
                mi_dlif[ch1, ch2] = mi(dlif_a, dlif_b, bins=bins)
    return mirror(mi_ldif), mirror(mi_dlif)

def auto_mi(epochs, total_offset, interest_region, n_chan, bins):
    ami_ldif = zeros((n_chan, total_offset))
    ami_dlif = zeros((n_chan, total_offset))
    for ch in range(n_chan):
        print ch
        for tau in range(total_offset):
            ldif_a = epochs[7, interest_region,     ch, :].flatten()
            ldif_b = epochs[7, interest_region+tau, ch, :].flatten()
            dlif_a = epochs[3, interest_region,     ch, :].flatten()
            dlif_b = epochs[3, interest_region+tau, ch, :].flatten()
            ami_ldif[ch, tau] = mi(ldif_a, ldif_b, bins=bins)
            ami_dlif[ch, tau] = mi(dlif_a, dlif_b, bins=bins)
    return ami_ldif, ami_dlif

if __name__ == "__main__":
    
    epochs = numpy.load("epochs.npy")
    
    import psyco
    histogram2d = psyco.proxy(numpy.histogram2d)
    
    # Transform the data to emphasize power. 
    # This is something to play around with.
    epochs = abs(epochs)**2 
    
    # Constants.
    bins = 64
    n_chan = 157
    total_offset = 150
    interest_region = arange(250, 450+1)
    
    print "Computing cross mutual information"
    mi_ldif, mi_dlif = cross_mi(epochs, interest_region, n_chan, bins)
    
    #print "Computing auto mutual information with time shifts"
    #ami_ldif, ami_dlif = auto_mi(epochs, total_offset, interest_region, n_chan, bins)
        
    savez("midata.npz", mi_ldif, mi_dlif)#, ami_ldif, ami_dlif)
        
    figure(1)
    imshow(mi_ldif, interpolation="nearest", origin="lower")
    colorbar()                                             
    title("Ungrammatical condition")                       
    xlabel("Channels")                                     
    ylabel("Channels")                                     
                                                           
    figure(2)                                              
    imshow(mi_dlif, interpolation="nearest", origin="lower")
    colorbar()
    title("Grammatical condition")
    xlabel("Channels")
    ylabel("Channels")
    
    #mi_means_ldif = zeros(n_chan)
    #mi_means_dlif = zeros(n_chan)
    #for ch in arange(n_chan):
    #    mi_means_ldif[ch] = mean(mi_ldif[ch, :])
    #    mi_means_dlif[ch] = mean(mi_dlif[ch, :])
    #
    #figure(3)
    #bar(range(n_chan), mi_means_ldif)
    #sorted_chans_ldif = argsort(mi_means_ldif)
    #best = sorted_chans_ldif[-1]
    #print "ldif", sorted_chans_ldif
    #
    #figure(4)
    #bar(range(n_chan), mi_means_dlif)
    #sorted_chans_dlif = argsort(mi_means_dlif)
    #best = sorted_chans_dlif[-1]
    #print "dlif", sorted_chans_dlif
    #
    #
    ##
 
            
    #figure(5)
    #plot(offsets, mi_phase_ldif, 'o', label="Ungrammatical")
    #plot(offsets, mi_phase_dlif, 'o', label="Grammatical")
    #axis('tight')
    #title("Time-shifted auto-mutual information")
    #xlabel("Time shift (ms)")
    #ylabel("Mutual information")
    #legend(loc=2)
    #
    #print "ldif = ", mi_ldif[best,:]
    #print "dlif = ", mi_dlif[best,:]
    
    #print "ldif = ", mi_ldif[19,:]
    #print "dlif = ", mi_dlif[19,:]
    
    show()
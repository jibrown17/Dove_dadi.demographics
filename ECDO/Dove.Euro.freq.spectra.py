##################################################################################################
# File: single_fs_from_nex.py
# For use in python v. 3.7
#
# Created by Joshua I. Brown in August, 2020
#
# Description: generate single-population allele frequency spectra from SNPs from Nexus dataset
#
# Usage: python fs_from_nex_GBS.py
#
# NOTE: This script depends on a modified version of the custom script "nex_mo.py" originally 
# written by R. Gutenkunst specifically for use with alignments in nexus format.
##################################################################################################



import dadi
import nex_mo



if __name__ == '__main__':
    pop_assignments = {
'DOV128a':'Euro',
'DOV128b':'Euro',
'DOV146a':'Euro',
'DOV146b':'Euro',
'DOV096a':'Euro',
'DOV096b':'Euro',
'DOV097a':'Euro',
'DOV097b':'Euro',
'DOV123a':'Euro',
'DOV123b':'Euro',
'DOV124a':'Euro',
'DOV124b':'Euro',
'DOV125a':'Euro',
'DOV125b':'Euro',
'DOV126a':'Euro',
'DOV126b':'Euro',
'DOV127a':'Euro',
'DOV127b':'Euro',
'DOV129a':'Euro',
'DOV129b':'Euro',
'DOV130a':'Euro',
'DOV130b':'Euro',
'DOV131a':'Euro',
'DOV131b':'Euro',
'DOV132a':'Euro',
'DOV132b':'Euro',
'DOV133a':'Euro',
'DOV133b':'Euro',
'DOV134a':'Euro',
'DOV134b':'Euro',
'DOV135a':'Euro',
'DOV135b':'Euro',
'DOV136a':'Euro',
'DOV136b':'Euro',
'DOV137a':'Euro',
'DOV137b':'Euro',
'DOV138a':'Euro',
'DOV138b':'Euro',
'DOV139a':'Euro',
'DOV139b':'Euro',
'DOV140a':'Euro',
'DOV140b':'Euro',
'DOV141a':'Euro',
'DOV141b':'Euro',
'DOV142a':'Euro',
'DOV142b':'Euro',
'DOV143a':'Euro',
'DOV143b':'Euro',
'DOV144a':'Euro',
'DOV144b':'Euro',
'DOV145a':'Euro',
'DOV145b':'Euro',
'DOV147a':'Euro',
'DOV147b':'Euro',
'DOV148a':'Euro',
'DOV148b':'Euro',
'DOV149a':'Euro',
'DOV149b':'Euro',
}
    


    dd = nex_mo.data_dict_from_file('Euro.nex', pop_assignments)
    
    
    
    fs_onepop = dadi.Spectrum.from_data_dict(dd, ['Euro'], [37], polarized = False)
    import pylab
    pylab.figure(figsize = (8,6))
    dadi.Plotting.plot_1d_fs(fs_onepop)
    pylab.savefig("Dove.Euro.spectra.XX.projection.pdf", bbox_inches = 'tight')
    pylab.savefig("Dove.Euro.spectra.XX.projection.png", dpi = 300, bbox_inches = 'tight')



fs_onepop.tofile('Dove.Euro.spectra.XX.projection.fs')
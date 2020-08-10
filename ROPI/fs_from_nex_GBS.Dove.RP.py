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
'DOV164a':'RP',
'DOV164b':'RP',
'DOV150a':'RP',
'DOV150b':'RP',
'DOV155a':'RP',
'DOV155b':'RP',
'DOV156a':'RP',
'DOV156b':'RP',
'DOV157a':'RP',
'DOV157b':'RP',
'DOV158a':'RP',
'DOV158b':'RP',
'DOV159a':'RP',
'DOV159b':'RP',
'DOV160a':'RP',
'DOV160b':'RP',
'DOV161a':'RP',
'DOV161b':'RP',
'DOV162a':'RP',
'DOV162b':'RP',
'DOV163a':'RP',
'DOV163b':'RP',
'DOV165a':'RP',
'DOV165b':'RP',
'DOV166a':'RP',
'DOV166b':'RP',
'DOV167a':'RP',
'DOV167b':'RP',
'DOV168a':'RP',
'DOV168b':'RP',
'DOV169a':'RP',
'DOV169b':'RP',
'DOV170a':'RP',
'DOV170b':'RP',
'DOV171a':'RP',
'DOV171b':'RP',
'DOV172a':'RP',
'DOV172b':'RP',
'DOV173a':'RP',
'DOV173b':'RP',
'DOV174a':'RP',
'DOV174b':'RP',
'DOV175a':'RP',
'DOV175b':'RP',
'DOV176a':'RP',
'DOV176b':'RP',
'DOV177a':'RP',
'DOV177b':'RP',
'DOV178a':'RP',
'DOV178b':'RP',
'DOV179a':'RP',
'DOV179b':'RP',
'DOV180a':'RP',
'DOV180b':'RP',
'DOV181a':'RP',
'DOV181b':'RP',
'DOV182a':'RP',
'DOV182b':'RP',
'DOV183a':'RP',
'DOV183b':'RP',
'DOV184a':'RP',
'DOV184b':'RP',
'DOV185a':'RP',
'DOV185b':'RP',
'DOV186a':'RP',
'DOV186b':'RP'
}



    dd = nex_mo.data_dict_from_file('Dove.RP.nex', pop_assignments)
    
    
    
    fs_onepop = dadi.Spectrum.from_data_dict(dd, ['RP'], [28], polarized = False)
    import pylab
    pylab.figure(figsize = (8,6))
    dadi.Plotting.plot_1d_fs(fs_onepop)
    pylab.savefig("Dove.RP.spectra.XX.projection.pdf", dpi = 300, bbox_inches = 'tight')
    pylab.savefig("Dove.RP.spectra.XX.projection.png", dpi = 300, bbox_inches = 'tight')



fs_onepop.tofile('Dove.RP.spectra.XX.projection.fs')
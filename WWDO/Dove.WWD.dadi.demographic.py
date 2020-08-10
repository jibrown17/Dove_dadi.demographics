##################################################################################################
# For use in python v. 3.7
#
# Originally created by Joshua I. Brown in August, 2020
#
# Description: Stepwise demographic modelling for site frequency spectrum data
#
#
# NOTE: This script depends on frequency spectrum file created by fs_from_nex_GBS.py
##################################################################################################



import dadi
import numpy
import sys
from numpy import array
from dadi import Demographics2D, Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import pylab



data = dadi.Spectrum.from_file('Dove.WWD.spectra.59.projection.fs')
ns = data.sample_sizes
data.mask_corners()
pts_1 = [100, 110, 120]



def psmc_100steps(params, ns, pts_1):

    nu99, nu98, nu97, nu96, nu95, nu94, nu93, nu92, nu91, nu90, nu89, nu88, nu87, nu86, nu85, nu84, nu83, nu82, nu81, nu80, nu79, nu78, nu77, nu76, nu75, nu74, nu73, nu72, nu71, nu70, nu69, nu68, nu67, nu66, nu65, nu64, nu63, nu62, nu61, nu60, nu59, nu58, nu57, nu56, nu55, nu54, nu53, nu52, nu51, nu50, nu49, nu48, nu47, nu46, nu45, nu44, nu43, nu42, nu41, nu40, nu39, nu38, nu37, nu36, nu35, nu34, nu33, nu32, nu31, nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nuA, T99, T98, T97, T96, T95, T94, T93, T92, T91, T90, T89, T88, T87, T86, T85, T84, T83, T82, T81, T80, T79, T78, T77, T76, T75, T74, T73, T72, T71, T70, T69, T68, T67, T66, T65, T64, T63, T62, T61, T60, T59, T58, T57, T56, T55, T54, T53, T52, T51, T50, T49, T48, T47, T46, T45, T44, T43, T42, T41, T40, T39, T38, T37, T36, T35, T34, T33, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0 = params
    


    xx = Numerics.default_grid(pts_1)
    # initial phi with ancestral pop: (nuA = 1)
    phi = PhiManip.phi_1D(xx,nu=nuA)
    # stays at nuA for T0 duration of time:
    phi = Integration.one_pop(phi,xx,T0,nuA)
    # followed by a number of time steps, with associated pop changes:
    phi = Integration.one_pop(phi, xx, T1, nu1)
    phi = Integration.one_pop(phi, xx, T2, nu2)
    phi = Integration.one_pop(phi, xx, T3, nu3)
    phi = Integration.one_pop(phi, xx, T4, nu4)
    phi = Integration.one_pop(phi, xx, T5, nu5)
    phi = Integration.one_pop(phi, xx, T6, nu6)
    phi = Integration.one_pop(phi, xx, T7, nu7)
    phi = Integration.one_pop(phi, xx, T8, nu8)
    phi = Integration.one_pop(phi, xx, T9, nu9)
    phi = Integration.one_pop(phi, xx, T10, nu10)
    phi = Integration.one_pop(phi, xx, T11, nu11)
    phi = Integration.one_pop(phi, xx, T12, nu12)
    phi = Integration.one_pop(phi, xx, T13, nu13)
    phi = Integration.one_pop(phi, xx, T14, nu14)
    phi = Integration.one_pop(phi, xx, T15, nu15)
    phi = Integration.one_pop(phi, xx, T16, nu16)
    phi = Integration.one_pop(phi, xx, T17, nu17)
    phi = Integration.one_pop(phi, xx, T18, nu18)
    phi = Integration.one_pop(phi, xx, T19, nu19)
    phi = Integration.one_pop(phi, xx, T20, nu20)
    phi = Integration.one_pop(phi, xx, T21, nu21)
    phi = Integration.one_pop(phi, xx, T22, nu22)
    phi = Integration.one_pop(phi, xx, T23, nu23)
    phi = Integration.one_pop(phi, xx, T24, nu24)
    phi = Integration.one_pop(phi, xx, T25, nu25)
    phi = Integration.one_pop(phi, xx, T26, nu26)
    phi = Integration.one_pop(phi, xx, T27, nu27)
    phi = Integration.one_pop(phi, xx, T28, nu28)
    phi = Integration.one_pop(phi, xx, T29, nu29)
    phi = Integration.one_pop(phi, xx, T30, nu30)
    phi = Integration.one_pop(phi, xx, T31, nu31)
    phi = Integration.one_pop(phi, xx, T32, nu32)
    phi = Integration.one_pop(phi, xx, T33, nu33)
    phi = Integration.one_pop(phi, xx, T34, nu34)
    phi = Integration.one_pop(phi, xx, T35, nu35)
    phi = Integration.one_pop(phi, xx, T36, nu36)
    phi = Integration.one_pop(phi, xx, T37, nu37)
    phi = Integration.one_pop(phi, xx, T38, nu38)
    phi = Integration.one_pop(phi, xx, T39, nu39)
    phi = Integration.one_pop(phi, xx, T40, nu40)
    phi = Integration.one_pop(phi, xx, T41, nu41)
    phi = Integration.one_pop(phi, xx, T42, nu42)
    phi = Integration.one_pop(phi, xx, T43, nu43)
    phi = Integration.one_pop(phi, xx, T44, nu44)
    phi = Integration.one_pop(phi, xx, T45, nu45)
    phi = Integration.one_pop(phi, xx, T46, nu46)
    phi = Integration.one_pop(phi, xx, T47, nu47)
    phi = Integration.one_pop(phi, xx, T48, nu48)
    phi = Integration.one_pop(phi, xx, T49, nu49)
    phi = Integration.one_pop(phi, xx, T50, nu50)
    phi = Integration.one_pop(phi, xx, T51, nu51)
    phi = Integration.one_pop(phi, xx, T52, nu52)
    phi = Integration.one_pop(phi, xx, T53, nu53)
    phi = Integration.one_pop(phi, xx, T54, nu54)
    phi = Integration.one_pop(phi, xx, T55, nu55)
    phi = Integration.one_pop(phi, xx, T56, nu56)
    phi = Integration.one_pop(phi, xx, T57, nu57)
    phi = Integration.one_pop(phi, xx, T58, nu58)
    phi = Integration.one_pop(phi, xx, T59, nu59)
    phi = Integration.one_pop(phi, xx, T60, nu60)
    phi = Integration.one_pop(phi, xx, T61, nu61)
    phi = Integration.one_pop(phi, xx, T62, nu62)
    phi = Integration.one_pop(phi, xx, T63, nu63)
    phi = Integration.one_pop(phi, xx, T64, nu64)
    phi = Integration.one_pop(phi, xx, T65, nu65)
    phi = Integration.one_pop(phi, xx, T66, nu66)
    phi = Integration.one_pop(phi, xx, T67, nu67)
    phi = Integration.one_pop(phi, xx, T68, nu68)
    phi = Integration.one_pop(phi, xx, T69, nu69)
    phi = Integration.one_pop(phi, xx, T70, nu70)
    phi = Integration.one_pop(phi, xx, T71, nu71)
    phi = Integration.one_pop(phi, xx, T72, nu72)
    phi = Integration.one_pop(phi, xx, T73, nu73)
    phi = Integration.one_pop(phi, xx, T74, nu74)
    phi = Integration.one_pop(phi, xx, T75, nu75)
    phi = Integration.one_pop(phi, xx, T76, nu76)
    phi = Integration.one_pop(phi, xx, T77, nu77)
    phi = Integration.one_pop(phi, xx, T78, nu78)
    phi = Integration.one_pop(phi, xx, T79, nu79)
    phi = Integration.one_pop(phi, xx, T80, nu80)
    phi = Integration.one_pop(phi, xx, T81, nu81)
    phi = Integration.one_pop(phi, xx, T82, nu82)
    phi = Integration.one_pop(phi, xx, T83, nu83)
    phi = Integration.one_pop(phi, xx, T84, nu84)
    phi = Integration.one_pop(phi, xx, T85, nu85)
    phi = Integration.one_pop(phi, xx, T86, nu86)
    phi = Integration.one_pop(phi, xx, T87, nu87)
    phi = Integration.one_pop(phi, xx, T88, nu88)
    phi = Integration.one_pop(phi, xx, T89, nu89)
    phi = Integration.one_pop(phi, xx, T90, nu90)
    phi = Integration.one_pop(phi, xx, T91, nu91)
    phi = Integration.one_pop(phi, xx, T92, nu92)
    phi = Integration.one_pop(phi, xx, T93, nu93)
    phi = Integration.one_pop(phi, xx, T94, nu94)
    phi = Integration.one_pop(phi, xx, T95, nu95)
    phi = Integration.one_pop(phi, xx, T96, nu96)
    phi = Integration.one_pop(phi, xx, T97, nu97)
    phi = Integration.one_pop(phi, xx, T98, nu98)
    phi = Integration.one_pop(phi, xx, T99, nu99)
    
    # get sfs:
    fs = Spectrum.from_phi(phi,ns,(xx,))
    return fs

func = psmc_100steps

params = array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017])

upper_bound = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15]
lower_bound = [1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1.00E-05, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07]


reps = 50
while reps > 0:

    func_ex = dadi.Numerics.make_extrap_log_func(func)
    p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=lower_bound, upper_bound=upper_bound)
    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_1, 
    lower_bound=lower_bound, upper_bound=upper_bound, verbose= 1)
    
    model = func_ex(popt, ns, pts_1)
    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    ll_opt = dadi.Inference.ll_multinom(model, data)


    outfile = "Output.Figures\Dove.WWD.demographic.resid."
    scaled_model = dadi.Inference.optimally_scaled_sfs(model, data)
    pylab.figure(figsize = (8,6))
    dadi.Plotting.plot_1d_comp_Poisson(scaled_model, data, fig_num = None, residual = 'Anscombe', plot_masked = False, show = False)
    pylab.savefig(outfile + str(reps) + ".pdf", dpi = 300, bbox_inches = 'tight')
    pylab.savefig(outfile + str(reps) + ".png", dpi = 300, bbox_inches = 'tight')

    print ('Optimized parameters', repr([theta, ll_opt, popt]))
    file = open("Dove.WWD.demographic.3mil.txt", "a")
    file.write("\nReplicate" + str(reps) + repr([theta, ll_opt, popt]))
    file.close()
    reps = reps - 1

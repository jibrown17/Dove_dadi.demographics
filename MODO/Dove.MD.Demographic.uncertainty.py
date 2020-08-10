##################################################################################################
# For use in python v. 3.7
#
# Originally created by Joshua I. Brown in August, 2020
#
# Description: FIM Godambe uncertainty analysis for Stepwise demographic modelling from site frequency spectrum data
#
#
# NOTE1: This script depends on frequency spectrum file created by fs_from_nex_GBS.py
# as well as the optimum parameters output from psmc_100steps model
#
# NOTE2: popt should be filled with 100 optimum parameter output from psmc_100steps modelling
#
# Note3: eps setting might need to be adjusted to maximize number of parameters returned
##################################################################################################



import dadi
from dadi import Misc
from numpy import array
from dadi import Demographics2D, Numerics, PhiManip, Integration, Spectrum, Plotting, Inference, Godambe
import pylab
import os

data = dadi.Spectrum.from_file('..\Dove.MD.spectra.61.projection.fs')
ns = data.sample_sizes
data.mask_corners()
pts_l = [70, 80, 90]


popt = array([2.96815611, 8.08181663, 1.00794642, 2.4937849, 1.3024127, 1.2521957, 0.96844368, 0.76829648, 0.044929, 0.93051195, 1.08403878, 0.35000964, 0.80929236, 1.20800958, 0.5164057, 0.51478456, 0.34672341, 0.26208969, 1.05350482, 0.86885817, 0.4906577, 1.07977762, 0.561707, 0.81686422, 0.73641625, 0.58154166, 1.86317774, 0.56834889, 1.78657548, 0.76110536, 1.75061795, 1.00363856, 0.3607688, 1.10117532, 0.93423354, 0.44015494, 0.46035059, 0.9069816, 0.67699126, 1.11167357, 0.5431997, 0.72601026, 1.25217817, 1.65263253, 0.61419044, 1.76335727, 1.77779789, 1.07524213, 0.47834043, 1.63455055, 1.67643482, 0.58881888, 1.27038124, 0.53932669, 0.52035346, 1.32332054, 1.66758414, 0.67216648, 1.01905396, 1.11516596, 0.48928674, 0.56886739, 1.32611131, 0.58464586, 0.64235066, 0.83221323, 1.06989439, 1.30926319, 0.76289588, 0.86336025, 1.21118653, 0.95357196, 1.01738049, 1.74928153, 0.55164122, 1.96324081, 0.78157661, 0.78255532, 1.12922869, 0.90628536, 0.70982824, 1.20085058, 0.78110149, 0.57716117, 0.53213967, 1.12608463, 0.61398131, 0.72779294, 0.77223211, 0.55077234, 0.68856367, 0.84598914, 1.39794878, 1.92286073, 0.78409454, 0.70244123, 1.67346127, 0.66771644, 1.43678233, 1.26537844, 0.08222437, 0.06654138, 0.06745076, 0.075441, 0.04737331, 0.08032472, 0.034285, 0.03026729, 0.13437665, 0.03931192, 0.03058096, 0.04070508, 0.03581052, 0.03645803, 0.07824934, 0.05578025, 0.10950835, 0.09973036, 0.05085921, 0.04775462, 0.03761175, 0.05736868, 0.03039343, 0.0616437, 0.03216039, 0.05085589, 0.03465843, 0.02883101, 0.03210428, 0.09206346, 0.07219586, 0.06039391, 0.10170144, 0.0364304, 0.06622764, 0.05243868, 0.09195893, 0.08431568, 0.02860828, 0.04088037, 0.05260896, 0.02591021, 0.06271884, 0.02856289, 0.03561092, 0.06630729, 0.0271773, 0.02901717, 0.0706367, 0.0259205, 0.02591877, 0.0303138, 0.04406768, 0.03078876, 0.09688871, 0.09034556, 0.0602777, 0.02608832, 0.07463813, 0.02621465, 0.07357389, 0.04595958, 0.02969308, 0.07903388, 0.03250043, 0.04595167, 0.05979371, 0.09602857, 0.0412124, 0.09257762, 0.03410599, 0.02830409, 0.06674192, 0.05496302, 0.06081634, 0.07150405, 0.04280067, 0.07887487, 0.02982208, 0.08846264, 0.03281631, 0.09496371, 0.03563491, 0.03234779, 0.02744988, 0.05862141, 0.02663753, 0.06653985, 0.03276419, 0.05431482, 0.03144328, 0.02520009, 0.06310801, 0.08975984, 0.08404353, 0.08042964, 0.06884432, 0.04564625, 0.08471321, 0.081185])


def psmc_100steps(params, ns, pts_l):

    nu99, nu98, nu97, nu96, nu95, nu94, nu93, nu92, nu91, nu90, nu89, nu88, nu87, nu86, nu85, nu84, nu83, nu82, nu81, nu80, nu79, nu78, nu77, nu76, nu75, nu74, nu73, nu72, nu71, nu70, nu69, nu68, nu67, nu66, nu65, nu64, nu63, nu62, nu61, nu60, nu59, nu58, nu57, nu56, nu55, nu54, nu53, nu52, nu51, nu50, nu49, nu48, nu47, nu46, nu45, nu44, nu43, nu42, nu41, nu40, nu39, nu38, nu37, nu36, nu35, nu34, nu33, nu32, nu31, nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nuA, T99, T98, T97, T96, T95, T94, T93, T92, T91, T90, T89, T88, T87, T86, T85, T84, T83, T82, T81, T80, T79, T78, T77, T76, T75, T74, T73, T72, T71, T70, T69, T68, T67, T66, T65, T64, T63, T62, T61, T60, T59, T58, T57, T56, T55, T54, T53, T52, T51, T50, T49, T48, T47, T46, T45, T44, T43, T42, T41, T40, T39, T38, T37, T36, T35, T34, T33, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0 = params
    


    xx = Numerics.default_grid(pts_l)
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
func_ex = dadi.Numerics.make_extrap_func(func)
model = func_ex(popt, ns, pts_l)


uncerts_fim = dadi.Godambe.FIM_uncert(func_ex, pts_l, popt, data, multinom=True, eps = 0.01)
file = open("Dove.MD.Deomgraphic.uncertaintymultinom.false.txt", "a")
file.write('\nEstimated parameter standard deviations from FIM 0.01: {0}'.format(uncerts_fim))
file.close()

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

data = dadi.Spectrum.from_file('..\Dove.RP.spectra.28.projection.fs')
ns = data.sample_sizes
data.mask_corners()
pts_l = [50, 60, 70]


popt = array([3.480855789, 2.873985657, 1.585800637, 1.858710614, 1.894249512, 2.315294579, 2.609922826, 3.245398423, 3.413463759, 3.715352479, 3.625198640, 4.086547009, 3.803064358, 4.001388391, 3.585272820, 2.853539608, 2.799530220, 2.401181599, 2.091806792, 1.95158604, 1.778461611, 1.635481615, 1.425818897, 1.321687683, 1.278282085, 1.232750805, 1.124840660, 1.053786238, 1.034159862, 0.986832459, 0.918969986, 0.777366910, 0.959394937, 0.939062678, 0.756239790, 0.858196627, 0.901421971, 0.811705389, 0.808305198, 0.881853419, 0.871947159, 0.705385579, 0.759343202, 0.885325443, 0.88892782, 0.869804533, 0.801306285, 0.820562577, 0.799215757, 0.887476275, 0.946734341, 0.822997645, 0.833190482, 0.899113395, 0.878605579, 0.76923142, 0.895112922, 0.751909312, 0.879698855, 0.943888312, 0.947928887, 0.910826452, 0.953940724, 0.955806673, 0.834342025, 0.867258069, 0.895959569, 0.922253542, 0.888079399, 0.910396245, 1.041646096, 0.864881616, 0.90710853, 0.918553579, 0.987214041, 0.941968783, 1.00714442, 0.996576929, 0.998896018, 1.121509978, 0.998941537, 0.941017286, 1.086337861, 1.005134866, 1.13311259, 0.954317996, 1.076697474, 1.068795191, 1.107989851, 1.041755763, 1.057945023, 1.196259966, 1.018041752, 1.052488244, 1.061760806, 1.097907906, 0.987215271, 1.051848006, 1.178997638, 0.294974415, 0.082587607, 0.034244271, 0.035595365, 0.031914651, 0.041717962, 0.040612282, 0.033738677, 0.034084432, 0.033199467, 0.031966217, 0.030864233, 0.029532416, 0.034087892, 0.034462781, 0.035463901, 0.024884478, 0.028735612, 0.025058619, 0.02506213, 0.022351306, 0.023022337, 0.021157847, 0.02030686, 0.020868394, 0.022555123, 0.020100697, 0.021099062, 0.021993316, 0.024724142, 0.019453262, 0.025536902, 0.026291816, 0.021025785, 0.023939057, 0.023768654, 0.019980565, 0.020598373, 0.024671614, 0.022208298, 0.022200809, 0.024387557, 0.023635965, 0.023013962, 0.02285928, 0.018834981, 0.020873966, 0.022958784, 0.024660908, 0.022537121, 0.021108005, 0.019831579, 0.01870309, 0.021840155, 0.020837561, 0.020675153, 0.022471133, 0.021837082, 0.021878558, 0.023405976, 0.02150501, 0.019615387, 0.020117187, 0.019126801, 0.018613185, 0.021680591, 0.022520529, 0.020733733, 0.019634038, 0.020564262, 0.02183455, 0.020282268, 0.021473674, 0.020002258, 0.019835366, 0.020836005, 0.019343708, 0.019788815, 0.020367828, 0.018866855, 0.019038347, 0.018067178, 0.019536308, 0.021121245, 0.02084206, 0.018370318, 0.019625015, 0.020417508, 0.020963028, 0.019391359, 0.019186811, 0.019983053, 0.019353382, 0.020643593, 0.021596631, 0.018576662, 0.020526466, 0.021457288, 0.020034738, 0.019579912, 0.020987947])


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
file = open("Dove.RP.Deomgraphic.uncertaintymultinom.false.txt", "a")
file.write('\nEstimated parameter standard deviations from FIM 0.01: {0}'.format(uncerts_fim))
file.close()

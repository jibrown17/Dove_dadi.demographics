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

data = dadi.Spectrum.from_file('Dove.Euro.spectra.37.projection.fs')
ns = data.sample_sizes
data.mask_corners()
pts_l = [50, 60, 70]


popt = array([1.40725005, 3.441410089, 3.290958281, 4.317148681, 4.315404266, 2.974332419, 2.511498343, 1.833090138, 1.589364203, 1.267292661, 0.96457261, 0.845428673, 0.77915858, 0.840041285, 0.7939469, 0.87987163, 0.924603875, 0.874404477, 0.969677705, 0.899909181, 0.872119562, 0.902562742, 1.020593663, 0.89684564, 0.88317703, 0.988100252, 1.116856614, 0.920791124, 0.871013791, 0.940709394, 1.031743994, 1.063494752, 1.041452899, 1.004997424, 0.938657229, 0.962137804, 1.041318701, 1.098895296, 1.050044155, 1.017258428, 1.056212013, 1.110897757, 1.0336126, 0.957308242, 1.092604363, 0.985214815, 1.019502186, 1.120661098, 1.084576087, 1.11527676, 1.041497059, 0.974914328, 1.062685326, 0.986729885, 1.03824969, 0.920432198, 1.182482942, 1.128851737, 1.079427624, 1.147005264, 1.034426405, 1.108100867, 1.120954014, 1.100577574, 1.113562809, 1.082749609, 1.10336734, 1.097610901, 0.963327184, 1.043893095, 1.054458057, 1.119908893, 1.033906928, 1.1251969, 0.971295578, 1.0409, 1.020225874, 1.080408596, 1.084000488, 1.140459721, 1.098151733, 1.042359642, 1.089591514, 1.050127366, 1.119862127, 1.101128757, 1.038379195, 0.997112921, 1.106606974, 1.003773498, 1.05964012, 1.052268829, 1.101078205, 1.020978883, 1.239615388, 0.966163555, 0.981634133, 1.082593337, 1.104689548, 1.060831872, 0.051150828, 0.067349241, 0.045940244, 0.048149175, 0.044635925, 0.028571658, 0.025400531, 0.02081164, 0.021376624, 0.022547615, 0.022324301, 0.024802798, 0.025027285, 0.022691996, 0.02238314, 0.023092704, 0.024447785, 0.02341983, 0.024845447, 0.023864991, 0.026499867, 0.024129836, 0.025493044, 0.026899597, 0.027681527, 0.025455379, 0.024850386, 0.024366498, 0.024596468, 0.027289406, 0.027147659, 0.024869259, 0.025014425, 0.02512804, 0.02417242, 0.027193783, 0.028745963, 0.025200196, 0.027648513, 0.028188523, 0.024712289, 0.026917022, 0.023511976, 0.024371179, 0.024638771, 0.026136917, 0.025913945, 0.025368969, 0.026373592, 0.024855381, 0.02836445, 0.025850902, 0.024739647, 0.025626913, 0.024614722, 0.027358872, 0.024853578, 0.027023172, 0.024681063, 0.024000603, 0.024090255, 0.024213303, 0.022075505, 0.02580331, 0.0260304, 0.029063004, 0.026062429, 0.025728925, 0.022400648, 0.025198499, 0.023411164, 0.024267576, 0.023316789, 0.025664002, 0.02738461, 0.027501974, 0.027621826, 0.027562613, 0.025453833, 0.025963885, 0.024644237, 0.025603974, 0.023563291, 0.02746057, 0.028258646, 0.023270464, 0.026458981, 0.028572472, 0.029077404, 0.023560791, 0.023838519, 0.024547119, 0.025300496, 0.025647043, 0.025613164, 0.025905781, 0.028053009, 0.026132659, 0.025833187, 0.027671429])


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
file = open("Dove.Euro.Deomgraphic.uncertaintymultinom.false.txt", "a")
file.write('\nEstimated parameter standard deviations from FIM 0.01: {0}'.format(uncerts_fim))
file.close()



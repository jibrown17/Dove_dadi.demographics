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

data = dadi.Spectrum.from_file('..\Dove.WWD.spectra.59.projection.fs')
ns = data.sample_sizes
data.mask_corners()
pts_l = [70, 80, 90]


popt = array([3.062671307, 3.014618746, 3.544226831, 3.058573124, 2.592078287, 2.17387838, 1.897196882, 1.652280528, 1.482993668, 1.383487423, 1.342953369, 1.282623093, 1.239798862, 1.114274013, 1.085777238, 1.006652794, 1.033032657, 1.070897154, 1.017461495, 1.032370551, 1.005156273, 0.983674287, 1.016761207, 0.981416856, 0.985864438, 1.035302269, 0.946156259, 0.86352937, 1.129731171, 0.966877714, 1.021034884, 1.002107189, 1.003111686, 0.988780401, 0.991100421, 1.074433056, 0.998628904, 0.960458004, 0.956958748, 0.980868002, 0.985583106, 0.970506007, 1.017239443, 1.061238634, 0.92373299, 1.024701845, 1.043181321, 0.973531793, 0.945075506, 1.04690614, 1.079650078, 0.978008538, 0.951855993, 1.036883519, 1.006447582, 0.91005365, 1.039219448, 1.040266099, 1.047538231, 0.975648969, 0.976507284, 1.063194457, 0.963577854, 0.892525829, 1.011424403, 1.007153483, 1.083877823, 1.092931415, 0.957935694, 1.043119036, 1.07028242, 0.984483865, 0.986167978, 0.912333549, 0.994604985, 1.029860677, 1.145072564, 0.973761149, 1.173890557, 1.096125621, 1.087456604, 0.968474168, 1.090984804, 0.973052938, 1.042112133, 0.961472148, 0.991460846, 1.011038036, 1.052005606, 1.104403966, 0.998003414, 1.103776756, 1.064449382, 1.081548759, 1.037432383, 1.009406054, 1.162944308, 1.089626518, 1.069633923, 0.542680154, 0.074631347, 0.067005377, 0.047263635, 0.030506046, 0.023371623, 0.022567091, 0.019792598, 0.019329389, 0.018913599, 0.019426437, 0.01785238, 0.019694926, 0.019214001, 0.018070218, 0.019611052, 0.020245831, 0.017988788, 0.018917813, 0.019059411, 0.018046711, 0.019113406, 0.017766285, 0.018469798, 0.018073254, 0.017531784, 0.016975488, 0.019472311, 0.018698715, 0.017348351, 0.018094732, 0.019825938, 0.017736828, 0.018678977, 0.017330822, 0.017369106, 0.017054071, 0.01711767, 0.018938056, 0.016756089, 0.017421785, 0.01732535, 0.01845543, 0.017690316, 0.01764808, 0.016355234, 0.020681393, 0.017892613, 0.018057709, 0.016255036, 0.017261433, 0.017399067, 0.017834089, 0.018862131, 0.017501265, 0.019203162, 0.018448087, 0.017356349, 0.018369111, 0.017300215, 0.018397856, 0.018291244, 0.017382372, 0.017949732, 0.018430406, 0.018405468, 0.018011955, 0.019113156, 0.018658146, 0.019194469, 0.016814043, 0.019644062, 0.018957308, 0.017255077, 0.016920141, 0.019632765, 0.016269025, 0.017059592, 0.01615029, 0.016587594, 0.016377016, 0.017211369, 0.018049531, 0.018071227, 0.018457596, 0.016783489, 0.017362155, 0.016878927, 0.019593636, 0.018065313, 0.016317725, 0.017858356, 0.018558923, 0.018743943, 0.019681532, 0.016033705, 0.016918051, 0.017452099, 0.018441493, 0.016696674, 0.019599972])


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
file = open("Dove.WWD.Deomgraphic.uncertaintymultinom.false.txt", "a")
file.write('\nEstimated parameter standard deviations from FIM 0.01: {0}'.format(uncerts_fim))
file.close()


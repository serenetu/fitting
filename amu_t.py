__author__ = 'SereneTu'
import numpy as np
# import sympy as sympy
from scipy.optimize import minimize
# from sympy import *
from scipy import integrate
from scipy.integrate import quadrature
from scipy.integrate import fixed_quad
import matplotlib.pyplot as plt
from scipy.integrate import quad
import random
import copy
import math
import sys
import src.func
import HISQ_DDS as cmunu_src
# import fit_from_cmatrix as amu_src
import functools


def log(function):
    @functools.wraps(function)
    def wrapper(self, *args, **kwargs):
        print 'Call %s::%s():' % (self.__class__.__name__, function.__name__)
        return function(self, *args, **kwargs)
    return wrapper


def make_jackknife_dic(dic):

    keys = dic.keys()
    jackknife_dic = {}
    for i in range(len(keys)):
        array = []
        new_configs = keys[:i] + keys[i+1:]

        for config in new_configs:
            array.append(dic[config])
        array = np.array(array)
        mean = np.mean(array, axis=0)
        jackknife_dic[keys[i]] = mean

    return jackknife_dic


def jackknife_avg_err_dic(jackknife_dic):
    mean, err = jackknife_avg_err_array(jackknife_dic.values())
    return mean, err


def jackknife_avg_err_array(jackknife_array):
    array = np.array(jackknife_array)
    mean = np.mean(array, axis=0)
    err = (len(jackknife_array) - 1.) ** (1. / 2.) * np.std(array, axis=0)
    return mean, err


class AmuT:

    def __init__(self):
        self.alpha
        self.ainv
        self.mu_lat
        self.fac
        self.T_MAX
        return

    def z_q2(self):
        mu = self.mu_lat
        return lambda q2: ((q2 ** 2. + 4. * mu ** 2. * q2) ** (1./2.) - q2) / (2. * mu ** 2. * q2)

    def k_q2(self):
        mu = self.mu_lat
        z = self.z_q2()
        return lambda q2: mu ** 2. * q2 * (z(q2)) ** 3. * (1. - q2 * z(q2)) / (1. + mu ** 2. * q2 * (z(q2)) ** 2.)

    def plt_k_q2(self):
        k = self.k_q2()
        x = np.linspace(0, 1, 50)
        plt.plot(x, k(x))
        return

    def pi_q2(self, cmunu):

        def func(q2):
            res = 0.
            for t in range(0, self.T_MAX / 2):
                tf = float(t)
                res_t = ((math.cos(q2 ** (1. / 2.) * tf) - 1.) / (2. * math.sin(q2 ** (1. / 2.) / 2.)) ** 2. + 1. / 2. * tf ** 2.) * cmunu[t]
                res += res_t
            for t in range(self.T_MAX / 2, self.T_MAX):
                tf = float(t)
                res_t = ((math.cos(q2 ** (1. / 2.) * (tf - self.T_MAX)) - 1.) / (2. * math.sin(q2 ** (1. / 2.) / 2.)) ** 2. + 1. / 2. * (tf - self.T_MAX) ** 2.) * cmunu[t]
                res += res_t
            return res

        return lambda q2: func(q2)

    def plt_pi_q2(self, cmunu):
        q2_array = [(2.*np.pi*i/self.T_MAX) ** 2. for i in range(1, self.T_MAX)]
        q2_array = np.array(q2_array)
        pi_q2 = self.pi_q2(cmunu)
        plt.plot(q2_array[:30] * self.ainv ** 2., pi_q2(q2_array)[:30], '.')
        return

    def func_f_q2(self, t):
        tf = float(t)
        return lambda q2: (math.cos(q2 ** (1. / 2.) * tf) - 1.) / (2. * math.sin(q2 ** (1. / 2.) / 2.)) ** 2. + 1. / 2. * tf ** 2.

    def get_window_t(self, t, t0=0.4, t1=1.0, delta=0.3): # t0 and t1 are in unit of fm
        gevinv_to_fm = 0.1975 # 1 GeV-1 = 0.1975 fm
        res = (1. + np.tanh((t * self.a * gevinv_to_fm - t0) / delta)) / 2.
        res -= (1. + np.tanh((t * self.a * gevinv_to_fm - t1) / delta)) / 2.
        return res

    def get_window_t_list(self, t0=0.4, t1=1.0, delta=0.3):
        self.window_tlist = []
        for tt in range(self.T_MAX / 2 + 1):
            self.window_tlist.append(self.get_window_t(tt, t0, t1, delta))
        return self.window_tlist

    def update_w_t_list_with_window_t_list(self):
        try:
            self.window_tlist
        except NameError:
            print 'Run Amu_T::get_window_t_list() First'
        else:
            for tt in range(self.T_MAX / 2 + 1):
                self.w_tlist[tt] = self.w_tlist[tt] * self.window_tlist[tt]
        return

    def get_w_t(self, t):
        tf = float(t)
        k = self.k_q2()
        f = self.func_f_q2(tf)
        # kf = lambda q2: k(q2) * f(q2)
        kf = lambda q2: k(q2) * f(q2)
        res = self.fac * np.array(quad(kf, 0., np.inf))
        return res[0]

    def get_w_t_list(self):
        self.w_tlist = []
        for tt in range(self.T_MAX / 2 + 1):
            self.w_tlist.append(self.get_w_t(tt))
        return self.w_tlist

    def get_wc_t_list(self, cmunu):
        try:
            self.w_tlist
        except NameError:
            print 'Run Amu_T::get_w_tlist() First'
        else:
            res = []
            fcmunu = self.fold_cmunu(cmunu)
            # for tt in range(self.T_MAX / 2 + 1):
            for tt in range(len(fcmunu)):
                res.append(self.w_tlist[tt] * fcmunu[tt])
            return np.array(res)
        return

    def get_wc_t_list_jackknife(self, cmunu_dic):
        wt_ct_dic = {}
        for config in cmunu_dic:
            wt_ct_dic[config] = self.get_wc_t_list(cmunu_dic[config])

        self.wt_ct_jackknife_dic = make_jackknife_dic(wt_ct_dic)
        return self.wt_ct_jackknife_dic

    def get_wc_t_list_jackknife_avg_err(self):
        try:
            self.wt_ct_jackknife_dic
        except NameError:
            print 'Run Amu_T::wt_ct_jackknife() First'
        else:
            mean, err = jackknife_avg_err_dic(self.wt_ct_jackknife_dic)
            return mean, err

    def fold_cmunu(self, cmunu):
        cmunu_fold = []
        len_fold = len(cmunu) / 2
        cmunu_fold.append(cmunu[0])
        for i in range(1, len_fold):
            cmunu_fold.append((cmunu[i]+cmunu[-1-i+1])/2.)
        cmunu_fold.append(cmunu[len_fold+1])
        '''
        cmunu_fold = []
        len_fold = len(cmunu) / 2
        for i in range(1, len_fold):
            cmunu_fold.append((cmunu[i]+cmunu[-1-i])/2.)
        cmunu_fold.append(cmunu[len_fold+1])
        '''
        return np.array(cmunu_fold)

    def amu_tcut_list(self, cmunu):
        wt_ct = self.get_wc_t_list(cmunu)
        amu_tcut = []
        for tcut in range(self.T_MAX / 2 + 1):
            amu_tcut.append(np.sum(wt_ct[:tcut]))
        return np.array(amu_tcut)

    def get_amu_tcut_list_jackknife(self, cmunu_dic):
        amu_tcut_dic = {}
        for config in cmunu_dic:
            amu_tcut_dic[config] = self.amu_tcut_list(cmunu_dic[config])

        self.amu_tcut_list_jackknife_dic = make_jackknife_dic(amu_tcut_dic)
        return self.amu_tcut_list_jackknife_dic

    def get_amu_tcut_list_jackknife_avg_err(self):
        try:
            self.amu_tcut_list_jackknife_dic
        except NameError:
            print 'Run Amu_T::amu_tcut_list_jackknife() First'
        else:
            mean, err = jackknife_avg_err_dic(self.amu_tcut_list_jackknife_dic)
            return mean, err

    @log
    def plt_wc_t_list_jackknife(self, label =''):
        # t_list = range(self.T_MAX / 2 + 1)
        wc_t, wc_t_err = self.get_wc_t_list_jackknife_avg_err()
        t_list = range(len(wc_t))
        plt.errorbar(t_list, wc_t, wc_t_err, fmt='.', elinewidth=1, label = label)
        return

    @log
    def plt_amu_tcut_list_jackknife(self, label=''):
        # t_list = range(self.T_MAX / 2 + 1)
        amu, amu_err = self.get_amu_tcut_list_jackknife_avg_err()
        t_list = range(len(amu))
        plt.errorbar(t_list, amu, amu_err, fmt='.', elinewidth=1, label = label)
        return


class AmuT4864(AmuT):

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        #self.ainv = 1.629278350515464  # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.a = 0.12224 / 0.1975      # GeV-1  # From Email
        self.ainv = 1. / self.a
        self.mu_lat = 0.1056583715 * self.a
        self.fac = 4. * self.alpha ** 2. * (5. / 9.) * 10. ** 10. * 2.
        self.T_MAX = 64
        return


class AmuT6496(AmuT):

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        self.a = 1. / 2.255453                # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.a = 0.08786 / 0.1975             # GeV-1   # From Email
        self.ainv = 1. / self.a
        self.mu_lat = 0.1056583715 * self.a
        self.fac = 4. * self.alpha ** 2. * (5./9.) * 10. ** 10. * 2.
        self.T_MAX = 96
        return


class AmuT96192(AmuT):

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        self.a = 1. / 2.255453                # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.a = 0.05662 / 0.1975             # GeV-1 (1 GeV-1 = 0.1975 fm) # From Email
        self.ainv = 1. / self.a
        self.mu_lat = 0.1056583715 * self.a
        self.fac = 4. * self.alpha ** 2. * (5./9.) * 10. ** 10. * 2.
        self.T_MAX = 192
        return


if __name__ == "__main__":

    # =================================
    # l96192 amu_t
    # =================================

    # Read Cmunu
    CMUNU_PATH = '/Users/tucheng/Desktop/Physics/research/hvp/analysis/data/HISQ_cmunu/l96192/'
    CMUNU_EXACT_PATH = CMUNU_PATH + 'l96_Exact_cmunu'
    CMUNU_SUB_PATH = CMUNU_PATH + 'l96_Sub_cmunu'
    CMUNU_AMA_PATH = CMUNU_PATH + 'l96_AMA_cmunu'
    CMUNU_LMA_PATH = CMUNU_PATH + 'l96_LMA_cmunu'
    CMUNU_LMASUB_PATH = CMUNU_PATH + 'l96_LMASUB_cmunu'

    readcmunu = cmunu_src.Do96192()

    readcmunu.E_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_EXACT_PATH)
    readcmunu.ES_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_SUB_PATH)
    readcmunu.A_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_AMA_PATH)
    readcmunu.L_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_LMA_PATH)
    readcmunu.LS_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_LMASUB_PATH)

    # Subtraction
    CONFIG_LIST_E = readcmunu.E_cmunu_config_dic.keys()[:]
    CONFIG_LIST_E.sort()
    CONFIG_LIST_E = CONFIG_LIST_E[:]
    print 'CONFIG_LIST E:', len(CONFIG_LIST_E), CONFIG_LIST_E
    CONFIG_LIST_ES = readcmunu.ES_cmunu_config_dic.keys()[:]
    CONFIG_LIST_ES.sort()
    CONFIG_LIST_ES = CONFIG_LIST_ES[:]
    print 'CONFIG_LIST ES:', len(CONFIG_LIST_ES), CONFIG_LIST_ES
    CONFIG_LIST_A = readcmunu.A_cmunu_config_dic.keys()[:]
    CONFIG_LIST_A.sort()
    CONFIG_LIST_A = CONFIG_LIST_A[:]
    print 'CONFIG_LIST A:', len(CONFIG_LIST_A), CONFIG_LIST_A
    CONFIG_LIST_L = readcmunu.L_cmunu_config_dic.keys()[:]
    CONFIG_LIST_L.sort()
    CONFIG_LIST_L = CONFIG_LIST_L[:]
    print 'CONFIG_LIST L:', len(CONFIG_LIST_L), CONFIG_LIST_L
    CONFIG_LIST_LS = readcmunu.LS_cmunu_config_dic.keys()[:]
    CONFIG_LIST_LS.sort()
    CONFIG_LIST_LS = CONFIG_LIST_LS[:]
    print 'CONFIG_LIST LS:', len(CONFIG_LIST_LS), CONFIG_LIST_LS

    CONFIG_LIST_FOR_CONTRACT = [504, 648, 696, 744, 792, 840, 888, 936, 984, 1032, 1128, 1176, 1224, 1272, 1368, 1512, 1560, 1608, 1656, 1704, 1752, 1800]

    readcmunu.cmunu_subtract(CONFIG_LIST_FOR_CONTRACT, 'Exact Sub AMA LMASUB LMA')

    # plot some cmunu
    for config in CONFIG_LIST_FOR_CONTRACT:
        readcmunu.plot_cmumu_t(readcmunu.cmunu_config_subtract_dic[config])
    plt.show()

    # l96 wc_t
    amu = AmuT96192()
    amu.get_w_t_list()
    amu.get_wc_t_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu.plt_wc_t_list_jackknife(label='l96: amu(t)')
    plt.legend(scatterpoints=1)
    plt.show()

    # l96 amu_tcut
    amu.get_amu_tcut_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l96 amu_tcut_avg', amu_tcut_avg
    print 'l96 amu_tcut_err', amu_tcut_err
    print 'l96 amu_tcut_avg at tcut = 50:', amu_tcut_avg[50]
    print 'l96 amu_tcut_err at tcut = 50:', amu_tcut_err[50]
    amu.plt_amu_tcut_list_jackknife(label='l96: amu')
    plt.legend(scatterpoints=1)
    plt.show()

    # l96 wc_t with window
    amu = AmuT96192()
    amu.get_w_t_list()
    amu.get_window_t_list()
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu.plt_wc_t_list_jackknife(label='l96: amu(t) With Window Method')
    plt.legend(scatterpoints=1)
    plt.show()

    # l96 amu_tcut with window
    amu.get_amu_tcut_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l96 amu_tcut_avg window', amu_tcut_avg
    print 'l96 amu_tcut_err window', amu_tcut_err
    print 'l96 amu_tcut_avg window at tcut = 50:', amu_tcut_avg[50]
    print 'l96 amu_tcut_err window at tcut = 50:', amu_tcut_err[50]
    amu.plt_amu_tcut_list_jackknife(label='l96: amu With Window Method')

    plt.legend(scatterpoints=1)
    plt.show()

    # l96 extend ama lmasub avg
    CONFIG_FOR_EXTEND = [1464]
    print 'CONFIG_FOR_EXTEND for ama and lmasub', str(len(CONFIG_FOR_EXTEND)), CONFIG_FOR_EXTEND
    CONFIG_FOR_AVG = [504, 648, 696, 744, 792, 840, 888, 936, 984, 1032, 1128, 1176, 1224, 1272, 1368, 1512, 1560, 1608, 1656, 1704, 1752, 1800]
    print 'CONFIG_FOR_AVG for ama and lmasub', str(len(CONFIG_FOR_AVG)), CONFIG_FOR_AVG
    # readcmunu.extend_avg_to_Exact(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    # readcmunu.extend_avg_to_Sub(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    readcmunu.extend_avg_to_AMA(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    readcmunu.extend_avg_to_LMASUB(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)

    readcmunu.cmunu_subtract(CONFIG_FOR_EXTEND + CONFIG_FOR_AVG, 'Exact Sub AMA LMASUB LMA')

    # plot cmunu after extend
    for config in CONFIG_FOR_EXTEND + CONFIG_FOR_AVG:
        readcmunu.plot_cmumu_t(readcmunu.cmunu_config_subtract_dic[config])
    plt.show()

    # l96 wc_t with window
    amu = AmuT96192()
    amu.get_w_t_list()
    amu.get_window_t_list(delta=0.3)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu.plt_wc_t_list_jackknife(label='l96: amu(t) With Window Method')
    plt.legend(scatterpoints=1)
    plt.show()

    # l96 amu_tcut with window
    amu.get_amu_tcut_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l96 amu_tcut_avg window', amu_tcut_avg
    print 'l96 amu_tcut_err window', amu_tcut_err
    print 'l96 amu_tcut_avg window at tcut = 50:', amu_tcut_avg[50]
    print 'l96 amu_tcut_err window at tcut = 50:', amu_tcut_err[50]
    amu.plt_amu_tcut_list_jackknife(label='l96: amu With Window Method')

    plt.legend(scatterpoints=1)
    plt.show()

    # =================================
    # l6496 amu_t
    # =================================

    # Read Cmunu
    CMUNU_PATH = '/Users/tucheng/Desktop/Physics/research/hvp/analysis/data/HISQ_cmunu/l6496/'

    CMUNU_EXACT_PATH = CMUNU_PATH + 'l64_Exact_cmunu'
    CMUNU_SUB_PATH = CMUNU_PATH + 'l64_Sub_cmunu'
    CMUNU_AMA_PATH = CMUNU_PATH + 'l64_AMA_cmunu'
    CMUNU_LMA_PATH = CMUNU_PATH + 'l64_LMA_cmunu'
    CMUNU_LMASUB_PATH = CMUNU_PATH + 'l64_LMASUB_cmunu'

    readcmunu = cmunu_src.Do6496()

    readcmunu.E_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_EXACT_PATH)
    readcmunu.ES_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_SUB_PATH)
    readcmunu.A_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_AMA_PATH)
    readcmunu.L_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_LMA_PATH)
    readcmunu.LS_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_LMASUB_PATH)

    # Subtraction
    CONFIG_LIST_E = readcmunu.E_cmunu_config_dic.keys()[:]
    CONFIG_LIST_E.sort()
    CONFIG_LIST_E = CONFIG_LIST_E[:]
    print 'CONFIG_LIST E:', len(CONFIG_LIST_E), CONFIG_LIST_E
    for config in CONFIG_LIST_E:
        readcmunu.plot_cmumu_t(readcmunu.E_cmunu_config_dic[config])
    plt.show()

    CONFIG_LIST_ES = readcmunu.ES_cmunu_config_dic.keys()[:]
    CONFIG_LIST_ES.sort()
    CONFIG_LIST_ES = CONFIG_LIST_ES[:]
    print 'CONFIG_LIST ES:', len(CONFIG_LIST_ES), CONFIG_LIST_ES
    for config in CONFIG_LIST_ES:
        readcmunu.plot_cmumu_t(readcmunu.ES_cmunu_config_dic[config])
    plt.show()

    CONFIG_LIST_A = readcmunu.A_cmunu_config_dic.keys()[:]
    CONFIG_LIST_A.sort()
    CONFIG_LIST_A = CONFIG_LIST_A[:]
    print 'CONFIG_LIST A:', len(CONFIG_LIST_A), CONFIG_LIST_A
    for config in CONFIG_LIST_A:
        readcmunu.plot_cmumu_t(readcmunu.A_cmunu_config_dic[config])
    plt.show()

    CONFIG_LIST_L = readcmunu.L_cmunu_config_dic.keys()[:]
    CONFIG_LIST_L.sort()
    CONFIG_LIST_L = CONFIG_LIST_L[:]
    print 'CONFIG_LIST L:', len(CONFIG_LIST_L), CONFIG_LIST_L
    for config in CONFIG_LIST_L:
        readcmunu.plot_cmumu_t(readcmunu.L_cmunu_config_dic[config])
    plt.show()

    CONFIG_LIST_LS = readcmunu.LS_cmunu_config_dic.keys()[:]
    CONFIG_LIST_LS.sort()
    CONFIG_LIST_LS = CONFIG_LIST_LS[:]
    print 'CONFIG_LIST LS:', len(CONFIG_LIST_LS), CONFIG_LIST_LS
    for config in CONFIG_LIST_LS:
        readcmunu.plot_cmumu_t(readcmunu.LS_cmunu_config_dic[config])
    plt.show()

    CONFIG_LIST = [720, 732, 744, 768, 792, 804, 816, 828, 840, 852, 864, 876, 888, 900, 912, 924, 936, 948, 960, 984, 996, 1008, 1020, 1032, 1044, 1056, 1080, 1092, 1104, 1116, 1128, 1140, 1152, 1164, 1176]
    readcmunu.cmunu_subtract(CONFIG_LIST, 'Exact Sub AMA LMASUB LMA')

    # plot some cmunu
    for config in CONFIG_LIST:
        readcmunu.plot_cmumu_t(readcmunu.cmunu_config_subtract_dic[config])
    plt.show()

    # l64 LMA wc_t
    amu = AmuT6496()
    amu.get_w_t_list()
    amu.get_wc_t_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu.plt_wc_t_list_jackknife('l64: Subtraction')
    plt.legend(scatterpoints=1)
    plt.show()

    # l64 LMA amu_tcut
    amu.get_amu_tcut_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l64 amu_tcut_avg', amu_tcut_avg
    print 'l64 amu_tcut_err', amu_tcut_err
    amu.plt_amu_tcut_list_jackknife(label='l64: amu ')

    plt.legend(scatterpoints=1)
    plt.show()

    # l64 wc_t with window
    amu = AmuT6496()
    amu.get_w_t_list()
    amu.get_window_t_list()
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu.plt_wc_t_list_jackknife('l64: Subtraction With Window')
    plt.legend(scatterpoints=1)
    plt.show()

    # l64 amu_tcut with window
    amu.get_amu_tcut_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l64 amu_tcut_avg window', amu_tcut_avg
    print 'l64 amu_tcut_err window', amu_tcut_err
    amu.plt_amu_tcut_list_jackknife(label='l64: amu window ')

    plt.legend(scatterpoints=1)
    plt.show()

    print '============================================================================================================'

    # l64 extend ama and lmasub avg
    CONFIG_FOR_AVG = [720, 732, 744, 768, 804, 816, 828, 840, 852, 864, 876, 888, 900, 912, 924, 936, 948, 960, 984, 996, 1008, 1020, 1032, 1044, 1056, 1080, 1092, 1104, 1116, 1128, 1140, 1152, 1164, 1176]
    CONFIG_FOR_EXTEND = [708, 756, 780, 792, 972, 1068]
    print 'CONFIG_LIST To AVG:', len(CONFIG_FOR_AVG), CONFIG_FOR_AVG
    print 'CONFIG_LIST To Extend:', len(CONFIG_FOR_EXTEND), CONFIG_FOR_EXTEND
    readcmunu.extend_avg_to_Exact(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    readcmunu.extend_avg_to_Sub(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    readcmunu.extend_avg_to_AMA(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    readcmunu.extend_avg_to_LMASUB(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)

    readcmunu.cmunu_subtract(CONFIG_FOR_EXTEND + CONFIG_FOR_AVG, 'Exact Sub AMA LMASUB LMA')
    # plot some cmunu
    for config in CONFIG_FOR_EXTEND + CONFIG_FOR_AVG:
        readcmunu.plot_cmumu_t(readcmunu.cmunu_config_subtract_dic[config])
    plt.show()

    # l64 wc_t with window
    amu = AmuT6496()
    amu.get_w_t_list()
    amu.get_window_t_list(delta=0.3)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu.plt_wc_t_list_jackknife('l64: Subtraction With Window')
    plt.legend(scatterpoints=1)
    plt.show()

    # l64 amu_tcut with window
    amu.get_amu_tcut_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l64 amu_tcut_avg window', amu_tcut_avg
    print 'l64 amu_tcut_err window', amu_tcut_err
    amu.plt_amu_tcut_list_jackknife(label='l64: amu window ')

    plt.legend(scatterpoints=1)
    plt.show()

    # =================================
    # l4864 amu_t
    # =================================

    # Read Cmunu
    CMUNU_PATH = '/Users/tucheng/Desktop/Physics/research/hvp/analysis/data/HISQ_cmunu/l4864/'

    CMUNU_EXACT_PATH = CMUNU_PATH + 'l48_Exact_cmunu'
    CMUNU_SUB_PATH = CMUNU_PATH + 'l48_Sub_cmunu'
    CMUNU_AMA_PATH = CMUNU_PATH + 'l48_AMA_cmunu'
    CMUNU_LMA_PATH = CMUNU_PATH + 'l48_LMA_cmunu'
    CMUNU_LMASUB_PATH = CMUNU_PATH + 'l48_LMASUB_cmunu'

    readcmunu = cmunu_src.Do4864()

    readcmunu.E_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_EXACT_PATH)
    readcmunu.ES_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_SUB_PATH)
    readcmunu.A_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_AMA_PATH)
    readcmunu.L_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_LMA_PATH)
    readcmunu.LS_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_LMASUB_PATH)

    # Subtraction
    CONFIG_LIST_E = readcmunu.E_cmunu_config_dic.keys()[:]
    CONFIG_LIST_E.sort()
    CONFIG_LIST_E = CONFIG_LIST_E[:]
    print 'CONFIG_LIST E:', CONFIG_LIST_E
    CONFIG_LIST_ES = readcmunu.ES_cmunu_config_dic.keys()[:]
    CONFIG_LIST_ES.sort()
    CONFIG_LIST_ES = CONFIG_LIST_ES[:]
    print 'CONFIG_LIST ES:', CONFIG_LIST_ES
    CONFIG_LIST_A = readcmunu.A_cmunu_config_dic.keys()[:]
    CONFIG_LIST_A.sort()
    CONFIG_LIST_A = CONFIG_LIST_A[:]
    print 'CONFIG_LIST A:', CONFIG_LIST_A
    CONFIG_LIST_L = readcmunu.L_cmunu_config_dic.keys()[:]
    CONFIG_LIST_L.sort()
    CONFIG_LIST_L = CONFIG_LIST_L[:]
    print 'CONFIG_LIST L:', CONFIG_LIST_L
    CONFIG_LIST_LS = readcmunu.LS_cmunu_config_dic.keys()[:]
    CONFIG_LIST_LS.sort()
    CONFIG_LIST_LS = CONFIG_LIST_LS[:]
    print 'CONFIG_LIST LS:', CONFIG_LIST_LS

    readcmunu.cmunu_subtract(CONFIG_LIST_E, 'Exact Sub AMA LMASUB LMA')

    # plot some cmunu
    for config in CONFIG_LIST_E:
        readcmunu.plot_cmumu_t(readcmunu.cmunu_config_subtract_dic[config])
    plt.show()

    # l4864 wc_t
    amu = AmuT4864()
    amu.get_w_t_list()
    amu.get_wc_t_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu.plt_wc_t_list_jackknife('l48: Subtraction')
    plt.legend(scatterpoints=1)
    plt.show()

    # l4864 amu_tcut
    amu.get_amu_tcut_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l48 amu_tcut_avg', amu_tcut_avg
    print 'l48 amu_tcut_err', amu_tcut_err
    amu.plt_amu_tcut_list_jackknife(label='l48: amu ')

    plt.legend(scatterpoints=1)
    plt.show()

    # l4864 wc_t with window_t
    amu = AmuT4864()
    amu.get_w_t_list()
    amu.get_window_t_list(delta=0.3)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu.plt_wc_t_list_jackknife('l48: Subtraction With Window')
    plt.legend(scatterpoints=1)
    plt.show()

    # l4864 amu_tcut
    amu.get_amu_tcut_list_jackknife(readcmunu.cmunu_config_subtract_dic)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l48 amu_tcut_avg with window', amu_tcut_avg
    print 'l48 amu_tcut_err with window', amu_tcut_err
    amu.plt_amu_tcut_list_jackknife(label='l48: amu with window ')

    plt.legend(scatterpoints=1)
    plt.show()

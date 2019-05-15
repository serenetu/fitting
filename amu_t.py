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
        self.a
        self.ainv
        self.mu_lat
        self.fac
        self.T_MAX
        self.L
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
        # return lambda q2: (math.cos(q2 ** (1. / 2.) * tf) - 1.) / q2 + 1. / 2. * tf ** 2.

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

    def get_wc_t_list_from_jk_cmunu(self, cmunu_jk_dict):
        self.wt_ct_jackknife_dic = {}
        for config in cmunu_jk_dict:
            self.wt_ct_jackknife_dic[config] = self.get_wc_t_list(cmunu_jk_dict[config])
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
        return np.array(cmunu_fold)

    def amu_tcut_list(self, cmunu):
        wt_ct = self.get_wc_t_list(cmunu)
        amu_tcut = []
        for tcut in range(self.T_MAX / 2 + 1):
            amu_tcut.append(np.sum(wt_ct[:tcut+1]))
        return np.array(amu_tcut)

    def get_amu_tcut_list_from_jk_cmunu(self, cmunu_jk_dict):
        self.amu_tcut_list_jackknife_dict = {}
        for config in cmunu_jk_dict:
            self.amu_tcut_list_jackknife_dict[config] = self.amu_tcut_list(cmunu_jk_dict[config])

        return self.amu_tcut_list_jackknife_dict

    def get_amu_jk(self):
        try:
            self.amu_tcut_list_jackknife_dict
        except NameError:
            print 'Run Amu_T::get_amu_tcut_list_from_jk_cmunu First'
            return
        amu_jk = {}
        for config in self.amu_tcut_list_jackknife_dict:
            amu_jk[config] = [self.amu_tcut_list_jackknife_dict[config][-1]]
        mean, err = jackknife_avg_err_dic(amu_jk)
        return mean[0], err[0]

    def get_amu_jk_bound(self, t0, t1):
        amu_jk = {}
        for config in self.amu_upper_bound_dict:
            amu_jk[config] = [np.mean(np.append(self.amu_upper_bound_dict[config][t0:t1 + 1],
                                                self.amu_lower_bound_dict[config][t0:t1 + 1]))]
        mean, err = jackknife_avg_err_dic(amu_jk)
        return mean[0], err[0]

    def get_amu_tcut_list_jackknife_avg_err(self):
        try:
            self.amu_tcut_list_jackknife_dict
        except NameError:
            print 'Run Amu_T::amu_tcut_list_jackknife() First'
        else:
            mean, err = jackknife_avg_err_dic(self.amu_tcut_list_jackknife_dict)
            return mean, err

    def get_wc_t_list_upper_bound(self, cmunu, t):
        try:
            self.w_tlist
        except NameError:
            print 'Run Amu_T::get_w_tlist() First'
            return
        res = []
        fcmunu = self.fold_cmunu(cmunu)
        e0 = 2. * math.sqrt(self.pi_lat ** 2. + (2. * math.pi / self.L) ** 2.)
        for tt in range(len(fcmunu)):
            if tt > t:
                fcmunu[tt] = fcmunu[t] * math.exp(-e0 * (tt - t))
            res.append(self.w_tlist[tt] * fcmunu[tt])
        return np.array(res)

    def get_wc_t_list_lower_bound(self, cmunu, t):
        try:
            self.w_tlist
        except NameError:
            print 'Run Amu_T::get_w_tlist() First'
            return
        res = []
        fcmunu = self.fold_cmunu(cmunu)
        for tt in range(len(fcmunu)):
            if tt > t:
                fcmunu[tt] = 0.
            res.append(self.w_tlist[tt] * fcmunu[tt])
        return np.array(res)

    def get_amu_list_upper_bound(self, cmunu):
        amu_tcut = []
        for tcut in range(self.T_MAX / 2 + 1):
            wt_ct = self.get_wc_t_list_upper_bound(cmunu, tcut)
            amu_tcut.append(np.sum(wt_ct))
        return np.array(amu_tcut)

    def get_amu_list_lower_bound(self, cmunu):
        amu_tcut = []
        for tcut in range(self.T_MAX / 2 + 1):
            wt_ct = self.get_wc_t_list_lower_bound(cmunu, tcut)
            amu_tcut.append(np.sum(wt_ct))
        return np.array(amu_tcut)

    def get_amu_list_upper_bound_jk(self, cmunu_jk_dict):
        self.amu_upper_bound_dict = {}
        for config in cmunu_jk_dict:
            self.amu_upper_bound_dict[config] = self.get_amu_list_upper_bound(cmunu_jk_dict[config])
        return self.amu_upper_bound_dict

    def get_amu_list_lower_bound_jk(self, cmunu_jk_dict):
        self.amu_lower_bound_dict = {}
        for config in cmunu_jk_dict:
            self.amu_lower_bound_dict[config] = self.get_amu_list_lower_bound(cmunu_jk_dict[config])
        return self.amu_lower_bound_dict

    def plt_amu_upper_bound_jk(self, label=''):
        amu, amu_err = jackknife_avg_err_dic(self.amu_upper_bound_dict)
        t_list = range(len(amu))
        plt.errorbar(t_list, amu, amu_err, fmt='x', capsize=2, elinewidth=1, label=label, color='black')
        return

    def plt_amu_lower_bound_jk(self, label=''):
        amu, amu_err = jackknife_avg_err_dic(self.amu_lower_bound_dict)
        t_list = range(len(amu))
        plt.errorbar(t_list, amu, amu_err, fmt='+', capsize=2, elinewidth=1, label=label, color='red')
        return

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
        # amu, amu_err = self.get_amu_tcut_list_jackknife_avg_err()
        amu, amu_err = jackknife_avg_err_dic(self.amu_tcut_list_jackknife_dict)
        t_list = range(len(amu))
        plt.errorbar(t_list, amu, amu_err, fmt='.', elinewidth=1, label = label)
        return


class AmuT4864(AmuT):

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        #self.ainv = 1.629278350515464  # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        # self.a = 0.12224 / 0.1975      # GeV-1  # From Email
        self.a = 0.12121 / 0.1975      # GeV-1  # From Email
        self.ainv = 1. / self.a
        self.mu_lat = 0.1056583715 * self.a
        self.fac = 4. * self.alpha ** 2. * (5. / 9.) * 10. ** 10. * 2.
        self.T_MAX = 64
        self.L = 48
        # self.pi_lat = 0.1349766 * self.a
        self.pi_lat = 3.91 / self.L
        return


class AmuT6496(AmuT):

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        # self.a = 1. / 2.255453                # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        # self.a = 0.08786 / 0.1975             # GeV-1   # From Email
        self.a = 0.08787 / 0.1975             # GeV-1   # From Email
        self.ainv = 1. / self.a
        self.mu_lat = 0.1056583715 * self.a
        self.fac = 4. * self.alpha ** 2. * (5./9.) * 10. ** 10. * 2.
        self.T_MAX = 96
        self.L = 64
        # self.pi_lat = 0.1349766 * self.a
        self.pi_lat = 3.66 / self.L
        return


class AmuT96192(AmuT):

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        # self.a = 1. / 2.255453                # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        # self.a = 0.05662 / 0.1975             # GeV-1 (1 GeV-1 = 0.1975 fm) # From Email
        self.a = 0.05684 / 0.1975
        self.ainv = 1. / self.a
        self.mu_lat = 0.1056583715 * self.a
        self.fac = 4. * self.alpha ** 2. * (5./9.) * 10. ** 10. * 2.
        self.T_MAX = 192
        self.L = 96
        # self.pi_lat = 0.1349766 * self.a
        self.pi_lat = 3.73 / self.L
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

    # Config List
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

    # l96 extend ama lmasub avg
    CONFIG_FOR_EXTEND = [1464]
    print 'CONFIG_FOR_EXTEND for ama and lmasub', str(len(CONFIG_FOR_EXTEND)), CONFIG_FOR_EXTEND
    CONFIG_FOR_AVG = [504, 648, 696, 744, 792, 840, 888, 936, 984, 1032, 1128, 1176, 1224, 1272, 1368, 1512, 1560, 1608, 1656, 1704, 1752, 1800]
    print 'CONFIG_FOR_AVG for ama and lmasub', str(len(CONFIG_FOR_AVG)), CONFIG_FOR_AVG
    readcmunu.compute_jk_exact(config_list='All')
    readcmunu.compute_jk_sub(config_list='All')
    readcmunu.compute_jk_ama(CONFIG_FOR_AVG)
    readcmunu.compute_jk_lmasub(CONFIG_FOR_AVG)
    readcmunu.compute_jk_lma(config_list='All')
    readcmunu.extend_avg_to_ama_jk(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    readcmunu.extend_avg_to_lmasub_jk(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)

    # E - S + A - LS + L
    readcmunu.cmunu_jk_subtract(CONFIG_FOR_EXTEND + CONFIG_FOR_AVG, 'Exact Sub AMA LMASUB LMA')

    # l96 wc_t with window
    amu = AmuT96192()
    amu.get_w_t_list()
    amu.get_window_t_list(t0=0.4, t1=1.0, delta=0.15)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l96: amu(t) With Window')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l96 amu_tcut with window
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l96 amu_tcut_avg window', amu_tcut_avg
    print 'l96 amu_tcut_err window', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l96: amu window ')
    # plt.legend(scatterpoints=1)
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l96 amu window 0.4 1.0 0.15:', amu_mean, amu_err)

    # l96 wc_t with window
    amu = AmuT96192()
    amu.get_w_t_list()
    amu.get_window_t_list(t0=0.4, t1=1.0, delta=0.3)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l96: amu(t) With Window')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l96 amu_tcut with window
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l96 amu_tcut_avg window', amu_tcut_avg
    print 'l96 amu_tcut_err window', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l96: amu window ')
    # plt.legend(scatterpoints=1)
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l96 amu window 0.4 1.0 0.3:', amu_mean, amu_err)

    # l96 wc_t with window
    amu = AmuT96192()
    amu.get_w_t_list()
    amu.get_window_t_list(t0=0.4, t1=1.3, delta=0.15)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l96: amu(t) With Window')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l96 amu_tcut with window
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l96 amu_tcut_avg window', amu_tcut_avg
    print 'l96 amu_tcut_err window', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l96: amu window ')
    # plt.legend(scatterpoints=1)
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l96 amu window 0.4 1.3 0.15:', amu_mean, amu_err)

    # l96 wc_t
    amu = AmuT96192()
    amu.get_w_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l96: jk Subtraction')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l96 amu_tcut
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l96 amu_tcut_avg', amu_tcut_avg
    print 'l96 amu_tcut_err', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l96: amu ')
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l96 amu:', amu_mean, amu_err)

    amu.get_amu_list_upper_bound_jk(readcmunu.cmunu_jk_subtract_dict)
    amu.get_amu_list_lower_bound_jk(readcmunu.cmunu_jk_subtract_dict)
    amu_mean, amu_err = amu.get_amu_jk_bound(45, 49)
    print('l96 amu bound:', amu_mean, amu_err)
    amu.plt_amu_lower_bound_jk(label='lower bound')
    amu.plt_amu_upper_bound_jk(label='upper bound')
    plt.ylim(bottom=0, top=1000)
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

    # Config List
    CONFIG_LIST_E = readcmunu.E_cmunu_config_dic.keys()[:]
    CONFIG_LIST_E.sort()
    print 'CONFIG_LIST E:', len(CONFIG_LIST_E), CONFIG_LIST_E
    # for config in CONFIG_LIST_E:
        # readcmunu.plot_cmumu_t(readcmunu.E_cmunu_config_dic[config])
    # plt.show()
    CONFIG_LIST_ES = readcmunu.ES_cmunu_config_dic.keys()[:]
    CONFIG_LIST_ES.sort()
    print 'CONFIG_LIST ES:', len(CONFIG_LIST_ES), CONFIG_LIST_ES
    # for config in CONFIG_LIST_ES:
        # readcmunu.plot_cmumu_t(readcmunu.ES_cmunu_config_dic[config])
    # plt.show()
    CONFIG_LIST_A = readcmunu.A_cmunu_config_dic.keys()[:]
    CONFIG_LIST_A.sort()
    print 'CONFIG_LIST A:', len(CONFIG_LIST_A), CONFIG_LIST_A
    # for config in CONFIG_LIST_A:
        # readcmunu.plot_cmumu_t(readcmunu.A_cmunu_config_dic[config])
    # plt.show()
    CONFIG_LIST_L = readcmunu.L_cmunu_config_dic.keys()[:]
    CONFIG_LIST_L.sort()
    print 'CONFIG_LIST L:', len(CONFIG_LIST_L), CONFIG_LIST_L
    # for config in CONFIG_LIST_L:
        # readcmunu.plot_cmumu_t(readcmunu.L_cmunu_config_dic[config])
    # plt.show()
    CONFIG_LIST_LS = readcmunu.LS_cmunu_config_dic.keys()[:]
    CONFIG_LIST_LS.sort()
    print 'CONFIG_LIST LS:', len(CONFIG_LIST_LS), CONFIG_LIST_LS
    # for config in CONFIG_LIST_LS:
        # readcmunu.plot_cmumu_t(readcmunu.LS_cmunu_config_dic[config])
    # plt.show()

    # l64 extend ama and lmasub avg
    CONFIG_FOR_AVG = [708, 720, 732, 744, 768, 792, 804, 816, 828, 840, 852, 864, 876, 888, 900, 912, 924, 936, 948, 960, 984, 996, 1008, 1020, 1032, 1044, 1056, 1080, 1092, 1104, 1116, 1128, 1140, 1152, 1164, 1176]
    CONFIG_FOR_EXTEND = [756, 780, 972, 1068]
    print 'CONFIG_LIST To AVG:', len(CONFIG_FOR_AVG), CONFIG_FOR_AVG
    print 'CONFIG_LIST To Extend:', len(CONFIG_FOR_EXTEND), CONFIG_FOR_EXTEND
    readcmunu.compute_jk_exact(CONFIG_FOR_AVG)
    readcmunu.compute_jk_sub(CONFIG_FOR_AVG)
    readcmunu.compute_jk_ama(CONFIG_FOR_AVG)
    readcmunu.compute_jk_lmasub(CONFIG_FOR_AVG)
    readcmunu.compute_jk_lma(config_list='All')
    readcmunu.extend_avg_to_exact_jk(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    readcmunu.extend_avg_to_sub_jk(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    readcmunu.extend_avg_to_ama_jk(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)
    readcmunu.extend_avg_to_lmasub_jk(CONFIG_FOR_AVG, CONFIG_FOR_EXTEND)

    # E - S + A - LS + L
    readcmunu.cmunu_jk_subtract(CONFIG_FOR_EXTEND + CONFIG_FOR_AVG, 'Exact Sub AMA LMASUB LMA')

    # plot some cmunu
    # for config in CONFIG_FOR_EXTEND + CONFIG_FOR_AVG:
        # readcmunu.plot_cmumu_t(readcmunu.cmunu_jk_subtract_dict[config])
    # plt.show()

    # l64 wc_t with window
    amu = AmuT6496()
    amu.get_w_t_list()
    amu.get_window_t_list(t0=0.4, t1=1.0, delta=0.15)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l64: Subtraction With Window')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l64 amu_tcut with window
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l64 amu_tcut_avg window', amu_tcut_avg
    print 'l64 amu_tcut_err window', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l64: amu window ')
    # plt.legend(scatterpoints=1)
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l64 amu window 0.4 1.0 0.15:', amu_mean, amu_err)

    # l64 wc_t with window
    amu = AmuT6496()
    amu.get_w_t_list()
    amu.get_window_t_list(t0=0.4, t1=1.0, delta=0.3)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l64: Subtraction With Window')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l64 amu_tcut with window
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l64 amu_tcut_avg window', amu_tcut_avg
    print 'l64 amu_tcut_err window', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l64: amu window ')
    # plt.legend(scatterpoints=1)
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l64 amu window 0.4 1.0 0.3:', amu_mean, amu_err)

    # l64 wc_t with window
    amu = AmuT6496()
    amu.get_w_t_list()
    amu.get_window_t_list(t0=0.4, t1=1.3, delta=0.15)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l64: Subtraction With Window')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l64 amu_tcut with window
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l64 amu_tcut_avg window', amu_tcut_avg
    print 'l64 amu_tcut_err window', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l64: amu window ')
    # plt.legend(scatterpoints=1)
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l64 amu window 0.4 1.3 0.15:', amu_mean, amu_err)

    # l64 wc_t
    amu = AmuT6496()
    amu.get_w_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l64: jk Subtraction')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l64 amu_tcut
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l64 amu_tcut_avg', amu_tcut_avg
    print 'l64 amu_tcut_err', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l64: amu ')
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l64 amu:', amu_mean, amu_err)

    amu.get_amu_list_upper_bound_jk(readcmunu.cmunu_jk_subtract_dict)
    amu.get_amu_list_lower_bound_jk(readcmunu.cmunu_jk_subtract_dict)
    amu_mean, amu_err = amu.get_amu_jk_bound(31, 36)
    print('l64 amu bound:', amu_mean, amu_err)
    amu.plt_amu_lower_bound_jk(label='lower bound')
    amu.plt_amu_upper_bound_jk(label='upper bound')
    plt.ylim(bottom=0, top=900)
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

    # Config List
    CONFIG_LIST_E = readcmunu.E_cmunu_config_dic.keys()[:]
    CONFIG_LIST_E.sort()
    print 'CONFIG_LIST E:', len(CONFIG_LIST_E), CONFIG_LIST_E
    CONFIG_LIST_ES = readcmunu.ES_cmunu_config_dic.keys()[:]
    CONFIG_LIST_ES.sort()
    print 'CONFIG_LIST ES:', len(CONFIG_LIST_ES), CONFIG_LIST_ES
    CONFIG_LIST_A = readcmunu.A_cmunu_config_dic.keys()[:]
    CONFIG_LIST_A.sort()
    print 'CONFIG_LIST A:', len(CONFIG_LIST_A), CONFIG_LIST_A
    CONFIG_LIST_L = readcmunu.L_cmunu_config_dic.keys()[:]
    CONFIG_LIST_L.sort()
    print 'CONFIG_LIST L:', len(CONFIG_LIST_L), CONFIG_LIST_L
    CONFIG_LIST_LS = readcmunu.LS_cmunu_config_dic.keys()[:]
    CONFIG_LIST_LS.sort()
    print 'CONFIG_LIST LS:', len(CONFIG_LIST_LS), CONFIG_LIST_LS

    # compute jk
    readcmunu.compute_jk_exact(config_list='All')
    readcmunu.compute_jk_sub(config_list='All')
    readcmunu.compute_jk_ama(config_list='All')
    readcmunu.compute_jk_lmasub(config_list='All')
    readcmunu.compute_jk_lma(config_list='All')

    # E - S + A - LS + L
    readcmunu.cmunu_jk_subtract(CONFIG_LIST_E, 'Exact Sub AMA LMASUB LMA')

    # plot jk cmunu
    # for config in CONFIG_LIST_E:
        # readcmunu.plot_cmumu_t(readcmunu.cmunu_jk_subtract_dict[config])
    # plt.show()

    # l4864 wc_t with window
    amu = AmuT4864()
    amu.get_w_t_list()
    amu.get_window_t_list(t0=0.4, t1=1.0, delta=0.15)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l48: jk Subtraction With Window')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l4864 amu_tcut with window
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l48 amu_tcut_avg with window', amu_tcut_avg
    print 'l48 amu_tcut_err with window', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l48: amu with window ')
    # plt.legend(scatterpoints=1)
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l48 amu window (0.4, 1.0, 0.15):', amu_mean, amu_err)

    # l4864 wc_t with window
    amu = AmuT4864()
    amu.get_w_t_list()
    amu.get_window_t_list(t0=0.4, t1=1.0, delta=0.3)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l48: jk Subtraction With Window')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l4864 amu_tcut with window
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l48 amu_tcut_avg with window', amu_tcut_avg
    print 'l48 amu_tcut_err with window', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l48: amu with window ')
    # plt.legend(scatterpoints=1)
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l48 amu window (0.4, 1.0, 0.3):', amu_mean, amu_err)

    # l4864 wc_t with window
    amu = AmuT4864()
    amu.get_w_t_list()
    amu.get_window_t_list(t0=0.4, t1=1.3, delta=0.15)
    amu.update_w_t_list_with_window_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l48: jk Subtraction With Window')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l4864 amu_tcut with window
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l48 amu_tcut_avg with window', amu_tcut_avg
    print 'l48 amu_tcut_err with window', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l48: amu with window ')
    # plt.legend(scatterpoints=1)
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l48 amu window (0.4, 1.3, 0.15):', amu_mean, amu_err)

    # l4864 wc_t
    amu = AmuT4864()
    amu.get_w_t_list()
    amu.get_wc_t_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    # amu.plt_wc_t_list_jackknife('l48: jk Subtraction')
    # plt.legend(scatterpoints=1)
    # plt.show()

    # l4864 amu_tcut
    amu.get_amu_tcut_list_from_jk_cmunu(readcmunu.cmunu_jk_subtract_dict)
    amu_tcut_avg, amu_tcut_err = amu.get_amu_tcut_list_jackknife_avg_err()
    print 'l48 amu_tcut_avg', amu_tcut_avg
    print 'l48 amu_tcut_err', amu_tcut_err
    # amu.plt_amu_tcut_list_jackknife(label='l48: amu ')
    # plt.show()
    amu_mean, amu_err = amu.get_amu_jk()
    print('l48 amu:', amu_mean, amu_err)

    amu.get_amu_list_upper_bound_jk(readcmunu.cmunu_jk_subtract_dict)
    amu.get_amu_list_lower_bound_jk(readcmunu.cmunu_jk_subtract_dict)
    amu_mean, amu_err = amu.get_amu_jk_bound(23, 26)
    print('l48 amu bound:', amu_mean, amu_err)
    amu.plt_amu_lower_bound_jk(label='lower bound')
    amu.plt_amu_upper_bound_jk(label='upper bound')
    plt.ylim(bottom=0, top=900)
    plt.legend(scatterpoints=1)
    plt.show()

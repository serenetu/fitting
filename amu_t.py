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
import fit_from_cmatrix as amu_src


def mean_dic(dic):
    array = []
    for key in dic:
        array.append(np.array(dic[key]))
    array = np.array(array)
    mean = np.mean(array, axis = 0)
    return mean


if __name__ == "__main__":

    '''
    #-----------> l4864 <--------------
    # >>>>> setting <<<<<<
    HISQ_PATH = '/Users/tucheng/Desktop/transfer/HISQ_extract/l4864/'
    PIQSQU_OUTPUT_FILE = '/Users/tucheng/Desktop/fitting/results/l48/l48-cmatrix'

    CMUNU_EXACT_PATH  = '/Users/tucheng/Desktop/fitting/results/l48/l48_Exact_cmunu'
    CMUNU_SUB_PATH    = '/Users/tucheng/Desktop/fitting/results/l48/l48_Sub_cmunu'
    CMUNU_AMA_PATH    = '/Users/tucheng/Desktop/fitting/results/l48/l48_AMA_cmunu'
    CMUNU_LMA_PATH    = '/Users/tucheng/Desktop/fitting/results/l48/l48_LMA_cmunu'
    CMUNU_LMASUB_PATH = '/Users/tucheng/Desktop/fitting/results/l48/l48_LMASUB_cmunu'

    LMA_LABEL = 'l4864f211b600m00184m0507m628a-LMA'
    LMASUB_LABEL = 'l4864f211b600m00184m0507m628a-LMASUB'
    EXACT_LABEL = 'l4864f211b600m00184m0507m628a-EXACT'
    SUB_LABEL = 'l4864f211b600m00184m0507m628a-SUB'
    AMA_LABEL = 'l4864f211b600m00184m0507m628a-AMA'

    readcmunu = cmunu_src.Do4864(HISQ_PATH)

    # >>>>> read cmunu <<<<<<

    readcmunu.E_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_EXACT_PATH)
    readcmunu.ES_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_SUB_PATH)
    readcmunu.A_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_AMA_PATH)
    readcmunu.L_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_LMA_PATH)
    readcmunu.LS_cmunu_config_dic = readcmunu.read_cmunu(CMUNU_LMASUB_PATH)

    CONFIG_LIST = readcmunu.L_cmunu_config_dic.keys()
    CONFIG_LIST.sort()
    print 'CONFIG_LIST:', CONFIG_LIST, 'num:', len(CONFIG_LIST)

    readcmunu.cmunu_subtract(CONFIG_LIST, 'Exact Sub AMA LMASUB LMA')
    readcmunu.cmunu_subtract(CONFIG_LIST, 'LMA')
    # plot some cmunu
    for config in CONFIG_LIST:
        readcmunu.plot_cmumu_t(readcmunu.cmunu_config_subtract_dic[config])
    plt.show()


    # LMA wc_t
    LMA_CMUNU_DIC = readcmunu.L_cmunu_config_dic

    amu_LMA_l48 = amu_src.AmuT_4864()

    amu_LMA_l48.get_w_t_list()
    print 'l48 wt', amu_LMA_l48.w_tlist
    amu_LMA_l48.get_wc_t_list_jackknife(LMA_CMUNU_DIC)
    amu_LMA_l48.plt_wc_t_list_jackknife('l48: LMA')
    #plt.show()


    # subtract Exact Sub AMA LMASUB LMA
    SUBTRACT_TYPE = 'Exact Sub AMA LMASUB LMA'
    readcmunu.cmunu_subtract(CONFIG_LIST, SUBTRACT_TYPE)

    EAL_CMUNU_DIC = readcmunu.cmunu_config_subtract_dic

    # Exact Sub AMA LMASUB LMA wc_t
    amu_EAL = amu_src.AmuT_4864()
    amu_EAL.get_w_t_list()
    amu_EAL.get_wc_t_list_jackknife(EAL_CMUNU_DIC)
    amu_EAL.plt_wc_t_list_jackknife('l48: EAL')

    plt.legend(scatterpoints=1)
    plt.show()



    # LMA amu_tcut
    amu_LMA_l48.get_amu_tcut_list_jackknife(LMA_CMUNU_DIC)
    LMA_amu_tcut_avg, LMA_amu_tcut_err = amu_LMA_l48.get_amu_tcut_list_jackknife_avg_err()
    print 'l48 LMA_amu_tcut_avg', LMA_amu_tcut_avg
    print 'l48 LMA_amu_tcut_err', LMA_amu_tcut_err
    amu_LMA_l48.plt_amu_tcut_list_jackknife(label='l48: amu LMA')

    # Exact Sub AMA LMASUB LMA amu_tcut
    amu_EAL.get_amu_tcut_list_jackknife(EAL_CMUNU_DIC)
    EAL_amu_tcut_avg, EAL_amu_tcut_err = amu_EAL.get_amu_tcut_list_jackknife_avg_err()
    print 'l48 EAL_amu_tcut_avg', EAL_amu_tcut_avg
    print 'l48 EAL_amu_tcut_err', EAL_amu_tcut_err
    amu_EAL.plt_amu_tcut_list_jackknife(label='l48: amu EAL')


    plt.legend(scatterpoints=1)
    plt.show()
    '''




    '''
    # subtract Exact Sub
    SUBTRACT_TYPE = 'Exact Sub'
    readcmunu.cmunu_subtract(CONFIG_LIST, SUBTRACT_TYPE)

    # subtract Exact Sub jackknife
    subtract_jack = cmunu_src.Cmunu_Jackknife(CONFIG_LIST, readcmunu.cmunu_config_subtract_dic)

    cmunu_jackknife_dic = subtract_jack.get_cmunu_jackknife()

    CMUNU_MEAN, CMUNU_ERR = subtract_jack.cmunu_avg_err()

    # subtract Exact Sub amu
    amu = amu_src.AmuT_4864()
    T_LIST = [i for i in range(64)]
    amu.get_w_t_list()
    amu_t = amu.get_wc_t_list(CMUNU_MEAN)
    amu_t_error = amu.get_wc_t_list(CMUNU_ERR)
    print 'amu(t) Exact Sub:', amu_t
    amu.plt_wc_t_list_jackknife('l48: Exact Sub')
    '''

    '''
    # =================================
    # l96192 amu_t
    # =================================

    # Read Cmunu
    HISQ_PATH = ''

    CMUNU_EXACT_PATH  = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l96192/l96_Exact_cmunu'
    CMUNU_SUB_PATH    = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l96192/l96_Sub_cmunu'
    CMUNU_AMA_PATH    = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l96192/l96_AMA_cmunu'
    CMUNU_LMA_PATH    = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l96192/l96_LMA_cmunu'
    CMUNU_LMASUB_PATH = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l96192/l96_LMASUB_cmunu'

    readcmunu = cmunu_src.Do96192(HISQ_PATH)

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
    amu = amu_src.AmuT_96192()
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
    amu = amu_src.AmuT_96192()
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
    amu = amu_src.AmuT_96192()
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
    exit()
    '''

    # =================================
    # l6496 amu_t
    # =================================

    # Read Cmunu
    HISQ_PATH = ''

    CMUNU_EXACT_PATH  = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l6496/l64_Exact_cmunu'
    CMUNU_SUB_PATH    = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l6496/l64_Sub_cmunu'
    CMUNU_AMA_PATH    = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l6496/l64_AMA_cmunu'
    CMUNU_LMA_PATH    = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l6496/l64_LMA_cmunu'
    CMUNU_LMASUB_PATH = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l6496/l64_LMASUB_cmunu'

    readcmunu = cmunu_src.Do6496(HISQ_PATH)

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
    amu = amu_src.AmuT_6496()
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
    amu = amu_src.AmuT_6496()
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
    # CONFIG_LIST = [720, 732, 744, 768, 792, 804, 816, 828, 840, 852, 864, 876, 888, 900, 912, 924, 936, 948, 960, 984, 996, 1008, 1020, 1032, 1044, 1056, 1080, 1092, 1104, 1116, 1128, 1140, 1152, 1164, 1176]
    CONFIG_LIST = [720, 732, 744, 768, 804, 816, 828, 840, 852, 864, 876, 888, 900, 912, 924, 936, 948, 960, 984, 996, 1008, 1020, 1032, 1044, 1056, 1080, 1092, 1104, 1116, 1128, 1140, 1152, 1164, 1176]
    readcmunu.cmunu_subtract(CONFIG_LIST, 'Exact Sub AMA LMASUB LMA')

    # l64 extend ama and lmasub avg
    CONFIG_FOR_AVG = [720, 732, 744, 768, 804, 816, 828, 840, 852, 864, 876, 888, 900, 912, 924, 936, 948, 960, 984, 996, 1008, 1020, 1032, 1044, 1056, 1080, 1092, 1104, 1116, 1128, 1140, 1152, 1164, 1176]
    CONFIG_FOR_EXTEND = [708, 756, 780, 792, 972, 1068]
    # CONFIG_FOR_AVG = [792, 804, 984, 996, 1008, 1020, 1032, 1044, 1056, 1080, 1092, 1104, 1116, 1128, 1140, 1152, 1164, 1176]
    # CONFIG_FOR_EXTEND = [708, 720, 732, 744, 756, 768, 780, 816, 828, 840, 852, 864, 876, 888, 900, 912, 924, 936, 948, 960, 972, 1068]
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
    #exit()

    # l64 wc_t with window
    amu = amu_src.AmuT_6496()
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

    '''
    # =================================
    # l4864 amu_t
    # =================================

    # Read Cmunu
    HISQ_PATH = ''

    CMUNU_EXACT_PATH  = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l4864/l48_Exact_cmunu'
    CMUNU_SUB_PATH    = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l4864/l48_Sub_cmunu'
    CMUNU_AMA_PATH    = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l4864/l48_AMA_cmunu'
    CMUNU_LMA_PATH    = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l4864/l48_LMA_cmunu'
    CMUNU_LMASUB_PATH = '/Users/tucheng/Desktop/Physics/research/hvp/HISQ_cmunu/l4864/l48_LMASUB_cmunu'

    readcmunu = cmunu_src.Do4864(HISQ_PATH)

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
    amu = amu_src.AmuT_4864()
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
    amu = amu_src.AmuT_4864()
    amu.get_w_t_list()
    amu.get_window_t_list()
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
    '''

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
import functools


def log(function):
    @functools.wraps(function)
    def wrapper(self, *args, **kwargs):
        print 'Call %s::%s():' % (self.__class__.__name__, function.__name__)
        return function(self, *args, **kwargs)
    return wrapper

#?????????
class Read_AVG:
    def __init__(self, path):
        data = open(path)
        data = data.readlines()
        self.qsqu_all = []
        self.pi_qsqu_all = []
        self.points_num_all = - 1

        for line in data:
            if ('c_matrix:' in line) or ('COVARIANCE:' in line):
                break
            self.points_num_all += 1 #Jump over the first line which is always the label called 'AVERAGES:'
            if self.points_num_all > 0:
                self.qsqu_all.append(float(data[self.points_num_all].split()[0]))
                self.pi_qsqu_all.append(float(data[self.points_num_all].split()[1]))
        self.qsqu_all = np.array(self.qsqu_all)
        self.pi_qsqu_all = np.array(self.pi_qsqu_all)

    def split_qsqu_pi_qsqu_avg_c_matrix(self, points_num):
        self.qsqu = np.hsplit(self.qsqu_all, np.array([points_num]))[0]
        self.pi_qsqu = np.hsplit(self.pi_qsqu_all, np.array([points_num]))[0]

    def plots(self, points_num, typ = '.b'):
        self.split_qsqu_pi_qsqu_avg_c_matrix(points_num)
        plt.errorbar(self.qsqu, self.pi_qsqu, fmt=typ, elinewidth=1)

class read_pi_qsqu_cmatrix_file:

    def __init__(self, path):
        data = open(path)
        data = data.readlines()
        self.qsqu_all = []
        self.pi_qsqu_all = []
        self.cmatrix_all = []
        self.points_num_all = - 1

        for line in data:
            if ('c_matrix:' in line) or ('COVARIANCE:' in line):
                break
            self.points_num_all += 1
            if self.points_num_all > 0:
                self.qsqu_all.append(float(data[self.points_num_all].split()[0]))
                self.pi_qsqu_all.append(float(data[self.points_num_all].split()[1]))
        for i in range(0, self.points_num_all):
            self.cmatrix_all.append([])
            for j in range(0, self.points_num_all):
                self.cmatrix_all[i].append(float(data[self.points_num_all + 2 + self.points_num_all * i + j].split()[2]))
        self.qsqu_all = np.array(self.qsqu_all)
        self.pi_qsqu_all = np.array(self.pi_qsqu_all)
        self.cmatrix_all = np.array(self.cmatrix_all)
        return

    def split_qsqu_pi_qsqu_avg_c_matrix(self, points_num):
        self.qsqu = np.hsplit(self.qsqu_all, np.array([points_num]))[0]
        self.pi_qsqu = np.hsplit(self.pi_qsqu_all, np.array([points_num]))[0]
        self.cmatrix = np.hsplit(self.cmatrix_all, np.array([points_num]))[0]
        self.cmatrix = np.vsplit(self.cmatrix, np.array([points_num]))[0]
        self.diagCmatrix = np.diag(np.diag(self.cmatrix))
        return

    def plots(self, points_num, typ = '.b'):
        self.split_qsqu_pi_qsqu_avg_c_matrix(points_num)
        plt.errorbar(self.qsqu, self.pi_qsqu, yerr=np.diag(self.cmatrix) ** (1.0 / 2.0), fmt=typ, elinewidth=1)
        #plt.show()
        return

class ReadPiQsqu:

    def __init__(self, path):
        data = open(path)
        data = data.readlines()
        self.qsqu_all = []
        self.pi_qsqu_all = []
        self.points_num_all = - 1

        for line in data:
            if ('c_matrix:' in line) or ('COVARIANCE:' in line):
                break
            self.points_num_all += 1
            if self.points_num_all > 0:
                self.qsqu_all.append(float(data[self.points_num_all].split()[0]))
                self.pi_qsqu_all.append(float(data[self.points_num_all].split()[1]))
        self.qsqu_all = np.array(self.qsqu_all)
        self.pi_qsqu_all = np.array(self.pi_qsqu_all)
        return

    def plots(self, points_num, typ = '.b'):
        plt.errorbar(self.qsqu_all[:points_num], self.pi_qsqu_all[:points_num], fmt=typ, elinewidth=1)
        #plt.show()
        return

def func_inverse(matrix, dis = 0):
    matrix_inverse = np.linalg.inv(matrix)
    if dis != 0:
        print
        print 'inverse'
        print matrix_inverse
        print
        print 'unit test'
        print matrix.dot(matrix_inverse)
        print
    return matrix_inverse

def func_pi(parameters, x, type):
    m = len(parameters)
    function_pi_sum_part = 0
    if type == 'NOT DDS':
        if m % 2 == 0:
            for i in range(1, (m - 2) / 2 + 1):
                function_pi_sum_part += parameters[i] / (parameters[i + (m - 2) / 2] + x)
            function_pi_sum = parameters[0] - x * (parameters[m - 1] + function_pi_sum_part)
            return function_pi_sum
        elif m % 2 == 1:
            for i in range(1, (m - 1) / 2 + 1):
                function_pi_sum_part += parameters[i] / (parameters[i + (m - 1) / 2] + x)
            function_pi_sum = parameters[0] - x * function_pi_sum_part
            return function_pi_sum
    if type == 'DDS':
        if m % 2 == 0:
            for i in range(0, m / 2):
                function_pi_sum_part += parameters[i] / (parameters[i + m / 2] + x)
            function_pi_sum = x * function_pi_sum_part
            return function_pi_sum
        elif m % 2 == 1:
            for i in range(0, (m - 1) / 2):
                function_pi_sum_part += parameters[i] / (parameters[i + (m - 1) / 2] + x)
            function_pi_sum = x * (parameters[m - 1] + function_pi_sum_part)
            return function_pi_sum

def chi_squ(parameters, qsqu, pi_qsqu_avg, c_inverse, type):
    d_data_theory = np.array([np.array(pi_qsqu_avg) - np.array(func_pi(parameters, qsqu, type))])
    return d_data_theory.dot(c_inverse.dot(d_data_theory.T))[0][0]

def func_pi_der(parameters, x):
    m = len(parameters)
    der = [0] * m
    x = np.array(x)
    if m % 2 == 0:
        der[0] = [1.0] * len(x)
        for i in range(1, (m - 2) / 2 + 1):
            der[i] = - x * 1.0 / (parameters[i + (m - 2) / 2] + x)
            der[i + (m - 2) / 2] = x * (parameters[i] / (parameters[i + (m - 2) / 2] + x) ** 2.0)
        der[m - 1] = - x
        return np.array(der)
    elif m % 2 == 1:
        der[0] = [1.0] * len(x)
        for i in range(1, (m - 1) / 2 + 1):
            der[i] = - x * 1.0 / (parameters[i + (m - 1) / 2] + x)
            der[i + (m - 1) / 2] = x * (parameters[i] / (parameters[i + (m - 1) / 2] + x) ** 2.0)
        return np.array(der)

def chi_squ_der(parameters, qsqu, pi_qsqu_avg, c_inverse):
    return (np.array(func_pi_der(parameters, qsqu)).dot(
        c_inverse.dot(np.array([np.array(func_pi(parameters, qsqu)) - np.array(pi_qsqu_avg)]).T))).T[0]

def hessian_matrix(parameters, x, c_inverse, dis = 0):
    hessian = 2.0 * np.array(func_pi_der(parameters, x)).dot(c_inverse.dot(np.array(func_pi_der(parameters, x)).T))
    hessian_inv = np.linalg.inv(hessian)
    delta_ab = 2.0 * hessian_inv
    uncertainty = np.diag(delta_ab) ** (1.0 / 2.0)
    if dis != 0:
        print
        print 'hessian_matrix:'
        print hessian
        print
        print 'hessian_matrix_inverse:'
        print hessian_inv
        print
        print 'unit test:'
        print hessian.dot(hessian_inv)
        print
        print 'delta_ab:'
        print delta_ab
        print
    return hessian, delta_ab, uncertainty

def z_func(x, mu):
    return ((x ** 2.0 + 4.0 * mu ** 2.0 * x) ** (1.0 / 2.0) - x) / (2.0 * mu ** 2.0 * x)

def f_func(x, mu):
    return (mu ** 2.0 * x * (z_func(x, mu)) ** 3.0 * (1.0 - x * z_func(x, mu))) / (1.0 + mu ** 2.0 * x * (z_func(x, mu)) ** 2.0)

def comb_func(parameters, x, alpha, mu, type):
    if type == 'NOT DDS':
        return 4.0 * alpha ** 2.0 * f_func(x, mu) * (parameters[0] - func_pi(parameters, x, type))
    if type == 'DDS':
        return 4.0 * alpha ** 2.0 * f_func(x, mu) * func_pi(parameters, x, type)


class fitting:

    def __init__(self, qsqu, pi_qsqu, c_matrix, parameters_value, type):
        self.c_matrix = c_matrix
        self.qsqu = qsqu
        self.pi_qsqu = pi_qsqu
        self.c_inverse = func_inverse(self.c_matrix)

        '''min nelder-mead'''
        self.res = minimize(chi_squ, parameters_value, args=(self.qsqu, self.pi_qsqu, self.c_inverse, type,), method='nelder-mead', options={'maxiter': 1000000, 'maxfev': 1000000})
        '''min BFGS'''
        #self.res = minimize(chi_squ, (self.res.x), args=(self.qsqu, self.pi_qsqu, self.c_inverse,), method='BFGS', jac=chi_squ_der, options={'maxiter': 100000, 'gtol': 10 ** (-7)})

        '''hessian and delta_ab'''
        self.hessian, self.delta_ab, self.uncertainty = hessian_matrix((self.res.x), self.qsqu, self.c_inverse)

    def fitting_output(self):

        '''print results'''
        print
        print '########################################################################################################'
        print 'fitting points:'
        print len(self.qsqu)
        print 'q^2:'
        print self.qsqu
        print 'pi(q^2):'
        print self.pi_qsqu
        print 'err of piq2'
        print np.diag(self.c_matrix) ** (1.0 / 2.0)
        print 'chi2:'
        print self.res.fun
        # print 'num_points - num_parameters = DOF:'
        dof = len(self.qsqu) - len(self.res.x)
        # print num_points,'-', len(delitem), '-', len(parameters), '=', dof
        print 'chi2/DOF:'
        print self.res.fun / dof
        print 'parameters:'
        print self.res.x
        print 'uncertainty:'
        print self.uncertainty

#????????
class fit6496:
    def __init__(self):
        self.mu = 0.1056583715
        self.ainv = 2.241748633879781 #From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.alpha = 1.0 / 137.035999074
        self.charge_cor = 1.0
        self.fitting_points = [4, 5, 6, 7, 8, 9, 10]
        self.delitem = []
        self.patht = '/Users/tucheng/Desktop/fitting/cmatrix/l6496/21configs_AS'
        self.parameters_guess = np.array([ 0.09495962,  0.64574499, 0])
        self.cmatrix_filet = read_pi_qsqu_cmatrix_file(self.patht)
        self.cmatrix_filet.qsqu_all = self.ainv ** 2.0 * self.cmatrix_filet.qsqu_all
        self.cmatrix_filet.plots(10, typ='.r')
        plt.show()
    def run(self):
        for points_num in self.fitting_points:
            self.cmatrix_filet.split_qsqu_pi_qsqu_avg_c_matrix(points_num)
            fittingfunc = fitting(self.cmatrix_filet.qsqu, self.cmatrix_filet.pi_qsqu, self.cmatrix_filet.cmatrix,
                                  self.parameters_guess, 'DDS')
            fittingfunc.fitting_output()
            print 'a_mu (from 0 to 0.1):'
            print (quad(lambda x: comb_func((fittingfunc.res.x), x, self.alpha, self.mu, 'DDS'), 0, 0.1))[0]
            print 'g-2 (from 0.1 to 4*ainv):'
            print (quad(lambda x: comb_func((fittingfunc.res.x), x,  self.alpha, self.mu, 'DDS'), 0.1, 4*self.ainv))[0]


class FitJackknife:

    def __init__(self, jackknife_path, jackknife_label):
        self.mu = 0.1056583715
        self.ainv = 2.241748633879781 # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.alpha = 1.0 / 137.035999074
        self.charge_cor = 5./9.
        self.parameters_guess = np.array([0.09495962, 0.64574499])
        self.delitem = []
        self.jackknife_path = jackknife_path
        self.jackknife_label = jackknife_label

        self.get_jackknife_config_list()
        self.get_jackknife_files_full_path_list()
        self.numJackknife = len(self.jackknife_config_list)


        #self.patht = '/Users/tucheng/Desktop/fitting/cmatrix/l6496/21configs_AS'
        #self.cmatrix_filet = read_cmatrix_file(self.patht)
        #self.cmatrix_filet.qsqu_all = self.ainv ** 2.0 * self.cmatrix_filet.qsqu_all
        #self.cmatrix_filet.plots(10, typ='.r')
        #plt.show()
        return

    def get_jackknife_config_list(self):
        self.jackknife_config_list = []
        fileList = (src.func.walkfiles(self.jackknife_path))[1]
        for file in fileList:
            if self.jackknife_label in file:
                self.jackknife_config_list.append(file.split(".")[-1])
        return

    def get_jackknife_files_full_path_list(self):
        self.jackknife_files_full_path_list = []
        num = len(self.jackknife_config_list)
        for i in range(num):
            self.jackknife_files_full_path_list.append(self.jackknife_path + '/' + self.jackknife_label + '.' + self.jackknife_config_list[i])
        return

    def read_jackknife_cmatrix(self):
        self.jackknife_cmatrix_list = [] # List of Pointers to the Class 'read_cmatrix_file' for each Jackknife
        pi_qsqu_jackknife = []
        for one_jackknife_path in self.jackknife_files_full_path_list:
            one_jackknife_cmatrix = read_pi_qsqu_cmatrix_file(one_jackknife_path)
            self.jackknife_cmatrix_list.append(one_jackknife_cmatrix)

        for i in range(self.numJackknife):
            pi_qsqu_jackknife.append((self.jackknife_cmatrix_list[i]).pi_qsqu_all)
        pi_qsqu_jackknife = np.array(pi_qsqu_jackknife)
        self.qsqu_in_GeV = self.ainv ** 2.0 * (self.jackknife_cmatrix_list[0]).qsqu_all

        self.pi_qsqu_avg = np.mean(pi_qsqu_jackknife, axis=0)
        print 'Jackknife Average of Pi(q^2):'
        print np.array(self.pi_qsqu_avg)
        self.pi_qsqu_error = (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(pi_qsqu_jackknife, axis=0)
        print 'Jackknife Error of Pi(q^2)'
        print np.array(self.pi_qsqu_error)
        return

    def plot_pi_qsqu(self, num_of_plots, typ):
        plt.errorbar(self.qsqu_in_GeV[:num_of_plots], self.pi_qsqu_avg[:num_of_plots],
                     yerr=self.pi_qsqu_error[:num_of_plots], fmt=typ, elinewidth=1)
        return

    def fit_jackknife(self, num_fit_points):
        print '========================================================================================================'
        print 'Jackknife Fit for ' + str(num_fit_points) + ' Points'
        a_mu_List = []
        paramList = []

        for i in range(self.numJackknife):
            (self.jackknife_cmatrix_list[i]).split_qsqu_pi_qsqu_avg_c_matrix(num_fit_points)
            fittingfunc = fitting(self.qsqu_in_GeV[:num_fit_points], (self.jackknife_cmatrix_list[i]).pi_qsqu,
                                  (self.jackknife_cmatrix_list[i]).cmatrix, self.parameters_guess, 'DDS')
            fittingfunc.fitting_output()

            paramList.append(fittingfunc.res.x)
            aJack = (quad(lambda x: comb_func((fittingfunc.res.x), x, self.alpha, self.mu, 'DDS'), 0, np.inf))[0]
            a_mu_List.append(aJack)
        a_mu_List = np.array(a_mu_List)
        paramList = np.array(paramList)
        amu = self.charge_cor * np.mean(a_mu_List)
        self.paramAvg = np.mean(paramList, axis = 0)
        amuErr = self.charge_cor * (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(a_mu_List)
        paramErr = (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(paramList, axis = 0)
        print 'amu list'
        print a_mu_List
        print 'amu jackknife'
        print amu
        print 'amu jackknife error'
        print amuErr
        print 'parameter jackknife'
        print self.paramAvg
        print 'parameter error'
        print paramErr
        return

    def fit_jackknife_uncorr(self, numFittingPoints):
        print '========================================================='
        print 'fitting for ' + str(numFittingPoints) + ' points'
        a_mu_List = []
        paramList = []

        for i in range(self.numJackknife):
            (self.jackknife_cmatrix_list[i]).split_qsqu_pi_qsqu_avg_c_matrix(numFittingPoints)
            #fittingfunc = fitting(self.QsquGeV[:numFittingPoints], (self.JackCmatrixList[i]).pi_qsqu,
                                  #(self.JackCmatrixList[i]).cmatrix, self.parameters_guess, 'DDS')
            fittingfunc = fitting(self.qsqu_in_GeV[:numFittingPoints], (self.jackknife_cmatrix_list[i]).pi_qsqu,
                                  (self.jackknife_cmatrix_list[i]).diagCmatrix, self.parameters_guess, 'DDS')
            fittingfunc.fitting_output()

            paramList.append(fittingfunc.res.x)
            aJack = (quad(lambda x: comb_func((fittingfunc.res.x), x, self.alpha, self.mu, 'DDS'), 0, 2))[0]
            a_mu_List.append(aJack)
        a_mu_List = np.array(a_mu_List)
        paramList = np.array(paramList)
        amu = np.mean(a_mu_List)
        self.paramAvg = np.mean(paramList, axis = 0)
        amuErr = (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(a_mu_List)
        paramErr = (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(paramList, axis = 0)
        print 'amu list'
        print a_mu_List
        print 'amu jackknife'
        print amu
        print 'amu jackknife error'
        print amuErr
        print 'parameter jackknife'
        print self.paramAvg
        print 'parameter error'
        print paramErr
        return

    def plot_jackknife_fit(self, end_point, tp ='-r'):
        x = np.linspace(0, end_point, 500 * end_point)
        plt.plot(x, np.array(func_pi(self.paramAvg, x, 'DDS')), tp)
        return


class Fit4864Jackknife(FitJackknife):

    def __init__(self, jackknife_path, jackknife_label):
        self.mu = 0.1056583715
        self.ainv = 1.629278350515464               # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.alpha = 1.0 / 137.035999074
        self.charge_cor = 5.0/9.0                   # (1/3)^2 + (2/3)^2 u+d
        self.parameters_guess = np.array([0.09495962, 0.64574499,0.0])
        self.delitem = []
        self.jackknife_path = jackknife_path
        self.jackknife_label = jackknife_label

        self.get_jackknife_config_list()
        self.get_jackknife_files_full_path_list()
        self.numJackknife = len(self.jackknife_config_list)
        return


class Amu_T():

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        self.ainv = 1.629278350515464  # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.mu_lat = 0.1056583715 / self.ainv
        self.fac = 4. * self.alpha ** 2. * (5./9.) * 2. # the last '2' is integral over t 0 to \[Infinity]
        self.T_MAX = 64
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


def make_jackknife_dic(dic):  # need to improve

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


class AmuT_4864(Amu_T):

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        #self.ainv = 1.629278350515464  # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.a = 0.12224 / 0.1975      # GeV-1  # From Email
        self.ainv = 1. / self.a
        self.mu_lat = 0.1056583715 * self.a
        self.fac = 4. * self.alpha ** 2. * (5. / 9.) * 10. ** 10. * 2.
        self.T_MAX = 64
        return


class AmuT_6496(Amu_T):

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        self.a = 1. / 2.255453                # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.a = 0.08786 / 0.1975             # GeV-1   # From Email
        self.ainv = 1. / self.a
        self.mu_lat = 0.1056583715 * self.a
        self.fac = 4. * self.alpha ** 2. * (5./9.) * 10. ** 10. * 2.
        self.T_MAX = 96
        return

class AmuT_96192(Amu_T):

    def __init__(self):
        self.alpha = 1.0 / 137.035999074
        self.a = 1. / 2.255453                # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.a = 0.05662 / 0.1975             # GeV-1 (1 GeV-1 = 0.1975 fm) # From Email
        self.ainv = 1. / self.a
        self.mu_lat = 0.1056583715 * self.a
        self.fac = 4. * self.alpha ** 2. * (5./9.) * 10. ** 10. * 2.
        self.T_MAX = 192
        return




class Fit6496Jackknife:

    def __init__(self, jackknife_path, jackknife_label):
        self.mu = 0.1056583715
        self.ainv = 2.241748633879781 # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.alpha = 1.0 / 137.035999074
        self.charge_cor = 1.0
        self.parameters_guess = np.array([0.09495962, 0.64574499])
        self.delitem = []
        self.jackknife_path = jackknife_path
        self.jackknife_label = jackknife_label

        self.get_jackknife_config_list()
        self.get_jackknife_files_full_path_list()
        self.numJackknife = len(self.jackknife_config_list)


        #self.patht = '/Users/tucheng/Desktop/fitting/cmatrix/l6496/21configs_AS'
        #self.cmatrix_filet = read_cmatrix_file(self.patht)
        #self.cmatrix_filet.qsqu_all = self.ainv ** 2.0 * self.cmatrix_filet.qsqu_all
        #self.cmatrix_filet.plots(10, typ='.r')
        #plt.show()
        return

    def get_jackknife_config_list(self):
        self.jackknife_config_list = []
        fileList = (src.func.walkfiles(self.jackknife_path))[1]
        for file in fileList:
            if self.jackknife_label in file:
                self.jackknife_config_list.append(file.split(".")[-1])
        return

    def get_jackknife_files_full_path_list(self):
        self.jackknife_files_full_path_list = []
        num = len(self.jackknife_config_list)
        for i in range(num):
            self.jackknife_files_full_path_list.append(self.jackknife_path + '/' + self.jackknife_label + '.' + self.jackknife_config_list[i])
        return

    def read_jackknife_cmatrix(self):
        self.jackknife_cmatrix_list = [] # List of Pointers to the Class 'read_cmatrix_file' for each Jackknife
        pi_qsqu_jackknife = []
        for one_jackknife_path in self.jackknife_files_full_path_list:
            one_jackknife_cmatrix = read_pi_qsqu_cmatrix_file(one_jackknife_path)
            self.jackknife_cmatrix_list.append(one_jackknife_cmatrix)

        for i in range(self.numJackknife):
            pi_qsqu_jackknife.append((self.jackknife_cmatrix_list[i]).pi_qsqu_all)
        pi_qsqu_jackknife = np.array(pi_qsqu_jackknife)
        self.qsqu_in_GeV = self.ainv ** 2.0 * (self.jackknife_cmatrix_list[0]).qsqu_all

        self.pi_qsqu_avg = np.mean(pi_qsqu_jackknife, axis=0)
        print 'Jackknife Average of Pi(q^2):'
        print np.array(self.pi_qsqu_avg)
        self.pi_qsqu_error = (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(pi_qsqu_jackknife, axis=0)
        print 'Jackknife Error of Pi(q^2)'
        print np.array(self.pi_qsqu_error)
        return

    def plot_pi_qsqu(self, num_of_plots, typ):
        plt.errorbar(self.qsqu_in_GeV[:num_of_plots], self.pi_qsqu_avg[:num_of_plots],
                     yerr=self.pi_qsqu_error[:num_of_plots], fmt=typ, elinewidth=1)
        return

    def fit_jackknife(self, num_fit_points):
        print '========================================================================================================'
        print 'Jackknife Fit for ' + str(num_fit_points) + ' Points'
        a_mu_List = []
        paramList = []

        for i in range(self.numJackknife):
            (self.jackknife_cmatrix_list[i]).split_qsqu_pi_qsqu_avg_c_matrix(num_fit_points)
            fittingfunc = fitting(self.qsqu_in_GeV[:num_fit_points], (self.jackknife_cmatrix_list[i]).pi_qsqu,
                                  (self.jackknife_cmatrix_list[i]).cmatrix, self.parameters_guess, 'DDS')
            fittingfunc.fitting_output()

            paramList.append(fittingfunc.res.x)
            aJack = (quad(lambda x: comb_func((fittingfunc.res.x), x, self.alpha, self.mu, 'DDS'), 0, 2))[0]
            a_mu_List.append(aJack)
        a_mu_List = np.array(a_mu_List)
        paramList = np.array(paramList)
        amu = np.mean(a_mu_List)
        self.paramAvg = np.mean(paramList, axis = 0)
        amuErr = (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(a_mu_List)
        paramErr = (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(paramList, axis = 0)
        print 'amu list'
        print a_mu_List
        print 'amu jackknife'
        print amu
        print 'amu jackknife error'
        print amuErr
        print 'parameter jackknife'
        print self.paramAvg
        print 'parameter error'
        print paramErr
        return

    def fittingJackknifeUncorr(self, numFittingPoints):
        print '========================================================='
        print 'fitting for ' + str(numFittingPoints) + ' points'
        a_mu_List = []
        paramList = []

        for i in range(self.numJackknife):
            (self.jackknife_cmatrix_list[i]).split_qsqu_pi_qsqu_avg_c_matrix(numFittingPoints)
            #fittingfunc = fitting(self.QsquGeV[:numFittingPoints], (self.JackCmatrixList[i]).pi_qsqu,
                                  #(self.JackCmatrixList[i]).cmatrix, self.parameters_guess, 'DDS')
            fittingfunc = fitting(self.qsqu_in_GeV[:numFittingPoints], (self.jackknife_cmatrix_list[i]).pi_qsqu,
                                  (self.jackknife_cmatrix_list[i]).diagCmatrix, self.parameters_guess, 'DDS')
            fittingfunc.fitting_output()

            paramList.append(fittingfunc.res.x)
            aJack = (quad(lambda x: comb_func((fittingfunc.res.x), x, self.alpha, self.mu, 'DDS'), 0, 2))[0]
            a_mu_List.append(aJack)
        a_mu_List = np.array(a_mu_List)
        paramList = np.array(paramList)
        amu = np.mean(a_mu_List)
        self.paramAvg = np.mean(paramList, axis = 0)
        amuErr = (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(a_mu_List)
        paramErr = (self.numJackknife - 1.0) ** (1.0 / 2.0) * np.std(paramList, axis = 0)
        print 'amu list'
        print a_mu_List
        print 'amu jackknife'
        print amu
        print 'amu jackknife error'
        print amuErr
        print 'parameter jackknife'
        print self.paramAvg
        print 'parameter error'
        print paramErr
        return

    def plot_jackknife_fit(self, endPoint, tp ='-r'):
        x = np.linspace(0, endPoint, 500*endPoint)
        plt.plot(x, np.array(func_pi(self.paramAvg, x, 'DDS')), tp)
        return


class Plot6496PiQsquFromSingleFile:

    def __init__(self, file_path):
        self.mu = 0.1056583715
        self.ainv = 2.241748633879781 # From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.alpha = 1.0 / 137.035999074
        self.charge_cor = 1.0
        self.delitem = []

        self.filepath = file_path
        return

    def plot(self, num_plot):
        self.pi_qsqu_file = ReadPiQsqu(self.filepath)
        self.pi_qsqu_file.qsqu_all = self.ainv ** 2.0 * self.pi_qsqu_file.qsqu_all
        self.pi_qsqu_file.plots(num_plot, typ='.r')
        #plt.show()
        return

#????????????????
class fit48:
    def __init__(self):
        self.mu = 0.1056583715
        self.ainv = 3.486342381540
        self.alpha = 1.0 / 137.035999074
        self.charge_cor = 1.0
        self.fitting_points = [4, 5, 6, 7, 8, 9, 10]
        self.delitem = []
        self.patht = '/Users/tucheng/Desktop/Fitting/results/48c fitting/c_48_t'
        self.parameters_guess = np.array([ 0.09495962,  0.64574499, 0])
        self.cmatrix_filet = read_pi_qsqu_cmatrix_file(self.patht)
        self.cmatrix_filet.qsqu_all = self.ainv ** 2.0 * self.cmatrix_filet.qsqu_all
        self.cmatrix_filet.plots(10, typ='.r')
        plt.show()
    def run(self):
        for points_num in self.fitting_points:
            self.cmatrix_filet.split_qsqu_pi_qsqu_avg_c_matrix(points_num)
            fittingfunc = fitting(self.cmatrix_filet.qsqu, self.cmatrix_filet.pi_qsqu, self.cmatrix_filet.cmatrix,
                                  self.parameters_guess, 'DDS')
            fittingfunc.fitting_output()
            print 'g-2 (from 0.1 to 4*ainv):'
            print (quad(lambda x: comb_func((fittingfunc.res.x), x,  self.alpha, self.mu, 'DDS'), 0.1, 4*self.ainv))[0]

#????????????????
class fit96:
    def __init__(self):
        self.mu = 0.1056583715 # in GeV units
        self.ainv = 3.29316455696 #From https://arxiv.org/abs/1212.4768 TABLE 1: M_Pi*N_s/(M_Pi*L)
        self.alpha = 1.0 / 137.035999074
        self.charge_cor = 1.0
        self.fitting_points = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        self.delitem = []
        self.patht = '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_t'
        self.parameters_guess = np.array([ 0.09495962,  0.64574499, 0])
        self.cmatrix_filet = read_pi_qsqu_cmatrix_file(self.patht)
        self.cmatrix_filet.qsqu_all = self.ainv ** 2.0 * self.cmatrix_filet.qsqu_all
        self.cmatrix_filet.plots(10, typ='.r')
        plt.show()
    def run(self):
        for points_num in self.fitting_points:
            self.cmatrix_filet.split_qsqu_pi_qsqu_avg_c_matrix(points_num)
            fittingfunc = fitting(self.cmatrix_filet.qsqu, self.cmatrix_filet.pi_qsqu, self.cmatrix_filet.cmatrix,
                                  self.parameters_guess, 'DDS')
            fittingfunc.fitting_output()
            print 'a_mu (from 0 to 0.1):'
            print (quad(lambda x: comb_func((fittingfunc.res.x), x, self.alpha, self.mu, 'DDS'), 0, 0.1))[0]
            print 'a_mu (from 0.1 to 4*ainv):'
            print (quad(lambda x: comb_func((fittingfunc.res.x), x,  self.alpha, self.mu, 'DDS'), 0.1, 4*self.ainv))[0]


if __name__ == "__main__":


    '''
    start = Fit4864Jackknife('/Users/tucheng/Desktop/fitting/results/l48/EAL_Jackknife/', 'jackknife-26configs_Exact Sub AMA LMASUB LMA')
    start.read_jackknife_cmatrix()
    start.plot_pi_qsqu(30, '.r')

    start.fit_jackknife(4)
    start.plot_jackknife_fit(2.5, '-b')
    #start.fit_jackknife_uncorr(7)
    #start.plot_jackknife_fit(2.5, '-g')
    plt.show()
    '''

    '''
    plottest = Plot6496PiQsquFromSingleFile('/Users/tucheng/Desktop/fitting/data/test/out')
    plottest.plot(10)
    plottest = Plot6496PiQsquFromSingleFile('/Users/tucheng/Desktop/fitting/data/test-fnal/out-fnal')
    plottest.plot(10)
    plt.show()
    
    start = Fit6496Jackknife('/Users/tucheng/Desktop/fitting/cmatrix/l6496/jackknife', 'EAL.jackknife-21')
    start.read_jackknife_cmatrix()
    start.plot_pi_qsqu(15, '.r')
    #plt.show()
    #start.fittingJackknife(7)
    #start.plotJackknifeFit(2.5,'-b')
    #start.fittingJackknifeUncorr(7)
    #start.plotJackknifeFit(2.5, '-g')
    #print start.qsqu_all
    #print start.pi_qsqu_all
    
    plt.show()
    #start.run()
    
    
    #start = plot6496_AVG()
    '''

    amu = AmuT_6496()
    print 'l6496 w(t):', amu.get_w_t_list([i for i in range(96)])
    amu = AmuT_4864()
    print 'l4864 w(t):', amu.get_w_t_list([i for i in range(64)])

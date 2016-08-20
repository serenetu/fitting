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

class read_cmatrix_file:
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


    def split_qsqu_pi_qsqu_avg_c_matrix(self, points_num):
        self.qsqu = np.hsplit(self.qsqu_all, np.array([points_num]))[0]
        self.pi_qsqu = np.hsplit(self.pi_qsqu_all, np.array([points_num]))[0]
        self.cmatrix = np.hsplit(self.cmatrix_all, np.array([points_num]))[0]
        self.cmatrix = np.vsplit(self.cmatrix, np.array([points_num]))[0]

    def plots(self, points_num, typ = '.b'):
        self.split_qsqu_pi_qsqu_avg_c_matrix(points_num)
        plt.errorbar(self.qsqu, self.pi_qsqu, yerr=np.diag(self.cmatrix) ** (1.0 / 2.0), fmt=typ, elinewidth=1)
        #plt.show()

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
        print '######################################################################################################'
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

class fit64:
    def __init__(self):
        self.mu = 0.1056583715
        self.ainv = 3.486342381540
        self.alpha = 1.0 / 137.035999074
        self.charge_cor = 1.0
        self.fitting_points = [4, 5, 6, 7, 8, 9, 10]
        self.delitem = []
        self.patht = '/Users/tucheng/Desktop/Fitting/results/64c fitting/c_64_t'
        self.parameters_guess = np.array([ 0.09495962,  0.64574499, 0])
        self.cmatrix_filet = read_cmatrix_file(self.patht)
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

start = fit64()
start.run()

'''
#fit
mu = 0.1056583715
alpha = 1.0 / 137.035999074
charge_cor = 1.0
# charge_cor = 1.0 / 9.0  # strange
# charge_cor = 5.0 / 9.0  # light

fitting_points = [4, 5, 6, 7]
delitem = []
#lasp = 64
#lasp = 48
patht = '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_t'
pathx = '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_x'
pathy = '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_y'
pathz = '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_z'

parameters_init = np.array([ 0.09495962,  0.64574499, 0])

cmatrix_filet = read_cmatrix_file(patht)
cmatrix_filex = read_cmatrix_file(pathx)
cmatrix_filey = read_cmatrix_file(pathy)
cmatrix_filez = read_cmatrix_file(pathz)
#print cmatrix_file.cmatrix_all


cmatrix_filex.plots(10)
cmatrix_filey.plots(10)
cmatrix_filez.plots(10)
cmatrix_filet.plots(20, typ = '.r')
plt.show()



for points_num in fitting_points:
    cmatrix_filex.split_qsqu_pi_qsqu_avg_c_matrix(points_num)
    fittingfunc = fitting(cmatrix_filex.qsqu, cmatrix_filex.pi_qsqu, cmatrix_filex.cmatrix, parameters_init, 'DDS')
    fittingfunc.fitting_output()
'''
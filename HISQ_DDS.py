__author__ = 'SereneTu'

import numpy as np
# import sympy as sympy
from scipy.optimize import minimize
# from sympy import *
from scipy import integrate
from scipy.integrate import quadrature
from scipy.integrate import fixed_quad
import matplotlib.pyplot as plt
import random
import copy
import math
import sys
import src.func
import os


def order_by_Column(matrix, column):
    index = np.argsort(matrix[:, column])
    matrix1 = matrix[index]
    return matrix1

def order_by_two_Column(matrix, main_column, sec_column):
    matrix1 = matrix.tolist
    matrix1.sort(key=lambda l: (l[main_column], l[sec_column]))
    return np.array(matrix1)

# Combine two vectors to a matrix
# each vector is located to one COLUMN
def Mix_Vec_Vec(vector1, vector2):
    Vec1 = np.array([vector1])
    Vec2 = np.array([vector2])
    matrix = np.concatenate((Vec1.T, Vec2.T), axis=1)
    return matrix

def Mix_Vec_Matx(vector, matrix):
    Vec = np.array([vector])
    matrix_mix = np.concatenate((Vec.T, matrix), axis=1)
    return matrix_mix


def Split_Vec_Vec(matrix):
    matrix_array = np.array(matrix)
    vector1, vector2 = np.hsplit(matrix_array, [1])
    vector1 = vector1.T
    vector2 = vector2.T
    return vector1[0], vector2[0]

def Split_Vec_Matx(matrix):
    matrix_array = np.array(matrix)
    vector, matrix1 = np.hsplit(matrix_array, [1])
    vector = vector.T
    return vector[0], matrix1


def Avg_Same(vector, matrix):
    line_num = len(vector)
    num_same = 1
    vector_same = []
    matrix_same = []
    vector_same.append(vector[0])
    matrix_same.append(matrix[0])
    ii = 0

    for i in range(1, line_num):
        if vector[i] == vector[i - 1]:
            matrix_same[ii] = np.array(matrix_same[ii]) + np.array(matrix[i])
            num_same += 1
        else:
            matrix_same[ii] = matrix_same[ii] / num_same
            vector_same.append(vector[i])
            matrix_same.append(matrix[i])
            num_same = 1
            ii += 1
    matrix_same[-1] = matrix_same[-1] / num_same
    return np.array(vector_same), np.array(matrix_same)

def SortAvg(Vec, Mtx):
    vector = np.array(Vec)
    matrix = np.array(Mtx)

    index = np.argsort(vector)

    vector = vector[index]
    matrix = matrix[index]

    line_num = len(vector)
    num_same = 1
    vector_same = []
    matrix_same = []
    vector_same.append(vector[0])
    matrix_same.append(matrix[0])

    ii = 0

    for i in range(1, line_num):
        if vector[i] == vector[i - 1]:
            matrix_same[ii] = np.array(matrix_same[ii]) + np.array(matrix[i])
            num_same += 1
        else:
            matrix_same[ii] = matrix_same[ii] / num_same
            vector_same.append(vector[i])
            matrix_same.append(matrix[i])
            num_same = 1
            ii += 1
    matrix_same[-1] = matrix_same[-1] / num_same
    return np.array(vector_same), np.array(matrix_same)

class read_dds_files:
    def __init__(self, HISQpath, filename, l, Type):
        self.HISQpath = HISQpath
        self.filename = filename
        self.l = l[0]
        if Type == 'TimeSlice':
            self.TMAX = l[-1]
        elif Type == 'SpacialSlice':
            self.TMAX = l[0]

    # Sum elements in vector2 which T are the same
    def Add_Cmumu(self, vector1, vector2):
        line_num = len(vector1)
        num_same = 1
        vector1_same = []
        vector2_same = []
        vector1_same.append(vector1[0])
        vector2_same.append(vector2[0])
        ii = 0

        for i in range(1, line_num):
            if vector1[i] == vector1[i - 1]:
                vector2_same[ii] = np.array(vector2_same[ii]) + np.array(vector2[i])
                num_same += 1
            else:
                vector1_same.append(vector1[i])
                vector2_same.append(vector2[i])
                num_same = 1
                ii += 1
        return np.array(vector1_same), np.array(vector2_same)

    def output_Missing(self, path):

        file = open(path, 'a')


        file.write('Missing a line in a file:\n')
        file.write(str(self.fullpath) + '\n')
        file.write('Line #\n')
        file.write(str(self.linenum))
        file.close()

    def read_file_samemu_AllInOne(self, path, mu, slice):

        '''read a specific file line by line

        Args:
            data: the data read by using open()
            munu: index of mu to find C_mumu

        Returns:
            An array of sum of C_mumu from t = 0 to MAX

        '''

        self.linenum, self.Exact_num, self.appx_num, self.appx_avg_num = 0, 0, 0, 0

        self.Cmunu_Rel_Exact = []
        self.Cmunu_Img_Exact = []
        self.t_Exact = []
        self.Cmunu_Rel_appx = []
        self.Cmunu_Img_appx = []
        self.t_appx = []
        self.Cmunu_Rel_appx_avg = []
        self.Cmunu_Img_appx_avg = []
        self.t_appx_avg = []

        self.t_Exact_Add, self.Cmunu_Rel_Exact_Add = [], []
        self.t_appx_Add, self.Cmunu_Rel_appx_Add = [], []
        self.t_appx_avg_Add, self.Cmunu_Rel_appx_avg_Add = [], []

        self.Cmunu = []

        last_t = self.TMAX - 1
        last_mu = 3
        last_nu = 3
        count_num = 0

        temp_t, temp_Cmunu_Rel= [], []

        data = open(path)

        for line in data.readlines():
            self.linenum += 1
            if slice in line:
                if not 'CONTACT' in line:
                    linesplit = line.split()
                    t_line = int(linesplit[8])
                    mu_line = int(linesplit[4])
                    nu_line = int(linesplit[6])

                    if (mu_line in mu) & (nu_line in mu) & (mu_line == nu_line):
                        '''
                        if ((self.TMAX + t_line - last_t) % self.TMAX) != 1:
                            self.output_Missing('/Users/tucheng/Desktop/missing')
                        '''

                        if (t_line < last_t) & ((mu_line < last_mu) or (nu_line < last_nu)):
                            count_num += 1

                            temp_t = np.array(temp_t)
                            temp_Cmunu_Rel = np.array(temp_Cmunu_Rel)
                            index = np.argsort(temp_t)
                            if len(index) != 0:

                                temp_t = temp_t[index]
                                temp_Cmunu_Rel = temp_Cmunu_Rel[index]
                                # print temp_t[3*67],temp_t[3*67+1],temp_t[3*67+2], temp_Cmunu_Rel[3*67], temp_Cmunu_Rel[3*67 +1], temp_Cmunu_Rel[3*67 +2]
                                temp_t, temp_Cmunu_Rel = Avg_Same(temp_t, temp_Cmunu_Rel)
                                #print temp_t
                                #print temp_Cmunu_Rel
                                #print len(temp_Cmunu_Rel)
                                #print




                            # Read Exact
                            if count_num in range(2, 10):
                                self.t_Exact = temp_t
                                self.Cmunu_Rel_Exact.append(temp_Cmunu_Rel)

                            # Read Appx
                            elif count_num in range(10, 18):
                                self.t_appx = temp_t
                                self.Cmunu_Rel_appx.append(temp_Cmunu_Rel)
                            # Read Appx_avg
                            elif count_num >= 18:
                                self.t_appx_avg = temp_t
                                self.Cmunu_Rel_appx_avg.append(temp_Cmunu_Rel)



                            temp_t, temp_Cmunu_Rel= [], []
                            temp_t.append(int((int(linesplit[8]) - int(linesplit[2]) + self.TMAX) % self.TMAX))
                            temp_Cmunu_Rel.append(float(linesplit[9]))

                            #print temp_Cmunu_Rel
                            #print


                            #print temp_t
                            #print line
                            #print


                        else:
                            temp_t.append(int((int(linesplit[8]) - int(linesplit[2]) + self.TMAX) % self.TMAX))
                            temp_Cmunu_Rel.append(float(linesplit[9]))

                            #print temp_Cmunu_Rel
                            #print

                            #print temp_t
                            #print line
                            #print


                        last_t = t_line
                        last_mu = mu_line
                        last_nu = nu_line

        data.close()

        temp_t = np.array(temp_t)
        temp_Cmunu_Rel = np.array(temp_Cmunu_Rel)
        index = np.argsort(temp_t)
        if len(index) != 0:
            temp_t = temp_t[index]
            temp_Cmunu_Rel = temp_Cmunu_Rel[index]
            temp_t, temp_Cmunu_Rel = Avg_Same(temp_t, temp_Cmunu_Rel)
        self.t_appx_avg = temp_t
        self.Cmunu_Rel_appx_avg.append(temp_Cmunu_Rel)



        print  'num of Exact, appx_Exact, appx:'
        print len(self.Cmunu_Rel_Exact), len(self.Cmunu_Rel_appx), len(self.Cmunu_Rel_appx_avg)
        self.Cmunu_Rel_Exact = np.average(self.Cmunu_Rel_Exact, axis=0)
        #print self.Cmunu_Rel_Exact
        self.Cmunu_Rel_appx = np.average(self.Cmunu_Rel_appx, axis=0)
        self.Cmunu_Rel_appx_avg = np.average(self.Cmunu_Rel_appx_avg, axis=0)

        self.Cmunu = self.Cmunu_Rel_Exact - self.Cmunu_Rel_appx + self.Cmunu_Rel_appx_avg

        if len(self.Cmunu) != self.TMAX:
            print 'warning'
            exit()

    def ReadFileSamemu(self, path, mu, slice):
        data = open(path)
        t = []
        real = []
        for line in data.readlines():
            if slice in line:
                if not 'CONTACT' in line:
                    linesplit = line.split()
                    mu_line = int(linesplit[4])
                    nu_line = int(linesplit[6])

                    if (mu_line in mu) & (nu_line in mu) & (mu_line == nu_line):
                        t.append(int((int(linesplit[8]) - int(linesplit[2]) + self.TMAX) % self.TMAX))
                        real.append(float(linesplit[9]))
        #print t
        #print real[0], real[0+192], real[192*2], real[192*3], real[192*4], real[192*5], real[192*6], real[192*7]
        difft, self.Cmunu = SortAvg(t, real)

    def print_raw(self):
        print '##############################################################################'
        # print 'Exact_num, appx_num, appx_avg_num:'
        # print self.Exact_num, self.appx_num, self.appx_avg_num
        print 't_Exact:'
        print self.t_Exact
        print 't_appx'
        print self.t_appx
        print 't_appx_avg'
        print self.t_appx_avg
        print '##############################################################################'

    def print_T_Cmunu(self):
        print '##############################################################################'
        print 'T'
        print self.t_Exact_Add
        print 'Cmumu'
        print self.Cmunu
        print '##############################################################################'

    def check_linenum(self, path, num):
        data = open(path)
        data = data.readlines()
        return True if (len(data) == num) else False

    def read_allconfig_AllInOne(self):
        self.Cmunu_Configs = []
        self.conf_num_list = []
        i = 0
        folderlist = (src.func.walkfiles(self.HISQpath, prt=0))[0]
        for self.folder in folderlist:
            if ('l' + str(self.l)) in self.folder:
                self.conf_num = self.folder.split(".")[-1]
                self.fullpath = self.HISQpath + '/' + self.folder + '/' + self.filename + '.' + self.conf_num
                if os.path.exists(self.fullpath):
                    if os.path.getsize(self.fullpath) != 0:
                        self.conf_num_list.append(int(self.conf_num))
                        print 'Read File: ' + self.fullpath

                        self.read_file_samemu_AllInOne(self.fullpath, [0, 1, 2], 'VEC-CORR')\

                        self.Cmunu_Configs.append(self.Cmunu)

                        i += 1

        self.Cmunu_Configs = np.array(self.Cmunu_Configs)
        # print self.Cmunu_Configs



        self.conf_num_list, self.Cmunu_Configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.Cmunu_Configs), 0))

        self.conf_num_list, self.Cmunu_Configs = Avg_Same(self.conf_num_list, self.Cmunu_Configs)

        return self.Cmunu_Configs

    def ReadAllconfig(self, VEC = 'VEC-CORR'):
        self.Cmunu_Configs = []
        self.conf_num_list = []
        i = 0
        folderlist = (src.func.walkfiles(self.HISQpath, prt=0))[0]
        for self.folder in folderlist:
            if ('l' + str(self.l)) in self.folder:
                conf_num = self.folder.split(".")[-1]
                self.fullpath = self.HISQpath + '/' + self.folder + '/' + self.filename + '.' + conf_num
                if os.path.exists(self.fullpath):
                    if os.path.getsize(self.fullpath) != 0:
                        self.conf_num_list.append(int(conf_num))
                        print 'Read File: ' + self.fullpath

                        self.ReadFileSamemu(self.fullpath, [0, 1, 2], VEC)

                        self.Cmunu_Configs.append(self.Cmunu)

                        i += 1

        self.Cmunu_Configs = np.array(self.Cmunu_Configs)

        self.conf_num_list, self.Cmunu_Configs = SortAvg(self.conf_num_list, self.Cmunu_Configs)

        return self.Cmunu_Configs

    def read_allconfig_tslices(self, Line_Right):
        self.Cmunu_Configs = []
        self.conf_num_list = []
        i = 0
        folderlist = (src.func.walkfiles(self.HISQpath, prt=0))[0]
        for self.folder in folderlist:
            if ('l' + str(self.l)) in self.folder:
                self.conf_num = self.folder.split(".")[-1]
                self.fullpath = self.HISQpath + '/' + self.folder + '/' + self.filename + '.' + self.conf_num
                if os.path.exists(self.fullpath):
                    if os.path.getsize(self.fullpath) != 0:
                        print 'Read File: ' + self.fullpath
                        if self.check_linenum(self.fullpath, Line_Right):
                            self.conf_num_list.append(int(self.conf_num))

                            self.read_file_samemu_AllInOne(self.fullpath, [0, 1, 2], 'VEC-CORRt')

                            self.Cmunu_Configs.append(self.Cmunu)
                        else:
                            print 'not enough measurement'

                        i += 1

        self.Cmunu_Configs = np.array(self.Cmunu_Configs)
        # print self.Cmunu_Configs


        #print len(self.conf_num_list), len(self.Cmunu_Configs)
        self.conf_num_list, self.Cmunu_Configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.Cmunu_Configs), 0))

        self.conf_num_list, self.Cmunu_Configs = Avg_Same(self.conf_num_list, self.Cmunu_Configs)

        return self.Cmunu_Configs

    def read_allconfig_xyzslices(self, Line_Right):
        self.Cmunu_Configs = []
        self.conf_num_list = []
        i = 0
        folderlist = (src.func.walkfiles(self.HISQpath, prt=0))[0]
        for self.folder in folderlist:
            if ('l' + str(self.l)) in self.folder:
                self.conf_num = self.folder.split(".")[-1]
                self.fullpath = self.HISQpath + '/' + self.folder + '/' + self.filename + '.' + self.conf_num
                if os.path.exists(self.fullpath):
                    if os.path.getsize(self.fullpath) != 0:
                        print 'Read File: ' + self.fullpath

                        if self.check_linenum(self.fullpath, Line_Right):


                            self.conf_num_list.append(int(self.conf_num))
                            self.read_file_samemu_AllInOne(self.fullpath, [1, 2, 3], 'VEC-CORRx')
                            self.Cmunu_Configs.append(self.Cmunu)


                            self.conf_num_list.append(int(self.conf_num))
                            self.read_file_samemu_AllInOne(self.fullpath, [0, 2, 3], 'VEC-CORRy')
                            self.Cmunu_Configs.append(self.Cmunu)


                            self.conf_num_list.append(int(self.conf_num))
                            self.read_file_samemu_AllInOne(self.fullpath, [0, 1, 3], 'VEC-CORRz')
                            self.Cmunu_Configs.append(self.Cmunu)


                        else:
                            print 'not enough measurement'


                        i += 1

        self.Cmunu_Configs = np.array(self.Cmunu_Configs)
        self.conf_num_list, self.Cmunu_Configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.Cmunu_Configs), 0))
        self.conf_num_list, self.Cmunu_Configs = Avg_Same(self.conf_num_list, self.Cmunu_Configs)

        return self.Cmunu_Configs

    def read_allconfig_xslices(self, Line_Right):
        self.Cmunu_Configs = []
        self.conf_num_list = []
        i = 0
        folderlist = (src.func.walkfiles(self.HISQpath, prt=0))[0]
        for self.folder in folderlist:
            if ('l' + str(self.l)) in self.folder:
                self.conf_num = self.folder.split(".")[-1]
                self.fullpath = self.HISQpath + '/' + self.folder + '/' + self.filename + '.' + self.conf_num
                if os.path.exists(self.fullpath):
                    if os.path.getsize(self.fullpath) != 0:
                        print 'Read File: ' + self.fullpath

                        if self.check_linenum(self.fullpath, Line_Right):


                            self.conf_num_list.append(int(self.conf_num))
                            self.read_file_samemu_AllInOne(self.fullpath, [1, 2, 3], 'VEC-CORRx')
                            self.Cmunu_Configs.append(self.Cmunu)


                        else:
                            print 'not enough measurement'


                        i += 1

        self.Cmunu_Configs = np.array(self.Cmunu_Configs)
        self.conf_num_list, self.Cmunu_Configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.Cmunu_Configs), 0))
        self.conf_num_list, self.Cmunu_Configs = Avg_Same(self.conf_num_list, self.Cmunu_Configs)

        return self.Cmunu_Configs

    def read_allconfig_yslices(self, Line_Right):
        self.Cmunu_Configs = []
        self.conf_num_list = []
        i = 0
        folderlist = (src.func.walkfiles(self.HISQpath, prt=0))[0]
        for self.folder in folderlist:
            if ('l' + str(self.l)) in self.folder:
                self.conf_num = self.folder.split(".")[-1]
                self.fullpath = self.HISQpath + '/' + self.folder + '/' + self.filename + '.' + self.conf_num
                if os.path.exists(self.fullpath):
                    if os.path.getsize(self.fullpath) != 0:
                        print 'Read File: ' + self.fullpath

                        if self.check_linenum(self.fullpath, Line_Right):

                            self.conf_num_list.append(int(self.conf_num))
                            self.read_file_samemu_AllInOne(self.fullpath, [0, 2, 3], 'VEC-CORRy')
                            self.Cmunu_Configs.append(self.Cmunu)


                        else:
                            print 'not enough measurement'


                        i += 1

        self.Cmunu_Configs = np.array(self.Cmunu_Configs)
        self.conf_num_list, self.Cmunu_Configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.Cmunu_Configs), 0))
        self.conf_num_list, self.Cmunu_Configs = Avg_Same(self.conf_num_list, self.Cmunu_Configs)

        return self.Cmunu_Configs

    def read_allconfig_zslices(self, Line_Right):
        self.Cmunu_Configs = []
        self.conf_num_list = []
        i = 0
        folderlist = (src.func.walkfiles(self.HISQpath, prt=0))[0]
        for self.folder in folderlist:
            if ('l' + str(self.l)) in self.folder:
                self.conf_num = self.folder.split(".")[-1]
                self.fullpath = self.HISQpath + '/' + self.folder + '/' + self.filename + '.' + self.conf_num
                if os.path.exists(self.fullpath):
                    if os.path.getsize(self.fullpath) != 0:
                        print 'Read File: ' + self.fullpath

                        if self.check_linenum(self.fullpath, Line_Right):

                            self.conf_num_list.append(int(self.conf_num))
                            self.read_file_samemu_AllInOne(self.fullpath, [0, 1, 3], 'VEC-CORRz')
                            self.Cmunu_Configs.append(self.Cmunu)


                        else:
                            print 'not enough measurement'


                        i += 1

        self.Cmunu_Configs = np.array(self.Cmunu_Configs)
        self.conf_num_list, self.Cmunu_Configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.Cmunu_Configs), 0))
        self.conf_num_list, self.Cmunu_Configs = Avg_Same(self.conf_num_list, self.Cmunu_Configs)

        return self.Cmunu_Configs

def Pi_q(q, Cmunu, tmax):
    Pi = 0
    #tmax = len(Cmunu)
    for t in range(0, tmax / 2):
        Pi += ((math.cos(q * t) - 1.0) / q ** 2.0 + 1.0 / 2.0 * t ** 2.0) * Cmunu[t]
    for t in range(tmax / 2, tmax):
        Pi += ((math.cos(q * (t - tmax)) - 1.0) / q ** 2.0 + 1.0 / 2.0 * (t - tmax) ** 2.0) * Cmunu[t]
    return Pi

class Make_Pi_matrix:
    def __init__(self, qsqu, Cmunu_Configs, tmax):
        config, t = Cmunu_Configs.shape
        lenqsqu = len(qsqu)
        self.matrix = []
        for config_index in range(0, config):
            self.matrix.append([])
            for n in range(0, lenqsqu):
                self.matrix[config_index].append(Pi_q((qsqu[n]) ** (1.0 / 2.0), Cmunu_Configs[config_index], tmax))

class Make_qsqu:
    def __init__(self, maxnx, maxny, maxnz, maxnt, LL):
        self.L = LL
        self.qsqu = []
        for nx in range(0, maxnx):
            for ny in range(0, maxny):
                for nz in range(0, maxnz):
                    for nt in range(0, maxnt):
                        self.qsqu.append(self.func_qsqu(nx, ny, nz, nt))
        self.same_qsqu()

    def func_qsqu(self, nx, ny, nz, nt):
        return (2.0 * math.sin(math.pi * nx / self.L[0])) ** 2.0 + \
               (2.0 * math.sin(math.pi * ny / self.L[1])) ** 2.0 + \
               (2.0 * math.sin(math.pi * nz / self.L[2])) ** 2.0 + \
               (2.0 * math.sin(math.pi * nt / self.L[3])) ** 2.0

    def same_qsqu(self):
        self.qsqu = np.array(self.qsqu)
        index = np.argsort(self.qsqu)
        self.qsqu = self.qsqu[index]
        line_num = len(self.qsqu)
        qsqu_same = []
        qsqu_same.append(self.qsqu[0])

        for i in range(1, line_num):
            if abs(self.qsqu[i] - self.qsqu[i - 1])>10**(-7):
                qsqu_same.append(self.qsqu[i])
        qsqu_same = np.array(qsqu_same)
        qsqu_same = np.delete(qsqu_same, 0, 0)
        self.qsqu = qsqu_same

class Make_CovMatx:
    def __init__(self, qsqu, matrix, path):
        self.qsqu = qsqu
        self.savepath = path
        self.matrix_t = np.array(matrix).T
        self.t, self.n = (self.matrix_t).shape
        self.matrix_t_ave = np.mean(self.matrix_t, axis=1)
        self.c_matrix = np.zeros((self.t, self.t))
        for i in range(0, self.t):
            for j in range(0, self.t):
                # c_matrix[i, j] = 1.0 / (n - 1.0) * 1.0 / n * sum((matrix_t[i, :] - matrix_t_ave[i]) * (matrix_t[j, :] - matrix_t_ave[j]))
                configsum = 0
                for k in range(0, self.n):
                    configsum += (self.matrix_t[i, k] - self.matrix_t_ave[i]) * (self.matrix_t[j, k] - self.matrix_t_ave[j])
                self.c_matrix[i, j] = 1.0 / (self.n - 1.0) * 1.0 / self.n * configsum
        #return self.c_matrix

    def output(self):
        file = open(self.savepath, 'w')
        file.write('AVERAGES:\n')
        for i in range(0, self.t):
            file.write(str(self.qsqu[i]) + ' ' + str(self.matrix_t_ave[i]) + '\n')
        file.write('c_matrix:' + str(self.n) + '\n')
        for i in range(0, self.t):
            for j in range(0, self.t):
                file.write(str(self.qsqu[i]) + ' ' + str(self.qsqu[j]) + ' ' + str(self.c_matrix[i][j]) + '\n')
        file.close()


'''
path = '/Volumes/Seagate Backup Plus Drive/lqcdproj/gMinus2/blum/HISQ'
# L = [64, 64, 64, 144]
# L = [48, 48, 48, 144]
L = [96, 96, 96, 192]

missing = open('/Users/tucheng/Desktop/missing_96', 'w')
read = read_dds_files(path, 'vec_t_ama', [0, 1, 2], [0, 1, 2], L, 'TimeSlice')
read.read_allconfig()
qsqu = Make_qsqu(1, 1, 1, 30, L)
pi_matrix = Make_Pi_matrix(qsqu.qsqu, read.Cmunu_Configs)
Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/c_64').output()
missing.close()
'''

class do_48:
    def __init__(self):
        self.path = '/Volumes/Seagate Backup Plus Drive/lqcdproj/gMinus2/blum/HISQ'
        self.L = [48, 48, 48, 144]

    def runTimeSlices(self):
        print 'Reading TimeSlice'
        read_Exact = read_dds_files(self.path, 'vec_Exact', self.L, 'TimeSlice')
        read_Exact.ReadAllconfig()
        read_Sub = read_dds_files(self.path, 'vec_Sub', self.L, 'TimeSlice')
        read_Sub.ReadAllconfig()
        read_AMA = read_dds_files(self.path, 'vec_AMA', self.L, 'TimeSlice')
        read_AMA.ReadAllconfig()
        Cmunu_Configs = read_Exact.Cmunu_Configs - read_Sub.Cmunu_Configs + read_AMA.Cmunu_Configs
        qsqu = Make_qsqu(1, 1, 1, 30, self.L)
        pi_matrix = Make_Pi_matrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/48c fitting/c_48_t').output()

class do_64:
    def __init__(self):
        self.path = '/Volumes/Seagate Backup Plus Drive/lqcdproj/gMinus2/blum/HISQ'
        self.L = [64, 64, 64, 144]

    def runTimeSlices(self):
        print 'Reading TimeSlice'
        read_t = read_dds_files(self.path, 'vec', self.L, 'TimeSlice')
        read_t.read_allconfig_AllInOne()
        qsqu = Make_qsqu(1, 1, 1, 30, self.L)
        pi_matrix = Make_Pi_matrix(qsqu.qsqu, read_t.Cmunu_Configs, self.L[3])
        Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/64c fitting/c_64_t').output()


class do_96:
    def __init__(self):
        self.path = '/Volumes/Seagate Backup Plus Drive/lqcdproj/gMinus2/blum/HISQ'
        #self.path = '/Users/tucheng/Desktop'
        self.L = [96, 96, 96, 192]
        self.outputfile = '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_t'

    def runTimeSlicesHalf(self):
        print 'Reading TimeSlice'
        read_t = read_dds_files(self.path, 'vec_ama_xyzt', self.L, 'TimeSlice')
        read_t.read_allconfig_tslices(952320)
        qsqu = Make_qsqu(1, 1, 1, 30, self.L)
        pi_matrix = Make_Pi_matrix(qsqu.qsqu, read_t.Cmunu_Configs, self.L[3] / 2)
        Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_thalf').output()

    def runTimeSlices(self):
        print 'Reading TimeSlice'
        read_t = read_dds_files(self.path, 'vec_ama_xyzt', self.L, 'TimeSlice')
        read_t.read_allconfig_tslices(952320)
        print read_t.Cmunu_Configs[0]
        qsqu = Make_qsqu(1, 1, 1, 30, self.L)
        pi_matrix = Make_Pi_matrix(qsqu.qsqu, read_t.Cmunu_Configs, self.L[3])
        Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_t').output()

    def runXYZSlices(self):
        print 'Reading SpacialSlice'
        read_xyz = read_dds_files(self.path, 'vec_ama_xyzt', self.L, 'SpacialSlice')
        read_xyz.read_allconfig_xyzslices(952320)
        qsqu = Make_qsqu(30, 1, 1, 1, self.L)
        pi_matrix = Make_Pi_matrix(qsqu.qsqu, read_xyz.Cmunu_Configs, self.L[0])
        Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_xyz').output()

    def runXSlices(self):
        print 'Reading X SpacialSlice'
        read_x = read_dds_files(self.path, 'vec_ama_xyzt', self.L, 'SpacialSlice')
        read_x.read_allconfig_xslices(952320)
        qsqu = Make_qsqu(30, 1, 1, 1, self.L)
        pi_matrix = Make_Pi_matrix(qsqu.qsqu, read_x.Cmunu_Configs, self.L[0])
        Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_x').output()

    def runYSlices(self):
        print 'Reading Y SpacialSlice'
        read_y = read_dds_files(self.path, 'vec_ama_xyzt', self.L, 'SpacialSlice')
        read_y.read_allconfig_yslices(952320)
        qsqu = Make_qsqu(30, 1, 1, 1, self.L)
        pi_matrix = Make_Pi_matrix(qsqu.qsqu, read_y.Cmunu_Configs, self.L[0])
        Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_y').output()

    def runZSlices(self):
        print 'Reading Y SpacialSlice'
        read_z = read_dds_files(self.path, 'vec_ama_xyzt', self.L, 'SpacialSlice')
        read_z.read_allconfig_zslices(952320)
        qsqu = Make_qsqu(30, 1, 1, 1, self.L)
        pi_matrix = Make_Pi_matrix(qsqu.qsqu, read_z.Cmunu_Configs, self.L[0])
        Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_z').output()

    def runTimeSlicesFromSeperateFiles(self):
        print 'Reading TimeSlice'

        read_Exact = read_dds_files(self.path, 'vec_Exact', self.L, 'TimeSlice')
        read_Exact.ReadAllconfig('VEC-CORRt')
        read_Sub = read_dds_files(self.path, 'vec_Sub', self.L, 'TimeSlice')
        read_Sub.ReadAllconfig('VEC-CORRt')
        read_AMA = read_dds_files(self.path, 'vec_AMA', self.L, 'TimeSlice')
        read_AMA.ReadAllconfig('VEC-CORRt')
        Cmunu_Configs = read_Exact.Cmunu_Configs - read_Sub.Cmunu_Configs + read_AMA.Cmunu_Configs
        # print Cmunu_Configs[0]
        qsqu = Make_qsqu(1, 1, 1, self.L[3], self.L)
        # print qsqu.qsqu
        pi_matrix = Make_Pi_matrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()


start = do_96()
#start.runTimeSlices()
start.runTimeSlicesFromSeperateFiles()

#do_64()


#start = do_96()
#start.runTimeSlicesHalf()

#start = do_64()
#start.runTimeSlices()

#start = do_48()
#start.runTimeSlices()
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

# ????????????
def order_by_Column(matrix, column):
    index = np.argsort(matrix[:, column])
    matrix1 = matrix[index]
    return matrix1

# ?????????????
def order_by_two_Column(matrix, main_column, sec_column):
    matrix1 = matrix.tolist
    matrix1.sort(key=lambda l: (l[main_column], l[sec_column]))
    return np.array(matrix1)

# ?????????????
def Mix_Vec_Vec(vector1, vector2):
    Vec1 = np.array([vector1])
    Vec2 = np.array([vector2])
    matrix = np.concatenate((Vec1.T, Vec2.T), axis=1)
    return matrix

# ?????????????
def Mix_Vec_Matx(vector, matrix):
    Vec = np.array([vector])
    matrix_mix = np.concatenate((Vec.T, matrix), axis=1)
    return matrix_mix

# ?????????????
def Split_Vec_Vec(matrix):
    matrix_array = np.array(matrix)
    vector1, vector2 = np.hsplit(matrix_array, [1])
    vector1 = vector1.T
    vector2 = vector2.T
    return vector1[0], vector2[0]

# ?????????????
def Split_Vec_Matx(matrix):
    matrix_array = np.array(matrix)
    vector, matrix1 = np.hsplit(matrix_array, [1])
    vector = vector.T
    return vector[0], matrix1

# ??????????????
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


def sort_avg(Vec, Mtx):
    """Sort and Merge Vector and Matrix Depend On Values in Vector

    First Sort Vector and Matrix Depend On Values in Vector
    Then Deduplicate Vector and Take Average Related Lines in Matrix

    :param Vec: Vector
    :param Mtx: Matrix
    :return:    Deduplicated Vector and Matrix in numpy.array Format
    """

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


def pi_q(q, Cmunu, tmax):
    Pi = 0
    for t in range(0, tmax / 2):
        Pi += ((math.cos(q * t) - 1.0) / q ** 2.0 + 1.0 / 2.0 * t ** 2.0) * Cmunu[t]
    for t in range(tmax / 2, tmax):
        Pi += ((math.cos(q * (t - tmax)) - 1.0) / q ** 2.0 + 1.0 / 2.0 * (t - tmax) ** 2.0) * Cmunu[t]
    return Pi


class ReadCorrelationFunction:

    def __init__(self, HISQpath, filename, l, Type):
        self.HISQpath = HISQpath
        self.filename = filename
        self.l = l[0]
        self.lt = l[3]
        if Type == 'TimeSlice':
            self.TMAX = l[-1]
        elif Type == 'SpacialSlice':
            self.TMAX = l[0]
        return

    def read_allconfigs_samemu(self, VEC='VEC-CORRt'):
        """Read C_munu From All Proper Configure Folders

        C_munu Have The Same Value of mu and nu
        C_munu is saved into self.cmunu_configs and self.cmunu_configs_dic

        :param VEC: Label of The Lines Which Will Be Read
        """
        self.cmunu_configs = []
        self.cmunu_configs_dic = {}
        self.config_list = []
        folderlist = (src.func.walkfiles(self.HISQpath, prt=0))[0]
        for folder in folderlist:
            if ('l' + str(self.l)) + str(self.lt) in folder:
                conf_num = folder.split(".")[-1]
                fullpath = self.HISQpath + '/' + folder + '/' + self.filename + '.' + conf_num
                if os.path.exists(fullpath):
                    if os.path.getsize(fullpath) != 0:
                        self.config_list.append(int(conf_num))
                        print 'Read File: ' + fullpath
                        self.read_samemu(fullpath, [0, 1, 2], VEC)
                        self.cmunu_configs.append(self.__onefile_cmunu)

        self.cmunu_configs = np.array(self.cmunu_configs)
        self.config_list, self.cmunu_configs = sort_avg(self.config_list, self.cmunu_configs)

        for i in range(len(self.config_list)):
            self.cmunu_configs_dic[self.config_list[i]] = np.array(self.cmunu_configs[i])

        return

    def read_samemu(self, file_full_path, mu, slice):
        """Read C_munu From One File

        C_munu Have The Same Value of mu and nu
        C_munu is saved into self.__onefile_cmunu
        
        :param file_full_path: Full Path of the File
        :param mu: Value of mu That Will Be Read
        :param slice: Label of The Lines That Will Be Read
        """
        data = open(file_full_path)
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
        difft, self.__onefile_cmunu = sort_avg(t, real)

#????????????????????????
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

        self.__onefile_cmunu = []

        last_t = self.TMAX - 1
        last_mu = 3
        last_nu = 3
        count_num = 0

        temp_t, temp_Cmunu_Rel = [], []

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
                                # print temp_t
                                # print temp_Cmunu_Rel
                                # print len(temp_Cmunu_Rel)
                                # print

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

                            temp_t, temp_Cmunu_Rel = [], []
                            temp_t.append(int((int(linesplit[8]) - int(linesplit[2]) + self.TMAX) % self.TMAX))
                            temp_Cmunu_Rel.append(float(linesplit[9]))

                            # print temp_Cmunu_Rel
                            # print


                            # print temp_t
                            # print line
                            # print


                        else:
                            temp_t.append(int((int(linesplit[8]) - int(linesplit[2]) + self.TMAX) % self.TMAX))
                            temp_Cmunu_Rel.append(float(linesplit[9]))

                            # print temp_Cmunu_Rel
                            # print

                            # print temp_t
                            # print line
                            # print

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
        # print self.Cmunu_Rel_Exact
        self.Cmunu_Rel_appx = np.average(self.Cmunu_Rel_appx, axis=0)
        self.Cmunu_Rel_appx_avg = np.average(self.Cmunu_Rel_appx_avg, axis=0)

        self.__onefile_cmunu = self.Cmunu_Rel_Exact - self.Cmunu_Rel_appx + self.Cmunu_Rel_appx_avg

        if len(self.__onefile_cmunu) != self.TMAX:
            print 'warning'
            exit()

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
        print self.__onefile_cmunu
        print '##############################################################################'

    def check_linenum(self, path, num):
        data = open(path)
        data = data.readlines()
        return True if (len(data) == num) else False

    def read_allconfig_AllInOne(self):
        self.cmunu_configs = []
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

                        self.read_file_samemu_AllInOne(self.fullpath, [0, 1, 2], 'VEC-CORR')
                        self.cmunu_configs.append(self.__onefile_cmunu)

                        i += 1

        self.cmunu_configs = np.array(self.cmunu_configs)
        # print self.Cmunu_Configs



        self.conf_num_list, self.cmunu_configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.cmunu_configs), 0))

        self.conf_num_list, self.cmunu_configs = Avg_Same(self.conf_num_list, self.cmunu_configs)

        return self.cmunu_configs

    def CheckLenCmunu(self, Len):
        config_num = len(self.cmunu_configs)
        WrongIndex = []
        for i in range(0, config_num):
            if len(self.cmunu_configs[i]) != Len:
                WrongIndex.append(i)
                print 'Config ' + str(self.conf_num_list[i]) + ' has not enough t for Cmunu'
        return WrongIndex

    def read_allconfig_tslices(self, Line_Right):
        self.cmunu_configs = []
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

                            self.cmunu_configs.append(self.__onefile_cmunu)
                        else:
                            print 'not enough measurement'

                        i += 1

        self.cmunu_configs = np.array(self.cmunu_configs)
        # print self.Cmunu_Configs


        # print len(self.conf_num_list), len(self.Cmunu_Configs)
        self.conf_num_list, self.cmunu_configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.cmunu_configs), 0))

        self.conf_num_list, self.cmunu_configs = Avg_Same(self.conf_num_list, self.cmunu_configs)

        return self.cmunu_configs

    def read_allconfig_xyzslices(self, Line_Right):
        self.cmunu_configs = []
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
                            self.cmunu_configs.append(self.__onefile_cmunu)

                            self.conf_num_list.append(int(self.conf_num))
                            self.read_file_samemu_AllInOne(self.fullpath, [0, 2, 3], 'VEC-CORRy')
                            self.cmunu_configs.append(self.__onefile_cmunu)

                            self.conf_num_list.append(int(self.conf_num))
                            self.read_file_samemu_AllInOne(self.fullpath, [0, 1, 3], 'VEC-CORRz')
                            self.cmunu_configs.append(self.__onefile_cmunu)


                        else:
                            print 'not enough measurement'

                        i += 1

        self.cmunu_configs = np.array(self.cmunu_configs)
        self.conf_num_list, self.cmunu_configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.cmunu_configs), 0))
        self.conf_num_list, self.cmunu_configs = Avg_Same(self.conf_num_list, self.cmunu_configs)

        return self.cmunu_configs

    def read_allconfig_xslices(self, Line_Right):
        self.cmunu_configs = []
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
                            self.cmunu_configs.append(self.__onefile_cmunu)


                        else:
                            print 'not enough measurement'

                        i += 1

        self.cmunu_configs = np.array(self.cmunu_configs)
        self.conf_num_list, self.cmunu_configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.cmunu_configs), 0))
        self.conf_num_list, self.cmunu_configs = Avg_Same(self.conf_num_list, self.cmunu_configs)

        return self.cmunu_configs

    def read_allconfig_yslices(self, Line_Right):
        self.cmunu_configs = []
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
                            self.cmunu_configs.append(self.__onefile_cmunu)


                        else:
                            print 'not enough measurement'

                        i += 1

        self.cmunu_configs = np.array(self.cmunu_configs)
        self.conf_num_list, self.cmunu_configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.cmunu_configs), 0))
        self.conf_num_list, self.cmunu_configs = Avg_Same(self.conf_num_list, self.cmunu_configs)

        return self.cmunu_configs

    def read_allconfig_zslices(self, Line_Right):
        self.cmunu_configs = []
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
                            self.cmunu_configs.append(self.__onefile_cmunu)


                        else:
                            print 'not enough measurement'

                        i += 1

        self.cmunu_configs = np.array(self.cmunu_configs)
        self.conf_num_list, self.cmunu_configs = Split_Vec_Matx(
            order_by_Column(Mix_Vec_Matx(self.conf_num_list, self.cmunu_configs), 0))
        self.conf_num_list, self.cmunu_configs = Avg_Same(self.conf_num_list, self.cmunu_configs)

        return self.cmunu_configs


class MakePiMatrix:
    def __init__(self, qsqu, Cmunu_Configs, tmax):
        config, t = Cmunu_Configs.shape
        lenqsqu = len(qsqu)
        self.matrix = []
        for config_index in range(0, config):
            self.matrix.append([])
            for n in range(0, lenqsqu):
                self.matrix[config_index].append(pi_q((qsqu[n]) ** (1.0 / 2.0), Cmunu_Configs[config_index], tmax))


class MakeQsqu:
    def __init__(self, maxnx, maxny, maxnz, maxnt, LL):
        self.L = LL
        self.qsqu = []
        for nx in range(0, maxnx):
            for ny in range(0, maxny):
                for nz in range(0, maxnz):
                    for nt in range(0, maxnt):
                        self.qsqu.append(self.func_qsqu(nx, ny, nz, nt))
        self.same_qsqu()
        return

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
            if abs(self.qsqu[i] - self.qsqu[i - 1]) > 10 ** (-7):
                qsqu_same.append(self.qsqu[i])
        qsqu_same = np.array(qsqu_same)
        qsqu_same = np.delete(qsqu_same, 0, 0)
        self.qsqu = qsqu_same
        return


class MakeCovarianceMatrix:

    def __init__(self, qsqu, matrix, path):
        self.qsqu = qsqu
        self.savepath = path
        self.matrix_t = np.array(matrix).T
        self.t, self.n = (self.matrix_t).shape
        self.matrix_t_ave = np.mean(self.matrix_t, axis=1)
        self.c_matrix = np.zeros((self.t, self.t))
        if self.n != 1:

            for i in range(0, self.t):
                for j in range(0, self.t):
                    configsum = 0
                    for k in range(0, self.n):
                        configsum += (self.matrix_t[i, k] - self.matrix_t_ave[i]) * (
                        self.matrix_t[j, k] - self.matrix_t_ave[j])
                    self.c_matrix[i, j] = 1.0 / (self.n - 1.0) * 1.0 / self.n * configsum
        return

    def output(self):
        file = open(self.savepath, 'w')

        # Write averages:
        file.write('AVERAGES:\n')
        for i in range(0, self.t):
            file.write(str(self.qsqu[i]) + ' ' + str(self.matrix_t_ave[i]) + '\n')

        # Write c_matrix
        if self.n != 1:
            file.write('c_matrix:' + str(self.n) + '\n')
            for i in range(0, self.t):
                for j in range(0, self.t):
                    file.write(str(self.qsqu[i]) + ' ' + str(self.qsqu[j]) + ' ' + str(self.c_matrix[i][j]) + '\n')

        file.close()
        return


class do_48:
    def __init__(self):
        self.path = '/Volumes/Seagate Backup Plus Drive/lqcdproj/gMinus2/blum/HISQ'
        self.L = [48, 48, 48, 144]

    def runTimeSlices(self):
        print 'Reading TimeSlice'
        read_Exact = ReadCorrelationFunction(self.path, 'vec_Exact', self.L, 'TimeSlice')
        read_Exact.read_allconfigs_samemu()
        read_Sub = ReadCorrelationFunction(self.path, 'vec_Sub', self.L, 'TimeSlice')
        read_Sub.read_allconfigs_samemu()
        read_AMA = ReadCorrelationFunction(self.path, 'vec_AMA', self.L, 'TimeSlice')
        read_AMA.read_allconfigs_samemu()
        Cmunu_Configs = read_Exact.cmunu_configs - read_Sub.cmunu_configs + read_AMA.cmunu_configs
        qsqu = MakeQsqu(1, 1, 1, 30, self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/48c fitting/c_48_t').output()


class do_6496:
    def __init__(self, data_path):
        self.path = data_path
        self.L = [64, 64, 64, 96]
        self.Vol = self.L[0] * self.L[1] * self.L[2] * self.L[3]
        self.numPointSrc = 2.0 * 2.0 * 2.0 * 2.0

        self.read_E = None
        self.read_ES = None
        self.read_A = None
        self.read_L = None
        self.read_LS = None

        self.CmunuFactor = 1.0 / 4.0
        return

    def read_allconfig_Exact_Sub_AMA_LMASUB_LMA(self, Ename, ESname, Aname, Lname, LSname):

        self.read_E = ReadCorrelationFunction(self.path, Ename, self.L, 'TimeSlice')
        self.read_E.read_allconfigs_samemu('VEC-CORRt')
        print '========================================================================================================'
        print 'Read Exact Config List'
        print (self.read_E.config_list).tolist()
        print '========================================================================================================'

        self.read_ES = ReadCorrelationFunction(self.path, ESname, self.L, 'TimeSlice')
        self.read_ES.read_allconfigs_samemu('VEC-CORRt')
        print '========================================================================================================'
        print 'Read Sub Config List'
        print (self.read_ES.config_list).tolist()
        print '========================================================================================================'

        self.read_A = ReadCorrelationFunction(self.path, Aname, self.L, 'TimeSlice')
        self.read_A.read_allconfigs_samemu('VEC-CORRt')
        print '========================================================================================================'
        print 'Read AMA Config List'
        print (self.read_A.config_list).tolist()
        print '========================================================================================================'

        self.read_L = ReadCorrelationFunction(self.path, Lname, self.L, 'TimeSlice')
        self.read_L.read_allconfigs_samemu('VEC-CORRt')
        for i in self.read_L.cmunu_configs_dic:
            self.read_L.cmunu_configs_dic[i] = self.read_L.cmunu_configs_dic[i] / self.Vol
        print '========================================================================================================'
        print 'Read LMA Config List'
        print (self.read_L.config_list).tolist()
        print '========================================================================================================'

        self.read_LS = ReadCorrelationFunction(self.path, LSname, self.L, 'TimeSlice')
        self.read_LS.read_allconfigs_samemu('VEC-CORRt')
        for i in self.read_LS.cmunu_configs_dic:
            self.read_LS.cmunu_configs_dic[i] = self.read_LS.cmunu_configs_dic[i] / self.numPointSrc
        print '========================================================================================================'
        print 'Read LMASUB Config List'
        print (self.read_LS.config_list).tolist()
        print '========================================================================================================'

        return

    def read_allconfig_LMA(self, Lname):

        self.read_L = ReadCorrelationFunction(self.path, Lname, self.L, 'TimeSlice')
        self.read_L.read_allconfigs_samemu('VEC-CORRt')
        for i in self.read_L.cmunu_configs_dic:
            self.read_L.cmunu_configs_dic[i] = self.read_L.cmunu_configs_dic[i] / self.Vol
        print '========================================================================================================'
        print 'Read LMA Config List'
        print (self.read_L.config_list).tolist()
        print '========================================================================================================'

        return

    def read_allconfig_Exact_Sub_AMA(self, Ename, ESname, Aname):

        self.read_E = ReadCorrelationFunction(self.path, Ename, self.L, 'TimeSlice')
        self.read_E.read_allconfigs_samemu('VEC-CORRt')
        print '========================================================================================================'
        print 'Read Exact Config List'
        print (self.read_E.config_list).tolist()
        print '========================================================================================================'

        self.read_ES = ReadCorrelationFunction(self.path, ESname, self.L, 'TimeSlice')
        self.read_ES.read_allconfigs_samemu('VEC-CORRt')
        print '========================================================================================================'
        print 'Read Sub Config List'
        print (self.read_ES.config_list).tolist()
        print '========================================================================================================'

        self.read_A = ReadCorrelationFunction(self.path, Aname, self.L, 'TimeSlice')
        self.read_A.read_allconfigs_samemu('VEC-CORRt')
        print '========================================================================================================'
        print 'AMA Config List'
        print (self.read_A.config_list).tolist()
        print '========================================================================================================'

        return

    def make_pi_qsqu_covmatx_Exact_Sub_AMA_LMASUB_LMA(self, Ename, ESname, Aname, Lname, LSname, configList, out_dir):

        print 'Make Pi(q^2) and Covariance Matrix For Exact Sub AMA LMASUB LMA'
        Cmunu_Configs = []

        for i in configList:
            Cmunu_Configs.append(self.read_E.cmunu_configs_dic[i] - self.read_ES.cmunu_configs_dic[i] +
                                 self.read_A.cmunu_configs_dic[i] - self.read_LS.cmunu_configs_dic[i] +
                                 self.read_L.cmunu_configs_dic[i])
        Cmunu_Configs = np.array(Cmunu_Configs)

        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs * self.CmunuFactor, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, out_dir).output()
        return

    def make_pi_qsqu_covmatx_Exact_Sub_AMA(self, Ename, ESname, Aname, configList, out_dir):
        print 'Make Pi(q^2) and Covariance Matrix For Exact Sub AMA'

        Cmunu_Configs = []

        for i in configList:
            Cmunu_Configs.append(self.read_E.cmunu_configs_dic[i] - self.read_ES.cmunu_configs_dic[i] +
                                 self.read_A.cmunu_configs_dic[i])
        Cmunu_Configs = np.array(Cmunu_Configs)

        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs * self.CmunuFactor, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, out_dir).output()
        return

    def make_pi_qsqu_covmatx_LMA(self, configList, out_dir):
        print 'Make Pi(q^2) and Covariance Matrix For LMA'

        Cmunu_Configs = []

        for i in configList:
            Cmunu_Configs.append(self.read_L.cmunu_configs_dic[i])
        Cmunu_Configs = np.array(Cmunu_Configs)

        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs * self.CmunuFactor, self.L[3])
        make_covmtrx = MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, out_dir)
        make_covmtrx.output()
        plt.errorbar(range(0, self.L[3]), self.read_L.cmunu_configs[0][0:], fmt='.', elinewidth=1)
        plt.yscale('log', nonposy='clip')
        #plt.show()
        return




# ?????????????????????
    def runTimeSlices(self):
        print 'Reading TimeSlice'
        read_t = ReadCorrelationFunction(self.path, 'vec', self.L, 'TimeSlice')
        read_t.read_allconfig_AllInOne()
        qsqu = MakeQsqu(1, 1, 1, 30, self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, read_t.cmunu_configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/64c fitting/c_64_t').output()

    def runTimeSlicesFromSeparateFiles_ESA(self):
        print 'Reading TimeSlice From Separate Exact Aub AMA files'

        read_Exact = ReadCorrelationFunction(self.path, 'vec_Exact', self.L, 'TimeSlice')
        read_Exact.read_allconfigs_samemu('VEC-CORRt')

        read_Sub = ReadCorrelationFunction(self.path, 'vec_Sub', self.L, 'TimeSlice')
        read_Sub.read_allconfigs_samemu('VEC-CORRt')

        read_AMA = ReadCorrelationFunction(self.path, 'vec_AMA', self.L, 'TimeSlice')
        read_AMA.read_allconfigs_samemu('VEC-CORRt')

        Cmunu_Configs = read_Exact.cmunu_configs - read_Sub.cmunu_configs + read_AMA.cmunu_configs

        print 'E-ES+AMA'
        print (read_Exact.cmunu_configs - read_Sub.cmunu_configs + read_AMA.cmunu_configs).tolist()

        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlicesFromSeparateFiles_EAL(self):
        print 'Reading TimeSlice From Separate Exact Sub AMA LMA LMASUB files'

        read_Exact = ReadCorrelationFunction(self.path, 'vec_Exact', self.L, 'TimeSlice')
        read_Exact.read_allconfigs_samemu('VEC-CORRt')

        read_Sub = ReadCorrelationFunction(self.path, 'vec_Sub', self.L, 'TimeSlice')
        read_Sub.read_allconfigs_samemu('VEC-CORRt')

        read_AMA = ReadCorrelationFunction(self.path, 'vec_AMA', self.L, 'TimeSlice')
        read_AMA.read_allconfigs_samemu('VEC-CORRt')

        read_LMASUB = ReadCorrelationFunction(self.path, 'vec_LMASUB', self.L, 'TimeSlice')
        read_LMASUB.read_allconfigs_samemu('VEC-CORRt')

        read_LMA = ReadCorrelationFunction(self.path, 'vec_LMA', self.L, 'TimeSlice')
        read_LMA.read_allconfigs_samemu('VEC-CORRt')
        read_LMA.cmunu_configs = read_LMA.cmunu_configs / (self.L[0] * self.L[1] * self.L[2] * self.L[3])

        Cmunu_Configs = read_Exact.cmunu_configs - read_Sub.cmunu_configs + read_AMA.cmunu_configs - read_LMASUB.cmunu_configs + read_LMA.cmunu_configs

        print 'E-ES'
        print (read_Exact.cmunu_configs - read_Sub.cmunu_configs).tolist()
        print 'AMA'
        print (read_AMA.cmunu_configs).tolist()
        print 'E-ES+AMA'
        print (read_Exact.cmunu_configs - read_Sub.cmunu_configs + read_AMA.cmunu_configs).tolist()
        print 'LMA'
        print (read_LMA.cmunu_configs).tolist()
        print 'LMA-LMASUB'
        print (read_LMA.cmunu_configs - read_LMASUB.cmunu_configs).tolist()
        print 'E-ES+AMA-LMASUB+LMA'
        print Cmunu_Configs.tolist()

        # plt.errorbar(range(0, self.L[3]), (read_Exact.Cmunu_Configs - read_Sub.Cmunu_Configs + read_AMA.Cmunu_Configs)[0][0:], fmt='.b', elinewidth=1)
        # plt.errorbar(range(0, self.L[3]), (read_AMA.Cmunu_Configs)[0][:], fmt='.r', elinewidth=1)
        # plt.errorbar(range(0, self.L[3]), Cmunu_Configs[0][0:], fmt='.g', elinewidth=1)
        plt.errorbar(range(0, self.L[3]), read_LMA.cmunu_configs[0][0:], fmt='.b', elinewidth=1)
        plt.errorbar(range(0, self.L[3]), read_LMASUB.cmunu_configs[0][0:], fmt='.r', elinewidth=1)
        plt.yscale('log', nonposy='clip')
        # plt.show()



        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlicesFromSeparateFiles_L_LS(self, Lname, LSname, Lplt, LSplt):
        print 'Reading TimeSlice From LMA and LMASUB files'

        read_L = ReadCorrelationFunction(self.path, Lname, self.L, 'TimeSlice')
        read_L.read_allconfigs_samemu('VEC-CORRt')
        read_L.cmunu_configs = read_L.cmunu_configs / (self.L[0] * self.L[1] * self.L[2] * self.L[3])

        read_LS = ReadCorrelationFunction(self.path, LSname, self.L, 'TimeSlice')
        read_LS.read_allconfigs_samemu('VEC-CORRt')
        print read_LS.cmunu_configs
        read_LS.cmunu_configs = read_LS.cmunu_configs / ((self.L[0] / 32.0) ** 3.0 * (self.L[3] / 48.0))

        Cmunu_Configs = read_L.cmunu_configs - read_LS.cmunu_configs

        print 'L'
        print (read_L.cmunu_configs).tolist()

        print 'LS'
        print (read_LS.cmunu_configs).tolist()

        print 'L-LS'
        print (read_L.cmunu_configs - read_LS.cmunu_configs).tolist()

        # plt.errorbar(range(0, self.L[3]), (read_Exact.Cmunu_Configs - read_Sub.Cmunu_Configs + read_AMA.Cmunu_Configs)[0][0:], fmt='.b', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), (read_AMA.Cmunu_Configs)[0][1:], fmt='.r', elinewidth=1)
        plt.errorbar(range(1, self.L[3]), read_L.cmunu_configs[0][1:], fmt=Lplt, elinewidth=1)
        plt.errorbar(range(1, self.L[3]), read_LS.cmunu_configs[0][1:], fmt=LSplt, elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), Cmunu_Configs[0][1:], fmt='.r', elinewidth=1)

        # plt.errorbar(range(1, self.L[3]), read_L.Cmunu_Configs[0][1:], fmt='.b', elinewidth=1)
        plt.show()

        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlicesFromSeparateFiles_LS(self):
        print 'Reading TimeSlice From LS files'

        read_LMASUB = ReadCorrelationFunction(self.path, 'vec_LS', self.L, 'TimeSlice')
        read_LMASUB.read_allconfigs_samemu('VEC-CORRt')
        read_LMASUB.cmunu_configs = read_LMASUB.cmunu_configs / ((self.L[0] / 32.0) ** 3.0 * self.L[3] / 48.0)

        Cmunu_Configs = read_LMASUB.cmunu_configs

        print 'LS'
        print (read_LMASUB.cmunu_configs).tolist()

        # plt.errorbar(range(0, self.L[3]), (read_Exact.Cmunu_Configs - read_Sub.Cmunu_Configs + read_AMA.Cmunu_Configs)[0][0:], fmt='.b', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), (read_AMA.Cmunu_Configs)[0][1:], fmt='.r', elinewidth=1)
        plt.errorbar(range(1, self.L[3]), read_LMASUB.cmunu_configs[0][1:], fmt='.r', elinewidth=1)

        # plt.errorbar(range(1, self.L[3]), read_LMA.Cmunu_Configs[0][1:], fmt='.b', elinewidth=1)
        # plt.show()




        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlicesFromSeparateFiles_L_LS_LOWMODEAPPROX(self, Lname, LWAPPRXname, Lplt, LWAPPRXplt):
        print 'Reading TimeSlice From LMA and LMASUB_LOWMODEAPPROX files'

        read_L = ReadCorrelationFunction(self.path, Lname, self.L, 'TimeSlice')
        read_L.read_allconfigs_samemu('VEC-CORRt')
        read_L.cmunu_configs = read_L.cmunu_configs / (self.L[0] * self.L[1] * self.L[2] * self.L[3])

        read_LS_LOWMODEAPPROX = ReadCorrelationFunction(self.path, LWAPPRXname, self.L, 'TimeSlice')
        read_LS_LOWMODEAPPROX.read_allconfigs_samemu('VEC-CORRt')
        read_LS_LOWMODEAPPROX.cmunu_configs = read_LS_LOWMODEAPPROX.cmunu_configs

        Cmunu_Configs = read_L.cmunu_configs - read_LS_LOWMODEAPPROX.cmunu_configs

        print 'L'
        print (read_L.cmunu_configs).tolist()

        print 'LS_LOWMODEAPPROX'
        print (read_LS_LOWMODEAPPROX.cmunu_configs).tolist()

        # plt.errorbar(range(0, self.L[3]), (read_Exact.Cmunu_Configs - read_Sub.Cmunu_Configs + read_AMA.Cmunu_Configs)[0][0:], fmt='.b', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), (read_AMA.Cmunu_Configs)[0][1:], fmt='.r', elinewidth=1)
        plt.errorbar(range(1, self.L[3]), read_L.cmunu_configs[0][1:], fmt=Lplt, elinewidth=1)
        plt.errorbar(range(1, self.L[3]), read_LS_LOWMODEAPPROX.cmunu_configs[0][1:], fmt=LWAPPRXplt, elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), Cmunu_Configs[0][1:], fmt='.r', elinewidth=1)

        # plt.errorbar(range(1, self.L[3]), read_L.Cmunu_Configs[0][1:], fmt='.b', elinewidth=1)
        # plt.ylim(-0.0004, 0.0010)
        # plt.show()




        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlices_EAL_dic_selectedT(self, configList, numT, fileDir):
        print 'Reading TimeSlice From Separate Exact Sub AMA LMA LMASUB files'

        Cmunu_Configs = []

        for i in configList:
            Cmunu_Configs.append(self.read_E.Cmunu_Configs_dic[i] - self.read_ES.Cmunu_Configs_dic[i] +
                                 self.read_A.Cmunu_Configs_dic[i] - self.read_LS.Cmunu_Configs_dic[i] +
                                 self.read_L.Cmunu_Configs_dic[i])
        Cmunu_Configs = np.array(Cmunu_Configs)

        # manage the middle Cmunu(T)
        for i in range(len(configList)):
            left = Cmunu_Configs[i][:(self.L[3] / 2 - numT / 2)].tolist()
            middle = np.zeros(numT).tolist()
            right = Cmunu_Configs[i][(self.L[3] / 2 + numT / 2):].tolist()
            Cmunu_Configs[i] = np.array(left + middle + right)

        Cmunu_Configs = np.array(Cmunu_Configs)

        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs * self.CmunuFactor, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, fileDir).output()
        return

    def runTimeSlices_LS_dic(self, LSname, configList, fileDir):
        print 'Reading TimeSlice From LMASUB files'

        Cmunu_Configs = []

        read_LS = ReadCorrelationFunction(self.path, LSname, self.L, 'TimeSlice')
        read_LS.read_allconfigs_samemu('VEC-CORRt')
        for i in read_LS.cmunu_configs_dic:
            # print read_LS.Cmunu_Configs_dic[i]
            # print self.numPointSrc
            read_LS.cmunu_configs_dic[i] = read_LS.cmunu_configs_dic[i] / self.numPointSrc

        # print 'configs List for LMASUB:'
        # print list(configList)

        for i in configList:
            Cmunu_Configs.append(read_LS.cmunu_configs_dic[i])
        Cmunu_Configs = np.array(Cmunu_Configs)

        # print 'L'
        # print (read_L.Cmunu_Configs ).tolist()

        # print 'LS'
        # print (Cmunu_Configs ).tolist()

        # print 'L-LS'
        # print (read_L.Cmunu_Configs - read_LS.Cmunu_Configs).tolist()


        # plt.errorbar(range(0, self.L[3]), (read_Exact.Cmunu_Configs - read_Sub.Cmunu_Configs + read_AMA.Cmunu_Configs)[0][0:], fmt='.b', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), (read_AMA.Cmunu_Configs)[0][1:], fmt='.r', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), read_L.Cmunu_Configs[0][1:], fmt=Lplt, elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), Cmunu_Configs[0][1:], fmt='.r', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), Cmunu_Configs[0][1:], fmt='.r', elinewidth=1)

        # plt.errorbar(range(1, self.L[3]), read_L.Cmunu_Configs[0][1:], fmt='.b', elinewidth=1)
        # plt.show()


        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs * self.CmunuFactor, self.L[3])
        self.outputfile = fileDir
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlices_EL_dic(self, Ename, ESname, Lname, configList, fileDir):
        print 'Reading TimeSlice From Separate Exact Sub LMA files'

        Cmunu_Configs = []

        read_E = ReadCorrelationFunction(self.path, Ename, self.L, 'TimeSlice')
        read_E.read_allconfigs_samemu('VEC-CORRt')

        read_ES = ReadCorrelationFunction(self.path, ESname, self.L, 'TimeSlice')
        read_ES.read_allconfigs_samemu('VEC-CORRt')

        read_L = ReadCorrelationFunction(self.path, Lname, self.L, 'TimeSlice')
        read_L.read_allconfigs_samemu('VEC-CORRt')
        for i in read_L.cmunu_configs_dic:
            read_L.cmunu_configs_dic[i] = read_L.cmunu_configs_dic[i] / self.Vol

        # print 'Configs Read:'
        # print list(configList)

        for i in configList:
            Cmunu_Configs.append(read_E.cmunu_configs_dic[i] - read_ES.cmunu_configs_dic[i] +
                                 read_L.cmunu_configs_dic[i])
        Cmunu_Configs = np.array(Cmunu_Configs)

        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs * self.CmunuFactor, self.L[3])
        self.outputfile = fileDir
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlices_AL_dic(self, Aname, Lname, LSname, configList, fileDir):
        print 'Reading TimeSlice From Separate AMA LMA LMASUB files'

        Cmunu_Configs = []

        read_A = ReadCorrelationFunction(self.path, Aname, self.L, 'TimeSlice')
        read_A.read_allconfigs_samemu('VEC-CORRt')

        read_L = ReadCorrelationFunction(self.path, Lname, self.L, 'TimeSlice')
        read_L.read_allconfigs_samemu('VEC-CORRt')
        for i in read_L.cmunu_configs_dic:
            read_L.cmunu_configs_dic[i] = read_L.cmunu_configs_dic[i] / self.Vol

        read_LS = ReadCorrelationFunction(self.path, LSname, self.L, 'TimeSlice')
        read_LS.read_allconfigs_samemu('VEC-CORRt')
        for i in read_L.cmunu_configs_dic:
            read_LS.cmunu_configs_dic[i] = read_LS.cmunu_configs_dic[i] / self.numPointSrc

        # print 'Configs Read:'
        # print list(configList)

        for i in configList:
            Cmunu_Configs.append(read_A.cmunu_configs_dic[i] - read_LS.cmunu_configs_dic[i] +
                                 read_L.cmunu_configs_dic[i])
        Cmunu_Configs = np.array(Cmunu_Configs)

        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs * self.CmunuFactor, self.L[3])
        self.outputfile = fileDir
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlices_AS_dic(self, Aname, LSname, configList, fileDir):
        print 'Reading TimeSlice From Separate AMA LMASUB files'

        Cmunu_Configs = []

        read_A = ReadCorrelationFunction(self.path, Aname, self.L, 'TimeSlice')
        read_A.read_allconfigs_samemu('VEC-CORRt')

        read_LS = ReadCorrelationFunction(self.path, LSname, self.L, 'TimeSlice')
        read_LS.read_allconfigs_samemu('VEC-CORRt')
        for i in read_LS.cmunu_configs_dic:
            read_LS.cmunu_configs_dic[i] = read_LS.cmunu_configs_dic[i] / self.numPointSrc

        print 'Configs Read:'
        print list(configList)

        for i in configList:
            Cmunu_Configs.append(read_A.cmunu_configs_dic[i] - read_LS.cmunu_configs_dic[i])
        Cmunu_Configs = np.array(Cmunu_Configs)

        # print 'L'
        # print (read_L.Cmunu_Configs ).tolist()

        # print 'LS'
        # print (read_LS.Cmunu_Configs ).tolist()

        # print 'L-LS'
        # print (read_L.Cmunu_Configs - read_LS.Cmunu_Configs).tolist()


        # plt.errorbar(range(0, self.L[3]), (read_Exact.Cmunu_Configs - read_Sub.Cmunu_Configs + read_AMA.Cmunu_Configs)[0][0:], fmt='.b', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), (read_AMA.Cmunu_Configs)[0][1:], fmt='.r', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), read_L.Cmunu_Configs[0][1:], fmt=Lplt, elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), read_LS.Cmunu_Configs[0][1:], fmt=LSplt, elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), Cmunu_Configs[0][1:], fmt='.r', elinewidth=1)

        # plt.errorbar(range(1, self.L[3]), read_L.Cmunu_Configs[0][1:], fmt='.b', elinewidth=1)
        # plt.show()


        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs * self.CmunuFactor, self.L[3])
        self.outputfile = fileDir
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlices_A_dic(self, Aname, configList, fileDir):
        print 'Reading TimeSlice From Separate AMA files'

        Cmunu_Configs = []

        read_A = ReadCorrelationFunction(self.path, Aname, self.L, 'TimeSlice')
        read_A.read_allconfigs_samemu('VEC-CORRt')

        # print 'Configs Read:'
        # print list(configList)

        for i in configList:
            Cmunu_Configs.append(read_A.cmunu_configs_dic[i])
        Cmunu_Configs = np.array(Cmunu_Configs)

        # print 'L'
        # print (read_L.Cmunu_Configs ).tolist()

        # print 'LS'
        # print (read_LS.Cmunu_Configs ).tolist()

        # print 'L-LS'
        # print (read_L.Cmunu_Configs - read_LS.Cmunu_Configs).tolist()


        # plt.errorbar(range(0, self.L[3]), (read_Exact.Cmunu_Configs - read_Sub.Cmunu_Configs + read_AMA.Cmunu_Configs)[0][0:], fmt='.b', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), (read_AMA.Cmunu_Configs)[0][1:], fmt='.r', elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), read_L.Cmunu_Configs[0][1:], fmt=Lplt, elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), read_LS.Cmunu_Configs[0][1:], fmt=LSplt, elinewidth=1)
        # plt.errorbar(range(1, self.L[3]), Cmunu_Configs[0][1:], fmt='.r', elinewidth=1)

        # plt.errorbar(range(1, self.L[3]), read_L.Cmunu_Configs[0][1:], fmt='.b', elinewidth=1)
        # plt.show()


        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs * self.CmunuFactor, self.L[3])
        self.outputfile = fileDir
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runJackknife_L(self, Lname, configList, outDir):
        print '==============================='
        print 'Run Jackknife for LMA'
        print 'Save to ' + str(outDir)

        lenConfig = len(configList)
        if lenConfig == []:
            print 'No Config List'
            return
        if lenConfig == 1:
            print 'Only One Configuration'
            return
        for i in range(lenConfig):
            print 'In Jackknife Config ' + str(configList[i])
            JackConfig = configList[:i] + configList[i + 1:]
            # print JackConfig
            fileDir = outDir + '/' + 'LMA.jackknife-' + str(lenConfig) + '.' + str(configList[i])
            self.make_pi_qsqu_covmatx_LMA(JackConfig, fileDir)
            print 'Creat file' + fileDir
            print
        return

    def runJackknife_LS(self, LSname, configList, outDir):
        print '==============================='
        print 'Run Jackknife for LMASUB'
        print 'Save to ' + str(outDir)

        lenConfig = len(configList)
        if lenConfig == []:
            print 'No Config List'
            return
        if lenConfig == 1:
            print 'Only One Configuration'
            return
        for i in range(lenConfig):
            print 'In Jackknife Config ' + str(configList[i])
            JackConfig = configList[:i] + configList[i + 1:]
            # print JackConfig
            fileDir = outDir + '/' + 'LMASUB.jackknife-' + str(lenConfig) + '.' + str(configList[i])
            self.runTimeSlices_LS_dic(LSname, JackConfig, fileDir)
            print 'Creat file' + fileDir
            print
        return

    def runJackknife_A(self, Aname, configList, outDir):
        print '==============================='
        print 'Run Jackknife for AMA'
        print 'Save to ' + str(outDir)

        lenConfig = len(configList)
        if lenConfig == []:
            print 'No Config List'
            return
        if lenConfig == 1:
            print 'Only One Configuration'
            return
        for i in range(lenConfig):
            print 'In Jackknife Config ' + str(configList[i])
            JackConfig = configList[:i] + configList[i + 1:]
            # print JackConfig
            fileDir = outDir + '/' + 'AMA.jackknife-' + str(lenConfig) + '.' + str(configList[i])
            self.runTimeSlices_A_dic(Aname, JackConfig, fileDir)
            print 'Creat file' + fileDir
            print
        return

    def runJackknife_AS(self, Aname, LSname, configList, outDir):
        print '==============================='
        print 'Run Jackknife for AMA and LMASUB'
        print 'Save to ' + str(outDir)

        lenConfig = len(configList)
        if lenConfig == []:
            print 'No Config List'
            return
        if lenConfig == 1:
            print 'Only One Configuration'
            return
        for i in range(lenConfig):
            print 'In Jackknife Config ' + str(configList[i])
            JackConfig = configList[:i] + configList[i + 1:]
            # print JackConfig
            fileDir = outDir + '/' + 'AS.jackknife-' + str(lenConfig) + '.' + str(configList[i])
            self.runTimeSlices_AS_dic(Aname, LSname, JackConfig, fileDir)
            print 'Creat file' + fileDir
            print
        return

    # clear
    def runJackknife_EAL(self, Ename, ESname, Aname, Lname, LSname, configList, outDir):
        lenConfig = len(configList)
        if lenConfig == []:
            print 'No Config List'
            return
        if lenConfig == 1:
            print 'Only One Configuration'
            return

        if (self.read_E is None or
                    self.read_ES is None or
                    self.read_A is None or
                    self.read_L is None or
                    self.read_LS is None):
            print 'Need to run readAllConfig_EAL() before'
            return

        print '======================================================='
        print 'Run Jackknife for Exact, Sub, AMA, LMA and LMASUB'
        print 'Save to ' + str(outDir)

        for i in range(lenConfig):
            print 'In Jackknife Config ' + str(configList[i])
            JackConfig = configList[:i] + configList[i + 1:]
            fileDir = outDir + '/' + 'EAL.jackknife-' + str(lenConfig) + '.' + str(configList[i])
            self.make_pi_qsqu_covmatx_Exact_Sub_AMA_LMASUB_LMA(Ename, ESname, Aname, Lname, LSname, JackConfig, fileDir)
            print 'Creat file' + fileDir
            print
        return

    def runJackknife_AL(self, Aname, Lname, LSname, configList, outDir):
        print '==============================='
        print 'Run Jackknife for AMA, LMA and LMASUB'
        print 'Save to ' + str(outDir)

        lenConfig = len(configList)
        if lenConfig == []:
            print 'No Config List'
            return
        if lenConfig == 1:
            print 'Only One Configuration'
            return
        for i in range(lenConfig):
            print 'In Jackknife Config ' + str(configList[i])
            JackConfig = configList[:i] + configList[i + 1:]
            # print JackConfig
            fileDir = outDir + '/' + 'AL.jackknife-' + str(lenConfig) + '.' + str(configList[i])
            self.runTimeSlices_AL_dic(Aname, Lname, LSname, JackConfig, fileDir)
            print 'Creat file' + fileDir
            print
        return

    def runJackknife_EL(self, Ename, ESname, Lname, configList, outDir):
        print '==============================='
        print 'Run Jackknife for Exact, Sub and LMA'
        print 'Save to ' + str(outDir)

        lenConfig = len(configList)
        if lenConfig == []:
            print 'No Config List'
            return
        if lenConfig == 1:
            print 'Only One Configuration'
            return
        for i in range(lenConfig):
            print 'In Jackknife Config ' + str(configList[i])
            JackConfig = configList[:i] + configList[i + 1:]
            # print JackConfig
            fileDir = outDir + '/' + 'EL.jackknife-' + str(lenConfig) + '.' + str(configList[i])
            self.runTimeSlices_EL_dic(Ename, ESname, Lname, JackConfig, fileDir)
            print 'Creat file' + fileDir
            print
        return

    # clear
    def runJackknife_EA(self, Ename, ESname, Aname, configList, outDir, outName):
        lenConfig = len(configList)
        if lenConfig == []:
            print 'No Config List'
            return
        if lenConfig == 1:
            print 'Only One Configuration'
            return

        if (self.read_E is None or
                    self.read_ES is None or
                    self.read_A is None):
            print 'Need to run readAllConfig_EA() before'
            return

        if os.path.isdir(outDir) == False:
            print('mkdir ' + outDir)
            os.mkdir(outDir, 0755)

        print '==============================='
        print 'Run Jackknife for Exact, Sub and AMA'
        print 'Save to ' + str(outDir)

        for i in range(lenConfig):
            print 'In Jackknife Config ' + str(configList[i])
            JackConfig = configList[:i] + configList[i + 1:]
            fileDir = outDir + '/' + outName + str(lenConfig) + 'configs' + '.' + str(configList[i])
            self.make_pi_qsqu_covmatx_Exact_Sub_AMA(Ename, ESname, Aname, JackConfig, fileDir)
            print 'Creat file' + fileDir
            print
        return

    # clear
    def runJackknife_EAL_selectedT(self, Ename, ESname, Aname, Lname, LSname, configList, numT, outDir):
        lenConfig = len(configList)
        if lenConfig == []:
            print 'No Config List'
            return
        if lenConfig == 1:
            print 'Only One Configuration'
            return

        if (self.read_E is None or
                    self.read_ES is None or
                    self.read_A is None or
                    self.read_L is None or
                    self.read_LS is None):
            print 'Need to run readAllConfig_EAL() before'
            return

        print '======================================================='
        print 'Run Jackknife for Exact, Sub, AMA, LMA and LMASUB'
        print 'Set Cmunu(T) to 0, which T is the middle ' + str(numT)
        print 'Save to ' + str(outDir)

        for i in range(lenConfig):
            print 'In Jackknife Config ' + str(configList[i])
            JackConfig = configList[:i] + configList[i + 1:]
            fileDir = outDir + '/' + 'EAL.jackknife-' + str(lenConfig) + '-numMiddleT' + str(numT) + '.' + str(
                configList[i])
            self.runTimeSlices_EAL_dic_selectedT(JackConfig, numT, fileDir)
            print 'Creat file' + fileDir
            print
        return


class do_96:
    def __init__(self):
        self.path = '/Volumes/Seagate Backup Plus Drive/lqcdproj/gMinus2/blum/HISQ'
        # self.path = '/Users/tucheng/Desktop'
        self.L = [96, 96, 96, 192]
        self.outputfile = '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_t'

    def runTimeSlicesHalf(self):
        print 'Reading TimeSlice'
        read_t = ReadCorrelationFunction(self.path, 'vec_ama_xyzt', self.L, 'TimeSlice')
        read_t.read_allconfig_tslices(952320)
        qsqu = MakeQsqu(1, 1, 1, 30, self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, read_t.cmunu_configs, self.L[3] / 2)
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix,
                     '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_thalf').output()

    def runTimeSlices(self):
        print 'Reading TimeSlice'
        read_t = ReadCorrelationFunction(self.path, 'vec_ama_xyzt', self.L, 'TimeSlice')
        read_t.read_allconfig_tslices(952320)
        print read_t.cmunu_configs[0]
        qsqu = MakeQsqu(1, 1, 1, 30, self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, read_t.cmunu_configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_t').output()

    def runXYZSlices(self):
        print 'Reading SpacialSlice'
        read_xyz = ReadCorrelationFunction(self.path, 'vec_ama_xyzt', self.L, 'SpacialSlice')
        read_xyz.read_allconfig_xyzslices(952320)
        qsqu = MakeQsqu(30, 1, 1, 1, self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, read_xyz.cmunu_configs, self.L[0])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix,
                     '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_xyz').output()

    def runXSlices(self):
        print 'Reading X SpacialSlice'
        read_x = ReadCorrelationFunction(self.path, 'vec_ama_xyzt', self.L, 'SpacialSlice')
        read_x.read_allconfig_xslices(952320)
        qsqu = MakeQsqu(30, 1, 1, 1, self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, read_x.cmunu_configs, self.L[0])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_x').output()

    def runYSlices(self):
        print 'Reading Y SpacialSlice'
        read_y = ReadCorrelationFunction(self.path, 'vec_ama_xyzt', self.L, 'SpacialSlice')
        read_y.read_allconfig_yslices(952320)
        qsqu = MakeQsqu(30, 1, 1, 1, self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, read_y.cmunu_configs, self.L[0])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_y').output()

    def runZSlices(self):
        print 'Reading Y SpacialSlice'
        read_z = ReadCorrelationFunction(self.path, 'vec_ama_xyzt', self.L, 'SpacialSlice')
        read_z.read_allconfig_zslices(952320)
        qsqu = MakeQsqu(30, 1, 1, 1, self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, read_z.cmunu_configs, self.L[0])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/96c fitting/c_96_z').output()

    def runTimeSlicesFromSeperateFiles(self):
        print 'Reading TimeSlice'

        read_Exact = ReadCorrelationFunction(self.path, 'vec_Exact', self.L, 'TimeSlice')
        read_Exact.read_allconfigs_samemu('VEC-CORRt')
        read_Sub = ReadCorrelationFunction(self.path, 'vec_Sub', self.L, 'TimeSlice')
        read_Sub.read_allconfigs_samemu('VEC-CORRt')
        read_AMA = ReadCorrelationFunction(self.path, 'vec_AMA', self.L, 'TimeSlice')
        read_AMA.read_allconfigs_samemu('VEC-CORRt')
        Cmunu_Configs = read_Exact.cmunu_configs - read_Sub.cmunu_configs + read_AMA.cmunu_configs
        # print Cmunu_Configs[0]
        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        # print qsqu.qsqu
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()


class do_4444:
    def __init__(self):
        self.path = '/Users/tucheng/Desktop/HISQ/test44'
        self.L = [4, 4, 4, 4]
        self.outputfile = '/Users/tucheng/Desktop/Fitting/results/4444c/test44'

    def runTimeSlices(self):
        print 'Reading TimeSlice'
        read_t = ReadCorrelationFunction(self.path, 'vec', self.L, 'TimeSlice')
        read_t.read_allconfig_AllInOne()
        qsqu = MakeQsqu(1, 1, 1, 30, self.L)
        pi_matrix = MakePiMatrix(qsqu.qsqu, read_t.cmunu_configs, self.L[3])
        # Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, '/Users/tucheng/Desktop/Fitting/results/64c fitting/c_64_t').output()

    def runTimeSlicesFromSeparateFiles(self):
        print 'Reading TimeSlice From Separate Exact Aub AMA files'

        read_Exact = ReadCorrelationFunction(self.path, 'vec_Exact', self.L, 'TimeSlice')
        read_Exact.read_allconfigs_samemu('VEC-CORRt')

        read_Sub = ReadCorrelationFunction(self.path, 'vec_Sub', self.L, 'TimeSlice')
        read_Sub.read_allconfigs_samemu('VEC-CORRt')

        read_AMA = ReadCorrelationFunction(self.path, 'vec_AMA', self.L, 'TimeSlice')
        read_AMA.read_allconfigs_samemu('VEC-CORRt')

        Cmunu_Configs = read_Exact.cmunu_configs - read_Sub.cmunu_configs + read_AMA.cmunu_configs
        # print Cmunu_Configs[0]
        qsqu = MakeQsqu(1, 1, 1, self.L[3], self.L)
        # print qsqu.qsqu
        pi_matrix = MakePiMatrix(qsqu.qsqu, Cmunu_Configs, self.L[3])
        MakeCovarianceMatrix(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlicesFromSingleFile_LMAcode(self):
        print 'Reading TimeSlice From LMAcode'

        read_LMA = ReadCorrelationFunction(self.path, 'test_2221_LMA_src1212', self.L, 'TimeSlice')
        read_LMA.read_allconfigs_samemu('VEC-CORRt')
        read_LMA.cmunu_configs = read_LMA.cmunu_configs / (self.L[0] * self.L[1] * self.L[2] * self.L[3])

        print 'LMA All Srcs'
        print read_LMA.cmunu_configs

        read_LMASUB = ReadCorrelationFunction(self.path, 'test_2221_LMASUB_src1212', self.L, 'TimeSlice')
        read_LMASUB.read_allconfigs_samemu('VEC-CORRt')
        read_LMASUB.cmunu_configs = read_LMASUB.cmunu_configs / (
        (self.L[0] - self.L[0] / 2.0) ** 3.0 * (self.L[3] - self.L[3] / 2.0))

        print 'LMA Sub Srcs'
        print read_LMASUB.cmunu_configs



        # qsqu = Make_qsqu(1, 1, 1, self.L[3], self.L)
        # pi_matrix = Make_Pi_matrix(qsqu.qsqu, read_LMA.Cmunu_Configs, self.L[3])
        # Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()

    def runTimeSlicesFromSingleFile_LOWMODEAPPROX(self):
        print 'Reading TimeSlice From LOWMODEAPPROX'

        read_One = ReadCorrelationFunction(self.path, 'LMA_eigen10_LOWMODEAPPROX_nocontact_2221', self.L, 'TimeSlice')
        read_One.read_allconfigs_samemu('VEC-CORRt')

        print 'LOWAPPROX All Srcs'
        print read_One.cmunu_configs

        read_One = ReadCorrelationFunction(self.path, 'test_2221_LMAAPPROX_src1212', self.L, 'TimeSlice')
        read_One.read_allconfigs_samemu('VEC-CORRt')

        print 'LOWAPPROX Sub Srcs'
        print read_One.cmunu_configs

        # qsqu = Make_qsqu(1, 1, 1, self.L[3], self.L)
        # pi_matrix = Make_Pi_matrix(qsqu.qsqu, read_One.Cmunu_Configs, self.L[3])
        # Make_CovMatx(qsqu.qsqu, pi_matrix.matrix, self.outputfile).output()


start = do_6496('/Users/tucheng/Desktop/fitting/data/test3000/')
start.read_allconfig_LMA('vec_LMA')
start.make_pi_qsqu_covmatx_LMA([720], '/Users/tucheng/Desktop/fitting/data/test3000/out3000')

#start1 = do_6496('/Users/tucheng/Desktop/fitting/data/test-fnal/')
#start1.read_allconfig_LMA('vec_LMA')
#start1.make_pi_qsqu_covmatx_LMA([504], '/Users/tucheng/Desktop/fitting/data/test-fnal/out-fnal')
plt.show()
# start.readAllConfig_EA('vec_Exact-strange', 'vec_Sub-strange', 'vec_AMA-strange')


# start.runTimeSlicesFromSeparateFiles_EAL()
# start.runTimeSlicesFromSeparateFiles_ESA()
# start.runTimeSlicesFromSeparateFiles_L_LS('vec_sLMA-1000evecs', 'vec_sLS-1000evecs', '.k', '.g')
# start.runTimeSlicesFromSeparateFiles_L_LS('vec_LMA-single', 'vec_LS-double_new', '.g', '.b')
# start.runTimeSlicesFromSeparateFiles_L_LS('vec_dLMA-1000evecs', 'vec_dLS-1000evecs', '.b', '.r')
# start.runTimeSlicesFromSeparateFiles_L_LS_LOWMODEAPPROX('vec_LMA-single', 'vec_LMAPPRX-double', '.g', '.r')
# start.runTimeSlicesFromSeparateFiles_L_LS_LOWMODEAPPROX('vec_LMA-single', 'vec_LMAPPRX-single', '.g', '.b')
# plt.yscale('log')
# start.runTimeSlices_EAL_dic('vec_Exact', 'vec_Sub', 'vec_AMA', 'vec_LMA', 'vec_LMASUB')
# start.runTimeSlices_L_dic('vec_LMA', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128])
# start.runTimeSlices_EL_dic('vec_Exact', 'vec_Sub','vec_LMA', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128])
# start.runTimeSlices_AL_dic('vec_AMA', 'vec_LMA', 'vec_LMASUB', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128])
# start.runTimeSlices_AS_dic('vec_AMA', 'vec_LMASUB', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128])
# start.runTimeSlices_LS_dic('vec_LMASUB', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128], '/Users/tucheng/Desktop/fitting/cmatrix/l6496/21configs_LS')
# start.runTimeSlices_A_dic('vec_AMA', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128])
# start.runTimeSlices_LS_dic('vec_LMASUB', [1092])
# start.runTimeSlices_L_dic('vec_LMA', [504], '/Users/tucheng/Desktop/fitting/data/test/out')

# plt.show()
# start.runTimeSlicesFromSeparateFiles_LS()
# start.runTimeSlices_LS_dic('vec_LMASUB', [1020], '/Users/tucheng/Desktop/fitting/cmatrix/l6496/test')

# 1092 1104 1164 372 1080 1140 1176 1020 540 1128

# start.runJackknife_L('vec_LMA', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128], '/Users/tucheng/Desktop/fitting/cmatrix/l6496/jackknife/')
# start.runJackknife_LS('vec_LMASUB', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128], '/Users/tucheng/Desktop/fitting/cmatrix/l6496/jackknife/')
# start.runJackknife_A('vec_AMA', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128], '/Users/tucheng/Desktop/fitting/cmatrix/l6496/jackknife/')
# start.runJackknife_AS('vec_AMA', 'vec_LMASUB', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128], '/Users/tucheng/Desktop/fitting/cmatrix/l6496/jackknife/')
# start.runJackknife_EAL('vec_Exact', 'vec_Sub', 'vec_AMA', 'vec_LMA', 'vec_LMASUB', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128], '/Users/tucheng/Desktop/fitting/cmatrix/l6496/jackknife/')
# start.runJackknife_EA('vec_Exact-strange', 'vec_Sub-strange', 'vec_AMA-strange',
# start.conf_num_list,
# '/Users/tucheng/Desktop/fitting/cmatrix/l6496/jackknife-strange/', 'EA.jackknife-strange.')

# start.runJackknife_EA('vec_Exact', 'vec_Sub', 'vec_AMA', [1152, 672, 204, 1092, 504, 1056, 1032, 1068, 1104, 684, 1008, 456, 1164, 372, 528, 1080, 1140, 1176, 1020, 540, 1128], '/Users/tucheng/Desktop/fitting/cmatrix/l6496/jackknife/')
# start.runTimeSlices_L_dic('vec_LMA', [204], '/Users/tucheng/Desktop/fitting/cmatrix/l6496/test')
# do_64()


# start = do_96()
# start.runTimeSlicesHalf()

# start = do_6496()
# start.runTimeSlices()

# start = do_48()
# start.runTimeSlices()

# start = do_4444()
# start.runTimeSlicesFromSingleFile_LMAcode()
# start.runTimeSlicesFromSingleFile_LOWMODEAPPROX()

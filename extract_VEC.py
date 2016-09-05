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
import shutil
import time

class mix_exact_sub_sloppy:

    '''
    Mix Exact, Sub, AMA terms into one file
    This is designed for 96192 lattice
    '''

    def __init__(self, path):

        '''
        Read in path and build a folder list under that path
        :param path: HISQ path
        '''

        self.path = path
        self.folderlist = (src.func.walkfiles(self.path, prt=0))[0]

    def share_exact(self):

        '''
        Because x0-t0 and x16-t24 folders share the same Exact_Sub values and these files just locate under x0-t0,
        so this function is to copy Exact_Sub to x16-t24
        The size of Exact_Sub should be > 10^6
        :return:
        '''

        for folder_x0t0 in self.folderlist:
            if ('l96' in folder_x0t0) & ('x0-t0' in folder_x0t0):
                exact_num = 0
                conf_num = folder_x0t0.split(".")[-1]
                for folder_x16t24 in self.folderlist:
                    if ('l96' in folder_x16t24) & ('x16-t24' in folder_x16t24) & (('.' + conf_num) in folder_x16t24):

                        filelist = (src.func.walkfiles(self.path + '/' + folder_x0t0, prt=0))[1]
                        for file in filelist:
                            if ('exact' in file) & (('.' + conf_num + '.') in file):
                                exact_size = os.path.getsize(self.path + '/' + folder_x0t0 + '/' + file)
                                if exact_size > 10.0 ** 6.0:
                                    exact_num = 1
                                    file_exact = file
                        if exact_num == 1:
                            print 'copy ' + self.path + '/' + folder_x0t0 + '/' + file_exact + \
                                  ' to ' + self.path + '/' + folder_x16t24 + '/'

                            shutil.copy(self.path + '/' + folder_x0t0 + '/' + file_exact,
                                        self.path + '/' + folder_x16t24 + '/')

    def Combine_96(self):

        '''
        Combine Exact_Sub and AMA files into 'vec_ama_xyzt' files.
        First check sizes of Exact_Sub (> 10^6) and times and sizes of AMA (> 10^7 and newest)
        and read in 'src_global_xyzt' lines because the values of them in below lines are wrong
        Then read and write the 'VEC' lines
        :return: NA
        '''

        wrong_config = []
        right_config_num = 0
        for folder in self.folderlist:
            if ('l96' in folder) & ('t' in folder):
                conf_num = folder.split(".")[-1]
                filelist = (src.func.walkfiles(self.path + '/' + folder, prt=0))[1]


                exact_num, sloppy_num = 0, 0
                sloppy_time_old = 0
                for file in filelist:
                    if ('exact' in file) & (('.' + conf_num + '.') in file):
                        exact_size = os.path.getsize(self.path + '/' + folder + '/' + file)
                        if exact_size > 10.0 ** 6.0:
                            exact_num = 1
                            file_exact = file

                    if ('sloppy' in file) & (('.' + conf_num + '.') in file):
                        sloppy_time = os.path.getmtime(self.path + '/' + folder + '/' + file)
                        sloppy_size = os.path.getsize(self.path + '/' + folder + '/' + file)
                        if (sloppy_size > 10.0 ** 7.0) & (sloppy_time > sloppy_time_old):
                            sloppy_num = 1
                            file_sloppy = file
                            sloppy_time_old = sloppy_time

                #print self.path + '/' + folder + '/'
                #print exact_num, sloppy_num

                if (exact_num == 1) & (sloppy_num == 1):
                    right_config_num += 1

                    file = open(self.path + '/' + folder + '/' + 'vec_ama_xyzt' + '.' + conf_num, 'w')

                    data_exact_read = open(self.path + '/' + folder + '/' + file_exact)
                    data_sloppy_read = open(self.path + '/' + folder + '/' + file_sloppy)

                    num_src_global_xyzt_exact = 0
                    num_src_global_xyzt_sloppy = 0

                    for lines in data_exact_read.readlines():
                        if 'src_global_xyzt' in lines:
                            src_global_xyzt = lines.split()
                            s_x = src_global_xyzt[1]
                            s_y = src_global_xyzt[2]
                            s_z = src_global_xyzt[3]
                            s_t = src_global_xyzt[4]
                            num_src_global_xyzt_exact += 1

                        if 'VEC' in lines:
                            linestr = lines.split()
                            if 'VEC-CORRt' in lines:
                                file.write(linestr[0] + ' ' + linestr[1] + ' ' + s_t + ' ' + linestr[3] + ' '
                                           + linestr[4] + ' ' + linestr[5] + ' ' + linestr[6] + ' ' + linestr[7]
                                           + ' ' + linestr[8] + ' ' + linestr[9] + ' ' + linestr[10] + '\n')
                            if 'VEC-CORRx' in lines:
                                file.write(linestr[0] + ' ' + linestr[1] + ' ' + s_x + ' ' + linestr[3] + ' '
                                           + linestr[4] + ' ' + linestr[5] + ' ' + linestr[6] + ' ' + linestr[7]
                                           + ' ' + linestr[8] + ' ' + linestr[9] + ' ' + linestr[10] + '\n')
                            if 'VEC-CORRy' in lines:
                                file.write(linestr[0] + ' ' + linestr[1] + ' ' + s_y + ' ' + linestr[3] + ' '
                                           + linestr[4] + ' ' + linestr[5] + ' ' + linestr[6] + ' ' + linestr[7]
                                           + ' ' + linestr[8] + ' ' + linestr[9] + ' ' + linestr[10] + '\n')
                            if 'VEC-CORRz' in lines:
                                file.write(linestr[0] + ' ' + linestr[1] + ' ' + s_z + ' ' + linestr[3] + ' '
                                           + linestr[4] + ' ' + linestr[5] + ' ' + linestr[6] + ' ' + linestr[7]
                                           + ' ' + linestr[8] + ' ' + linestr[9] + ' ' + linestr[10] + '\n')

                    for lines in data_sloppy_read.readlines():
                        if 'src_global_xyzt' in lines:
                            src_global_xyzt = lines.split()
                            s_x = src_global_xyzt[1]
                            s_y = src_global_xyzt[2]
                            s_z = src_global_xyzt[3]
                            s_t = src_global_xyzt[4]
                            num_src_global_xyzt_sloppy += 1
                        if 'VEC' in lines:
                            linestr = lines.split()
                            if 'VEC-CORRt' in lines:
                                file.write(linestr[0] + ' ' + linestr[1] + ' ' + s_t + ' ' + linestr[3] + ' '
                                           + linestr[4] + ' ' + linestr[5] + ' ' + linestr[6] + ' ' + linestr[7]
                                           + ' ' + linestr[8] + ' ' + linestr[9] + ' ' + linestr[10] + '\n')
                            if 'VEC-CORRx' in lines:
                                file.write(linestr[0] + ' ' + linestr[1] + ' ' + s_x + ' ' + linestr[3] + ' '
                                           + linestr[4] + ' ' + linestr[5] + ' ' + linestr[6] + ' ' + linestr[7]
                                           + ' ' + linestr[8] + ' ' + linestr[9] + ' ' + linestr[10] + '\n')
                            if 'VEC-CORRy' in lines:
                                file.write(linestr[0] + ' ' + linestr[1] + ' ' + s_y + ' ' + linestr[3] + ' '
                                           + linestr[4] + ' ' + linestr[5] + ' ' + linestr[6] + ' ' + linestr[7]
                                           + ' ' + linestr[8] + ' ' + linestr[9] + ' ' + linestr[10] + '\n')
                            if 'VEC-CORRz' in lines:
                                file.write(linestr[0] + ' ' + linestr[1] + ' ' + s_z + ' ' + linestr[3] + ' '
                                           + linestr[4] + ' ' + linestr[5] + ' ' + linestr[6] + ' ' + linestr[7]
                                           + ' ' + linestr[8] + ' ' + linestr[9] + ' ' + linestr[10] + '\n')

                    file.close()
                    print 'create ' + self.path + '/' + folder + '/' + 'vec_ama_xyzt' + '.' + conf_num
                    print 'num_src_global_xyzt_exact: ' + str(num_src_global_xyzt_exact)
                    print 'num_src_global_xyzt_sloppy: ' + str(num_src_global_xyzt_sloppy)
                    data_exact_read.close()
                    data_sloppy_read.close()

                else:
                    wrong_config.append(folder)

        print 'wrong_config:'
        for wrong in wrong_config:
            print wrong

class Extract_VEC:

    '''
    Extract 'VEC' lines from 'out-ama' files to 'vec' files
    !!Attention!!: those 'out-ama' files should contain Exact and Sub and AMA in order (usually 8 & 8 & rest)
    '''

    def __init__(self, path, ensemble):

        '''
        Read in path and build a folder list under that path
        And set the ensemble label ('l64', 'l48'...)
        :param path: HISQ path
        :param ensemble: 'l64', 'l48'...
        '''

        self.path = path
        self.folderlist = (src.func.walkfiles(self.path, prt=0))[0]
        self.ensemble = ensemble

    def run(self):

        '''
        First check the sizes and the times of 'out-ama' (> 10^7 and newest)
        Read and write 'VEC' lines
        !!Attention Again!!: those 'out-ama' files should contain Exact and Sub and AMA in order (usually 8 & 8 & rest)
        :return: NA
        '''

        wrong_config = []
        right_config_num = 0
        for folder in self.folderlist:
            if (self.ensemble in folder) & ('t' in folder):
                conf_num = folder.split(".")[-1]
                filelist = (src.func.walkfiles(self.path + '/' + folder, prt=0))[1]

                out_ama_num = 0
                ama_time_old = 0
                for file in filelist:

                    if ('out-ama' in file) & (('.' + conf_num + '.') in file):
                        out_ama_time = os.path.getmtime(self.path + '/' + folder + '/' + file)
                        out_ama_size = os.path.getsize(self.path + '/' + folder + '/' + file)
                        if (out_ama_size > 10.0 ** 7.0) & (out_ama_time > ama_time_old):
                            out_ama_num = 1
                            file_out_ama = file
                            ama_time_old = out_ama_time

                if (out_ama_num == 1):
                    right_config_num += 1

                    file = open(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num, 'w')

                    data_out_ama_read = open(self.path + '/' + folder + '/' + file_out_ama)

                    num_src_global_xyzt_out_ama = 0

                    for lines in data_out_ama_read.readlines():
                        if 'src_global_xyzt' in lines:
                            num_src_global_xyzt_out_ama += 1
                        if 'VEC' in lines:
                            file.write(lines)

                    file.close()
                    print 'create ' + self.path + '/' + folder + '/' + 'vec' + '.' + conf_num
                    print 'num_src_global_xyzt_out_ama: ' + str(num_src_global_xyzt_out_ama)
                    data_out_ama_read.close()

                else:
                    wrong_config.append(folder)

        print 'wrong_config:'
        for wrong in wrong_config:
            print wrong

    def check_Size_Date(self, folder, size):

        '''
        Check Size and Compare Date of the files under the folder
        :param folder:
        :param size:
        :return: True/False: have/no such file, and give file name if True
        '''

        conf_num = folder.split(".")[-1]
        filelist = (src.func.walkfiles(self.path + '/' + folder, prt=0))[1]
        ama_time_old = 0
        out_ama_num = 0
        for file in filelist:
            if ('out-ama' in file) & (('.' + conf_num + '.') in file):
                out_ama_time = os.path.getmtime(self.path + '/' + folder + '/' + file)
                out_ama_size = os.path.getsize(self.path + '/' + folder + '/' + file)
                if (out_ama_size > size) & (out_ama_time > ama_time_old):
                    out_ama_num = 1
                    file_out_ama = file
                    ama_time_old = out_ama_time
        if (out_ama_num == 1):
            return True, file_out_ama
        else:
            return False, False

    def check_Num_Size(self, folder, name):

        '''
        Check the Number of files named by 'name' under the folder
        :param folder:
        :param name:
        :return:
        '''

        num1 = 0
        num2 = 0
        filename = []
        filelist = (src.func.walkfiles(self.path + '/' + folder, prt=0))[1]
        for file in filelist:
            if name in file:
                out_ama_size = os.path.getsize(self.path + '/' + folder + '/' + file)
                if ((out_ama_size > 10.0 ** 6.0) and (out_ama_size < 10.0 ** 7.0)):
                    num1 +=1
                    filename.append(file)

        for file in filelist:
            if name in file:
                out_ama_size = os.path.getsize(self.path + '/' + folder + '/' + file)
                if (out_ama_size > 10.0 ** 7.0):
                    num2 += 1
                    filename.append(file)
        return num1, num2, filename

    def check_Complete(self, FullFilePath):

        '''
        Check if the file contain all Exact Sub and AMA values
        BY SIMPLY CHECK THE FIRST TWO 'src_global_xyzt' LINES. ONLY t CHANGE IN FIRST TWO LINES WILL BE ACCEPTED
        :return:
        '''

        file = open(FullFilePath)
        linenum = 0
        s_x = 0
        s_y = 0
        s_z = 0
        s_t = 0
        for line in file:
            if 'src_global_xyzt' in line:
                s_x_last = s_x
                s_y_last = s_y
                s_z_last = s_z
                s_t_last = s_t
                linenum += 1
                src_global_xyzt = line.split()
                s_x = src_global_xyzt[1]
                s_y = src_global_xyzt[2]
                s_z = src_global_xyzt[3]
                s_t = src_global_xyzt[4]
                if linenum == 1:
                    s_x_last = s_x
                    s_y_last = s_y
                    s_z_last = s_z
                    s_t_last = s_t
                if linenum == 2:
                    if (s_x_last == s_x) & (s_y_last == s_y) & (s_z_last == s_z):
                        return True
                    else:
                        return False

    def run_check_Size_Date_Complete(self):

        '''
        First check Size_Date_Complete
        Read and write 'VEC' lines
        :return: NA
        '''

        right_config_num = 0
        NotCompleteList = []
        for folder in self.folderlist:
            if (self.ensemble in folder) & ('t' in folder):
                checkres1 = self.check_Size_Date(folder, 10.0 ** 7.0)

                if checkres1[0]:

                    checkres2 = self.check_Complete(self.path + '/' + folder + '/' + checkres1[1])

                    if checkres2:

                        conf_num = folder.split(".")[-1]
                        right_config_num += 1

                        file = open(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num, 'w')

                        data_out_ama_read = open(self.path + '/' + folder + '/' + checkres1[1])

                        num_src_global_xyzt_out_ama = 0

                        for lines in data_out_ama_read.readlines():
                            if 'src_global_xyzt' in lines:
                                num_src_global_xyzt_out_ama += 1
                            if 'VEC' in lines:
                                file.write(lines)

                        file.close()
                        print 'create ' + self.path + '/' + folder + '/' + 'vec' + '.' + conf_num
                        print 'num_src_global_xyzt_out_ama: ' + str(num_src_global_xyzt_out_ama)
                        data_out_ama_read.close()
                    else:
                        NotCompleteList.append(self.path + '/' + folder + '/' + checkres1[1])
                        conf_num = folder.split(".")[-1]
                        file = open(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num, 'w')
                        file.close()
        #print out the 'Not Complete' List
        print 'Not Complete:'
        for list in NotCompleteList:
            print list

    def run_Num_Size(self):

        '''
        First Check the number of the 'out-ama' files
        In l48, there are two situations (actually three, the last one is total no file):
        1. Exact_Sub and AMA are separated (Size of E_S: 10^6-10^7, Size of AMA: >10^7)
           I also need to check E_S is really E_S and AMA is really AMA by using check_Complete()
        2. All of them are in one file (Size: >10^7)
           The same, need to check_Complete()
        :return:
        '''

        NotCompleteList = []
        for folder in self.folderlist:
            if (self.ensemble in folder) & ('t' in folder):
                check_Num_Size = self.check_Num_Size(folder, 'out-ama')
                print folder
                print check_Num_Size

                # All in one file
                if (check_Num_Size[0] == 0) and (check_Num_Size[1] == 1):
                    if self.check_Complete(self.path + '/' + folder + '/' + check_Num_Size[2][0]):
                        num_src_global_xyzt_out_ama = 0
                        conf_num = folder.split(".")[-1]
                        fileread = open(self.path + '/' + folder + '/' + check_Num_Size[2][0])
                        filewrite = open(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num, 'w')
                        for lines in fileread.readlines():
                            if 'src_global_xyzt' in lines:
                                num_src_global_xyzt_out_ama += 1
                                filewrite.write(lines)
                            if 'VEC' in lines:
                                filewrite.write(lines)
                        fileread.close()
                        filewrite.close()
                        print 'create ' + self.path + '/' + folder + '/' + 'vec' + '.' + conf_num
                        print 'num_src_global_xyzt_out_ama: ' + str(num_src_global_xyzt_out_ama)
                    else:
                        NotCompleteList.append(self.path + '/' + folder + '/')
                        conf_num = folder.split(".")[-1]
                        file = open(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num, 'w')
                        file.close()

                # Exact_Sub and AMA are separated
                elif (check_Num_Size[0] == 1) and (check_Num_Size[1] == 1):
                    if self.check_Complete(self.path + '/' + folder + '/' + check_Num_Size[2][0]) and \
                            (not self.check_Complete(self.path + '/' + folder + '/' + check_Num_Size[2][1])):
                        conf_num = folder.split(".")[-1]
                        num_src_global_xyzt_out_ama = 0
                        fileread1 = open(self.path + '/' + folder + '/' + check_Num_Size[2][0])
                        filewrite = open(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num, 'w')
                        for lines in fileread1.readlines():
                            if 'src_global_xyzt' in lines:
                                num_src_global_xyzt_out_ama += 1
                                filewrite.write(lines)
                            if 'VEC' in lines:
                                filewrite.write(lines)
                        fileread1.close()
                        fileread2 = open(self.path + '/' + folder + '/' + check_Num_Size[2][1])
                        for lines in fileread2.readlines():
                            if 'src_global_xyzt' in lines:
                                num_src_global_xyzt_out_ama += 1
                                filewrite.write(lines)
                            if 'VEC' in lines:
                                filewrite.write(lines)
                        fileread2.close()
                        filewrite.close()
                        print 'create ' + self.path + '/' + folder + '/' + 'vec' + '.' + conf_num
                        print 'num_src_global_xyzt_out_ama: ' + str(num_src_global_xyzt_out_ama)
                    else:
                        NotCompleteList.append(self.path + '/' + folder + '/')
                        conf_num = folder.split(".")[-1]
                        file = open(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num, 'w')
                        file.close()
                else:
                    NotCompleteList.append(self.path + '/' + folder + '/')
                    conf_num = folder.split(".")[-1]
                    file = open(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num, 'w')
                    file.close()

        # print out the 'Not Complete' List
        print 'Not Complete:'
        for list in NotCompleteList:
            print list

    def Separate_E_S_A(self):

        '''
        Read in the files that contain all three Exact Sub AMA
        And Write them into three files individually
        :return:
        '''

        print 'Separate Exact Sub AMA:'
        for folder in self.folderlist:
            if (self.ensemble in folder) & ('t' in folder):
                conf_num = folder.split(".")[-1]
                if os.path.isfile(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num):
                    if os.path.getsize(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num) != 0:
                        vec_Exact = open(self.path + '/' + folder + '/' + 'vec_Exact' + '.' + conf_num, 'w')
                        vec_Sub = open(self.path + '/' + folder + '/' + 'vec_Sub' + '.' + conf_num, 'w')
                        vec_AMA = open(self.path + '/' + folder + '/' + 'vec_AMA' + '.' + conf_num, 'w')
                        FileRead = open(self.path + '/' + folder + '/' + 'vec' + '.' + conf_num, 'r')
                        AllLines = FileRead.readlines()
                        TotalLineNum = len(AllLines)
                        LineNum = -1
                        src_linenum = [] # The line number of the line with 'src_global_xyzt: 0 0 0 0'
                        for lines in AllLines:
                            LineNum += 1
                            if 'src_global_xyzt: 0 0 0 0' in lines:
                                src_linenum.append(LineNum)
                        if len(src_linenum) == 3:
                            for n in range(src_linenum[0], src_linenum[1]):
                                if 'VEC' in AllLines[n]:
                                    vec_Exact.write(AllLines[n])
                            vec_Exact.close()
                            print 'create ' + self.path + '/' + folder + '/' + 'vec_Exact' + '.' + conf_num
                            for n in range(src_linenum[1], src_linenum[2]):
                                if 'VEC' in AllLines[n]:
                                    vec_Sub.write(AllLines[n])
                            vec_Sub.close()
                            print 'create ' + self.path + '/' + folder + '/' + 'vec_Sub' + '.' + conf_num
                            for n in range(src_linenum[2], TotalLineNum):
                                if 'VEC' in AllLines[n]:
                                    vec_AMA.write(AllLines[n])
                            vec_AMA.close()
                            print 'create ' + self.path + '/' + folder + '/' + 'vec_AMA' + '.' + conf_num
                            FileRead.close()
                        else:
                            print 'wrong config:'
                            print self.path + '/' + folder + '/' + 'vec' + '.' + conf_num







'''
Mix_96 = mix_exact_sub_sloppy('/Volumes/Seagate Backup Plus Drive/lqcdproj/gMinus2/blum/HISQ/')
Mix_96.share_exact()
Mix_96.Combine_96()
'''

'''
l64 = Extract_VEC('/Volumes/Seagate Backup Plus Drive/lqcdproj/gMinus2/blum/HISQ/', 'l64')
l64.run()
'''

l48 = Extract_VEC('/Volumes/Seagate Backup Plus Drive/lqcdproj/gMinus2/blum/HISQ/', 'l48')
#l48.run_check_Size_Date_Complete()
l48.run_Num_Size()
l48.Separate_E_S_A()

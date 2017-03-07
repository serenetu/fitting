import numpy as np
from scipy.optimize import minimize
from scipy import integrate
from scipy.integrate import quadrature
from scipy.integrate import fixed_quad
import matplotlib.pyplot as plt
import random
import copy
import math
import sys
import os
import shutil


def walkfiles(path, prt = 0):
    dir_name = []
    file_names = []
    for dirpath, dirnames, filenames in os.walk(path):
        dir_name.append(dirnames)
        file_names.append(filenames)
    if dir_name == []:
        dir_name.append([])
    if file_names == []:
        file_names.append([])
    if prt == 1:
        print dir_name[0]
        print file_names[0]
    return dir_name[0], file_names[0]


def sort_config(filename, prt = 0):
    name_splt = []
    for i in filename:
        splt = i.split(".")
        splt[1] = int(splt[1])
        name_splt.append(splt)
    name_sort = sorted(name_splt, key=lambda name_splt : name_splt[1])
    for i in range (0, len(name_sort)):
        name_sort[i] = '.'.join(str(x) for x in name_sort[i])
    if prt == 1:
        print name_sort
    return name_sort


def read_confignm_AMA_Exact_SUB(L, path, prt = 0):
    hvp_nm, hvpAMA_nm, hvpExact_nm, hvpSUB_nm = [], [], [], []
    if int(L) == 48:
        for filenames in walkfiles(path, prt = 0)[1]:
            if 'l48144f211b672m0048m024m286a' in filenames:
                name_config = int(filenames.split(".")[-1])
                if 'hvpAMA' in filenames:
                    hvpAMA_nm.append(name_config)
                if 'hvpExact' in filenames:
                    hvpExact_nm.append(name_config)
                if 'hvpSUB' in filenames:
                    hvpSUB_nm.append(name_config)
        hvp_nm = sorted(list(set(hvpAMA_nm) & set(hvpExact_nm) & set(hvpSUB_nm)))
        if prt == 1:
            print 'configurations list:'
            print hvp_nm
            print 'number of configurations:'
            print len(hvp_nm)
            print
        return hvp_nm, sorted(hvpAMA_nm), sorted(hvpExact_nm), sorted(hvpSUB_nm)
    elif int(L) == 64:
        for filenames in walkfiles(path, prt = 0)[1]:
            if 'l64144f211b672m0024m024m286a' in filenames:
                name_config = int(filenames.split(".")[-1])
                if 'hvpAMA' in filenames:
                    hvpAMA_nm.append(name_config)
                if 'hvpExact' in filenames:
                    hvpExact_nm.append(name_config)
                if 'hvpSUB' in filenames:
                    hvpSUB_nm.append(name_config)
        hvp_nm = sorted(list(set(hvpAMA_nm) & set(hvpExact_nm) & set(hvpSUB_nm)))
        if prt == 1:
            print 'configurations list:'
            print hvp_nm
            print 'number of configurations:'
            print len(hvp_nm)
            print
        return hvp_nm, sorted(hvpAMA_nm), sorted(hvpExact_nm), sorted(hvpSUB_nm)
    elif int(L) == 96:
        for filenames in walkfiles(path, prt = 0)[1]:
            if 'l96192f211b672m0008m022m260a' in filenames:
                name_config = int(filenames.split(".")[-1])
                if 'hvpAMA' in filenames:
                    hvpAMA_nm.append(name_config)
                if 'hvpExact' in filenames:
                    hvpExact_nm.append(name_config)
                if 'hvpSUB' in filenames:
                    hvpSUB_nm.append(name_config)
        hvp_nm = sorted(list(set(hvpAMA_nm) & set(hvpExact_nm) & set(hvpSUB_nm)))
        if prt == 1:
            print 'configurations list:'
            print hvp_nm
            print 'number of configurations:'
            print len(hvp_nm)
            print 'configurations list of hvpAMA:'
            print sorted(hvpAMA_nm)
            print 'configurations list of hvpExact:'
            print sorted(hvpExact_nm)
            print 'configurations list of hvpSUB:'
            print sorted(hvpSUB_nm)
        return hvp_nm, sorted(hvpAMA_nm), sorted(hvpExact_nm), sorted(hvpSUB_nm)



def cp_rn(HISQ_PATH, L, confignm, HISQ_RN, prt = 0):
    hvpAMA_nm, hvpExact_nm, hvpSUB_nm = [], [], []
    if L == 48:
        for name in confignm:
            if 'l48144f211b672m0048m024m286a' in name:
                name_1 = name.split("-")[0]
                name_config = name.split(".")[-1]
                for root, dirs, files in os.walk(HISQ_PATH + name):
                    if files != []:
                        for hvp_nm in files:
                            shutil.copy(HISQ_PATH + name + "/" + hvp_nm, HISQ_RN)
                            os.rename(HISQ_RN + hvp_nm, HISQ_RN + name_1 + "-" + hvp_nm + "." + name_config)
                            if 'hvpAMA' in hvp_nm:
                                hvpAMA_nm.append(name_config)
                                if prt == 1 :
                                    print "copy hvpAMA, config " + name_config
                            if 'hvpExact' in hvp_nm:
                                hvpExact_nm.append(name_config)
                                if prt == 1 :
                                    print "copy hvpExact, config " + name_config
                            if 'hvpSUB' in hvp_nm:
                                hvpSUB_nm.append(name_config)
                                if prt == 1 :
                                    print "copy hvpSUB, config " + name_config
    elif L == 64:
        for name in confignm:
            if 'l64144f211b672m0024m024m286a' in name:
                name_1 = name.split(".")[0]
                name_config = name.split(".")[-1]
                for root, dirs, files in os.walk(HISQ_PATH + name):
                    if files != []:
                        for hvp_nm in files:
                            shutil.copy(HISQ_PATH + name + "/" + hvp_nm, HISQ_RN)
                            os.rename(HISQ_RN + hvp_nm, HISQ_RN + name_1 + "-" + hvp_nm + "." + name_config)
                            if 'hvpAMA' in hvp_nm:
                                hvpAMA_nm.append(name_config)
                                if prt == 1 :
                                    print "copy hvpAMA, config " + name_config
                            if 'hvpExact' in hvp_nm:
                                hvpExact_nm.append(name_config)
                                if prt == 1 :
                                    print "copy hvpExact, config " + name_config
                            if 'hvpSUB' in hvp_nm:
                                hvpSUB_nm.append(name_config)
                                if prt == 1 :
                                    print "copy hvpSUB, config " + name_config
    elif L == 96:
        for name in confignm:
            if 'l96192f211b672m0008m022m260a' in name:
                name_1 = name.split(".")[0]
                name_config = name.split(".")[-1]
                for root, dirs, files in os.walk(HISQ_PATH + name):
                    if files != []:
                        for hvp_nm in files:
                            shutil.copy(HISQ_PATH + name + "/" + hvp_nm, HISQ_RN)
                            os.rename(HISQ_RN + hvp_nm, HISQ_RN + name_1 + "-" + hvp_nm + "." + name_config)
                            if 'hvpAMA' in hvp_nm:
                                hvpAMA_nm.append(name_config)
                                if prt == 1 :
                                    print "copy hvpAMA, config " + name_config
                            if 'hvpExact' in hvp_nm:
                                hvpExact_nm.append(name_config)
                                if prt == 1 :
                                    print "copy hvpExact, config " + name_config
                            if 'hvpSUB' in hvp_nm:
                                hvpSUB_nm.append(name_config)
                                if prt == 1 :
                                    print "copy hvpSUB, config " + name_config
    hvp_nm = list(set(hvpAMA_nm) & set(hvpExact_nm) & set(hvpSUB_nm))
    if prt == 1:
        print len(hvp_nm)
    return sorted(hvp_nm), sorted(hvpAMA_nm), sorted(hvpExact_nm), sorted(hvpSUB_nm)


def all_mode_averaging(L, path, hvp_confignm, outpath, prt = 0):
    u = []
    v = []
    nx = []
    ny = []
    nz = []
    nt = []
    pi_uv = []
    pi_uv_AMA = []
    pi_uv_Exact = []
    pi_uv_SUB = []

    i = 0
    if int(L) == 48:
        for config_nm in hvp_confignm:
            pi_uv_AMA.append([])
            pi_uv_Exact.append([])
            pi_uv_SUB.append([])
            pi_uv.append([])

            data = open(path + 'l48144f211b672m0048m024m286a-hvpAMA.2000.dat.0.' + str(config_nm))
            for line in data.readlines():
                data_fix = line
                data_fix = data_fix.split()
                if i == 0:
                    u.append(int(data_fix[1]))
                    v.append(int(data_fix[2]))
                    nx.append(int(data_fix[4]))
                    ny.append(int(data_fix[5]))
                    nz.append(int(data_fix[6]))
                    nt.append(int(data_fix[7]))
                pi_uv_AMA[i].append(float(data_fix[8]))

            data = open(path + 'l48144f211b672m0048m024m286a-hvpExact.2000.dat.0.' + str(config_nm))
            for line in data.readlines():
                data_fix = line
                data_fix = data_fix.split()
                pi_uv_Exact[i].append(float(data_fix[8]))

            data = open(path + 'l48144f211b672m0048m024m286a-hvpSUB.2000.dat.0.' + str(config_nm))
            for line in data.readlines():
                data_fix = line
                data_fix = data_fix.split()
                pi_uv_SUB[i].append(float(data_fix[8]))

            pi_uv[i] = np.array(pi_uv_AMA[i]) + np.array(pi_uv_Exact[i]) - np.array(pi_uv_SUB[i])

            outdata = open(outpath + 'l48144f211b672m0048m024m286a.2000.dat.0.' + str(config_nm), 'w')
            for ii in range(0, len(u)):
                outdata.write('untVACPOL0 ' + str(u[ii]) + ' ' + str(v[ii]) + ' mom ' + str(nx[ii]) + ' ' + str(ny[ii]) + ' ' + str(nz[ii]) + ' ' + str(nt[ii]) + ' ' + str(pi_uv[i][ii]) + '\n')
            outdata.close()

            i += 1
        if prt == 1:
            print 'mu:'
            print u
            print 'nu:'
            print v
            print 'nx:'
            print nx
            print 'ny:'
            print ny
            print 'nz:'
            print nz
            print 'nt:'
            print nt
            print 'Pi_munu:'
            print np.array(pi_uv)
        return u, v, nx, ny, nz, nt, np.array(pi_uv)


def all_mode_averaging_64(L, path, hvp_confignm, outpath, prt = 0):
    u = []
    v = []
    nx = []
    ny = []
    nz = []
    nt = []
    pi_uv = []
    pi_uv_AMA = []
    pi_uv_Exact = []
    pi_uv_SUB = []

    i = 0
    if int(L) == 64:
        for config_nm in hvp_confignm:
            pi_uv_AMA.append([])
            pi_uv_Exact.append([])
            pi_uv_SUB.append([])
            pi_uv.append([])
            pi_uv_AMA_temp = []
            pi_uv_Exact_temp = []
            pi_uv_SUB_temp = []

            for filenames in walkfiles(path, prt = 0)[1]:
                filenames_config = filenames.split(".")[-1]
                if ('hvpAMA' in filenames) & (str(config_nm) in filenames_config):

                    pi_uv_AMA_temp.append([])
                    data = open(path + filenames)
                    for line in data.readlines():
                        data_fix = line
                        data_fix = data_fix.split()
                        if i == 0:
                            u.append(int(data_fix[1]))
                            v.append(int(data_fix[2]))
                            nx.append(int(data_fix[4]))
                            ny.append(int(data_fix[5]))
                            nz.append(int(data_fix[6]))
                            nt.append(int(data_fix[7]))
                        pi_uv_AMA_temp[-1].append(float(data_fix[8]))
                    i = 1
                elif ('hvpExact' in filenames) & (str(config_nm) in filenames_config):
                    pi_uv_Exact_temp.append([])
                    data = open(path + filenames)
                    for line in data.readlines():
                        data_fix = line
                        data_fix = data_fix.split()

                        pi_uv_Exact_temp[-1].append(float(data_fix[8]))
                elif ('hvpSUB' in filenames) & (str(config_nm) in filenames_config):
                    pi_uv_SUB_temp.append([])
                    data = open(path + filenames)
                    for line in data.readlines():
                        data_fix = line
                        data_fix = data_fix.split()

                        pi_uv_SUB_temp[-1].append(float(data_fix[8]))


            pi_uv_AMA[-1] = np.average(pi_uv_AMA_temp, axis=0)
            pi_uv_Exact[-1] = np.average(pi_uv_Exact_temp, axis=0)
            pi_uv_SUB[-1] = np.average(pi_uv_SUB_temp, axis=0)


            pi_uv[-1] = np.array(pi_uv_AMA[-1]) + np.array(pi_uv_Exact[-1]) - np.array(pi_uv_SUB[-1])
            print len(u)
            print len(pi_uv[-1])

            outdata = open(outpath + 'l64144f211b672m0024m024m286a.3000.dat.0.' + str(config_nm), 'w')
            for ii in range(0, len(u)):
                outdata.write('untVACPOL0 ' + str(u[ii]) + ' ' + str(v[ii]) + ' mom ' + str(nx[ii]) + ' ' + str(ny[ii]) + ' ' + str(nz[ii]) + ' ' + str(nt[ii]) + ' ' + str(pi_uv[-1][ii]) + '\n')
            outdata.close()

            i += 1
        if prt == 1:
            print 'mu:'
            print u
            print 'nu:'
            print v
            print 'nx:'
            print nx
            print 'ny:'
            print ny
            print 'nz:'
            print nz
            print 'nt:'
            print nt
            print 'Pi_munu:'
            print np.array(pi_uv)
        return u, v, nx, ny, nz, nt, np.array(pi_uv)


def all_mode_averaging_96(L, path, hvp_confignm, outpath, prt = 0):
    u = []
    v = []
    nx = []
    ny = []
    nz = []
    nt = []
    pi_uv = []
    pi_uv_AMA = []
    pi_uv_Exact = []
    pi_uv_SUB = []

    i = 0
    if int(L) == 96:
        for config_nm in hvp_confignm:
            pi_uv_AMA.append([])
            pi_uv_Exact.append([])
            pi_uv_SUB.append([])
            pi_uv.append([])
            pi_uv_AMA_temp = []
            pi_uv_Exact_temp = []
            pi_uv_SUB_temp = []

            for filenames in walkfiles(path, prt = 0)[1]:
                filenames_config = filenames.split(".")[-1]
                if ('hvpAMA' in filenames) & (str(config_nm) in filenames_config):

                    pi_uv_AMA_temp.append([])
                    data = open(path + filenames)
                    for line in data.readlines():
                        data_fix = line
                        data_fix = data_fix.split()
                        if i == 0:
                            u.append(int(data_fix[1]))
                            v.append(int(data_fix[2]))
                            nx.append(int(data_fix[4]))
                            ny.append(int(data_fix[5]))
                            nz.append(int(data_fix[6]))
                            nt.append(int(data_fix[7]))
                        pi_uv_AMA_temp[-1].append(float(data_fix[8]))
                    i = 1
                elif ('hvpExact' in filenames) & (str(config_nm) in filenames_config):
                    pi_uv_Exact_temp.append([])
                    data = open(path + filenames)
                    for line in data.readlines():
                        data_fix = line
                        data_fix = data_fix.split()

                        pi_uv_Exact_temp[-1].append(float(data_fix[8]))
                elif ('hvpSUB' in filenames) & (str(config_nm) in filenames_config):
                    pi_uv_SUB_temp.append([])
                    data = open(path + filenames)
                    for line in data.readlines():
                        data_fix = line
                        data_fix = data_fix.split()

                        pi_uv_SUB_temp[-1].append(float(data_fix[8]))


            pi_uv_AMA[-1] = np.average(pi_uv_AMA_temp, axis=0)
            pi_uv_Exact[-1] = np.average(pi_uv_Exact_temp, axis=0)
            pi_uv_SUB[-1] = np.average(pi_uv_SUB_temp, axis=0)


            pi_uv[-1] = np.array(pi_uv_AMA[-1]) + np.array(pi_uv_Exact[-1]) - np.array(pi_uv_SUB[-1])
            print len(u)
            print len(pi_uv[-1])

            outdata = open(outpath + 'l96192f211b672m0008m022m260a.2000.dat.0.' + str(config_nm), 'w')
            for ii in range(0, len(u)):
                outdata.write('untVACPOL0 ' + str(u[ii]) + ' ' + str(v[ii]) + ' mom ' + str(nx[ii]) + ' ' + str(ny[ii]) + ' ' + str(nz[ii]) + ' ' + str(nt[ii]) + ' ' + str(pi_uv[-1][ii]) + '\n')
            outdata.close()

            i += 1
        if prt == 1:
            print 'mu:'
            print u
            print 'nu:'
            print v
            print 'nx:'
            print nx
            print 'ny:'
            print ny
            print 'nz:'
            print nz
            print 'nt:'
            print nt
            print 'Pi_munu:'
            print np.array(pi_uv)
        return u, v, nx, ny, nz, nt, np.array(pi_uv)



def read_config(L, path, prt = 0):
    u, v, nx, ny, nz, nt, pi_uv = [], [], [], [], [], [], []
    i = 0
    for filenames in walkfiles(path, prt = 0)[1]:
        if (L == 96) & ('l96192f211b672m0008m022m260a' in filenames):
            data = open(path + filenames)
            pi_uv.append([])
            for line in data.readlines():
                data_fix = line
                data_fix = data_fix.split()
                if i == 0:
                    u.append(int(data_fix[1]))
                    v.append(int(data_fix[2]))
                    nx.append(int(data_fix[4]))
                    ny.append(int(data_fix[5]))
                    nz.append(int(data_fix[6]))
                    nt.append(int(data_fix[7]))
                pi_uv[i].append(float(data_fix[8]))
            i += 1
        elif (L == 64) & ('l64144f211b672m0024m024m286a' in filenames):
            data = open(path + filenames)
            pi_uv.append([])
            for line in data.readlines():
                data_fix = line
                data_fix = data_fix.split()
                if i == 0:
                    u.append(int(data_fix[1]))
                    v.append(int(data_fix[2]))
                    nx.append(int(data_fix[4]))
                    ny.append(int(data_fix[5]))
                    nz.append(int(data_fix[6]))
                    nt.append(int(data_fix[7]))
                pi_uv[i].append(float(data_fix[8]))
            i += 1
        elif (L == 48) & ('l48144f211b672m0048m024m286a' in filenames):
            data = open(path + filenames)
            pi_uv.append([])
            for line in data.readlines():
                data_fix = line
                data_fix = data_fix.split()
                if i == 0:
                    u.append(int(data_fix[1]))
                    v.append(int(data_fix[2]))
                    nx.append(int(data_fix[4]))
                    ny.append(int(data_fix[5]))
                    nz.append(int(data_fix[6]))
                    nt.append(int(data_fix[7]))
                pi_uv[i].append(float(data_fix[8]))
            i += 1
    if prt == 1:
        print 'mu:'
        print np.array(u)
        print 'nu:'
        print np.array(v)
        print 'nx:'
        print np.array(nx)
        print 'ny:'
        print np.array(ny)
        print 'nz:'
        print np.array(nz)
        print 'nt:'
        print np.array(nt)
        print 'Pi_munu:'
        print np.array(pi_uv)
    return u, v, nx, ny, nz, nt, np.array(pi_uv)


def func_q(n, L):
    return 2.0 * math.sin(math.pi * n / L)
    # return 2.0*math.pi*n/L


def func_qsqu(nx, ny, nz, nt, L):
    return (2.0 * math.sin(math.pi * nx / L[0])) ** 2.0 + (2.0 * math.sin(math.pi * ny / L[1])) ** 2.0 + (2.0 * math.sin(math.pi * nz / L[2])) ** 2.0 + (2.0 * math.sin(math.pi * nt / L[3])) ** 2.0
    # return (2.0*math.pi)**2.0*((nx/L[0])**2.0+(ny/L[1])**2.0+(nz/L[2])**2.0+(nt/L[3])**2.0)


def func_delta(u, v):
    if u == v:
        return 1.0
    else:
        return 0.0


def func_qsqu_original_pi_qsqu_original(u, v, nx, ny, nz, nt, pi_uv, L, ainv, prt = 0):
    ii = 0
    nu = 0
    lu = 0
    nv = 0
    lv = 0
    index = []
    qsqu = []
    pi_qsqu = []
    line_num = len(u)
    for i in range(0, line_num):
        if u[i] == 0:
            nu = nx[i]
            lu = L[0]
        elif u[i] == 1:
            nu = ny[i]
            lu = L[1]
        elif u[i] == 2:
            nu = nz[i]
            lu = L[2]
        elif u[i] == 3:
            nu = nt[i]
            lu = L[3]
        if v[i] == 0:
            nv = nx[i]
            lv = L[0]
        elif v[i] == 1:
            nv = ny[i]
            lv = L[1]
        elif v[i] == 2:
            nv = nz[i]
            lv = L[2]
        elif v[i] == 3:
            nv = nt[i]
            lv = L[3]
        temp_qsqu_original = func_qsqu(nx[i], ny[i], nz[i], nt[i], L)
        temp = temp_qsqu_original * func_delta(u[i], v[i]) - func_q(nu, lu) * func_q(nv, lv)

        if abs(temp) > 10.0 ** (-5.0):
            qsqu.append(temp_qsqu_original)
            pi_qsqu.append([])
            index.append(i)
            for j in range(0, len(pi_uv)):
                pi_qsqu[ii].append(pi_uv[j][i] / temp / 4.0)
            ii += 1
    pi_qsqu = (np.array(pi_qsqu)).T
    if prt == 1:
        print 'original Q^2:'
        print qsqu
        print 'original Pi(Q^2):'
        print pi_qsqu
        print 'lenth of original Q^2:'
        print len(qsqu)
        print
    return ainv ** 2.0 * np.array(qsqu), pi_qsqu


def func_qsqu_original_pi_qsqu_original_spltq2(u, v, nx, ny, nz, nt, pi_uv, L, ainv, prt = 0):
    ii = 0
    nu = 0
    lu = 0
    nv = 0
    lv = 0
    index = []
    qsqu = []
    pi_qsqu = []
    line_num = len(u)
    for i in range(0, line_num):
        if u[i] == 0:
            nu = nx[i]
            lu = L[0]
        elif u[i] == 1:
            nu = ny[i]
            lu = L[1]
        elif u[i] == 2:
            nu = nz[i]
            lu = L[2]
        elif u[i] == 3:
            nu = nt[i]
            lu = L[3]
        if v[i] == 0:
            nv = nx[i]
            lv = L[0]
        elif v[i] == 1:
            nv = ny[i]
            lv = L[1]
        elif v[i] == 2:
            nv = nz[i]
            lv = L[2]
        elif v[i] == 3:
            nv = nt[i]
            lv = L[3]
        temp_qsqu_original = func_qsqu(nx[i], ny[i], nz[i], nt[i], L)
        temp = temp_qsqu_original * func_delta(u[i], v[i]) - func_q(nu, lu) * func_q(nv, lv)

        if abs(temp) > 10.0 ** (-5.0):
            temp_qsqu_original = ainv ** 2.0 * temp_qsqu_original
            if nt[i] != 0:
                temp_qsqu_original = - temp_qsqu_original
            qsqu.append(temp_qsqu_original)
            pi_qsqu.append([])
            index.append(i)
            for j in range(0, len(pi_uv)):
                pi_qsqu[ii].append(pi_uv[j][i] / temp / 4.0)
            ii += 1
    pi_qsqu = (np.array(pi_qsqu)).T
    if prt == 1:
        print 'original Q^2:'
        print qsqu
        print 'original Pi(Q^2):'
        print pi_qsqu
        print 'lenth of original Q^2:'
        print len(qsqu)
        print
    return np.array(qsqu), pi_qsqu


def func_irrep_extract(u, v, nx, ny, nz, nt, pi_uv, uv, prt = 0):
    u_new = []
    v_new = []
    nx_new = []
    ny_new = []
    nz_new = []
    nt_new = []
    pi_uv_new = []
    for i in range (0, len(u)):
        if (u[i] == uv[0]) & (v[i] == uv[1]):
            u_new.append(u[i])
            v_new.append(v[i])
            nx_new.append(nx[i])
            ny_new.append(ny[i])
            nz_new.append(nz[i])
            nt_new.append(nt[i])
            pi_uv_new.append((pi_uv.T)[i])
    pi_uv_new = (np.array(pi_uv_new)).T
    if prt == 1:
        print 'mu = ', uv[0], 'nu = ', uv[1]
        print 'lenth:'
        print len(u_new)
        print 'mu:'
        print u_new
        print 'nu'
        print v_new
        print 'nx'
        print nx_new
        print 'ny'
        print ny_new
        print 'nz'
        print nz_new
        print 'nt'
        print nt_new
        print 'Pi_', uv[0], uv[1]
        print pi_uv_new
    return u_new, v_new, nx_new, ny_new, nz_new, nt_new, pi_uv_new


def func_mix_v_m(vector, matrix):
    v = np.array([vector])
    matrix_array = np.array(matrix)
    return np.concatenate((v, matrix_array), axis=0)


def func_order_byrow(matrix):
    matrix_sort = np.array(sorted(matrix.T, key=lambda row : row[0]))
    return matrix_sort.T


def func_split_v_m(matrix):
    matrix_array = np.array(matrix)
    vector, matrix_seperated = np.vsplit(matrix_array, [1])
    return vector[0], matrix_seperated


def func_same_qsqu(vector1, matrix, prt = 0):
    matrix_T = matrix.T
    line_num = len(vector1)
    num_same = 1
    vector1_same = []
    matrix_same = []
    vector1_same.append(vector1[0])
    matrix_same.append([])
    matrix_same[0] = matrix_T[0]
    ii = 1
    for i in range(1, line_num):
        if vector1[i] == vector1[i - 1]:
            matrix_same[ii - 1] = np.add(matrix_same[ii - 1], matrix_T[i])
            num_same += 1
        else:
            vector1_same.append(vector1[i])
            matrix_same[ii - 1] = matrix_same[ii - 1] / num_same
            matrix_same[ii - 1] = matrix_same[ii - 1].tolist()
            matrix_same.append([])
            matrix_same[ii] = matrix_T[i]
            num_same = 1
            ii += 1
    matrix_same[ii - 1] = matrix_same[ii - 1] / num_same
    matrix_same[ii - 1] = matrix_same[ii - 1].tolist()
    matrix_same = (np.array(matrix_same)).T
    if prt == 1:
        print 'Q^2:'
        print vector1_same
        print 'lenth of Q^2:'
        print len(vector1_same)
        print 'Pi(Q^2)'
        print matrix_same
        print 'number of configurations:'
        print len(matrix_same)
        print
    return vector1_same, matrix_same


def c(matrix):
    matrix_t = matrix.T
    t, n = matrix_t.shape
    matrix_t_ave = np.mean(matrix_t, axis = 1)
    c_matrix = np.zeros((t, t))
    for i in range(0, t):
        for j in range(0, t):
            # c_matrix[i, j] = 1.0 / (n - 1.0) * 1.0 / n * sum((matrix_t[i, :] - matrix_t_ave[i]) * (matrix_t[j, :] - matrix_t_ave[j]))
            configsum = 0
            for k in range(0, n):
                configsum += (matrix_t[i, k] - matrix_t_ave[i]) * (matrix_t[j, k] - matrix_t_ave[j])
            c_matrix[i, j] = 1.0 / (n - 1.0) * 1.0 / n * configsum
    return c_matrix


def func_pi(parameters, x):
    m = len(parameters)
    function_pi_sum_part = 0
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


def chi_squ(parameters, qsqu, pi_qsqu_avg, c_inverse):
    d_data_theory = np.array([np.array(pi_qsqu_avg) - np.array(func_pi(parameters, qsqu))])
    return d_data_theory.dot(c_inverse.dot(d_data_theory.T))[0][0]


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



def get_pi_ii_qsqu_pi_qsqu(u, v, nx, ny, nz, nt, pi_uv, L, ainv):
    pi_ii_org_qsqu, pi_ii_org_pi_qsqu = [], []
    for i in [0, 1, 2]:
        pi_ii_u_temp, pi_ii_v_temp, pi_ii_nx_temp, pi_ii_ny_temp, pi_ii_nz_temp, pi_ii_nt_temp, pi_ii_pi_uv_temp = func_irrep_extract(u, v, nx, ny, nz, nt, pi_uv, [i, i])
        pi_ii_org_qsqu_temp, pi_ii_org_pi_qsqu_temp = func_qsqu_original_pi_qsqu_original(pi_ii_u_temp, pi_ii_v_temp, pi_ii_nx_temp, pi_ii_ny_temp, pi_ii_nz_temp, pi_ii_nt_temp, pi_ii_pi_uv_temp, L, ainv)
        if pi_ii_org_pi_qsqu == []:
            pi_ii_org_qsqu = pi_ii_org_qsqu_temp
            pi_ii_org_pi_qsqu = pi_ii_org_pi_qsqu_temp
        else:
            pi_ii_org_qsqu = np.concatenate((pi_ii_org_qsqu, pi_ii_org_qsqu_temp), axis=1)
            pi_ii_org_pi_qsqu = np.concatenate((pi_ii_org_pi_qsqu, pi_ii_org_pi_qsqu_temp), axis=1)
    pi_ii_qsqu_original_pi_qsqu_original_mix = func_mix_v_m(pi_ii_org_qsqu, pi_ii_org_pi_qsqu)
    pi_ii_qsqu_original_pi_qsqu_original_mix_ordered = func_order_byrow(pi_ii_qsqu_original_pi_qsqu_original_mix)
    pi_ii_qsqu_ordered, pi_ii_pi_qsqu_ordered = func_split_v_m(pi_ii_qsqu_original_pi_qsqu_original_mix_ordered)
    pi_ii_qsqu, pi_ii_pi_qsqu = func_same_qsqu(pi_ii_qsqu_ordered, pi_ii_pi_qsqu_ordered)
    return pi_ii_qsqu, pi_ii_pi_qsqu


def get_pi_ii_qsqu_pi_qsqu_spltq2(u, v, nx, ny, nz, nt, pi_uv, L, ainv):
    pi_ii_org_qsqu, pi_ii_org_pi_qsqu = [], []
    for i in [0, 1, 2]:
        pi_ii_u_temp, pi_ii_v_temp, pi_ii_nx_temp, pi_ii_ny_temp, pi_ii_nz_temp, pi_ii_nt_temp, pi_ii_pi_uv_temp = func_irrep_extract(u, v, nx, ny, nz, nt, pi_uv, [i, i])
        pi_ii_org_qsqu_temp, pi_ii_org_pi_qsqu_temp = func_qsqu_original_pi_qsqu_original_spltq2(pi_ii_u_temp, pi_ii_v_temp, pi_ii_nx_temp, pi_ii_ny_temp, pi_ii_nz_temp, pi_ii_nt_temp, pi_ii_pi_uv_temp, L, ainv)
        if pi_ii_org_pi_qsqu == []:
            pi_ii_org_qsqu = pi_ii_org_qsqu_temp
            pi_ii_org_pi_qsqu = pi_ii_org_pi_qsqu_temp
        else:
            pi_ii_org_qsqu = np.concatenate((pi_ii_org_qsqu, pi_ii_org_qsqu_temp), axis=1)
            pi_ii_org_pi_qsqu = np.concatenate((pi_ii_org_pi_qsqu, pi_ii_org_pi_qsqu_temp), axis=1)
    pi_ii_qsqu_original_pi_qsqu_original_mix = func_mix_v_m(pi_ii_org_qsqu, pi_ii_org_pi_qsqu)
    pi_ii_qsqu_original_pi_qsqu_original_mix_ordered = func_order_byrow(pi_ii_qsqu_original_pi_qsqu_original_mix)
    pi_ii_qsqu_ordered, pi_ii_pi_qsqu_ordered = func_split_v_m(pi_ii_qsqu_original_pi_qsqu_original_mix_ordered)
    pi_ii_qsqu, pi_ii_pi_qsqu = func_same_qsqu(pi_ii_qsqu_ordered, pi_ii_pi_qsqu_ordered)

    return_pi_qsqu(pi_ii_qsqu)
    pi_ii_qsqu_pi_qsqu_mix = func_mix_v_m(pi_ii_qsqu, pi_ii_pi_qsqu)
    pi_ii_qsqu_pi_qsqu_mix_ordered = func_order_byrow(pi_ii_qsqu_pi_qsqu_mix)
    pi_ii_qsqu, pi_ii_pi_qsqu = func_split_v_m(pi_ii_qsqu_pi_qsqu_mix_ordered)



    return pi_ii_qsqu, pi_ii_pi_qsqu


def return_pi_qsqu(pi_qsqu):
    for i in range(0, len(pi_qsqu)):
        if pi_qsqu[i] < 0.0:
            pi_qsqu[i] = -pi_qsqu[i]


def get_pi_e_qsqu_pi_qsqu(u, v, nx, ny, nz, nt, pi_uv, L, ainv):
    pi_ii_org_qsqu, pi_ii_org_pi_qsqu = [], []
    for i in [0, 1, 2]:
        pi_ii_u_temp, pi_ii_v_temp, pi_ii_nx_temp, pi_ii_ny_temp, pi_ii_nz_temp, pi_ii_nt_temp, pi_ii_pi_uv_temp = func_irrep_extract(u, v, nx, ny, nz, nt, pi_uv, [i, i])
        pi_ii_org_qsqu_temp, pi_ii_org_pi_qsqu_temp = func_qsqu_original_pi_qsqu_original(pi_ii_u_temp, pi_ii_v_temp, pi_ii_nx_temp, pi_ii_ny_temp, pi_ii_nz_temp, pi_ii_nt_temp, pi_ii_pi_uv_temp, L, ainv)
        if pi_ii_org_pi_qsqu == []:
            pi_ii_org_qsqu = pi_ii_org_qsqu_temp
            pi_ii_org_pi_qsqu = pi_ii_org_pi_qsqu_temp
        else:
            pi_ii_org_qsqu = np.concatenate((pi_ii_org_qsqu, pi_ii_org_qsqu_temp), axis=1)
            pi_ii_org_pi_qsqu = np.concatenate((pi_ii_org_pi_qsqu, pi_ii_org_pi_qsqu_temp), axis=1)
    pi_ii_qsqu_original_pi_qsqu_original_mix = func_mix_v_m(pi_ii_org_qsqu, pi_ii_org_pi_qsqu)
    pi_ii_qsqu_original_pi_qsqu_original_mix_ordered = func_order_byrow(pi_ii_qsqu_original_pi_qsqu_original_mix)
    pi_ii_qsqu_ordered, pi_ii_pi_qsqu_ordered = func_split_v_m(pi_ii_qsqu_original_pi_qsqu_original_mix_ordered)
    pi_ii_qsqu, pi_ii_pi_qsqu = func_same_qsqu(pi_ii_qsqu_ordered, pi_ii_pi_qsqu_ordered)
    return pi_ii_qsqu, pi_ii_pi_qsqu



def get_pi00_q2_piq2(u, v, nt, nx, ny, nz, pi_uv, L, ainv):
    pi_00_u, pi_00_v, pi_00_nx, pi_00_ny, pi_00_nz, pi_00_nt, pi_00_pi_uv = func_irrep_extract(u, v, nx, ny, nz, nt, pi_uv, [3, 3], prt = 0)
    pi_00_org_qsqu, pi_00_org_pi_qsqu = func_qsqu_original_pi_qsqu_original(pi_00_u, pi_00_v, pi_00_nx, pi_00_ny, pi_00_nz, pi_00_nt, pi_00_pi_uv, L, ainv, prt = 0)
    pi_00_qsqu_original_pi_qsqu_original_mix = func_mix_v_m(pi_00_org_qsqu, pi_00_org_pi_qsqu)
    pi_00_qsqu_original_pi_qsqu_original_mix_ordered = func_order_byrow(pi_00_qsqu_original_pi_qsqu_original_mix)
    pi_00_qsqu_ordered, pi_00_pi_qsqu_ordered = func_split_v_m(pi_00_qsqu_original_pi_qsqu_original_mix_ordered)
    pi_00_qsqu, pi_00_pi_qsqu = func_same_qsqu(pi_00_qsqu_ordered, pi_00_pi_qsqu_ordered, prt = 0)
    return pi_00_qsqu, pi_00_pi_qsqu




def get_pi_4i_qsqu_pi_qsqu(u, v, nx, ny, nz, nt, pi_uv, L, ainv):
    pi_4i_org_qsqu, pi_4i_org_pi_qsqu = [], []
    for i in [0, 1, 2]:
        pi_4i_u_temp, pi_4i_v_temp, pi_4i_nx_temp, pi_4i_ny_temp, pi_4i_nz_temp, pi_4i_nt_temp, pi_4i_pi_uv_temp = func_irrep_extract(u, v, nx, ny, nz, nt, pi_uv, [3, i])
        pi_4i_org_qsqu_temp, pi_4i_org_pi_qsqu_temp = func_qsqu_original_pi_qsqu_original(pi_4i_u_temp, pi_4i_v_temp, pi_4i_nx_temp, pi_4i_ny_temp, pi_4i_nz_temp, pi_4i_nt_temp, pi_4i_pi_uv_temp, L, ainv)
        if pi_4i_org_pi_qsqu == []:
            pi_4i_org_qsqu = pi_4i_org_qsqu_temp
            pi_4i_org_pi_qsqu = pi_4i_org_pi_qsqu_temp
        else:
            pi_4i_org_qsqu = np.concatenate((pi_4i_org_qsqu, pi_4i_org_qsqu_temp), axis=1)
            pi_4i_org_pi_qsqu = np.concatenate((pi_4i_org_pi_qsqu, pi_4i_org_pi_qsqu_temp), axis=1)
    for i in [0, 1, 2]:
        pi_4i_u_temp, pi_4i_v_temp, pi_4i_nx_temp, pi_4i_ny_temp, pi_4i_nz_temp, pi_4i_nt_temp, pi_4i_pi_uv_temp = func_irrep_extract(u, v, nx, ny, nz, nt, pi_uv, [i, 3])
        pi_4i_org_qsqu_temp, pi_4i_org_pi_qsqu_temp = func_qsqu_original_pi_qsqu_original(pi_4i_u_temp, pi_4i_v_temp, pi_4i_nx_temp, pi_4i_ny_temp, pi_4i_nz_temp, pi_4i_nt_temp, pi_4i_pi_uv_temp, L, ainv)
        if pi_4i_org_pi_qsqu == []:
            pi_4i_org_qsqu = pi_4i_org_qsqu_temp
            pi_4i_org_pi_qsqu = pi_4i_org_pi_qsqu_temp
        else:
            pi_4i_org_qsqu = np.concatenate((pi_4i_org_qsqu, pi_4i_org_qsqu_temp), axis=1)
            pi_4i_org_pi_qsqu = np.concatenate((pi_4i_org_pi_qsqu, pi_4i_org_pi_qsqu_temp), axis=1)
    pi_4i_qsqu_original_pi_qsqu_original_mix = func_mix_v_m(pi_4i_org_qsqu, pi_4i_org_pi_qsqu)
    pi_4i_qsqu_original_pi_qsqu_original_mix_ordered = func_order_byrow(pi_4i_qsqu_original_pi_qsqu_original_mix)
    pi_4i_qsqu_ordered, pi_4i_pi_qsqu_ordered = func_split_v_m(pi_4i_qsqu_original_pi_qsqu_original_mix_ordered)
    pi_4i_qsqu, pi_4i_pi_qsqu = func_same_qsqu(pi_4i_qsqu_ordered, pi_4i_pi_qsqu_ordered)
    return pi_4i_qsqu, pi_4i_pi_qsqu



def get_pi_ij_qsqu_pi_qsqu(u, v, nx, ny, nz, nt, pi_uv, L, ainv):
    org_qsqu, org_pi_qsqu = [], []
    for i in [0, 1, 2]:
        sub = list(set([0, 1, 2]) - set([i]))
        for j in sub:
            u_temp, v_temp, nx_temp, ny_temp, nz_temp, nt_temp, pi_uv_temp = func_irrep_extract(u, v, nx, ny, nz, nt, pi_uv, [i, j])
            org_qsqu_temp, org_pi_qsqu_temp = func_qsqu_original_pi_qsqu_original(u_temp, v_temp, nx_temp, ny_temp, nz_temp, nt_temp, pi_uv_temp, L, ainv)
            if org_pi_qsqu == []:
                org_qsqu = org_qsqu_temp
                org_pi_qsqu = org_pi_qsqu_temp
            else:
                org_qsqu = np.concatenate((org_qsqu, org_qsqu_temp), axis=1)
                org_pi_qsqu = np.concatenate((org_pi_qsqu, org_pi_qsqu_temp), axis=1)
    qsqu_original_pi_qsqu_original_mix = func_mix_v_m(org_qsqu, org_pi_qsqu)
    qsqu_original_pi_qsqu_original_mix_ordered = func_order_byrow(qsqu_original_pi_qsqu_original_mix)
    qsqu_ordered, pi_qsqu_ordered = func_split_v_m(qsqu_original_pi_qsqu_original_mix_ordered)
    qsqu, pi_qsqu = func_same_qsqu(qsqu_ordered, pi_qsqu_ordered)
    return qsqu, pi_qsqu



def get_pi_qsqu_avg_cmatrix_output(path, L, irrep, qsqu, pi_qsqu, pointsnm, jkblk = 0, prt = 0):
    if jkblk == 0:
        pi_qsqu_avg = np.mean(pi_qsqu, axis=0)
        c_matrix = c(pi_qsqu)
        write_q2_piq2_cmatrix_file(path, L, irrep, qsqu, pi_qsqu_avg, c_matrix, pointsnm, jkblk = 0)
        if prt == 1:
            print 'Pi(Q^2)_avg for' + irrep + ':'
            print pi_qsqu_avg, len(pi_qsqu_avg)
            print 'lenth of Pi(Q^2)_avg'
            print 'c_matrix for' + irrep + ':'
            print c_matrix
            print
        return pi_qsqu_avg, c_matrix
    else:
        for i in range(0, len(pi_qsqu)):
            pi_qsqu_del = np.delete(pi_qsqu, i, 0)
            pi_qsqu_avg_del = np.mean(pi_qsqu_del, axis=0)
            c_matrix_del = c(pi_qsqu_del)
            write_q2_piq2_cmatrix_file(path, L, irrep + '_jkblk' + str(i), qsqu, pi_qsqu_avg_del, c_matrix_del, pointsnm, jkblk = 1)



def write_q2_piq2_cmatrix_file(path, L, name, qsqu, pi_qsqu_avg, c_matrix, pointsnm, jkblk = 0):
    if L == 48:
        file = open(path + 'HISQ_l48144f211b672m0048m024m286a_' + name, 'w')
    elif L == 64:
        file = open(path + 'HISQ_l64144f211b672m0024m024m286a_' + name, 'w')
    elif L == 96:
        file = open(path + 'HISQ_l96192f211b672m0008m022m260a_' + name, 'w')
    file.write('AVERAGES:\n')
    for i in range(0, pointsnm):
        file.write(str(qsqu[i]) + ' ' + str(pi_qsqu_avg[i]) + '\n')
    file.write('c_matrix:\n')
    for i in range(0, pointsnm):
        for j in range (0, pointsnm):
            file.write(str(qsqu[i]) + ' ' + str(qsqu[j]) + ' ' + str(c_matrix[i][j]) + '\n')
    file.close()
    if L == 48:
        print 'output' + ' ' + 'HISQ_l48144f211b672m0048m024m286a_' + name + ' successfully!'
    if L == 64:
        print 'output' + ' ' + 'HISQ_l64144f211b672m0024m024m286a_' + name + ' successfully!'



def get_pi00_sub_q2_piq2(u, v, nt, nx, ny, nz, pi_uv, L, ainv):



    pi_00_u, pi_00_v, pi_00_nx, pi_00_ny, pi_00_nz, pi_00_nt, pi_00_pi_uv = func_irrep_extract(u, v, nx, ny, nz, nt, pi_uv, [3, 3], prt = 0)
    pi_00_org_qsqu, pi_00_org_pi_qsqu = func_qsqu_original_pi_qsqu_original(pi_00_u, pi_00_v, pi_00_nx, pi_00_ny, pi_00_nz, pi_00_nt, pi_00_pi_uv, L, ainv, prt = 0)
    pi_00_qsqu_original_pi_qsqu_original_mix = func_mix_v_m(pi_00_org_qsqu, pi_00_org_pi_qsqu)
    pi_00_qsqu_original_pi_qsqu_original_mix_ordered = func_order_byrow(pi_00_qsqu_original_pi_qsqu_original_mix)
    pi_00_qsqu_ordered, pi_00_pi_qsqu_ordered = func_split_v_m(pi_00_qsqu_original_pi_qsqu_original_mix_ordered)
    pi_00_qsqu, pi_00_pi_qsqu = func_same_qsqu(pi_00_qsqu_ordered, pi_00_pi_qsqu_ordered, prt = 0)
    return pi_00_qsqu, pi_00_pi_qsqu



def subtract(u, v, nx, ny, nz, nt, pi_uv, L):
    print pi_uv
    pi_uv_T = pi_uv.T
    confignm = len(pi_uv_T[0])
    pi_uv_T_test = []
    uu = []
    vv = []
    nxx = []
    nyy = []
    nzz = []
    ntt = []
    pi_uuvv_T = []
    for i in range(1, len(u)/16):

        for uuu in range(0, 4):
            for vvv in range(0, 4):
                uu.append(uuu)
                vv.append(vvv)
                nxx.append(nx[i * 16 + uuu * 4 + vvv])
                nyy.append(ny[i * 16 + uuu * 4 + vvv])
                nzz.append(nz[i * 16 + uuu * 4 + vvv])
                ntt.append(nt[i * 16 + uuu * 4 + vvv])
                pi_uuvv_T_temp = np.zeros(confignm)
                for k in range(0, 4):
                    for l in range(0, 4):

                        pi_uv_T_test.append(pi_uv_T[i * 16 + k * 4 + l][0])
                        #print 'pi_uuvv_T_temp'
                        #print pi_uuvv_T_temp
                        #print 'sub'
                        #print transproj(uuu, k, nx[i * 16], ny[i * 16], nz[i * 16], nt[i * 16], L)
                        #print transproj(uuu, k, nx[i * 16], ny[i * 16], nz[i * 16], nt[i * 16], L) * (np.array(pi_uv_T[i * 16 + k * 4 + l]) - np.array(pi_uv_T[k * 4 + l])) * transproj(l, vvv, nx[i * 16], ny[i * 16], nz[i * 16], nt[i * 16], L)
                        pi_uuvv_T_temp += transproj(uuu, k, nx[i * 16], ny[i * 16], nz[i * 16], nt[i * 16], L) * (np.array(pi_uv_T[i * 16 + k * 4 + l]) - np.array(pi_uv_T[k * 4 + l])) * transproj(l, vvv, nx[i * 16], ny[i * 16], nz[i * 16], nt[i * 16], L)
                pi_uuvv_T.append(pi_uuvv_T_temp)
    pi_uuvv = (np.array(pi_uuvv_T)).T



    return uu, vv, nxx, nyy, nzz, ntt, pi_uuvv



def subtract_no_p(u, v, nx, ny, nz, nt, pi_uv, L):
    pi_uv_T = pi_uv.T
    pi_uv_T_test = []
    uu = []
    vv = []
    nxx = []
    nyy = []
    nzz = []
    ntt = []
    pi_uuvv_T = []
    for i in range(1, len(u)/16):

        for uuu in range(0, 4):
            for vvv in range(0, 4):
                uu.append(uuu)
                vv.append(vvv)
                nxx.append(nx[i * 16 + uuu * 4 + vvv])
                nyy.append(ny[i * 16 + uuu * 4 + vvv])
                nzz.append(nz[i * 16 + uuu * 4 + vvv])
                ntt.append(nt[i * 16 + uuu * 4 + vvv])
                # print np.array(pi_uv_T[i * 16 + uuu * 4 + vvv]) - np.array(pi_uv_T[uuu * 4 + vvv])
                pi_uuvv_T.append(np.array(pi_uv_T[i * 16 + uuu * 4 + vvv]) - np.array(pi_uv_T[uuu * 4 + vvv]))
    pi_uuvv = (np.array(pi_uuvv_T)).T
    return uu, vv, nxx, nyy, nzz, ntt, pi_uuvv



def transproj(u, v, nx, ny, nz, nt, L):

    if u == 0:
        nu = nx
        lu = L[0]
    elif u == 1:
        nu = ny
        lu = L[1]
    elif u == 2:
        nu = nz
        lu = L[2]
    elif u == 3:
        nu = nt
        lu = L[3]
    if v == 0:
        nv = nx
        lv = L[0]
    elif v == 1:
        nv = ny
        lv = L[1]
    elif v == 2:
        nv = nz
        lv = L[2]
    elif v == 3:
        nv = nt
        lv = L[3]

    return func_delta(u, v) - func_q(nu, lu) * func_q(nv, lv) / func_qsqu(nx, ny, nz, nt, L)

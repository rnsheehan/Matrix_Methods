# Import libraries
# You should try an import the bare minimum of modules
import sys # access system routines
import os
import glob
import re

import math
import scipy
import numpy as np
import matplotlib.pyplot as plt

# add path to our file
sys.path.append('c:/Users/Robert/Programming/Python/Common/')
sys.path.append('c:/Users/Robert/Programming/Python/Plotting/')

import Common
import Plotting

MOD_NAME_STR = "Plots" # use this in exception handling messages

def iface_r_t():
    # make a plot of the computed interface reflection / transmission curves
    # R. Sheehan 15 - 6 - 2020

    FUNC_NAME = ".iface_r_t()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        filename = "Air_to_Glass_R_T.txt"
        #filename = "Glass_to_Air_R_T.txt"
            
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', skiprows = 4, unpack = True)

            hv_data = []; labels = []; marks = [];
            hv_data.append([data[0], data[1]]); labels.append('$r_{TE}$'); marks.append(Plotting.labs_lins[0]); 
            hv_data.append([data[0], data[3]]); labels.append('$r_{TM}$'); marks.append(Plotting.labs_dashed[0]); 
            hv_data.append([data[0], data[2]]); labels.append('$t_{TE}$'); marks.append(Plotting.labs_lins[1]);
            hv_data.append([data[0], data[4]]); labels.append('$t_{TM}$'); marks.append(Plotting.labs_dashed[1]); 
                            
            # make the plot of the data set
            args = Plotting.plot_arg_multiple()

            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'Input Angle (rad)'
            args.y_label = 'Reflection / Transmission'
            args.fig_name = filename.replace('.txt','')

            Plotting.plot_multiple_curves(hv_data, args)

            del hv_data; del labels; del marks; 
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def layer_r_t():
    # make a plot of the computed dielectric layer reflection / transmission curves
    # R. Sheehan 17 - 6 - 2020

    FUNC_NAME = ".layer_r_t()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        #filename = "Air_Silicon_Silica_R_T.txt"
        #filename = "Air_Silicon_Silica_R_T_Alt.txt"
        #filename = "Air_Silica_Silicon_R_T.txt"
        filename = "Air_Silica_Silicon_R_T_Alt.txt"
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)

            lscale = (1.0/1550.0)

            hv_data = []; labels = []; marks = [];
            hv_data.append([data[0]*lscale, data[1]]); labels.append('$R_{TE}$'); marks.append(Plotting.labs_lins[0]); 
            hv_data.append([data[0]*lscale, data[3]]); labels.append('$R_{TM}$'); marks.append(Plotting.labs_dashed[0]); 
            hv_data.append([data[0]*lscale, data[2]]); labels.append('$T_{TE}$'); marks.append(Plotting.labs_lins[1]);
            hv_data.append([data[0]*lscale, data[4]]); labels.append('$T_{TM}$'); marks.append(Plotting.labs_dashed[1]);            
                            
            # make the plot of the data set
            args = Plotting.plot_arg_multiple()

            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'Layer Thickness / Wavelength'
            args.y_label = 'Reflectivity / Transmissivity'
            args.plt_range = [data[0][0]*lscale, data[0][-1]*lscale, 0.0, 1.0]
            args.fig_name = filename.replace('.txt','')

            Plotting.plot_multiple_curves(hv_data, args)

            del hv_data; del labels; del marks; 
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def layer_r_t_alt_compare():
    # make a plot of the computed dielectric layer reflection / transmission curves
    # R. Sheehan 17 - 6 - 2020

    FUNC_NAME = ".layer_r_t()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        #filename = "Air_Silicon_Silica_R_T.txt"
        #filename = "Air_Silicon_Silica_R_T_Alt.txt"
        filename1 = "Air_Silica_Silicon_R_T.txt"
        filename2 = "Air_Silica_Silicon_R_T_Alt.txt"
        if glob.glob(filename1) and glob.glob(filename2):
            # import the dataset
            data1 = np.loadtxt(filename1, delimiter = ',', unpack = True)
            data2 = np.loadtxt(filename2, delimiter = ',', unpack = True)

            lscale = (1.0/1550.0)

            hv_data = []; labels = []; marks = [];
            hv_data.append([data1[0]*lscale, data1[1]]); labels.append('$R_{TE}$'); marks.append(Plotting.labs_lins[0]); 
            hv_data.append([data1[0]*lscale, data1[2]]); labels.append('$T_{TE}$'); marks.append(Plotting.labs_lins[1]);
            hv_data.append([data2[0]*lscale, data2[1]]); labels.append('$R_{TE}^{alt}$'); marks.append(Plotting.labs_dashed[0]); 
            hv_data.append([data2[0]*lscale, data2[2]]); labels.append('$T_{TE}^{alt}$'); marks.append(Plotting.labs_dashed[1]);            
                            
            # make the plot of the data set
            args = Plotting.plot_arg_multiple()

            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'Layer Thickness / Wavelength'
            args.y_label = 'Reflectivity / Transmissivity'
            args.plt_range = [data1[0][0]*lscale, data1[0][-1]*lscale, 0.0, 1.0]
            args.fig_name = filename1.replace('.txt','') + '_Compar'

            Plotting.plot_multiple_curves(hv_data, args)

            del hv_data; del labels; del marks; 
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def layer_AR():
    # make a plot of the computed dielectric layer AR coating reflection curve
    # R. Sheehan 17 - 6 - 2020

    FUNC_NAME = ".layer_AR()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        filename = "Air_MgF2_Glass.txt"
        #filename = "Air_SiO2_SiN.txt"
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)

            hv_data = []; labels = []; marks = [];
            #hv_data.append([data[0], data[1]]); labels.append('T = 267.24 nm'); marks.append(Plotting.labs_lins[0]);
            hv_data.append([data[0], data[1]]); labels.append('T = 101.86 nm'); marks.append(Plotting.labs_lins[0]); 
            hv_data.append([data[0], data[3]]); labels.append('2 T'); marks.append(Plotting.labs_lins[1]); 
            #hv_data.append([data[0], data[2]]); labels.append('3 T'); marks.append(Plotting.labs_lins[2]);
            #hv_data.append([data[0], data[4]]); labels.append('4 T'); marks.append(Plotting.labs_lins[3]);
            hv_data.append([data[0], data[5]]); labels.append('5 T'); marks.append(Plotting.labs_lins[4]); 
                            
            # make the plot of the data set
            args = Plotting.plot_arg_multiple()

            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'Wavelength nm'
            args.y_label = 'Reflectivity'
            args.plt_range = [data[0][0], data[0][-1], 0.0, 0.05]
            args.fig_name = filename.replace('.txt','')

            Plotting.plot_multiple_curves(hv_data, args)

            del hv_data; del labels; del marks; 
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def layer_HR():
    # make a plot of the computed dielectric layer HR coating reflection curve
    # R. Sheehan 17 - 6 - 2020

    FUNC_NAME = ".layer_HR()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        filename = "HR_Coating_15.txt"
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)

            hv_data = []; labels = []; marks = [];
            hv_data.append([data[0], data[1]]); labels.append('R'); marks.append(Plotting.labs_lins[0]);
            hv_data.append([data[0], data[2]]); labels.append('T'); marks.append(Plotting.labs_lins[1]); 
                            
            # make the plot of the data set
            args = Plotting.plot_arg_multiple()

            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'Wavelength nm'
            args.y_label = 'Reflectivity / Transmissivity'
            args.plt_range = [data[0][0], data[0][-1], 0.0, 1.0]
            args.fig_name = filename.replace('.txt','')

            Plotting.plot_multiple_curves(hv_data, args)

            del hv_data; del labels; del marks; 
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def layer_HR_2():
    # make a plot of the computed dielectric layer HR reflection curves
    # R. Sheehan 19 - 6 - 2020

    FUNC_NAME = ".layer_HR_2()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        hv_data = []; labels = []; marks = [];

        files = ["HR_Coating_7.txt","HR_Coating_11.txt","HR_Coating_15.txt"]

        count = 7;
        i=0
        for f in files:
            if glob.glob(f):
                data = np.loadtxt(f, delimiter = ',', unpack = True)
                hv_data.append([data[0], data[1]]); labels.append('Layers = %(v1)d'%{"v1":count}); marks.append(Plotting.labs_lins[i])
            count = count + 4
            i = i + 1 

        # make the plot of the data set
        args = Plotting.plot_arg_multiple()

        args.loud = True
        args.crv_lab_list = labels
        args.mrk_list = marks
        args.x_label = 'Wavelength nm'
        args.y_label = 'Reflectivity'
        args.plt_range = [data[0][0], data[0][-1], 0.0, 1.0]
        args.fig_name = 'HR_Coating'

        Plotting.plot_multiple_curves(hv_data, args)

        del hv_data; del labels; del marks; 
       
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def layer_BP():
    # make a plot of the computed dielectric BP filter transmission curve
    # R. Sheehan 22 - 6 - 2020

    FUNC_NAME = ".layer_BP()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        filename = "BP_Filter_3.txt"
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)

            hv_data = []; labels = []; marks = [];
            hv_data.append([data[0], data[1]]); labels.append('R'); marks.append(Plotting.labs_lins[0]);
            hv_data.append([data[0], data[2]]); labels.append('T'); marks.append(Plotting.labs_lins[1]); 
                            
            # make the plot of the data set
            args = Plotting.plot_arg_multiple()

            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'Wavelength nm'
            args.y_label = 'Reflectivity / Transmissivity'
            args.plt_range = [data[0][0], data[0][-1], 0.0, 1.0]
            args.fig_name = filename.replace('.txt','')

            Plotting.plot_multiple_curves(hv_data, args)

            del hv_data; del labels; del marks; 
            
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def layer_BP_2():
    # make a plot of the computed dielectric layer HR reflection curves
    # R. Sheehan 19 - 6 - 2020

    FUNC_NAME = ".layer_BP_2()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        hv_data = []; labels = []; marks = [];

        files = ["BP_Filter_1.txt","BP_Filter_2.txt","BP_Filter_3.txt","BP_Filter_4.txt"]

        count = 1;
        i=0
        for f in files:
            if glob.glob(f):
                data = np.loadtxt(f, delimiter = ',', unpack = True)
                hv_data.append([data[0], data[2]]); labels.append('Layers = %(v1)d'%{"v1":count}); marks.append(Plotting.labs_lins[i])
            count = count + 1
            i = i + 1 

        # make the plot of the data set
        args = Plotting.plot_arg_multiple()

        args.loud = True
        args.crv_lab_list = labels
        args.mrk_list = marks
        args.x_label = 'Wavelength nm'
        args.y_label = 'Transmissivity'
        args.plt_range = [data[0][0], data[0][-1], 0.0, 1.0]
        args.fig_name = 'BP_Filter'

        Plotting.plot_multiple_curves(hv_data, args)

        del hv_data; del labels; del marks; 
       
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def layer_BP_BW():
    # make a plot of the computed dielectric BP filter transmission curve bandwidth
    # R. Sheehan 22 - 6 - 2020

    FUNC_NAME = ".layer_BP_BW()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        # the data
        x = [1, 2, 3, 4, 5]
        y = [151.70, 22.48, 3.86, 0.67, 0.12]
                        
        # make the plot of the data set
        args = Plotting.plot_arg_single()

        args.loud = True
        #args.curve_label = labels
        args.marker = Plotting.labs_mrk_only[3]
        args.x_label = 'No. Layer Pairs'
        args.y_label = 'Transmission BW (nm)'
        #args.plt_range = [data[0][0], data[0][-1], 0.0, 1.0]
        args.fig_name = 'BP_Filter_BW'

        Plotting.plot_single_curve(x, y, args)
        
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def main():
    pass

if __name__ == '__main__':
    main()

    pwd = os.getcwd() # get current working directory

    #print(pwd)
    
    #iface_r_t()

    #layer_r_t()

    layer_r_t_alt_compare()

    #layer_AR()

    #layer_HR_2()

    #layer_BP()

    #layer_BP_2()

    #layer_BP_BW()

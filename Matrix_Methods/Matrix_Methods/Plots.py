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
        filename = "Air_Silicon_Silica_R_T.txt"
        #filename = "Air_Silica_Silicon_R_T.txt"
        if glob.glob(filename):
            # import the dataset
            data = np.loadtxt(filename, delimiter = ',', unpack = True)

            hv_data = []; labels = []; marks = [];
            hv_data.append([data[0]/1550.0, data[1]]); labels.append('$R_{TE}$'); marks.append(Plotting.labs_lins[0]); 
            hv_data.append([data[0]/1550.0, data[3]]); labels.append('$R_{TM}$'); marks.append(Plotting.labs_dashed[0]); 
            hv_data.append([data[0]/1550.0, data[2]]); labels.append('$T_{TE}$'); marks.append(Plotting.labs_lins[1]);
            hv_data.append([data[0]/1550.0, data[4]]); labels.append('$T_{TM}$'); marks.append(Plotting.labs_dashed[1]); 
                            
            # make the plot of the data set
            args = Plotting.plot_arg_multiple()

            args.loud = True
            args.crv_lab_list = labels
            args.mrk_list = marks
            args.x_label = 'Layer Thickness / Wavelength'
            args.y_label = 'Reflectivity / Transmissivity'
            args.plt_range = [data[0][0]/1550.0, data[0][-1]/1550.0, 0.0, 1.0]
            args.fig_name = filename.replace('.txt','')

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
            hv_data.append([data[0], data[1]]); labels.append('T = 287.04 nm'); marks.append(Plotting.labs_lins[0]); 
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
    # make a plot of the computed dielectric layer AR coating reflection curve
    # R. Sheehan 17 - 6 - 2020

    FUNC_NAME = ".layer_HR()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        filename = "HR_Coating.txt"
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

def main():
    pass

if __name__ == '__main__':
    main()

    pwd = os.getcwd() # get current working directory

    #print(pwd)
    
    #iface_r_t()

    #layer_r_t()

    #layer_AR()

    layer_HR()

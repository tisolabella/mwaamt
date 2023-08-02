#####################################################################
# IMPORTS
#####################################################################

import json, sys, inspect, datetime, warnings, pickle
import numpy as np
from csv import writer
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
from pathlib import Path

from mwaamt.datahandle import data_read
from mwaamt.funk import singlefit, doublefit
from mwaamt.funk import set_resolution, get_chisq
from mwaamt.funk import average, stddev, weighted_average






#---- Wrapper function

def run(configuration_file_path):

    #------ Get the starting date and time
    start_time = datetime.datetime.now()

    warnings.filterwarnings("ignore")






    #####################################################################
    # CONSTANTS & PRINTOUT STRINGS
    #####################################################################
    MISSING_KEYWORD = "ERROR: keyword not present in the configuration file"







    #####################################################################
    # GET THE INPUT
    #####################################################################


    # Input config file as CLI argument
    print("\n---> Opening configuration file...\n")
    try:
        with open(configuration_file_path, 'r') as f:
            cfg = json.load(f)
    except Exception as e:
        print(e)

    # Read the data with the method from datahandle.py
    data = data_read(configuration_file_path)







    #####################################################################
    # PRELIMINARY FITTING AND CORRELATIONS
    #####################################################################


    #----------- Check for booked presets:
    try:
        levo_booked = True if 'levoglucosan' in cfg['presets'] else False
    except NameError as e:
        print(MISSING_KEYWORD, e)
    # Fill in with other presets
    ##### TODO implement user-defined fits and correlations (for the GUI)

    # The valus in the cfg file are the default to be used if there is
    # no preliminary analysis
    alpha_BC, alpha_FF, alpha_WB = cfg['alpha BC'], cfg['alpha FF'], cfg['alpha WB']
    best_alpha_BC, best_alpha_FF, best_alpha_WB = alpha_BC, alpha_FF, alpha_WB

    #----------- Carry out the LEVOGLUCOSAN CORRELATION ADJUSTMENT
    if levo_booked:
        print(f'---> Performing correlation maximisation with'
                ' "levoglucosan" preset...\n')
        # Set two parameters
        for iteration_number in range(cfg['iterations']):
            print(f"\t*Iteration {iteration_number + 1}")
            # Set the search resolution
            alpha_BC_set, alpha_FF_set, alpha_WB_set = set_resolution(best_alpha_BC, 
                    best_alpha_FF, best_alpha_WB, iteration=iteration_number)
            #------ Do the alpha_BC iteration -------------
            # List for the correlations
            R_2_alpha_BC = []
            BC_correlation_pairs = {}
            for alpha_BC in alpha_BC_set:
                # List to save the BrC(lambda_short) and levo
                BrC_set, levo_set = [], []
                # Fix the alpha BC value in a new function
                def typefit(x, A, B, alpha_BrC):
                    return doublefit(x, A, alpha_BC, B, alpha_BrC)
                # Iterate over all the samples
                for sample in data:
                    prop = sample.properties
                    # Do fit
                    try:
                        fitres = curve_fit(typefit, prop.wavelength,
                                prop.abs, p0=(1e3, 1e10, 3),
                                bounds=([1, 1, 1], [1e15, 1e15, 10]),
                                sigma=prop.u_abs)
                        B = fitres[0][1]
                        alpha_BrC = fitres[0][2]
                        # Apportion BrC at the shortest wavelength
                        lambda_short =  min(prop.wavelength) 
                        BrC = B * lambda_short ** ( - alpha_BrC) 
                        BrC_set.append(BrC)
                    except Exception as e:
                        print(f'FIT ERROR for ALPHA_BC: {e}')
                        BrC_set.append(0)
                    levo_set.append(prop.Levoglucosan)
                # Calculate regression and append R^2
                try:
                    regression_res = linregress(levo_set, y=BrC_set)
                    R_2_alpha_BC.append(regression_res.rvalue ** 2)
                    BC_correlation_pairs[alpha_BC] = regression_res.rvalue ** 2
                    # For the article
                    #print(f'alpha_BC {alpha_BC}: slope {regression_res.slope} inter {regression_res.intercept} R2 {regression_res.rvalue ** 2}')
                except Exception as e:
                    print(f'REGRESSION ERROR for ALPHA BC: {e}')
                    R_2_alpha_BC.append(0)
                    BC_correlation_pairs[alpha_BC] = 0
            max_R2_index = R_2_alpha_BC.index(max(R_2_alpha_BC))
            # For stability, use the best only if its significantly
            # better than the previous best:
            if R_2_alpha_BC[max_R2_index] - BC_correlation_pairs[best_alpha_BC] < cfg['threshold']:
                pass # Do not change the value
            else:
                best_alpha_BC = alpha_BC_set[max_R2_index]
            #For the article
            #if iteration_number == cfg['iterations'] - 1: # Final lap
                #print(f'alpha_BC {alpha_BC_set}\nR2 {R_2_alpha_BC}')

            #------ Do the alpha_FF iteration -------------
            # List for the correlations
            R_2_alpha_FF = []
            FF_correlation_pairs = {}
            for alpha_FF in alpha_FF_set:
                # List to save the BC_WB(lambda_long) and levo
                BC_WB_set, levo_set = [], []
                # Fix the alpha BC value in a new function 
                def typefit(x, A, B, alpha_BrC):
                    return doublefit(x, A, best_alpha_BC, B, alpha_BrC)
                # Fix the alpha FF and alpha WB values in a new function 
                def sourcefit(x, A_p, B_p):
                    return doublefit(x, A_p, alpha_FF, B_p, alpha_WB)
                # Iterate over all the samples
                for sample in data:
                    prop = sample.properties
                    # Do fit
                    try:
                        typefitres = curve_fit(typefit, prop.wavelength,
                                prop.abs, p0=(1e3, 1e10, 3),
                                bounds=([1, 1, 1], [1e15, 1e15, 10]),
                                sigma=prop.u_abs)
                        sourcefitres = curve_fit(sourcefit, prop.wavelength,
                                prop.abs, p0=(1e3, 1e10),
                                bounds=([1, 1], [1e15, 1e15]),
                                sigma=prop.u_abs)
                        A = typefitres[0][0]
                        A_p = sourcefitres[0][0]
                        # Apportion BC WB at the longest wavelength
                        lambda_long =  max(prop.wavelength) 
                        BC_WB = (A - A_p) * lambda_long ** ( - best_alpha_BC) 
                        BC_WB_set.append(BC_WB)
                    except Exception as e:
                        print(f'FIT ERROR for ALPHA_FF: {e}')
                        BC_WB_set.append(0)
                    levo_set.append(prop.Levoglucosan)
                # Calculate regression and append R^2
                try:
                    regression_res = linregress(levo_set, y=BC_WB_set)
                    R_2_alpha_FF.append(regression_res.rvalue ** 2)
                    FF_correlation_pairs[alpha_FF] = regression_res.rvalue ** 2
                except Exception as e:
                    print(f'REGRESSION ERROR for ALPHA FF: {e}')
                    R_2_alpha_FF.append(0)
                    FF_correlation_pairs[alpha_FF] = 0
            max_R2_index = R_2_alpha_FF.index(max(R_2_alpha_FF))
            # For stability, use the best only if its significantly
            # better than the previous best:
            if R_2_alpha_FF[max_R2_index] - FF_correlation_pairs[best_alpha_FF] < cfg['threshold']:
                pass # Do not change the value
            else:
                best_alpha_FF = alpha_FF_set[max_R2_index]
            #For the article
            #if iteration_number == cfg['iterations'] - 1: # Final lap
                #print(f'alpha_FF {alpha_FF_set}\nR2 {R_2_alpha_FF}')

            #------ Do the alpha_WB iteration -------------
            # List for the correlations
            R_2_alpha_WB = []
            WB_correlation_pairs = {}
            for alpha_WB in alpha_WB_set:
                # List to save the BC_WB(lambda_long) and levo
                BC_WB_set, levo_set = [], []
                # Fix the alpha BC value in a new function 
                def typefit(x, A, B, alpha_BrC):
                    return doublefit(x, A, best_alpha_BC, B, alpha_BrC)
                # Fix the alpha FF and alpha WB values in a new function 
                def sourcefit(x, A_p, B_p):
                    return doublefit(x, A_p, best_alpha_FF, B_p, alpha_WB)
                # Iterate over all the samples
                for sample in data:
                    prop = sample.properties
                    # Do fit
                    try:
                        typefitres = curve_fit(typefit, prop.wavelength,
                                prop.abs, p0=(1e3, 1e10, 3),
                                bounds=([1, 1, 1], [1e15, 1e15, 10]),
                                sigma=prop.u_abs)
                        sourcefitres = curve_fit(sourcefit, prop.wavelength,
                                prop.abs, p0=(1e3, 1e10),
                                bounds=([1, 1], [1e15, 1e15]),
                                sigma=prop.u_abs)
                        A = typefitres[0][0]
                        A_p = sourcefitres[0][0]
                        # Apportion BC WB at the longest wavelength
                        lambda_long =  max(prop.wavelength) 
                        BC_WB = (A - A_p) * lambda_long ** ( - best_alpha_BC) 
                        BC_WB_set.append(BC_WB)
                    except Exception as e:
                        print(f'FIT ERROR for ALPHA_FF: {e}')
                        BC_WB_set.append(0)
                    levo_set.append(prop.Levoglucosan)
                # Calculate regression and append R^2
                try:
                    regression_res = linregress(levo_set, y=BC_WB_set)
                    R_2_alpha_WB.append(regression_res.rvalue ** 2)
                    WB_correlation_pairs[alpha_WB] = regression_res.rvalue ** 2
                    # For the article
                    #print(f'alpha_WB {alpha_WB}: slope {regression_res.slope} inter {regression_res.intercept} R2 {regression_res.rvalue ** 2}')
                except Exception as e:
                    print(f'REGRESSION ERROR for ALPHA WB: {e}')
                    R_2_alpha_WB.append(0)
                    WB_correlation_pairs[alpha_WB] = 0
            max_R2_index = R_2_alpha_WB.index(max(R_2_alpha_WB))
            # For stability, use the best only if its significantly
            # better than the previous best:
            if R_2_alpha_WB[max_R2_index] - WB_correlation_pairs[best_alpha_WB] < cfg['threshold']:
                pass # Do not change the value
            else:
                best_alpha_WB = alpha_WB_set[max_R2_index]
            #For the article
            #if iteration_number == cfg['iterations'] - 1: # Final lap
                #print(f'alpha_WB {alpha_WB_set}\nR2 {R_2_alpha_WB}')

            #--- Save plots for the correlation vs parameter values
            try:
                if cfg['plots']:
                    print('\t\tsaving plots')
                    fig, ax = plt.subplots()
                    ax.plot(best_alpha_BC, BC_correlation_pairs[best_alpha_BC], 'rx', label='best')
                    ax.plot(alpha_BC_set, R_2_alpha_BC, 'k-', label=r'$\alpha_{BC}$')
                    ax.plot(best_alpha_FF, FF_correlation_pairs[best_alpha_FF], 'rx')
                    ax.plot(alpha_FF_set, R_2_alpha_FF, 'r-', label=r'$\alpha_{FF}$')
                    ax2 = ax.twiny()
                    ax2.plot(best_alpha_WB, WB_correlation_pairs[best_alpha_WB], 'rx')
                    ax2.plot(alpha_WB_set, R_2_alpha_WB, 'b-', label=r'$\alpha_{WB}$')
                    ax.grid(alpha=0.3)
                    ax.set_xlabel(r'Parameter value for $\alpha_{BC}$ and $\alpha_{FF}$')
                    ax2.set_xlabel(r'Parameter value for $\alpha_{WB}$')
                    ax.set_ylabel(r'Levoglucosan analysis $R^2$')
                    fig.legend(bbox_to_anchor=(0.15, 0.15, 0.15, 0.25))
                    plt.tight_layout()
                    direc = cfg['working directory'] + f'plots/preplots/'
                    try:
                        plt.savefig(direc + f'/iter{iteration_number + 1}.png', dpi = 300)
                    except FileNotFoundError as fnfe:
                        Path(direc).mkdir(parents=True, exist_ok=True)
                        plt.savefig(direc + f'/iter{iteration_number + 1}.png', dpi = 300)
                    plt.close()
            except KeyError as ke:
                print(MISSING_KEYWORD, ke)

        #--- From now on, use alpha_XX for the best values
        #--- or the default values if no optimizatios was done
        alpha_BC, alpha_FF, alpha_WB = best_alpha_BC, best_alpha_FF, best_alpha_WB
        print(f"\nThe best parameters for the the correlation with levoglucosan are "
                f"(alpha_BC = {round(alpha_BC, 2)}, alpha_FF = {round(alpha_FF, 2)}, alpha_WB = {round(alpha_WB, 2)}).\n")







    #################################################################
    # ALPHA BrC variation with ALPHA BC (swipe)
    #################################################################

    try:
        if cfg['alpha bc swipe']:
            print("---> Performing alpha bc swipe...\n")
            sw_alpha_BC_set = np.linspace(0.8, 1.2, 20).tolist()
            sw_alpha_BrC_set = []
            sw_alpha_BrC_stddev_set = []
            for abc in sw_alpha_BC_set:
                tmp_abrc_list = []
                for sample in data:
                    prp = sample.properties
                    def typefit(x, A, B, alpha_BrC):
                        """Fix alpha BC to the best value"""
                        return doublefit(x, A, abc, B, alpha_BrC)
                    type_fitres = curve_fit(typefit, prp.wavelength,
                            prp.abs, p0=(1e3, 1e10, 3),
                            bounds=([1, 1, 1], [1e15, 1e15, 10]),
                            sigma=prp.u_abs)
                    # Save alpha_BrC
                    tmp_abrc_list.append(type_fitres[0][2])
                sw_alpha_BrC_set.append(average(tmp_abrc_list))
                sw_alpha_BrC_stddev_set.append(stddev(tmp_abrc_list))
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)

    # For the article
    #print(f'avg {sw_alpha_BrC_set} \nstdev {sw_alpha_BrC_stddev_set}')




    ##################################################################
    # FITTING PROCEDURE WITH OPTIMIZED PARAMETERS
    ##################################################################

    print(f'---> Fitting the experimental data...\n')

    failed_fit_count = 0
    failed_fit = []

    for sample in data:

        #--------- AAE fit
        prp = sample.properties
        try:
            aae_fitres = curve_fit(singlefit, prp.wavelength,
                    prp.abs, p0=(1e5, 1), bounds=([1,0.5], [1e15, 3]),
                    sigma=prp.u_abs)
            # Save the exponent
            prp.aae = aae_fitres[0][1]
            prp.u_aae = np.sqrt(aae_fitres[1][1][1])
            #--------- Second AAE fit fixing the exp for lower uncertainty
            def singlefit_fix(x, scale):
                return singlefit(x, scale, prp.aae)
            second_aae_fitres = curve_fit(singlefit_fix, prp.wavelength,
                    prp.abs, p0=(1e5), bounds=(1, 1e15),
                    sigma=prp.u_abs)
            # Save the scale
            prp.scale = second_aae_fitres[0][0]
            prp.u_scale = np.sqrt(second_aae_fitres[1][0][0])
            # Calculates the chisquared for this fit
            expected = [singlefit_fix(x, prp.scale) for x in prp.wavelength]
            ndf = len(prp.wavelength) - 2
            tmp, prp.red_chisq_aae_fit = get_chisq(prp.abs, expected,
                    prp.u_abs, ndf)
        except RuntimeError as re:
            print(FIT_ERROR, re)
            failed_fit_count += 1
            failed_fit.append(['aae', f'sample {sample.name}'])

        #--------- Components fit
        def typefit(x, A, B, alpha_BrC):
            """Fix alpha BC to the best value"""
            return doublefit(x, A, alpha_BC, B, alpha_BrC)
        try:
            type_fitres = curve_fit(typefit, prp.wavelength,
                    prp.abs, p0=(1e3, 1e10, 3),
                    bounds=([1, 1, 1], [1e15, 1e15, 10]),
                    sigma=prp.u_abs)
            # Save alpha_BrC
            prp.alpha_brc = type_fitres[0][2]
            prp.u_alpha_brc = np.sqrt(type_fitres[1][2][2])
            #------ Improved components fit
            def typefit_fix(x, A, B):
                """Fix alpha BrC for uncertainty improvement"""
                return typefit(x, A, B, prp.alpha_brc)
            second_type_fitres = curve_fit(typefit_fix, prp.wavelength,
                    prp.abs, p0=(1e3, 1e10),
                    bounds=([1, 1], [1e15, 1e15]),
                    sigma=prp.u_abs)
            prp.A = type_fitres[0][0]
            prp.u_A = np.sqrt(type_fitres[1][0][0])
            prp.B = type_fitres[0][1]
            prp.u_B = np.sqrt(type_fitres[1][1][1])
            # Calculates the chisquared for this fit
            expected = [typefit_fix(x, prp.A, prp.B) for x in prp.wavelength]
            ndf = len(prp.wavelength) - 3
            tmp, prp.red_chisq_type_fit = get_chisq(prp.abs, expected,
                    prp.u_abs, ndf)
        except RuntimeError as re:
            print(FIT_ERROR, re)
            failed_fit_count += 1
            failed_fit.append(['type', f'sample {sample.name}'])

        #---- Source fit
        def sourcefit(x, A_p, B_p):
            """Fix alpha_FF and alpha_WB to the best values"""
            return doublefit(x, A_p, alpha_FF, B_p, alpha_WB)
        try:
            source_fitres = curve_fit(sourcefit, prp.wavelength,
                    prp.abs, p0=(1e3, 1e10),
                    bounds=([1, 1], [1e15, 1e15]),
                    sigma=prp.u_abs)
            prp.A_p = source_fitres[0][0]
            prp.u_A_p = np.sqrt(source_fitres[1][0][0])
            prp.B_p = source_fitres[0][1]
            prp.u_B_p = np.sqrt(source_fitres[1][1][1])
            # Calculates the chisquared for this fit
            expected = [sourcefit(x, prp.A_p, prp.B_p) for x in prp.wavelength]
            ndf = len(prp.wavelength) - 2
            tmp, prp.red_chisq_source_fit = get_chisq(prp.abs, expected,
                    prp.u_abs, ndf)
        except RuntimeError as re:
            print(FIT_ERROR, re)
            failed_fit_count += 1
            failed_fit.append(['source', f'sample {sample.name}'])






    ####################################################################
    # OPTICAL APPORTIONMENT
    ####################################################################

    print("---> Performing optical apportionment...\n")
    for sample in data:
        prp = sample.properties
        A = prp.A
        B = prp.B
        A_p = prp.A_p
        B_p = prp.B_p
        alpha_BrC = prp.alpha_brc
        prp.bc_wb, prp.bc_wb_frac = [], []
        prp.bc_ff, prp.bc_ff_frac = [], []
        prp.brc, prp.brc_frac = [], []
        for w, a in zip(prp.wavelength, prp.abs):
            #--- BC wood burning
            value = (A - A_p) / (w ** alpha_BC)
            prp.bc_wb.append(value)
            prp.bc_wb_frac.append(value / a)
            #--- BC fossil fuel
            value = A_p / (w ** alpha_FF)
            prp.bc_ff.append(value)
            prp.bc_ff_frac.append(value / a)
            #--- BrC 
            value = B / (w ** alpha_BrC)
            prp.brc.append(value)
            prp.brc_frac.append(value / a)

                 



    ###################################################################
    # MASS APPORTIONMENT
    ###################################################################

    # Check whether k1 and k2 are provided 
    if cfg['k1'] == 0 or cfg['k2'] == 0:
        do_fit = True
    else:
        k1 = cfg['k1']
        k2 = cfg['k2']
        do_fit = False

    print("---> Performing mass apportionment...\n")
    #--- EC apportionment
    for sample in data:
        prp = sample.properties
        # Get the longest and shortest wavelengths
        # in case the list is not in order
        lambda_long = max(prp.wavelength)
        i_l = prp.wavelength.index(lambda_long)
        # Apportion EC
        prp.ec_ff = prp.ec * (prp.bc_ff[i_l] / (prp.abs[i_l] - prp.brc[i_l]))
        prp.ec_wb = prp.ec * (prp.bc_wb[i_l] / (prp.abs[i_l] - prp.brc[i_l]))
    #--- OC apportionment
    # Create and populate the lists for the fits
    k1_list, k2_list = [], []
    k1_x, k1_y = [], []
    k2_x, k2_y = [], []
    for sample in data:
        prp = sample.properties
        # Get the longest and shortest wavelengths
        # in case the list is not in order
        lambda_long = max(prp.wavelength)
        i_l = prp.wavelength.index(lambda_long)
        if prp.aae < cfg['AAE high'] and prp.aae > cfg['AAE low']:
            k1_list.append(sample)
            k1_x.append(prp.bc_ff[i_l])
            k1_y.append(prp.oc)
        else:
            k2_list.append(sample)
    # Perform fit if needed
    if do_fit:
        fit_1 = linregress(k1_x, y=k1_y)
        k1 = fit_1.slope
        u_k1 = fit_1.stderr
        int1 = fit_1.intercept
    # Calculate OC_FF for all samples using the k1 just found
    for sample in data:
        prp = sample.properties
        lambda_long = max(prp.wavelength)
        i_l = prp.wavelength.index(lambda_long)
        prp.oc_ff = k1 * prp.bc_ff[i_l]
    # Create the x and y lists for the second regression
    for sample in k2_list:
        prp = sample.properties
        lambda_short = min(prp.wavelength)
        i_s = prp.wavelength.index(lambda_short)
        k2_x.append(prp.brc[i_s])
        k2_y.append(prp.oc - prp.oc_ff)
    # Perform fit if needed
    if do_fit:
        fit_2 = linregress(k2_x, y=k2_y)
        k2 = fit_2.slope
        u_k2 = fit_2.stderr
        int2 = fit_2.intercept
    # Calculate OC_WB and OC_NC for all samples using the k1 just found
    for sample in data:
        prp = sample.properties
        lambda_short = min(prp.wavelength)
        i_s = prp.wavelength.index(lambda_short)
        prp.oc_wb = k2 * prp.brc[i_s]
        #input(f'name {sample.name} lambda_short {lambda_short} brc {prp.brc[i_s]} k2 {k2} oc_wb {prp.oc_wb}')
        prp.oc_nc = prp.oc - prp.oc_ff - prp.oc_wb




    ###################################################################
    # FIT PARAMETERS WRITEOUT
    ###################################################################

    try:
        out_path = cfg['working directory'] + 'fitres.csv'
        print(f"---> Writing fit results to {out_path}\n")
        with open(out_path, 'w') as f:
            writa = writer(f)
            header = ['Name', 'Scale', 'eScale', 'AAE', 
                    'eAAE', 'Red_chisq', 'A', 'eA', 'B', 
                    'eB', 'alpha_BrC', 'ealpha_BrC', 'Red_chisq', 
                    'A\'', 'eA\'', 'B\'', 'eB\'', 'Red_chisq']
            writa.writerow(header)
            for sample in data:
                prp = sample.properties
                linetowrite = [sample.name, prp.scale, prp.u_scale, 
                        prp.aae, prp.u_aae, prp.red_chisq_aae_fit, 
                        prp.A, prp.u_A, prp.B, prp.u_B, prp.alpha_brc,
                        prp.u_alpha_brc, prp.red_chisq_type_fit, 
                        prp.A_p, prp.u_A_p, prp.B_p, prp.u_B_p, 
                        prp.red_chisq_source_fit]
                writa.writerow(linetowrite) 
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)
    except FileNotFoundError as fnfe:
        Path(cfg['working directory']).mkdir(parents=True, exist_ok=True)







    ###################################################################
    # OPTICAL APPORTIONMENT RESULTS WRITEOUT
    ###################################################################

    # Use the first sample to write header
    header_1 = ['',]
    for w in data[0].properties.wavelength:
        header_1.append(f'{w} nm')
        header_1.append('')
        header_1.append('')
    header_1 = header_1 + header_1[1:]
    header_2 = ['Name',]
    for i in range(len(data[0].properties.wavelength)):
        header_2 += ['BC_FF_frac', 'BC_WB_frac', 'BrC_frac']
    for i in range(len(data[0].properties.wavelength)):
        header_2 += ['BC_FF', 'BC_WB', 'BrC']

    try:
        out_path = cfg['working directory'] + 'appres.csv'
        print(f"---> Writing optical apportionment"
                f" results to {out_path}\n")
        with open(out_path, 'w') as f:
            writa = writer(f)
            writa.writerow(header_1)
            writa.writerow(header_2)
            for sample in data:
                prp = sample.properties
                line_to_write = [sample.name,]
                zippo = zip(prp.bc_ff_frac, prp.bc_wb_frac, prp.brc_frac)
                for bcff_frac, bcwb_frac, brc_frac in zippo:
                    line_to_write.append(bcff_frac)
                    line_to_write.append(bcwb_frac)
                    line_to_write.append(brc_frac)
                zippo = zip(prp.bc_ff, prp.bc_wb, prp.brc)
                for bcff, bcwb, brc in zippo:
                    line_to_write.append(bcff)
                    line_to_write.append(bcwb)
                    line_to_write.append(brc)
                writa.writerow(line_to_write)
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)
    except FileNotFoundError as fnfe:
        Path(cfg['working directory']).mkdir(parents=True, exist_ok=True)




    ###################################################################
    # MASS APPORTIONMENT RESULTS WRITEOUT
    ###################################################################

    # Use the first sample to write header
    header = ['Name', 'EC_FF', 'EC_WB', 'OC_FF', 'OC_WB', 'OC_NC']

    try:
        out_path = cfg['working directory'] + 'mappres.csv'
        print(f"---> Writing mass apportionment"
                f" results to {out_path}\n")
        with open(out_path, 'w') as f:
            writa = writer(f)
            writa.writerow(header)
            for sample in data:
                prp = sample.properties
                line_to_write = [sample.name,]
                line_to_write.append(prp.ec_ff)
                line_to_write.append(prp.ec_wb)
                line_to_write.append(prp.oc_ff)
                line_to_write.append(prp.oc_wb)
                line_to_write.append(prp.oc_nc)
                writa.writerow(line_to_write)
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)
    except FileNotFoundError as fnfe:
        Path(cfg['working directory']).mkdir(parents=True, exist_ok=True)



    #####################################################################
    # SAVE PLOTS
    #####################################################################

    #---- Save fit plots plots for the two fits if booked
    try:
        if cfg['plots']:
            print(f'---> Saving fit plots in {cfg["working directory"]}' + f'plots/fitplots/' + '\n') 
            for sample in data:
                prp = sample.properties
                # Data points
                plt.errorbar(prp.wavelength, prp.abs, 
                      xerr=prp.u_wavelength, yerr=prp.u_abs,
                      fmt='.k', elinewidth=0.8, markersize=1.2,
                      label=f'data {sample.name}')
                x = np.linspace(prp.wavelength[0], prp.wavelength[-1], 500)
                # Type fit
                ytype = typefit_fix(x, prp.A, prp.B)
                plt.plot(x, ytype, 'r', 
                      label='Component fit (BC + BrC)')
                ybc = singlefit(x, prp.A, alpha_BC)
                plt.plot(x, ybc, 'k', linewidth=0.5,
                      label='BC contribution')
                ybrc = singlefit(x, prp.B, prp.alpha_brc)
                plt.plot(x, ybrc, 'y', linewidth=0.5,
                        label='BrC contribution')
                # Source fit
                ysource = sourcefit(x, prp.A_p, prp.B_p)
                plt.plot(x, ysource, 'b', linestyle='dashed', 
                      label='Source fit (FF + WB)')
                yff = singlefit(x, prp.A_p, alpha_FF)
                plt.plot(x, yff, 'm', linewidth=0.5, linestyle='dashed',
                      label='FF contribution')
                ywb = singlefit(x, prp.B_p, alpha_WB)
                plt.plot(x, ywb, 'g', linewidth=0.5, linestyle='dashed',
                        label='WB contribution')
                plt.xlabel('Wavelength [nm]')
                plt.ylabel(r'Absorption coefficient, $b_{abs}$  ' + '[Mm'+ r'$^{-1}$' + ']') if prp.data_type == 'Babs' else plt.ylabel('100 ABS')
                plt.grid(alpha=0.3)
                plt.legend()
                plt.tight_layout()
                direc = cfg['working directory'] + f'plots/fitplots/'
                try:
                    plt.savefig(direc+f'{sample.name}.png', dpi = 300)
                except FileNotFoundError as fnfe:
                    Path(direc).mkdir(parents=True, exist_ok=True)
                    plt.savefig(direc+f'{sample.name}.png', dpi = 300)
                plt.close()
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)

    #----- Save plots for the variation of alpha_brc if booked
    try:
        if cfg['plots']:
            print(f'---> Saving alpha plot in {cfg["working directory"]}' + f'plots/fitplots/' + '\n') 
            # Lists to store values for plotting
            alpha, error, names = [], [], []
            for sample in data:
                prp = sample.properties
                alpha.append(prp.alpha_brc)
                error.append(prp.u_alpha_brc)
                names.append(sample.name)
            fig, ax = plt.subplots() 
            ax.errorbar(names, alpha, xerr=None, yerr=error, fmt='.r')
            ax.set_xticklabels(names, rotation=75)
            ax.grid(alpha=0.3)
            ax.set_ylabel(r'$\alpha_{BrC}$')
            plt.tight_layout()
            direc = cfg['working directory'] + f'plots/fitplots/'
            try:
                plt.savefig(direc+f'alpha_BrC.png', dpi = 300)
            except FileNotFoundError as fnfe:
                Path(direc).mkdir(parents=True, exist_ok=True)
                plt.savefig(direc+f'alpha_BrC.png', dpi = 300)
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)


    #----- Save plots for the swipe of alpha_bc if booked
    try:
        if cfg['alpha bc swipe']:
            print(f'---> Saving swipe plot in {cfg["working directory"]}' + f'plots/fitplots/' + '\n') 
            fig, ax = plt.subplots() 
            y = sw_alpha_BrC_set
            u_y = sw_alpha_BrC_stddev_set
            u_up, u_down = [a+b for (a,b) in zip(y, u_y)], [a-b for (a,b) in zip(y, u_y)]
            ax.plot(sw_alpha_BC_set, y, '-r', label='average')
            ax.fill_between(sw_alpha_BC_set, u_up, u_down, color='red', 
                    alpha=0.3, label=r'$1\sigma$ band')
            ax.grid(alpha=0.3)
            ax.set_ylabel(r'$\alpha_{BrC}$ averaged over all samples')
            ax.set_xlabel(r'$\alpha_{BC}$')
            ax.legend()
            plt.tight_layout()
            direc = cfg['working directory'] + f'plots/fitplots/'
            try:
                plt.savefig(direc+f'swipe.png', dpi = 300)
            except FileNotFoundError as fnfe:
                Path(direc).mkdir(parents=True, exist_ok=True)
                plt.savefig(direc+f'swipe.png', dpi = 300)
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)


    #----- Save plots for the optical apportionment if booked
    try:
        if cfg['plots']:
            print(f'---> Saving optical apportionment plots in {cfg["working directory"]}' + f'plots/appoplots/' + '\n') 
            # Lists to store values for plotting
            bc_ff_short, bc_wb_short, brc_short, names = [], [], [], []
            bc_ff_long, bc_wb_long, brc_long = [], [], []
            # Shortest and longest wavelength indices
            # (in case the lambda are not in ascending order)
            lambda_short = min(data[0].properties.wavelength)
            i_short = data[0].properties.wavelength.index(lambda_short)
            lambda_long = max(data[0].properties.wavelength)
            i_long = data[0].properties.wavelength.index(lambda_long)
            for sample in data:
                prp = sample.properties
                brc_short.append(prp.brc[i_short])
                bc_wb_short.append(prp.bc_wb[i_short])
                bc_ff_short.append(prp.bc_ff[i_short])
                bc_ff_long.append(prp.bc_ff[i_long])
                bc_wb_long.append(prp.bc_wb[i_long])
                brc_long.append(prp.brc[i_long])
                names.append(sample.name)
            fig1, ax1 = plt.subplots() # For short lambda
            fig2, ax2 = plt.subplots() # For long lambda
            ax1.plot(names, bc_ff_short, '-g', label=r'BC$_{FF}$'+ f'@ {lambda_short} nm')
            ax1.plot(names, bc_ff_short, '.g')
            ax1.plot(names, bc_wb_short, '-b', label=r'BC$_{WB}$'+ f'@ {lambda_short} nm')
            ax1.plot(names, bc_wb_short, '.b')
            ax1.plot(names, brc_short, '-r', label=r'BrC'+ f'@ {lambda_short} nm')
            ax1.plot(names, brc_short, '.r')
            ax1.set_xticklabels(names, rotation=75)
            ax1.grid(alpha=0.3)
            ax1.set_ylabel(r'Absorption coefficient, $b_{abs}$  ' + '[Mm'+ r'$^{-1}$' + ']') if prp.data_type == 'Babs' else plt.ylabel('100 ABS')
            ax1.legend()
            ax2.plot(names, bc_ff_long, '-g', label=r'BC$_{FF}$'+ f'@ {lambda_long} nm')
            ax2.plot(names, bc_ff_long, '.g')
            ax2.plot(names, bc_wb_long, '-b', label=r'BC$_{WB}$'+ f'@ {lambda_long} nm')
            ax2.plot(names, bc_wb_long, '.b')
            ax2.plot(names, brc_long, '-r', label=r'BrC'+ f'@ {lambda_long} nm')
            ax2.plot(names, brc_long, '.r')
            ax2.set_xticklabels(names, rotation=75)
            ax2.grid(alpha=0.3)
            ax2.set_ylabel(r'Absorption coefficient, $b_{abs}$  ' + '[Mm'+ r'$^{-1}$' + ']') if prp.data_type == 'Babs' else plt.ylabel('100 ABS')
            ax2.legend()
            fig1.tight_layout()
            fig2.tight_layout()
            direc = cfg['working directory'] + f'plots/appoplots/'
            try:
                fig1.savefig(direc+f'short_lambda.png', dpi = 300)
                fig2.savefig(direc+f'long_lambda.png', dpi = 300)
            except FileNotFoundError as fnfe:
                Path(direc).mkdir(parents=True, exist_ok=True)
                fig1.savefig(direc+f'short_lambda.png', dpi = 300)
                fig2.savefig(direc+f'long_lambda.png', dpi = 300)
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)


    #--- Save k1 and k2 regression plots if the fit is booked
    print('---> Saving mass apportionment plots...\n')
    if cfg['plots'] and do_fit:
        def line(x, m, q):
            return m * x + q
        x1 = np.linspace(min(k1_x), max(k1_x), 100)
        line1 = line(x1, k1, int1)
        plt.plot(k1_x, k1_y, '.r')
        plt.plot(x1, line1, '--k',)
        plt.ylabel(r'OC concentration   [$\mu$g/m$^3$]') 
        plt.xlabel(r'$b_{abs}^{BC,FF}$'+f'@{lambda_long} nm    [Mm$^{-1}$]')
        plt.grid(alpha=0.3)
        direc = cfg['working directory'] + f'plots/mappoplots/'
        try:
            plt.savefig(direc+'k1_fit.jpg', dpi = 300)
        except FileNotFoundError as fnfe:
            Path(direc).mkdir(parents=True, exist_ok=True)
            plt.savefig(direc+'k1_fit.jpg', dpi = 300)
        plt.close()
        x2 = np.linspace(min(k2_x), max(k2_x), 200)
        line2 = line(x2, k2, int2)
        plt.plot(k2_x, k2_y, '.r')
        plt.plot(x2, line2, '--k')
        plt.ylabel(r'OC - OC$_{FF}$ concentration   [$\mu$g/m$^3$]') 
        plt.xlabel(r'$b_{abs}^{BrC}$'+f'@{lambda_short} nm   [Mm$^{-1}$]')
        plt.grid(alpha=0.3)
        direc = cfg['working directory'] + f'plots/mappoplots/'
        plt.savefig(direc+'k2_fit.png', dpi = 300)
        plt.close()


    #----- Save time series plots for the mass apportionment if booked
    try:
        if cfg['plots']:
            print(f'---> Saving time series mass apportionment plots in {cfg["working directory"]}' + f'plots/mappoplots/' + '\n') 
            # Lists to store values for plotting
            oc_ff, oc_wb, oc_nc, ec_ff, ec_wb, names = [], [], [], [], [], []
            for sample in data:
                prp = sample.properties
                names.append(sample.name)
                oc_ff.append(prp.oc_ff)
                oc_wb.append(prp.oc_wb)
                oc_nc.append(prp.oc_nc)
                ec_ff.append(prp.ec_ff)
                ec_wb.append(prp.ec_wb)
            fig1, ax1 = plt.subplots() 
            ax1.plot(names, oc_ff, '-g', label=r'OC$_{FF}$')
            ax1.plot(names, oc_ff, '.g')
            ax1.plot(names, oc_wb, '-b', label=r'OC$_{WB}$')
            ax1.plot(names, oc_wb, '.b')
            ax1.plot(names, oc_nc, '-r', label=r'OC$_{NC}$')
            ax1.plot(names, oc_nc, '.r')
            ax1.plot(names, ec_ff, '-k', label=r'EC$_{FF}$')
            ax1.plot(names, ec_ff, '.k')
            ax1.plot(names, ec_wb, '-m', label=r'EC$_{WB}$')
            ax1.plot(names, ec_wb, '.m')
            ax1.set_xticklabels(names, rotation=75)
            ax1.grid(alpha=0.3)
            ax1.set_ylabel(r'Concentration ' + r'[$\mu$g cm$^{-3}$]')
            ax1.legend()
            fig1.tight_layout()
            direc = cfg['working directory'] + f'plots/mappoplots/'
            try:
                fig1.savefig(direc+f'time_series.png', dpi = 300)
            except FileNotFoundError as fnfe:
                Path(direc).mkdir(parents=True, exist_ok=True)
                fig1.savefig(direc+f'time_series.png', dpi = 300)
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)



    #####################################################################
    # WRITE LOG FILE
    #####################################################################

    #----- Get the time of the analysis end
    end_time = datetime.datetime.now()

    #----- Prepare the lines to write
    date_start_line = f'Start time:\t{start_time.day}/{start_time.month}/{start_time.year}, {start_time.hour}:{start_time.minute}:{start_time.second}\n'
    date_end_line = f'End time:\t{end_time.day}/{end_time.month}/{end_time.year}, {end_time.hour}:{end_time.minute}:{end_time.second}\n'
    input_file_line = 'Input file:\t' + cfg['input file'] + '\n'
    output_folder_line = 'Output folder:\t' + cfg['working directory'] + '\n'
    presets_line = 'Booked presets:\t' + str(cfg['presets']) + '\n'
    best_par_line = f'Best parameters: \n\talpha_BC = {best_alpha_BC} \n\talpha_FF = {best_alpha_FF} \n\talpha_WB = {best_alpha_WB}\n'
    saved_par_plots_line = f"Parameter optimization plots in:\t{cfg['working directory']}plots/preplots/\n" if cfg['plots'] and len(cfg['presets']) > 0 else f"Parameter optimization plots not saved\n"
    failed_fit_count_line = f'N° failed fits:\t{failed_fit_count}\n'
    failed_fit_line = f'Failed fits:\t{failed_fit}\n'
    # Get a list to do statistics on alpha brown
    list_for_abrc = [d.properties.alpha_brc for d in data]
    list_for_uabrc = [d.properties.u_alpha_brc for d in data]
    avg_alpha, stddev_alpha = average( list_for_abrc), stddev(list_for_abrc)
    alpha_mean_line = f"Average alpha_BrC:\t {round(avg_alpha, 7)}\n"
    alpha_stddev_line = f"Uncertainty (stdev) on alpha_BrC:\t {round(stddev_alpha, 7)}\n"
    saved_fit_plots_line = f"Fit plots in:\t{cfg['working directory']}plots/fitplots/\n" if cfg['plots'] else f"Fit plots not saved\n"
    saved_appo_plots_line = f"Optical apportionment plots in:\t{cfg['working directory']}plots/appoplots/\n" if cfg['plots'] else f"Optical apportionment plots not saved\n"
    k1_line = f"k1:\t {round(k1, 3)}\t\tR^2:\t{round(fit_1.rvalue ** 2 ,3)}\n" if do_fit else f"k1:\t {round(k1, 3)}\n"
    k2_line = f"k2:\t {round(k2, 3)}\t\tR^2:\t{round(fit_2.rvalue ** 2 ,3)}\n" if do_fit else f"k2:\t {round(k2, 3)}\n"
    saved_reg_plots_line = f"Linear regression plots in:\t{cfg['working directory']}plots/mappoplots/\n" if cfg['plots'] and do_fit else f"Linear regression plots not saved\n"
    try:
        with open(cfg['working directory'] + 'log.txt', 'w') as f:
            print(f'---> Writing log file in {cfg["working directory"]}' + '\n') 
            f.write('---------- GENERAL\n')
            f.write(date_start_line)
            f.write(date_end_line)
            f.write(input_file_line)
            f.write(output_folder_line)
            f.write('\n---------- PREPROCESSING\n')
            f.write(presets_line)
            f.write(best_par_line)
            f.write(saved_par_plots_line)
            f.write('\n---------- FIT\n')
            f.write(failed_fit_count_line)
            f.write(failed_fit_line)
            f.write(alpha_mean_line)
            f.write(alpha_stddev_line)
            f.write(saved_fit_plots_line)
            f.write('\n---------- OPTICAL APPORTIONMENT\n')
            f.write(saved_appo_plots_line)
            f.write('\n---------- MASS APPORTIONMENT\n')
            f.write(k1_line)
            f.write(k2_line)
            f.write(saved_reg_plots_line)
    except KeyError as ke:
        print(MISSING_KEYWORD, ke)
    except FileNotFoundError as fnfe:
        Path(cfg['working directory']).mkdir(parents=True, exist_ok=True)



    ##################################################################
    # SAVE PICKLE
    ##################################################################
    print(f"---> Saving .pkl data file for internal use in {cfg['working directory']}\n")
    try:
        with open(cfg['working directory'] + 'data.pkl', 'wb') as f:
            pickle.dump(data, f, -1)
    except FileNotFoundError as fnfe:
        Path(cfg['working directory']).mkdir(parents=True, exist_ok=True)


    print('*** DONE ***\n\n')
    










#----- Execution function
def cli_script():
    run(sys.argv[1])

import sys
from os.path import normpath as path
from shutil import copy2
try:
    import matplotlib.pyplot as plt
    import pandas as pd
    from scipy.interpolate import interp1d
    from scipy.optimize import curve_fit
    import numpy as np
    from scipy.special import exp10
except ModuleNotFoundError:
    print('Error loading required Python packages, this script requires: Numpy, Scipy, Pandas, and Matplotlib.\nSee tests.md for details.')
    sys.exit(1)
from test_units.Data import Data
import test_units.plot_settings


def fit_RC(dat_sim):
    coeffs = np.polyfit(dat_sim.x[1:], np.log(dat_sim.y[1:]), 1)
    RC_time, offset = -1 / coeffs[0], coeffs
    dat_fit = Data(dat_sim.x[1:], np.exp(np.polyval(coeffs, dat_sim.x[1:])))
    return RC_time, coeffs, dat_fit


def fit_exp(dat_sim, p_guess=None):
    def exp(x, Tau, f0, b):
        return f0*np.exp(-x/Tau) + b
    coeffs, pcov = curve_fit(exp, dat_sim.x, dat_sim.y, p0=p_guess)

    dat_fit = Data(dat_sim.x, exp(dat_sim.x, coeffs[0], coeffs[1], coeffs[2]))
    return coeffs, dat_fit


def fit_exp_no_offset(dat_sim, p_guess=None, x_offset=0):
    def exp(x, Tau, f0):
        return f0*np.exp(-(x-x_offset)/Tau) + 0
    coeffs, pcov = curve_fit(exp, dat_sim.x, dat_sim.y, p0=p_guess)

    dat_fit = Data(dat_sim.x, exp(dat_sim.x, coeffs[0], coeffs[1]))
    return coeffs, dat_fit


def calc_rmse(y_sim, y_test):
    """Calculates the root mean square error. Make sure the x corresponds exactly
    between the two arrays."""
    mse = (np.square(y_sim - y_test).mean())
    rmse = mse ** 0.5
    return rmse


def calc_selected_rmse(dat_sim, dat_test, log_x=False, log_y=False, x_min=None, x_max=None):

    if log_x == False:
        x_sim, x_test = dat_sim.x, dat_test.x
    else:
        x_sim, x_test = dat_sim.log_x, dat_test.log_x
        if x_min != None and x_max != None:
            x_max_log = np.log10(max([abs(x_min), abs(x_max)]))
            x_min_log = np.log10(min([abs(x_min), abs(x_max)]))
            x_min, x_max = x_min_log, x_max_log
        elif x_min != None:
            x_min = np.log10(abs(x_min))
        elif x_max != None:
            x_max = np.log10(abs(x_max))

    if log_y == False:
        y_sim, y_test = dat_sim.y, dat_test.y
    else:
        y_sim, y_test = dat_sim.log_y, dat_test.log_y

    y_interp = interp1d(x_test, y_test, fill_value="extrapolate")

    if x_min != None:
        idxs = np.nonzero(np.array(x_sim) > x_min)[0]
        x_sim = x_sim[idxs]
        y_sim = y_sim[idxs]

    if x_max != None:
        idxs = np.nonzero(np.array(x_sim) < x_max)[0]
        x_sim = x_sim[idxs]
        y_sim = y_sim[idxs]

    y_test_int = y_interp(x_sim)
    rmse = calc_rmse(y_sim, y_test_int)

    if log_x == True:
        x_sim = exp10(x_sim)

    if log_y == True:
        y_test_int = exp10(y_test_int)

    selected_pts = Data(x_sim, y_test_int)
    return rmse, selected_pts


def get_data_from_sim(simulation, test_idx, x='Vext', y='Jext', output_file='jv'):
    test_dev_par_path = 'test_{:n}/device_parameters.txt'.format(test_idx)
    test_trap_path = 'test_{:n}/traps.txt'.format(test_idx)
    try:
        copy2(path(test_trap_path), path(simulation.work_dir))
    except FileNotFoundError:
        pass
    copy2(path(test_dev_par_path), path(simulation.work_dir))
    if simulation.code_name.lower() == 'zimt':
        test_tvg_path = 'test_{:n}/tVG.txt'.format(test_idx)
        copy2(path(test_tvg_path), simulation.work_dir)
    simulation.run()
    data = simulation.return_out_dic()[output_file]
    sim_data = Data(data[x], data[y])

    return sim_data


def get_data_from_file(test_idx, x_test, y_test):
    test_data_path = 'test_{:n}/test{:n}.dat'
    test_data = pd.read_csv(path(test_data_path.format(
        test_idx, test_idx)), delim_whitespace=True)
    test_data = Data(test_data[x_test], test_data[y_test])
    return test_data


def plot_test_results(test_idx, dat_sim, ax_lab, dat_tests=[], log_x=False, log_y=False, dat_err=[]):
    figure_path = 'test_{:n}/test_{:n}_result.png'
    fig = plt.figure()
    plt.plot(dat_sim.x, dat_sim.y, 'o', ms=4, color='black',
             markerfacecolor='white', markeredgewidth=0.8, label='Simulation')

    for label, dat_test in dat_tests:
        plt.plot(dat_test.x, dat_test.y, label=label)

    for error, dat_section in dat_err:
        plt.plot(dat_section.x, dat_section.y, 'x',
                 label='RMSE={:.1e}'.format(error))

    if log_x:
        plt.xscale('log')
    if log_y:
        plt.yscale('log')
    plt.legend()
    plt.xlabel(ax_lab['x'])
    plt.ylabel(ax_lab['y'])
    plt.tight_layout()
    plt.savefig(path(figure_path.format(test_idx, test_idx)))
    print('Test {:n} finished'.format(test_idx))


def print_test_results(test_results):
    all_tests_passed = True
    for result in test_results:
        if result['passed'] == False:
            print('Test {:n} FAILED!'.format(result['test nr']))
            all_tests_passed = False
    if all_tests_passed:
        print('All tests passed successfully!')

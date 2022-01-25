import sys
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt 
except ModuleNotFoundError:
    print('Error loading required Python packages, this script requires: Numpy, Scipy, Pandas, and Matplotlib.\nSee tests.md for details.')
    sys.exit(1)
    
fig_size_cm = [10, 7]
mpl.rcParams['savefig.dpi'] = 600

font_weight = 'normal' #'bold'
fs_small = 6
fs_medium = 7
fs_large = 9

tick_padding = 3.4
line_thickness = 0.7# default = 1.1

def get_line_thickness():
    return line_thickness
# print(plt.rcParams.keys())
def cm2inch(value_cm):
    return value_cm / 2.54

#plt.rc('savefig', orientation='landscape', format='pdf')
plt.rc('figure', figsize=[cm2inch(fig_size_cm[0]), cm2inch(fig_size_cm[1])]) #a4 size
plt.rc('font', weight=font_weight, size=fs_medium)          # controls default text sizes

plt.rc('axes', linewidth=line_thickness, titlesize=fs_large)     # fontsize of the axes title
plt.rc('axes', labelweight=font_weight, labelsize=fs_medium)    # fontsize of the x and y labels
plt.rc('xtick', top=True, bottom=True, direction='in', labelsize=fs_small)    # fontsize of the tick labels
plt.rc('xtick.major', pad=tick_padding, size=line_thickness * 4, width=line_thickness)
plt.rc('xtick.minor', pad=tick_padding, size=line_thickness * 2.5, width=line_thickness)
plt.rc('ytick', left=True, right=True, direction='in', labelsize=fs_small)    # fontsize of the tick labels
plt.rc('ytick.major', pad=tick_padding, size=line_thickness * 4, width=line_thickness)
plt.rc('ytick.minor', pad=tick_padding, size=line_thickness * 2.5, width=line_thickness)
plt.rc('legend', fontsize=fs_small)    # legend fontsize
plt.rc('figure', titlesize=fs_large)  # fontsize of the figure title
plt.rc('lines', markersize=4*line_thickness, linewidth=line_thickness*1.5)

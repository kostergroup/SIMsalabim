"""
This class runs simsalabim, taking command line arguments as dictionary.
It then runs simsalabim and saves the output in dictionary format without
creating output files in order to save my ssd (to my knowledge).

Created by: Marten Koopmans
"""
import subprocess
import copy
import os
import sys

try:
    import pandas as pd
except ModuleNotFoundError:
    print('Error loading required Python packages, this script requires: Numpy, Scipy, Pandas, and Matplotlib.\nSee tests.md for details.')
    sys.exit(1)

class SIMsalabim():
    """Runs simsalabim from input dictionary, gathers output and returns
    the output in a dictionary."""
    def __init__(self, inp_dict=None, file_nr=0, include_scpars=True, include_tj=True, include_var=True, include_jv=True, code_name='SimSS', work_dir='../SimSS', output_sub_folder='.temp', show_term_output=False):
        self.input_dic = inp_dict
        self.code_name = code_name
        self.set_work_dir(work_dir, output_sub_folder)
        self.work_dir = os.path.normpath(work_dir)
        self.windows = not (os.name == 'posix')
        self.show_term_output = show_term_output
        self.output_file_ext = str(file_nr) + '_file.dat'
        self.output_dic = {}
        self.out_fil_nam = {}

        if code_name == 'SimSS' or code_name == 'SIMsalabim':
            if include_jv == True:
                self.output_dic['jv'] = None
                self.out_fil_nam['-JVFile'] = output_sub_folder +'/JV' + str(file_nr) + '_file.dat'

            if include_scpars == True:
                self.output_dic['scpars'] = None
                self.out_fil_nam['-scParsFile'] = output_sub_folder + '/SCpars' + str(file_nr) + '_file.dat'

        if code_name == 'ZimT':
            if include_tj == True:
                self.output_dic['tj'] = None
                self.out_fil_nam['-tJFile'] = output_sub_folder + '/tj' + str(file_nr) + '_file.dat'
                
        if include_var == True:
            self.output_dic['var'] = None
            self.out_fil_nam['-varFile'] = output_sub_folder +'/var' + str(file_nr) + '_file.dat'

        self.make_command_list()


    def set_input_dic(self, input_dictionary):
        """feed an input dictionary to the simulation."""
        self.input_dic = input_dictionary


    def make_command_list(self):
        """Makes the command list to run SIMsalabim."""
        if self.windows:
            executable = self.code_name.lower() + '.exe'
        else:
            executable = './' + self.code_name.lower()
        command_list = [executable]
        if self.input_dic:
            for input_par in self.input_dic:
                command_list.append('-' + input_par)
                command_list.append(str(self.input_dic[input_par]))

        # add the object specefik output file names to the cmd list.
        for command in self.out_fil_nam:
            command_list.append(command)
            command_list.append(os.path.normpath(self.out_fil_nam[command]))

        self.cmd_list = command_list


    def compile_code(self):
        if self.show_term_output == True:
            output_direct = None
        else:
            output_direct = subprocess.DEVNULL
        try:
            subprocess.check_call(['fpc', self.code_name.lower()], encoding='utf8', stdout=output_direct, cwd=self.work_dir, shell=self.windows)
        except subprocess.CalledProcessError:
            print(self.work_dir)
            print('Code \'{}\' failed to compile!'.format(self.code_name))
            sys.exit(1)


        
    def run(self):
        """runs simulation..."""
        output = None
        try:
            output = subprocess.run(self.cmd_list, encoding='utf8', capture_output=True, cwd=self.work_dir, shell=self.windows)
            return(output)
        except subprocess.CalledProcessError as err:
            print('{} exited with non-zero status.'.format(self.code_name))

    def set_work_dir(self, working_directory, output_subdir):
        work_path = os.path.dirname(working_directory)
        if not os.path.exists(work_path):
            raise EOFError('SIMsalabim path does not exist: \'{}\''.format(working_directory))

        self.work_dir = os.path.normpath(working_directory)
        sub_path = os.path.join(working_directory, output_subdir)
        if not os.path.exists(sub_path):
            os.makedirs(sub_path)


    def return_out_dic(self):
        """Return the output (files) of SIMsalabim in dictionary form."""
        for file in self.out_fil_nam:
            try:
                file_path = os.path.join(self.work_dir, self.out_fil_nam[file])
                
                with open(file_path) as out_fil:
                    dic_entry = file[1:].split('File')[0].lower()

                    data = pd.read_csv(out_fil, delim_whitespace=True)
                    if data.empty:
                        data = None
                    self.output_dic[dic_entry] = data

            except FileNotFoundError:
                pass
        return copy.deepcopy(self.output_dic)     


    def run_return_output(self):
        """This is the function to run a SIMsalabim simulation and return the output to outside the object."""
        self.make_command_list()
        self.run_simsalabim()
        return self.return_out_dic()

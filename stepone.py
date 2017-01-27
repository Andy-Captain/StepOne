import tkinter as tk
import tkinter.filedialog
import csv
import os
import re
import pandas as pd
import numpy as np
import shutil
import zipfile as zip

__author__ = 'Lukas Jaworski'
__version__ = '0.6.7'
# TODO: Add a hell of a lot more documentation. Also should I mirror some functions in R?
# TODO: Make a version with a lot less hand holding and assume the user knows what they're doing (remove setters and
#       getters) <- Most of this code should be identical in both versions, one will just be intended for less
#       experienced to novice programmers. This other version should run much faster, but since I am not actually a CS
#       some of the finer points of optimization, especially inter language ones elude me.


class StepOne(object):
    """
    A class that works with the output files of the StepOne/StepOnePlus(TM) Real-Time PCR system.
    It works with the scientific computing environment in Python and is helpful if you have Anaconda installed as there
    are many libraries that are required to run this class.
    """
    # TODO: Add method that pulls from multicomponent data and the meltcurve to generate your own meltcurves.
    #       This can be done because the cycles 0-39 hold the amplification data and the rest holds the melt curve data
    #       The meltcurve_results holds the temperatures while the multicompenent_data holds the raw floresence values.
    #       The raw florescence values have no scale or units and are referenced to the ROX. <- optional.
    # TODO: Add own exponentional alogrithm. <-This is very optional.
    def __init__(self, working_directory=None):
        """
        Initializes a bunch of internal variables and gets a working directory in which its functions will be run.

        :param working_directory: String containing a valid directory.
        """
        self.__directory_attempts = 3
        self.__labels = ('Sample', 'Detector', 'Cq', 'Contaminated')
        self.__directory = ''
        self.__file_list = []
        self.__dataframe = pd.DataFrame()
        self.__sample_ids = {}
        self.__temp = 'temp'
        self.__sample_id_filename = 'sampleids.csv'
        self.__controls = ['H2O']
        self.__lastfile = ''
        if working_directory is None:
            self.directory_dialog(working_directory)
        else:
            self.dir = working_directory

    @property
    def dir(self):
        """
        Returns the currently stored working directory.

        :return: String containing the stored working directory.
        """
        return self.__directory

    @dir.setter
    def dir(self, directory):
        """
        Takes a string and uses it as the path to the current working directory.

        :param directory: String containing a valid directory.
        :return: Internally set working directory to input path.
        """
        try:
            assert isinstance(directory, str)
        except AssertionError:
            raise TypeError('Arguments must be strings.')
        if not os.path.isdir(directory):
            raise NotADirectoryError('Not a valid directory!')
        else:
            self.__directory = directory

    @dir.deleter
    def dir(self):
        """
        Sets the directory to a slash (root directory)

        :return: Internally sets directory path to root.
        """
        self.__directory = os.path.abspath(os.sep)

    @property
    def df(self):
        """
        Returns the current dataframe stored within the StepOne object.

        :return: Internally stored dataframe.
        """
        return self.__dataframe

    @df.setter
    def df(self, dataframe):
        """
        When setting the dataframe manually this method tries to check that the column headings of the incoming
        dataframe match the internally set column titles and warn you if there is a any mismatch, it will not tell you
        which is mismatched and how many, that is up to the user to check, conversely the internally stored labels can
        be changed using the set_df_labels method. Conversely the dataframe is still a pandas object and all panads
        functions will work on it, the resulting dataframe may need to be imported.

        :param dataframe: Import a Pandas Dataframe into the StepOne object.
        :return: Internally store input Dataframe.
        """
        try:
            assert isinstance(dataframe, pd.DataFrame)
        except:
            raise TypeError('Need a dataframe')
        if dataframe.empty:
            raise Warning('Trying to import an empty dataframe, if you want to delete the current dataframe use del')
        else:
            self.__dataframe = dataframe

    @df.deleter
    def df(self):
        """
        The deleter simple sets the dataframe to an empty dataframe.

        :return: Set internal Dataframe to the empty Dataframe.
        """
        self.__dataframe = pd.DataFrame()

    @property
    def controls(self):
        """
        Returns the current list of controls.

        :return: List of current controls.
        """
        return self.__controls

    @controls.setter
    def controls(self, controls):
        """
        Manually sets the internal list of controls. If a single value is used then a simple string can be entered.
        The list of controls is used by the scrub controls method to remove the controls from the dataframe as there is
        no interest in analyzing them.

        :param controls: A list-like object containing the IDs(Names) of the controls in the plate
        :return: Store the list of controls internally.
        """
        if isinstance(controls, str):
            controls = [controls]
        controls = list(controls)
        try:
            for control in controls:
                assert isinstance(control, str)
        except AssertionError:
            raise ValueError('Need a list of strings of the filenames, even fore just one file')
        else:
            controls = [control.strip() for control in controls]
            self.__controls = controls

    @controls.deleter
    def controls(self):
        """
        The deleter sets the controls list to the default of a single string list of 'H2O'.

        :return: Internally set the controls to ['H2O'].
        """
        self.__controls = ['H2O']

    # TODO Change the implementation of the labels function and the methods it depends on. Maybe remove it altogether
    #  and pass it as a default parameter
    @property
    def df_labels(self):
        """
        The labels used in the Pandas' dataframe object which is stored within this object.

        :return: Dataframe labels
        """
        return self.__labels

    @df_labels.setter
    def df_labels(self, parameters):
        """
        Set the dataframe labels directly by inputing a tuple-like (gets converted to tuple) object.
        Order MATTERS and should be ('Sample', 'Detector', 'Cq', 'Contaminated'). Place an empty string if you don't
        plan on using a label.

        :param parameters: Tuple-like object containing the four labels for ('Sample', 'Detector', 'Cq', 'Contaminated')
        :return: Internally stores the labels that will be used in the dataframe for these critical parameters
        """
        parameters = tuple(parameters)
        if len(parameters) != 4:
            raise ValueError('Incorrect number of parameters, needs 4 ("Detector", "ID", "Cq", "Contaminated")')
        for value in parameters:
            try:
                assert isinstance(value, str)
            except AssertionError:
                raise TypeError('Tuple values must be strings.')
        self.__labels = parameters

    @df_labels.deleter
    def df_labels(self):
        """
        The deleter resets the major data labels to their default values.

        :return: Internally stored labels are set to ('Sample', 'Detector', 'Cq', 'Contaminated').
        """
        self.__labels = ('Sample', 'Detector', 'Cq', 'Contaminated')

    @property
    def eds_files(self):
        """
        Displays the internally stored .eds file list. This list is the list of files that this object works on.

        :return: Internally stored .eds file list.
        """
        return self.__file_list

    @eds_files.setter
    def eds_files(self, file_list):
        """
        Takes a list like object and checks each item for the proper file ending then stores the list internally.
        This is used to set the list of files that will be analyzed manually. There is a method that scans the working
        directory and gets a list of all of the files there. This method will also replace the file list contained here,
        so caustion is advised.

        :param file_list: List-like object containing strings or a string. All strings must end in .eds
        :return: Internally stores the list of files inside the object
        """
        if isinstance(file_list, str):
            file_list = [file_list]
        file_list = list(file_list)
        try:
            for file in file_list:
                assert isinstance(file, str)
        except AssertionError:
            raise ValueError('Need a list of strings of the filenames, even fore just one file')
        else:
            file_list = [file.strip() for file in file_list]
            for file in file_list:
                if not file.endswith('.eds'):
                    raise ReferenceError('Only .eds files are allowed')
            self.__file_list = file_list

    @eds_files.deleter
    def eds_files(self):
        """
        The deleter sets the internal file list variable to the empty list.

        :return: Internal file list becomes the empty list.
        """
        self.__file_list = []

    @property
    def __last_file(self):
        """
        Returns the last file that was processed by the object.

        :return: String containing the last file that this object processed.
        """
        return self.__lastfile

    @__last_file.setter
    def __last_file(self, filename):
        """
        Stores the last file processed. This is for checking file endings and for error logging.

        :param filename: String containing a filename with a .eds ending.
        :return: Internally stores the name in a lastfile variable.
        """
        filename = str(filename)
        filename = filename.strip()
        if not filename.endswith('.eds'):
            raise ReferenceError('Only .eds files are allowed')
        else:
            self.__lastfile = filename

    @__last_file.deleter
    def __last_file(self):
        """
        Gets rid of the name of the last file that was processed and sets it to the empty string.

        :return: sets the internal lastfile variable to the empty string.
        """
        self.__lastfile = ''

    @property
    def sample_ids_file(self):
        """
        Returns the name of the file that is used convert from sample ids to full experimental parameters.
        This is done because the StepOne has limited space for sample ids which would not allow all experimental
        groupings to be displayed. This is the file used to fill in that information. This file is completely optional.
        The program assumes this file is in the working directory.

        :return: Internal name of sample id file.
        """
        return self.__sample_id_filename

    @sample_ids_file.setter
    def sample_ids_file(self, filename):
        """
        Sets the name of the file used to convert sample ids to experimental parameters.

        :param filename: String containing the filename of the sample ids to parameters file.
        :return: Internally store filename as a string.
        """
        try:
            assert isinstance(filename, str)
        except AssertionError:
            raise TypeError('Needs a string filename.')
        filename = filename.strip()
        if not (filename.endswith('.txt') or filename.endswith('.csv')):
            raise ReferenceError('Need a txt or csv file that can be read by the csv library')
        self.__sample_id_filename = filename

    @sample_ids_file.deleter
    def sample_ids_file(self):
        """
        Sets the sample id filename to the default of 'sampleids.csv'

        :return: Set internal sample id filename variable to 'sampleids.csv
        """
        self.__sample_id_filename = 'sampleids.csv'

    # TODO: Delete this if individual textfiles can be extracted.
    @property
    def temp(self):
        """
        Returns the name of the subdirectory where the extracted files are placed.

        :return: String containing subdirectory name.
        """
        return self.__temp

    @temp.setter
    def temp(self, subdirectory):
        """
        Set the name of the subdirectory in which the extracted files will be placed.

        :param subdirectory: String with name of subdirectory, no preceding slash.
        :return: Internally stores string containing subdirectory name.
        """
        try:
            assert isinstance(subdirectory, str)
        except AssertionError:
            raise TypeError('Needs a string subdirectory path (no preceding slash).')
        subdirectory = subdirectory.strip().lstrip(os.sep)
        self.__temp = subdirectory

    @temp.deleter
    def temp(self):
        """
        Sets the subdirectory name to the default of 'temp'

        :return: Internallys sets subdirectory folder name to 'temp'
        """
        self.__temp = 'temp'

    @property
    def sample_ids(self):
        """
        A dictionary containing the sample ids as the keys and the experimental parameters as the values. These values
        should later be used in filtering the samples for statistical analysis.

        :return: Dictionary containing the sample ids as the keys and the experimental parameters as the values.
        """
        return self.__sample_ids

    @sample_ids.setter
    def sample_ids(self, ids):
        """
        Uses this dictionary to add sample parameters using the simplified id numbers/names/values used to in creating
        the runfile and adding all of the experimental parameters of interest, uses the builtin

        :param ids: Dictionary-like with ID(key):sample parameters(value)<-should be string with sample parameters being
                                                                           separated by a consistent separator.
        :return: Internally store sample id dictionary.
        """
        self.__sample_ids = dict(ids)

    @sample_ids.deleter
    def sample_ids(self):
        """

        :return: Internally sets the sample id dictionary to the empty dictionary.
        """
        self.__sample_ids = {}

    def set_df_labels(self, detector=None, identifier=None, cq=None, contamination=None):
        """
        Use this function to set the dataframe labels, in a parameter-wise way. detector is the column label for the
        genes used in the study. identifier is column label for sample names used in the study. cq is the Cq/Ct number.
        contamination is is the label for if the sample is found to be contaminated using the flag_contamination
        method. DEFAULT LABELS:('Sample', 'Detector', 'Cq', 'Contaminated')

        :param detector: String containing the dataframe label for the Detectors/Genes of interest.
        :param identifier: String containing the dataframe label for the sample ids/names of the experiment.
        :param cq: String containing the dataframe label for the Cq/Ct values for the experiment.
        :param contamination: String containing the dataframe label indicating whether the sample was contaminated.
        :return: Internally set the dataframe labels.
        """
        if detector is not None:
            try:
                assert isinstance(detector, str)
            except AssertionError:
                raise TypeError('Arguments must be strings.')
            detector = detector
        else:
            detector = self.df_labels[1]
        if identifier is not None:
            try:
                assert isinstance(identifier, str)
            except AssertionError:
                raise TypeError('Arguments must be strings.')
            identifier = identifier
        else:
            identifier = self.df_labels[0]
        if cq is not None:
            try:
                assert isinstance(cq, str)
            except AssertionError:
                raise TypeError('Arguments must be strings.')
            cq = cq
        else:
            cq = self.df_labels[2]
        if contamination is not None:
            try:
                assert isinstance(contamination, str)
            except AssertionError:
                raise TypeError('Arguments must be strings.')
            contamination = contamination
        else:
            contamination = self.df_labels[3]
        self.df_labels = (detector, identifier, cq, contamination)

    def directory_dialog(self, directory=None):
        """
        Opens a file dialog where the user selects the directory they want to work from or conversely the working
        directory can be set manually using a string containing the directory. The directory can also be manually set by
        simply setting the <StepOne>.dir to the desired path (i.e. myStepOne.dir = 'C:/Users/Guest/Python')

        :param directory: String containing a valid directory.
        :return: Internally sets path of working directory.
        """
        self.__directory_attempts -= 1
        if not self.__directory_attempts:
            raise RecursionError('Too many directory attempts, crashing program...')
        if directory is None:
            root = tk.Tk()
            self.dir = tk.filedialog.askdirectory(parent=root, title='Please select a directory')
            root.destroy()
        else:
            try:
                self.dir = directory
            except (TypeError, NotADirectoryError):
                print('directory lookup failed, need a valid string representing a directory.')
                self.directory_dialog(None)
        self.__directory_attempts = 3

    def get_sample_ids(self, sample_id_filename=None):
        """
        Makes a dictionary object containing the sample IDs as the keys and the sample parameters as values. The file in
        question needs to be a csv with the sample IDs in the first column and the description parameters in the second.

        :param sample_id_filename: filename (as a string) DEFAULT: 'sampleids.csv'
        :return: Internally set filename.
        """
        if sample_id_filename is None:
            sample_id_filename = self.sample_ids_file
        else:
            sample_id_filename = sample_id_filename
        filepath = self.dir + '/' + sample_id_filename
        sample_id_dict = {}
        with open(filepath, 'r') as sample_id_file_object:
            sample_id_file_reader = csv.reader(sample_id_file_object)
            for line in sample_id_file_reader:
                sample_id_dict[line[0]] = line[1]
            sample_id_file_object.close()
        self.sample_ids = sample_id_dict

    def get_eds_files(self, directory=None):
        """
        Makes a list of all of the files with the .eds ending in the current directory

        :param directory: String containing a directory.
        :return: Internally store a list of strings containing the names of all the .eds files in directory.
                 It uses the internally stored working directory if none is provided.
        """
        if directory is None:
            directory = self.dir
        else:
            directory = directory
        file_list = []
        filename_regex = re.compile(r'(.*)(\.eds)$')
        for filename in os.listdir(directory):
            regexed_filename_parts = filename_regex.search(filename)
            if not (regexed_filename_parts is None):
                file_list.append(filename)
        self.eds_files = file_list

    @staticmethod
    def flag_contamination(dataframe, neg_cont='H2O', gene_list=None, thresh=37.0, diff=10.0, labels=('Sample',
                                                                                                      'Detector',
                                                                                                      'Cq',
                                                                                                      'Contaminated')):
        """
        This method is used to flag contamination in the PCR by flagging the samples in which the negative control
        showed amplification. This method needs to be called ONCE pcr plate at a time (one plate per input dataframe)
        otherwise controls from other plates may faux contaminate the current data.

        :param dataframe: A Pandas dataframe
        :param neg_cont: String containing the sample name of the negative control.
        :param gene_list: List-like object which contains the list of genes which you want to flag, if set to None it
                          will simply pull all the unique names found in the dataframe's Detector column.
        :param thresh: Float containing the threshold value (Cq value) below which the negative control is considered
                       contaminated.
        :param diff: Float which is closest (in terms of Cq value) a non negative control sample can get to a negative
                     control which is above thresh but below 40
        :param labels: Quartet (four string) tuple which holds the dataframe column labels for the (Detector/Gene,
                       Sample Name, Cq, is contaminated) respectively, these need to be set correctly or this method
                       won't run.
        :return: Dataframe with a column indicating whether the sample was contaminated.
        """
        try:
            assert isinstance(dataframe, pd.DataFrame)
        except AssertionError:
            raise TypeError('dataframe needs to be a pandas DataFrame object!')
        if gene_list is None:
            gene_list = set(dataframe.index.get_level_values(1))
        else:
            gene_list = gene_list
            if isinstance(gene_list, str):
                gene_list = [gene_list]
        gene_list = list(gene_list)
        try:
            assert isinstance(neg_cont, str)
        except AssertionError:
            raise TypeError('neg_cont needs to be a string indicating the ids of the negative control wells.')
        try:
            assert isinstance(thresh, float)
        except AssertionError:
            raise TypeError('thresh needs to be a float.')
        try:
            assert isinstance(diff, float)
        except AssertionError:
            raise TypeError('diff needs to be a float.')
        dataframe.sort_index(inplace=True)
        idx = pd.IndexSlice
        neg_cont_only = dataframe.loc[idx[neg_cont, :], :]
        no_neg_cont = dataframe.drop(neg_cont, level=0)
        if neg_cont_only.empty:
            raise ValueError('No negative controls found!')
        if no_neg_cont.empty:
            raise ValueError('No samples found!')
        for gene in gene_list:
            cq_value = list(neg_cont_only.loc[idx[:, gene], labels[2]])
            if cq_value:
                thresh_breached = cq_value[0] - thresh < 0.0001
                diff_thresh = (cq_value[0] - no_neg_cont.loc[idx[:, gene], labels[2]].max(axis=0)) < diff
                undetected = 40.0 - cq_value[0] < 0.0001
                if thresh_breached or (diff_thresh and not undetected):
                    dataframe.loc[idx[:, gene], labels[3]] = True
        dataframe.loc[~dataframe[labels[3]].isin([True]), labels[3]] = False
        return dataframe

    @staticmethod
    def scrub_contamination(dataframe, contam='Contaminated'):
        """

        :param dataframe:
        :param contam:
        :return:
        """
        try:
            assert isinstance(dataframe, pd.DataFrame)
        except AssertionError:
            raise TypeError('dataframe needs to be a pandas DataFrame object!')
        clean = dataframe.loc[dataframe[contam] == False]
        dirty = dataframe.loc[dataframe[contam] == True]
        return clean, dirty

    @staticmethod
    def scrub_controls(dataframe, controls='H2O'):
        """

        :param dataframe:
        :param controls:
        :return:
        """
        idx = pd.IndexSlice
        dirty = pd.DataFrame()
        try:
            assert isinstance(dataframe, pd.DataFrame)
        except AssertionError:
            raise TypeError('dataframe needs to be a pandas DataFrame object!')
        if isinstance(controls, str):
            controls = [controls]
        controls = list(controls)
        for cont in controls:
            try:
                assert isinstance(cont, str)
            except AssertionError:
                raise TypeError('Controls must be strings!')
            else:
                dirty = dirty.append(dataframe.sort_index().loc[idx[controls, :], :])
                dataframe = dataframe.drop(controls, level=0)
        return dataframe, dirty

    def extract_ct_values(self, filename, directory=None, dh='append', f_contam=True,
                          contam_params=('H2O', 37.0, 10.0)):
        """

        :param filename:
        :param directory:
        :param dh:
        :param f_contam:
        :param contam_params:
        :return:
        """

        temp_df = pd.DataFrame()
        dh_options = ('append', 'return', 'replace')
        try:
            assert isinstance(dh, str)
        except AssertionError:
            raise TypeError('Need a string keyword.')
        if dh not in dh_options:
            raise ValueError('Not a valid dh selection, options are: "append", "replace", "return"')
        resultfile_subdirectory = '/apldbio/sds/analysis_result.txt'
        if directory is None:
            directory = self.dir
        else:
            directory = directory
        try:
            assert isinstance(directory, str)
        except AssertionError:
            raise TypeError('Need a string to a directory.')
        if not os.path.isdir(directory):
            raise NotADirectoryError('%s is not a valid directory' % directory)
        assert isinstance(f_contam, bool)
        if f_contam is True:
            cp1 = str(contam_params[0])
            cp2 = float(contam_params[1])
            cp3 = float(contam_params[2])
            contam_params = (cp1, cp2, cp3)
        try:
            with zip.ZipFile(directory + '/' + filename) as edsfile:
                with edsfile.open('apldbio/sds/analysis_result.txt') as Workfile:
                    for each_line in Workfile:
                        each_line_parts = each_line.decode('utf-8').split('\t')
                        try:
                            int(each_line_parts[0])
                        except ValueError:
                            continue
                        else:
                            if 0 <= int(each_line_parts[0]) <= 383:
                                data = [[each_line_parts[1].strip().replace(' ', '_'),
                                         each_line_parts[2].strip().replace(' ', '_'),
                                         float(each_line_parts[4])]]
                                temp_df = temp_df.append(pd.DataFrame(data, columns=[self.df_labels[0],
                                                                                     self.df_labels[1],
                                                                                     self.df_labels[2]]),
                                                         ignore_index=True)
                    Workfile.close()
                edsfile.close()
            temp_df.set_index(list(self.df_labels[0:2]), inplace=True)
            if temp_df.empty:
                raise Exception('No data found! Choose a different file...')
            if f_contam:
                temp_df = self.flag_contamination(temp_df, neg_cont=contam_params[0], thresh=contam_params[1],
                                                  diff=contam_params[2], labels=self.df_labels)
            else:
                self.__last_file = filename
            if dh is 'append':
                self.merge_frames(temp_df)
            elif dh is 'return':
                return temp_df
            elif dh is 'replace':
                self.df = temp_df
            else:
                raise Exception('Not an acceptable selection')
        except FileNotFoundError:
            print(filename)
            print("This file doesn't exist!")

    def add_melt_temps(self, filename, directory=None, dh='append', f_mult_tm=True):
        """

        :param filename:
        :param directory:
        :param dh:
        :param f_mult_tm:
        :return:
        """

        temp_df = pd.DataFrame()
        dh_options = ('append', 'return', 'replace')
        try:
            assert isinstance(dh, str)
        except AssertionError:
            raise TypeError('Need a string keyword.')
        if dh not in dh_options:
            raise ValueError('Not a valid dh selection, options are: "append", "replace", "return"')
        if directory is None:
            directory = self.dir
        else:
            directory = directory
        try:
            assert isinstance(directory, str)
        except AssertionError:
            raise TypeError('Need a string to a directory.')
        if not os.path.isdir(directory):
            raise NotADirectoryError('%s is not a valid directory' % directory)
        try:
            assert isinstance(f_mult_tm, bool)
        except AssertionError:
            raise TypeError('f_mult_tm must be a bool.')
        try:
            with zip.ZipFile(directory + '/' + filename) as edsfile:
                with edsfile.open('apldbio/sds/meltcuve_result.txt') as Workfile:
                    for each_line in Workfile:
                        each_line_parts = each_line.decode('utf-8').split('\t')
                        try:
                            int(each_line_parts[0])
                        except ValueError:
                            continue
                        else:
                            try:
                                float(each_line_parts[4])
                            except ValueError:
                                melt_temps = each_line_parts[4].split(',')
                                melt_temps = [float(temp.strip()) for temp in melt_temps if temp.strip()]
                            else:
                                if each_line_parts[4] is '\n':
                                    continue
                                melt_temps = [float(each_line_parts[4])]
                            if 0 <= int(each_line_parts[0]) <= 383:
                                data_id = [[each_line_parts[1].strip().replace(' ', '_'),
                                            each_line_parts[2].strip().replace(' ', '_')]]

                                temp_df = temp_df.append(pd.DataFrame(data_id, columns=[self.df_labels[0],
                                                                                        self.df_labels[1]]),
                                                         ignore_index=True)
                                for index, temp in enumerate(melt_temps, start=1):
                                    melt_header = 'Melt_Temp_' + str(index)
                                    if melt_header not in temp_df.columns:
                                        temp_df[melt_header] = np.nan
                                    temp_df.loc[(temp_df[self.df_labels[0]] == data_id[0][0]) &
                                                (temp_df[self.df_labels[1]] == data_id[0][1]), melt_header] = temp
                    Workfile.close()
                edsfile.close()
            temp_df.set_index(list(self.df_labels[0:2]), inplace=True)
            if temp_df.empty:
                raise Exception('No data found! Choose a different file...')
            if f_mult_tm:
                pass
            else:
                self.__last_file = filename
            if dh is 'append':
                self.merge_frames(temp_df)
            elif dh is 'return':
                return temp_df
            elif dh is 'replace':
                self.df = temp_df
            else:
                raise Exception('Not an acceptable selection')
        except FileNotFoundError:
            print(filename)
            print('This file does not have any data, going to the next one...')

    def merge_frames(self, dataframe, controls=None):
        """
        Merges dataframes based on ID and Detector.

        :param dataframe: The dataframe that will be merged with the internally stored dataframe.
        :param controls:   A list like object providing a list of the ID of the control wells.
        :return: Updates the internal dataframe with the new values.
        """
        if controls is None:
            controls = set(self.controls)
        else:
            controls = set(controls)
        try:
            assert isinstance(dataframe, pd.DataFrame)
        except AssertionError:
            raise TypeError('Need a pandas DataFrame.')
        if not self.df.empty:
            if (set(dataframe.index.get_level_values(0)) - controls)\
                    .intersection(set(self.df.index.get_level_values(0)) - controls):
                if set(self.df.columns.values).intersection(dataframe.columns.values):
                    self.df.update(dataframe)
                else:
                    self.df = pd.concat([self.df, dataframe], axis=1)
            else:
                self.df = pd.concat([dataframe, self.df])
        else:
            self.df = pd.concat([dataframe, self.df])

    def _remove_tempfiles(self):
        """

        :return:
        """
        iterations = 0
        while True:  # I had to add this because the virus scanner kept it from being removed otherwise without an error
            try:
                shutil.rmtree(self.dir + '/' + self.temp)
                break
            except WindowsError:
                if iterations > 10000:  # any longer than this, then the problem is more than just the virus scanner
                    raise
                iterations += 1

    @staticmethod
    def unzip_file(directory, filename, unzip_folder):
        """
        Takes the eds file[filename], with the extension included, in the directory[directory] and unzips the
        contents of the file into a folder[unzip_folder] within the same directory.

        :param directory:
        :param filename:
        :param unzip_folder:
        :return:
        """
        assert isinstance(filename, str)
        if not filename.endswith('.eds'):
            raise Exception('File ending is invalid, need an eds file.')
        assert isinstance(directory, str)
        if not os.path.isdir(directory):
            raise Exception('Not a valid directory')
        assert isinstance(unzip_folder, str)
        zip_ref = zip.ZipFile(directory + '/' + filename, 'r')  # read when doing actual files
        zip_ref.extractall(directory + '/' + unzip_folder)
        zip_ref.close()

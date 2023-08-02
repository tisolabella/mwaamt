# Imports
import sys
import json
import pandas as pd
import numpy as np
from m2too.property import Sample
from pprint import pprint

# Constants
DEFAULT_ABS_UNC = .1

# Class definition
class DataPrepper():
    """
    A class that reads data from a variety of formats and 
    prepares it for usage within the m2too package.
    The standard way to usage of the class is indirectly through
    the dataprep.prepare function that wraps the DataPrepper class
    functionalities.
    """
    def __init__(self, config=None):
        """infile can be None if the data is read directly from the terminal
        (NOT YET IMPLEMENTED)
        """
        self.config = config
        self.file = self.config['input file']

    def get_data(self):
        self.get_data_from_file()

    def get_data_from_stdin(self):
        """
        To implement once I figure out a way to turn input strings
        into list (simple overcast doesn't work, it turns the string
        into a list of characters, e.g. [1,2]= ['[','1',',','2',']']
        """
        pass

    def get_data_from_file(self):
        if self.file[-3:] == 'csv':
            self.open_csv()
        elif self.file[-4:] == 'json':
            self.open_json()

    def open_csv(self):
        '''
        The format of the .csv input file has to be
        Name,Wavelength,ABS,EC,OC              (<--- header row)
        N---,W [nm]----,A--,ec,oc [ug/cm3]     (<---data rows)
        
        Potentially implement a routine in which it is possible
        to read an additional column, named errABS, for the
        uncertainty on the ABS value.
        '''
        rawdata = pd.read_csv(self.file)

        # Get the three basic analysis components
        try:
            names = [x for x in rawdata['Name']]
            wavelengths = [float(x) for x in rawdata['Wavelength']]
            ABS = [float(x) for x in rawdata['ABS']]
        except KeyError as ke:
            print(f"KeyError: {ke} \nPossible solution: the input file has to be formatted  as follows: a .csv file with the following template: \n  Name \tWavelength \tABS \tuABS \tEC \tOC \n  --- \t --- \t --- \t --- \t ---\t --- \n  where the Wavelength is in nm, and EC and OC are in ug/cm3. \n  uABS is optional, as are EC and OC.")
        
        # Check if ABS uncertainty is provided
        uncertainty = True if 'uABS' in rawdata.keys() else False
        if uncertainty:
            u_ABS = [float(x) for x in rawdata['uABS']]

        # Check if EC and OC are provided
        mass_appo = True if 'EC' in rawdata.keys() and 'OC' in rawdata.keys() else False
        # Check if a mass apportionment is required 
        try:
            if self.config['mass_appo'] and not mass_appo:
                print("WARNING: the mass_appo flag has been set to True but the EC and OC values were not provided in the input  file. The mass apportionment will not be performed.")
            if not self.config['mass_appo'] and mass_appo:
                print("WARNING: EC and OC were provided but the mass_appo flag has not been set to True. Mass apportionment will not be performed.")
                mass_appo = False 
        except KeyError as ke:
            print(f"KeyError: {ke} is not present in the configuration file.")

        # Unique names
        names = np.unique(np.array(names))
        names.tolist()
        # Create the Sample objects
        samples = [Sample(name) for name in names]
        for name in names:
            # Scan for wavelengths and ABS
            # Temporary lists to populate from file
            tmp_wlength, tmp_u_wlength, tmp_abs, tmp_u_abs, tmp_ec, tmp_oc = [], [], [], [], [], [] 
            for n, w, a in zip(rawdata['Name'], rawdata['Wavelength'], rawdata['ABS']):
                if n == name:
                    tmp_wlength.append(float(w))
                    tmp_u_wlength.append(5.0)
                    tmp_abs.append(float(a))
                    if not uncertainty:
                        tmp_u_abs.append(DEFAULT_ABS_UNC * float(a))
            if uncertainty:
                for n, ua in zip(rawdata['Name'], rawdata['uABS']):
                    if n == name:
                        tmp_u_abs.append(float(ua))
            if mass_appo:
                for n, ec, oc in zip(rawdata['Name'], rawdata['EC'], rawdata['OC']):
                    if n == name:
                        tmp_ec.append(float(ec))
                        tmp_oc.append(float(oc))
            for sample in samples:
                if sample.name == name:
                    sample.properties.wavelength = tmp_wlength
                    sample.properties.u_wavelength = tmp_u_wlength
                    sample.properties.abs = tmp_abs
                    sample.properties.u_abs = tmp_u_abs
                    if mass_appo:
                        tmp_ec = tmp_ec[0] # It's all the same values
                        tmp_oc = tmp_oc[0]
                        sample.properties.ec = tmp_ec
                        sample.properties.oc = tmp_oc
        self.data = samples


    def open_json(self):
        """To implement"""
        pass 


# Function definitions
def prepare(config_file):
    """
    Wrapper for the functionality of the DataPrepper class
    """
    if type(config_file) == str:
        with open(config_file) as f:
            config = json.load(f)
    elif type(config_file) == dict:
            config = config_file
    else:
        print("Input format not recognized")

    dp = DataPrepper(config = config)
    dp.get_data()
    return dp.data

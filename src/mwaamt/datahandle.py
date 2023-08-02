###########################################################
# IMPORTS
###########################################################

import sys
import json
import pandas as pd
import numpy as np
from pprint import pprint


###########################################################
# NUMEERICAL CONSTANTS
###########################################################
DEFAULT_ABS_UNCERTAINTY = .05



###########################################################
# PRINTOUT STRINGS
###########################################################
INPUT_FILE_KEY_NOT_PRESENT = "CONFIGURATION FILE ERROR: the key 'input file' was not found in the configuration file"
INPUT_FILE_NO_ABS = "INPUT FILE ERROR: the input file does not contain either a column named ABS nor a column named Babs"
MASS_APPO_NO_EC_OC = "WARNING: the mass_appo flag has been set to True but the EC and OC values were not provided in the input  file. The mass apportionment will not be performed."
MASS_APPO_NO_FLAG = "WARNING: EC and OC were provided but the mass_appo flag has not been set to True. Mass apportionment will not be performed."



###########################################################
# CLASS DEFINITIONS
###########################################################

class Sample():
    """Represents an individual sample under analysis"""
    def __init__(self, name):
        self.name = name
        self.properties = Property()
    def __str__(self):
        # Name string
        ns = f'Sample object - Name {self.name}'
        # Properties string
        ps = str(list(self.properties.__dict__.keys()))
        string_to_print = ns + '\nProperties - ' + ps +'\n'
        return string_to_print

class Property():
    """Container object for all the numerical properties of the sample."""
    def __init__(self):
        """The properties are:

        # Input data
        data_type (string, either "ABS" or "Babs")
        wavelength (list, nm)
        u_wavelength (list, nm, corresponding to wavelength)
        abs (list, corresponding to wavelength)
        u_abs (list, corresponding to wavelength)
        ec (float, ug/cm3)
        oc (float, ug/cm3)

        # From the fit
        aae (float)
        u_aae (float)
        scale (float)
        u_scale (float)
        chisq_aae_fit (float)
        a (float)
        u_a (float)
        b (float)
        u_b (float)
        alpha_brc (float)
        u_alpha_brc (float)
        chisq_type_fit (float)
        a_p (float)
        u_a_p (float)
        b_p (float)
        u_b_p (float)
        chisq_source_fit (float)

        # From the optical apportionment
        brc (list, corresponding to wavelength)
        brc_fraction (list, corresponding to wavelength)
        bc_ff (list, corresponding to wavelength)
        bc_ff_fraction (list, corresponding to wavelength)
        bc_wb (list, corresponding to wavelength)
        bc_wb_fraction (list, corresponding to wavelength)

        # From the mass apportionment
        ec_ff (float)
        ec_ff_fraction (float)
        ec_wb (float)
        ec_wb_fraction (float)
        oc_ff (float)
        oc_ff_fraction (float)
        oc_wb (float)
        oc_wb_fraction (float)
        oc_nc (float)
        oc_nc_fraction (float)

        """
        pass





###########################################################
# DATA READING SECTION
###########################################################

# Function to call for reading the data
def data_read(configuration_file):
    """Read the data from an input file.

    Parameters:
        configuration_file : str,
            the path to the JSON configuration file which contains, 
            among other information, the input data file path under the 
            key 'input file'.

    Returns:
        data : sequence,
            a list of Sample objects, each one corresponding to a physical
            sample for analysis.
    """
    # Open the input data file
    with open(configuration_file, 'r') as f:
        config = json.load(f)
    try:
        input_file = config['input file']
        # For additional measurements on data, e.g. levoglucosan
        additional_measurements = config['additional']
        data_type = config['data type']
        rawdata = pd.read_csv(input_file)
    except KeyError as e:
        print(e)
        print(INPUT_FILE_KEY_NOT_PRESENT)
    # Read the names
    try:
        names = [x for x in rawdata['Name']]
    except NameError as e:
        print(e)
    # Read the wavelength directly from the input file header
    wlength, u_wlength = [], []
    for k in rawdata.keys():
        try:
            wlength.append(int(k))
            u_wlength.append(config['wavelength error'])
        except ValueError:
            pass
    # Check if a mass apportionment is required
    mass_appo_requested = True if 'EC' in rawdata.keys() and 'OC' in rawdata.keys() else False
    # Create the list of sample objects. This will be the output
    data = [Sample(name) for name in names]
    # Fill the Sample.property 
    for sample in data:
        sample.properties.wavelength = [x for x in wlength]
        sample.properties.u_wavelength = [x for x in u_wlength]
        sample.properties.abs, sample.properties.u_abs = [], []
        sample.properties.data_type = data_type
        for wl in wlength:
            # Fill with the absorption values
            for NAME, ABS in zip(rawdata['Name'], rawdata[str(wl)]):
                if NAME == sample.name:
                    sample.properties.abs.append(ABS)
                    sample.properties.u_abs.append(DEFAULT_ABS_UNCERTAINTY * ABS)
        if mass_appo_requested:
            # Fill with carbon concentrations if needed
            for NAME, EC, OC in zip(rawdata['Name'], rawdata['EC'], rawdata['OC']):
                if NAME == sample.name:
                    sample.properties.ec = EC
                    sample.properties.oc = OC
        if len(additional_measurements) > 0:
            for ad_meas_name in additional_measurements:
                for NAME, AD_MEAS in zip(rawdata['Name'], rawdata[ad_meas_name]):
                    if NAME == sample.name:
                        setattr(sample.properties, ad_meas_name, AD_MEAS)
    return data


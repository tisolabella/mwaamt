# MWAAMT - MWAA Model Toolkit
Source code for the toolkit presented in [1]. 

## Installation instructions
* Clone the git repository;
* `cd` to the root `python/` directory;
* run `python3 -m pip install .`

## How to use
After editing the configuration file to suit the specific analysis needs (see below for details), run `mmt <path-to-config-file>`.

## The configuration file
The various options to run the MWAA model analysis are contained in a JSON file, usually named config.json. A template for the configuration file can be found in `python/examples/`. The file format needs to follow the JSON structure. 
Here is an explanation of the various options and parameters. For more information on the various parameters, see [1].
* `input file`: `string`, path to the data file;
* `working directory`: `string`, path to the directory where result files and plots will be saved;
* `alpha bc swipe`: `boolean`, books the $\alpha_{BC}$ swipe analysis;
* `mass appo`: `boolean`, books the mass apportionment procedure;
* `plots`: `boolean`, enables plot saving;
* `data type`: `Babs` for optical absorption coefficient in Mm$^{-1}$, `ABS` for raw absorbance;
* `additional`: `list of strings`, indicates the independent measurements available in addition to the optical data. Currently supports only `["Levoglucosan"]`. If an empty list, the preprocessing step will not be performed;
* `wavelength error`: `float`, the instrument-dependent uncertainty on the wavelength;
* `presets`: `list of strings`, indicates what preset to use for the preprocessing step. Currently supports only `["levoglucosan"]`;
* `iterations`: `int`, the number of iterations in the preprocessing step;
* `threshold`: `float`, the R$^2$ threshold value for the preprocessing step;
* `alpha BC`, `alpha FF`, `alpha WB`: `float`, the starting values for the three AAE values. If the preprocessing step has not been booked, these will be the fixed values for the absorption exponents;
* `AAE low` and `AAE high`: `float`, the AAE boundary values for the particulate to be considered EC in the mass apportionment procedure (ignored if `k1` and `k2` are assigned a value);
* `k1` and `k2`: `float`, the linear regression coefficients for the mass apportionment procedure; if either of these is assigned a `0`, the fitting procedure is booked, otherwise the provided values will be used for mass apportionment.  

## References
[1] Isolabella, T., Bernardoni, V., Bigi, A., Brunoldi, M., Mazzei, F., Parodi, F., Prati, P., Vernocchi, V., and Massabò, D.:_A new software toolkit for optical apportionment of carbonaceous aerosol_, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-1936, 2023. 

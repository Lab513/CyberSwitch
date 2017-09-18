This repository contains the code and data that were used in the paper:
Lugagne, J.-B., Carillo, S. S., Kirch, M., Köhler, A., Batt, G., & Hersen, P. (2017). Balancing a genetic toggle switch by real-time feedback control and periodic forcing. Nature Communications.

All experimental and simulated data is provided in the Data folder.
The parameters that were with our model are saved in the two following files:
parameters1.mat: parameter values that result from fitting to calibration data.
parameters2.mat: parameter values that result from fitting to characterization data and control experiment (includes manual adjustments).

The folder ParameterSearch contains the necessary code for fitting our model to the experiemntal calibration data provided in the data folder.
The scrip main_param_search.m executes the whole workflow for fitting the data and displaying the fitting results.

The folder ModelSimulation contains all necessary functions for simulating dynamic experiments with our model.
The two main functions are generate_data.m and generate_data_explicit_input.m

The folder Control contains all functions and objects necessary for generating control simulations. (With the functions in ModelSimulation/).
The main function is mainControl.m

The Utilities folder contains a collection of auxiliary functions that are used in other parts of the code.

You can re-plot all the graphs of the paper with the script PaperFigures.m (I recommend executing it section by section, some parts of the script take time.)
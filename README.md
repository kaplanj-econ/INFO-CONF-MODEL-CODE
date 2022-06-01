# INFO-CONF-MODEL-CODE
C++ code for simulating HLB ABM economic model
The files in the repository will allow the user to simulate the HLB ABM economic model.

ADDITIONAL LIBRARIES REQUIRED
=================================
Please install the boost and cereal libraries into the headers folder before compilation.

HOW TO USE
=================================
Once the dependencies are installed, use the makefile to compile the program. The program can then be run from the command line. The program can be executed on its own, in which the default config files located at configs/bioConfig.json and configs/econConfig.json are used for the biology config file and the economic config file, respectively. To provide your own config files with custom parameters, run the program as so:

[executable name] [path to econ config] [path to bio config]

CONTACT
==================================
Please direct any questions about this program to Jonathan Kaplan (kaplanj@csus.edu).

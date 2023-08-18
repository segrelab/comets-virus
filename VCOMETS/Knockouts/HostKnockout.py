import cobra
import matplotlib.pyplot as plt
import pandas as pd

import sys
sys.path.append('/projectnb/cometsfba/pythonlibs/lib/python3.8/site-packages')

import cometspy_virus_test as c
    
covid_cobra = cobra.io.read_sbml_model('../Models/iAB_AMO1410_SARS-CoV-2.xml')

covid_cobra.objective = 3392

#----------------------#

output = open("HostKnockout.txt", "w")

for r in range(len(covid_cobra.reactions)):
    rxn = covid_cobra.reactions[r]
    
    lb = rxn.lower_bound
    rb = rxn.upper_bound
    
    rxn.lower_bound = rxn.upper_bound = 0
    
    solution = cobra.flux_analysis.pfba(covid_cobra)
    
    output.write(str(solution.fluxes['BIOMASS_mac']) + '\n')
    
    rxn.lower_bound = lb
    rxn.upper_bound = rb
    
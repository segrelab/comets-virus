import cobra
import matplotlib.pyplot as plt
import pandas as pd

import sys
sys.path.append('/projectnb/cometsfba/pythonlibs/lib/python3.8/site-packages')

import cometspy_virus_test as c
    
covid_cobra = cobra.io.read_sbml_model('../Models/iAB_AMO1410_SARS-CoV-2.xml')
vbof = covid_cobra.reactions[3393]
    
#Note that the coefficient file leaves the ATP coefficient as 0
coefficientFile = open('../Models/Coefficients.txt')
new_coefficients = eval(coefficientFile.readline())
old_coefficients = eval(coefficientFile.readline())

for m in vbof.metabolites:
    vbof.add_metabolites({m: new_coefficients[str(m)] - old_coefficients[str(m)]})
        
#Add the coefficent for ATP
atp_c = covid_cobra.metabolites[225]
vbof.add_metabolites({atp_c: -23.063351651177445})

#Define the metabolites for lipids in the viral membrane
sphmyln = covid_cobra.metabolites[1535]
pchol = covid_cobra.metabolites[1528]

#Add these lipids to the VBOF (both primary and secondary) according to calculated coefficients
vbof.add_metabolites({sphmyln: -0.20233837286384723, pchol: -0.20233837286384723})



#----------------------#

output = open("VirusKnockout.txt", "w")

for r in range(len(covid_cobra.reactions)):
    rxn = covid_cobra.reactions[r]
    
    lb = rxn.lower_bound
    rb = rxn.upper_bound
    
    rxn.lower_bound = rxn.upper_bound = 0
    
    solution = cobra.flux_analysis.pfba(covid_cobra)
    
    output.write(str(solution.fluxes['VBOF']) + '\n')
    
    rxn.lower_bound = lb
    rxn.upper_bound = rb
    
#Runner B has dynamic host biomass for virus biomass calculations

import cobra
import matplotlib.pyplot as plt
import pandas as pd

#Import COMETS from the right location
import sys
sys.path.append('/projectnb/cometsfba/pythonlibs/lib/python3.8/site-packages')

import cometspy_virus_test as c



#Minimum ATP requirement thresholds to test
thresholds = [1]

#World media glucose concentrations to test
#glcConcentrations = [100, 500, 1000e9] #in mols/cm^3
glcConcentrations = [500] #in mols/cm^3

#Maximum Uptake Rate
uptakes = [1]

#Virus resource partition ratios to test
partitions = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]

#Option to keep metabolite concentration constant
for staticMetabolites in [True]:

    for uptake in uptakes:
        for glcConc in glcConcentrations: 
            for vmin in thresholds:

                #New output file for each ATP threshold and glucose concentration
                if (not staticMetabolites):
                    output = open("ResourcePartitionC|uptake=" + str(uptake) + "|vmin=" + str(vmin) + "|Glc=" + str(glcConc) + ".csv", "w")
                else:
                    output = open("ResourcePartitionC|static|uptake=" + str(uptake) + "|vmin=" + str(vmin) + "|Glc=" + str(glcConc) + ".csv", "w")

                output.write('Glucose Concentration,ATP Threshold,Virus Partition Ratio,Host Lifespan,Virus Mass,Virus Rate\n') #Write header


                for r in partitions:

                    vResource = r
                    cResource = 1 - r

                    #Load the model file into COBRA
                    covid_cobra = cobra.io.read_sbml_model('../../Models/iAB_AMO1410_SARS-CoV-2.xml')


                    #Define ATP for virion production (ATPv)
                    atp_v = cobra.Metabolite(
                        'atp_v',
                        formula='C10H12N5O13P3',
                        name='ATP C10H12N5O13P3',
                        compartment='v')

                    #Define ATP for host cell maintenance (ATPh)
                    atp_h = cobra.Metabolite(
                        'atp_h',
                        formula='C10H12N5O13P3',
                        name='ATP C10H12N5O13P3',
                        compartment='h')

                    #Store the current VBOF as the primary one
                    vbof = covid_cobra.reactions[3393]

                    #Define the secondary VBOF
                    vbofb = cobra.Reaction('VBOFb')
                    vbofb.name = 'VBOF secondary'
                    vbofb.lower_bound = 0
                    vbofb.upper_bound = 1000

                    #Change the coefficients (stored in the file "Coefficients.txt") of old metabolites in the VBOF
                    #First line of file is a dictionary of new coefficients for each reaction
                    #Second line of file is a dictionary of old coefficients for each reaction

                    #Note that the new VBOF does not use any atp_c (cytosolic ATP) anymore so the coefficient for it is 0
                    coefficientFile = open('../../Models/Coefficients.txt')
                    new_coefficients = eval(coefficientFile.readline())
                    old_coefficients = eval(coefficientFile.readline())

                    for m in vbof.metabolites:
                        vbof.add_metabolites({m: new_coefficients[str(m)] - old_coefficients[str(m)]})
                        vbofb.add_metabolites({m: new_coefficients[str(m)]})

                    #Add the coefficents for the newly partitioned ATP (ATPv and ATPh)
                    vbof.add_metabolites({atp_v: -23.063351651177445}) #Standard VBOF consumes ATPv
                    vbofb.add_metabolites({atp_h: -23.063351651177445}) #Secondary VBOF consumes ATPh

                    #Define the metabolites for lipids in the viral membrane
                    sphmyln = covid_cobra.metabolites[1535]
                    pchol = covid_cobra.metabolites[1528]

                    #Add these lipids to the VBOF (both primary and secondary) according to calculated coefficients
                    vbof.add_metabolites({sphmyln: -0.20233837286384723, pchol: -0.20233837286384723})
                    vbofb.add_metabolites({sphmyln: -0.20233837286384723, pchol: -0.20233837286384723})

                    #Add the secondary VBOF to the model
                    covid_cobra.add_reactions([vbofb])

                    #Define ATP Maintenance Reaction (not in the model)
                    atpm = cobra.Reaction('ATPM')
                    atpm.name = 'ATP Maintenance'
                    atpm.lower_bound = 0.  # This is the default
                    atpm.upper_bound = 1000.  # This is the default

                    h2o_c = covid_cobra.metabolites[496]
                    adp_c = covid_cobra.metabolites[169]
                    h_c = covid_cobra.metabolites[494]
                    pi_c = covid_cobra.metabolites[702]

                    atpm.add_metabolites({ #Note that ATPM consumes ATPh (host ATP)
                        atp_h: -1,
                        h2o_c: -1,
                        adp_c: 1,
                        h_c: 1,
                        pi_c: 1
                    })

                    covid_cobra.add_reactions([atpm])

                    #Define the ATP conversion reaction for ATPv --> ATPh
                    atpvconversion = cobra.Reaction('ATPVConversion')
                    atpvconversion.name = 'ATP_v Conversion'
                    atpvconversion.lower_bound = 0.  # This is the default
                    atpvconversion.upper_bound = 1000.  # This is the default

                    atpvconversion.add_metabolites({
                        atp_v: -1,
                        atp_h: 1
                    })

                    covid_cobra.add_reactions([atpvconversion])

                    #Define the ATP partition reaction for ATPc --> ATPv + ATPh
                    atppartition = cobra.Reaction('ATPpartition')
                    atppartition.name = 'ATP Partitioning'
                    atppartition.lower_bound = 0.  # This is the default
                    atppartition.upper_bound = 1000.  # This is the default

                    atp_c = covid_cobra.metabolites[225]

                    atppartition.add_metabolites({
                        atp_c: -1,
                        atp_v: vResource,
                        atp_h: cResource
                    })
                    covid_cobra.add_reactions([atppartition])


                    #Load model into COMETS
                    covid = c.model(covid_cobra)
                    covid.initial_pop = [0, 0, 5e-6]

                    #Set objective weights to preserve a hierarchy for the reactions
                    covid.change_objective('VBOFb', 1)
                    covid.change_objective('ATPM', 1000)
                    covid.change_objective('VBOF', 10e6)

                    #Set the biomass to host biomass and adjust the ATP maintenance requirement threshold
                    covid.change_biomass('BIOMASS_mac')
                    covid.change_maintenance('ATPM', vmin)
                    covid.change_bounds('ATPM', 0, vmin)

                    if (uptake != 1):
                        #covid.change_bounds('EX_arg_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        covid.change_bounds('EX_glc_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_gln_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_his_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_ile_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_leu_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_lys_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_met_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_o2_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_ocdca_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_ocdcea_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_phe_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_pi_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_pyr_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_thr_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_trp_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_ttdca_LPAREN_e_RPAREN_', uptake, 1000.0)
                        #covid.change_bounds('EX_val_DASH_L_LPAREN_e_RPAREN_', uptake, 1000.0)


                    #Use parsimonious FBA
                    covid.change_objective_style("MAX_OBJECTIVE_MIN_TOTAL")

                    #Set the world media
                    test_tube = c.layout()
                    test_tube.add_model(covid)

                    if (not staticMetabolites):
                        test_tube.set_specific_metabolite('arg_DASH_L_e', 19.9e-9)
                        test_tube.set_specific_metabolite('glc_DASH_D_e', glcConc * 1e-9) #500e-9
                        test_tube.set_specific_metabolite('gln_DASH_L_e', 0)
                        test_tube.set_specific_metabolite('his_DASH_L_e', 10e-9)
                        test_tube.set_specific_metabolite('ile_DASH_L_e', 40.05e-9)
                        test_tube.set_specific_metabolite('leu_DASH_L_e', 40.05e-9)
                        test_tube.set_specific_metabolite('lys_DASH_L_e', 39.9e-9)
                        test_tube.set_specific_metabolite('met_DASH_L_e', 10.05e-9)
                        test_tube.set_specific_metabolite('o2_e', 1000)
                        test_tube.set_specific_metabolite('h2o_e', 1000)
                        #test_tube.set_specific_metabolite('ocdcea_e', 1000)#
                        test_tube.set_specific_metabolite('phe_DASH_L_e', 20e-9)
                        test_tube.set_specific_metabolite('pi_e', 45e-9)
                        test_tube.set_specific_metabolite('pyr_e', 0)
                        test_tube.set_specific_metabolite('thr_DASH_L_e', 39.9e-9)
                        test_tube.set_specific_metabolite('trp_DASH_L_e', 3.92e-9)
                        #test_tube.set_specific_metabolite('ttdca_e', 1000)#
                        test_tube.set_specific_metabolite('val_DASH_L_e', 40.15e-9)
                        #test_tube.set_specific_metabolite('Tyr_DASH_ggn_c', 1000)#
                        #test_tube.set_specific_metabolite('ocdca_e', 1000)#

                    else:
                        test_tube.set_specific_static_at_location('arg_DASH_L_e', (0, 0), 19.9e-9)
                        test_tube.set_specific_static_at_location('glc_DASH_D_e', (0, 0), glcConc * 1e-9) #500e-9
                        test_tube.set_specific_static_at_location('gln_DASH_L_e', (0, 0), 0)
                        test_tube.set_specific_static_at_location('his_DASH_L_e', (0, 0), 10e-9)
                        test_tube.set_specific_static_at_location('ile_DASH_L_e', (0, 0), 40.05e-9)
                        test_tube.set_specific_static_at_location('leu_DASH_L_e', (0, 0), 40.05e-9)
                        test_tube.set_specific_static_at_location('lys_DASH_L_e', (0, 0), 39.9e-9)
                        test_tube.set_specific_static_at_location('met_DASH_L_e', (0, 0), 10.05e-9)
                        test_tube.set_specific_static_at_location('o2_e', (0, 0), 1000)
                        test_tube.set_specific_static_at_location('h2o_e', (0, 0), 1000)
                        #test_tube.set_specific_static_at_location('ocdcea_e', (0, 0), 1000)#
                        test_tube.set_specific_static_at_location('phe_DASH_L_e', (0, 0), 20e-9)
                        test_tube.set_specific_static_at_location('pi_e', (0, 0), 45e-9)
                        test_tube.set_specific_static_at_location('pyr_e', (0, 0), 0)
                        test_tube.set_specific_static_at_location('thr_DASH_L_e', (0, 0), 39.9e-9)
                        test_tube.set_specific_static_at_location('trp_DASH_L_e', (0, 0), 3.92e-9)
                        #test_tube.set_specific_static_at_location('ttdca_e', (0, 0), 1000)#
                        test_tube.set_specific_static_at_location('val_DASH_L_e', (0, 0), 40.15e-9)
                        #test_tube.set_specific_static_at_location('Tyr_DASH_ggn_c', (0, 0), 1000)#
                        #test_tube.set_specific_static_at_location('ocdca_e', (0, 0), 1000)#


                    # Set the parameters
                    sim_params = c.params()

                    sim_params.set_param('defaultVmax', 18.5)
                    sim_params.set_param('defaultKm', 0.000015)
                    sim_params.set_param('maxCycles', 1000)
                    sim_params.set_param('timeStep', 0.01)
                    sim_params.set_param('spaceWidth', 0.01)
                    sim_params.set_param('maxSpaceBiomass', 10)
                    sim_params.set_param('minSpaceBiomass', 1e-11)
                    sim_params.set_param('writeMediaLog', True)
                    sim_params.set_param('writeFluxLog', True)
                    sim_params.set_param('FluxLogRate', 1)

                    #Run the experiment
                    experiment = c.comets(test_tube, sim_params)
                    try:
                        experiment.run()
                    except:
                        continue

                    #Compute the length of the simulation ()
                    lifespan = len(experiment.total_biomass)

                    #Virus flux is the sum of the primary and secondary VBOF fluxes
                    experiment.fluxes_by_species['macrophage_SARS_CoV_2']['VBOF'] += \
                        experiment.fluxes_by_species['macrophage_SARS_CoV_2']['VBOFb']
                    
                    initialVirusRate = experiment.fluxes_by_species['macrophage_SARS_CoV_2']['VBOF'][0]

                    #Compute the total virus biomass given the flux at each cycle
                    deltaT = 0.01
                    totalVirusBiomass = 0

                    for i in range(0, len(experiment.fluxes_by_species['macrophage_SARS_CoV_2']['VBOF'])):
                        t = i * deltaT
                        totalVirusBiomass = (totalVirusBiomass) + \
                                (experiment.fluxes_by_species['macrophage_SARS_CoV_2']['VBOF'][i] * deltaT * \
                                 experiment.total_biomass['macrophage_SARS_CoV_2'][i])


                    #Write data to output file
                    output.write(str(glcConc) + ", " + str(vmin) + ", " + str(vResource) + ", " + 
                                 str(lifespan) + ', ' + str(totalVirusBiomass) + ", " + str(initialVirusRate) + '\n')
                    output.flush()




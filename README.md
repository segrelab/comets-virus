# comets-virus

## Figures 4-5
FluxComparer compares the fluxes between the model with and without lipids. It uses ReactionLabels' data to categorize reactions and generate Figures 4 and 5.

## Figure 6
VirusRunner is the main notebook that sets up a SARS-CoV-2 infected macrophage experiment (including partitioning ATP and enabling host maintenance). It graphs the virus growth and host death over time. DynamicModelGrapher synthesizes data generated from VirusRunner with different partition ratios.

## Figure 7
PartitionRatioRunnerC compares generates data about virus biomass and host lifespan based on different partition ratios and a set host maintenance. GrapherC reads and plots the data.

## Figure 8
PartitionRatioRunnerB compares generates data about virus biomass and host lifespan based on different partition ratios and different host maintenance values. GrapherB reads and plots the data.

## Knockouts
HostKnockout and VirusKnockout generate data about host/virus flux for each gene knockout. KnockoutComparer then determines potential drug targets based on those data. PartialKnockout tries 10% and 50% knockout on those potential drug targets to further test their efficacy.

# EcoFish


![What is this](Images/fish2.png)

Image Citation: https://gen-fish.ca/uses-and-limitations-of-environmental-dna-edna-in-fisheries-management/

This Repository contains code used for analysis of Data and visualizations of Environmental DNA measurements taken from EcoFish Ltd datasets.

Four main set of codes will be available.

The first set of code is an analysis of a density related experiment.
The second set of code is the analysis of a dilution related experiment.
The third set of code is from the field where we apply covariate analysis.
Finally, we include code for the testing and outlier search for several eDNA related primers.


# Density Experiment 
The goal was to investigate the relationship between
Transformed CT values ('TCT') and Coho salmon 'density' (or the similar measurement
of biomass). CT is the Cycle Threshold and we define TCT=50-CT. The
experiment consisted of manipulating juvenile Coho Salmon densities in treatments of
1, 2, 4, 8, 16, 32, and 65 fish in replicated 10,000L tanks.

![What is this](Images/tctdensity.png)



# Dilution Experiment

In this experiment, three juvenile Coho were allowed to acclimate to four tanks (tanks 19, 20, 21 and 24) and were subsequently removed. After the Coho were removed, water was allowed to flow out of each tank at a known rate as it was replaced, or diluted, by hatchery water. Measurements were obtained at several intervals of varying levels of flow. Moreover, several control samples were taken from the hatchery kitchen sink and from the hatchery pond.

![What is this](Images/TCTflow.png)

# Field Study

To study how eDNA analysis can be used to detect Coho in the wild, several field studies were conducted. In particular,
four streams in British Columbia were studied (steams AAA, BBB, CCC and DDD).
Water samples were collected and associated environmental and physical covariates
were recorded. The main fish of interest was Coho Salmon (Oncorhynchus kisutch),
although we also considered Cutthroat Trout (Oncorhynchus clarkii ) and Rainbow
Trout (Oncorhynchus mykiss). EcoFish personnel also took biomass measurements
of each of the species caught using electrofishing.

![What is this](Images/pcaimage.png)

# Distribution Checks

Analysis of several eDNA related qPCR targets. Checking for consistency and outliers under several distributions.
![What is this](Images/distributions.png)


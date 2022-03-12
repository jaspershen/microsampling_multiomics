# Load data and helper functions -- required for all scripts
source("scripts/load-data.R")

# HR and Temp (Figure 1)
source("scripts/figure1.R") # TODO: I can't run this -- missing file "Users/jessilyn/Desktop/framework_paper/Figure1/Slide 2/slide2_C_participant_data_summary.csv"

# Make tables for supplementary materials
# source("scripts/tables.R")

# Population models and CCA (CCA goes to Figure 1 and )
source("scripts/population-models.R")
source("scripts/group-comparison.R")
source("scripts/cca.R")

# Plot summary metrics of features
source("scripts/figure2F.R") # TODO: I can't run this -- missing file "20180622/20180621_DayPrior_noDemog_RF_Features.csv"
source("scripts/figure3.R") # TODO: I can't run this, errors
source("scripts/figure3sup.R")

# Create ranked list of clinical laboratory tests by the correlation coefficients between observed and predicted values
source("scripts/individual-models.R")

# Individual models (Figure 5)
source("scripts/figure5.R") # TODO: I can't run this: missing files framework_timecourse/*

# Supplementary figures
source("scripts/supplementary-figures.R")

# Supplementary tables
source("scripts/supplementary-tables.R")

# Not in paper
source("scripts/not-in-paper.R")

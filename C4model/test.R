# Test the C4model

# Load packages
library(R.utils)

# Source model code and functions
source("C4model.R")
sourceDirectory("functions/", modifiedOnly = FALSE)

# Test model
C4model()

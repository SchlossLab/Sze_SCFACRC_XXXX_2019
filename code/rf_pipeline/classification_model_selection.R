# Author: Begum Topcuoglu
# Date: 2019-01-14
######################################################################
# Description:
# This function defines defines:
#     1. Tuning budget as a grid the classification methods chosen
#     2. Cross-validation method
#     3. Caret name for the classification method chosen
######################################################################

######################################################################
# Dependencies and Outputs:
#    Filename to put to function:
#     "Random_Forest"


# Usage:
# Call as source when using the function. The function is:
#   tuning_grid()

# Output:
#  List of:
#     1. Tuning budget as a grid the classification methods chosen
#     2. Cross-validation method
#     3. Caret name for the classification method chosen
######################################################################


######################################################################
#------------------------- DEFINE FUNCTION -------------------#
######################################################################

tuning_grid <- function(model, n_features){

  # Cross-validation method
  cv <- trainControl(method="repeatedcv",
                     repeats = 100,
                     number=5,
                     returnResamp="final",
                     classProbs=TRUE,
                     summaryFunction=twoClassSummary,
                     indexFinal=NULL,
                     savePredictions = TRUE)

  # Grid and caret method defined for random forest classification model
  if(model=="Random_Forest"){
		if(n_features > 1000) {powers <- 5^(0:10)
		} else if(n_features > 30){powers <- 4^(0:10)
		} else if(n_features > 15){powers <- 2^(0:10)
		}	else if(n_features > 7)	{ powers <- seq(1, 15, by=2)
		} else { powers <- 1:7 }

		powers <- powers[powers <= n_features]
    grid <-  expand.grid(mtry = powers)
    method = "rf"
  }
  else {
    print("Model not available")
  }
  params <- list(grid, method, cv)
  return(params)
}

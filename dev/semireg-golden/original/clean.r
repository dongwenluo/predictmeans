########## R-function: clean ##########

# Removes missing vales (NA's)
# from an object.

# Last changed: 13 NOV 2003

clean <- function(x)
   return(x[is.na(x)==0])

######### End of clean ##########



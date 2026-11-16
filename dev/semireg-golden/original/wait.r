########## R-function: wait ##########

# Causes the computer to wait
# before continuing.

# Last changed: 13 NOV 2003

   wait <- function()
   {
      cat("Hit Enter to continue\n")
      ans <- readline()
      invisible()
   }

######### End of wait ##########


# Radionuclide class definition.
# Started 9/2/2010 by Jared Nance (nancejk@uw.edu)

setClass(
         # This is the class definition of a radionuclide.
         Class="radionuclide",

         # A radionuclide has a name (chemical name, such as
         # Pb210), a half life, and some other characteristics.
         representation=representation(
           # Isotope name
           name="character",
           # Half life in seconds
           hl="numeric"
           ),

         # The default isotope is lead 210
         prototype=prototype(
           # Our generic isotope is lead-210
           name="Pb210",
           # With a half life of 22.3 years
           hl=22.3*31556926
           ),

         # Validity function.  Should be called at construction
         # and any time the object is modified!
         validity=function(object)
         {
           # No negative half lives.
           if( object@hl <= 0 ) {
             stop("Half life cannot be negative!")
           }

           # OK
           return(TRUE)
         }

         # end of definition
         )

# Initializer method.  verifies that the half life is a positive number,
# etc
setMethod(
          f="initialize",
          signature="radionuclide",
          definition=function(.Object,name="Pb210",hl=22.3*31556926)
          {
            # Set the fields
            .Object@name <- name
            .Object@hl <- hl

            # Return the object
            validObject(.Object)
            return(.Object)
          }
          )

### Radioactive decay.  This just computes exp(-t/hl).

# Generic def
setGeneric(
           name="decayFactor",
           def=function(.Object,...)
           {
             standardGeneric("decayFactor")
           }
           )

# Real function def
setMethod(
          f="decayFactor",
          signature="radionuclide",
          # Time here is in seconds
          definition=function(.Object,time=1)
          {
            res <- numeric()
            for( element in 1:length(time) ) {
              res <- c(res, exp(-log(2)*time[element]/.Object@hl))
            }
             
            return(res)
          }
          )

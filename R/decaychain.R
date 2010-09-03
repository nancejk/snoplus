# Decay Chain class and method definitions.
# Started 9/2/2010 by Jared Nance (nancejk@uw.edu)

# Needs the radionuclide definition
source('radionuclide.R')

setClass(
         # Class name is decaychain
         Class="decaychain",

         # Each decay chain has a name, a list of
         # nuclei in the chain, and an initial
         # activity list.
         representation=representation(
           name="character",
           chain="list",
           activity0 = "numeric"
           ),

         # Validating function.
         validity=function(object)
         {
           # Each isotope in the chain must have an initial
           # activity for t=0
           if( length(object@chain) != length(object@activity0) ) {
             cat(length(object@chain))
             cat(length(object@activity0))
             stop("All isotopes must have an initial activity specified!")
           }

           # Each element of the chain also needs to be
           # a radionuclide, otherwise that doesn't make much sense
           for( element in object@chain ) {
             if( class(element) != "radionuclide" ) {
               stop("Some element of the decay chain is not a radionuclide.")
             }
           }

           # No negative activities are permitted of course
           for( element in object@activity0 ) {
             if( element < 0 ) {
               stop("Initial activities must be positive numbers.")
             }
           }

           # OK
           return(TRUE)
         }
         )

### Initializer method.  Takes a list of radioisotopes and
### initial activities and builds a decay chain out of it.
setMethod(
          f="initialize",
          signature="decaychain",
          definition=function(.Object,isotope1,isotope2,activity)
          {
            # Assign slots
            .Object@chain[1] <- isotope1
            .Object@chain[2] <- isotope2
            .Object@activity0 <- activity

            # Validate
            validObject(.Object)

            # Finish
            return(.Object)
          }
          )

### Calculate the activity of the elements of the chain
### at time t in seconds.  The initial activity of the
### species will be used for t=0.

# Generic method def
setGeneric(
           name="calcActivities",
           def=function(.Object,...) {
             standardGeneric("calcActivities")
           }
           )

# Real method def
setMethod(
          f="calcActivities",
          signature="decaychain",
          definition=function(.Object,t=1.0)
          {
            isotope1 <- .Object@chain[[1]]
            isotope2 <- .Object@chain[[2]]
            ac1 <- .Object@activity0[1]
            ac2 <- .Object@activity0[2]
            # Result will be the same length as the initial
            # activities
            res <- data.frame(numeric(), numeric())
            names(res) <- c(isotope1@name, isotope2@name)

            # The parent of the n-th isotope is the n-1th
            # isotope.  The 1st isotope is the parent of the
            # entire chain.
            for( element in 1:length(t) ) {
              # Parent decays away exponentially
              term1 <- ac1*decayFactor(isotope1,t[element])

              # Child builds up due to parent and dies off exponentially
              ratio <- isotope2@hl/(isotope2@hl - isotope1@hl)
              sum <- decayFactor(isotope2,t[element]) - decayFactor(isotope1,t[element])
              term2 <- ac2*decayFactor(isotope2,t[element]) + ac1*ratio*sum
              res[element,] <- c(term1, term2)
            }
            
            # Return
            return(res)
          }
          )

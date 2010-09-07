# Decay Chain class and method definitions.
# Started 9/2/2010 by Jared Nance (nancejk@uw.edu)

# Needs the radionuclide definition
source('radionuclide.R')

# The class 'decaychain'.  This is really meant to be
# an abstract class definition, in that it only does
# things that are independent of the number of elements
# in the chain, and the subclasses below (twoelementchain,
# threeelementchain, etc etc) should be used for actual
# calculations.
setClass(
         # Class name is decaychain
         Class="decaychain",

         # Each decay chain has a name, a list of
         # nuclei in the chain, and an initial
         # activity list.
         representation=representation(
           name="character",
           chain="list",
           activity0 = "numeric",
           numisotopes = "numeric"
           ),

         # Validating function.
         validity=function(object)
         {
           # Each isotope in the chain must have an initial
           # activity for t=0
           if( length(object@chain) != length(object@activity0) ) {
             stop("All isotopes must have an initial activity specified!")
           }

           # Make sure the length of the chain and the number of
           # elements is consistent
           if( length(object@chain) != object@numisotopes ) {
             stop("Length of chain not equal to number of elements somehow.")
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


# Generic method def for calculating chain activities.
setGeneric(
           name="calcActivities",
           def=function(.Object,...) {
             standardGeneric("calcActivities")
           }
           )

#####################################
### Two element chain definition. ###
#####################################
setClass(
         Class="twoelementchain",
         contains="decaychain"
         )

### Initialize method for two element chain
### to build, call new("twoelementchain", isotope1, isotope2, c(ac1, ac2))
setMethod(
          f="initialize",
          signature="twoelementchain",
          definition=function(.Object,isotope1,isotope2,activity)
          {
            # Assign slots
            .Object@chain[1] <- isotope1
            .Object@chain[2] <- isotope2
            .Object@activity0 <- activity
            .Object@numisotopes <- 2

            # Validate
            validObject(.Object)

            # Finish
            return(.Object)
          }
          )

# Calculates the activities in a two element chain for
# any time t in seconds
setMethod(
          f="calcActivities",
          signature="twoelementchain",
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

###########################
### Three element chain ###
###########################
setClass(
         Class="threeelementchain",
         contains="decaychain"
         )

### Initialization for three element chain
setMethod(
          f="initialize",
          signature="threeelementchain",
          definition=function(.Object,isotope1,isotope2,isotope3,activity)
          {
            # Assign slots
            .Object@chain[1] <- isotope1
            .Object@chain[2] <- isotope2
            .Object@chain[3] <- isotope3
            .Object@activity0 <- activity
            .Object@numisotopes <- 3

            # Validate
            validObject(.Object)

            # OK
            return(.Object)
          }
          )

### Calculation of activity for three element chain for any time t in seconds.
### A nasty business.
setMethod(
          f="calcActivities",
          signature="threeelementchain",
          definition=function(.Object,t=1.0)
          {
            i1 <- .Object@chain[[1]]
            i2 <- .Object@chain[[2]]
            i3 <- .Object@chain[[3]]
            a1 <- .Object@activity0[1]
            a2 <- .Object@activity0[2]
            a3 <- .Object@activity0[3]

            # Result will be the same length as the initial
            # activities
            res <- data.frame(numeric(), numeric(), numeric())
            names(res) <- c(i1@name, i2@name, i3@name)

            # GO
            for( el in 1:length(t) ) {
              ### Primary decay
              # Radioactive decay of primary
              a1p <- i1@hl*a1*decayFactor(i1,t[el])

              ### Secondary decay
              # Radioactive decay of intermediate species
              term1 <- a2*decayFactor(i2,t[el])

              # Child builds up due to parent and dies off exponentially
              ratio <- i2@hl/(i2@hl - i1@hl)
              sum <- decayFactor(i2,t[el]) - decayFactor(i1,t[el])
              term2 <- a1*ratio*sum
              a2p <- i2@hl*(term1 + term2)

              ### Ternary decay
              # First term  is just radioactive decay
              term1 <- a3*decayFactor(i3,t[el])

              # Second term is interplay between elements 2 and 3
              ratio <- i3@hl/(i3@hl-i2@hl)
              sum <- decayFactor(i3,t[el]) - decayFactor(i2,t[el])
              term2 <- a2*ratio*sum

              # Third term is due to primary element decay
              ratio <- i2@hl*i3@hl/(i1@hl - i2@hl)
              sum1 <- (decayFactor(i2,t[el]) - decayFactor(i3,t[el]))/(i2@hl-i3@hl)
              sum2 <- (decayFactor(i3,t[el]) - decayFactor(i1,t[el]))/(i3@hl-i1@hl)
              term3 <- a1 * ratio * (sum1 + sum2)
              a3p <- i3@hl*(term1 + term2 + term3)

              # Assign row
              res[el,] <- c(a1p,a2p,a3p)
            }

            # Return
            return(res)
          }
          )

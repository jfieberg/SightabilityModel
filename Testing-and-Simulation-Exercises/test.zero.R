# Test that values of 0 are handled properly, especially for cases where a block has 
# no group seen and so we must insert a record with all zeros in that block so that
# that block is counted

library(SightabilityModel)


# Create a survey with two strata, but in one of the strata, all of the group values are zero

Stratum.info.csv <- textConnection(
"Stratum, Stratum.Area, Stratum.Blocks
A, 1000,  100
B, 1000,  100
C, 1000,  100")

Stratum.info <- read.csv(Stratum.info.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)

Block.info.csv <- textConnection(
"Block.ID,  Block.Area
1,     10
2,     10
3,     10
4,     10
5,     10")

Block.info <- read.csv(Block.info.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)

Group.info.csv <- textConnection(
"Block.ID, Stratum, Bulls, Cows, Total
1,  A, 1, 2, 3
1,  A, 0, 0, 0
2,  A, 3, 6, 9
2,  A, 0, 0, 0
3,  B, 1, 2, 3
3,  B, 0, 0, 0
4,  B, 0, 0, 0
5,  C, 1, 2, 3
5,  C, 0, 0, 0")

Group.info <- read.csv(Group.info.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)

options(survey.lonely.psu="average")

# Test the MoosePopR() function
# We see that the abundance for stratum B is half that of stratum C because
# stratum B has 2 blocks surveyed (one with no groups)
# If you just deleted block 4, you would get a positively biased estimate.
MoosePopR( survey.data       = Group.info,
           survey.block.area = Block.info,
           stratum.data      = Stratum.info,
           abundance=~Total,
           survey.lonely.psu="remove")

MoosePopR( survey.data       = Group.info,
           survey.block.area = Block.info,
           stratum.data      = Stratum.info,
           density=~Total,
           survey.lonely.psu="remove")

# unfortunately, we cannot deal with lonely PSU in a ratio estimator
# because it is not clear how to set the covariance and coef methods for this
# wrong object
#MoosePopR( survey.data       = Group.info,
#           survey.block.area = Block.info,
#           stratum.data      = Stratum.info,
#           numerator=~Bulls, denominator=~Cows,
#          survey.lonely.psu="remove")
#
# returns
#Error in .fun(piece, ...) : 
#  Sorry, unable to deal with lonely PSU when finding a ratio estimate 

# Test the SightabilityPopR() function
# The abundance in stratum B should be 1/2 of that from stratum A.
# We can estimate the variance due to sampling for strata with lonely PSU - hot sure how that is possible
# except is must come from multiple groups within the same block so likely is an 
# underestimate. The default sighting adjustment is 100% visibility with very small standard error
# so these results should be very similar to MoosePopR() above.
SightabilityPopR( survey.data       = Group.info,
           survey.block.area = Block.info,
           stratum.data      = Stratum.info,
           abundance=~Total)

SightabilityPopR( survey.data       = Group.info,
           survey.block.area = Block.info,
           stratum.data      = Stratum.info,
           density=~Total)

SightabilityPopR( survey.data       = Group.info,
           survey.block.area = Block.info,
           stratum.data      = Stratum.info,
           numerator=~Bulls, denominator=~Cows)


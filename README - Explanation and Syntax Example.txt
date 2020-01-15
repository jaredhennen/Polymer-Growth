These programs were used for a project with the goal of finding the distribution of
polymer chain lengths (end to end distance) given some polymer consisting of two
monomers grown in 2 dimensions. The monomers consists of two bonds of equal length
with a fixed bond angle, where the bond angle is the only parameter that differs
between each monomer. Two solutions were found, a first approach where overlapping
bonds was allowed (as this can be most directly extended to 3D growth) and a second
approach requiring no overlap.


To create a computational distribution for overlap allowed polymer growth and compare
it to the analytical distribution run the following code in MATLAB

>> simulateandcomparetodist(10000,8,2,110/180*pi,128/180*pi,1)

To create a computational distribution of polymer growth when overlap is NOT allowed
and compare it to the analytical distribution for overlap allowed runthe following
code in MATLAB

>> simulateandcomparenooverlapwithoverlap(100,4,1,110/180*pi,128/180*pi,1)

Those two programs should provide the syntax and comments necessary to run and
understand subfunctions

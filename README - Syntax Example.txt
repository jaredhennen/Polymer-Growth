To create a computational distribution for overlap allowed polymer growth and compare
it to the analytical distribution (eg slide 25) run the following code in MATLAB

>> simulateandcomparetodist(10000,8,2,110/180*pi,128/180*pi,1)

To create a computational distribution for overlap NOT allowed polymer growth and
compare it to the analytical distribution for overlap allowed (eg slide 29) run the
following code in MATLAB

>> simulateandcomparenooverlapwithoverlap(100,4,1,110/180*pi,128/180*pi,1)

Those two programs should provide the syntax and comments necessary to run and
understand subfunctions
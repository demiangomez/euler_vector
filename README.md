# Euler vector calculation: A class to calculate and plot Euler vectors of rotation 
This is the code related to the paper "Linear system approach to Euler vector inversion" from Smalley, Gómez, and Bevis.
The only dependency is m_map: https://www.eoas.ubc.ca/~rich/map.html

As input, the class can take as first argument:

  a) text file,
  
  b) a matrix with the same format as the file, see below,
  
  c) a folder with json files from Parallel.GAMIT.
  
 A second argument (optional) provides the list of stations to use to 
 determine the vector. Finally, a name for the vector (also optional)
 
 If using a text file, it needs to have the following columns 
  (in order, table headers are ignored):
 stnm   x   y z  n_vel  e_vel    OR
 stnm lat lon 0  n_vel  e_vel    
 
 stnm       : name of the station
 x y z      : meters
 lat lon    : decimal degrees
 n_vel e_vel: meters/year (do NOT use mm/yr)

 The Euler vector will be calculated using all stations or using a subset
 determined by the either a cell array with the stations names
 or another text file with a list of station names.
 The Euler vector is determined using the technique explained in Smalley
 et al. (xxxx) and a robost least squares algorithm 
 that relies on a goodness of fit test to remove outliers. This is why
 you should not use mm/yr instead of m/year. The fit starts with very
 loose weights for each station (equivalent to 1 m/yr) and it reweighs
 the dataset on several iterations removing outliers.

 EXAMPLES:
 using a list of json files from Parallel.GAMIT
   evect = euler_vector('../Solutions', 'href.txt');
 using a text file as input
   evect = euler_vector('sam.dat', 'href.txt');
 or just a single argument to use all the station in the file
   evect = euler_vector('sam.dat');
 the list of sites can also be provided using a cell array
   evect = euler_vector('sam.dat', {'stn1', 'stn2', 'stn3' ...});
 plot the result and histograms of residuals
   evect.plot_result(250);
 to predict a station's velocity, use
   vel = evect.compute_vel(x)
 where x accepts the same inputs as the object
 if the prediction and a plot is desired, then use
   vel = evect.plot_prediction(files, 100, true);
 where 100 is the scale for the velocity vectors, and true determines
 if error ellipses should be plotted or not

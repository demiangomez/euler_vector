clear
close all

% read the locations where we want to estimate V
sg = readtable('sg.dat');
lla = [table2array(sg(:, 2:3)) zeros(height(sg), 1) table2array(sg(:, 4:5))];

% now using the dat file, selecting stations from href.txt
epole1 = euler_pole('sam.dat', 'href.txt', 'dat with HREF');

epole1.plot_result(250);

% compute and plot the velocities
v2 = epole1.plot_prediction(lla, 100, true);

% now using the dat file,  but without specifying a list of stations
epole2 = euler_pole('sam.dat', [], 'data without HREF');

epole2.plot_result(250);
v3 = epole2.plot_prediction(lla, 100, true);

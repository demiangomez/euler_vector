clear
epole1 = euler_pole('example.dat', [], 'dat without HREF');

fprintf('(without HREF) wx, wy, wz (mas/yr): %.3f %.3f %.3f\n', epole1.pole_X * 1e-9 * 180/pi * 3600 * 1000)
fprintf('(without HREF) sigma wx, wy, wz (mas/yr): %.3f %.3f %.3f\n', sqrt(diag(epole1.cov_xyz)) * 1e-9 * 180/pi * 3600 * 1000)
fprintf('(without HREF) pole position: %.3f %.3f %.3f\n', epole1.pole)
fprintf('(without HREF) pole sigmas  : %.3f %.3f %.3f\n', epole1.pole_sigma_lat, epole1.pole_sigma_lon, epole1.pole_sigma_rot)

% create a custom made plot
v_scale = 250;
lon_lim = [-142.74      -24.952];
lat_lim = [-56.694       10.252];

m = [0.09 0.03];

figure(10)
clf

subplot_tight(2,3,1, m)
m_proj('Mercator','lon', lon_lim, 'lat', lat_lim)
hold on
% data vectors
fr(1) = m_quiver([epole1.stations(:).lon]', [epole1.stations(:).lat]', ...
    [epole1.stations(:).ve]' .* v_scale, [epole1.stations(:).vn]' .* v_scale, 0,'b');

% modeled
fr(2) = m_quiver([epole1.stations(:).lon]', [epole1.stations(:).lat]', ...
    [epole1.vm(:).ve] .* v_scale, [epole1.vm(:).vn] .* v_scale, 0, 'r');

% plot the location of the pole
fr(3) = m_plot(epole1.pole(2), epole1.pole(1), 'or','MarkerFaceColor','r');

% plot ellipse
epole1.plot_ellipse(epole1.cov_lla, epole1.pole, 0, 0, 30)

grid on
m_coast('color',[0 .6 0]);
m_grid('tickdir','out', 'fontsize', 12);

legend(fr, 'Data', 'Model', 'Euler pole', 'location', 'southwest')
title('(a) Map of velocities and Euler pole (without href)')

% histogram of errors in N
subplot_tight(2,3,2, m)
histogram([epole1.residuals(:).vn]' * 1000)
title('(b) North residuals')
grid on
ylabel('Frequency')

% histogram of errors in E
subplot_tight(2,3,3, m)
histogram([epole1.residuals(:).ve]' * 1000)
title('(c) East residuals')
set(gca,'yticklabel',{[]})
grid on

% ==================================
% redo but now using the HREF file
% ==================================

epole1 = euler_pole('example.dat', 'href.txt', 'dat with HREF');

fprintf('(with HREF) wx, wy, wz (mas/yr): %.3f %.3f %.3f\n', epole1.pole_X * 1e-9 * 180/pi * 3600 * 1000)
fprintf('(with HREF) sigma wx, wy, wz (mas/yr): %.3f %.3f %.3f\n', sqrt(diag(epole1.cov_xyz)) * 1e-9 * 180/pi * 3600 * 1000)
fprintf('(with HREF) pole position: %.3f %.3f %.3f\n', epole1.pole)
fprintf('(with HREF) pole sigmas  : %.3f %.3f %.3f\n', epole1.pole_sigma_lat, epole1.pole_sigma_lon, epole1.pole_sigma_rot)

subplot_tight(2,3,4, m)
m_proj('Mercator','lon', lon_lim, 'lat', lat_lim)
hold on
% data vectors
fr(1) = m_quiver([epole1.stations(:).lon]', [epole1.stations(:).lat]', ...
    [epole1.stations(:).ve]' .* v_scale, [epole1.stations(:).vn]' .* v_scale, 0,'b');

% modeled
fr(2) = m_quiver([epole1.stations(:).lon]', [epole1.stations(:).lat]', ...
    [epole1.vm(:).ve] .* v_scale, [epole1.vm(:).vn] .* v_scale, 0, 'r');

% plot the location of the pole
fr(3) = m_plot(epole1.pole(2), epole1.pole(1), 'or','MarkerFaceColor','r');

% plot ellipse
epole1.plot_ellipse(epole1.cov_lla, epole1.pole, 0, 0, 30)

grid on
m_coast('color',[0 .6 0]);
m_grid('tickdir','out', 'fontsize', 12);

legend(fr, 'Data', 'Model', 'Euler pole', 'location', 'southwest')
title('(d) Map of velocities and Euler pole (with href)')

% histogram of errors in N
subplot_tight(2,3,5, m)
histogram([epole1.residuals(:).vn]' * 1000)
title('(e) North residuals')
grid on
ylabel('Frequency')
xlabel('[mm/yr]')

% histogram of errors in E
subplot_tight(2,3,6, m)
histogram([epole1.residuals(:).ve]' * 1000)
title('(f) East residuals')
set(gca,'yticklabel',{[]})
grid on
xlabel('[mm/yr]')

set(gcf,'color','w')
exportgraphics(gcf,'paper/fig_3.png','Resolution',300);
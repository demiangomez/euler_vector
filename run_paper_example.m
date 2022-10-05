clear
evect1 = euler_vector('example.dat', [], 'dat without HREF');

fprintf('(without HREF) wx, wy, wz (mas/yr): %.3f %.3f %.3f\n', evect1.vector_X * 1e-9 * 180/pi * 3600 * 1000)
fprintf('(without HREF) sigma wx, wy, wz (mas/yr): %.3f %.3f %.3f\n', sqrt(diag(evect1.cov_xyz)) * 1e-9 * 180/pi * 3600 * 1000)
fprintf('(without HREF) pole position: %.3f %.3f %.3f\n', evect1.pole)
fprintf('(without HREF) pole sigmas  : %.3f %.3f %.3f\n', evect1.pole_sigma_lat, evect1.pole_sigma_lon, evect1.pole_sigma_rot)

% create a custom made plot
v_scale = 200;
lon_lim_b = [-142.74  -24.952];
lat_lim   = [-56.694   10.252];

lon_lim_s = [-84      -24.952];

m = [0.06 0.05];
r = 4;
c = 6;

figs = reshape((1:r*c)', c, r)';

figure(10)
clf
f = figs(1:3, 1:3);
subplot_tight(r,c,f(:), m)
m_proj('Mercator','lon', lon_lim_s, 'lat', lat_lim)
hold on
% data vectors
fr(1) = m_quiver([evect1.stations(:).lon]', [evect1.stations(:).lat]', ...
    [evect1.stations(:).ve]' .* v_scale, [evect1.stations(:).vn]' .* v_scale, 0,'b');

% modeled
fr(2) = m_quiver([evect1.stations(:).lon]', [evect1.stations(:).lat]', ...
    [evect1.vm(:).ve] .* v_scale, [evect1.vm(:).vn] .* v_scale, 0, 'r');

grid on
m_coast('color',[0 .6 0], 'linewidth', 1.5);
m_grid('tickdir','out');

legend(fr, 'Data', 'Model', 'location', 'northeast')

title('(A-i) Velocities and Euler pole (without href)')

f = figs(4, 1);
% histogram of errors in N
subplot_tight(r,c,f,m)
histogram([evect1.residuals(:).vn]' * 1000)
title('(A-ii) N residuals')
grid on
ylabel('Frequency')
xlabel('[mm/yr]')
    
% histogram of errors in E
subplot_tight(r,c,figs(4,2), m)
histogram([evect1.residuals(:).ve]' * 1000)
title('(A-iii) E residuals')
%ylabel('Frequency')
xlabel('[mm/yr]')
grid on

% save to plot later
evectA.cov_lla  = evect1.cov_lla;
evectA.pole     = evect1.pole;
evectA.stations = evect1.stations;
evectA.vm       = evect1.vm;

% ==================================
% redo but now using the HREF file
% ==================================

evect1 = euler_vector('example.dat', 'href.txt', 'dat with HREF');

fprintf('(with HREF) wx, wy, wz (mas/yr): %.3f %.3f %.3f\n', evect1.vector_X * 1e-9 * 180/pi * 3600 * 1000)
fprintf('(with HREF) sigma wx, wy, wz (mas/yr): %.3f %.3f %.3f\n', sqrt(diag(evect1.cov_xyz)) * 1e-9 * 180/pi * 3600 * 1000)
fprintf('(with HREF) pole position: %.3f %.3f %.3f\n', evect1.pole)
fprintf('(with HREF) pole sigmas  : %.3f %.3f %.3f\n', evect1.pole_sigma_lat, evect1.pole_sigma_lon, evect1.pole_sigma_rot)

f = figs(1:3, 4:end);
subplot_tight(r,c,f(:), m)
m_proj('Mercator','lon', lon_lim_s, 'lat', lat_lim)
hold on
% data vectors
fr(1) = m_quiver([evect1.stations(:).lon]', [evect1.stations(:).lat]', ...
    [evect1.stations(:).ve]' .* v_scale, [evect1.stations(:).vn]' .* v_scale, 0,'b');

% modeled
fr(2) = m_quiver([evect1.stations(:).lon]', [evect1.stations(:).lat]', ...
    [evect1.vm(:).ve] .* v_scale, [evect1.vm(:).vn] .* v_scale, 0, 'r');

grid on
m_coast('color',[0 .6 0], 'linewidth', 1.5);
m_grid('tickdir','out');

legend(fr, 'Data', 'Model', 'location', 'northeast')

title('(B-i) Velocities and Euler pole (with href)')

% ascencion island inset
axes('Position',[.50 .32 .72 .17])
m_proj('Mercator','lon', [-14.9, -14.0], 'lat', [-8.2 -7.0])
hold on
sc = 3;
fr(1) = m_quiver([evect1.stations(:).lon]', [evect1.stations(:).lat]', ...
    [evect1.stations(:).ve]' .* v_scale./sc, [evect1.stations(:).vn]' .* v_scale./sc, 0,'b');
% modeled
fr(2) = m_quiver([evect1.stations(:).lon]', [evect1.stations(:).lat]', ...
    [evect1.vm(:).ve] .* v_scale./sc, [evect1.vm(:).vn] .* v_scale./sc, 0, 'r');
m_coast('color',[0 .6 0], 'linewidth', 1.5);
m_grid('tickdir','out');

f = figs(4, 5);
% histogram of errors in N
subplot_tight(r,c,f,m)
histogram([evect1.residuals(:).vn]' * 1000)
title('(B-ii) N residuals')
grid on
ylabel('Frequency')
xlabel('[mm/yr]')
    
% histogram of errors in E
subplot_tight(r,c,figs(4,6), m)
histogram([evect1.residuals(:).ve]' * 1000)
title('(B-iii) E residuals')
%ylabel('Frequency')
xlabel('[mm/yr]')
grid on

set(gcf,'color','w')


f = figs(4, 3:4);
subplot_tight(r,c,f, m)
m_proj('Mercator','lon', [evect1.pole(2) evect1.pole(2)] + [-5 +5], ...
                  'lat', [evect1.pole(1) evect1.pole(1)] + [-1.5 +1.5])
hold on
% plot ellipse
evect1.plot_ellipse(evectA.cov_lla, evectA.pole, 0, 0, 10)
% plot the location of the pole
h(1) = m_plot(evectA.pole(2), evectA.pole(1), 'or','MarkerFaceColor','r');

% plot ellipse
evect1.plot_ellipse(evect1.cov_lla, evect1.pole, 0, 0, 10, [0 .6 0])
% plot the location of the pole
h(2) = m_plot(evect1.pole(2), evect1.pole(1), 'o', 'color', [0 .6 0],'MarkerFaceColor',[0 .6 0]);

% convert ITRF solution to lat lon
X = [-0.0751   -0.0835   -0.0389] * (1e-6 * pi/180);
[pp(1), pp(2), pp(3)] = cart2euler(X(1), X(2), X(3));
pp = [pp(1:2) pp(3)] * 180/pi;

h(3) = m_plot(pp(2), pp(1), 'xb','MarkerFaceColor','b', 'LineWidth', 2);

grid on
m_coast('color',[0 .6 0]);
m_grid('tickdir','out');

legend(h, 'Euler pole (A)', 'Euler pole (B)', 'Altamimi et al. (2017)', 'location', 'southoutside')

title('(c) Euler poles and error ellipses')


% ascencion island inset
axes('Position',[.025 .32 .72 .17])
m_proj('Mercator','lon', [-14.9, -14.0], 'lat', [-8.2 -7.0])
hold on
sc = 3;
fr(1) = m_quiver([evectA.stations(:).lon]', [evectA.stations(:).lat]', ...
    [evectA.stations(:).ve]' .* v_scale./sc, [evectA.stations(:).vn]' .* v_scale./sc, 0,'b');
% modeled
fr(2) = m_quiver([evectA.stations(:).lon]', [evectA.stations(:).lat]', ...
    [evectA.vm(:).ve] .* v_scale./sc, [evectA.vm(:).vn] .* v_scale./sc, 0, 'r');
m_coast('color',[0 .6 0], 'linewidth', 1.5);
m_grid('tickdir','out');

set(gcf, 'position', [2188          11        1160         962])

%exportgraphics(gcf,'paper/Fig_3.png','Resolution',300);
%exportgraphics(gcf,'paper/Fig_3.eps','Resolution',300);


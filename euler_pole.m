classdef euler_pole
    % euler_pole: class to calculate and plot an Euler pole of rotation 
    % the only dependency is m_map: https://www.eoas.ubc.ca/~rich/map.html
    % As input, the class can take as first argument:
    %  a) text file, 
    %  b) a matrix with the same format as the file, see below
    %  c) a folder with json files from Parallel.GAMIT.
    %  
    % A second argument (optional) provides the list of stations to use to 
    % determine the pole. Finally, a name for the pole (also optional)
    % 
    % If using a text file, it needs to have the following columns 
    %  (in order, table headers are ignored):
    % stnm   x   y z  n_vel  e_vel    OR
    % stnm lat lon 0  n_vel  e_vel    
    % 
    % stnm       : name of the station
    % x y z      : meters
    % lat lon    : decimal degrees
    % n_vel e_vel: meters/year (do NOT use mm/yr)
    %
    % The pole will be calculated using all stations or using a subset
    % determined by the either a cell array with the stations names
    % or another text file with a list of station names.
    % The Euler pole is determined using the technique explained in Smalley
    % et al. (xxxx) and a robost least squares algorithm 
    % that relies on a goodness of fit test to remove outliers. This is why
    % you should not use mm/yr instead of m/year. The fit starts with very
    % loose weights for each station (equivalent to 1 m/yr) and it reweighs
    % the dataset on several iterations removing outliers.
    %
    % EXAMPLES:
    % using a list of json files from Parallel.GAMIT
    %   epole = euler_pole('../Solutions', 'href.txt');
    % using a text file as input
    %   epole = euler_pole('sam.dat', 'href.txt');
    % or just a single argument to use all the station in the file
    %   epole = euler_pole('sam.dat');
    % the list of sites can also be provided using a cell array
    %   epole = euler_pole('sam.dat', {'stn1', 'stn2', 'stn3' ...});
    % plot the result and histograms of residuals
    %   epole.plot_result(250);
    % to predict a station's velocity, use
    %   vel = epole.compute_vel(x)
    % where x accepts the same inputs as the object
    % if the prediction and a plot is desired, then use
    %   vel = epole.plot_prediction(files, 100, true);
    % where 100 is the scale for the velocity vectors, and true determines
    % if error ellipses should be plotted or not
    
    properties
        pole
        pole_X
        pole_sigma_lat
        pole_sigma_lon
        pole_sigma_rot
        cov_xyz
        cov_lla
        stations
        json_files
        A
        residuals
        dof
        n
        vm
        rms
        wrms
        wrms_n
        wrms_e
        pole_name % better to reference the pole itself, not the figure
    end
    
    methods
        function self = euler_pole(varargin)
            %EULER_POLE calculates an Euler pole based on a given input
            %   input_files can be
            %
            
            self.stations = self.init_stations_struct();
            self.json_files = [];
            self.A = [];
            self.residuals = [];
            self.dof = NaN;
            self.vm = [];
            self.rms = [];
            self.wrms = [];
            
            if nargin == 3
                self.pole_name = varargin{3};
            else
                self.pole_name = '(default)';
            end

            if nargin == 0
                error('Invalid number of arguments')
            end
            
            stations_input = varargin{1};
            pole_stations = {};
            
            if nargin > 1
                % more then one input argument
                if ~isempty(varargin{2})
                    % if second argument is not empty, we have a list of
                    % stations to build the HREF
                    pole_stations = varargin{2};
                    
                    if ~iscell(pole_stations)
                        if isfile(pole_stations)
                            % file with station list, load it
                            pole_stations = table2cell(readtable(pole_stations));
                        elseif ischar(pole_stations)
                            pole_stations = cellstr(pole_stations);
                        end
                    end
                end
            end
            
            % determine the input type
            if isfolder(stations_input)
                % a directory, search for json files
                self.json_files = dir(fullfile(stations_input, '*.json'));
                j = 1;
                for i = 1:length(self.json_files)
                    % if stations exists in pole_stations or pole_stations
                    % is empty, read the json file
                    if or(sum(~cellfun(@isempty, regexp(lower(self.json_files(i).name), lower(pole_stations)))), isempty(pole_stations))
                        
                        stn = self.load_json(fullfile(self.json_files(i).folder, self.json_files(i).name));
                        % make sure no NaNs make it in here
                        if ~isnan(stn.ve)
                            self.stations(j) = stn;
                            j = j + 1;
                        end
                    end
                end
            elseif isfile(stations_input)
                % it's a file, load the values
                data = readtable(stations_input);
                
                j = 1;
                for i = 1:height(data)
                    if or(sum(~cellfun(@isempty, regexp(lower(table2cell(data(i, 1))), lower(pole_stations)))), isempty(pole_stations))
                        
                        stn = self.load_table_row(data(i, :));
                        % make sure no NaNs make it in here
                        if ~isnan(stn.ve)
                            self.stations(j) = stn;
                            j = j + 1;
                        end
                    end
                end
            end
            
            self.n = size(self.stations, 2);
            
            % now calculate the pole
            self = self.calculate_pole();
        end
        
        function self = calculate_pole(self)
            % Weigther Euler Pole calculation
            % x: n x 3 matrix with station positions in ECEF coordinates
            % v: n x 3 matrix with station velocities in NEU (can also be just NE)

            self.A = self.make_design([[self.stations(:).x]' [self.stations(:).y]' [self.stations(:).z]']);

             % degrees of freedom of the system
            self.dof = size(self.A, 1) - size(self.A, 2);

            % chi-squared values for 95% confidence and the DOF of this system
            X1 = chi2inv(1-0.05/2, self.dof);
            X2 = chi2inv(0.05/2, self.dof);

            % observations vector
            L = [self.stations(:).vn self.stations(:).ve]';

            % weigh each station acording to longitude. This is because in the west 
            % side of Argentina the velocities are affected by plate boundary deformation =>
            % assume a linear relationship east-west for weights

            pass = false;
            factor = 1;
            limit = 2.5;

            % I am assuming Vx, Vy and Vz have the same weight.
            P = eye(self.n * 2);
            
            iter = 1;
            % this is a very simplistic cycle. On a real problem you wouldn't do this.
            while and(~pass, iter < 10)

                % invert for the model parameters
                X = (self.A' * P * self.A) \ self.A' * P * L;

                % residuals
                V = self.A * X - L;

                % square root of the variance of unit weight
                So = sqrt(V' * P * V / self.dof);

                % estimated chi
                chi = So ^ 2 * self.dof;

                % find the a priori sigma for the observations
                factor = factor .* So;
                % normalize residuals by factor (sigma)
                s = abs(V ./ factor);

                if chi < X2 || chi > X1
                    % if it falls in here it's because it didn't pass the Chi2 test
                    % this is the simplest assumption: reweight using So
                    % downweight outliers exponentially
                    f = ones(size(V));

                    f(s > limit) = 10.^(limit - s(s > limit));
                    % replace with eps: don't let the problem become
                    % unstable!
                    f(f < eps) = eps;
                    % recompute the weights
                    P = diag((f./factor).^2);
                else
                    pass = true;
                    break
                end
                iter = iter + 1;
            end

            if ~pass
                warning(['Goodness of fit test for Euler pole ' self.pole_name ' did not converge!'])
            end

            % modeled velocities
            tvm = self.A * X;
            tvm = reshape(tvm, [self.n 2]);
            self.vm.vn = tvm(:, 1);
            self.vm.ve = tvm(:, 2);
            
            % save the residuals of the pole calculation
            tv = reshape(V, [self.n 2]);
            self.residuals.vn = tv(:, 1);
            self.residuals.ve = tv(:, 2);

            % compute the covariance matrix (we will use it later to
            % compute the uncertainties of the predicted velocities)
            self.cov_xyz = So.^2.*inv(self.A' * P * self.A);

            % RMS error of model
            self.rms = sqrt(sum((V(:)).^2)./(self.n * 2 - 1));

            % WRMS error of model
            % variances updated by the LS estimation
            w=diag(P);
            self.wrms = sqrt(sum(w(:).*(V(:)).^2)./sum(w));
            self.wrms_n = sqrt(sum(w(1:end/2).*(V(1:end/2)).^2)./sum(w(1:end/2)));
            self.wrms_e = sqrt(sum(w(end/2+1:end).*(V(end/2+1:end)).^2)./sum(w(end/2+1:end)));

            % save the pole as is
            self.pole_X = X;
            
            % convert euler pole to lat lon position
            [pp(1), pp(2), pp(3)] = cart2euler(X(1), X(2), X(3));
            
            % now scale result to 1e-9 (see design matrix) and Myr
            pp = [pp(1:2) pp(3) * 1e-9 * 1e6] * 180/pi;
            
            self.pole = pp;

            % calculate the sigmas of the pole (as rotation, lat, lon
            % uncertainties)
            mT_ = sqrt(sum(X.^2));         % rotation angular velocity
            mXY = sqrt(X(1).^2 + X(2).^2); % terms for the partials
            pXZ = X(1) * X(3);             % terms for the partials
            pYZ = X(2) * X(3);             % terms for the partials

            % Propagate the Euler pole covariance 
            % See eq (15) from EPC (Goudarzi et al 2013)
            G = [         X(1)/mT_            X(2)/mT_         X(2)/mT_;
                 -1/mT_.^2*pXZ/mXY   -1/mT_.^2*pYZ/mXY     1/mT_.^2*mXY;
                    -X(2)/(mXY.^2)       X(1)/(mXY.^2)               0];
            
            % propagate
            self.cov_lla = G * self.cov_xyz * G';

            % get the sigmas in degrees
            ps = sqrt(diag(self.cov_lla)) * 180/pi;
            self.pole_sigma_lat = ps(2);
            self.pole_sigma_lon = ps(3);

            % save the rotation sigma in deg/Myr
            self.pole_sigma_rot = ps(1) * 1e-9 * 1e6;
        end
        
        function plot_result(self, v_scale)
            % plot the results from the inversion
            figure('Name', self.pole_name);
            clf
            subplot(2,2,[1 3])
            m_proj('Mercator','lon', [min([self.stations(:).lon self.pole(2)])-10 ...
                                      max([self.stations(:).lon self.pole(2)])+10], ...
                              'lat', [min([self.stations(:).lat self.pole(1)])- 5 ...
                                      max([self.stations(:).lat self.pole(1)])+ 5])

            fr(1) = m_quiver([self.stations(:).lon]', [self.stations(:).lat]', ...
                [self.stations(:).ve]' .* v_scale, [self.stations(:).vn]' .* v_scale, 0,'b');
            
            hold on
            
            fr(2) = m_quiver([self.stations(:).lon]', [self.stations(:).lat]', ...
                [self.vm(:).ve] .* v_scale, [self.vm(:).vn] .* v_scale, 0, 'r');
            
            % plot the location of the pole
            fr(3) = m_plot(self.pole(2), self.pole(1), 'or','MarkerFaceColor','r');
            
            % plot the error ellipse of the pole
            % uses zero velocity so that the ellipse is on top of the pole
            self.plot_ellipse(self.cov_lla, self.pole, 0, 0, 30)

            axis equal
            m_coast('color',[0 .6 0]);
            title('(a) Map of velocities and Euler pole')
            grid on
            m_grid('tickdir','out', 'fontsize', 12);
            legend(fr, 'Data', 'Model', 'Euler pole')

            % histogram of errors in N
            subplot(2,2,2)
            histfit([self.residuals(:).vn]' * 1000)
            title('(b) North Histogram')
            grid on
            
            % histogram of errors in E
            subplot(2,2,4)
            histfit([self.residuals(:).ve]' * 1000)
            title('(c) East Histogram')
            grid on
            xlabel('[mm/yr]')
            
            % plot the RMS and WRMS
            %subplot(3,2,6)
            %stem([1 2],[self.rms self.wrms]*1000)
            %xlim([0 3])
            %title('RMS and WRMS')
            %ylabel('Error [mm]')
            %set(gca,'xtick',[1 2],'xticklabel',{'RMS', 'WRMS'});
            %grid on
        end
        
        function vel = compute_vel(self, x)
            % compute the velocities at sites given in x
            % x can be lla or xyz
            
            % pass through function to determine if input is json o array
            x = self.determine_input(x);
            
            Ai = self.make_design(x);
            
            tvm = Ai * self.pole_X;
            tvm = reshape(tvm, [size(x, 1) 2]);
            
            % propagate the uncertainties
            sig = Ai * self.cov_xyz * Ai';
            % sig is a matrix that has a first block for the N uncertainty
            % and a second block for the E uncertainty
            sign = sig(1:size(x, 1), 1:size(x, 1));
            sige = sig(size(x, 1)+1:end, size(x, 1)+1:end);
            % the NE covariances of each predicted velocity are in the
            % diagonal of the covariance block
            cov = sig(1:size(x, 1), size(x, 1)+1:end);
            
            vel.vm.vn = tvm(:, 1);
            vel.vm.ve = tvm(:, 2);
            
            % populate the covariances into the output structure
            vel.vm.vn_sig = sqrt(diag(sign));
            vel.vm.ve_sig = sqrt(diag(sige));
            vel.vm.ne_cov = diag(cov);

            % store the covariance matrix for each prediciton in a 
            % 3D matrix: each slice is the matrix for each prediction
            vel.vm.cov(1,1,:) = diag(sign);
            vel.vm.cov(2,2,:) = diag(sige);
            vel.vm.cov(2,1,:) = diag(cov);
            vel.vm.cov(1,2,:) = diag(cov);

            if size(x, 2) >= 5
                vel.residuals.vn = x(:, 4) - vel.vm.vn;
                vel.residuals.ve = x(:, 5) - vel.vm.ve;
                vel.vn = x(:, 4);
                vel.ve = x(:, 5);
            else
                vel.residuals.vn = nan(size(x, 1), 1);
                vel.residuals.ve = nan(size(x, 1), 1);
                vel.vn = nan(size(x, 1), 1);
                vel.ve = nan(size(x, 1), 1);
            end
        end
        
        function vel = plot_prediction(self, x, v_scale, plot_ellipses)
            % compute and plot the predicted velocities at given sites
            % x can be lla or xyz (n x 3) or a cell array with json files
            % if cols in x > 3, then vn ad ve are assumed to be passed to
            % in which case, residuals will be plotted too
            % v_scale is to scale the vectors
            
            if nargin < 4
                plot_ellipses = false;
            end

            % pass through function to determine if input is json o array
            x = self.determine_input(x);
            
            vel = self.compute_vel(x);
            
            [~, lla] = self.xyz_lla(x);
            
            lla(:, 1:2) = lla(:, 1:2) * 180/pi;
            
            figure('Name', self.pole_name)
            clf
            m_proj('Mercator','lon', [min(lla(:,2))- 5 max(lla(:,2))+ 5], ...
                              'lat', [min(lla(:,1))- 2 max(lla(:,1))+ 2])
            
            m_quiver(lla(:, 2), lla(:, 1), ...
                     vel.ve * v_scale, vel.vn * v_scale, 0, 'b')
            hold on
            m_quiver(lla(:, 2), lla(:, 1), ...
                     vel.vm.ve .* v_scale, vel.vm.vn .* v_scale, 0, 'r')
            
            m_quiver(lla(:, 2), lla(:, 1), ...
                     vel.residuals.ve .* v_scale, ...
                     vel.residuals.vn .* v_scale, 0, 'k')
            
            axis equal
            m_coast('color',[0 .6 0]);
            if all(~isnan(vel.vn))
                title('Model: red; Data: blue')
            else
                title('Model: red')
            end
            xlabel('lon')
            ylabel('lat')
            grid on
            m_grid('tickdir','out', 'fontsize', 12);
            
            if plot_ellipses
                for i = 1:size(lla,1)
                    % plot the ellipse
                    self.plot_ellipse(vel.vm.cov(:,:,i), lla(i,:), ...
                        vel.vm.vn(i), vel.vm.ve(i), v_scale)
                end
            end
        end

        function plot_ellipse(self, cov, lla, vn, ve, scale)
            % function to plot the error ellipses
            % it assumes the user is sending a 2x2 covariance matrix
            % the matrix should have the shape [varN  covNE]
            %                                  [covEN varE ]
            % lla, ve, ve, and scale are used to place the ellipse at the
            % end of the quiver (and scaled by scale * 3)

            [eivec, eival] = eig(cov);
            
            % from matlab's documentation:
            % By default eig does not always return the eigenvalues and 
            % eigenvectors in sorted order.
            [eival, ind] = sort(diag(eival));
            
            % sort the eigenvectors
            eivec = eivec(ind,ind);

            max_evc = eivec(:, end);
            
            % Get the largest eigenvalue
            max_evl = eival(end);
            
            % Get the smallest eigenvector and eigenvalue
            min_evl = eival(1);
            
            % Calculate the angle between the N-axis and the 
            % largest eigenvector
            angle = atan2(max_evc(2), max_evc(1));
            
            % Make sure the angle is between 0 and 2pi
            if angle < 0
                angle = angle + 2*pi;
            end
            
            % Get the 95% confidence interval error ellipse using the
            % Fisher distribution
            chisquare_val = sqrt(2*finv(0.95,2,self.dof));

            gtheta = linspace(0, 2*pi);

            a      = chisquare_val*sqrt(max_evl);
            b      = chisquare_val*sqrt(min_evl);
            
            % determine the scale to apply to the ellipse so that it is
            % proportional to the vector plotted by m_quiver
            % this piece of code was extracted from m_quiver
            [XN,YN]=m_ll2xy([lla(:,2) lla(:,2)]', ...
                [lla(:,1) lla(:,1)+.001]','clip','off');
            [XE,YE]=m_ll2xy([lla(:,2) ...
                lla(:,2)+(.001)./cos(lla(:,1)*pi/180)]', ...
                [lla(:,1) lla(:,1)]','clip','off');
            
            mU=a.*reshape(diff(XE),size(lla(:,1)))*1000 + ...
                b.*reshape(diff(XN),size(lla(:,1)))*1000;
            mV=a.*reshape(diff(YE),size(lla(:,1)))*1000 + ...
                b.*reshape(diff(YN),size(lla(:,1)))*1000;

            mVe = ve.*reshape(diff(XE),size(lla(:,1)))*1000 + ...
                vn.*reshape(diff(XN),size(lla(:,1)))*1000;
            mVn = ve.*reshape(diff(YE),size(lla(:,1)))*1000 + ...
                vn.*reshape(diff(YN),size(lla(:,1)))*1000;

            % the ellipse in x and y coordinates 
            % apply an addition scale factor to the size of the ellipse
            ellipse_N  = mU*cos(gtheta) * scale * 3;
            ellipse_E  = mV*sin(gtheta) * scale * 3;
            
            % Rotation matrix
            R = [ cos(angle) -sin(angle)
                  sin(angle)  cos(angle) ];
            
            %let's rotate the ellipse to some angle phi
            r_ellipse = [ellipse_N
                         ellipse_E]' * R;
            
            % determine the lat lon in the plot's coordinate system
            [x, y] = m_ll2xy(lla(:,2), lla(:,1));

            % Draw the ellipse
            plot(r_ellipse(:,1) + x + mVe * scale, ...
                 r_ellipse(:,2) + y + mVn * scale, '-r')
        end
    end
    
    methods(Static)
        function stn_struct = init_stations_struct()
            stn_struct = struct('name', '', 'lat', NaN, 'lon', NaN, ...
                'x', NaN, 'y', NaN, 'z', NaN, 'vn', NaN, 've', NaN);
        end
        
        function stn_struct = load_json(file)
            % open a json file from Parallel.GAMIT
            ts = jsondecode(fileread(file));
            
            stn_struct.name = lower(ts.Station);
            stn_struct.lat = ts.lat;
            stn_struct.lon = ts.lon;
            stn_struct.x = ts.ref_x;
            stn_struct.y = ts.ref_y;
            stn_struct.z = ts.ref_z;
            
            if ismember('Polynomial', fieldnames(ts))
                neu = ts.Polynomial.params(:,2);
                
                stn_struct.vn = neu(1);
                stn_struct.ve = neu(2);
            else
                stn_struct.vn = NaN;
                stn_struct.ve = NaN;
            end
        end
        
        function stn_struct = load_table_row(data)

            t = table2struct(data);
            % get the name of the fields to make this function independent
            % of headers
            fields = fieldnames(t);
            
            stn_struct.name = lower(t.(fields{1}));
            
            p = [t.(fields{2}) t.(fields{3}) t.(fields{4})];

            % detect the type of input
            [x, lla] = euler_pole.xyz_lla(p);
            
            % lat lon input
            stn_struct.lat = lla(:,1) * 180/pi;
            stn_struct.lon = lla(:,2) * 180/pi;
            
            % ecef input
            stn_struct.x = x(:,1);
            stn_struct.y = x(:,2);
            stn_struct.z = x(:,3);
            
            % velocity of stations
            stn_struct.vn = t.(fields{5});
            stn_struct.ve = t.(fields{6});
        end
        
        function [x, lla] = xyz_lla(lla_or_x)
            % this function return both lla and xyz given any of both
            n = sqrt(sum(lla_or_x.^2, 2));
            
            if n(1) > 6000e3
                % vector in xyz, get lla
                lla = euler_pole.x2g(lla_or_x');
                lla = lla';
                x = lla_or_x;
            else
                % discards the height
                lla = lla_or_x(:, 1:2) * pi/180;
                % + 6370e3 because heights are given and g2x wants radius
                x = euler_pole.g2x([lla(:, 1:2)'; lla_or_x(:, 3)' + 6370e3]);
                x = x';
            end    
        end
        
        function x = determine_input(x)
            % determine the type of input an reshape it to match our format
            if iscell(x)
                % a cell array with json files!
                for i = 1:length(x)
                    stn(i) = euler_pole.load_json(x{i});
                end
                x = [stn(:).lat; 
                     stn(:).lon; 
                     zeros(1, i); 
                     stn(:).vn; 
                     stn(:).ve]';
            end
        end
        
        function A = make_design(x)
            % make design can accept any vector in either lat lon or XYZ
            % dimenensions must be n x 3
            
            [~, lla] = euler_pole.xyz_lla(x);
            
            n = size(lla,1);
            
            slat = sin(lla(:,1));
            slon = sin(lla(:,2));
            clat = cos(lla(:,1));
            clon = cos(lla(:,2));

            slatslon = slat.*slon;
            slatclon = slat.*clon;
            
            Re=6378137;
            
            % build the design matrix
            A=Re*[slon -clon zeros(n,1); -slatclon -slatslon clat]/1e9;
        end

        function Vg = x2g(Vx)
            %x2g transforms vector(s) V from cartesian coordinates (x,y,z)
            %    to spherical coordinates (lat,lon,r).

            %    SEND: Vx = vector V stated in global cartesian form (x,y,z)
            % RETURNS: Vg = vector V stated in spherical form (lat,lon,r)

            %  SHAPES:  If transforming one vector at a time then Vx can be a
            %  a column or row matrix, and Vg will have same shape. If you are 
            %  transforming m vectors at one time then the individual vectors
            %  should be organized as the columns of Vx. In this case both Vx
            %  and Vg will be 3 by m matrices. 
            %   NOTES:  (lat,lon) returned in radians.
            % VERSION: 1.0              (17 Sept '91)            Mike Bevis

            [a, b]=size(Vx);
            if a == 1         % if Vx is a row matrix (row vector)
              Vx=Vx';         % then reshape it into a column matrix
              Vg=Vx ;         % copy shape of Vx
            end
            x=Vx(1,:); y=Vx(2,:); z=Vx(3,:);
            Vg(1,:)= atan2(z, sqrt(x.^2 + y.^2));
            Vg(2,:)= atan2(y, x);
            Vg(3,:)= sqrt(x.^2 + y.^2 +z.^2);
            if a == 1          % if Vx was sent as a row matrix
              Vx=Vx'; Vg=Vg';  % then turn Vx and Vg back into row matrices
            end
        end
        
        function Vx = g2x(Vg)
            %g2x spherical geographical coordinates to global cartesian coordinates.
            %    Transforms vector(s) V from geographical coordinates 
            %    (lat,lon,r) to cartesian coordinates (x,y,z)
            %
            %    SEND: Vg = vector(s) V stated in geographical form (lat,lon,r)
            % RETURNS: Vx = vector(s) V stated in global cartesian form (x,y,z)
            %
            %  SHAPES:  If transforming one vector at a time then Vg can be a
            %  a column or row matrix, and Vx will have same shape. If you are 
            %  transforming m vectors at one time then the individual vectors
            %  should be organized as the columns of Vg. In this case both Vg
            %  and Vx will be 3 by m matrices. 
            %   NOTES: (lat,lon) must be given in radians
            % VERSION: 2.0                (1 Oct '91)            Mike Bevis


            [a,b]=size(Vg);
            if a == 1         % if Vg is a row matrix (row vector)
              Vg=Vg';         % then reshape it into a column matrix
            end
            sla=sin(Vg(1,:));  cla=cos(Vg(1,:)); 
            slo=sin(Vg(2,:));  clo=cos(Vg(2,:));
            r=Vg(3,:);
            Vx(1,:)=r.*cla.*clo; Vx(2,:)=r.*cla.*slo; Vx(3,:)=r.*sla;
            if a == 1          % if Vg was sent as a row matrix
              Vg=Vg'; Vx=Vx';  % then turn Vg and Vx back into row matrices
            end
        end
    end
end

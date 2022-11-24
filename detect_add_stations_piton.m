%% 
% detect_add_stations.m
%
% Calculates detectability of different station combinations using a known
% seismic network by adding stations to the network
%
%
% SECTION 1: CALCULATING VARIABLES, SARA AND OUTPUT
%    Part 1) Loading station location, creating other variables
%            (contains values/data that user must input)
%    Part 3) Creating grid of points
%            (contains values/data that user must input)
%    Part 3) Creating grid of added stations
%            (contains values/data that user must input)
%    Part 4) Creating the names of the station pairs ratios
%    Part 5) Calculating ampltitude ratios
%    Part 6) Calculating difference in amplitude ratios
%    Part 7) Calculating detectability
%    Part 8) Calculating volumes 
%    Part 9) Calculating centroids and stds of each combination
%    Part 10) Exporting output data as a table
%
% SECTION 2: PLOTTING OUTPUT
%    Part 11) Plotting 3d grid of points
%    Part 12) Plotting 2d detection grid with scaled volume and scaled
%             centroid shift (relies on importing data from
%             detect_positions.m)
%
%
% EXTERNAL FUNCTIONS NEEDED TO RUN ALL CODE (most availiable from MATHWORKS file exchange):
% deg2utm.m (Rafael Palacios, 2006)
% utm2deg.m (Rafael Palacios, 2006)
% distinguishable_colors.m (Timothy E. Holy, 2010-2011)
% GetSRTMData.m (Sebastian Hölz, 2004)
% axesLabelsAlign3D.m (Matthew Arthington, 2010)
% colorcet.m (Peter Kovesi, 2015)
%
% Katie James, 2020

 close all
 clear all

%% SECTION 1: CALCULATING VARIABLES
 
%% PART 1) LOADING STATION LOCATION DATA AND CREATING OTHER VARIABLES
% This section contains variables that the user must input

%% 1a) Loading in station data and extracting station positions in decimal degrees

% Load in station position data
load Piton_stations.mat % currently as txt data
data = PdF_st_location; % input_data = cell array containing data
Ns = size(data,1); % length of input data/ number of stations

% Extracting longitude, latitudes & elevations in decimal degrees
lat = cell2mat(data(:,2)); % Latitude column of data
lon = cell2mat(data(:,3)); % Longitude column
elev = cell2mat(data(:,6)); % elevations

[utme_o,utmn_o] = deg2utm(lat,lon); % station locations of unedited network
sta = [utme_o utmn_o elev];

max_elev = 0; % max elevation of volcano,used to close isosurface, should be set to zero
sum_elev = 2632; % elevation of summit of volcano

summit_lon = 55.708889;
summit_lat = -21.2425; % coordinates of summmit of volcano
[utme_s,utmn_s] = deg2utm(summit_lat,summit_lon); % summit in utm

%% 1b) Extracting station names from input data

% Extracting station names from input data
sta_name = (data(:,1)); % column of data containing station names
sta_char = 3; % input the number of characters in the station name (all stations names must have the same no of characters)

z_value = 8000; % Arbitrary elevation value, used for plotting station names on relief map
               % - just has to be larger than any of the actual elevations

vector_latlon = [51.7, 51.9,-177.5,-178.3]; % vector of form 
% [latstart, latend, lonstart, lonend], denotes geographic locations of 
% seismic stations, (is needed to create DEM map - not used in grid
% creation)

%% 1c) Loading topography of volcano from srtm data
% Topography is used to define z location of grid of added stations

% insert function that downloads srtm data
data = readhgt(vector_latlon); % data = struct that gives volcano topography

%% 1d) Inputs for calulating B (used in amplitude ratio equation)

Q = 170; % quality factor (orig 170)
f  =10; % Mean value for the frequency used  (orig 10)
beta = 1700; % shear wave velocity
n= 1; % for body waves

B= pi*f./(Q*beta); % B = isotropic medium constant

%% 1e) Picking location that is used in amplitude ratio graph

lon_index = 2; % pick index for longitude (CONVERT TO LOCATION)
lat_index = 2; % pick index for latitude
UTMzone = '40 K'; % must define UTM zone this point is in (for amp loc)

%% 1f) Pick direction of migration, threshold and no of station pairs

mig_dir = 'z'; % indicate the direction of migration (x, y or z);
gs = 1000; % minimum value is 100 m
start_depth = -10000; % gives the starting depth, max depth is 20 km, must be negative,  min is summit elevation + 1000m.

x_horz_gs = 500; % minimum value is 2000m  
y_horz_gs = 500; % minimum value is 2000m  

% define threshold based off of absolute difference in amplitude ratios
thres = 0.1;
st_detect = Ns+1; % surface encloses all gridpoints with st_detect and greater number of station pairs, must be integer

%% PART 2) CREATING GRID THAT CAN BE USED FOR DETECTABILITY
% This section contains variables that the user must input

%% 2a) Create 3d grid

% centroid of stations
cen_station = [mean(sta(:,1)) mean(sta(:,2)) mean(sta(:,3))];

detect_grid = []; % creates a grid of points, migration will be from these points and used in the calculation for detectability

% Creation of grid uses latitude and longitude. USER MUST INPUT LAT AND LON START AND END
% VALUES FOR GRID IN UTM

% utm (m) for latitude and longititude
x0 = cen_station(2)-20e3; % latitude inputs (start), 
x00 = cen_station(2)+20e3 ; % (end)
y0 =  cen_station(1)-20e3 ; % longitude inputs (start)
y00 = cen_station(1)+ 20e3; % (end), 

dx = 500; % spacing between latitude points in m
dy = 500; % spacing between longitude points in m 

% convert slice of lat and lon points 
utmy_slice =[x0-y_horz_gs:dx:x00];
utmx_slice =[y0-x_horz_gs:dy:y00];
lon_slice ={}; lat_slice={};

for i = 1:length(utmx_slice)
    for j=1:length(utmy_slice)
      [lat1,lon1]= utm2deg(utmx_slice(i),utmy_slice(j),UTMzone);
      lon_slice{i,j} = lon1;
      lat_slice{i,j} = lat1;              
    end
end

% convert grid dimensions back to dec deg
[latstart, lonstart]= utm2deg(y0-x_horz_gs,x0-y_horz_gs,UTMzone);
[latend, lonend]= utm2deg(y00,x00,UTMzone);
latdiff = (latend-latstart)/ceil((size(lat_slice,2)-1)/1);
londiff = (lonend-lonstart)/ceil((size(lat_slice,1)-1)/1);

latgrid = [latstart:latdiff:latend];
longrid = [lonstart:londiff:lonend];


% define elevation values in m
out_srtm=[];
for i = 1:length(longrid)
    for j = 1:length(latgrid)
        %data = readhgt(latgrid(j),longrid(i)); % download required srtm tile/s
        out = GetSRTMData([longrid(i) latgrid(j)]); % out gives x,y,z for one point
        out_srtm{i,j} = out;
        out_z(i,j) = out(3);
    end
end

% define elevation values in m
depth_int = 50; % depth interval - must include max_elev
min_depth = max_elev+(sum_elev+1000); % min depth for grid
depth_idx = gs/depth_int; % alters code so uses user set grid spacing gs

x_idx = x_horz_gs/dy;
y_idx = y_horz_gs/dx;

max_depth = (start_depth) - gs; % max depth for grid, with correction for user set grid spacing

layer_mat = max_depth:depth_int:min_depth; % gives 1d matrix of elevation values
z_len = length(layer_mat); % length of elevation data

% use meshgrid to create longitude, latitude and depth vectors
[ymesh,xmesh, zmesh] = meshgrid(y0-x_horz_gs:dy:y00,x0-y_horz_gs:dx:x00,max_depth:depth_int:min_depth); 

utmew = ymesh; % rename lon vector 
utmns = xmesh; % lat vector
z_3d = zmesh; % depth vector
utm_lenns = size(utmns,1); % Length of utm ns data
utm_lenew = size(utmew,2);  % Length of utm ew data

% interpolating topography to make detection surface smooth
latdiff = (x00-x0)/ceil((size(lat_slice,2)-1)/1);
londiff = (y00-y0)/ceil((size(lat_slice,1)-1)/1);
outlon = [y0:londiff:y00];
outlat=[x0:latdiff:x00];

%[outlon outlat] = deg2utm(out_y,out_x);
for i = 1:length(outlon)
    for j = 1:length(outlat)
        outlon_all(i,j) = outlon(i);
        outlat_all(i,j) = outlat(j);
    end
end

Vq = interp2(outlat_all,outlon_all,out_z,utmns(:,:,1),utmew(:,:,1)); %Vq is interpolated topography data


%% 2b) Finding location of point used in seismic amplitude ratio graph

% UTM location of location used in amplitude ratio plot
point_lon = utmew(1,lon_index,1); % longitude value for location in utm
point_lat = utmns(lat_index,1,1); % latitude value for location in utm

% conversion of gridpoints from utm to decimal degrees (used only for DEM
% map)
[point_lat_deg, point_lon_deg] = utm2deg(point_lon, point_lat, UTMzone);


%% Part 3) Creating grid stations will be added too
% This section contains variables that the user must input

%% 3a) creating grid of stations to be added

addstat = 1; % No. of stations to add on each iteration
Ns_add = Ns+addstat; % New total number of seismic stations

% 2d grid (x and y change, z is kept constant)
% has to be within limits of the srtm data (otherwise topography cannot be
% used for added stations)
statx0 = cen_station(1)-8e3;  % longitude inputs
statx00 = cen_station(1)+7e3; 
staty0 = cen_station(2)-8e3; % latitude inputs
staty00 = cen_station(2)+8e3; 

stat_dx = 5e3; % longitude spacing
stat_dy = 4e3; % latitude spacing

[ystat, xstat] = meshgrid(staty0:stat_dy:staty00,statx0:stat_dx:statx00); % creating grid (2d)

g_len_ew = size(ystat,1);
g_len_ns = size(ystat,2);
g_len_tot =  g_len_ns*g_len_ew;

ystat_vec = reshape(ystat,g_len_tot,[]); % 1d form of ystat
xstat_vec = reshape(xstat,g_len_tot,[]); % 1d form of xstat

add_name = 'ADD';
sta_name{Ns_add,:} = add_name; % adds test station name to original dataset


% use srtm elevations for locations of new stations
 xstat_all = []; ystat_all = [];
for i = 1:length(ystat_vec) 
    [ystat_deg, xstat_deg] = utm2deg(xstat_vec(i),ystat_vec(i),UTMzone);
    ystat_all = [ystat_all; ystat_deg];
    xstat_all = [xstat_all; xstat_deg];
end

out_srtm2 =[];
for i = 1:length(ystat_vec)
    out = GetSRTMData([xstat_all(i) ystat_all(i)]);
    out_srtm2 = [out_srtm2; out];
end

z_add = out_srtm2(:,3);

%sets all added station elevations to zero
% for i = 1:g_len_tot
%     z_add(i) = 0;
% end

%% Part 4) Creating the names of the station pairs ratios (e.g. station a / station b)

ratio_data = {}; % contains names of station pairs

    for st1 = 1:Ns_add
        for st2 = 1:Ns_add
            if st2 > st1
                ratio_pair = strcat(sta_name(st2,:),'/',sta_name(st1,:));
                ratio_data = [ratio_data; ratio_pair];             
            end
        end
    end

Ns2 = size(ratio_data,1);
ratio_data_mat = cell2mat(ratio_data); 

%% 5) Calculating amplitude and amplitude ratios
% amplitude is now not attenuated - each location has a source amplitude

%% 5a) Calculating amplitude for each station

final_all_d = {};
final_all_amp = {};
src_amp = 1; % amplitude of events
sta_cell = {};

% for loop to calculate distance and amplitude between each grid point and station
for g1 = 1:g_len_tot
        sta(Ns_add,1) = xstat_vec(g1); % sta gives coords of all stations (original+added)
        sta(Ns_add,2) = ystat_vec(g1);
        sta(Ns_add,3) = z_add(g1);
        sta_cell = [sta_cell; sta];
    
    for st = 1:Ns_add
             all_amp = [];
             all_d = [];

        for k = 1:utm_lenns
            for i = 1:utm_lenew

                d = sqrt((utmew(k,i,1)-sta(st,1)).^2+(utmns(k,i,1)-sta(st,2)).^2+(z_3d(k,i,:)-sta(st,3)).^2); 
                amp = src_amp.*(exp(-B*d)./(d.^n));
                all_amp = [all_amp; amp];
                all_d = [all_d; d];

               end
            end

        all_d2 = reshape(all_d,utm_lenew,utm_lenns,z_len); %reshape to right shape
        all_amp2 = reshape(all_amp,utm_lenew,utm_lenns,z_len);

        all_d2 = 1./permute(all_d2,[2 1 3]); %transpose and flip standard form
        all_amp2 = 1./permute(all_amp2,[2 1 3]);
       % now amplitude ratios are in correct format

        final_all_d = [final_all_d all_d2];
        final_all_amp = [final_all_amp all_amp2];


    end
end

dist = final_all_d; % cell containing distances between eahc grid point and the staion
amp = final_all_amp; % cell containing amplitude for each grid point each station

amp = reshape(amp,Ns_add,(g_len_tot)); % each column is a diff combo
dist = reshape(dist,Ns_add,(g_len_tot));
    
    
%% 5b) Calculating amplitude ratios for each station pair

% for loop to calculate amplitude ratios
amp_ratio = {}; % cell containing amplitude ratios

for clen = 1:(g_len_tot)
    amp_use = amp(:,clen);
    
    for st1 = 1:Ns_add
        for st2 = 1:Ns_add
            if st2 > st1 
                  amp_pair = log10(amp_use{st1,:}./amp_use{st2,:}); % amplitude ratio without
                  amp_ratio = [amp_ratio; amp_pair];
            end 
        end
    end   
end

amp_ratio = reshape(amp_ratio, Ns2, (g_len_tot));

amp_ratio_all = amp_ratio;
ratio_data_all = ratio_data;

%% 6) Calculating differences between amp ratios for a migration

%% 6a) Calculating amplitude ratio differences

ratio_diff = {}; % gives absolute difference for amp ratios for each station pair
statpair = [];
diff =[];

ratio_diff_per = {}; 
per_change = []; % gives percentage change for amp ratios for each station pair

for stn = 1:g_len_tot
for i =1:Ns2
    statpair = cell2mat(amp_ratio_all(i,stn));
                
    for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m =1:z_len
                
                if mig_dir == 'x'
                     if j > x_idx   
                   diff(j,k,m) = statpair(j,k,m)-statpair(j-x_idx,k,m); % x direction , e-w 
                    end
                   
                else if mig_dir == 'y'
                        if k>y_idx
                        diff(j,k,m) = statpair(j,k,m)-statpair(j,k-y_idx,m); 
                        end
                        
                else if mig_dir == 'z'
                        if m >depth_idx
                        diff(j,k,m) = (statpair(j,k,m)-statpair(j,k,m-depth_idx));
                        end      

                end
                end
                end
                end
                end
                end
    ratio_diff{i,stn} = abs(diff); 
end
end

%creating new station name array
for i = 1:Ns2
    ratio_name_sum(i,:) = ratio_data_all(i,:);
end


%% 7) Detection capability and detection plot for every station combination
%% 7a) calculating detection capability for all stations combined

statpair = [];
detect_cell ={}; % contains all detection matrices
detect_final={}; 

for stn = 1:g_len_tot
    for i = 1:Ns2 % loop that calculates how many stations can detect the event
        statpair = ratio_diff{i,stn};
        detect = zeros(utm_lenns,utm_lenew,z_len); 
           % detection matrix for 1 station combination
    
    for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len-1
                if  (abs(statpair(j,k,m)) >= thres)
                    detect(j,k,m) = detect(j,k,m) + 1;

                else
                    detect(j,k,m) = 0;
                end
            end
        end
    end
   detect_cell{i,stn} = detect;
    end
end

for stn = 1:g_len_tot
    detects = zeros(utm_lenns,utm_lenew,z_len); % one detection matrix 
for i = 1:Ns2 % loop to perform Boolean sum operation
    for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len-1
                
                    if (detect_cell{i,stn}(j,k,m) == 1) %|| (detect_cell{Ns2+i,stn}(j,k,m) == 1)
                         detects(j,k,m) = detects(j,k,m) + 1;
                   % else
                       % detects(j,k,m) = detects(j,k,m);
                    end
                end
                
            end
    end
   
end
 detect_final{stn} = detects;
end

for stn=1:g_len_tot
    detects  = detect_final{stn};

      for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len

                  if zmesh(j,k,m)>= Vq(j,k)
                    detects(j,k,m) = 0;
                else
                    continue
                end
            end
     end
  end

detect_final{stn} = detects;
end

detect_final_topo ={};

for stn=1:g_len_tot
    detects  = detect_final{stn};

      for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len
                
                  if zmesh(j,k,m)>= (Vq(j,k)-50)
                    detects(j,k,m) = 0;
                else
                    continue
                end
            end
     end
  end

detect_final_topo{stn} = detects;
end


%% Part 8) getting volumes of surface - to define which combination has the best detectability

volume_all = []; % contains volumes
vol_new =[];

 for i = 1:g_len_tot
    detect = detect_final{i};
%     [r,c,z] = ind2sub(size(detect),find(detect>= st_detect)); % r,c,z corresponds to row column and depth positions
%     div = size(r,1); % gives size to divide by for centroid
%     vol = div*abs(dx)*abs(dy)*depth_int ; % calculates volume as m^3
%     vol = ceil(vol/10e9); % round and change into km^3
%     volume_all(i) = vol;
    vol_new(i) = sum(detect_final_topo(:));
 end

for i = 1:g_len_tot % add column with names of stations next too and detect grid
        volume{i,1} = volume_all(i);
        volume{i,2} =  xstat_vec(i);% change to station position
        volume{i,3} =  ystat_vec(i);
        volume{i,4} = detect_final{i};
        volume{i,8} =  detect_final_topo{i};
        volume{i,9} = z_add(i);
        volume{i,10} = vol_new(i);
        volume{i,11} = amp_ratio_all{i};
        volume{i,12} = ratio_diff{i};
end

%% Part 9) Calculating 3D centroid positions, distance of centroid from summit and the standard deviation

 for i = 1:g_len_tot
   detect = volume{i,4}; 
    
    r = []; c = []; z =[];
    [r,c,z] = ind2sub(size(detect),find(detect>= st_detect)); % r,c,z corresponds to row column and depth positions corresponds
    
    v_value_x = []; v_value_y = [];  v_value_z = [];
    for j = 1:size(r,1)
        v_value_x(j) = utmew(r(j),c(j),z(j));
        v_value_y(j) = utmns(r(j),c(j),z(j));
        v_value_z(j) = z_3d(r(j),c(j),z(j));
    end

    % location of centroid
    div = size(r,1);
    cen(i,1) = (sum(v_value_x))/div;
    cen(i,2) = (sum(v_value_y))/div;
    cen(i,3) = (sum(v_value_z))/div;
    volume{i,5} = cen(i,:);
    
    % distance of centroid from summit
    volume{i,6} = sqrt( ((utme_s/1000)-(cen(i,1)/1000))^2 + ((utmn_s/1000)-(cen(i,2)/1000))^2 + ((sum_elev/1000)-(cen(i,3)/1000))^2 );
    
    % standard deviation on centroid
    cen_std(i,1) = std(v_value_x(:));
    cen_std(i,2) = std(v_value_y(:));
    cen_std(i,3) = std(v_value_z(:));
    volume{i,7} = cen_std(i,:);
    
 end

%% Part 10) Convert output data cell to a table

detect_table = cell2table(volume,'VariableNames',{'Volume' 'EW_m loc' 'Ns_m loc'...
    'Detection','Centroid_Location','Distance_centroid_to_summit','SD_of_centroid'...
    'Detection with Topo' 'Elevation of added stations','new volume','amp_ratio_all','ratio_diff'});

save('add_piton.mat','volume','-v7.3');

%% SECTION II) PLOTTING
%
% Plot empty 3d grid
% Plot 2d detection grids for added stations
% Plot isosurfaces

% Not yet edited
% Plot DEM map
% Plot amp ratios
% Plot amp ratio diff
% plot isosurface

%% Part 11) Plot 3d grid with station locations

figure % plotting 3d grid and station locations together
hold on
grid on

for i = 1:z_len % plots gridpoints in utm;
    gr = plot3(utmew(:,:,1)/1000,utmns(:,:,1)/1000,z_3d(:,:,i)/1000,'*','markersize',4,'color',[0.8 0.8 0.8]); 
end

for i = 1:g_len_tot % plots added station grid- in utm
     ads = plot3(xstat_vec(i)/1000,ystat_vec(i)/1000,z_add(i)/1000,'.k','markersize',12,'marker','^','markerfacecolor',[0.9 0.5 0.5]); % in utm
end

for i = 1:Ns % plots station location - in utm
     triangle = plot3(utme_o(i)/1000,utmn_o(i)/1000,elev(i)/1000,'.k','markersize',12,'marker','^','markerfacecolor',[0.9 0.9 0.2]); % in utm
     text(utme_o(i)/1000,utmn_o(i)/1000,elev(i)/1000,sta_name(i,:),'fontsize',11,'fontweight','bold');
end

view(3)
title('3D Grid with Station Locations','fontsize',16,'fontweight','bold')
zlabel('Depth, km','fontsize',14,'fontweight','bold') % creating axis labels
xlabel('Easting, km','fontsize',14,'fontweight','bold')
ylabel('Northing, km','fontsize',14,'fontweight','bold')

leg = legend([gr(1) triangle ads], 'Gridpoints', 'Seismic Stations','Grid of Added Stations');
set(leg,'fontsize',14);
axis equal

%% Part 12) Detection Plot

% Plots detection for each added station. Detection is compared to
% detection calculation of original network (from detect_positions.m)

% Volume and radial distance and direction of new centroid is comapared to the
% original seismic network centroid. All variables except the number of stations
% must be identical between detect_add_stations.m and detect_positions.m

% DETECTION SHOULD ALWAYS INCREASE

%% 12a) calculating variables needed for plot

orig_data = load('detect_piton.mat'); % load in and extract synthetic data from real seismic network
orig_data = (struct2cell(orig_data));
orig_data = orig_data{1,1};
orig_detect = orig_data{1,19};
orig_vol = sum(orig_detect(:));

for i = 1:g_len_tot 
     orig= volume{i,8};
     origsum = sum(orig(:));
     vol_scale(i) = ((origsum-orig_vol)/orig_vol)*100;
end

% extract original location of centroid
orig_cen = orig_data{1,2};

% considers if x, y, z and whether is N,S,E,W, up or down when plotting
% errorbar
for i = 1:g_len_tot
    xdist(i) = (cen(i,1)-orig_cen(1,1))/abs(orig_cen(1,1));
    ydist(i) = (cen(i,2)-orig_cen(1,2))/abs(orig_cen(1,2));
    zdist(i) = (cen(i,3)-orig_cen(1,3))/abs(orig_cen(1,3));
    rad_dist(i) = sqrt((cen(i,1)-orig_cen(1,1))^2 + (cen(i,2)-orig_cen(1,2))^2 + (cen(i,3)-orig_cen(1,3))^2); % in m
   
    xdist_adj(i) = (cen(i,1)-orig_cen(1,1)); % abs distance between centroids (x)
    ydist_opp(i) = (cen(i,2)-orig_cen(1,2)); % (y)
    zdist_hyp(i) = (cen(i,3)-orig_cen(1,3)); % (z)
end

%% 12b) creating color scheme based off of depth of new centroid

% if centroid is originally negative, and centroid moves deeper then
% zdist_hyp will be negative. The deeper the new centroid is, the larger the
% distance will be (i.e. more negative)

%% LIGHTER DEEPER DEPTHS color scheme
% the deeper the new centroid, the lighter the arrow

min_zdist = min(zdist_hyp);
for i = 1:length(zdist_hyp)
    cr_depth(i,1) = (zdist_hyp(i)/min_zdist) ; % deeper will be lighter
    cr_depth(i,2) = (zdist_hyp(i)/min_zdist) ;
    cr_depth(i,3) = (zdist_hyp(i)/min_zdist) ;
%     if cr_depth(i,1) < 0
%         cr_depth(i,1) = abs(cr_depth(i,1));
%         cr_depth(i,2) = abs(cr_depth(i,1));
%         cr_depth(i,3) = abs(cr_depth(i,1)); % above 0m depths will be black
%     %end
%     else if cr_depth(i,1) >= 0
%         cr_depth(i,1) = 0.89;
%         cr_depth(i,2) = 0.89;
%         cr_depth(i,3) = 0.89; % avoids invisble arrows
%         end
%     end
        
end

[cr_depth2, ind]= (sort(cr_depth(:,1),'descend'));

xstat_vec=(xstat_vec(ind));
ystat_vec=(ystat_vec(ind));
z_add=(z_add(ind));
xdist_adj = (xdist_adj(ind));
ydist_opp = (ydist_opp(ind));
zdist_hyp = (zdist_hyp(ind));

cm=colorcet('L17');

cind = round(size(cm,1)/length(ind),1,'significant');

%% 12c) Plotting vertical slice (2D plane view) detection figure

figure
hold on
grid on
view(2)

[m,c] = contour3(utmew(:,:,1)/1000,utmns(:,:,1)/1000,Vq/1000,...
    [1e-10 1e-10],'color',[0.4 0.4 0.4]);%,'ShowText','on',... % elevation contour

for i = 1:g_len_tot % plot arrows and volume
    if ismissing(cr_depth2(i)) == 0 
       qu = quiver3(xstat_vec(i)/1000,ystat_vec(i)/1000,z_add(i)/1000,xdist_adj(i)/1000,ydist_opp(i)/1000,zdist_hyp(i)/1000,'color',cm(i*cind,:),'maxheadsize',1.2,'linewidth',2.5,'autoscalefactor',2);
       ad = scatter3(xstat_vec(i)/1000,ystat_vec(i)/1000,z_add(i)/1000,90,vol_scale(i),'filled','marker','o','Markeredgecolor','k');
    end
end

for i = 1:Ns % plots station location - in utm
     triangle = plot3(utme_o(i)/1000,utmn_o(i)/1000,elev(i)/1000,'.k','markersize',12,'marker','^','markerfacecolor',[0.5 0.8 0.5]); % in utm
     text((utme_o(i)/1000)+0.4,(utmn_o(i)/1000)+0.4,elev(i)/1000,sta_name(i,:),'fontsize',13);
end

% add reference from original data
ori = plot3((orig_cen(1)/1000),orig_cen(2)/1000,orig_cen(3)/1000,'marker','p','markerfacecolor',[0.2 0.9 0.9],'markeredgecolor','k','markersize',16,'linestyle','none');

% plot summit
smit = plot3((utme_s/1000),utmn_s/1000,sum_elev/1000,'marker','p','markerfacecolor',[0.9 0.2 0.9],'markeredgecolor','k','markersize',16,'linestyle','none');

xlabel('Easting, km','fontsize',14,'fontweight','bold');
ylabel('Northing, km','fontsize',14,'fontweight','bold');
zlabel('Depth, km','fontsize',14,'fontweight','bold');
title({'Summary of Detection compared'; 'to Original Network'},'fontsize', 16);

colorcet('L4') %L6 or L4 works well
%caxis([0 max(vol_scale)]);

hcb = colorbar; % volume colorbar
set(hcb,'fontsize',12);
title(hcb,{'% difference in volume'}','fontsize',14);

l = legend([ad triangle qu ori smit],'Added station','Original seismic network','Centroid shift, in degrees','Original centroid location','Summit location');
set(l,'location','westoutside','fontsize',16)

axis equal
ylim([staty0/1000-10 staty00/1000+10]);
xlim([statx0/1000-10 statx00/1000+10]);
set(gca,'fontsize',14);

%% 12d) Creating arrows (Centroid shift) colorbar

figure
hold on
for i = 1:g_len_tot % plot arrows and volume
    qu = quiver3(xstat_vec(i)/1000,ystat_vec(i)/1000,z_add(i)/1000,xdist_adj(i)/1000,ydist_opp(i)/1000,zdist_hyp(i)/1000,'color',cm(i*cind,:),'maxheadsize',1.2,'linewidth',2.5,'autoscalefactor',1);
end

[sort_cr,ind_cr] = sort(cr_depth,'descend'); % make a gradient color scheme

colormap(colorcet('L17');
acb = colorbar;
title(acb, 'Centroid Shift, km','fontsize',14);

tick_cen = ((0:0.1:1)*min_zdist)/1000;
tick_cen = round(tick_cen,2);
acb.TickLabels = num2cell(tick_cen);



% END ---------------------------------------------------------------------
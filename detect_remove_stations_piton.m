%% 
% detect_positions_stations.m
%
% Calculates detectability of different station combinations using a known
% seismic network
%
%
% SECTION 1: CALCULATING VARIABLES, SARA AND OUTPUT
%    Part 1) Loading station location, creating other variables, picking
%            number of stations to look at
%            (contains values/data that user must input)
%    Part 2) Creating grid 
%            (contains values/data that user must input)
%    Part 3) Creating the names of the station pairs ratios
%    Part 4) Calculating ampltitude ratios
%    Part 5) Calculating difference in amplitude ratios
%    Part 6) Calculating detectability
%    Part 7) Calculating volumes
%    Part 8) Calculating centroids and stds of each combination
%    Part 9) Exporting output data as table
%
% SECTION 2: PLOTTING OUTPUT
%    Part 10) Plotting 3d grid of points
%    Part 11) 2D detection grid - with scaled volume and scaled
%             centroid shift(relies on importing data from
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

%% 1b) Extracting station names from input data

% Load in station position data
load Piton_stations.mat % currently as txt data
data = PdF_st_location; % input_data = cell array containing data
Ns = size(data,1); % length of input data/ number of stations

% Extracting longitude, latitudes & elevations in decimal degrees
lat = cell2mat(data(:,2)); % Latitude column of data
lon = cell2mat(data(:,3)); % Longitude column
elev = cell2mat(data(:,6)); % elevations

[utme_o,utmn_o] = deg2utm(lat,lon); % station locations of unedited network

max_elev = 0; % max elevation of volcano,used to close isosurface, should be set to zero
sum_elev = 2632; % elevation of summit of volcano

summit_lon = 55.708889;
summit_lat = -21.2425; % coordinates of summmit of volcano
[utme_s,utmn_s] = deg2utm(summit_lat,summit_lon); % summit in utm


%% 1b) Extracting station names from input data

% Extracting station names from input data
sta_name = char(data(:,1)); % column of data containing station names
sta_char =3; % input the number of characters in the station name (all stations names must have the same no of characters)

%% 1c) Extracting  the number of stations to include in calculations

% each possible combination must be input by user
sta_name2 = nchoosek(1:Ns,2); % getting indexes for each possible station name combination
sta_name3 = nchoosek(1:Ns,3); % 3 stations
sta_name4 = nchoosek(1:Ns,4); % 4 stations
sta_name5 = nchoosek(1:Ns,5); % 5 stations etc.
sta_name6 = nchoosek(1:Ns,6);
sta_name7 = nchoosek(1:Ns,7);
sta_name8 = nchoosek(1:Ns,8);
sta_name9 = nchoosek(1:Ns,9);
sta_name10 = nchoosek(1:Ns,10);

% PICK NUMBER OF STATIONS
sta_name_use = sta_name9; % just swap for number of stations you are interested in

%% 1d) Extracting data associated with chosen stations

Nsr = size(sta_name_use,2); % number of stations now using
Nsrem = Ns-Nsr; % number of stations remaining in the network
comblen = size(sta_name_use,1); % number of possible combinations for sta_name_use number of stations

sta_name_cell = {}; % stations not removed
sta_rem_cell = {}; % names of removed stations
lon_cell = {}; % contains lon of remaining stations
lat_cell = {}; % contains lat of remaining stations
elev_cell = {}; % contains elevations of remaining stations
rem_loc = {}; % contains location of removed stations for each combination

for i = 1:comblen
    
    st_example = sta_name(sta_name_use(i,:),:); % getting each possible combination of station names
    rem_idx = (setxor(sta_name_use(i,:),1:Ns));
    sta_rem =  sta_name(rem_idx,:);
    sta_rem_cell = [sta_rem_cell; sta_rem]; % lists station removed for each combo
    sta_name_cell = [sta_name_cell; st_example]; % putting each station name combination into a cell
    
    for q = 1:Nsrem        
        lon_rem(q) = lon(rem_idx);
        lat_rem(q) = lat(rem_idx);
        elev_rem(q) = elev(rem_idx);
    end
    
    for j = 1:Nsr
        lon_one(j) = lon(sta_name_use(i,j));
        lat_one(j) = lat(sta_name_use(i,j));
        elev_one(j) = elev(sta_name_use(i,j));
    end
    
    lon_cell = [lon_cell; lon_one]; % getting location of stations for each possble combo in a cell
    lat_cell = [lat_cell; lat_one];
    elev_cell =[elev_cell; elev_one];
    
    rem_loc{i,1} = lon_rem;
    rem_loc{i,2} = lat_rem;
    rem_loc{i,3} = elev_rem;
     
end

%% 1e) Calculating positions of chosen stations in UTM and size of DEM map

sta_cell = {}; % al possible station combination locations in UTM (not removed stations)

% convert the not removed station positions from  decimal degrees to utm 
for i = 1:comblen
    
    lat_in = cell2mat(lat_cell(i));
    lon_in = cell2mat(lon_cell(i));
    elev_in = reshape(cell2mat(elev_cell(i)),Nsr,1);
    
    [utme,utmn] = deg2utm(lat_in,lon_in); % external function that calculates UTM from decimal degree
    
    sta = [utme utmn elev_in]; % vector containing utm locations of stations
    sta_cell = [sta_cell; sta];
    
end

% converts the removed station positions from decimal degrees to utm
sta_rem_loc = {}; % contains utm location of removed stations for each station combination

for i = 1:comblen
    lon_use = rem_loc{i,1};
    lat_use = rem_loc{i,2};
    elev_use = rem_loc{i,3};
    
   [utme_rem,utmn_rem] = deg2utm(lat_use,lon_use);
   sta_rem =[utme_rem utmn_rem elev_use];
   sta_rem_loc = [sta_rem_loc; sta_rem];
   
end

z_value = 8000; % Arbitrary elevation value, used for plotting station names on relief map
               % - just has to be larger than any of the actual elevations

vector_latlon = [51.7, 51.9,-177.5,-178.3]; % vector of form 
% [latstart, latend, lonstart, lonend], denotes geographic locations of 
% seismic stations, (is needed to create DEM map - not used in grid
% creation)

%% 1f) Inputs for calulating B (used in amplitude ratio equation)

Q = 170; % quality factor (orig 170)
f  =10; % Mean value for the frequency used  (orig 10)
beta = 1700; % shear wave velocity
n= 1; % for body waves

B= pi*f./(Q*beta); % B = isotropic medium constant

%% 1g) Picking location that is used in amplitude ratio graph

lon_index = 6; % pick index for longitude (CONVERT TO LOCATION)
lat_index = 7; % pick index for latitude
UTMzone = '40 K'; % must define UTM zone this point is in (for amp loc)

%% 1h) Pick direction of migration, grid spacing and threshold

mig_dir = 'z'; % indicate the direction of migration (x, y or z);
gs = 1000; % minimum value is 100 m
start_depth = -10000; % gives the starting depth, max depth is 20 km, must be negative,  min is summit elevation + 1000m.

x_horz_gs = 500; % minimum value is 2000m  
y_horz_gs = 500; % minimum value is 2000m  

% define threshold based off of absolute differencein ampltitude ratios
thres = 0.1;
st_detect = Ns-1; % surface encloses all gridpoints with st_detect and greater number of station pairs, must be integer

%% PART 2) CREATING GRID THAT CAN BE USED FOR DETECTABILITY
% This section contains variables that the user must input

%% 2a) Create 3d grid

detect_grid = []; % creates a grid of points, migration will be from these points and used in the calculation for detectability

% centroid of stations
cen_station = [mean(sta(:,1)) mean(sta(:,2)) mean(sta(:,3))];

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
latdiff = (latend-latstart)/ceil((size(lat_slice,2)/1)-1);
londiff = (lonend-lonstart)/ceil((size(lat_slice,1)/1)-1);

latgrid = [latstart:latdiff:latend];
longrid = [lonstart:londiff:lonend];


% define elevation values in m
out_srtm=[];
for i = 1:length(longrid)
    for j = 1:length(latgrid)
        out = GetSRTMData([longrid(i) latgrid(j)]); % out gives x,y,z for one point
        out_srtm{i,j} = out;
        out_z(i,j) = out(3);
    end
end

% define elevation values in m
depth_int = 50; % depth interval
min_depth = max_elev+(sum_elev+1000); % min depth for grid
depth_idx = gs/depth_int; % alters code so uses user set grid spacing gs

x_idx = x_horz_gs/dy;
y_idx = y_horz_gs/dx;

max_depth = start_depth-gs; % max depth for grid

layer_mat = max_depth:depth_int:min_depth; % gives 1d matrix of elevation values
z_len = length(layer_mat); % length of elevation data

% use meshgrid to create longitude, latitude and depth vectors
[ymesh,xmesh, zmesh] = meshgrid((y0-x_horz_gs):dy:y00,(x0-y_horz_gs):dx:x00,max_depth:depth_int:min_depth); 

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


%% Part 3) Creating the names of the station pairs ratios (e.g. station a / station b)

ratio_data = {}; % contains names of station pairs

for stn = 1:comblen
    st_data = cell2mat(sta_name_cell(stn));
 
    for st1 = 1:Nsr
        for st2 = 1:Nsr
            if st2 > st1
                ratio_pair = strcat(st_data(st2,:),'/',st_data(st1,:));
                ratio_data = [ratio_data; ratio_pair];
              
            end
        end

    end
end

Ns2 = size(nchoosek(1:Nsr,2),1); % length of ratio_data

ratio_data = reshape(ratio_data,Ns2,comblen); % char containing names of the station pairs ratios
% each column corresponds to a different station combination
ratio_data_mat = cell2mat(ratio_data);

%% 4) Calculating amplitude and amplitude ratios
% amplitude is not attenuated - each location has a source amplitude

%% 4a) Calculating amplitude for each station

final_all_d = {};
final_all_amp = {};
src_amp =1; % amplitude of events

% for loop to calculate distance and amplitude between each grid point and station
for clen = 1:comblen
    sta =  sta_cell{clen};
    
    for st = 1:Nsr 

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

amp = reshape(amp,Nsr,comblen); % each column is a diff combo
dist = reshape(dist,Nsr,comblen);
    
    
%% 4b) Calculating amplitude ratios for each station pair

% for loop to calculate amplitude ratios
amp_ratio = {}; % cell containing amplitude ratios

for clen = 1:comblen
    amp_use = amp(:,clen);
    
    for st1 = 1:Nsr 
        for st2 = 1:Nsr
            if st2 > st1 
                  amp_pair = log10(amp_use{st1,:}./amp_use{st2,:}); % amplitude ratio without
                  amp_ratio = [amp_ratio; amp_pair];
            end 
        end
    end   
end

amp_ratio = reshape(amp_ratio, Ns2, comblen);

amp_ratio_all = amp_ratio;
ratio_data_all = ratio_data;

%% 5) Calculating differences between amp ratios for a migration

%% 5a) Calculating amplitude ratio differences

ratio_diff = {}; % gives absolute difference for amp ratios for each station pair
statpair = [];
diff =[];


for stn = 1:comblen
for i =1:Ns2
    statpair = (amp_ratio_all{i,stn});
    
    for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m =1:z_len

                if mig_dir == 'x'
                    if j>x_idx
                    diff(j,k,m) = abs(statpair(j,k,m)-statpair(j-x_idx,k,m)); % x direction , e-w 
                    end
                   
                else if mig_dir == 'y'
                        if k>y_idx
                        diff(j,k,m) = abs(statpair(j,k,m)-statpair(j,k-y_idx,m)); % y direction, n-s
                        end
                        
                else if mig_dir == 'z'
                        if m>depth_idx
                        diff(j,k,m) = abs(statpair(j,k,m)-statpair(j,k,m-depth_idx));
                        end      
                        
               end
        end
    end
               
            end
        end
    end
    ratio_diff{i,stn} =abs(diff); 
end
end

ratio_diff_ind = (ratio_diff); % rename

for i = 1:Ns2
    ratio_name_sum(i,:) = ratio_data_all(i,:);
end



%% 6) Detection capability and detection plot for every station combination
%% 6a) calculating detection capability for all stations combined

% calculating detection
statpair = [];
detect_cell ={}; % cell containing all detection matrices
detect_final ={};

for stn = 1:comblen
    for i = 1:Ns2 % loop that calculates how many stations can detect the event
        statpair = ratio_diff_ind{i,stn};
 
    
    for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len-1
                if  (abs(statpair(j,k,m)) >= thres)
                    detect_cell{i,stn}(j,k,m) = detect_cell{i,stn}(j,k,m) + 1;

                 else
                     detect_cell{i,stn}(j,k,m) = 0;
                end
            end
        end
    end
end

end



for stn = 1:comblen
    detects = zeros(utm_lenns,utm_lenew,z_len);

for i = 1:Ns2 % loop to perform Boolean sum operation
    for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len-1
                
                    if (detect_cell{i,stn}(j,k,m) == 1) 
                         detects(j,k,m) = detects(j,k,m) + 1;

                    end
                end
                
            end
        end
end
detect_final{stn} = detects;
end

detect_final_topo = {};

for stn = 1:comblen
    detects = detect_final{stn};
      for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len
                
                if zmesh(j,k,m)> (Vq(j,k)-50)
                    detects(j,k,m) = 0;
                else
                    continue
                end
            end
        end
      end
detect_final_topo{stn} = detects;
end

%% Part 7) getting volumes of surface - to define which combination has the best detectability

volume_all = []; % cell containing volumes
vol_new =[];

 for i = 1:comblen
    detect = detect_final{i};
%     [r,c,z] = ind2sub(size(detect),find(detect>= st_detect)); % r,c,z corresponds to row column and depth positions
%     div = size(r,1); % gives size to divide by for centroid
%     vol = div*abs(dx)*abs(dy)*depth_int ; % calculates volume as m^3
%     vol = ceil(vol/10e9); % round and change into km^3
%     volume_all(i) = vol;
    vol_new(i) = sum(detect_final_topo(:));
 end


for i = 1:comblen 
    volume{i,1} = volume_all(i);
    volume{i,2} = (sta_name_cell(i)); % remaining stations
    volume{i,3} = sta_rem_cell(i); % stations removed
    volume{i,4} = sta_rem_loc(i,:); % location of removed station
    volume{i,5} = (detect_final_topo{i});
    volume{i,9}= detect_final{i};
    volume{i,10} = vol_new(i);
    volume{i,11} = amp_ratio_all{i};
    volume{i,12} = ratio_diff{i};
end

%% Part 8) Calculating 3D centroid location

 for i = 1:comblen
    detect = volume{i,5}; % extracts volumes
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
    cen(i,1) = (sum(v_value_x(:)))/div;
    cen(i,2) = (sum(v_value_y(:)))/div;
    cen(i,3) = (sum(v_value_z(:)))/div;
    volume{i,6} = cen(i,:);
    
    % distance of centroid from summit
    volume{i,7} = sqrt( ((utme_s/1000)-(cen(i,1)/1000))^2 + ((utmn_s/1000)-(cen(i,2)/1000))^2 +...
    ((sum_elev/1000)-(cen(i,3)/1000))^2 );
    
    % standard deviation on centroid
    cen_std(i,1) = std(v_value_x(:));
    cen_std(i,2) = std(v_value_y(:));
    cen_std(i,3) = std(v_value_z(:));
    volume{i,8} = cen_std(i,:);
    
 end

%% Part 9) convert output data to table

detect_table = cell2table(volume,'VariableNames',{'Volume' 'Stations_not_removed'...
    ' Stations_removed'  'Loc_of_removed_stations' 'Detection_topo_removed' 'Centroid_Location'...
    'Distance_centroid_to_summit' 'SD_of_centroid' 'Detection','new volume','amp_ratio_all','ratio_diff'});

save('remove_piton.mat','volume','-v7.3');

%% SECTION II) PLOTTING
%
% Plot empty 3d grid
% Plot DEM map
% Plot amp ratios
% Plot amp ratio diff
% plot isosurface

%% Part 10) plot 3d grid

figure % plotting 3d grid and station locations together
hold on
grid on

for i = 1:z_len % plots gridpoints in utm;
    gr = plot3(utmew(:,:,1)/1000,utmns(:,:,1)/1000,z_3d(:,:,i)/1000,'*','markersize',4,'color',[0.8 0.8 0.8]); 
end

for i = 1:Ns % plots station location - in utm
     triangle = plot3(utme_o(i)/1000,utmn_o(i)/1000,elev(i)/1000,'.k','markersize',12,'marker','^','markerfacecolor',[0.5 0.9 0.5]); % in utm
     text(utme_o(i)/1000,utmn_o(i)/1000,elev(i)/1000,sta_name(i,:),'fontsize',10);
end

view(3)
title('3D grid with station locations')
zlabel('Depth, km') % creating axis labels
xlabel('E-W, km')
ylabel('N-S, km')
legend([gr(1) triangle], 'Gridpoints', 'Seismic Stations');

%% Part 11) 2D plot showing the effect of removing just 1 station from network
% Plots detection when ONE station is removed - mute if removing more than one. Detection is compared to
% detection calculation of original network (from detect_positions.m)

% Volume and radial distance of centroid from the detectability volume and
% radial direction of the same centroid are comapared to data from the
% original seismic network. All variables except the number of stations
% must be identical between detect_removed_stations.m and detect_positions.m.

% DETECTION SHOULD ALWAYS DECREASE

%% 11a) Calculating variables needed for 2d plot

% load in and extract synthetic data from real seismic network
orig_data = load('detect_piton.mat'); 
orig_data = (struct2cell(orig_data));
orig_data = orig_data{1,1};
orig_sd = orig_data{1,4}; % orig sd
orig_sd_x = orig_sd(1);
orig_sd_y = orig_sd(2);
orig_sd_z = orig_sd(3);
orig_cen = orig_data{1,2}; % original location of centroid
orig_detect = orig_data{1,19};
orig_vol = sum(orig_detect(:));

for i = 1:comblen
      orig= volume{i,5};
      vol_scale(i) = ((sum(orig(:))-orig_vol)/orig_vol)*100;
end

% considers if x, y, z and whether is N,S,E,W, up or down when plotting
% errorbar
for i = 1:comblen
    xdist(i) = (cen(i,1)-orig_cen(1,1))/abs(orig_cen(1,1));
    ydist(i) = (cen(i,2)-orig_cen(1,2))/abs(orig_cen(1,2));
    zdist(i) = (cen(i,3)-orig_cen(1,3))/abs(orig_cen(1,3));
    rad_dist(i) = sqrt((cen(i,1)-orig_cen(1,1))^2 + (cen(i,2)-orig_cen(1,2))^2 + (cen(i,3)-orig_cen(1,3))^2); % in m
    
    xdist_adj(i) = (cen(i,1)-orig_cen(1,1));% abs distance between centroids (x)
    ydist_opp(i) = (cen(i,2)-orig_cen(1,2)); % (y)
    zdist_hyp(i) = (cen(i,3)-orig_cen(1,3)); % (z)
    mag(i) = sqrt(xdist_adj(i)^2 + ydist_opp(i)^2 +  zdist_hyp(i)^2);
    
end

%% 11b) Create color scheme for plots

%% Lighter arrows for larger depths 

max_zdist = max(zdist_hyp);
min_zdist = min(zdist_hyp);
for i = 1:length(zdist_hyp)
    cr_depth(i,1) = (zdist_hyp(i)/max_zdist) ; % deeper will be lighter
    cr_depth(i,2) = (zdist_hyp(i)/max_zdist) ;
    cr_depth(i,3) = (zdist_hyp(i)/max_zdist) ;
%     if cr_depth(i,1) > 0
%         cr_depth(i,1) = 0;
%         cr_depth(i,2) = 0;
%         cr_depth(i,3) = 0; % above 0m depths will be black i.e if moved upwards
% %     end
%      if cr_depth(i,1) >= 0.92
%         cr_depth(i,1) = 0.91;
%         cr_depth(i,2) = 0.91;
%         cr_depth(i,3) = 0.91; % avoids invisible arrows
%     end
        
end

[cr_depth2, ind]= (sort(cr_depth(:,1),'descend'));

sta_rem_loc = cell2mat(sta_rem_loc);
xdist_adj = (xdist_adj(ind));
ydist_opp = (ydist_opp(ind));
zdist_hyp = (zdist_hyp(ind));

cm = colorcet('L17');

cind = round(size(cm,1)/length(ind),1,'significant')-10;

%% 11c) Plotting vertical slice (2D plane view) detection figure

figure
hold on
grid on
axis square
view(2)

[m,c] = contour3(utmew(:,:,1)/1000,utmns(:,:,1)/1000,Vq/1000,...
    [1e-10 1e-10],'color',[0.4 0.4 0.4]);%,'ShowText','on',... % elevation contour

 cc = distinguishable_colors(comblen);

for i = 1:comblen % plot arrows and volume    ;
    st_rem = sta_rem_loc(i,:);
    qu(i) = quiver3(st_rem(1)/1000,st_rem(2)/1000,st_rem(3)/1000,xdist_adj(i)/1000,ydist_opp(i)/1000,zdist_hyp(i)/1000,'color',cm(i*cind,:),'maxheadsize',1.2,'linewidth',2.5,'autoscalefactor',1.5);
    ad(i) = scatter3(st_rem(1)/1000,st_rem(2)/1000,st_rem(3)/1000,140,vol_scale(i),'filled','marker','^','Markeredgecolor','k');
    text((st_rem(1)/1000)+0.2, st_rem(2)/1000+0.1,st_rem(3)/1000+0.1,sta_rem_cell{i,:},'fontsize',15,'fontweight','bold');

    % plots surface contour - mute if not wanted
%     stname(i,:) = cell2mat(volume{i,3});
%     detect_grid = (volume{i,5}); 
%     p = patch(isosurface(utmew/1000,utmns/1000,z_3d/1000, detect_grid, st_detect));
%     %set(p(i), 'FaceColor', [0.3 0.3 0.3] , 'EdgeColor','none','FaceAlpha',0.05);
%     p2 = patch(isosurface(utmew/1000,utmns/1000,z_3d/1000, orig_detect, st_detect));
%     set(p2, 'FaceColor', [0.3 0.3 0.3] , 'EdgeColor','none','FaceAlpha',0.05)
%     set(p, 'FaceColor', cc(i,:) , 'EdgeColor','none','FaceAlpha',0.05);
%     view(3)

end

detect_grid = (volume{1,5}); 
p1 = patch(isosurface(utmew/1000,utmns/1000,z_3d/1000, detect_grid, st_detect));
p2 = patch(isosurface(utmew/1000,utmns/1000,z_3d/1000, orig_detect, st_detect));
set(p2, 'FaceColor', [0.3 0.3 0.3] , 'EdgeColor','none','FaceAlpha',0.05);
set(p1, 'FaceColor', [0.47 0.67 0.19] , 'EdgeColor','none','FaceAlpha',0.1);

% add reference from original data
ori = plot3((orig_cen(1)/1000),orig_cen(2)/1000,orig_cen(3)/1000,'marker','p','markerfacecolor',[0.2 0.9 0.9],'markeredgecolor','k','markersize',18,'linestyle','none');

% plot summit
smit = plot3((utme_s/1000),utmn_s/1000,sum_elev/1000,'marker','p','markerfacecolor',[0.9 0.2 0.9],'markeredgecolor','k','markersize',18,'linestyle','none');

xlabel('Easting, km','fontsize',14,'fontweight','bold');
ylabel('Northing, km','fontsize',14,'fontweight','bold');
zlabel('Depth, km','fontsize',14,'fontweight','bold');
title('Summary of Detection compared to Original Network','fontsize', 16);

colormap(colorcet('L4','reverse',true))

hcb = colorbar; % colorbar for volume
set(hcb,'fontsize',12);
title(hcb,{'% difference in volume'}','fontsize',14);

l = legend([ad(1) qu(1) ori smit],'Removed station','Centroid shift, in degrees','Original centroid location','Summit location');
set(l,'location','westoutside','fontsize',16)

axis equal

xlim([(min(utme_o)/1000)-1 (max(utme_o)/1000)+1]);
ylim([(min(utmn_o)/1000)-1 (max(utmn_o)/1000)+1]);
set(gca,'fontsize',14);

%% 11d) Creating arrows (Centroid shift) colorbar

figure
hold on
for i = 1:comblen % plot arrows and volume
    qu = quiver3(st_rem(1)/1000,st_rem(2)/1000,st_rem(3)/1000,xdist_adj(i)/1000,ydist_opp(i)/1000,zdist_hyp(i)/1000,'color',cm(i*cind,:),'maxheadsize',1.2,'linewidth',2.5,'autoscalefactor',1);
end

[sort_cr,ind_cr] = sort(cr_depth,'ascend'); % make a gradient color scheme
colormap(colorcet('L17');
acb = colorbar;
set(acb,'fontsize',12);
title(acb, 'Centroid Shift, km','fontsize',14,'fontweight','bold');

tick_cen = ((0:0.1:1)*max_zdist)/1000;
tick_cen = round(tick_cen,2);
acb.TickLabels = num2cell(tick_cen);


% END ---------------------------------------------------------------------
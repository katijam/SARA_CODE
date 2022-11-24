%% 
% detect_positions.m
%
% script that gives idea of detection capability of seismic network for
% various migration scenarios
%
% this codes define the detection surface in terms of the number of station
% pairs
%
% THIS CODE USES MIGRATION FROM EVERY POSITION WITHIN A GRID (NOT RELATED
% TO THE EVENTS)
% No migration is needed
%
% If the difference between the ratios for a specified window is small (without noise) -
% then there is a low detection capability, large difference indicates a
% good detectability.
%
% Part 1) Loading station data and creating variables
%         (contains values/data that user must input)
% Part 2) Creating grid of points
% Part 3) DEM map plot
% Part 4) Calculating amplitude and amp ratios
% Part 5) Calculating differences between amp ratios for a migration
% Part 6) Finding detection capability and detection plot (points and surface) for EVERY station combination
% Part 7) Calculating volume, centroid and std on centroid of detectability
%         surface. Also saves output data as 'orig_data.mat'
% Part 9) Isosurface (detection surface) plots (with slices)
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


max_elev = 0; % max elevation of volcano,used to close isosurface, should be set to zero
sum_elev = 2632; % elevation of summit of volcano

summit_lon = 55.708889;
summit_lat = -21.2425; % coordinates of summmit of volcano
[utme_s,utmn_s] = deg2utm(summit_lat,summit_lon); % summit in utm


%% 1b) Extracting station names from input data

% Extracting station names from input data
sta_name = (data(:,1)); % column of data containing station names
sta_char =3; % input the number of characters in the station name (all stations names must have the same no of characters)


%% 1c) Calculating positions of stations in UTM and size of DEM map

% convert station positions from  decimal degrees to utm 
[utme,utmn] = deg2utm(lat,lon); % external function that calculates UTM from decimal degree
elev = reshape(elev,Ns,1);
sta = [utme utmn elev]; % vector containing utm locations of stations

z_value = 8000; % Arbitrary elevation value, used for plotting station names on relief map
               % - just has to be larger than any of the actual elevations

vector_latlon = [-21.15 -21.4 55.65 55.8]; 
% [latstart, latend, lonstart, lonend], denotes geographic locations of 
% seismic stations, (is needed to create DEM map - not used in grid
% creation)

%% 1d) Creating the names of the station pairs ratios (e.g. station a / station b)

ratio_data = {}; % contains the name of the station pairs

for st1 = 1:Ns
    for st2 = 1:Ns
        if st2 > st1
            ratio_pair = strcat(sta_name(st2,:),'/',sta_name(st1,:));
            ratio_data = [ratio_data; ratio_pair];
        end
    end
end

ratio_data = cell2mat(ratio_data); % char containing names of the station pairs ratios
Ns2 = size(ratio_data,1); % number of station pair ratios

%% 1e) Inputs for calulating B (used in amplitude ratio equation)

Q = 170; % quality factor (orig 170)
f  = 10; % Mean value for the frequency used  (orig 10)
beta = 1700; % shear wave velocity
n= 1; % for body waves

B= pi*f./(Q*beta); % B = isotropic medium constant

%% 1f) Picking location that is used in amplitude ratio graph

lon_index = 2; % pick index for longitude 
lat_index = 2; % pick index for latitude 
UTMzone = '40 K'; % must define UTM zone this point is in (for amp loc), in terms 'no.no. letter'

%% 1g) Pick direction of migration, grid spacing and threshold

mig_dir = 'z'; % indicate the direction of migration (x, y or z);
gs = 1000; % minimum value is 100 m
start_depth = -10000; % gives the starting depth
x_horz_gs = 500; 
y_horz_gs = 500; 

% Define a detectability threshold and number of station pairs needed for
% event to be detectable and interstaion distance
thres = 0.1; % define threshold based off of absolute difference in amplitude ratios - FIXED
st_detect = Ns; % surface encloses all gridpoints with st_detect and greater number of station pairs, must be integer

%% PART 2) CREATING GRID THAT CAN BE USED FOR DETECTABILITY
% This section contains variables that the user must input

%% 2a) Create 3d grid

% centroid of stations
cen_station = [mean(sta(:,1)) mean(sta(:,2)) mean(sta(:,3))];
detect_grid = []; % creates a grid of points, migration will be from these points and used in the calculation for detectability

% Creation of grid uses latitude and longitude. USER MUST INPUT LAT AND LON START AND END
% VALUES FOR GRID IN UTM

% utm (m) for latitude and longititude
x0 = cen_station(2)-20e3; % latitude inputs (start), 20 is standard now
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
h = waitbar(0,'Please wait...');
for i = 1:length(longrid)
    for j = 1:length(latgrid)
        
        waitbar(i/length(longrid),h,['Progess = ',num2str(i), '/', num2str(longrid)]);
        
        out = GetSRTMData([longrid(i) latgrid(j)]); % out gives x,y,z for one point
        out_z(i,j) =  out(3);
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


for i = 1:length(outlon)
    for j = 1:length(outlat)
        outlon_all(i,j) = outlon(i);
        outlat_all(i,j) = outlat(j);
    end
end

close(h);

Vq = interp2(outlat_all,outlon_all,out_z,utmns(:,:,1),utmew(:,:,1),'spline'); %Vq is interpolated topography data

%% 2b) Finding location of point used in seismic amplitude ratio graph

% UTM location of location used in amplitude ratio plot
point_lon = utmew(1,lon_index,1); % longitude value for location in utm
point_lat = utmns(lat_index,1,1); % latitude value for location in utm

% conversion of gridpoints from utm to decimal degrees (used only for DEM
% map)
[point_lat_deg, point_lon_deg] = utm2deg(point_lon, point_lat, UTMzone);

%% 2c) Plot 3d grid of points

figure; % plotting 3d grid and station locations together
hold on
grid on

for i = 1:z_len % plots gridpoints in utm;
gr = plot3(utmew(:,:,1)/1000,utmns(:,:,1)/1000,z_3d(:,:,i)/1000,'*','markersize',2,'color',[0.8 0.8 0.8]); 
end

for i = 1:Ns % plots station location - in utm
     triangle = plot3(utme(i)/1000,utmn(i)/1000,elev(i)/1000,'.k','markersize',12,'marker','^','markerfacecolor',[[0.5 0.8 0.5]]); % in utm
     text(((utme(i)/1000)+0.7),utmn(i)/1000,elev(i)/1000,sta_name(i,:),'fontsize',11,'fontweight','bold');
end

view(3)
zlabel('Depth, km','fontsize',14)%,'fontweight','bold') % creating axis labels
xlabel('Easting, km','fontsize',14)%,'fontweight','bold')
ylabel('Northing, km','fontsize',14)%,'fontweight','bold')

xlim([y0/1000 y00/1000]);
ylim([x0/1000 x00/1000]);
axis equal
zlim([max_depth/1000 (min_depth/1000)+2]);

set(gca,'fontsize',14);
axesLabelsAlign3D()
%set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])

%% 4) Calculating amplitude and amplitude ratios
% Amplitude is not attenuated - each location has a source amplitude

%% 4a) Calculating amplitude for each station

final_all_d = {};
final_all_amp = {};
src_amp =1; % amplitude of events

% for loop to calculate distance and amplitude between each grid point and station
for st = 1:Ns 
    
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

dist = final_all_d; % cell containing distances between eahc grid point and the staion
amp = final_all_amp; % cell containing amplitude for each grid point each station

%% 4b) Calculating amplitude ratios for each station pair

% for loop to calculate amplitude ratios
amp_ratio = {}; % cell containing amplitude ratios
for st1 = 1:Ns 
    for st2 = 1:Ns
        if st2 > st1     
              amp_pair = log10(amp{:,st1}./amp{:,st2}); % amplitude ratio
              amp_ratio = [amp_ratio; amp_pair]; % amplitude ratios in cell form
              
        end 
    end
end

% combining together to consider every amplitude ratio possibility
amp_ratio_all = {}; % cell for amplitude ratio values for each cell for all st pair combos

amp_ratio_all = amp_ratio;
ratio_data_all = ratio_data;

%% 4c) Plotting all possible amplitude ratio combination for each station pair

% this amplitude graph uses the location defined by point_lon and point_lat
% to create a figure showing a 2d slice of amplitude ratios

figure('Renderer', 'Painters')
ss = []; % matrices for plotting
e=[];
cc = distinguishable_colors(2*Ns2); % creating colors for lines

for i = 1:Ns2 % for loop to plot amp ratios for each station pair
     amp_ratio_mat = amp_ratio_all{i};  
     
    xp=[];
    yp=[];
    
         for j = 1:z_len             
             clear h
             xp(j) = layer_mat(j);
             yp(j) = abs(amp_ratio_mat(lon_index,lat_index,j));             
         end
         
     h = semilogy(xp/1000,yp); % semi log y graph
     set(h,'color',cc(i,:),'marker','.','markersize',10,'linewidth',0.75,'linestyle','-')
     hold on
     grid on
     e = [e; h];
     ss = [ss; i];
    
end
hold on

% [leg,icons,p,txt2] = legend(e,ratio_data_all(ss,:),'FontSize',5);
% set(icons(length(ss)+2:2:end),'markersize',14);
% set(leg,'location','eastoutside')
% title(leg,'Station Pair','Visible','on','fontsize',12,'fontweight','bold');
% txt = ['Amplitude Ratio Slice at Longitude: ', num2str((point_lon)/1000),' km & Latitude: ', num2str((point_lat)/1000), ' km'];
% title(txt,'fontsize',16,'fontweight','bold') 
% ylabel('Amplitude Ratio, Arbitrary Unit','fontsize',14,'fontweight','bold') 
% xlabel('Depth, km','fontsize',14,'fontweight','bold')
% title('LAR');

ylabel('LAR','fontsize',14)%,'fontweight','bold') 
xlabel('Depth, km','fontsize',14)%,'fontweight','bold')
xlim([-10 2.6])

set(gca,'fontsize',14);

%% 5) Calculating differences between amp ratios for a migration

%% 5a) Calculating amplitude ratio differences

ratio_diff = {}; % gives absolute difference for amp ratios for each station pair
statpair = [];
diff =[];

ratio_diff_per = {}; 
per_change = []; % gives percentage change for amp ratios for each station pair

for i =1:Ns2
    statpair = amp_ratio_all{i};
    
    for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m =1:z_len

                if mig_dir == 'x'
                    if j>x_idx
                    diff(j,k,m) = statpair(j,k,m)-statpair(j-x_idx,k,(m));
                    end
                   
                else if mig_dir == 'y'
                         if k>y_idx
                        diff(j,k,m) = statpair(j,k,m)-statpair(j,k-y_idx,(m));
                        end
                        
                else if mig_dir == 'z'
                        if m >depth_idx
                        diff(j,k,m) = statpair(j,k,m)-statpair(j,k,(m-depth_idx));
                        end      
                        
                end
                end
               
            end
        end
        end
    end
    ratio_diff{i} = abs(diff); 
end

ratio_diff_ind = (ratio_diff); % rename
    
%creating new station name array
for i = 1:Ns2
    ratio_name_sum(i,:) = ratio_data_all(i,:);
end

%% 5b) plotting ratio diff for a specfic point

cc = distinguishable_colors(Ns2);

%this amplitude graph uses point_lon and point_lat to create a figure of
%the amplitude ratio difference
figure
%subplot(2,2,3);
ss = [];
e=[];

for i = 1:Ns2 % for loop to plot amp ratios for each station pair
       
     ratio_diff_mat = ratio_diff_ind{i}; 
     xp=[];
     yp=[];
    
         for j = 1:z_len             
             clear h
             xp(j) = layer_mat(j);
             yp(j) = abs(ratio_diff_mat(lon_index,lat_index,j));                   
         end
              
     h = semilogy(xp/1000,yp);
     set(h,'color',cc(i,:),'marker','.','markersize',10,'linewidth',0.75,'linestyle','-')
     hold on
     grid on
     e = [e; h];
     ss = [ss; i];
    
end

tline = yline(thres,'--','linewidth',3.5,'color',[1 0 0]);

% [leg,icons,p,txt2] = legend([e],ratio_name_sum(ss,:),'FontSize',10);
% set(icons(length(ss)+2:2:end),'markersize',14)
% set(leg,'location','eastoutside')  
% title(leg,'Station Pair','Visible','on','fontsize',14,'fontweight','bold');
% txt = ['Amplitude Ratio Difference Slice at Longitude: ', num2str((point_lon)/1000),' km and Latitude: ', num2str((point_lat)/1000),' km'];
% title(txt,'fontsize',16,'fontweight','bold') 
% ylabel('Amplitude Ratio Difference, Arbitrary Unit','fontsize',14,'fontweight','bold') 
% xlabel('Depth, km','fontsize',14,'fontweight','bold')
% title('LAR Difference');

ylabel('LAR Difference','fontsize',14)%,'fontweight','bold') 
xlabel('Depth, km','fontsize',14)%,'fontweight','bold')

set(gca,'fontsize',14);
xlim([-10 2.6])

%% 6) Detection capability and detection plot for every station combination
%% 6a) calculating detection capability for all stations combined

% calculating detection
detects = zeros(utm_lenns,utm_lenew,z_len); % contains final detection values
statpair = [];
detect_cell ={};

for i = 1:Ns2 % reduce number of pairs back to Ns2, for pairs 0-Ns2
    statpair = ratio_diff_ind{i};
    detect_ind = zeros(utm_lenns,utm_lenew,z_len);
        
    for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len-1
                if  (abs(statpair(j,k,m)) >= thres)
                    detect_ind(j,k,m) =  1;

                else
                    detect_ind(j,k,m) = 0;
            end
        end
    end
    end
    detect_cell = [detect_cell; detect_ind];
end
        
for i = 1:Ns2 % loop to calculate detection matrix
    for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len-1
                
                    if (detect_cell{i}(j,k,m) == 1)
                         detects(j,k,m) = detects(j,k,m) + 1;

                    end
                end
                
            end
        end
end


detects_final = detects;

%subtract elevations
for i=1:Ns2
      for j = 1:utm_lenns
        for k = 1:utm_lenew
            for m = 1:z_len
                
                if zmesh(j,k,m)> (Vq(j,k)-50)
                    detects_final(j,k,m) = 0;
                else
                    continue
                end
            end
     end
  end
end

%% 7) Calculating volume, centroid and std on centroid

orig_data = {}; % cell containing data output 

vol = []; % volume of detection
% [r,c,z] = ind2sub(size(detects_final),find(detects_final>= st_detect)); % r,c,z corresponds to row column and depth positions
% div = size(r,1); % gives size to divide by for centroid
% vol = div*abs(dx)*abs(dy)*depth_int ; % calculates volume as m^3
% vol = ceil(vol/10e9); % round and change into km^3
orig_data{1} = vol;

for j = 1:div
    v_value_x(j) = utmew(r(j),c(j),z(j));
    v_value_y(j) = utmns(r(j),c(j),z(j));
    v_value_z(j) = z_3d(r(j),c(j),z(j));
end
    
% location of centroid
cen(1) = (sum(v_value_x(:)))/div;
cen(2) = (sum(v_value_y(:)))/div;
cen(3) = (sum(v_value_z(:)))/div;
orig_data{2} = cen;

% distance of centroid from summit
orig_data{3} = sqrt( (utme_s/1000)-cen(1)^2 + ((utmn_s/1000)-cen(2))^2 + ((sum_elev/1000)-cen(3))^2 );

% standard deviation on centroid
cen_std(1) = std(v_value_x(:));
cen_std(2) = std(v_value_y(:));
cen_std(3) = std(v_value_z(:));

orig_data{4} = cen_std;
orig_data{5} = detects;
orig_data{6} = st_detect;

orig_data{7} = utmew;
orig_data{8} = utmns;
orig_data{9} = z_3d;
orig_data{10} = Vq;
orig_data{11} = detect_cell;
orig_data{12} = Ns2;
orig_data{13} = amp_ratio_all;
orig_data{14} = ratio_diff_ind;
orig_data{15} = layer_mat;
orig_data{16} =outlon_all;
orig_data{17} = outlat_all;
orig_data{18} = out_z;
orig_data{19} = detects_final;

orig_data{2,1} = 'Volume';
orig_data{2,2} = 'Centroid';
orig_data{2,3} = 'Centroid to Summit';
orig_data{2,4} = 'STD centroid';
orig_data{2,5} = 'Detection Coordinates';
orig_data{2,6} = 'Station Threshold';

orig_data{2,7} = 'x coords';
orig_data{2,8} = 'y coords';
orig_data{2,9} = 'z coords';
orig_data{2,10} = 'Topography';
orig_data{2,11} = 'detection of all station pairs';
orig_data{2,12} = 'number of station pairs';
orig_data{2,13} = 'Amplitude Ratios';
orig_data{2,14} = 'Amp Ratio Diff';
orig_data{2,15} = 'Layer Vector';
orig_data{2,16} = 'longrid';
orig_data{2,17} = 'latgrid';
orig_data{2,18} = 'elevation';
orig_data{2,19} = 'detection elevation removed';

save('detect_piton.mat', 'orig_data','-v7.3'); % saves orig_data as a mat file


%% Part 8) Plotting detection capability for visualisation of capability of network in 3d

%% 8a) Adding isosurface (detection surface) to detection grid 

figure('renderer','opengl')
hold on
grid on

p = patch(isosurface(utmew/1000,utmns/1000,z_3d/1000, detects_final, st_detect));

view(3)
set(p, 'FaceColor', [0.45 0.55 0.9], 'EdgeColor','none','FaceAlpha',0.9);

camlight left; 
lighting gouraud
caxis([st_detect-1 Ns2])
txt = ['DCV'] ;

%additional stuff for plot if plot is in separate figure to detection grid

for i = 1:Ns % plots station location - in utm
     triangle = plot3(utme(i)/1000,utmn(i)/1000,elev(i)/1000,'.k','markersize',12,'marker','^','markerfacecolor',[0.5 0.8 0.5]); % in utm
     text(((utme(i)/1000)+0.5),utmn(i)/1000,elev(i)/1000,sta_name(i,:),'fontsize',8,...
        'fontweight','bold','backgroundcolor','w','margin',0.5)
end

%smit = plot3((utme_s/1000),utmn_s/1000,sum_elev/1000,'marker','p','markerfacecolor',[0.9 0.2 0.9],'markeredgecolor','k','markersize',16,'linestyle','none');

zlabel('Depth, km','fontsize',14)%,'fontweight','bold') % creating axis labels
xlabel('Easting, km','fontsize',14)%,'fontweight','bold')
ylabel('Northing, km','fontsize',14)%,'fontweight','bold')

axis equal
xlim([min(p.XData(:))-5 max(p.XData(:))+5 ]);
ylim([min(p.YData(:))-5 max(p.YData(:))+5 ]);
zlim([(min(p.ZData(:)))-3 (max(p.ZData(:)))+3]);
view(3)

% add on topography contour if desired
%[m,c] = contour3(utmew(:,:,1)/1000,utmns(:,:,1)/1000,Vq/1000,[0.5, 1,1.5 2]);%,'ShowText','on',...
    %'labelspacing',180);
% colormap(parula);
% caxis([0 2])
% set(c,'linewidth',1);

lgd = legend([p triangle],txt, 'Station Locations');
set(lgd,'fontsize',20,'location','eastoutside');


set(gca,'fontsize',14);
axesLabelsAlign3D()
%set(gca,'xtick',[],'ytick',[])
set(gca,'xticklabel',[],'yticklabel',[])

a=get(gca);
a.axes.SortMethod='ChildOrder';

% END ---------------------------------------------------------------------

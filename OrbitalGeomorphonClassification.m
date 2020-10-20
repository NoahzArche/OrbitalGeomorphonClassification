% Geomorphon classification for orbital DEMs (elevation in m)

% DEM in TIFF format
t = Tiff(filename);
dem = read(t);
filename = 'Orbital DEM';


%% Configuration Parameters 

res = 1;                        % resolution in meters (HiRISE orbital DEMs 1m)

% Geomorphon parameters in cells
window_size = 701;
skip = 80;
tdegree = 1.5;                  % flatness threshold in degrees
t = tdegree * pi/180;           % convert threshold to radians
scanning_size = size(dem,1);  



%%

% The DEL w/e/n/s/nw/ne/sw/se matrices contain -1's, 0's, and 1's
DELTA_W_matrix = nan(size(dem)); 
DELTA_E_matrix = nan(size(DELTA_W_matrix));
DELTA_N_matrix = nan(size(DELTA_W_matrix));
DELTA_S_matrix = nan(size(DELTA_W_matrix));
DELTA_NW_matrix = nan(size(DELTA_W_matrix));
DELTA_NE_matrix = nan(size(DELTA_W_matrix));
DELTA_SW_matrix = nan(size(DELTA_W_matrix));
DELTA_SE_matrix = nan(size(DELTA_W_matrix));

% index of neighboring cells which are taken into account (ignoring skipped cells)
i = 1 : (window_size-1)/2 * (1/skip);                
% distance from central cell to all the neighboring cells in meters
distance(i) = res * skip * i;                  
distance_diag = sqrt(2) * (res * skip * i); 


%%  Ternary Pattern Computation

% Scanning through each cell in the DEM 
for row = (window_size + 1)/2 : size(dem,1) - (window_size - 1)/2   
    for col = (window_size + 1)/2 : size(dem,2) - (window_size - 1)/2
                          
        % Elevations (in m) of all neighboring cells relative to central cell
        east_relative_elevation = double(dem( row, (col + skip) : skip :...
            col + (window_size - 1)/2 )) - double(dem(row,col));    
        west_relative_elevation = double(dem( row, (col - skip) : -skip :...
            col - (window_size - 1)/2 )) - double(dem(row,col));    
        north_relative_elevation = double(dem( (row - skip) : -skip :...
            row - (window_size - 1)/2, col )) - double(dem(row,col));
        south_relative_elevation = double(dem( (row + skip) : skip: row +...
            (window_size - 1)/2, col )) - double(dem(row,col));               
        nw_relative_elevation = double(dem( (col - 1 - skip) * scanning_size +...
            row - skip :...
            - (scanning_size + 1) * skip : (col - (window_size-1)/2 - 1) *...
            scanning_size )) - double(dem(row,col));        
        sw_relative_elevation = double(dem( (col - 1 - skip) * scanning_size +...
            row + skip : - (scanning_size - 1) * skip :...
            (col - (window_size-1)/2 -1) * scanning_size )) - double(dem(row,col));      
        ne_relative_elevation = double(dem( (col + skip -1) * scanning_size +...
            (row - skip) : (scanning_size - 1) * skip :...
            (col + (window_size-1)/2) * scanning_size )) - double(dem(row,col));
        se_relative_elevation = double(dem( (col + skip -1) * scanning_size +...
            (row + skip) : (scanning_size + 1) * skip :...
            (col + (window_size-1)/2) * scanning_size )) - double(dem(row,col));
        
        % Elevation angles in each lookup direction
        east_angles = atan( east_relative_elevation./ distance );
        west_angles = atan( west_relative_elevation./ distance );
        north_angles = atan( north_relative_elevation'./ distance );
        south_angles = atan( south_relative_elevation'./ distance );        
        nw_angles = atan( double(nw_relative_elevation)./ distance_diag );
        se_angles = atan( double(se_relative_elevation)./ distance_diag );
        sw_angles = atan( double(sw_relative_elevation)./ distance_diag );
        ne_angles = atan( double(ne_relative_elevation)./ distance_diag );
        
        
        
        % Assigning ternary pattern values
        % Starting with EAST direction (lookup direction 1)
        % beta and del denote the maximum and minimum elevation angles, respectively 
        beta_e_value = max(east_angles);             
        del_e_value = min(east_angles);         
        % Distinction of cases for each direction, resulting in the ternary
        % pattern value stored in the DEL matrices
        if (del_e_value + beta_e_value) > t
            
            DELTA_E_matrix(row,col) = 1;  
            
        elseif abs(del_e_value + beta_e_value) < t
            
            DELTA_E_matrix(row,col) = 0;  
            
        else    
            
            DELTA_E_matrix(row,col) = -1;
            
        end
        
        
        % This process is repeated for all the other lookup directions
               
        
        % WEST direction (lookup direction 2)                     
        beta_w_value = max(west_angles);           
        del_w_value = min(west_angles);         
        if (del_w_value + beta_w_value) > t
            
            DELTA_W_matrix(row,col) = 1;
            
        elseif abs(del_w_value + beta_w_value) < t
            
            DELTA_W_matrix(row,col) = 0;
            
        else 
            
            DELTA_W_matrix(row,col) = -1;
            
        end
              
        
        % NORTH direction (lookup direction 3)        
        beta_north_value = max(north_angles);   
        del_north_value = min(north_angles);         
        if (del_north_value + beta_north_value) > t
            
            DELTA_N_matrix(row,col) = 1;
            
        elseif abs(del_north_value + beta_north_value) < t
            
            DELTA_N_matrix(row,col) = 0;
            
        else  
            
            DELTA_N_matrix(row,col) = -1;
            
        end        
            
        
        % SOUTH direction (lookup direction 4)                             
        beta_south_value = max(south_angles);  
        del_south_value = min(south_angles);       
        if (del_south_value + beta_south_value) > t
            
            DELTA_S_matrix(row,col) = 1;
            
        elseif abs(del_south_value + beta_south_value) < t
            
            DELTA_S_matrix(row,col) = 0;
            
        else   
            
            DELTA_S_matrix(row,col) = -1;
            
        end  
             
        
        % SOUTH-EAST direction (lookup direction 5)        
        beta_se_value = max(se_angles);         
        del_se_value = min(se_angles);       
        if (del_se_value + beta_se_value) > t
            
            DELTA_SE_matrix(row,col) = 1;
            
        elseif abs(del_se_value + beta_se_value) < t
            
            DELTA_SE_matrix(row,col) = 0;
            
        else    
            
            DELTA_SE_matrix(row,col) = -1;
            
        end
              
        
        % NORTH-EAST direction (lookup direction 6)        
        beta_ne_value = max(ne_angles);        
        del_ne_value = min(ne_angles);       
        if (del_ne_value + beta_ne_value) > t
            
            DELTA_NE_matrix(row,col) = 1;
            
        elseif abs(del_ne_value + beta_ne_value) < t
            
            DELTA_NE_matrix(row,col) = 0;
            
        else           
            
            DELTA_NE_matrix(row,col) = -1;
            
        end
              
        
        % SOUTH-WEST direction (lookup direction 7)        
        beta_sw_value = max(sw_angles);        
        del_sw_value = min(sw_angles);         
        if (del_sw_value + beta_sw_value) > t
            
            DELTA_SW_matrix(row,col) = 1;
            
        elseif abs(del_sw_value + beta_sw_value) < t
            
            DELTA_SW_matrix(row,col) = 0;
            
        else     
            
            DELTA_SW_matrix(row,col) = -1;
            
        end
              
        
        % NORTH-WEST direction (lookup direction 8)        
        beta_nw_value = max(nw_angles);       
        del_nw_value = min(nw_angles);       
        if (del_nw_value + beta_nw_value) > t
            
            DELTA_NW_matrix(row,col) = 1;
            
        elseif abs(del_nw_value + beta_nw_value) < t
            
            DELTA_NW_matrix(row,col) = 0;
            
        else   
            
            DELTA_NW_matrix(row,col) = -1;
            
        end
               
    end
end

% Concatenating DELTA array (3 dimensions)
DELTA_array = cat(3, DELTA_N_matrix, DELTA_NE_matrix, DELTA_E_matrix,...
    DELTA_SE_matrix, DELTA_S_matrix, DELTA_SW_matrix, DELTA_W_matrix, DELTA_NW_matrix);   


%% CLASSIFYING GEOMORPHONS

case_matrix = nan(size(DELTA_array,1),size(DELTA_array,2));
DELTA = permute(DELTA_array,[3,2,1]);        

% Geomorphons of the 10 main landforms
flat =      [0,0,0,0,0,0,0,0]';
peak =      [-1,-1,-1,-1,-1,-1,-1,-1]';
ridge =     [-1,-1,0,-1,-1,-1,0,-1]';
shoulder =  [-1,-1,0,0,0,0,0,-1]';
spur =      [-1,-1,-1,1,1,1,-1,-1]';
slope =     [1,1,0,-1,-1,-1,0,1]';
pit =       [1,1,1,1,1,1,1,1]';
valley =    [1,1,0,1,1,1,0,1]';
footslope = [1,1,0,0,0,0,0,1]';
hollow =    [1,1,1,-1,-1,-1,1,1]';

geom_lists = [flat, peak, ridge, shoulder, spur, ...
    slope, pit, valley, footslope, hollow];
geoms = permute(repmat(geom_lists, [1,1,8]), [1,3,2]);
num_geoms = size(geom_lists, 2);

sorted_geoms = cell(num_geoms + 1, 1); 
for i = 1:length(sorted_geoms)
    sorted_geoms{i} = double.empty(8,0);
end

for j = 1:size(DELTA,2) 
    for n = 1:size(DELTA,3) 
        shifted_vecs = [
            circshift(DELTA(:,j,n), 0),...
            circshift(DELTA(:,j,n), 1),...
            circshift(DELTA(:,j,n), 2),...
            circshift(DELTA(:,j,n), 3),...
            circshift(DELTA(:,j,n), 4),...
            circshift(DELTA(:,j,n), 5),...
            circshift(DELTA(:,j,n), 6),...
            circshift(DELTA(:,j,n), 7)];

        % Sum of squared differences
        delta = sum((geoms - shifted_vecs).^2, 1);
        % Best match for each geomorphon
        match_min = permute(min(delta, [], 2), [1,3,2]);

        [val, idx] = min(match_min);
        if length(find(match_min == val)) == 1
            % Classify as unique geomorphon
            sorted_geoms{idx}(:, end+1) = DELTA(:,j,n);
            case_matrix(n,j) = idx-1;

        else
            % Classify as ambiguous geomorphon (unsorted)
            sorted_geoms{end}(:, end+1) = DELTA(:,j,n);        
            case_matrix(n,j) = 10;        
        end    
    end
end

% Display results
fprintf('%.2f of random patterns classified\n', ...
    100 * (1 - size(sorted_geoms{end}, 2) / (size(dem,1) * size(dem,2))^2))


%% COLORMAP

% Defining colormap to visualize the different types of geomorphons
% (the colors are defined by specifying a three-column matrix of RGB 
% triplets where each row defines one color):
% 
% flat:     white       [1 1 1]         case 0
% peak:     red         [.6 0 0]        case 1
% ridge:    light red   [.9 0 0]        case 2
% shoulder: orange      [.9 .5 0]       case 3
% spur:     green       [.8 .6 .5]      case 4
% slope:    yellow      [.9 .9 .3]      case 5
% pit:      black       [0 0 0]         case 6
% valley:   dark blue   [.1 .2 .4]      case 7
% footslope:green       [.1 .6 .5]      case 8
% hollow:   light blue  [.8 .9 1]       case 9      
% unsorted: gray        [0.5 0.5 0.5]   case 10

% Colormap matrix defined as 11x3 matrix
colors = [1 1 1 ; .6 0 0 ; .9 0 0 ; .9 .5 0 ; .8 .6 .5 ; .9 .9 .3 ; 0 0 0 ;...
    .1 .2 .4 ; .1 .6 .5 ; .8 .9 1 ; 0.5 0.5 0.5];


%% GEOMORPHOMETRIC MAP

% Geomorphomteric map size reduced in size (reduced by window size in both
% dimensions)
case_matrix_red = case_matrix((window_size + 1)/2 : size(DELTA_W_matrix, 1) -...
    (window_size - 1)/2, (window_size + 1)/2 : size(DELTA_W_matrix, 2) -...
    (window_size - 1)/2); 

figure('Name', [filename ' tdeg=' num2str(tdegree) ' Window size='...
    num2str(window_size) ' Skip=' num2str(skip)]); 
ax = axes();
imagesc(case_matrix_red);
colormap(colors);
colorTitleHandle = get(colorbar,'Title');   
titleString = 'Type of Geomorphon';
set(colorTitleHandle ,'String',titleString);

c=colorbar('TickLabels',{'Flat','Peak','Ridge','Shoulder','Spur','Slope',...
    'Pit','Valley','Footslope','Hollow','Unsorted'});
c.Ticks = [.5:(10/11):10.5];
ax.CLim = [-0.5,11.5];
caxis auto;
xlabel('x direction [cells]');
ylabel('y direction [cells]');
%set(gca, 'YDir', 'normal');
title([filename ' tdeg=' num2str(tdegree) ' Window size=' num2str(window_size)...
    ' Skip=' num2str(skip)]);
axis image;


%% DEM FOR COMPARISON

% Reducing DEM to same region as depicted in geomorphometric map
dem_red = dem((window_size+1)/2:size(DELTA_W_matrix,1)-(window_size-1)/2,...
    (window_size+1)/2:size(DELTA_W_matrix,2)-(window_size-1)/2);  

figure('Name',[filename ' original DEM']);
surf(dem_red, 'LineStyle', 'none');
xlabel('x direction [cells]');
ylabel('y direction [cells]');
zlabel('z elevation [cells]');



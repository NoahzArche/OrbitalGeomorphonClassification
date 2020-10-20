function GeomorphonAsRockDetector(NCfile)

% Reading NavCam DEM and reducing it to densly covered area

%t = Tiff('N_L000_1830_XYZ066ORR_S_0952_45RNGM1.TIF');
%sol = 'Sol0532';
t = Tiff(NCfile);
sol = NCfile(end-32:end-29);
dem = read(t);

% reducing DEM to region closest to rover (20x20 m^2 region centering
% rover)
dem = dem(3500:end-3501,3500:end-3501);

% mask to find most dense regions in NavCam DEM
mask = dem;
mask(mask ~= 0) = 1;
mask_opened = imopen(mask, strel('square', 150));
mask_opened = imopen(mask_opened, strel('square', 100));
mask_opened = bwareafilt(logical(mask_opened),...
    [50000 size(dem,1)*size(dem,2)]);

% rows and cols of where heightmap actually start
[n_row,n_col] = find(mask_opened ~= 0);   
min_row = min(n_row);
max_row = max(n_row);
min_col = min(n_col);
max_col = max(n_col);
dem = dem(min(n_row):max(n_row), min(n_col):max(n_col));

dem(dem == 0) = nan;


% shifting coordinates for accurate localisation of evaluated terrain
% regions
shift_row = (min_row + max_row - 1)/2 - 1000;
shift_col = (min_col + max_col - 1)/2 - 1000;


%% Configuration Parameters

filename = sol;

res = 0.01;                     % resolution in meters (for NavCam DEMs)
window_size = 101;
skip = 3;
tdegree = 10;                   % elevation flatness threshold in degrees
t = tdegree * pi/180;           % convert threshold to radians
scanning_size = size(dem,1);  




%%

% The DEL w/e/n/s/nw/ne/sw/se matrices contain the elements of the ternary
% patterns (i.e. -1's, 0's, and 1's)
DELTA_W_matrix = nan(size(dem)); 
DELTA_E_matrix = nan(size(DELTA_W_matrix));
DELTA_N_matrix = nan(size(DELTA_W_matrix));
DELTA_S_matrix = nan(size(DELTA_W_matrix));
DELTA_NW_matrix = nan(size(DELTA_W_matrix));
DELTA_NE_matrix = nan(size(DELTA_W_matrix));
DELTA_SW_matrix = nan(size(DELTA_W_matrix));
DELTA_SE_matrix = nan(size(DELTA_W_matrix));

% Matrix storing the number of extreme elevation angles for each cell in
% the DEM
DELTA_EXTREME = zeros(size(DELTA_W_matrix));
% extreme flatness threshold
t_EXTREME_deg = 60;
t_EXTREME = t_EXTREME_deg * pi/180; 


% Distance from central cell to all the neighboring cells
% index of neighboring cells which are taken into account (ignoring skipped ones)
i = 1 : (window_size-1)/2 * (1/skip);                
% distance in m
distance(i) = res * skip * i;                  
distance_diag = sqrt(2) * (res * skip * i); 



%%  Ternary Pattern Computation

% going through each pixel in the NavCam DEM 
for row = (window_size + 1)/2 : size(dem,1) - (window_size - 1)/2   
    for col = (window_size + 1)/2 : size(dem,2) - (window_size - 1)/2
        
        if ~ isnan(dem(row,col))
        
        % elevations of all the neighboring cells relative to central cell
        east_relative_elevation = double(dem( row, (col + skip) : skip :...
            col + (window_size - 1)/2 )) - double(dem(row,col));    
        west_relative_elevation = double(dem( row, (col - skip) : -skip :...
            col - (window_size - 1)/2 )) - double(dem(row,col));    
        north_relative_elevation = double(dem( (row - skip) : -skip :...
            row - (window_size - 1)/2, col )) - double(dem(row,col));
        south_relative_elevation = double(dem( (row + skip) : skip:...
            row + (window_size - 1)/2, col )) - double(dem(row,col));               
        nw_relative_elevation = double(dem( (col - 1 - skip) *...
            scanning_size + (row - skip) : - (scanning_size + 1) * skip :...
            (col - (window_size-1)/2 - 1) * scanning_size ))...
            - double(dem(row,col));        
        sw_relative_elevation = double(dem( (col - 1 - skip) *...
            scanning_size + (row + skip) :...
            - (scanning_size - 1) * skip : (col - (window_size-1)/2 -1) *...
            scanning_size +1 )) - double(dem(row,col));  
        
        ne_relative_elevation = double(dem( (col + skip -1) *...
            scanning_size + (row - skip) :...
            (scanning_size - 1) * skip : (col + (window_size-1)/2) *...
            scanning_size -1)) - double(dem(row,col));
        se_relative_elevation = double(dem( (col + skip -1) *...
            scanning_size + (row + skip) :(scanning_size + 1) * skip :...
            (col + (window_size-1)/2) * scanning_size ))...
            - double(dem(row,col));
        
        % claculate elevation angles in each direction
        east_angles = atan( east_relative_elevation./ distance );
        west_angles = atan( west_relative_elevation./ distance );
        north_angles = atan( north_relative_elevation'./ distance );
        south_angles = atan( south_relative_elevation'./ distance );        
        nw_angles = atan( double(nw_relative_elevation)./ distance_diag );
        se_angles = atan( double(se_relative_elevation)./ distance_diag );
        sw_angles = atan( double(sw_relative_elevation(1:...
            size(distance_diag,2)))./ distance_diag );
        ne_angles = atan( double(ne_relative_elevation(1:...
            size(distance_diag,2)))./ distance_diag );
        
        
        % Calculating elevation angles
        % starting with EAST direction (lookup direction 1)
        % beta and del denote max and min elevation angles, respectively 
        beta_e_value = max(east_angles);             
        del_e_value = min(east_angles);         
        % Distinction of cases for each direction, resulting in the ternary
        % pattern value stored in the DEL matrices
        if (del_e_value + beta_e_value) > t
            
            DELTA_E_matrix(row,col) = 1;  
            
            if (del_e_value + beta_e_value) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1;
                
            end
            
        elseif abs(del_e_value + beta_e_value) < t
            
            DELTA_E_matrix(row,col) = 0;        
            
        else    
            
            DELTA_E_matrix(row,col) = -1;
            
            if (abs(del_e_value + beta_e_value)) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1;
                
            end
            
        end
        
        
        % This process is repeated for all the other (7) lookup directions
                
        % WEST direction (lookup direction 2)                     
        beta_w_value = max(west_angles);           
        del_w_value = min(west_angles);         
        if (del_w_value + beta_w_value) > t
            
            DELTA_W_matrix(row,col) = 1;
            
            if (del_w_value + beta_w_value) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        elseif abs(del_w_value + beta_w_value) < t
            
            DELTA_W_matrix(row,col) = 0;
            
        else 
            
            DELTA_W_matrix(row,col) = -1;
            
            if (abs(del_w_value + beta_w_value)) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        end
                
        
        % NORTH direction (lookup direction 3)        
        beta_north_value = max(north_angles);   
        del_north_value = min(north_angles);         
        if (del_north_value + beta_north_value) > t
            
            DELTA_N_matrix(row,col) = 1;
            
            if (del_north_value + beta_north_value) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        elseif abs(del_north_value + beta_north_value) < t
            
            DELTA_N_matrix(row,col) = 0;
            
        else  
            
            DELTA_N_matrix(row,col) = -1;
            
            if (abs(del_north_value + beta_north_value)) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        end        
                     
        
        % SOUTH direction (lookup direction 4)                             
        beta_south_value = max(south_angles);  
        del_south_value = min(south_angles);       
        if (del_south_value + beta_south_value) > t
            
            DELTA_S_matrix(row,col) = 1;
            
            if (del_south_value + beta_south_value) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        elseif abs(del_south_value + beta_south_value) < t
            
            DELTA_S_matrix(row,col) = 0;
            
        else   
            
            DELTA_S_matrix(row,col) = -1;
            
            if (abs(del_south_value + beta_south_value)) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        end  
                
        
        % SOUTH-EAST direction (lookup direction 5)        
        beta_se_value = max(se_angles);         
        del_se_value = min(se_angles);       
        if (del_se_value + beta_se_value) > t
            
            DELTA_SE_matrix(row,col) = 1;
            
            if (del_se_value + beta_se_value) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        elseif abs(del_se_value + beta_se_value) < t
            
            DELTA_SE_matrix(row,col) = 0;
            
        else    
            
            DELTA_SE_matrix(row,col) = -1;
            
            if (abs(del_se_value + beta_se_value)) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        end
                
        
        % NORTH-EAST direction (lookup direction 6)        
        beta_ne_value = max(ne_angles);        
        del_ne_value = min(ne_angles);       
        if (del_ne_value + beta_ne_value) > t
            
            DELTA_NE_matrix(row,col) = 1;
            
            if (del_ne_value + beta_ne_value) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        elseif abs(del_ne_value + beta_ne_value) < t
            
            DELTA_NE_matrix(row,col) = 0;
            
        else           
            
            DELTA_NE_matrix(row,col) = -1;
            
            if (abs(del_ne_value + beta_ne_value)) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
        end
                     
        
        % SOUTH-WEST direction (lookup direction 7)        
        beta_sw_value = max(sw_angles);        
        del_sw_value = min(sw_angles);         
        if (del_sw_value + beta_sw_value) > t
            
            DELTA_SW_matrix(row,col) = 1;
            
            if (del_sw_value + beta_sw_value) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        elseif abs(del_sw_value + beta_sw_value) < t
            
            DELTA_SW_matrix(row,col) = 0;
            
        else     
            
            DELTA_SW_matrix(row,col) = -1;
            
            if (abs(del_sw_value + beta_sw_value)) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
        end
                
        
        % NORTH-WEST direction (lookup direction 8)        
        beta_nw_value = max(nw_angles);       
        del_nw_value = min(nw_angles);       
        if (del_nw_value + beta_nw_value) > t
            
            DELTA_NW_matrix(row,col) = 1;
            
            if (del_nw_value + beta_nw_value) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
            
        elseif abs(del_nw_value + beta_nw_value) < t
            
            DELTA_NW_matrix(row,col) = 0;
            
        else   
            
            DELTA_NW_matrix(row,col) = -1;
            
            if (abs(del_nw_value + beta_nw_value)) > t_EXTREME
                
                DELTA_EXTREME(row,col) = 1 + DELTA_EXTREME(row,col);
                
            end
        end
        
        end
    end
end

% Computing DELTA array (3 dimensions)
DELTA_array = cat(3, DELTA_N_matrix, DELTA_NE_matrix, DELTA_E_matrix,...
    DELTA_SE_matrix, DELTA_S_matrix, DELTA_SW_matrix, DELTA_W_matrix,...
    DELTA_NW_matrix);   

case_matrix = nan(size(DELTA_array,1),size(DELTA_array,2));
DELTA = permute(DELTA_array,[3,2,1]);        


%% CLASSIFYING TERNARY PATTERNS INTO GEOMORPHON CLASSES

% Class 1 patterns
class1_1 =  [-1,-1,-1,-1,-1,-1,-1,-1]'; % ideal
class1_2 =  [-1,-1,0,-1,-1,-1,-1,-1]';  % 1 zero
class1_3 =  [-1,-1,-1,-1,0,0,-1,-1]';   % 2 zeros

% Class 2 patterns
class2_1 =  [-1,-1,0,-1,-1,-1,0,-1]';   % ideal
class2_2 =  [-1,-1,0,-1,-1,0,-1,-1]';   % closer
class2_3 =  [-1,-1,0,-1,0,-1,-1,-1]';   % even closer
class2_4 =  [-1,-1,0,-1,-1,-1,0,0]';    % two zeros opposite of one 1
class2_5 =  [-1,-1,0,-1,-1,0,0,-1]';    % two zeros opposite of one 1
class2_6 =  [-1,0,0,-1,-1,0,0,-1]';     % two zeros opposite of two zeros
class2_7 =  [-1,0,-1,0,-1,0,-1,0]';     % a zero in each direction

% Class 3 patterns
class3_1 =  [-1,-1,0,0,0,0,0,-1]';      % 5 zeros
class3_2 =  [-1,-1,0,0,0,0,-1,-1]';     % 4 zeros
class3_3 =  [-1,-1,0,0,0,-1,-1,-1]';    % 3 zeros
class3_4 =  [-1,-1,0,0,-1,0,-1,-1]';    % 3 zeros with -1 between
class3_5 =  [-1,-1,0,-1,0,0,-1,-1]';    % 3 zeros with -1 between
class3_6 =  [-1,-1,0,0,0,-1,0,-1]';     % 4 zeros with -1 between
class3_7 =  [-1,-1,0,-1,0,0,0,-1]';     % 4 zeros with -1 between
class3_8 =  [-1,-1,0,0,-1,0,0,-1]';     % 4 zeros with -1 between
class3_9 =  [-1,0,0,-1,0,0,-1,0]';      % 4 zeros with -1 between

% Class 4 patterns
class4_1 =  [-1,-1,-1,1,-1,-1,-1,-1]';  % 1 ones
class4_2 =  [-1,-1,-1,1,1,-1,-1,-1]';   % 2 ones
class4_3 =  [-1,0,1,-1,-1,-1,-1,-1]';   % 1 one 1 zero
class4_4 =  [-1,-1,1,0,-1,-1,-1,-1]';   % 1 one 1 zero
class4_5 =  [-1,0,1,0,-1,-1,-1,-1]';    % 1 one 2 zeros
class4_6 =  [-1,-1,-1,1,1,1,-1,-1]';    % 3 ones

% Flat and sloped
flat =      [0,0,0,0,0,0,0,0]';         % ideal
slope =     [1,1,0,-1,-1,-1,0,1]';      % ideal
slope2 =    [1,1,0,-1,-1,0,0,1]';       % 2 zeros opposite 
slope3 =    [1,1,0,-1,-1,-1,0,0]';      % 2 zeros opposite

cases = [5, 5, 5, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4,...
    2, 2, 2, 2, 2, 2, 2, 1, 1, 1];

geom_lists = [slope, slope2, slope3, flat,...
    class3_1, class3_2, class3_3, class3_4, class3_5, class3_6,...
    class3_7, class3_8, class3_9,...
    class4_1, class4_2, class4_3, class4_4, class4_5, class4_6,...
    class2_1, class2_2, class2_3, class2_4, class2_5, class2_6, class2_7,...
    class1_1, class1_2, class1_3];
geoms = permute(repmat(geom_lists, [1,1,8]), [1,3,2]);
num_geoms = size(geom_lists, 2);

sorted_geoms = cell(num_geoms +1, 1); 
for g = 1:length(sorted_geoms)
    sorted_geoms{g} = double.empty(8,0);
end

for j = 1:size(DELTA,2) 
    for n = 1:size(DELTA,3) 
        if ~ isnan(dem(n,j))

        shifted_vecs = [
            circshift(DELTA(:,j,n), 0),...
            circshift(DELTA(:,j,n), 1),...
            circshift(DELTA(:,j,n), 2),...
            circshift(DELTA(:,j,n), 3),...
            circshift(DELTA(:,j,n), 4),...
            circshift(DELTA(:,j,n), 5),...
            circshift(DELTA(:,j,n), 6),...
            circshift(DELTA(:,j,n), 7)];

        % calculating sum of squared differences
        delta = sum(bsxfun(@minus,geoms,shifted_vecs).^2,1);
        % calculating best match for each ternary pattern
        match_min = permute(min(delta, [], 2), [1,3,2]);

        [val, idx] = min(match_min);
        if length(find(match_min == val)) == 1

            % categorizing as one of the geomorphon classes
            sorted_geoms{idx}(:, end+1) = DELTA(:,j,n);
            case_matrix(n,j) = cases(idx);

        else

            % categorizing as ambiguous (unsorted) ternary pattern
            sorted_geoms{end}(:, end+1) = DELTA(:,j,n);        
            case_matrix(n,j) = 10;        
        end 
        else

            case_matrix(n,j) = 10; 

        end
    end
end

% display results
disp(NCfile)
fprintf('%.2f of ternary patterns classified\n', ...
    100 * (1 - size(sorted_geoms{end}, 2) / (size(dem,1) * size(dem,2))^2))


%% Geomorphon-Class Map

% Visulazing the classified geomorphon classes in an intermediate 
% (geomorphometric) map
% The colors are defined by specifying a three-column matrix of RGB 
% triplets where each row defines one color:

% Class 1:          dark red    [.6 0 0]    case 1
% Class 2:          light red   [.9 0 0]    case 2
% Class 3:          orange      [.9 .5 0]   case 3
% Class 4:          coral       [.8 .6 .5]  case 4
% Flat (ground):    white       [1 1 1]     case 0
% Slope (ground):   yellow      [.9 .9 .3]  case 5
% Unsorted:         gray        [.5 .5 .5]  case 10

% Colormap matrix defined as 10x3 matrix
colors = [1 1 1 ; .6 0 0 ; .9 0 0 ; .9 .5 0 ; .8 .6 .5 ; .9 .9 .3 ; ...
    0 0 0 ; .1 .2 .4 ; .1 .6 .5 ; .8 .9 1 ; 0.5 0.5 0.5];

% assigning the nodata voids color black for visualisation
case_matrix(isnan(dem)) = 6;    
case_matrix_red = case_matrix((window_size + 1)/2 : size(case_matrix, 1) -...
    (window_size - 1)/2, (window_size + 1)/2:size(case_matrix, 2) -...
    (window_size - 1)/2); 

% figure('Name', [filename ' tdeg=' num2str(tdegree) ' Window size='...
%       num2str(window_size) ' Skip=' num2str(skip)]); 
% ax = axes();
% imagesc(case_matrix_red);
% colormap(colors);
% colorTitleHandle = get(colorbar,'Title');   
% titleString = 'Type of Geomorphon';
% set(colorTitleHandle ,'String',titleString);
% c=colorbar('TickLabels',{'Flat', 'Peak', 'Ridge', 'Shoulder','Spur',...
%     'Slope', 'Nodata regions', ' ', ' ', ' ', 'Unsorted'});
% c.Ticks = [.5:(10/11):10.5];
% ax.CLim = [-0.5,11.5];
% caxis auto;
% xlabel('x direction [cells]');
% ylabel('y direction [cells]');
% set(gca, 'YDir', 'normal');
% axis image;
% title([filename ' tdeg=' num2str(tdegree) ' Window size='...
%     num2str(window_size) ' Skip=' num2str(skip)]);

%% DEM FOR COMPARISON

dem_red = dem((window_size + 1)/2:size(dem, 1) - (window_size - 1)/2,...
    (window_size + 1)/2:size(dem, 2) - (window_size - 1)/2);  

% [xgrid_red, ygrid_red] = meshgrid( 0:res:(size(dem_red,1)-1) * res,...
%     0:res:(size(dem_red,2)-1) * res); 
% figure('Name',[filename ' original DEM']);
% surf(dem_red, 'LineStyle', 'none');
% xlabel('x direction [m]');
% ylabel('y direction [m]');
% zlabel('z elevation [m]');


%% Saving data



% saving variables into mat file
save(['OutputSol' sol], 'case_matrix', 'dem', 'DELTA_EXTREME', 'min_row',...
    'max_row', 'min_col', 'max_col', 'shift_row', 'shift_col')

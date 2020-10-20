
% Post-processing script for geomorphometric maps
% This script can be run after the geomorphon classification algorithm 
% to classify the ambiguous (unsorted) geomorphons in the resulting 
% (intermediate) geomorphometric map. 


% Binarizing unsorted regions (case 10 cells)
unsorted_geomap = case_matrix_red;
unsorted_geomap(unsorted_geomap ~= 10 ) = 0;
unsorted_geomap(unsorted_geomap ==10 ) = 1;  


% Displaying only the unsorted regions in the geomorphometric map
% figure('Name', 'Unsorted geomorphons binarized');
% imshow(unsorted_geomap);
% hold on;
[B_uns] = bwboundaries(unsorted_geomap,'noholes');
% colors_b=['b' 'g' 'r' 'c' 'm' 'y'];
% for k=1:length(B_uns)
%   boundary_uns = B_uns{k};                        % number of boundary
%   cidx_uns = mod(k,length(colors_b))+1;           % random color choice
%   plot(boundary_uns(:,2), boundary_uns(:,1), colors_b(cidx_uns),...
%       'LineWidth',1);
% end
% set(gca, 'YDir', 'normal');                       % to display correctly


%% Splitting the unsorted regions into sub-regions 

% Creating sub-regions by morphologically opening the unsorted regions
se = strel('cube', 17);
unsorted_subregions = imopen(unsorted_geomap, se); 

[B_unsorted_subregions, L_unsorted_subregions, N_unsorted_subregions] = ...
    bwboundaries(unsorted_subregions, 'noholes');


% figure('Name', 'Morphologically opened unsorted geomorphons');
% imshow(unsorted_subregions);
% hold on;
% for k=1:length(B_unsorted_subregions)
%   boundary_uns_sub = B_unsorted_subregions{k};    % number of boundary
%   cidx_uns_sub = mod(k,length(colors_b))+1;       % random coloring
%   plot(boundary_uns_sub(:,2), boundary_uns_sub(:,1),...
%       colors_b(cidx_uns_sub),'LineWidth',1);
% end
% set(gca, 'YDir', 'normal');                       % to display correctly


% Binary image of remaining unsorted regions (diff)
unsorted_diff = unsorted_geomap - unsorted_subregions;    

[B_unsorted_diff, L_unsorted_diff, N_unsorted_diff] = ...
    bwboundaries(unsorted_diff, 'noholes');


% Creating regions for the remaining unsorted geomorphons
% figure('Name', [filename ' Remaining unsorted geomorphons']);
% imshow(unsorted_diff);
% hold on;
% for k=1:length(B_unsorted_diff)
%   boundary_uns_diff = B_unsorted_diff{k};         % number of boundary
%   cidx_uns_rdiff = mod(k,length(colors_b))+1;     % random color choice
%   plot(boundary_uns_diff(:,2), boundary_uns_diff(:,1),...
%       colors_b(cidx_uns_rdiff),'LineWidth',1);
% end
% set(gca, 'YDir', 'normal'); 



%% CLASSIFYING THE UNSORTED-REGIONS

% Sorting the remaining (smaller) regions first

% structuring elementused for dilation
se = strel('cube', 3);                                   
case_numberchange = case_matrix_red;
% flat regions are now numbered 11 instead of 0 for computation
case_numberchange(case_numberchange == 0) = 11;    
% unsorted cases are now numbered 0 (will be changed back at the end)
case_numberchange(case_numberchange == 10) = 0;         
case_matrix_ext = case_numberchange;

% going through each of the smaller unsorted regions
for diff = 1 : N_unsorted_diff                               
    
    L_copy = L_unsorted_diff;  
    % zeros except for where unsorted geomorphon is (labeled with ones 
    % (numbered) for first one, 2's for second one etc.)
    L_copy(L_copy ~= diff) = 0;                         
    L_copy(L_copy == diff) = 1;    
    % the region of the first unsorted geomorphon is slightly increased
    L_dilate = imdilate(L_copy, se);    
    % binary image of outer boundary of the first unsorted geomorphon
    L_diff = L_dilate - L_copy;     
    
    % find case values within the boundary of the unsorted geomorphon
    local_case_matrix = L_diff .* case_numberchange;    
    local_case_matrix(local_case_matrix == 0) = nan;    
    
    % most common case within boundary
    case_number = mode(local_case_matrix,'all');         
    
    % the unsorted geomorphon region is now classified, matrix contains the 
    % newly sorted geomorphon (rest 0)
    classified_region = case_number * L_copy;           
                  
    case_matrix_ext = case_matrix_ext + classified_region; 
  
end

final_case_matrix = case_matrix_ext;



% Sorting the (larger) sub-regions (same process as above)

% going through each of the larger unsorted regions
for sub = 1 : N_unsorted_subregions                                 
    
    L_copy = L_unsorted_subregions;   
    % binarizing unsorted geomorphon region 
    L_copy(L_copy ~= sub) = 0;                          
    L_copy(L_copy == sub) = 1;
    % the region of the first unsorted geomorphon is slightly increased
    L_dilate = imdilate(L_copy, se);                    
    % binary image of outer boundary of the first unsorted geomorphon 
    L_diff = L_dilate - L_copy;                             
    
    % find case values within the boundary of the unsorted geomorphon
    local_case_matrix = L_diff .* case_matrix_ext;      
    local_case_matrix(local_case_matrix == 0) = nan; 
    
    % most common case within boundary 
    case_number = mode(local_case_matrix, 'all');       
   
    % the unsorted geomorphon region is now classified, matrix contains the 
    % newly sorted geomorphon (rest 0)
    classified_region = case_number * L_copy;           
                  
    final_case_matrix = final_case_matrix + classified_region; 
  
end

% return to initial case numbering
final_case_matrix(final_case_matrix == 11) = 0;         
final_case_matrix(final_case_matrix == 10) = nan;


%% FINAL GEOMORPHOMETRIC MAP (NO UNSORTED REGIONS)

% post-processed geomorphometric map (fully sorted)
colors1 = [1 1 1;.6 0 0;.9 0 0;.9 .5 0;.8 .6 .5;.9 .9 .3;0 0 0;.1 .2 .4;.1 .6 .5;.8 .9 1];
figure('Name', ['Final Geomorphometric Map tdeg=' num2str(tdegree) ' Window size=' num2str(window_size) ' Skip=' num2str(skip)]); 
ax1 = axes();
imagesc(final_case_matrix);
colormap(colors1);
colorTitleHandle = get(colorbar,'Title');               % setting a title for the colormap
titleString = 'Type of Geomorphon';
set(colorTitleHandle ,'String',titleString);
c=colorbar('TickLabels',{'Flat','Peak','Ridge','Shoulder','Spur','Slope',...
    'Pit','Valley','Footslope','Hollow'});
c.Ticks = [.5:.9:9.5];
ax1.CLim = [-0.5,10.5];
caxis auto;
xlabel('x direction [cells]');
ylabel('y direction [cells]');
%set(gca, 'YDir', 'normal');
title(['Final Geomorphometric Map tdeg=' num2str(tdegree) ' Window size=' num2str(window_size) ' Skip=' num2str(skip)]);
axis image;


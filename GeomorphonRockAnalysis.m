% Processing detected geomorphon class regions


GeomapOutputfiles = dir(['Output' '*.mat']);
number_of_files = length(GeomapOutputfiles);

dlmwrite(['NCLoopGeomapAnalysis' num2str(date) '.csv'],...
    {'Looping through sols'}, 'delimiter', '');
dlmwrite(['NCLoopGeomapAnalysis' num2str(date) '.csv'],...
    date, 'delimiter', '', '-append');


%% Configuration Parameters

window_size = 101; 
res = 0.01;



%% Processing Classified Ternary Patterns

% looping through every 
for i = 1:number_of_files
    
    GeomapOutputfile = GeomapOutputfiles(i).name;
    solnumber = GeomapOutputfile(end-10:end-4); 
    
    try 
        
        load(GeomapOutputfile);
        % contains 'case_matrix', 'dem', 'min_row', 'max_row', 'min_col',
        % 'max_col', 'shift_row', 'shift_col'
        
        case_matrix_red = case_matrix((window_size + 1)/2 : size(dem, 1)-...
            (window_size-1)/2, (window_size + 1)/2 : size(dem, 2) -...
            (window_size - 1)/2); 
        dem_red = dem((window_size + 1)/2 : size(dem, 1) - (window_size - 1)/2,...
             (window_size + 1)/2 : size(dem, 2) - (window_size - 1)/2);
        DELTA_EXTREME_red = DELTA_EXTREME((window_size + 1)/2 : size(dem, 1)-...
             (window_size - 1)/2, (window_size + 1)/2 : size(dem, 2) -...
             (window_size - 1)/2);
       
         
        %% Masking to remove NavCam nodata voids and their immediate surroundings
        
        % inverted mask marking nodata voids
        dem_mask_inv = zeros(size(dem_red));        
        dem_mask_inv(isnan(dem_red)) = 1;       
        % dilating nodata voids with different sized structing elements
        % dependent on the nodata void's filled area
        dem_mask_inv_20 = bwareafilt(logical(dem_mask_inv),[1 20]);
        dem_mask_inv_100 = bwareafilt(logical(dem_mask_inv),[21 100]);
        dem_mask_inv_100 = imdilate(dem_mask_inv_100, strel('disk', 3));
        dem_mask_inv_1000 = bwareafilt(logical(dem_mask_inv),[101 1000]);
        dem_mask_inv_1000 = imdilate(dem_mask_inv_1000, strel('disk', 9));
        dem_mask_inv_max = bwareafilt(logical(dem_mask_inv),[1001 4000000]);
        dem_mask_inv_max = imdilate(dem_mask_inv_max, strel('disk', 12)); 
        
        dem_mask_inv = dem_mask_inv + dem_mask_inv_20 + dem_mask_inv_100 +...
            dem_mask_inv_1000 + dem_mask_inv_max;        
        dem_mask_inv = bwareafilt(~dem_mask_inv,...
            [1500 size(dem_red,1) * size(dem_red,2)]) == 0;
        dem_mask_inv = imclose(dem_mask_inv, strel('disk', 10));
        dem_red(dem_mask_inv == 1) = nan;
        case_matrix_red(dem_mask_inv == 1) = 6;
        
        
        
        %% Hierarchical Attachment Process of Geomorphon Classes   
        
        % starting with geomorphon class 1 and 2 cell regions
        binary_rocks = case_matrix_red;
        % binarizing class 1 and class 2 geomorphons
        binary_rocks(binary_rocks == 2) = 1;
        binary_rocks(binary_rocks ~= 1) = 0;
        binary_rocks = imclose(binary_rocks,strel('disk', 1));
        binary_rocks = imfill(binary_rocks, 'holes'); 
                        
        binary_rocks_c = zeros(size(case_matrix_red));
        % dilating class 1 cell regions
        binary_rocks_c(case_matrix_red == 1) = 1;
        binary_rocks_c = imdilate(binary_rocks_c, strel('disk', 1));
        % adding class 2 cell regions
        binary_rocks_c(case_matrix_red == 2) = 1;
        binary_rocks_c = imclose(binary_rocks_c, strel('disk', 2));       
         
        
        % attaching cells of geomorphon classes 1, 2, 3, and 4
        classes_1_2_3_4 = case_matrix_red; 
        classes_1_2_3_4(classes_1_2_3_4 == 2) = 1; % class 2
        classes_1_2_3_4(classes_1_2_3_4 == 3) = 1; % class 3
        classes_1_2_3_4(classes_1_2_3_4 == 4) = 1; % class 4
        % binarizing
        classes_1_2_3_4(classes_1_2_3_4 ~= 1) = 0;
        classes_1_2_3_4 = imerode(classes_1_2_3_4, strel('disk', 1));
        classes_1_2_3_4 = imfill(classes_1_2_3_4, 'holes');
        
        % going through every class 1, 2, 3, 4 region and checking
        % if it contains a class 1, 2 region (previously marked as rock
        % region)
        [L_classes_1_2_3_4,n] = bwlabel(classes_1_2_3_4);        
        L_copy = L_classes_1_2_3_4;        
        L_copy(binary_rocks_c == 0) = 0;
        regions = unique(L_copy);
        for i = 2 : size(regions, 1)
            binary_rocks(L_classes_1_2_3_4 == regions(i)) = 1;
        end
        
        
        % Binary rocks have now been augmented by the attachment of 
        % geomorphon class 3 and 4 cells lying immediately surrounding
        % the previously marked rock cells (geomorphon class 1 and 2 cells)
        
        
        % attaching immediately surrounding unsorted cells
        cases_unsorted = case_matrix_red;
        cases_unsorted(cases_unsorted == 2) = 1;   
        cases_unsorted(cases_unsorted == 3) = 1;
        cases_unsorted(cases_unsorted == 4) = 1;
        cases_unsorted(cases_unsorted == 10) = 1;
        cases_unsorted(cases_unsorted ~= 1) = 0;
               
        % imclose the unsorted regions
        % attach them to the binary rock image regions
        cases_unsorted = imerode(cases_unsorted,strel('disk', 3));
        cases_unsorted = imclose(cases_unsorted,strel('disk', 3));
        
        % going through every unsorted region and checking if it lies
        % immediately surrounding a rock region
        [L_unsorted,n_unsorted] = bwlabel(cases_unsorted);       
        L_unsorted_copy = L_unsorted;        
        L_unsorted_copy(binary_rocks == 0) = 0;       
        regions_unsorted = unique(L_unsorted_copy);
        for i = 2 : size(regions_unsorted, 1)
            binary_rocks(L_unsorted == regions_unsorted(i)) = 1;
        end
        
        
        % Now adding the surrounding regions with extreme elevation angles
        % keeping the cells with at least 3 extreme elevation angles (i.e.
        % exceeding the extreme flatness threshold)
        DELTA_EXTREME_red(dem_mask_inv == 1) = 0;
        DELTA_EXTREME_red(DELTA_EXTREME_red < 3) = 0;
        DELTA_EXTREME_red(DELTA_EXTREME_red ~= 0) = 1;
        DELTA_EXTREME_red = imdilate(DELTA_EXTREME_red,strel('disk', 1));
        [L_extreme,n_extreme] = bwlabel(DELTA_EXTREME_red);
        % going through each region
        L_extreme_copy = L_extreme;        
        L_extreme_copy(binary_rocks == 0) = 0;      
        regions_extreme = unique(L_extreme_copy);
        for i = 2 : size(regions_extreme, 1)
            binary_rocks(L_extreme == regions_extreme(i)) = 1;
        end        

                
        % remove holes surroundings of isolated cells
        binary_rocks = imclose(binary_rocks, strel('disk',2)); 
        % fill regions between peak and ridge regions 
        binary_rocks = imfill(binary_rocks,'holes');
        
        % removing rocks containing less than 5 cells 
        binary_rocks(dem_mask_inv == 1) = 0;
        binary_rocks = bwareafilt(logical(binary_rocks),[5 10000]);
        
        
        
        
        %% Rock Stats
        
        % downsampling to 10 cm for the gradient computation and
        % normalizing the slope (i.e. [])
        dx = diff(dem(1:10:end,:),1,1)./.1;
        dy = diff(dem(:,1:10:end),1,2)./.1;
        slope = sqrt(mean(dx(:),'omitnan')^2 + mean(dy(:),'omitnan')^2);
                 
        percent_rocks = nnz(binary_rocks == 1)/(nnz(binary_rocks == 0) -...
            nnz(isnan(dem_red)));
        CC = bwconncomp(binary_rocks);
        n_rocks = CC.NumObjects;        
        stats = regionprops(CC, dem_red, 'MajorAxisLength',...
            'MinorAxisLength','EquivDiameter', 'FilledArea');
        stats_heights = regionprops(bwconncomp(imdilate(binary_rocks,...
            strel('disk', 5))), dem_red, 'MaxIntensity', 'MinIntensity');
        majoraxisl = cat(1, stats.MajorAxisLength);
        minoraxisl = cat(1, stats.MinorAxisLength);
        equivd = cat(1, stats.EquivDiameter);
        maxheight = cat(1, stats_heights.MaxIntensity);
        minheight = cat(1, stats_heights.MinIntensity);
        heights = maxheight - minheight;
        rock_area = cat(1, stats.FilledArea);
        % percent of cells nans
        percent_shadow = nnz(isnan(dem_red))/(size(dem_red,1) * size(dem_red,2));
        
        save([solnumber 'GeomapAnalysis.mat'], 'binary_rocks', 'dem_red',...
            'percent_rocks', 'percent_shadow', 'n_rocks', 'minoraxisl',...
            'majoraxisl', 'heights', 'rock_area', 'equivd', 'slope');
        
        disp(GeomapOutputfile)
        
        
        
    catch ME
        
        dlmwrite(['NCLoopGeomapAnalysis' num2str(date) '.csv'], 'failed NavCam mosaic DEM Geomap Rock Detection: ', 'delimiter', '', '-append');
        dlmwrite(['NCLoopGeomapAnalysis' num2str(date) '.csv'], GeomapOutputfile, 'delimiter', '', '-append');
        dlmwrite(['NCLoopGeomapAnalysis' num2str(date) '.csv'], ME.message, 'delimiter', '', '-append');
        
    end
end

% Generating Automated Reports Summarizing NavCam Geomap Rock results


GeomapAnalysisfiles = dir('*GeomapAnalysis.mat');
number_of_files = length(GeomapAnalysisfiles);

sols = [GeomapAnalysisfiles(1).name(4:5) '00'];
sols = str2double(sols);
Dataset = ['Dataset Sols ' num2str(sols) ' - ' num2str(sols + 100)];

window_size = 101;   
res = 0.01;
skip = 3;
tdegree = 10;

load('/Users/noahzr/Desktop/sols/GeomapSolRegions/Sorted_case_matrices.mat');
cases = eval(['cases' num2str(sols) 'sorted']);
colors = [1 1 1 ; .6 0 0 ; .9 0 0 ; .9 .5 0 ; .8 .6 .5 ; .9 .9 .3 ; 0 0 0 ;...
    .1 .2 .4 ; .1 .6 .5 ; .8 .9 1 ; 0.5 0.5 0.5];

percent_rocks_arr = zeros(size(1,number_of_files));
percent_shadow_arr = zeros(size(1,number_of_files));
n_rocks_arr = zeros(size(1,number_of_files));
%evaluatedarea_arr = zeros(size(1,number_of_files));
area_evaluated_abs_arr = zeros(size(1,number_of_files));
heights_arr = cell(1,number_of_files);
majoraxisl_arr = cell(1,number_of_files);
minoraxisl_arr = cell(1,number_of_files);
equivd_arr = cell(1,number_of_files);
q = zeros(size(1,number_of_files));
D = logspace(-2,.3);
cfa_value_arr = cell(1,number_of_files);
rock_area_arr = cell(1,number_of_files);
skipped_data = zeros(1,number_of_files);
slope_arr = zeros(size(1,number_of_files));

for i = 1:number_of_files
    
    GeomapAnalysisfile = GeomapAnalysisfiles(i).name;
    solnumber = GeomapAnalysisfile(4:7); %GeomapAnalysisfile(1:4);
    load(['OutputSol' solnumber '.mat'])
    load(GeomapAnalysisfile);
    % contains 'Binary_rocks', 'dem_red', 'percent_rocks', 'percent_shadow', 
    % 'n_rocks', 'minoraxisl', 'majoraxisl', 'heights', 'rock_area', 
    % 'equivd', 'slope'
    

% figure
% imshow(Binary_rocks)
% set(gca, 'YDir', 'normal');
% title(i);
% figure
% surf(dem_red,'LineStyle','none')
% title(i);

    if (nnz(~isnan(dem_red(:)))*0.0001 > 50) && slope <0.2 %&&  %&& (percent_shadow < 0.4) 
        
        % number of rock cells per non rock cells (nans excluded)
        percent_rocks_arr(i) = percent_rocks; 
        % number of nans in rectangle dem_red region
        percent_shadow_arr(i) = percent_shadow; 
        n_rocks_arr(i) = n_rocks;
        heights_arr(i) = mat2cell(heights, size(heights,1));
        %evaluated_area = nnz(~isnan(dem_red(:)))/(size(dem_red,1) *...
        %   size(dem_red,2));
        %evaluatedarea_arr(i) = evaluated_area;  % percentage
        area_evaluated_abs_arr(i) = nnz(~isnan(dem_red(:)));
        majoraxisl_arr(i) = mat2cell(majoraxisl, size(majoraxisl,1));
        minoraxisl_arr(i) = mat2cell(minoraxisl, size(minoraxisl,1)); 
        equivd_arr(i) = mat2cell(equivd, size(equivd,1));
        rock_area_arr(i) = mat2cell(rock_area, size(rock_area,1));
        slope_arr(i) = slope;
        
        q = 1.79 + 0.152/percent_rocks;
        % cumulative fractional area covered by rocks of diameter D or larger
        cfa_value = percent_rocks .* exp(-q * D);
        cfa_value_arr(i) = mat2cell(cfa_value, size(cfa_value,1));
               
    else
        skipped_data(i) = i;
    end
    
end



%% generating report of statistics


import mlreportgen.report.*
import mlreportgen.dom.*


rpt = Report(['RockAnalysisReportSols' num2str(sols)],'pdf');
open(rpt);
tp = TitlePage;
tp.Title = 'Rock Analysis Results';
tp.Subtitle = Dataset;
add(rpt,tp);

add(rpt, TableOfContents);


sec1 = Section;
sec1.Title = 'Dataset Information';

para1 = Paragraph(['Rocks were detected using the Geomorphon algorithm. '...
    'The window size was set to '...
    num2str(window_size) ...
    ', the number of skipped cells to ' ...
    num2str(skip) ...
    ', and the flatness threshold to ' ... 
    num2str(tdegree) ...
    ' deg.']);
para1.FontSize = '10';
add(sec1, para1);

if isempty(skipped_data(skipped_data ~= 0)) == 0
    para4 = Paragraph(['The NavCam DEMs were limited to the most densely covered '...
        'regions (dense meaning only a few small nodata regions). No surface '...
        'interpolation or other elevation data was added to the heightmaps. '...
        'Regions <= 50 m^2 were ignored (hence the missing data values for data points '...
        num2str(skipped_data(skipped_data ~= 0)) ').']);
else
    para4 = Paragraph(['The NavCam DEMs were limited to the most densely covered '...
        'regions (dense meaning only a few small nodata regions). No surface '...
        'interpolation or other elevation data was added to the heightmaps.)']);
end
para4.FontSize = '10';
add(sec1, para4);

para2 = Paragraph(['A total of ' num2str(number_of_files)...
    ' NavCam mosaics were evaluated.']);
para2.FontSize = '10';
add(sec1,para2);


p = Paragraph([num2str(1) ' = ' 'Sol ' GeomapAnalysisfiles(1).name(4:7)...
    ' // ' ' ' ' ' ' ' ' ' ' ']);
for i = 2:number_of_files
    %p = Paragraph([num2str(i) ' - ' GeomapAnalysisfiles(i).name(4:7) '...
    %   ' num2str(i+1) ' - ' GeomapAnalysisfiles(i+1).name(4:7)]);
    %p.Style = [p.Style {FontFamily('Menlo'), FontSize('8')}];
    append(p,[num2str(i) ' = ' 'Sol ' GeomapAnalysisfiles(i).name(4:7)...
        ' // ' ' ' ' ' ' ' ' ' ' ']);
end
%append(p, mlreportgen.dom.PageBreak)
p.FontFamilyName = 'Menlo';
p.FontSize = '8';
add(sec1,p);



subsec1 = Section;
subsec1.Title = ...
    'Size of Evaluated Area - Added for Statistical Significance of Results';
f1 = figure('visible','off','Position',[30 30 800 250]);
hold on;
for i = 1:number_of_files
    if skipped_data(i) == 0 && i<=max(size(cases))
        h = bar(i,area_evaluated_abs_arr(i)/10000);
        set(h,'FaceColor',colors(cases(i)+1,:));
        hold on;
    end
end
xlabel(Dataset)
ylabel('Area [m^2]')
title(['Size of Evaluated Area for ' Dataset])
add(subsec1, Figure(f1));
add(sec1, subsec1);

subsec2 = Section;
subsec2.Title = ...
    'Total Number of Detected Rocks - Added for Statistical Significance of Results';
f2 = figure('visible','off','Position', [30 30 800 250]);
hold on;
for i = 1:number_of_files
    if skipped_data(i) == 0 && i<=max(size(cases))
        h = bar(i,n_rocks_arr(i));
        set(h,'FaceColor',colors(cases(i)+1,:));
        hold on;
    end
end
xlabel(Dataset)
ylabel('Number of rocks')
title(['Total Number of Detected Rocks for ' Dataset])
add(subsec2, Figure(f2));
add(sec1, subsec2);

subsec2_1 = Section;
subsec2_1.Title = 'Nodata surface area percentage';
f1 = figure('visible','off','Position',[50 50 1200 650]);
hold on;
for i = 1:number_of_files
    if skipped_data(i) == 0 && i<=max(size(cases))
        h = bar(i,percent_shadow_arr(i) * 100);
        set(h,'FaceColor',colors(cases(i)+1,:));
        hold on;
    end
end
xlabel(Dataset)
ylabel('Percentage [%]')
title(['Nodata surface area percentage for ' Dataset])
add(subsec2_1, Figure(f1));
add(sec1, subsec2_1);

% added slope
subsec2_2 = Section;
subsec2_2.Title = 'Average slope of evaluated terrain region';
f1 = figure('visible','off','Position',[50 50 1200 650]);
hold on;
for i = 1:number_of_files
    if skipped_data(i) == 0 && i<=max(size(cases))
        h = bar(i,slope_arr(i));
        set(h,'FaceColor',colors(cases(i)+1,:));
        hold on;
    end
end
xlabel(Dataset)
ylabel('Slope []')
title(['Average slope for ' Dataset])
add(subsec2_2, Figure(f1));
add(sec1, subsec2_2);

add(rpt, sec1);




sec2 = Section;
sec2.Title = 'Results';

subsec21 = Section;
subsec21.Title = 'Number of Rocks Per Surface Area';
f3 = figure('visible','off','Position',[50 50 1200 650]);
hold on;
for i = 1:number_of_files
    if skipped_data(i) == 0 && i<=max(size(cases))
        h = bar(i,area_evaluated_abs_arr(i).\ n_rocks_arr(i) * 10000);
        set(h,'FaceColor',colors(cases(i)+1,:));
        hold on;
    end
end
xlabel(Dataset)
ylabel('Number of rocks per m^2')
title(['Number of Rocks Per Surface Area for ' Dataset])
add(subsec21, Figure(f3));
add(sec2, subsec21);

subsec22 = Section;
subsec22.Title = 'Percentage of Surface Area Covered by Rocks';
f4 = figure('visible','off','Position',[50 50 1200 650]);
hold on;
for i = 1:number_of_files
    if skipped_data(i) == 0 && i<=max(size(cases))
        h = bar(i,percent_rocks_arr(i) * 100);
        set(h,'FaceColor',colors(cases(i)+1,:));
        hold on;
    end
end
% rock cells per non rock cells
xlabel(Dataset)
ylabel('Percentage [%]')
title(['Percentage of Surface Area Covered by Rocks for ' Dataset])
add(subsec22, Figure(f4));
add(sec2, subsec22);

subsec23 = Section;
subsec23.Title = 'Rock Size Distributions';
smaller15_1 = zeros(1,number_of_files);
smaller50_1 = zeros(1,number_of_files);
smaller50_2 = zeros(1,number_of_files);
larger50_1 = zeros(1,number_of_files);
for i = 1 : number_of_files
    if skipped_data(i) == 0
        smaller15_1(i) = nnz(minoraxisl_arr{i}<5)/n_rocks_arr(i);
        smaller50_1(i) = nnz(minoraxisl_arr{i}<10)/n_rocks_arr(i) -...
            smaller15_1(i);
        smaller50_2(i) = nnz(minoraxisl_arr{i}<15)/n_rocks_arr(i) -...
            smaller15_1(i)-smaller50_1(i);
        larger50_1(i) = nnz(minoraxisl_arr{i}>15)/n_rocks_arr(i);
    end
end
f5 = figure('visible','off','Position',[50 50 1200 650]);
bar(cat(1, smaller15_1 *100, smaller50_1 *100,smaller50_2 *100,...
    larger50_1 *100)')
legend('< 5cm','5 - 10cm','10 - 15cm','> 15cm')
legend('show');
ylabel('Percentage [%] of rocks')
xlabel(Dataset)
title(['Minoraxis length measurements for ' Dataset])
add(subsec23, Figure(f5));



smaller15_1 = zeros(1,number_of_files);
smaller50_1 = zeros(1,number_of_files);
smaller50_2 = zeros(1,number_of_files);
larger50_1 = zeros(1,number_of_files);
for i = 1 : number_of_files
    if skipped_data(i) == 0
        smaller15_1(i) = nnz(equivd_arr{i}<5)/n_rocks_arr(i);
        smaller50_1(i) = nnz(equivd_arr{i}<10)/n_rocks_arr(i) -...
            smaller15_1(i);
        smaller50_2(i) = nnz(equivd_arr{i}<15)/n_rocks_arr(i) -...
            smaller15_1(i)-smaller50_1(i);
        larger50_1(i) = nnz(equivd_arr{i}>15)/n_rocks_arr(i);
    end
end
f5_1 = figure('visible','off','Position',[50 50 1200 650]);
bar(cat(1, smaller15_1 *100, smaller50_1 *100, smaller50_2 *100,...
    larger50_1 *100)')
legend('< 5cm','5 - 10cm', '10 - 15cm','> 15cm')
legend('show');
ylabel('Percentage [%] of rocks')
xlabel(Dataset)
title(['Equivalent rock diameters for ' Dataset])
add(subsec23, Figure(f5_1));
add(sec2, subsec23);



subsec23_2 = Section;
subsec23_2.Title = 'Rock Area Distributions';
smaller5_1 = zeros(1,number_of_files);
smaller15_1 = zeros(1,number_of_files);
smaller25_1 = zeros(1,number_of_files);
larger50_1 = zeros(1,number_of_files);
for i = 1 : number_of_files
    if skipped_data(i) == 0
        smaller5_1(i) = nnz(rock_area_arr{i}<25)/n_rocks_arr(i);
        smaller15_1(i) = nnz(rock_area_arr{i}<100)/n_rocks_arr(i) -...
            smaller5_1(i);
        smaller25_1(i) = nnz(rock_area_arr{i}<200)/n_rocks_arr(i) -...
            smaller15_1(i) - smaller5_1(i);
        larger50_1(i) = nnz(rock_area_arr{i}>200)/n_rocks_arr(i);
    end
end
f5_2 = figure('visible','off','Position',[50 50 1200 650]);
bar(cat(1, smaller5_1 *100, smaller15_1 *100, smaller25_1 *100,...
    larger50_1 *100)')
legend('< 25cm^2','25 - 100cm^2', '100 - 200cm^2', '> 200cm^2')
legend('show');
ylabel('Percentage [%] of rocks')
xlabel(Dataset)
title(['Rock Area Distributions for ' Dataset])
add(subsec23_2, Figure(f5_2));
add(sec2, subsec23_2);






subsec24 = Section;
subsec24.Title = 'Rock Height Distributions';
smaller5_1 = zeros(1,number_of_files);
smaller15_1 = zeros(1,number_of_files);
smaller25_1 = zeros(1,number_of_files);
larger50_1 = zeros(1,number_of_files);
for i = 1 : number_of_files
    if skipped_data(i) == 0
        smaller5_1(i) = nnz(heights_arr{i}<0.05)/nnz(heights_arr{i});
        smaller15_1(i) = nnz(heights_arr{i}<0.10)/nnz(heights_arr{i}) -...
            smaller5_1(i);
        smaller25_1(i) = nnz(heights_arr{i}<0.15)/nnz(heights_arr{i}) -...
            smaller15_1(i) - smaller5_1(i);
        larger50_1(i) = nnz(heights_arr{i}>0.15)/nnz(heights_arr{i});
    end
end
f6 = figure('visible','off','Position',[50 50 1200 650]);
bar(cat(1, smaller5_1 *100, smaller15_1 *100, smaller25_1 *100,...
    larger50_1 *100)')
legend('< 5cm','5 - 10cm', '10 - 15cm', '> 15cm')
legend('show');
ylabel('Percentage [%] of rocks')
xlabel(Dataset)
title(['Rock Height Distributions for ' Dataset])
add(subsec24, Figure(f6));
add(sec2, subsec24);



% CFA plot
subsec25 = Section;
subsec25.Title = 'Cumulative Fractional Area';
f6 = figure('visible','off','Position',[50 50 1200 650]);
for i = 1:round(number_of_files/2)
    if skipped_data(i) == 0
        if size(cfa_value_arr{i},2) ~= 0
            loglog(D,cfa_value_arr{i}, 'DisplayName',...
                GeomapAnalysisfiles(i).name(1:7))        
            hold on;
        end
    end
end
xlabel('Rock Diameter [m]')
ylabel('cumulative fractional area')
title(['Cumulate Fractional Area for ' Dataset ' Part 1']) 
legend('show');
add(subsec25, Figure(f6));

f7 = figure('visible','off','Position',[50 50 1200 650]);
for i = round(number_of_files/2)+1 : number_of_files
    if skipped_data(i) == 0
        if size(cfa_value_arr{i},2) ~= 0
            loglog(D,cfa_value_arr{i}, 'DisplayName',...
                GeomapAnalysisfiles(i).name(1:7))        
            hold on;
        end
    end
end
xlabel('Rock Diameter [m]')
ylabel('Cumulative fractional area')
title(['Cumulate Fractional Area for ' Dataset ' Part 2']) 
legend('show');
add(subsec25, Figure(f7));
add(sec2, subsec25);

add(rpt, sec2);

close(rpt)

rptview(rpt)



%% Combining with orbital analysis
load('/Users/noahzr/Desktop/sols/GeomapSolRegions/Sorted_case_matrices.mat');
% contains 'final_case_matrix500','cases500','cases500sorted', .. etc. 
load('/Users/noahzr/Desktop/sols/orbitalpart_case_matrix.mat');
geoms = ({'Flat','Peak','Ridge','Shoulder','Spur','Slope',...
    'Pit','Valley','Footslope','Hollow','Unsorted'});


cases_unique = unique(cases);       % cases covered in dataset
rocks_per_a = (area_evaluated_abs_arr.\ n_rocks_arr * 10000);   
rocks_per_a_mean = zeros(length(cases_unique));
rocks_per_a_std = zeros(length(cases_unique));

figure('Name', ['Number of Rocks Per Surface Area for ' Dataset]);
hold on;
for i = 1:length(cases_unique)
    if i<=max(size(cases))
        
    idx = (cases == cases_unique(i));
    idx(skipped_data ~= 0) = 0;
    rocks_per_a_mean(i) = mean(rocks_per_a(idx(1:length(rocks_per_a))),...
        'omitnan');
    rocks_per_a_std(i) = std(rocks_per_a(idx(1:length(rocks_per_a))),...
        'omitnan');
    h = bar(i,rocks_per_a_mean(i), 'FaceColor',colors(cases_unique(i)+1,:),...
        'EdgeColor','k','LineWidth',1,'BarWidth',0.5);
    
    hold on;
    
    h1 = bar(i, rocks_per_a_std(i),'FaceColor',colors(cases_unique(i)+1,:),...
        'FaceAlpha',0.5,'EdgeColor','k','LineWidth',5);
    
    hold on;
    
    text(h1.XEndPoints, h1.YEndPoints, string(nnz(idx)),...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    legend('mean','std');
    
    end
end
legend('show');
ylabel('Number of rocks per m^2')
xlabel(['Covered landforms in ' Dataset])
xticklabels(geoms(cases_unique+1))
xticks(1:length(cases_unique));
title(['Number of Rocks Per Surface Area for ' Dataset])



percent_rocks_mean = zeros(length(cases_unique));
percent_rocks_std = zeros(length(cases_unique));

figure('Name', ['Percentage of Surface Area Covered by Rocks for ' Dataset]);
hold on;
for i = 1:length(cases_unique)
    
    idx = (cases == cases_unique(i));
    idx(skipped_data ~= 0) = 0;
    percent_rocks_mean(i) = mean(percent_rocks_arr(idx(1:...
        length(percent_rocks_arr)))*100, 'omitnan');
    percent_rocks_std(i) = std(percent_rocks_arr(idx(1:...
        length(percent_rocks_arr)))*100, 'omitnan');
    h = bar(i,percent_rocks_mean(i), 'FaceColor',...
        colors(cases_unique(i)+1,:), 'EdgeColor','k', 'LineWidth', 1,...
        'BarWidth', 0.5);
    
    hold on;
    h1 = bar(i, percent_rocks_std(i),'FaceColor', colors(cases_unique(i)+1,:),...
        'FaceAlpha',0.5,'EdgeColor','k','LineWidth',5);
    
    hold on;
    text(h1.XEndPoints, h1.YEndPoints, string(nnz(idx)),'HorizontalAlignment',...
        'center','VerticalAlignment','bottom');
    legend('mean','std');
    
end
% rock cells per non rock cells
xlabel(['Covered landforms in ' Dataset])
ylabel('Percentage [%]')
xticklabels(geoms(cases_unique+1))
xticks(1:length(cases_unique));
title(['Percentage of Surface Area Covered by Rocks for ' Dataset])



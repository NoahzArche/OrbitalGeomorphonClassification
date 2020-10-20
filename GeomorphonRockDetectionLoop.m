% running for all NC files

NCfiles = dir('../*.TIF');
number_of_files = length(NCfiles);

% output file for loop
dlmwrite(['NCLoop' num2str(date) '.csv'], {'Looping through sols'},...
    'delimiter', '');
dlmwrite(['NCLoop' num2str(date) '.csv'], date, 'delimiter', '', '-append');


for i = 1:number_of_files
    
    NCfile = NCfiles(i).name;
    NCfile = ['../' NCfile];
    try
        GeomorphonAsRockDetector4(NCfile);
        
    catch ME
        dlmwrite(['NCLoop' num2str(date) '.csv'],...
            'failed NavCam mosaic DEMs: ', 'delimiter', '', '-append');
        dlmwrite(['NCLoop' num2str(date) '.csv'], NCfiles(i).name,...
            'delimiter', '', '-append');
        dlmwrite(['NCLoop' num2str(date) '.csv'], ME.message, 'delimiter',...
            '', '-append');

    end
end

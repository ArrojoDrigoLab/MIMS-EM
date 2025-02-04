function [pathstruc] = CreateMaskDir_AllMice_Final(folderpath)
% Input is the folder containing the folders for the mice, "folderpath".
% *NOTE* the mice folders must be the only files/folders in the folder
% containing the "mouse" for this version. FURTHERMORE, this version looks
% for the files to use all the way in the "Final_Outputs" folder, as
% opposed to the folder nested one level above containing the
% "Final_Outputs" folder

% Output is a structure. First structures are mice, then hour, then cell.
% Inside the cell structures, is a mat containing the file paths.


folderconts = deblank(string(ls(folderpath)));
folderconts(folderconts == '.' | folderconts == '..') = [];

%Changes needed to work in CodeOcean
% folderconts = string([]);
% predir = dir(folderpath);
% for pd = 1:size(predir,1)
%     folderconts = [folderconts; string(predir(pd).name)];
% end
% folderconts(folderconts == '.' | folderconts == '..') = [];

%ensure only folders containing 'mouse' are matched 
folderlogic = arrayfun(@(x) regexpi(x,'mouse'),folderconts,'UniformOutput',false);
folderlogic = ~cellfun(@isempty,folderlogic);


folderconts = folderconts(folderlogic);

pathstruc = struct();
for p = 1:numel(folderconts)
    pathstruc.(folderconts(p)) = struct();
    mousepath = folderpath+filesep+folderconts(p);
    
    mouseconts = deblank(string(ls(mousepath)));
    mouseconts(mouseconts == '.' | mouseconts == '..') = [];
    
    %Changes for CodeOcean
%     mouseconts = string([]);
%     predir = dir(mousepath);
%     for pd = 1:size(predir,1)
%         mouseconts = [mouseconts;string(predir(pd).name)];
%     end
%     mouseconts(mouseconts == '.' | mouseconts == '..') = [];
    
    for r = 1:numel(mouseconts)
        pathstruc.(folderconts(p)).(mouseconts(r)) = struct(); 
        hourpath = mousepath + filesep + mouseconts(r);
        
        hourconts = deblank(string(ls(hourpath)));
        hourconts(hourconts == '.' | hourconts == '..') = [];
        
        %Changes for CodeOcean
%         hourconts = string([]);
%         predir = dir(hourpath);
%         for pd = 1:size(predir,1)
%             hourconts = [hourconts;string(predir(pd).name)];
%         end
%         hourconts(hourconts == '.' | hourconts == '..') = [];
        
        for s = 1:numel(hourconts)
            pathstruc.(folderconts(p)).(mouseconts(r)).(hourconts(s)) = struct();
            cellpath = hourpath + filesep + hourconts(s);
            cellfinoutpath = cellpath + filesep + 'Final_Outputs';
            
            cellfinoutconts = deblank(string(ls(cellfinoutpath)));
            cellfinoutconts(cellfinoutconts == '.' | cellfinoutconts == '..') = [];
            
            %Changes for CodeOcean
%             cellfinoutconts = string([]);
%             predir = dir(cellfinoutpath);
%             for pd = 1:size(predir,1)
%                 cellfinoutconts = [cellfinoutconts;string(predir(pd).name)];
%             end
%             cellfinoutconts(cellfinoutconts == '.' | cellfinoutconts == '..') = [];
            
            
            pathstruc.(folderconts(p)).(mouseconts(r)).(hourconts(s)) = cellfinoutpath + filesep + cellfinoutconts;
            
        end
    end
end





end


clearvars
% Script requires a mito mask and a c13 mask at the least
foldpath = "D:\Aliyah\Final_Spatial_Liver_Cache";

imborderthresh = 5;

%Note that the cell masks in Mouse 1 are INVERTED. Treat as such for mouse
%1 codes.




%Changes from v4 to v5:


%In v4 and earlier, contact pixels were treated as image objects and thus were
%listed via their centroids. In this version, they are listed by all their
%individual pixel coordinates instead of combined into a centroid coordinate



%Also cleans up varnames and sumvarnames. No longer stows C13 and X,Y into
%them.

%This version expects an imported cytosol mask


%%
orgnames = {'Mito', 'ER', 'LD', 'Gly', 'Nuc'};
for on1 = 1:length(orgnames)
    for on2 = 1:length(orgnames)
        regstruc.mouse_1.(strcat('Contact',orgnames{on1},'_',orgnames{on2})) = strcat('Contacts_',orgnames{on1},'_to_',orgnames{on2});
        regstruc.mouse_1.(strcat('SumContact',orgnames{on1},'_',orgnames{on2})) = strcat('Summary_Table_Contact_Test_',orgnames{on1},'_to_',orgnames{on2});
    end
end

coninputs = fieldnames(regstruc.mouse_1);

regstruc.mouse_2 = regstruc.mouse_1;


%Organelles

regstruc.mouse_1.Mito = 'MITO_processed';
regstruc.mouse_1.ER = 'EndoR_processed_Final';
regstruc.mouse_1.Nucleus = 'NUCLEUS_liver';
regstruc.mouse_1.LD = 'LD_processed';
regstruc.mouse_1.Gly = 'Glycogen_processed';
regstruc.mouse_1.Cell = 'CELL_EM_Liver';
regstruc.mouse_1.C13 = 'X13C';
regstruc.mouse_1.Cytosol = 'Cytosol_processed_EM';

regstruc.mouse_2.Mito = 'Mito_processed';
regstruc.mouse_2.ER = regstruc.mouse_1.ER;
regstruc.mouse_2.Nucleus = 'nuclei.tif';
regstruc.mouse_2.LD = regstruc.mouse_1.LD;
regstruc.mouse_2.Gly = regstruc.mouse_1.Gly;
regstruc.mouse_2.Cell = 'Mouse.\d+.Liver.\d+hr.ROI.\d+.{0,4}tif';
regstruc.mouse_2.C13 = '13C';
regstruc.mouse_2.Cytosol = 'Cytosol_processed_EM';

regstruc.mouse_3 = regstruc.mouse_2;

orginputs = {'Mito','ER','Nucleus','LD','Gly','Cell','Cytosol','C13'};


%if for some reason you still cant match mito to er types, workaround is to
%use the the contact images and equalities 




%%

maskpath = CreateMaskDirectory_AllMice(foldpath);
struc3 = strel("square",3);
strucd = strel("diamond", 1);

errorlog = string([]);
mnames = fieldnames(maskpath);
for mnum = 3:numel(mnames) %editted
    hnames = fieldnames(maskpath.(mnames{mnum}));
    for hnum = 2:numel(hnames) %edit
        cnames = fieldnames(maskpath.(mnames{mnum}).(hnames{hnum}));
        parfor cnum = 1:numel(cnames) %par this level
%             try
                pathstruc = struct();
                
                workcell = maskpath.(mnames{mnum}).(hnames{hnum}).(cnames{cnum});
                if ~isempty(workcell)
                    for k = 1:numel(workcell)
                        if regexpi(workcell{k},'Final_Outputs') 
                            finalfold = workcell{k};
                        else
                            for oname = 1:numel(orginputs)
                                if regexpi(workcell{k},regstruc.(mnames{mnum}).(orginputs{oname}))
                                    pathstruc.(orginputs{oname}) = workcell{k};
                                end
                            end
                        end
                    end

                    finalfoldconts = deblank(string(ls(finalfold)));
                    finalfoldconts( finalfoldconts == '.' | finalfoldconts == '..') = [];

                    for k = 1:numel(finalfoldconts)
                        for conname = 1:numel(coninputs)
                            if regexpi(finalfoldconts{k},regstruc.(mnames{mnum}).(coninputs{conname}))
                                pathstruc.(coninputs{conname}) = strcat(finalfold, filesep, finalfoldconts{k});
                            end
                        end
                    end
                    
                    %Label nonexistant masks with 'none' in pathstruc
                    fullinputs = [orginputs coninputs'];
                    for checkiter = 1:numel(fullinputs)
                        try
                            pathstruc.(fullinputs{checkiter});
                        catch
                            pathstruc.(fullinputs{checkiter}) = 'none';
                        end
                    end

                    imstruc = struct();
                    %Preimport masks
                    for k = 1:numel(orginputs)
                        if strcmp(pathstruc.(orginputs{k}), 'none') == 0
                            if strcmp(orginputs{k},'C13')
                                imstruc.(orginputs{k}) = imread(pathstruc.(orginputs{k}));
                            else
                                imstruc.(orginputs{k}) = logical(imread(pathstruc.(orginputs{k})));
                            end
                        end
                    end
                    
                    
                    % Add check that all masks are the same size (compares the ith and ith+1 entry for size)
                    breakflag = 0;
                    checkfields = fieldnames(imstruc);
                    for mci = 1 : numel(checkfields)-1
                        if ~all(size(imstruc.(checkfields{mci})) == size(imstruc.(checkfields{mci+1}))) %if masks are not same size, all is false
                            breakflag = 1;
                            continue
                        end
                    end
                    if breakflag == 1
                        errorlog = [errorlog;[checkfields{mci},' and ',checkfields{mci+1},' masks not same size ',mnames{mnum},' ',hnames{hnum},' ',cnames{cnum}]];
                        continue
                    end
                    
                    
                    
                    sumvarnames = [];
                    
                    %Final seems unused, potential to remove
                    
                    %Use C13 for coordinate gridding, set background levels
                    if strcmp(pathstruc.C13,'none') == 0
                        %set background levels to 102
                        imstruc.C13(imstruc.C13 < 102) = 102;

                        %Apply mean filter, filter size 17
                        imstruc.C13 = imboxfilt(imstruc.C13,17);

                        
                    end
                    

                    varnames = string([]);
                    subimage = [];

                    

                    



                    %Import Contact Images Here

                    contactcount = 0;
                    ObjCell = {};
                    conarr = fieldnames(pathstruc);
                    %Looks for text starting with "Contact", followed by any non
                    %whitespace characters
                    conlogic = cellfun(@(x) regexpi(x,'^Contact.*'),conarr,'UniformOutput',false);
                    conlogic = cellfun(@(x) ~isempty(x),conlogic, 'UniformOutput',false);
                    conarr = conarr(cell2mat(conlogic));

                    for connum = 1:numel(conarr)
                        if ~strcmp(pathstruc.(conarr{connum}),'none')
                            tempcIM = imread(pathstruc.(conarr{connum}));
                            tempcARR = cast(tempcIM(:),'double');
                            tempstow = SquidPreprocStandard(tempcIM,imstruc.C13,imborderthresh);
                            if ~isempty(tempstow)
                                ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                            else
                                tempstow{1} = {};
                                ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                            end
                            varnames = cat(2,varnames,conarr{connum});
                            contactcount = contactcount+1;
                        end
                    end



                    %Main Batch Starts Here
                    
                    
                    %First Pass Subtraction
                    SIvarnamesFP = string([]);
                    
                    
                    if isfield(imstruc,'Cell')
                        imstruc.Cell(imstruc.Cell > 0 ) = 1;
                        imstruc.Cell = logical(bwfill(imstruc.Cell, 'holes'));
                        subimage(:,:,1) = ~imstruc.Cell;
                        SIvarnamesFP = [SIvarnamesFP,'Cell']; 
                    end
                    if isfield(imstruc,'Nucleus')
                        tempsub = imstruc.Nucleus;
                        for subid = 1:size(subimage,3)
                            tempsub(tempsub > 0 & subimage(:,:,subid) > 0) = 0;
                        end
                        subimage = cat(3,subimage,tempsub);
                        imstruc.Nucleus = tempsub;
                        SIvarnamesFP = [SIvarnamesFP,'Nucleus'];
                    else
                        tempsub = zeros(size(imstruc.C13));
                        imstruc.Nucleus = tempsub;
                        subimage = cat(3,subimage,tempsub);
                        SIvarnamesFP = [SIvarnamesFP,'Nucleus'];
                    end
                    if isfield(imstruc,'Mito')
                        tempsub = imstruc.Mito;
                        for subid = 1:size(subimage,3)
                            tempsub(tempsub > 0 & subimage(:,:,subid) > 0) = 0;
                        end
                        %fill holes
                        tempsub = bwfill(tempsub,'holes');
                        subimage = cat(3,subimage,tempsub);
                        imstruc.Mito = tempsub;
                        SIvarnamesFP = [SIvarnamesFP,'Mito'];
                    end
                    %Dont rewrite imstruc after this, will save in second pass
                    %for the ones not saved in the first pass
                    if isfield(imstruc,'ER')
                        tempsub = imstruc.ER;
                        %Erode tempsub 4 times
                        for eroid = 1:4
                            tempsub = imerode(tempsub,struc3);
                        end
                        %Subtract
                        for subid = 1:size(subimage,3)
                            tempsub(tempsub > 0 & subimage(:,:,subid) > 0) = 0;
                        end
                        %remove isolated pixels
                        tempsub = bwmorph(tempsub,'clean');
                        subimage = cat(3,subimage,tempsub);
                        SIvarnamesFP = [SIvarnamesFP, 'ER'];
                    end
                    if isfield(imstruc,'Gly')
                        tempsub = imstruc.Gly;

                        for subid = 1:size(subimage,3)
                            tempsub(tempsub > 0 & subimage(:,:,subid) > 0) = 0;
                        end
                        subimage = cat(3,subimage,tempsub);
                        SIvarnamesFP = [SIvarnamesFP,'Gly'];
                    end
                    if isfield(imstruc,'LD')
                        tempsub = imstruc.LD;

                        %dilate twice, ask 4bit or 8bit connectivity
                        for dilid = 1:2
                            tempsub = imdilate(tempsub, strucd);
                        end

                        for subid = 1:size(subimage,3)
                            tempsub(tempsub > 0 & subimage(:,:,subid) > 0) = 0;
                        end
                        tempsub = bwfill(tempsub,'holes');
                        subimage = cat(3,subimage,tempsub);
                        SIvarnamesFP = [SIvarnamesFP,'LD'];
                    end
                    %Treat Cytosol Mask using Organelle masks
                    % (((Cyt - Nuc) - Mito ) - LD)

                    if isfield(imstruc,'Cytosol')
                        imstruc.Cytosol(imstruc.Cytosol > 0 ) = 1;

                        if isfield(imstruc,'Nucleus')
                            imstruc.Cytosol(imstruc.Cytosol > 0 & imstruc.Nucleus > 0) = 0;
                        end
                        if isfield(imstruc,'Mito')
                            imstruc.Cytosol(imstruc.Cytosol > 0 & imstruc.Mito > 0) = 0;
                        end
                        if isfield(imstruc,'LD')
                            imstruc.Cytosol(imstruc.Cytosol > 0 & imstruc.LD > 0) = 0;
                        end

                    else
                        errorlog = [errorlog;['Cytosol Mask Missing ',mnames{mnum},hnames{hnum},cnames{cnum}]];
                        continue
                    end
                    
                    %Reset subimage for second pass
                
                
                    %Remove parts of subimage that need to be retreated as a
                    %result of cyto mask being upstream (anything after ER)

                    %Find position of ER in subimage using SIvarnamesFP (sub 1,
                    %index to keep)



                    subimage = subimage(:,:,1:(find(SIvarnamesFP == 'ER') - 1) );


                    %Second Subtraction Pass
                    
                    
                    
                    if isfield(imstruc,'Cell')
                        %Can just remove Cell periphery if not being used for now.

                        %IF RE-ADDING CELL PERIPH, NOTE THAT NEED TO REMOVE PIXELS
                        %ALONG IMAGE BORDER if they are not true cell periphery

                        %Previously had used PYTHON for that part.

        %                 imCell = bwmorph(imCell, 'remove');
        %                 Cellarray = cast(imCell(:),'double');

        %                 cellid = find(Cellarray);
        %                 tempstow = {};
        %                 for t = 1:numel(cellid)
        %                    Cellarray(cellid(t)) = t;
        %                    tempstow{t} = [Final(cellid(t),[1 2]), C13array(cellid(t))];
        %                 end
        %                 if ~isempty(tempstow) 
        %                     ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
        %                 end
        %                 if isempty(tempstow)
        %                    tempstow{1} = {};
        %                    ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
        %                 end
        %                 varnames=cat(2,varnames,'Cell_Periphery');
                    end
                    if isfield(imstruc,'Nucleus')
                        imstruc.Nucleus = bwmorph(imstruc.Nucleus, 'remove');
                        imstruc.Nucleus = bwmorph(imstruc.Nucleus, 'skeleton', Inf);
                        imstruc.Nucleus = bwmorph(imstruc.Nucleus, 'clean');
                        %Creates one cell for nuclear envelope coords and
                        %C13 val
                        prestownuc = SquidPreprocNoThresh(imstruc.Nucleus,imstruc.C13);
                        if ~isempty(prestownuc)
                            %Convert the single cell into one cell per
                            %coord
                            for nucid = 1:size(prestownuc{1},1)
                                tempstow{nucid} = prestownuc{1}(nucid,:);
                            end
                            ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        else
                           tempstow{1} = {};
                           ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        end
                        varnames=cat(2,varnames,"Nuc");
                    else
                        tempstow{1} = {};
                        ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        varnames=cat(2,varnames,"Nuc");
                    end
                    if isfield(imstruc,'Mito')        
                        tempstow = SquidPreprocStandard(imstruc.Mito,imstruc.C13,imborderthresh);
                        if ~isempty(tempstow) 
                            ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        else
                           tempstow{1} = {};
                           ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        end
                        varnames=cat(2,varnames,"Mito");
                    end
                    if isfield(imstruc,'ER')
                        tempsub = imstruc.ER;

                        %Cytosol Subtraction
                        tempsub(tempsub > 0 & imstruc.Cytosol > 0) = 0;


                        %Erode tempsub 4 times
                        for eroid = 1:4
                            tempsub = imerode(tempsub,struc3);
                        end


                        %Subtract
                        for subid = 1:size(subimage,3)
                            tempsub(tempsub > 0 & subimage(:,:,subid) > 0) = 0;
                        end
                        %remove isolated pixels
                        tempsub = bwmorph(tempsub,'clean');
                        subimage = cat(3,subimage,tempsub);
                        imstruc.ER = tempsub;
                        tempstow = SquidPreprocStandard(imstruc.ER,imstruc.C13,imborderthresh);


                        if ~isempty(tempstow) 
                            ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        else 
                           tempstow{1} = {};
                           ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        end
                        varnames = [varnames, "ER"];
                    end
                    if isfield(imstruc,'Gly')
                        tempsub = imstruc.Gly;
                        %Cytosol Subtraction
                        tempsub(tempsub > 0 & imstruc.Cytosol > 0) = 0;

                        for subid = 1:size(subimage,3)
                            tempsub(tempsub > 0 & subimage(:,:,subid) > 0) = 0;
                        end
                        subimage = cat(3,subimage,tempsub);
                        imstruc.Gly = tempsub;
                        tempstow = SquidPreprocStandard(imstruc.Gly,imstruc.C13,imborderthresh);

                        if ~isempty(tempstow) 
                            ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        else 
                           tempstow{1} = {};
                           ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        end
                        varnames = [varnames, "Gly"];
                    end
                    if isfield(imstruc,'LD')
                        tempsub = imstruc.LD;

                        %dilate twice, ask 4bit or 8bit connectivity
                        for dilid = 1:2
                            tempsub = imdilate(tempsub, strucd);
                        end

                        for subid = 1:size(subimage,3)
                            tempsub(tempsub > 0 & subimage(:,:,subid) > 0) = 0;
                        end
                        tempsub = bwfill(tempsub,'holes');
                        subimage = cat(3,subimage,tempsub);
                        imstruc.LD = tempsub;
                        tempstow = SquidPreprocStandard(imstruc.LD,imstruc.C13,imborderthresh);

                        if ~isempty(tempstow) 
                            ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                            varnames = [varnames, "LD"];
                        else 
                           tempstow{1} = {};
                           ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                           varnames = [varnames, "LD"];
                        end
                    end


                    


                    %Import Contact Summary Spreadsheets for Organelles
                    xlsSumContact = struct();
                    for on1 = 1:length(orgnames)
                        for on2 = 1:length(orgnames)
                            if strcmp(pathstruc.(strcat('SumContact',orgnames{on1},'_',orgnames{on2})),'none') == 0
                                xlsSumContact.(strcat(orgnames{on1},'_',orgnames{on2})) = readtable(pathstruc.(strcat('SumContact',orgnames{on1},'_',orgnames{on2})));
                                sumvarnames = [sumvarnames, strcat(orgnames{on1}, '_',orgnames(on2))]; %think about how this is used before, now modified 
                            end
                        end
                    end

                    interstruc = struct();
                    primarylist = cell(0,0);
                    %to determine the primary organelle, find the first
                    %underscore and take all the chars up to but not
                    %including that index


                    for intiter = 1:length(sumvarnames)
                        tempinter = sumvarnames{intiter};
                        USid = regexpi(tempinter,'_','start');
                        primorg = tempinter(1:USid-1);
                        secorg = tempinter(USid+1:end);
                        if isempty(ismember(primarylist,primorg)) | ~ismember(primarylist,primorg)
                            
                            interstruc.(primorg).Varnames = string(secorg);
                            primarylist = [primarylist,primorg];
                            interstruc.(primorg).Data = [xlsSumContact.(sumvarnames{intiter}).X, xlsSumContact.(sumvarnames{intiter}).Y, xlsSumContact.(sumvarnames{intiter}).Pixel_Contact_Count];

                        else
                            
                            interstruc.(primorg).Varnames = [interstruc.(primorg).Varnames, string(secorg)]; %fixed labeling here
                            [~,tempid] = ismember([interstruc.(primorg).Data(:,1),interstruc.(primorg).Data(:,2)],[xlsSumContact.(sumvarnames{intiter}).X,xlsSumContact.(sumvarnames{intiter}).Y],'rows');

                            matcharr = zeros(size(interstruc.(primorg).Data,1),1);
                            for eliter = 1:size(interstruc.(primorg).Data,1)
                                matchid = find(tempid == eliter);
                                matcharr(eliter) = xlsSumContact.(sumvarnames{intiter}).Pixel_Contact_Count(matchid);
                            end

                            interstruc.(primorg).Data = [interstruc.(primorg).Data,matcharr];
                        end

                    end

                    ptemp = fieldnames(interstruc);
                    for primiter = 1:length(ptemp)
                        interstruc.(ptemp{primiter}).Data(:,3:end) = double(interstruc.(ptemp{primiter}).Data(:,3:end) > 0);
                    end
                    
                    


                    strcat("Identifying ",mnames{mnum}," ",hnames{hnum}," ",cnames{cnum})



                    %Threshold by pixel size
                    Mitoobjid = find(strcmp(varnames,"Mito"));
                    ERobjid = find(strcmp(varnames,"ER"));
                    Glyobjid = find(strcmp(varnames,"Gly"));
                    LDobjid = find(strcmp(varnames,"LD"));

                    for k = 1:size(ObjCell,2)
                        if isempty(ObjCell{Mitoobjid,k}) == 0 && size(ObjCell{Mitoobjid,k},1) < 5000  
                            ObjCell{Mitoobjid,k} = [];
                        end
                        if isempty(ObjCell{ERobjid,k}) == 0 && size(ObjCell{ERobjid,k},1) < 150 
                            ObjCell{ERobjid,k} = [];
                        end
                        if isempty(ObjCell{Glyobjid,k}) == 0 && size(ObjCell{Glyobjid,k},1) < 500
                            ObjCell{Glyobjid,k} = [];
                        end
                        if isempty(ObjCell{LDobjid,k}) == 0 && size(ObjCell{LDobjid,k},1) < 5000
                            ObjCell{LDobjid,k} = [];
                        end
                    end







                    CentroidArr = [];
                    for x = 1:size(ObjCell,1)
                        for y = 1:size(ObjCell,2)
                            
                            %IF WANT TO SPEED UP:
                            %add another if statement for if within contactcount,
                            %treat cells differently on ObjCell, as well as nucleus and cell periph.
                            %Will speed up code significantly if so. Store all
                            %nucleus and cell periph pixels in a single cell
                            %instead.
                            if ~isempty(ObjCell{x,y})
                                temparr = [];
                                if x <= contactcount %if Cell is a "contact" cell 
                                    for z = 1:size(ObjCell{x,y},1)
                                        centroidX = ObjCell{x,y}(z,1);
                                        centroidY = ObjCell{x,y}(z,2);
                                        idrow = zeros(1,size(varnames,2)); %Creates spacing for organelle binary identifier
                                        idrow(1,x) = 1; %Labels organelle identifier
                                        AvgC13 = ObjCell{x,y}(z,3);
                                        PixelCount = 1;
                                        perim = 1;
                                        CurveScore = 0;


                                        %C13 Binning
                                        C13Bin = BinC13(AvgC13);
                                        
                                        OrgPlusC13Bin = strcat(varnames(x),"_",C13Bin);
                                        
                                        TempArr = [round(centroidX), round(centroidY), 1, idrow, AvgC13, PixelCount, perim, CurveScore, varnames(x), C13Bin, OrgPlusC13Bin];
                                        CentroidArr = [CentroidArr; TempArr];
                                    end
                                else
                                    centroidX = sum(ObjCell{x,y}(:,1))/size(ObjCell{x,y},1);
                                    centroidY = sum(ObjCell{x,y}(:,2))/size(ObjCell{x,y},1);
                                    idrow = zeros(1,size(varnames,2)); %Creates spacing for organelle binary identifier
                                    idrow(1,x) = 1;
                                    AvgC13 = sum(ObjCell{x,y}(:,3))/size(ObjCell{x,y},1);
                                    PixelCount = size(ObjCell{x,y},1);
                                    tempspace = logical(zeros(size(imstruc.Mito,1),size(imstruc.Mito,2)));
                                    for z = 1:size(ObjCell{x,y},1)
                                        tempspace(ObjCell{x,y}(z,1),ObjCell{x,y}(z,2)) = 1;
                                    end
                                    tempspace = logical(tempspace);
                                    tempspace = bwmorph(tempspace, 'remove');
                                    perim = numel(find(tempspace));
                                    CurveScore = (4*pi*PixelCount)/perim.^2;

                                    %C13 Binning
                                    C13Bin = BinC13(AvgC13);


                                    %The if statement below is added to label types of
                                    %mitochondria based on contacts

                                    %The if statement below is added to label types of
                                    %orgs based on contacts

                                    porglist = fieldnames(interstruc);
                                    if ~isempty(find(ismember(orgnames,varnames(x)))) & ~isempty(find(ismember(porglist,varnames(x))))
                                        [~,OrgMatchID] = ismember([round(centroidX),round(centroidY)],interstruc.(varnames{x}).Data(:,1:2),"rows");

                                        OrgConString = "";
                                        if OrgMatchID ~= 0
                                            BoolContactsArr = logical(interstruc.(varnames{x}).Data(OrgMatchID,3:end));

                                            for w = 1:size(BoolContactsArr,2)
                                                if BoolContactsArr(w)
                                                    OrgConString = strcat(OrgConString,"_",interstruc.(varnames{x}).Varnames{w});

                                                end
                                            end

                                        end

                                        if strcmp(OrgConString,"")
                                            OrgPlusC13Bin = strcat(varnames(x),"_",C13Bin);
                                            TempArr = [round(centroidX), round(centroidY), 1, idrow, AvgC13, PixelCount, perim, CurveScore, varnames(x), C13Bin, OrgPlusC13Bin];

                                        else
                                            OrgPlusC13Bin = strcat(varnames(x),OrgConString,"_",C13Bin);
                                            TempArr = [round(centroidX), round(centroidY), 1, idrow, AvgC13, PixelCount, perim, CurveScore, strcat(varnames(x),OrgConString), C13Bin, OrgPlusC13Bin];

                                        end
                                    else
                                        OrgPlusC13Bin = strcat(varnames(x),"_",C13Bin);
                                        TempArr = [round(centroidX), round(centroidY), 1, idrow, AvgC13, PixelCount, perim, CurveScore, varnames(x), C13Bin, OrgPlusC13Bin];
                                    end


                                    CentroidArr = [CentroidArr; TempArr];
                                    
                                end

                            end
                         end

                   end







                CentroidTable = array2table(CentroidArr,"VariableNames",["X", "Y", "Z",varnames(1,:),"Avg_C13", "PixelCount", "Perimeter_Pixel_Count","Curvature_Score", "Organelle", "C13_Bin", "Organelle_and_C13_Bin"]);










                % get uniques of each column, match to varnames index
                % or you can do it generically, then concat the number after

                % basic way: if organelle slot is empty and indicator is nonzero

                % alt way try first: for each row of finaldata, check each column, if nonzero store col number and ref varnames













                fullpath = strcat(foldpath,filesep,mnames{mnum},filesep,hnames{hnum},filesep,cnames{cnum},filesep,'Final_Outputs',filesep, 'Data_for_Squidpy_Centroids_', mnames{mnum},"_",hnames{hnum},"_",cnames{cnum}, '.csv');
                fullpath2 = strcat("D:\Aliyah\Final_Spatial_Liver_SummaryGraphs\Squidpy_Outputs",filesep,'Data_for_Squidpy_Centroids_', mnames{mnum},"_",hnames{hnum},"_",cnames{cnum}, '.csv');
                writetable(CentroidTable,fullpath,'WriteMode','overwrite');
                writetable(CentroidTable,fullpath2);

                 
               

                end
                
%             catch
%                 error(['Error in File ' strcat('mnum',num2str(mnum),' roinum',num2str(hnum),' cnum',num2str(cnum))])
%             end
        end
    end
end

display(errorlog)


%% Local Functions

function Binout = BinC13(AC13)
    if AC13 < 120
        Binout = "Low";
    end
    if AC13 >= 120 & AC13 < 140
       Binout = "MedLow"; 
    end
    if AC13 >= 140 & AC13 < 160
        Binout = "Med";
    end
    if AC13 >= 160 & AC13 <180
        Binout = "MedHigh";
    end
    if AC13 >= 180 & AC13 < 200
        Binout = "High";
    end
    if AC13 >= 200
        Binout = "VeryHigh";
    end

end
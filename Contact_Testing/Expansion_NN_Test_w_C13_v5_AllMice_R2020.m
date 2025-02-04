clearvars

%Set folder path, should contain "Mouse" folders
mousepath = "D:\Aliyah\Final_Spatial_Liver_Cache\Test_Set";

%Set Regular expressions used to match to organelle files
regstruc.mouse_3.Mito = 'Mito_processed';
regstruc.mouse_3.ER = 'EndoR_Processed_Final';
regstruc.mouse_3.Nucleus = 'nuclei.tif';
regstruc.mouse_3.LD = 'LD_processed';
regstruc.mouse_3.Gly = 'Glycogen_processed';
regstruc.mouse_3.C13 = '13C';
regstruc.mouse_3.Cell = 'Mouse.\d+.Liver.\d+hr.ROI.\d+.{0,4}tif';
regstruc.mouse_3.Cytosol = 'Cytosol_processed_EM';

regstruc.mouse_2 = regstruc.mouse_3;

regstruc.mouse_1.Mito = 'MITO_processed';
regstruc.mouse_1.ER = 'EndoR_Processed_Final';
regstruc.mouse_1.Nucleus = 'NUCLEUS_liver';
regstruc.mouse_1.LD = 'LD_processed';
regstruc.mouse_1.Gly = 'Glycogen_processed';
regstruc.mouse_1.C13 = 'X13C';
regstruc.mouse_1.Cell = 'CELL_EM_Liver';
regstruc.mouse_1.Cytosol = 'Cytosol_processed_EM';



%Still need to eventually rewrite, Preproc Mask Saving should be done AFTER
%thresholding by image object size. For now, perform on image J, or use
%some 2Dbwopen tool of some sort

%Make sure the image using to search is labeled as "I"
%image being searched is labeled "SearchImage"
%and C13 mask is labeled "imstruc.C13"

%Images being compared must be the same size

%Changes from v4 -> v5, different preprocessing/image subtraction pathway(cytosol mask added).

%Also, masks are now preimported instead of progressively imported. They
%are now stored in "imstruc", rather than individual variables. Instead of
%using strcmp with pathstruc to verify if a mask is found, will use isfield
%on imstruc.

%NOTE: CELL MASK IS INVERTED FOR MOUSE 1(dont see it inverted? perhaps imageJ issue)! TREAT AS SUCH. MICE 2 AND 3 USE
%EM MASKS INSTEAD.

%Added check to ensure that all masks are the same size. Also added an
%errorlog that will list all cells that had to be skipped with the
%corresponding reasoning


maskimport = {'Mito','LD','Gly','ER','Nucleus','Cell','Cytosol','C13'}; %Mask name keys to help with calls and organization

%%
%Define variables here
maskorder = {'Mito','LD','Gly','ER','Nucleus'}; %order of masks to run

maskexcept = {'Nucleus'}; %mask indexes to include in analysis but not run test on.
                            %Can be more than one mask, make a cell array
                            %of strings. Must be included in maskorder.
                            
pixthresharr = [5000,5000,500,500,500,500,0,0]; % number of pixels to threshold objects by in maskimport (at the minimum)
                                                % A 0 indicates no
                                                % thresholding necessary.
                                                % Objects which size is
                                                % below this size
                                                % will be removed

imborderthresh = 5; %If within this # of pixels to image edge (normal from the image edge)
                    % the object will be deleted
                    
threshold =[1 2]; % Minimum boundary is 1. 1 indicates direct contact.
                    %Sets range of distances (in pixel lengths) that will pass the
                    %contact test. For now, this code can only support
                    %whole numbers as values for threshold, may add support
                    %for decimals if needed in the future


sepcontactlength = 4; %In pixels, max distance between perimeter pixels that pass the
                      %contact test to be considered the same contact site
                      %(currently unused, support for number of contact
                      %sides 
                      
ContactMask = true; %returns image masks of contact sites if true

homotypictest = true; %determines if homotypic test is to be run, will run if true (an organelle searching for neighbors of the same organelle class) 

homotypicall = true; %run homotypic test on all organelles? (unused for now, have not programmed false case)

homotypicnum = []; %organelle names to run the homotypic test on IF homotypicall IS FALSE.
%  ^^^ didnt program this yet, can ignore 

savemaskopt = true; %saves preprocessed masks to provided location below if true

masksavepath = "D:\Aliyah\Final_Spatial_Liver_Cache\Preprocessed_Masks"; %path to save masks if savemaskopt true

%% Create file directory for images and Structure elements


maskpath = CreateMaskDirectory_AllMice(mousepath); %Create Structure for file directory

%Initiialize kernels used in convolutions
struc3 = strel("square",3);
strucd = strel("diamond", 1);

%% Main Loop

errorlog = string([]); %Create empty string arr for user made error messages
mnames = fieldnames(maskpath);
for mnum = 1:numel(mnames) %iter over mice, parallel loop goes here
    hnames = fieldnames(maskpath.(mnames{mnum}));
    for hnum = 1:numel(hnames) %iter over hours
        cnames = fieldnames(maskpath.(mnames{mnum}).(hnames{hnum}));
        for cnum = 1:numel(cnames) %iter over cells
%             try
                pathstruc = struct(); %Create empty structure for image file paths
                workcell = maskpath.(mnames{mnum}).(hnames{hnum}).(cnames{cnum}); %folder contents from desired cell
                
                %fill pathstruc
                for k = 1:numel(workcell)
                    for m = 1:numel(maskimport)
                        if regexpi(workcell{k},regstruc.(strcat('mouse_',num2str(mnum))).(maskimport{m}))
                            pathstruc.(maskimport{m}) = workcell{k};
                        end

                    end
                end


                %any image key not matched to an image in the folder directory is filled with 'none' 
                for k = 1:numel(maskimport)
                    try
                        pathstruc.(maskimport{k});
                    catch
                        pathstruc.(maskimport{k}) = 'none';
                    end
                end


                imstruc = struct();
                %Import masks
                for k = 1:numel(maskimport)
                    if strcmp(pathstruc.(maskimport{k}), 'none') == 0
                        if strcmp(maskimport{k},'C13')
                            imstruc.(maskimport{k}) = imread(pathstruc.(maskimport{k}));
                        else
                            imstruc.(maskimport{k}) = logical(imread(pathstruc.(maskimport{k})));
                        end
                    end
                end
                
                % Add check that all masks are the same size (compares the ith and ith+1 entry for size)
                breakflag = 0;
                checkfields = fieldnames(imstruc);
                for mci = 1 : numel(checkfields)-1
                    if ~all(size(imstruc.(checkfields{mci})) == size(imstruc.(checkfields{mci+1}))) %if masks are not same size, all() is false
                        breakflag = 1;
                        continue
                    end
                end
                if breakflag == 1
                    %Create error message and skip this image set if all
                    %images not same size
                    errorlog = [errorlog;[checkfields{mci},' and ',checkfields{mci+1},' masks not same size ',mnames{mnum},' ',hnames{hnum},' ',cnames{cnum}]];
                    continue
                end
                

                ObjCell = {}; %ObjCell is in [X,Y]

                %C13
                if strcmp(pathstruc.C13,'none') == 0
                    %Set background c13 levels to 102
                    imstruc.C13(imstruc.C13 < 102) = 102;                

                    %Apply mean filter, filter size 17x17
                    imstruc.C13 = imboxfilt(imstruc.C13,17);
                end

                subimage = [];
                varnames = [];
                
                
                
                %First Pass Subtraction
                SIvarnamesFP = string([]); %Create an array to store varnames for first pass of subimage
                
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

                    %dilate twice, ask 4bit or (8bit)(using 8 here) connectivity
                    for dilid = 1:2
                        tempsub = imdilate(tempsub, struc3);
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

                %Skip Cell, no new steps needed
                             
                %Nucleus Treatment, doesnt need if statement because
                %earlier if statement
                tempstow = SquidPreprocNoThresh(imstruc.Nucleus,imstruc.C13);
                if ~isempty(tempstow) 
                    ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                    varnames = [varnames,"Nucleus"];
                else 
                   tempstow{1} = {};
                   ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                   varnames = [varnames,"Nucleus"];
                end   
                %Mito 
                if isfield(imstruc,'Mito')
                    tempstow = SquidPreprocStandard(imstruc.Mito,imstruc.C13,imborderthresh);
                    %Reconstruct for mask stow if savemaskopt true
                    if savemaskopt
                        reconim = false(size(imstruc.C13));
                        for reciter = 1:numel(tempstow)
                            for rowit = 1:size(tempstow{reciter},1)
                                reconim(tempstow{reciter}(rowit,2),tempstow{reciter}(rowit,1)) = 1;
                            end
                        end
                        imstruc.Mito = reconim;
                    end

                    if ~isempty(tempstow) 
                        ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        varnames = [varnames, "Mito"];
                    else 
                       tempstow{1} = {};
                       ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                       varnames = [varnames, "Mito"];
                    end
                end
                %ER
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
                    if savemaskopt
                        reconim = logical(zeros(size(imstruc.C13)));
                        for reciter = 1:numel(tempstow)
                            for rowit = 1:size(tempstow{reciter},1)
                                reconim(tempstow{reciter}(rowit,2),tempstow{reciter}(rowit,1)) = 1;
                            end
                        end
                        imstruc.ER = reconim;
                    end

                    if ~isempty(tempstow) 
                        ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        varnames = [varnames, "ER"];
                    else 
                       tempstow{1} = {};
                       ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                       varnames = [varnames, "ER"];
                    end
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
                    if savemaskopt
                        reconim = logical(zeros(size(imstruc.C13)));
                        for reciter = 1:numel(tempstow)
                            for rowit = 1:size(tempstow{reciter},1)
                                reconim(tempstow{reciter}(rowit,2),tempstow{reciter}(rowit,1)) = 1;
                            end
                        end
                        imstruc.Gly = reconim;
                    end

                    if ~isempty(tempstow) 
                        ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        varnames = [varnames, "Gly"];
                    else 
                       tempstow{1} = {};
                       ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                       varnames = [varnames, "Gly"];
                    end
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
                    if savemaskopt
                        reconim = logical(zeros(size(imstruc.C13)));
                        for reciter = 1:numel(tempstow)
                            for rowit = 1:size(tempstow{reciter},1)
                                reconim(tempstow{reciter}(rowit,2),tempstow{reciter}(rowit,1)) = 1;
                            end
                        end
                        imstruc.LD = reconim;
                    end

                    if ~isempty(tempstow) 
                        ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                        varnames = [varnames, "LD"];
                    else 
                       tempstow{1} = {};
                       ObjCell = CellArrayCatUneq(ObjCell,tempstow,1);
                       varnames = [varnames, "LD"];
                    end
                end

                % Note, varnames maps the organelles from ObjCell
                % Also, objects on image edge are already removed
                % at this point









                %Should perform object size thresholding before the
                %reconstruction
                
                %Threshold images by size
                pixthreshdict = containers.Map(maskimport,pixthresharr);

                for orgiter = 1:size(ObjCell,1)
                    for coliter = 1:size(ObjCell,2)
                        if ~isempty(ObjCell{orgiter,coliter}) && size(ObjCell{orgiter,coliter},1) < pixthreshdict(varnames{orgiter})
                            ObjCell{orgiter,coliter} = [];
                        end
                    end
                end



                % put main loop here
                strcat("Identifying ", mnames{mnum}," ",hnames{hnum}, " ", cnames{cnum})

                %for main loop, Less efficient but easiest way is to load in the rows of the
                %ObjCell needed. Since the input of the contact test function
                %is two images, will have to convert objcell back into images,
                %stow into a structure, then use the string keys to call each
                %image when needed according to test inputs. Then call
                %regionprops and bwconncomp.


                %Reconstruct images in a struc
                postimstruc = struct();

                for orgiter = 1:size(ObjCell,1)
                     postimstow = zeros(size(imstruc.C13));
                    for coliter = 1:size(ObjCell,2)
                        if ~isempty(ObjCell{orgiter,coliter})
                            for pixiter = 1:size(ObjCell{orgiter,coliter},1)
                                postimstow(ObjCell{orgiter,coliter}(pixiter,2),ObjCell{orgiter,coliter}(pixiter,1)) = 1;
                            end
                        end
                    end
                    postimstruc.(varnames(orgiter)) = postimstow;
                end

                %SAVE MASKS HERE IF savemaskopt == true

                if savemaskopt
                        savefields = fieldnames(postimstruc);
                        %save reconstructed masks after thresholding by size
                        %and by dist to boundary
                        for saveiter = 1:numel(savefields)
                            removepart = pathstruc.(savefields{saveiter});
                            if strcmp(removepart,'none') == 0
                                pathend = regexpi(removepart,'Final_Spatial_Liver_Cache\\','end');
                                fullsavem = masksavepath + removepart(pathend:end);
                                imwrite(im2uint16(postimstruc.(savefields{saveiter})), fullsavem,'WriteMode','overwrite');
                            end
                        end
                        %Save c13 mask
                        c13remove = pathstruc.C13;
                        pathendc13 = regexpi(c13remove,'Final_Spatial_Liver_Cache\\','end');
                        fullsavec13 = masksavepath + c13remove(pathendc13:end);
                        imwrite(im2uint16(imstruc.C13), fullsavec13, 'WriteMode', 'overwrite');

                        %copy over original EM mask
                        cellsave = imread(pathstruc.Cell);
                        EMremove = pathstruc.Cell;
                        pathendEM = regexpi(EMremove,'Final_Spatial_Liver_Cache\\','end');
                        fullsaveEM = masksavepath + EMremove(pathendEM:end);
                        imwrite(im2uint16(cellsave),fullsaveEM,'WriteMode', 'overwrite');
                    end


                %Set up for loops to iterate over proper organelles
                for maskiter1 = 1:numel(maskorder)
                    for maskiter2 = 1:numel(maskorder)

                        if ~ismember(maskorder(maskiter1),maskexcept) %test if current mask is in maskexcept
                            %if not in maskexcept:

                            if strcmp(maskorder{maskiter1},maskorder{maskiter2}) & homotypictest == 1  
                                %if homotypic and homotypictest == true:

                                I = postimstruc.(maskorder{maskiter1});
                                SearchImage = postimstruc.(maskorder{maskiter2});


                                CCstruc = bwconncomp(I);
                                imObjs = regionprops(CCstruc,"Area", "Centroid", "PixelList"); %cant use perimeter, gives distance instead of pixels

                                %Performing the contact test below

                                %foundcell = DilationContactTest(imObjs,SearchImage, threshold); 
                                foundcell = DistContactTestHomotypic(imObjs, SearchImage, threshold); %CHECK THIS CODE ************
                                perimarr = PerimPixCount(imObjs,size(I)); 


                                %Delete Single Pixel Contact Sites
                                foundcell = IsoPixelContactDelete(foundcell, size(I));



                                %Create contactsarr from old code, phase out eventually
                                contactsarr = [];
                                for s = 1:size(foundcell,2)
                                    if ~isempty(foundcell{s})
                                        contactsarr = [contactsarr, size(foundcell{s},1)]; 
                                    else
                                        contactsarr = [contactsarr, 0];
                                    end
                                end


                                %numcontactsarr = ContactSiteCount(foundcell, sepcontactlength);




                                %Array for percent of perimeter is in contact with search object
                                percentcontactarr=[];

                                %iterate over contactsitecell, ~ifempty, count the rows,
                                %add, reset every new column

                                for n = 1:size(contactsarr,2)

                                   if contactsarr(n) >0
                                      percentcontact = contactsarr(n)./ perimarr(n);
                                   else
                                       percentcontact = 0;
                                   end
                                   percentcontactarr = [percentcontactarr; percentcontact];
                                end






                                %Create image masks for contact sites, if ContactMask = true
                                idx = 0;
                                ContactsiteImage = logical(zeros(size(I,1),size(I,2)));
                                if ContactMask == 1;
                                    for p = 1:size(foundcell,1)
                                        for r = 1:size(foundcell,2)
                                           if ~isempty(foundcell{p,r})
                                                for s = 1:size(foundcell{p,r})
                                                    ContactsiteImage(foundcell{p,r}(s,1),foundcell{p,r}(s,2)) = 1;

                                                end
                                           end
                                        end                   
                                    end
                                end

                                ContactsiteImage = uint16(ContactsiteImage);


                                %Avg C13 and Object Size for Sumtable.
                                %Also assigns centroid of object, want in [x,y]
                                C13AvgArr = [];
                                ObjSizeArr = [];
                                CentroidArr = [];
                                for n = 1:size(imObjs,1)               
                                    C13Plist = imObjs(n).PixelList; %[x,y]=[c,r]
                                    C13R = C13Plist(:,2);
                                    C13C = C13Plist(:,1);
                                    %Find Centroid here, note, this is inverted relative to
                                    %preprocessing. Should be column then row (x,y)
                                    CentroidPlist = imObjs(n).Centroid; %[x,y]
                                    RCentroid = CentroidPlist(:,2);
                                    CCentroid = CentroidPlist(:,1);
                                    C13Sum = 0;
                                    for m = 1:size(C13R,1) 
                                       C13Sum =  double(imstruc.C13(C13R(m,1),C13C(m,1))) + C13Sum;
                                    end
                                    C13Avg = C13Sum./size(C13R,1);
                                    ObjSize = imObjs(n).Area;
                                    ObjSizeArr = [ObjSizeArr;ObjSize];
                                    C13AvgArr = [C13AvgArr ; C13Avg];
                                    CentroidArr = [CentroidArr; RCentroid, CCentroid];
                                end


                                % Save Contactsiteimage as a file, along with the summary table
                                if ~isempty(CentroidArr)
                                    %Creates summary table to save for the test

                                    %  normal one v
                                    %Sumtable = [[1:size(imObjs,1)]', round(CentroidArr(:,2)),round(CentroidArr(:,1)), ones(size(CentroidArr,1),1), perimarr, ObjSizeArr, contactsarr', numcontactsarr, C13AvgArr, percentcontactarr];
                                    % run this when numcontactsarr isnt ready
                                    Sumtable = [[1:size(imObjs,1)]', round(CentroidArr(:,2)),round(CentroidArr(:,1)), ones(size(CentroidArr,1),1), perimarr, ObjSizeArr, contactsarr', C13AvgArr, percentcontactarr];
                                    Sumtable = array2table(Sumtable);

                                    %  normal one v
                                    %Sumtable.Properties.VariableNames = ["Object", "X","Y", "Z", "Length_of_Perimeter_in_Pixels","Object_Size_in_Pixels", "Pixel_Contact_Count", "Contact_Site_Count","Avg_C13_per_Obj" , "Percent_Perimeter_Contact"];
                                    % run this when numcontactsarr isnt ready
                                    Sumtable.Properties.VariableNames = ["Object", "X","Y", "Z", "Length_of_Perimeter_in_Pixels","Object_Size_in_Pixels", "Pixel_Contact_Count","Avg_C13_per_Obj" , "Percent_Perimeter_Contact"];


                                    strcat(string(numel(find(contactsarr))),' Passed the Test ', maskorder{maskiter1}, ' to ',maskorder{maskiter2})







                                    % end main loop here





                                    fullpath = strcat(mousepath,filesep,mnames{mnum},filesep,hnames{hnum},filesep,cnames{cnum},filesep,'Final_Outputs', ...
                                    '\Summary_Table_Contact_Test_',maskorder{maskiter1}, '_to_', maskorder{maskiter2},'_', cnames{cnum}, '.csv');
                                    writetable(Sumtable,fullpath,'WriteMode','overwrite');

                                    if ContactMask == 1
                                        fullpath2 = strcat(mousepath,filesep,mnames{mnum},filesep,hnames{hnum},filesep,cnames{cnum},filesep,'Final_Outputs',...
                                        '\Contacts_',maskorder{maskiter1}, '_to_', maskorder{maskiter2}, '_', cnames{cnum}, '.tiff');
                                        imwrite(ContactsiteImage, fullpath2,'WriteMode','overwrite');
                                    end

                                end




                            else
                                %if heterotypic:

                                %Obtain image classification for original im
                                I = postimstruc.(maskorder{maskiter1});
                                SearchImage = postimstruc.(maskorder{maskiter2});


                                CCstruc = bwconncomp(I);
                                imObjs = regionprops(CCstruc,"Area", "Centroid", "PixelList"); %cant use perimeter, gives distance instead of pixels


                                %Perform actual test here
                                foundcell = DistContactTest(imObjs, SearchImage, threshold); %CHECK THIS CODE ***********************
                                perimarr = PerimPixCount(imObjs,size(I));


                                %Delete Single Pixel Contact Sites
                                foundcell = IsoPixelContactDelete(foundcell, size(I));



                                %Create contactsarr from old code, phase out eventually
                                contactsarr = [];
                                for s = 1:size(foundcell,2)
                                    if ~isempty(foundcell{s})
                                        contactsarr = [contactsarr, size(foundcell{s},1)]; 
                                    else
                                        contactsarr = [contactsarr, 0];
                                    end
                                end


                                %numcontactsarr = ContactSiteCount(foundcell, sepcontactlength);

                                % Generate Image Masks and Contact Statistics



                                %Array for percent of perimeter is in contact with search object
                                percentcontactarr=[];

                                %iterate over contactsitecell, ~ifempty, count the rows,
                                %add, reset every new column

                                for n = 1:size(contactsarr,2)

                                   if contactsarr(n) >0
                                      percentcontact = contactsarr(n)./ perimarr(n);
                                   else
                                       percentcontact = 0;
                                   end
                                   percentcontactarr = [percentcontactarr; percentcontact];
                                end

                                %Create image masks for contact sites, if ContactMask = true
                                idx = 0;
                                ContactsiteImage = logical(zeros(size(I,1),size(I,2)));
                                if ContactMask == 1;
                                    for p = 1:size(foundcell,1)
                                        for r = 1:size(foundcell,2)
                                           if ~isempty(foundcell{p,r})
                                                for s = 1:size(foundcell{p,r})
                                                    ContactsiteImage(foundcell{p,r}(s,1),foundcell{p,r}(s,2)) = 1;

                                                end
                                           end
                                        end                   
                                    end
                                end

                                ContactsiteImage = uint16(ContactsiteImage);


                                %Avg C13 and Object Size for Sumtable.
                                %Also assigns centroid of object, want in [x,y]
                                C13AvgArr = [];
                                ObjSizeArr = [];
                                CentroidArr = [];
                                for n = 1:size(imObjs,1)               
                                    C13Plist = imObjs(n).PixelList; %[x,y]=[c,r]
                                    C13R = C13Plist(:,2);
                                    C13C = C13Plist(:,1);
                                    %Find Centroid here, note, this is inverted relative to
                                    %preprocessing. Should be column then row (x,y)
                                    CentroidPlist = imObjs(n).Centroid; %[x,y]
                                    RCentroid = CentroidPlist(:,2);
                                    CCentroid = CentroidPlist(:,1);
                                    C13Sum = 0;
                                    for m = 1:size(C13R,1) 
                                       C13Sum =  double(imstruc.C13(C13R(m,1),C13C(m,1))) + C13Sum;
                                    end
                                    C13Avg = C13Sum./size(C13R,1);
                                    ObjSize = imObjs(n).Area;
                                    ObjSizeArr = [ObjSizeArr;ObjSize];
                                    C13AvgArr = [C13AvgArr ; C13Avg];
                                    CentroidArr = [CentroidArr; RCentroid, CCentroid];
                                end

                                % Save Contactsiteimage as a file, along with the summary table
                                if ~isempty(CentroidArr)
                                    %Creates summary table to save for the test

                                    %  normal one v
                                    %Sumtable = [[1:size(imObjs,1)]', round(CentroidArr(:,2)),round(CentroidArr(:,1)), ones(size(CentroidArr,1),1), perimarr, ObjSizeArr, contactsarr', numcontactsarr, C13AvgArr, percentcontactarr];
                                    % run this when numcontactsarr isnt ready
                                    Sumtable = [[1:size(imObjs,1)]', round(CentroidArr(:,2)),round(CentroidArr(:,1)), ones(size(CentroidArr,1),1), perimarr, ObjSizeArr, contactsarr', C13AvgArr, percentcontactarr];
                                    Sumtable = array2table(Sumtable);

                                    
                                    Sumtable.Properties.VariableNames = ["Object", "X","Y", "Z", "Length_of_Perimeter_in_Pixels","Object_Size_in_Pixels", "Pixel_Contact_Count","Avg_C13_per_Obj" , "Percent_Perimeter_Contact"];


                                    strcat(string(numel(find(contactsarr)))," Passed the Test ", strcat(maskorder{maskiter1}," to ",maskorder{maskiter2}))

                                    % end main loop here
                                    fullpath = strcat(mousepath,filesep,mnames{mnum},filesep,hnames{hnum},filesep,cnames{cnum},filesep,'Final_Outputs', ...
                                    '\Summary_Table_Contact_Test_',maskorder{maskiter1}, '_to_', maskorder{maskiter2},'_', cnames{cnum}, '.csv');
                                    writetable(Sumtable,fullpath,'WriteMode','overwrite');

                                    if ContactMask == 1
                                        fullpath2 = strcat(mousepath,filesep,mnames{mnum},filesep,hnames{hnum},filesep,cnames{cnum},filesep,'Final_Outputs',...
                                        '\Contacts_',maskorder{maskiter1}, '_to_', maskorder{maskiter2}, '_', cnames{cnum}, '.tiff');
                                        imwrite(ContactsiteImage, fullpath2,'WriteMode','overwrite');
                                    end
                                end


                            end
                        end
                    end
                end
%             catch
%                 error(['Error in File ' strcat('mnum',num2str(mnum),' hnum',num2str(hnum),' cnum',num2str(cnum))])
%                 
%             end
        end
    end
end


display(errorlog)









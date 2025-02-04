%The function folder needs to be added to the Matlab path; it contains many
%custom functions shared across the multiple scripts, or alternatively they
%can be added as local functions at the bottom of the scripts that call
%them.

%The structure setup folder contains example files for creating the nested
%file organizational structure that is used as the folder path input in 
%the scripts contained in the Contact_Testing, Demographics_Summary, and
%Centroid_Preprocessing scripts. Although these aren't necessary, as we
%manually created some of the nested file organizational structures as
%well, and the scripts we used to create these structures varied greatly
%depending on the output of the deep learning segmentations.

%The ideal nested file organizational structure depends varies between
%datasets, as well as the preprocessing steps, leading to multiple versions
%of the code depending on which dataset is being run in the scripts contained
%in the Contact_Testing, Demographics_Summary, and Centroid_Preprocessing
%folders, as well as the custom functions. The sublabels between
%the different datasets also vary as mentioned, requiring different regular
%expression patterns to extract the sublabels.

%For example, the scripts containing the suffix "AllMice" use the 3 sub
%labels "mouse", "hour", and "cell" in that order. Thus, the nested file
%organizational structure would be:
% 1. A folder with the dataset title, containing:
% 2. The folders for each mouse, in our case, "mouse_1", "mouse_2", and
% "mouse_3". These would contain:
% 3. The folders for each hour timepoint, in our case "hour_1" and
% "hour_4", which would contain:
% 4. The folders for each cell at that timepoint, labeled "cell_x", where x
% is a whole, nonzero number. The cell folders would contain:
% 5. The confidence interval thresholded organelle masks, obtained from our
% deep learning segmentation pipelines. Each masks is binary, where the
% value is True where the organelle is expected to be find according to the
% selected confidence interval. The isotope mask and EM image (if available)
% are both Uint16. Finally, the cell folder (or whatever the final sublabel
% is depending on the dataset) would then contain a folder titled:
% 6. "Final_Outputs", which would eventually be propogated with the outputs
% from the scripts.

%This is just an example from our Spatial_Liver mice dataset. Other
%datasets can have more or less than 3 sublabels, and they may have
%different names, such as "ROI", "Area", "Diet", etc. . .

%For these scripts, the path listed in the first few lines of the code
%should be the absolute file path of the folder titled after the dataset as
%described above, containing the first sublabel folders. Furthermore,
%the available organelle masks vary between datasets, and cause differences
%in the codes between versions. For simplicity sake, in all the scripts in
%the "Contact_Testing" folder, the variable "savemaskopt" should be set to
%false, and you may disregard the "savemaskpath" variable. All other
%variables defined in that script section should be left untouched to
%repliccate our search criteria.

%The order we ran the scripts for each respective dataset would be as
%follows:
% 1. The script contained in the "Contact_Testing" folder for a given
% dataset suffix
% 2. The script contained in the "Centroid_Preprocessing" folder for the
% same suffix
% 3. The script contained in the "Demographics_Summaries" folder for the
% same suffix

% although steps 2 and 3 are interchangeable.



%Finally, included in the Concatenations_and_Sorting folder is a number of
%scripts we used to concatenate the outputs and consolidate/alter files 
%within the nested file organizational structure, although these are not
%necessary and were mainly used for readability.


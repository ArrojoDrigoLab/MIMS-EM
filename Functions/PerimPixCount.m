function [pcarr] = PerimPixCount(ObjStruc, sizearr)

%This function takes the structure output from regionprops and returns an
%array with scalar entries representing the # of perimeter pixels for the
%respective image object

%Inputs:

% ObjStruc: Structure array [n x 1], with 1 field at the minimum,
%PixelList. Obtained from running bwconncomp and regionprops on the
%original image.

% sizearr: [1x2] array with the size dimensions of the original image

%Outputs: 

%pcarr: array, [n x 1] each entry is # of perimeter pixels per image obj
kern8p = [-1 -1 -1;
        -1 8 -1;
        -1 -1 -1];

pcarr = [];




for p = 1:size(ObjStruc,1)
    %Create list of pix
    tempim = logical(zeros(sizearr(1),sizearr(2)));
    tplist = ObjStruc(p).PixelList; % [x,y]
    tplist = [tplist(:,2),tplist(:,1)]; %[r,c]
    
    %skeletonize, generate list again
    for r = 1:size(tplist,1)
        tempim(tplist(r,1),tplist(r,2)) = 1;
    end
    tempskel = imbinarize(conv2(tempim, kern8p, 'same'));
    perimc = numel(find(tempskel));
    
    
    pcarr = [pcarr;perimc];
    
end

end


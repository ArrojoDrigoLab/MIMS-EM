function [dfcell] = DistContactTest(ObjStruc,SearchI,thresholdarr)
% Inputs: DistContactTest(ObjStruc, SearchI,thresholdarr);

%ObjStruc: Structure array [n x 1], with 1 field at the  minimum,
%PixelList. Obtained from running bwconncomp and regionprops on the
%original image.

%SearchI: The image to be searched on using the original image; already
%should be preprocessed at this point.

%thresholdarr: an array [1 x 2] containing the dilation search area, or
%dilation iterations



% Output: dfcell, a cell array containing the pixels on the perimeter of the
% search objects that passed the test. Each cell corresponds to a different
% image object


kern8p = [-1 -1 -1;
        -1 8 -1;
        -1 -1 -1];

dfcell = cell(1,size(ObjStruc,1));

boxor = [thresholdarr(2)+1,thresholdarr(2)+1];

for p = 1:size(ObjStruc,1)
    % Create list of pix
    tempim = logical(zeros(size(SearchI)));
    tplist = ObjStruc(p).PixelList; % [x,y]
    tplist = [tplist(:,2),tplist(:,1)]; %[r,c]
    
    %skeletonize, generate list again
    for r = 1:size(tplist,1)
        tempim(tplist(r,1),tplist(r,2)) = 1;
    end
    tempskel = imbinarize(conv2(tempim, kern8p, 'same'));
    [tsr,tsc] = find(tempskel);
    tplist = [tsr,tsc];

    

    for r = 1:size(tplist,1)
        %Create box centered around each pixel
        
        sbox = SearchI(tplist(r,1)-thresholdarr(2):tplist(r,1)+thresholdarr(2),tplist(r,2)-thresholdarr(2):tplist(r,2)+thresholdarr(2));
        [sr,sc] = find(sbox);
        if ~isempty(sr)
            for s = 1:size(sr,1)
                if norm([boxor(1)-sr(s),boxor(2)-sc(s)]) >= thresholdarr(1) && norm([boxor(1)-sr(s),boxor(2)-sc(s)]) <= thresholdarr(2)
                    
                    %store
                    dfcell{p} = [dfcell{p}; tplist(r,[1 2]) ];
                    break
                end
    
    
            end
        end

    end




end
    






end


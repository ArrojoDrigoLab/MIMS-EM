function [outcell] = SquidPreprocStandard(immask,c13mask, imborderthresh)
%SquidPreProcStandard creates a cellarray 1 x n where each of the n cells is the pixel coords of each unique image object 
%If all objects are removed, returns an empty cell arr 

%inputs:
%immask,image mask to be processed
%c13 mask, c13 readout
%imborderthresh, objects containing points within this threshold (distance)
%to the boundary of the image will be removed

%output format: [1 x n] cell array 





immask(immask > 0) = 1;
imstruc = bwconncomp(immask); %pixels in [r,c]
imobjs = regionprops(imstruc, 'PixelList'); %pixels in [x,y] = [c,r]


% remove objects on the edge of the image 

deladj = 0;
for p = 1:size(imobjs,1)
    pixlist = imobjs(p-deladj).PixelList;
    %Dimension 1, (x),(c)
    borpixl = pixlist(:,1) > imborderthresh & pixlist(:,1) < (imstruc.ImageSize(2) - imborderthresh);
    if ~isempty(find(~borpixl,1))
        imobjs(p-deladj) = [];
        deladj = deladj + 1;
    else
        %Dimension 2 (y),(r), only if Dim 1 is within bounds
        borpixl = pixlist(:,2) > imborderthresh & pixlist(:,2) < (imstruc.ImageSize(1) - imborderthresh);
        if ~isempty(find(~borpixl,1))
           imobjs(p-deladj) = [];
           deladj = deladj + 1;
        end

    end

end

outcell = {};

try
    for p = 1:size(imobjs,1)
       outcell{1,p} = imobjs(p).PixelList; % [x,y] = [c,r] = dim2,dim1 
       c13stow = [];
       for r = 1:size(outcell{1,p},1)
          c13stow(r,1) = c13mask(outcell{1,p}(r,2),outcell{1,p}(r,1)); 
       end
       outcell{1,p} = [outcell{1,p},c13stow];
    end
catch
    outcell = {};
end





end


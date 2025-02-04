function [carrout] = CellArrayCatUneq(carr1,carr2,dim)
%CellArrayCatUneq 
%   Compares the size of two unequal cell arrays, modifies, and then concatenates
%   them together

%Inputs:
        % carr1, cell array 1, first to be concatenated

        % carr 2, cell array 2, second in order. If inputs are
        % (carr1,carr2,~), order of cat is [carr1;carr2](dim == 1) or
        % [carr1,carr2] (dim == 2)
        % depending on dimensions
        
        % dim, dimension of concatenation, (either 1 or 2) refer to carr 2 description
        
%Outputs:

        %carrout, [n x m] cell array, where n and m are determined by the
        %cell array inputs 
        

        
if dim == 1
    carr1size = size(carr1);
    carr2size = size(carr2);
    
    if carr1size(2) == carr2size(2)
        carrout = [carr1;carr2]; 
    end
    
    
    if carr1size(2) > carr2size(2)
        for p = 1:carr2size(1)
            for r = carr2size(2)+1 : carr1size(2)
                carr2{p,r} = {};
            end
        end
        
        carrout = [carr1;carr2];
    end
    
    if carr1size(2) < carr2size(2)
        for p = 1:carr1size(1)
           for r = carr1size(2) + 1 : carr2size(2)
                carr1{p,r} = {}; 
           end
        end
        
        carrout = [carr1;carr2];
    end


end




if dim == 2
    carr1size = size(carr1);
    carr2size = size(carr2);
    
    if carr1size(1) == carr2size(1)
        carrout = [carr1,carr2]; 
    end
    
    if carr1size(1) > carr2size(1)
        for p = 1:carr2size(2)
            for r = carr2size(1) + 1 : carr1size(1)
               carr2{r,p} = {}; 
            end
        end
        
        carrout = [carr1,carr2];
        
    end
    
    
    if carr1size(1) < carr2size(1)
        for p = 1:carr1size(2)
            for r = carr1size(1) + 1 : carr2size(1)
               carr1{r,p} = {}; 
            end
        end
        
        carrout = [carr1,carr2];
    end
    
    
end
        



        
if dim ~= 1 && dim ~= 2
    error("Dimension must have a value of 1 or 2") 
    
end


end


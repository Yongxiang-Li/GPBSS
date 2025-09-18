function [obj,x] = Fun_Hills(x,addNoise)
%--------------------------------------------------------
% The Hills problem
% range[0,100]x[0,100], global = [90,90]
% fmax = 20
%--------------------------------------------------------
    if nargin < 1
        x = nan;    obj = nan;
        return; 
    end
    if nargin < 2, addNoise = false; end

    design_space=[0,0;100,100];    x = design_space(1,:) + diff(design_space).*x;
        
    obj = 10*(sin(0.05*pi.*x(:,1))).^6./(2.^(2*((x(:,1)-90)/80).^2))+...
    10*(sin(0.05*pi.*x(:,2))).^6./(2.^(2*((x(:,2)-90)/80).^2));

    if addNoise 
        obj = obj+normrnd(0,0.5,[size(x,1),1]);
    end
end














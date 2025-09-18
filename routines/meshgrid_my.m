function [ X ] = meshgrid_my(X1, X2)
    [X1, X2] = meshgrid(X1, X2);
    X = [X1(:) X2(:)];
end


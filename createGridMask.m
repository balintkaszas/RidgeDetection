function mask = createGridMask(grid, shape, arbitraryPoints)
%%Creates a sparse matrix, given a uniform 2d grid and the coordinates of
%%given points. The resulting matrix is a binary mask, that has an entry 1
%%in the location, where a point falls close to a gridpoint.
%grid is N by 2, arbitraryPoints is m by 2.
% shape is the shape of the grid, N = prod(shape);
k = dsearchn(grid, arbitraryPoints); %returns linear indices, at which grid(k,:) ~= arbitraryPoints(j,:) for some j.
zeroMatrix = zeros(size(grid));
zeroMatrix(k,:) = 1;
mask = reshape(zeroMatrix(:,1), shape);
end 


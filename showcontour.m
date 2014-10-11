function showcontour(src, fignum)
%
%src = 'bk'

name_fits = strcat(src, '.fits');
img = fitsread(name_fits);

[ny, nx] = size(img);

[maxflux, ry, rx] = arr_max(img);

figure(fignum);
v = [-0.01 0.01 0.02 0.04 0.08 0.16 0.32 0.64] * maxflux;
[C, h] = contour(img, v);
%clabel(C, h, 'FontSize', 12);
colorbar();
axis square;

end

function [maxval, row, col] = arr_max(arr)
    [maxval, maxloc] = max(arr(:));
    [row, col] = ind2sub(size(arr), maxloc);
end

function [minval, row, col] = arr_min(arr)
    [minval, minloc] = min(arr(:));
    [row, col] = ind2sub(size(arr), minloc);
end
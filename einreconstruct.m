function einreconstruct
%
%bin_name = 'ein_0.05_10000.bin';
%bin_name = 'ein_0.01_50000.bin';
%bin_name = 'ein_0.002_100000.bin';
%bin_name = 'ein_0.002_200000.bin';
bin_name = 'ein_0.001_200000.bin';

fid = fopen(bin_name, 'r');
ng = fread(fid, [1,1], 'int')
niter = fread(fid, [1,1], 'int')
gain = fread(fid, [1,1], 'double')
flux = fread(fid, [1, niter], 'double');
arx = fread(fid, [1, niter], 'int16');
ary = fread(fid, [1, niter], 'int16');
res = fread(fid, [ng, ng], 'double');
fclose(fid);

ng4 = ng / 4;

img = zeros(ng, ng);

nb = 6;
gaux = 2;
gauy = 2;
pixw = 1;
patch = zeros(2*nb+1, 2*nb+1);

nb1 = nb + 1;
gaux2 = gaux^2;
gauy2 = gauy^2;

for j = -nb:nb
    for i = -nb:nb
        rx2 = (i * pixw)^2;
        ry2 = (j * pixw)^2;
        tmpx = exp(-0.5 * rx2 / gaux2);
        tmpy = exp(-0.5 * ry2 / gauy2);
        patch(nb1 + j, nb1 + i) = tmpx * tmpy;
    end
end

for i = 1:length(flux)    
    for gj = -nb:nb
        for gi = -nb:nb
            imgi = gi + arx(i);
            imgj = gj + ary(i);
            if imgi > 1 && imgi <= ng && imgj > 1 && imgj <= ng
                img(imgj, imgi) = img(imgj, imgi) + flux(i) * patch(gj + nb1, gi + nb1);
            end
        end
    end
end

img = img + res;

figure(101);
imagesc(img(ng4+1:ng4*3, ng4+1:ng4*3));
axis image;
colormap(gray);
colorbar();

[maxflux, ry, rx] = arr_max(img);

figure(11);
v = [-0.01 0.01 0.02 0.04 0.08 0.16 0.32 0.64] * maxflux;
[C, h] = contour(img(ng4+1:ng4*3, ng4+1:ng4*3), v);
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

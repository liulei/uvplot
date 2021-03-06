function einsdi()
% plot (l,m) figure
%

ng = 512;

% c = 3E8;
% %freq = 1.5352E10;
% freq = 1.0;
% res_mas = 0.1; % in mas
% res_rad = res_mas * 1E-3 / 3600. / 180. * pi;

umax = 2.0 * pi;
vmax = umax;

ngh = ng / 2;
uinc = umax / ngh;
vinc = uinc;
ulimit = uinc * ng / 4;
vlimit = ulimit;

target = 'ein_center';
uvname = strcat(target, '.uv');
pngname = strcat(target, '_dirt_', num2str(ng), '.png');

offset = 3;

arr = importdata(uvname);

u = arr(:, 1);
v = arr(:, 2);

vis = complex(arr(:, offset + 1), arr(:, offset + 2));
weight = arr(:, offset + 3);

visarr_r = zeros(ng, ng);
visarr_c = zeros(ng, ng);
visarr = complex(visarr_r, visarr_c);

mask = zeros(ng, ng);

% beamarr_r = zeros(ng, ng);
% beamarr_c = zeros(ng, ng);
% beamarr = complex(beamarr_r, beamarr_c);

%idu = floor((u - uleft) / du) + 1;
%idv = floor((v - vleft) / dv) + 1;

for i=1:length(u)    
    idu = floor(u(i) / uinc + 0.5);
    if(idu < 0)
        idu = idu + ng;
    end
    idu = idu + 1;
    
    idv = floor(v(i) / vinc + 0.5);
    if(idv < 0)
        idv = idv + ng;
    end
    idv = idv + 1;
    
    visarr(idv, idu) = visarr(idv, idu) + vis(i);
    mask(idv, idu) = 1.0;
end

% for i = 1:length(u)
%     
%     idu = floor(u(i) / du)
%     wu = (u(i) - idu * du) / du;
%     idv = floor(v(i) / )
% 
% end

%figure(2);
%imagesc(abs(visarr));
%colormap(jet);
%colorbar();

%figure(3);
%imagesc(abs(beamarr));
%colormap(jet);
%colorbar();

dirt_img = ifft2(visarr);
dirt_beam = ifft2(mask);

% dirt_img = fftshift(dirt_img);
% dirt_beam = fftshift(dirt_beam);
% dirt_img = flipud(dirt_img);
% dirt_beam = flipud(dirt_beam);

dirt = real(dirt_img);
[minval, row, col] = arr_min(dirt);
dirt = dirt - minval;

ng4 = ng / 4;
figure(4);
%figure('NumberTitle', 'off', 'Name', 'Dirty image');
%imagesc(real(flipud(dirt_img)));
imagesc(arr_cen(dirt, ng4));
axis image;
colormap(gray);
colorbar();
%print(gcf, '-dpng', pngname);

figure(5);
%figure('NumberTitle', 'off', 'Name', 'Dirty beam');
imagesc(real(arr_cen(dirt_beam, ng4)));
%imagesc(real(flipud(dirt_beam)));
axis image;
colormap(gray);
colorbar();

res  = dirt;
beam = real(dirt_beam);
den = zeros(ng, ng);

gain = 0.005;
trim_contour = 0.9;
niter = 1000;

cln = zeros(ng, ng);

for k = 1:niter
    
    resc = arr_cen(res, ng4);
    
    [peak, row, col] = arr_max(resc);
    thresh = trim_contour * peak;
    
    for j = 1:ng
        for i = 1:ng
            if res(j, i) > thresh
                cln(j, i) = res(j, i);
            else
                cln(j, i) = 0.0;
            end
        end
    end
    
    peak_gain = peak * gain;
    cln = cln * gain;
    den = den + cln;
    
%     den_fft = fft2(den) .* mask;
%     product = real(ifft2(den_fft));
%     res = dirt - product;
    
    cln_fft = fft2(cln) .* mask;
    product = real(ifft2(cln_fft));
    [maxval, row, col] = arr_max(arr_cen(product, ng4));
    
    product = product / maxval * peak_gain;
    res = res - product;
    
    [maxval, row, col] = arr_max(arr_cen(res, ng4));
    [minval, row, col] = arr_min(arr_cen(res, ng4));
    if (maxval - minval) / maxval < 0.01
        break;
    end
    
    if mod(k, 20) == 0
        
        fprintf('Iteration %d, peak %f\n', k, peak);
        
        figure(6);
        imagesc(arr_cen(res, ng4));
        axis image;
        colormap(gray);
        colorbar();
        
        figure(7);
        imagesc(arr_cen(den, ng4));
        axis image;
        colormap(gray);
        colorbar();
        
        figure(8);
        imagesc(arr_cen(product, ng4));
        axis image;
        colormap(gray);
        colorbar();        
        
        figure(9);
        imagesc(arr_cen(cln, ng4));
        axis image;
        colormap(gray);
        colorbar();
        
    end 
end

end

function [cen] = arr_cen(arr, ng4)
    tmp = flipud(fftshift(arr));
    cen = tmp(ng4+1:ng4*3, ng4+1:ng4*3);
end

function [maxval, row, col] = arr_max(arr)
    [maxval, maxloc] = max(arr(:));
    [row, col] = ind2sub(size(arr), maxloc);
end

function [minval, row, col] = arr_min(arr)
    [minval, minloc] = min(arr(:));
    [row, col] = ind2sub(size(arr), minloc);
end


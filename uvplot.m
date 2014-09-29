function uvplot()
% plot (l,m) figure
%

ng = 512;

target = 'bk';
uvname = strcat(target, '.uv');
pngname = strcat(target, '_dirt_', num2str(ng), '.png');

offset = 3;

arr = importdata(uvname);

u = arr(:, 1);
v = arr(:, 2);
vis = complex(arr(:, offset + 1), arr(:, offset + 2));
weight = arr(:, offset + 3);

mu = max(u);
mv = max(v);
maxuv = max(mu, mv) * 1.00001;

maxuv = maxuv * 4.6;

fsize = 17;
figure(1);
h = gca;
set(h, 'FontSize', fsize);
set(findall(h, 'type', 'text'), 'FontSize', fsize);
plot(u, v, 'ko', 'MarkerEdgeColor', 'None', 'MarkerFaceColor', 'k', 'MarkerSize', 2);
xlim([-maxuv, maxuv]);
ylim([-maxuv, maxuv]);
axis square;
xlabel('u');
ylabel('v');


uleft = -maxuv ;
vleft = -maxuv ;
uright = -uleft;
vright = -vleft;
du = (uright - uleft) / ng;
dv = (vright - vleft) / ng;

visarr_r = zeros(ng, ng);
visarr_c = zeros(ng, ng);
visarr = complex(visarr_r, visarr_c);

beamarr_r = zeros(ng, ng);
beamarr_c = zeros(ng, ng);
beamarr = complex(beamarr_r, beamarr_c);

%idu = floor((u - uleft) / du) + 1;
%idv = floor((v - vleft) / dv) + 1;

for i=1:length(u)
    
    idu = floor(u(i) / du);
    if(idu < 0)
        idu = idu + ng;
    end
    idu = idu + 1;
    
    idv = floor(v(i) / dv);
    if(idv < 0)
        idv = idv + ng;
    end
    idv = idv + 1;
    
    visarr(idv, idu) = visarr(idv, idu) + vis(i);
    beamarr(idv, idu) = 1.0;
end

%figure(2);
%imagesc(abs(visarr));
%colormap(jet);
%colorbar();

%figure(3);
%imagesc(abs(beamarr));
%colormap(jet);
%colorbar();

dirt_img = ifft2(visarr);
dirt_beam = ifft2(beamarr);

dirt_img = fftshift(dirt_img);
dirt_beam = fftshift(dirt_beam);

figure(4);
%imagesc(real(flipud(dirt_img)));
imagesc(real(dirt_img));
axis image;
colormap(gray);
colorbar();
print(gcf, '-dpng', pngname);

figure(5);
imagesc(real(dirt_beam));
%imagesc(real(flipud(dirt_beam)));
axis image;
colormap(gray);
colorbar();

res = real(dirt_img);
beam = real(dirt_beam);
[bmax, by, bx] = arr_max(beam);
beam = beam / bmax;

gain = 0.001;

niter = 1700;
flux = zeros(niter);
ary = zeros(niter, 'int32');
arx = zeros(niter, 'int32');

[maxflux, ry, rx] = arr_max(res);
[minflux, ry, rx] = arr_min(res);
fprintf('Before clean: %.4f --- %.4f\n', minflux, maxflux);

sum = 0.0;
for i = 1:niter
    [maxflux, ry, rx] = arr_max(res);
    %[minval, row, col] = arr_min(dirt_img)
    
    rmv = gain * maxflux * beam;
    bry = by - ry;
    brx = bx - rx;
    
    for bj = 1:ng
        for bi = 1:ng
             rj = bj + bry;
             ri = bi + brx;
             if rj > 0 & rj <= ng & ri > 0 & ri <= ng
                 res(rj, ri) = res(rj, ri) - rmv(bj, bi);
             end %if
        end % bi
    end % bj
    
    flux(i) = maxflux * gain;
    ary(i) = ry;
    arx(i) = rx;
    
    sum = sum + flux(i);
    if mod(i, 50) == 0
        [maxflux, ry, rx] = arr_max(res);
        [minflux, ry, rx] = arr_min(res);
        fprintf('Iteration %5d, flux cleaned: %.4f, %.4f --- %.4f\n', i, sum, minflux, maxflux);
        figure(6);
        imagesc(res);
        axis image;
        colormap(gray);
        colorbar();
    end
end
    
bw = 0.46;
pw = 0.1;

img = zeros(ng, ng);
irmax = 2 * ceil(bw / pix);

gauss = zeros(irmax * 2 + 1, irmax * 2 + 1);
for j = -irmax:irmax
    for i = -irmax:irmax
        rx = i * pw;
        ry = j * pw;
        tmp = (rx / bw)^2 + (ry / bw)^2;
        gauss(j + irmax + 1, i + irmax + 1) = exp(-())

for i = 1:length(flux)
    
    y0 = ary(i) - irmax;
    if y0 < 1
        y0 = 1;
    end
    y1 = ary(i) + irmax;
    if y1 > ng
        y1 = ng;
    end
    
    x0 = arx(i) - irmax;
    if x0 < 1
        x0 = 1;
    end
    x1 = arx(i) + irmax;
    if x1 > ng
        x1 = ng;
    end
    
    for j = y0:y1
        for i = x0:x1
            

end


function [maxval, row, col] = arr_max(arr)
    [maxval, maxloc] = max(arr(:));
    [row, col] = ind2sub(size(arr), maxloc);
end

function [minval, row, col] = arr_min(arr)
    [minval, minloc] = min(arr(:));
    [row, col] = ind2sub(size(arr), minloc);
end


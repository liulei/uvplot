function uvplot()
% plot (l,m) figure
%

ng = 512;

c = 3E8;
%freq = 1.5352E10;
freq = 1.0;
res_mas = 0.1; % in mas

res_rad = res_mas * 1E-3 / 3600. / 180. * pi;
uinc = 1. / (res_rad * ng);
vinc = uinc;
ulimit = uinc * ng / 4;
vlimit = ulimit;
fprintf('Required uv: %f, res: %f\n', ulimit, res_mas);

target = 'bk';
uvname = strcat(target, '.uv');
pngname = strcat(target, '_dirt_', num2str(ng), '.png');

offset = 3 + 9;

arr = importdata(uvname);

u = arr(:, 1);
v = arr(:, 2);

vis = complex(arr(:, offset + 1), arr(:, offset + 2));
weight = arr(:, offset + 3);

for offset = 6:3:12
    vis = vis + complex(arr(:, offset + 1), arr(:, offset + 2));
    weitht = weight + arr(:, offset + 3);
end

u = u * freq;
v = v * freq;

maxuv = max(u.^2 + v.^2);
maxuv = sqrt(maxuv);
minres = 1.0 / maxuv * 180. / pi * 3600. * 1000.;
fprintf('Provided max uv: %f, min res: %f\n', maxuv, minres);

%maxuv = maxuv * 4.2;

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


%uleft = -maxuv ;
%vleft = -maxuv ;
uleft = -ulimit * 2.0;
vleft = -vlimit * 2.0;
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
dirt_beam = ifft2(beamarr);

dirt_img = fftshift(dirt_img);
dirt_beam = fftshift(dirt_beam);

dirt_img = flipud(dirt_img);
dirt_beam = flipud(dirt_beam);

ng4 = ng / 4;
figure(4);
%imagesc(real(flipud(dirt_img)));
imagesc(real(dirt_img(ng4 + 1: ng4 * 3, ng4 + 1: ng4 * 3)));
axis image;
colormap(gray);
colorbar();
print(gcf, '-dpng', pngname);

figure(5);
imagesc(real(dirt_beam(ng4 + 1: ng4 * 3, ng4 + 1: ng4 * 3)));
%imagesc(real(flipud(dirt_beam)));
axis image;
colormap(gray);
colorbar();

res = real(dirt_img);
beam = real(dirt_beam);
[bmax, by, bx] = arr_max(beam)
beam = beam / bmax;

nb = 15;
xarr = bx - nb: bx + nb;
yarr = by - nb: by + nb;
beam(by, xarr)
figure(100);
plot(xarr, beam(by, xarr), 'r-');
hold on;
plot(yarr, beam(yarr, bx), 'b-');


%beam(by - 5: by + 5, bx - 5: bx + 5)

gain = 0.001;

niter = 20000;
flux = zeros(1, niter);
ary = zeros(1, niter, 'int16');
arx = zeros(1, niter, 'int16');

[maxflux, ry, rx] = arr_max(res);
[minflux, ry, rx] = arr_min(res);
fprintf('Before clean: %.4f --- %.4f\n', minflux, maxflux);

sum = 0.0;
for i = 1:niter
    [maxflux, ry, rx] = arr_max(res);
    %[minval, row, col] = arr_min(dirt_img)
    
    rmv = gain * maxflux * beam;
    rby = ry - by;
    rbx = rx - bx;
    
    for bj = 1:ng
        for bi = 1:ng
             rj = bj + rby;
             ri = bi + rbx;
             if rj > 0 && rj <= ng && ri > 0 && ri <= ng
                 res(rj, ri) = res(rj, ri) - rmv(bj, bi);
             end %if
        end % bi
    end % bj
    
    flux(i) = maxflux * gain;
    ary(i) = ry;
    arx(i) = rx;
    
    sum = sum + flux(i);
    if mod(i, 100) == 0
        [maxflux, ry, rx] = arr_max(res);
        [minflux, ry, rx] = arr_min(res);
        fprintf('Iteration %5d, flux cleaned: %.4f, %.4f --- %.4f\n', i, sum, minflux, maxflux);
    end
end

figure(6);
imagesc(res(ng4+1:ng4*3, ng4+1:ng4*3));
axis image;
colormap(gray);
colorbar();

img = zeros(ng, ng);

% nb = 8;
% %patch = zeros(2*nb+1, 2*nb+1);
% patch = beam(by-nb:by+nb, bx-nb:bx+nb);
% nb1 = nb+1;
% for j = -nb: nb
%     for i = -nb: nb
%         if patch(nb1+j, nb1+i) < 0.0
%             patch(nb1+j, nb1+i) = 0.0;
%         end
%     end
% end

nb = 25;
gaux = 0.4;
gauy = 1.0;
pixw = 0.1;
patch = zeros(2*nb+1, 2*nb+1);

nb1 = nb + 1;
gaux2 = gaux^2;
gauy2 = gauy^2;

for j = -nb:nb
    for i = -nb:nb
        rx2 = (i * pixw)^2;
        ry2 = (j * pixw)^2;
        tmpx = exp(-1.99 * rx2 / gaux2);
        tmpy = exp(-1.99 * ry2 / gauy2);
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

figure(7);
imagesc(img(ng4+1:ng4*3, ng4+1:ng4*3));
axis image;
colormap(gray);
colorbar();

[maxflux, ry, rx] = arr_max(img);

figure(8);
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


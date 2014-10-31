function uvcount()
% plot (l,m) figure
%

ng = 512;

c = 3E8;
%freq = 1.5352E10;
freq = 1.0;
res_mas = 0.1; % in mas

res_rad = res_mas * 1E-3 / 3600. / 180. * pi;
umax = 1.0 / res_rad / 2.0;
%umax = pi;

uinc = umax / ng * 2;
vinc = uinc;
ulimit = uinc * ng / 4;
vlimit = ulimit;
fprintf('Required uv: %f, res: %f\n', ulimit, res_mas);

src = 'bk';
%src = '16B';
%src = 'ein';
uvname = strcat(src, '.uv');
%pngname = strcat(target, '_dirt_', num2str(ng), '.png');

offset = 3;

arr = importdata(uvname);

u = arr(:, 1);
v = arr(:, 2);

vis = complex(arr(:, offset + 1), arr(:, offset + 2));
weight = arr(:, offset + 3);

% for offset = 6:3:12
%     vis = vis + complex(arr(:, offset + 1), arr(:, offset + 2));
%     weitht = weight + arr(:, offset + 3);
% end

u = u * freq;
v = v * freq;

maxuv = max(u.^2 + v.^2);
maxuv = sqrt(maxuv);
minres = 1.0 / maxuv * 180. / pi * 3600. * 1000.;
fprintf('Provided max uv: %f, min res: %f\n', maxuv, minres);

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

nmeas = length(u);
for i=1:nmeas   
    idu = floor(u(i) / du + 0.5);
    if(idu < 0)
        idu = idu + ng;
    end
    idu = idu + 1;
    
    idv = floor(v(i) / dv + 0.5);
    if(idv < 0)
        idv = idv + ng;
    end
    idv = idv + 1;
    
    visarr(idv, idu) = visarr(idv, idu) + vis(i);
    beamarr(idv, idu) = beamarr(idv, idu) + 1.0;
end

figure(2);
imagesc(abs(fftshift(beamarr)))
colormap(jet);
colorbar();

dirt_img = ifft2(visarr) * ng * ng / nmeas;
dirt_img = fftshift(dirt_img);
dirt_img = flipud(dirt_img);

dirt_beam = ifft2(beamarr);
dirt_beam = fftshift(dirt_beam);
dirt_beam = flipud(dirt_beam);

ng4 = ng / 4;
figure(4);
%imagesc(real(flipud(dirt_img)));
imagesc(real(dirt_img(ng4 + 1: ng4 * 3, ng4 + 1: ng4 * 3)));
axis image;
colormap(jet);
colorbar();

figure(5);
imagesc(real(dirt_beam(ng4 + 1: ng4 * 3, ng4 + 1: ng4 * 3)));
%imagesc(real(flipud(dirt_beam)));
axis image;
colormap(jet);
colorbar();


res = real(dirt_img);
beam = real(dirt_beam);
[imgmin, ~, ~] = arr_min(res);
[imgmax, ~, ~] = arr_max(res);
[bmax, by, bx] = arr_max(beam);
beam = beam / bmax;

fprintf('imgmax: %f, imgmin: %f, imgmax/imgmin: %f\n', ...
    imgmax, imgmin, imgmax / imgmin);

return;
end

function [maxval, row, col] = arr_max(arr)
    [maxval, maxloc] = max(arr(:));
    [row, col] = ind2sub(size(arr), maxloc);
end

function [minval, row, col] = arr_min(arr)
    [minval, minloc] = min(arr(:));
    [row, col] = ind2sub(size(arr), minloc);
end

function uvplot()
% plot (l,m) figure
%

offset = 3;

arr = importdata('out.dat');
u = vertcat(arr(:, 1), -arr(:, 1));
v = vertcat(arr(:, 2), -arr(:, 2));

tmp = complex(arr(:, offset + 1), arr(:, offset + 2));
vis = vertcat(tmp, conj(tmp));
weight = vertcat(arr(:, offset + 3), arr(:, offset + 3));

mu = max(u);
mv = max(v);
maxuv = max(mu, mv) * 1.05;

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

ng = 256;
uleft = -maxuv * 2;
vleft = -maxuv * 2;
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

idu = floor((u - uleft) / du) + 1;
idv = floor((v - vleft) / dv) + 1;

for i=1:length(idu)
    visarr(idv(i), idu(i)) = visarr(idv(i), idu(i)) + vis(i);
    beamarr(idv(i), idu(i)) = 1.0;
end

figure(2);
imagesc(abs(visarr));
colormap(jet);
colorbar();

figure(3);
imagesc(abs(beamarr));
colormap(jet);
colorbar();

dirt_img = fftshift(fft2(visarr));
beam_img = fftshift(fft2(beamarr));

figure(4);
imagesc(abs(dirt_img));
colormap(gray);
colorbar();
print(gcf, '-dpng', 'dirt.png');

figure(5);
imagesc(abs(beam_img));
colormap(gray);
colorbar();
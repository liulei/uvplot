function uvplot()
% plot (l,m) figure
%

ng = 256;

target = 'bk';
uvname = strcat(target, '.uv');
pngname = strcat(target, '_dirt_', num2str(ng), '.png');

offset = 3;

arr = importdata(uvname);

u = arr(:, 1);
v = arr(:, 2);
vis = complex(arr(:, offset + 1), arr(:, offset + 2));
weight = arr(:, offset + 3);

%u = vertcat(arr(:, 1), -arr(:, 1));
%v = vertcat(arr(:, 2), -arr(:, 2));
%tmp = complex(arr(:, offset + 1), arr(:, offset + 2));
%vis = vertcat(tmp, conj(tmp));
%weight = vertcat(arr(:, offset + 3), arr(:, offset + 3));

mu = max(u);
mv = max(v);
maxuv = max(mu, mv) * 1.00001;

maxuv = maxuv * 2;

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
beam_img = ifft2(beamarr);

dirt_img = fftshift(dirt_img);
beam_img = fftshift(beam_img);

figure(4);
imagesc(real(flipud(dirt_img)));
%imagesc(real(dirt_img));
axis image;
colormap(gray);
colorbar();
print(gcf, '-dpng', pngname);

%figure(5);
%imagesc(abs(beam_img));
%colormap(gray);
%colorbar();

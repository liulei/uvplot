function prepareuv
%
ng = 512;
source = 'Einstein';
fits_name = strcat(source, '.fits');
data = fitsread('Einstein.fits');

img = zeros(ng, ng);
ng4 = ng / 4;
img(ng4+1:ng4*3, ng4+1:ng4*3) = data;

figure(1);
imagesc(flipud(img(ng4+1:ng4*3, ng4+1:ng4*3)));
axis image;
colormap(gray);
colorbar();

img = fftshift(img);
uvarr = fft2(img) / ng;

arr = importdata('cont_sim4.vis');
u = arr(:, 2);
v = arr(:, 3);

ulimit = pi;
vlimit = pi;
uinc = 4.0 * pi / ng;
vinc = uinc;
len = floor(length(u) * 0.2);
vis = complex(zeros(1,len), zeros(1, len));
visarr = complex(zeros(ng, ng), zeros(ng, ng));
beamarr = complex(zeros(ng, ng), zeros(ng, ng));
fid = fopen('ein_matlab.uv', 'w');
for i=1:len
    idu = floor(u(i) / uinc + 0.5) + 1;
    if idu < 1
        idu = idu + ng;
    end
    if idu > ng
        idu = idu - ng;
    end
    idv = floor(v(i) / vinc + 0.5) + 1;
    if idv < 1
        idv = idv + ng;
    end
    if idv > ng
        idv = idv - ng;
    end
    vis(i) = uvarr(idv, idu);
    fprintf(fid, '%f %f %f %f %f %f\n', ...
        u(i), v(i), 0, real(vis(i)), imag(vis(i)), 1.0);
    
    visarr(idv, idu) = uvarr(idv, idu);
    beamarr(idv, idu) = 1.0;
end
fclose(fid);

figure(2);
imagesc(abs(flipud(visarr)));
axis image;
colormap(gray);
colorbar();

img_dirt = ifft2(visarr);
img_dirt = fftshift(img_dirt);
img_dirt = real(img_dirt);
%img_dirt = img_dirt(ng4+1:ng4*3, ng4+1:ng4*3);
figure(401);
imagesc(flipud(img_dirt(ng4+1:ng4*3, ng4+1:ng4*3)));
%imagesc(flipud(img_dirt));
axis image;
colormap(gray);
colorbar();

beam_dirt = ifft2(beamarr);
beam_dirt = fftshift(beam_dirt);
beam_dirt = real(beam_dirt);
figure(501)
imagesc(flipud(beam_dirt(ng4+1:ng4*3, ng4+1:ng4*3)));
%imagesc(flipud(img_dirt));
axis image;
colormap(gray);
colorbar();

end

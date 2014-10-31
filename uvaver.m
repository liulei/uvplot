function uvaver()
% average visibility, reduce number of values in each grid
%

ng = 512;
fac = 8;
src = 'bk';
%src = '16B';
%src = 'ein';
uvname = strcat(src, '.uv');
ng = ng * fac;

c = 3E8;
res_mas = 0.1; % in mas
res_rad = res_mas * 1E-3 / 3600. / 180. * pi;
umax = 1.0 / res_rad / 2.0;
%umax = pi;

uinc = umax / ng * 2 / fac;
vinc = uinc;
ulimit = uinc * ng / 4;
vlimit = ulimit;
fprintf('Required uv: %f, res: %f\n', ulimit, res_mas);

offset = 3;

arr = importdata(uvname);

u = arr(:, 1);
v = arr(:, 2);

vis = complex(arr(:, offset + 1), arr(:, offset + 2));
weight = arr(:, offset + 3);

maxuv = max(u.^2 + v.^2);
maxuv = sqrt(maxuv);
minres = 1.0 / maxuv * 180. / pi * 3600. * 1000.;
fprintf('Provided max uv: %f, min res: %f\n', maxuv, minres);

ugrid = zeros(ng, ng);
vgrid = zeros(ng, ng);

visgrid = complex(zeros(ng, ng), zeros(ng, ng));
beamgrid = zeros(ng, ng);

nmeas = length(u);
for i=1:nmeas   
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
    beamarr(idv, idu) = beamarr(idv, idu) + 1.0;
end

end

%%
% This experiment does fiber tracking before and after smoothing
%
% Gunnar Atli Sigurdsson, 2013


dodec=load('directions253.txt','-ascii');
mix=ReadXml('dwi-b1500decimated.xml_clone_basisMixtures.xml'));
ind=ReadXml('dwi-b1500decimated.xml_clone_basisIndecies.xml'));
dir = 5;

mix = mix(:,end:-1:1, :, :);
ind = ind(:,end:-1:1, :, :);
mask = ReadXml('masktest1.xml');
mask = ~mask(:,end:-1:1, :);

dodec(:, 2) = -dodec(:, 2);


Nx = size(mix,1);
Ny = size(mix,2);
Nz = size(mix,3);
Nw = size(mix,4);

%% clean data

weights = zeros(Nx, Ny, Nz, Nw);
dirs = zeros(Nx, Ny, Nz, Nw, 3);

for a = 1:Nx
    disp(a/Nx)
    for b = 1:Ny
        for c = 1:Nz
        w = squeeze(mix(a,b,c,:));

        if all(isnan(w))
            continue; %skip voxel
        end
        
        ind1 = squeeze(ind(a,b,c,:));
        ind1(isnan(w) | ind1 == -1) = 0;
        w(isnan(w)) = 0;
        X = dodec(ind1+1, :);
        
        dirs(a,b,c,:,:) = X;
        weights(a,b,c,:) = w;
        
        end
    end	
end


%% apply white matter mask

weights = weights.*repmat(mask, [1 1 1 Nw]);

%% visualize

figure(1); clf;
slice = 2;
visualizeField(weights, dirs, slice);
%set(findall(0,'Type','Line'),'LineWidth',2);

%% get seeds

seedsimage=ReadXml('spatial_positions.xml');
[xs ys] = find(seedsimage(:,end:-1:1,2));
seeds = [xs ys 2*ones(size(xs))];


%% simple fiber tracking (follow max, interative steps)

angle = 70/180*pi;
minW = eps;
step = 0.1;

points = [];
clear fibers;
for i=1:size(seeds,1) % for all seeds
    fiber = [];
    for minus = 1:2
        p = seeds(i,:)';
        fiber = flipud(fiber);
        x = p(1); y = p(2); z = 2;
        [w,wi] = max(shiftdim(weights(x,y,z,:)));
        dir = (-1)^minus*shiftdim(dirs(x,y,z,wi,:));
        olddir = dir;
        try
            while w > minW % && abs(dot(olddir,dir)) > cos(angle)
                olddir = dir;
                p = p + step*dir;
                x = round(p(1)); y = round(p(2)); z = 2;
                [~,wi] = max(abs(shiftdim(dirs(x,y,z,:,:))*olddir));
                w = weights(x,y,z,wi);
                dir = shiftdim(dirs(x,y,z,wi,:));
                if norm(dir(1:2)) < step
                    dir = olddir;
                end
                if dot(dir,olddir) < 0
                    dir = -dir; % go forward
                end
                fiber = [fiber; p'];
            end
        catch err
            disp(err)
        end
    end
    fibers{i} = fiber;
    points = [points; fiber];
end

%% plot fibers

figure(1);
try
    delete(lines)
catch err
end
lines = plot(points(:,1), -points(:,2), 'b.', 'MarkerSize', 1);


%% smooth data

iters = 5;
smoothdirs = collectionsmooth(weights,dirs,[],pi/4,iters);

%% smooth data, recalculating coefficients at each step

iters = 10;
%smoothdirs = dirs;
for i = 1:iters
    smoothdirs = collectionsmooth(weights,smoothdirs,[],pi/4,3);
end

%% visualize smooth data

figure(2); clf;
slice = 2;
visualizeField(weights, smoothdirs, slice);


%% simple fiber tracking on smooth data (follow max, interative steps)

angle = 70/180*pi;
minW = eps;
step = 0.1;

points2 = [];
clear fibers2;
for i=1:size(seeds,1) % for all seeds
    fiber = [];
    for minus = 1:2
        p = seeds(i,:)';
        fiber = flipud(fiber);
        x = p(1); y = p(2); z = 2;
        [w,wi] = max(shiftdim(weights(x,y,z,:)));
        dir = (-1)^minus*shiftdim(smoothdirs(x,y,z,wi,:));
        olddir = dir;
        try
            while w > minW % && abs(dot(olddir,dir)) > cos(angle)
                olddir = dir;
                p = p + step*dir;
                x = round(p(1)); y = round(p(2)); z = 2;
                [~,wi] = max(abs(shiftdim(smoothdirs(x,y,z,:,:))*olddir));
                w = weights(x,y,z,wi);
                dir = shiftdim(smoothdirs(x,y,z,wi,:));
                if norm(dir(1:2)) < step
                    dir = olddir;
                end
                if dot(dir,olddir) < 0
                    dir = -dir; % go forward
                end
                fiber = [fiber; p'];
            end
        catch err
            disp(err)
        end
    end
    fibers2{i} = fiber;
    points2 = [points2; fiber];
end

%% plot smooth fibers

figure(2);
try
    delete(smoothlines)
catch err
end
smoothlines = plot(points2(:,1), -points2(:,2), 'b.', 'MarkerSize', 1);





%% smoothing on an ideal phantom
%
% Gunnar Atli Sigurdsson, 2013

Nx = 10;
Ny = 10;
Nz = 3;
Nw = 4;

%% generate data
% data is composed of two equally weighted directions corrupted by noise
% data is normalized to the unit sphere

weights = zeros(Nx, Ny, Nz, Nw);
rng('default');
 
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            r = rand;
            if rand > 0.5
                weights(i,j,k,1) = 0.5; % no split
            else
                weights(i,j,k,[1 2]) = 0.25; %split
            end
            
            if rand > 0.5
                weights(i,j,k,3) = 0.5; % no split
            else
                weights(i,j,k,[3 4]) = 0.25; %split
            end
        end
    end
end

%weights = weights + randn(size(weights))/5;
%weights(weights < 0) = 0;

dirs = zeros(Nx, Ny, Nz, Nw, 3);
dirs(:,:,:,[1 2],1) = 1;
dirs(:,:,:,[3 4],2) = 1;
dirs = dirs + randn(size(dirs))/4;
dirsnorm = sqrt(sum(dirs.^2, 5));
dirs(:,:,:,:,1) = dirs(:,:,:,:,1)./dirsnorm;
dirs(:,:,:,:,2) = dirs(:,:,:,:,2)./dirsnorm;
dirs(:,:,:,:,3) = dirs(:,:,:,:,3)./dirsnorm;

%% smooth

iters = 10; %50
slice = 2;
smoothdirs = collectionsmoothSlice(weights,dirs,[],[],iters,slice);

%% visualize original

figure(1)
clf;
visualizeField(weights, dirs, slice);
%axis([1 size(weights,1) -size(weights,2) -1])

% visualize smooth
figure(2)
clf;
visualizeField(weights, smoothdirs, slice);
%axis([1 size(weights,1) -size(weights,2) -1])


%% stats

cart2sphv = @(X) cart2sph(X(:,1), X(:,2), X(:,3));

tmp = dirs(2:end-1, 2:end-1, slice, :, :);
tmp = shiftdim(tmp, 4);
tmp = tmp(:, :)';

[t,~,~] = cart2sphv(tmp);
t = t.*(t >= 0) + (t + pi).*(t < 0);

figure(3)
subplot(2,1,1)
hist(t, 30)
axis tight
xlabel('Angle')
ylabel('Frequency')
title('No smoothing')


tmp = smoothdirs(2:end-1, 2:end-1, slice, :, :);
tmp = shiftdim(tmp, 4);
tmp = tmp(:, :)';

[t,~,~] = cart2sphv(tmp);
t = t.*(t >= 0) + (t + pi).*(t < 0);

subplot(2,1,2)
hist(t, 30)
axis tight
xlabel('Angle')
ylabel('Frequency')
title('With smoothing')

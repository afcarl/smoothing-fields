%%
% Loads phantoms from JIST outputs
% and displays them.
% Smoothing is performed for comparison
%
% Gunnar Atli Sigurdsson 2013

% raw data:

%% Load CFARI reconstructed phantom

%filepath2 = '/home/gunnar/JISTProjects/CFARIsmoothing/totalcrossingoutDispersionandConst/exp-0000/'; %noisy version
filepath2 = '/home/gunnar/JISTProjects/SPIE2014/outputs/lownoisevaryingphantom/exp-0000/';

mix=ReadXml(strcat(filepath2,'exp-0000-DAAA/Output_calc.xml_clone_basisMixtures.xml')); %was DAAA
ind=ReadXml(strcat(filepath2,'exp-0000-DAAA/Output_calc.xml_clone_basisIndecies.xml'));

dir = 5;

filepath = '/home/gunnar/Documents/MATLAB/cfariinterp/';
dodec=load([filepath 'directions253.txt'],'-ascii');

Nx = size(mix,1);
Ny = size(mix,2);
Nz = size(mix,3);
Nw = size(mix,4);

% clean data
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

%% visualize phantom after CFARI reconstruction

figure(10); clf;
slice = 5;
visualizeField(weights, dirs, slice);
%axis([1 Nx/2 -Ny/2 -1])


%% load raw data, do gaussian smoothing on individual channels

%[vol, hdr] = ReadXml('/home/gunnar/Documents/MATLAB/smoothingSPIE/varyingphantom.xml'); %more noise
[vol, hdr] = ReadXml('/home/gunnar/JISTProjects/SPIE2014/outputs/lownoisevaryingphantom/exp-0000/exp-0000-BA/Create_a_Crossing_Fiber_Pattern/Output_calc.xml');

sigma = 2;
kernel = fspecial('gaussian', [10 1], sigma);
smooth1 = zeros(size(vol));

for i=1:size(vol,4)
   tmp = squeeze(vol(:,:,:,i));
   tmp = convn(tmp, kernel, 'same');
   tmp = convn(tmp, shiftdim(kernel, -1), 'same');
   %tmp = convn(tmp, shiftdim(kernel, -2), 'same'); % only in plane
    
   smooth1(:,:,:,i) = tmp;
end

writeXml(smooth1, hdr, '/home/gunnar/JISTProjects/SPIE2014/input/gauss2lowvary');


%% Load Gaussian output

%filepath2 = '/home/gunnar/JISTProjects/SPIE2014/outputs/gauss4reconstructed/exp-0000/exp-0000-BAAA/';
%filepath2 = '/home/gunnar/JISTProjects/SPIE2014/outputs/gauss6reconstructed/exp-0000/exp-0000-BAAA/';
filepath2 = '/home/gunnar/JISTProjects/SPIE2014/outputs/gauss4lowvarreconstructed/exp-0000/exp-0000-BAAA/';

mix=ReadXml(strcat(filepath2,'gauss4lowvary.xml_clone_basisMixtures.xml'));
ind=ReadXml(strcat(filepath2,'gauss4lowvary.xml_clone_basisIndecies.xml'));

filepath = '/home/gunnar/Documents/MATLAB/cfariinterp/';
dodec=load([filepath 'directions253.txt'],'-ascii');

Nx = size(mix,1);
Ny = size(mix,2);
Nz = size(mix,3);
Nw = size(mix,4);

% clean data
gweights = zeros(Nx, Ny, Nz, Nw);
gdirs = zeros(Nx, Ny, Nz, Nw, 3);

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
        
        gdirs(a,b,c,:,:) = X;
        gweights(a,b,c,:) = w;
        
        end
    end	
end

%% visualize Gaussian

figure(5); clf;
slice = 5;
visualizeField(gweights, gdirs, slice);
axis([1 Nx/2 -Ny/2 -1])
print -dpng -r300 varygauss4

%% smooth data from cfari reconstructed (not gaussian)

iters = 10;
slice = 5;
smoothdirs = collectionsmoothSliceXY(weights,dirs,[],pi/4,iters,slice);

%% visualize after smoothing

figure(3); clf
slice = 5;
visualizeField(weights, smoothdirs, slice)
axis([1 Nx/2 -Ny/2 -1])







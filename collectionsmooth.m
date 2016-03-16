function [newDirs] = collectionsmooth(fractions, dirs, kernel, cutoff, iters)
% Smooths the orientations in 'dirs'
% This methods smooths the whole volume and is extremely slow in MATLAB
% please look at the implementation in JAVA Imaging Science Toolkit
%
% This method calculates the weights at each iteration
% 
% from:
% 'Smoothing Fields of Weighted Collections with Applications to Diffusion MRI Processing'
% SPIE Medical Imaging 2014
%
% Gunnar Atli Sigurdsson, 2013

Nx = size(fractions,1);
Ny = size(fractions,2);
Nz = size(fractions,3);
Nw = size(fractions,4);

fprintf('\n\n       ');
coeffsX = zeros(Nx, Ny, Nz, Nw^2);
for a = 1:Nx-1
    fprintf('\b\b\b\b\b\b\bX: %1.2f', a/Nx);
    for b = 1:Ny
        for c = 1:Nz
            wA = squeeze(fractions(a,b,c,:));
            wB = squeeze(fractions(a+1,b,c,:));

            if any(isnan(wA)) || sum(wA) < eps
                continue;
            end
            
            if any(isnan(wB)) || sum(wB) < eps
                continue;
            end

            XA = squeeze(dirs(a,b,c,:,:));
            XB = squeeze(dirs(a+1,b,c,:,:));

            [K, ~] = coeffs( wA, wB, XA, XB );
            coeffsX(a,b,c,:) = K;
        end
    end	
end

coeffsY = zeros(Nx, Ny, Nz, Nw^2);
for a = 1:Nx
    fprintf('\b\b\b\b\b\b\bY: %1.2f', a/Nx);
    for b = 1:Ny-1
        for c = 1:Nz
            wA = squeeze(fractions(a,b,c,:));
            wB = squeeze(fractions(a,b+1,c,:));

            if any(isnan(wA)) || sum(wA) < eps
                continue;
            end
            
            if any(isnan(wB)) || sum(wB) < eps
                continue;
            end

            XA = squeeze(dirs(a,b,c,:,:));
            XB = squeeze(dirs(a,b+1,c,:,:));

            [K, ~] = coeffs( wA, wB, XA, XB );
            coeffsY(a,b,c,:) = K;
        end
    end	
end

coeffsZ = zeros(Nx, Ny, Nz, Nw^2);
for a = 1:Nx
    fprintf('\b\b\b\b\b\b\bZ: %1.2f', a/Nx);
    for b = 1:Ny
        for c = 1:Nz-1
            wA = squeeze(fractions(a,b,c,:));
            wB = squeeze(fractions(a,b,c+1,:));

            if any(isnan(wA)) || sum(wA) < eps
                continue;
            end
            
            if any(isnan(wB)) || sum(wB) < eps
                continue;
            end

            XA = squeeze(dirs(a,b,c,:,:));
            XB = squeeze(dirs(a,b,c+1,:,:));

            [K, ~] = coeffs( wA, wB, XA, XB );
            coeffsZ(a,b,c,:) = K;
        end
    end	
end




%% smooth

newDirs = dirs; %reset

fprintf('\n\n       ')
for i=1:iters
fprintf('\b\b\b\b\b\b\bI: %1.2f', i/iters);

% smooth in x dir
oldDirs = newDirs; % copy constructor
for a = 2:Nx-1
    for b = 1:Ny
        for c = 1:Nz
            XA = squeeze(oldDirs(a-1,b,c,:,:));
            XB = squeeze(oldDirs(a,b,c,:,:));
            XC = squeeze(oldDirs(a+1,b,c,:,:));

            KAB = squeeze(coeffsX(a-1,b,c,:));
            KBC = squeeze(coeffsX(a,b,c,:));

            newDirs(a,b,c,:,:) = smooth1D( KAB, KBC, XA, XB, XC, kernel, cutoff);
        end
    end
end


% smooth in y dir
oldDirs = newDirs; % copy constructor
for a = 1:Nx
    for b = 2:Ny-1
        for c = 1:Nz
            XA = squeeze(oldDirs(a,b-1,c,:,:));
            XB = squeeze(oldDirs(a,b,c,:,:));
            XC = squeeze(oldDirs(a,b+1,c,:,:));

            KAB = squeeze(coeffsY(a,b-1,c,:));
            KBC = squeeze(coeffsY(a,b,c,:));

            newDirs(a,b,c,:,:) = smooth1D( KAB, KBC, XA, XB, XC, kernel, cutoff);
        end
    end
end

% smooth in z dir
oldDirs = newDirs; % copy constructor
for a = 1:Nx
    for b = 1:Ny
        for c = 2:Nz-1
            XA = squeeze(oldDirs(a,b,c-1,:,:));
            XB = squeeze(oldDirs(a,b,c,:,:));
            XC = squeeze(oldDirs(a,b,c+1,:,:));

            KAB = squeeze(coeffsZ(a,b,c-1,:));
            KBC = squeeze(coeffsZ(a,b,c,:));
            
            newDirs(a,b,c,:,:) = smooth1D( KAB, KBC, XA, XB, XC, kernel, cutoff);
        end
    end
end


end



end

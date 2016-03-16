%% script to smooth data in measurement domain
%
% Gunnar Atli Sigurdsson 2014


[vol, hdr] = ReadXml('/home/gunnar/JISTProjects/SPIE2014/input/varyphantom.xml');

sigma = 1;
kernel = fspecial('gaussian', [10 1], sigma);
smooth1 = zeros(size(vol));

for i=1:size(vol,4)
   tmp = squeeze(vol(:,:,:,i));
   tmp = convn(tmp, kernel, 'same');
   tmp = convn(tmp, shiftdim(kernel, -1), 'same');
   %tmp = convn(tmp, shiftdim(kernel, -2), 'same'); %only xy
    
   smooth1(:,:,:,i) = tmp;
end

writeXml(smooth1, hdr, '/home/gunnar/JISTProjects/SPIE2014/outputs/20140111_varyphantom/gauss/varyphantomG1');

sigma = 2;
kernel = fspecial('gaussian', [10 1], sigma);
smooth1 = zeros(size(vol));

for i=1:size(vol,4)
   tmp = squeeze(vol(:,:,:,i));
   tmp = convn(tmp, kernel, 'same');
   tmp = convn(tmp, shiftdim(kernel, -1), 'same');
   %tmp = convn(tmp, shiftdim(kernel, -2), 'same'); %only xy
    
   smooth1(:,:,:,i) = tmp;
end

writeXml(smooth1, hdr, '/home/gunnar/JISTProjects/SPIE2014/outputs/20140111_varyphantom/gauss/varyphantomG2');

sigma = 4;
kernel = fspecial('gaussian', [10 1], sigma);
smooth1 = zeros(size(vol));

for i=1:size(vol,4)
   tmp = squeeze(vol(:,:,:,i));
   tmp = convn(tmp, kernel, 'same');
   tmp = convn(tmp, shiftdim(kernel, -1), 'same');
   %tmp = convn(tmp, shiftdim(kernel, -2), 'same'); %only xy
    
   smooth1(:,:,:,i) = tmp;
end

writeXml(smooth1, hdr, '/home/gunnar/JISTProjects/SPIE2014/outputs/20140111_varyphantom/gauss/varyphantomG4');

sigma = 8;
kernel = fspecial('gaussian', [10 1], sigma);
smooth1 = zeros(size(vol));

for i=1:size(vol,4)
   tmp = squeeze(vol(:,:,:,i));
   tmp = convn(tmp, kernel, 'same');
   tmp = convn(tmp, shiftdim(kernel, -1), 'same');
   %tmp = convn(tmp, shiftdim(kernel, -2), 'same'); %only xy
    
   smooth1(:,:,:,i) = tmp;
end

writeXml(smooth1, hdr, '/home/gunnar/JISTProjects/SPIE2014/outputs/20140111_varyphantom/gauss/varyphantomG8');

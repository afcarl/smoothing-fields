function [ newXB ] = smooth1D( KAB, KBC, XA, XB, XC, W, cutoff)
%[ newXB ] = smooth1D( KAB, KBC, XA, XB, XC, W, cutoff)
%   KAB is the weights between voxel A and voxel B
%   KBC is the weights between voxel A and voxel B
%   XA are the orientations of voxel A
%   XB are the orientations of voxel B
%   XC are the orientations of voxel C
%   W is an optional voxel kernel
%   cutoff is an optional threshold in radians that specifies a distance
%   beyond which to ignore weights.

if isempty(W)
    W = [1 2 1]/4;
end

cutoff = cos(cutoff);
if isempty(cutoff)
    cutoff = inf; % no cutoff
end

% smooth each orientation of voxel B
N = size(XB,1);
newXB = zeros(size(XB));
for i = 1:N
    % XB(i) <-- weighted average of all directions related to XB(i)
    this = XB(i,:)';
    
    %from voxel A, (coeffs are from A to B)
    w0 = KAB(i:N:end);
    
    %from voxel B, (coeffs are from B to C)
    w2 = KBC((i-1)*N+1:i*N);
    
    % cutoff directions in A
    dist = abs(sum(XA*this,2));
    w0 = w0.*(dist < cutoff);
    
    % cutoff directions in B
    dist = abs(sum(XC*this,2));
    w2 = w2.*(dist < cutoff);
    
    X = [XA; this'; XC];

    % weigh by weight mask (sum cancels out fractions)
    w = [W(1)*w0/nansum(w0); W(2); W(3)*w2/nansum(w2)];
    newXB(i,:) = simplesphmeanw(X,w,XB(i,:));
    
    if sum(w) < eps || any(isnan(newXB(i,:)))
        newXB(i,:) = XB(i,:);
    end

end

end


function [X1] = simplesphmeanw(X, w, est)
    signs = sign(X*shiftdim(est));
    X = X.*[signs signs signs];

    X1 = sum(X.*[w w w],1);
    X1 = X1/norm(X1);
end


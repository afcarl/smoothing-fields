function [ K, cost ] = coeffs( wA, wB, XA, XB )
%% coeffs Calculates coefficients
%   wA the fractions of the start voxel
%   wB the fractions of the end voxel
%   XA the orientations in the start voxel
%   XB the orientations in the end voxel
%   K(1) is the weight from orientation 1 of A
%   to orientation 1 of B.
%   K(2) is the cost from orientation 1 of A
%   to orientation 2 of B.
%
%   Gunnar Atli Sigurdsson, 2013

Ng = size(wB,1);
Nf = size(wA,1);

wA = wA/sum(wA); %normalize for linprog
wB = wB/sum(wB);

F = zeros(Ng, Ng*Nf);
for i=1:Nf
    F(:, (i-1)*Ng+1:i*Ng) = eye(Ng);
end

S = zeros(Ng, Ng*Nf);
iS = 1;
for i=1:Ng
    S(i,iS:iS+Nf-1) = ones(size(wA'));
    iS = iS+Nf;
end

D = zeros(Nf*Ng,1);
iD = 1;
for i=1:Nf
    for j=1:Ng
        vec = XB(j,:) - XA(i,:);
        vec2 = -XB(j,:) - XA(i,:);   
        k = min(norm(vec), norm(vec2)); %chord
        D(iD) = 2*asin(k/2); %arclength
        D(iD) = D(iD)^2;
        iD = iD + 1;
    end
end

%Fs = sparse(F);
%Ss = sparse(S);

K = linprog(D,S,wA,F,wB,zeros(size(D)),[],[], ...
    optimset('Display', 'off', 'LargeScale', 'off', 'simplex', 'on'));

%K = linprog(D,[],[],[F; S],[wB; wA],zeros(size(D)),[],[], ...
%    optimset('Display', 'off', 'LargeScale', 'off', 'simplex', 'on'));

cost = D'*K;

end


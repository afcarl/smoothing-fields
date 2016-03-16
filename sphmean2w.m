function [ meanX, err, var ] = sphmean2w( X, w )
%SPHMEAN Calculates the mean on a sphere
%   Assumes X = -X (Direction)
%   w should sum up to 1 and be non-negative

% % step 0, normalize w
% wsum = sum(w);
% if wsum < eps
%     w = 0;
% else
%     w = w/wsum;
% end

% step 1, find theta
cart2sphv = @(X) cart2sph(X(:,1), X(:,2), X(:,3));

% resolve ambiguity, wrap semicircle
[t0,p0,~] = cart2sphv(X);
t = 2*t0;
[X0,Y0,Z0] = sph2cart(t,zeros(size(t)),ones(size(t)));
%meanXY = sum([X0.*w Y0.*w Z0.*w],1)'/sum(w);
meanXY = sum([X0.*w Y0.*w Z0.*w],1)';
meanXY = meanXY/norm(meanXY,2);

% step 2, find phi
p = 2*p0.*(((abs(t0) < pi/2)*2) - 1);
[X1,Y1,Z1] = sph2cart(p,zeros(size(t)),ones(size(t)));
%meanXZ = sum([X1.*w Y1.*w Z1.*w],1)'/sum(w);
meanXZ = sum([X1.*w Y1.*w Z1.*w],1)';
meanXZ = meanXZ/norm(meanXZ,2);

% reassemble
[t0,~,~] = cart2sphv(meanXY');
[p0,~,~] = cart2sphv(meanXZ');

t = t0/2;
p = p0/2;
[X0,Y0,Z0] = sph2cart(t,p,1);
meanX = [X0 Y0 Z0]';
meanX = meanX/norm(meanX,2);

% calc error
[t0,p0,~] = cart2sphv(meanX');

err = 0;
n = size(X,1);
%for j=1:n
%    [t,p,~] = cart2sphv([X(j,:); -X(j,:)]);
%    errt = pi/2 - abs(abs(mod(t,2*pi) - mod(t0,2*pi)) - pi/2);
%    errp = pi/2 - abs(abs(mod(p,2*pi) - mod(p0,2*pi)) - pi/2);
%    err0 = min(sqrt(errt.^2+errp.^2));
%    err = err + err0;
%end

%rough error test (euclidean distance)
ntest = ones(n,1)*meanX';
err = sum(sqrt(min(sum((X - ntest).^2, 2), sum((-X - ntest).^2, 2))));
var = sum((min(sum((X - ntest).^2, 2), sum((-X - ntest).^2, 2))));

end


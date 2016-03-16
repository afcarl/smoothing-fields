function [] = visualizeField(fractions, dirs, slice)
%
% Visualizes a field of orientations
%
% Gunnar Atli Sigurdsson, 2013

hold all;
scale = 1;
for i=1:size(fractions,1)
    for j=1:size(fractions,2)
        m = fractions(i,j,slice,:);
        m = m/sum(m(isfinite(m)));
        for k=1:size(dirs,4)
            if(fractions(i,j,slice,k)>=0)
                % dim 1 is right to left
                % dim 2 is anterior to posterior
                % dim 3 is inferior to superior
                h=plot(i+scale*[-1 1]*dirs(i,j,slice,k,1)*m(k),...
                    -j-scale*[-1 1]*dirs(i,j,slice,k,2)*m(k),'k');
                vv = dirs(i,j,slice,k,:);
                [th,ph,r]=cart2sph(vv(1),vv(2),vv(3));
                [x,y,z]=sph2cart(th,ph,1);
                vv=abs([x y z]);
                set(h,'color',squeeze(min(1,max(0,abs(vv)))))
            end
        end
    end
end
axis equal tight;
axis off; 
set(gcf,'color','white')

end

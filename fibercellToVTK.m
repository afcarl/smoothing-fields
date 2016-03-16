function [ ] = fibercellToVTK( fibers, filename )
% writes all fibers in a cell to a VTK file
% fibers is a cell where fibers{i} is an array of points
% 
% Gunnar Atli Sigurdsson

N = numel(fibers);

pointCount = 0;
for i=1:N
    pointCount = pointCount + size(fibers{i},1);
end

% Write Header
fid = fopen(filename, 'w');
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'Curves\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET POLYDATA\n');
fprintf(fid, 'POINTS %d float\n', pointCount);

% Write points
for i=1:N
    fiber = fibers{i};
    for j=1:size(fiber, 1)
        fprintf(fid, '%f %f %f\n', fiber(j,:));
    end
end

% Write indices
fprintf(fid, 'LINES %d %d\n', N, N+pointCount);
ind = 0;
for i=1:N
    fprintf(fid, '%d', size(fibers{i},1));
    for j=1:size(fibers{i},1)
        fprintf(fid, ' %d', ind);
        ind = ind+1;
    end
    fprintf(fid, '\n');
end

% Write Lines
fprintf(fid, '\nCELL_DATA %d\n', N);
fprintf(fid, 'SCALARS Length float %d\n', 1); 
fprintf(fid, 'LOOKUP_TABLE default\n');

for i=1:N
    fprintf(fid, '%f\n', size(fibers{i},1)); 
end


fclose(fid);

end


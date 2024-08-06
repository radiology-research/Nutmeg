function [] = om_save_tri_flip(filename,points,faces,nrmls)

% OM_SAVE_TRI   Save .tri file
%
%   SYNTAX
%       [] = OM_SAVE_TRI(FILENAME,POINTS,FACES,NORMALS)
%
%   Created by Alexandre Gramfort on 2008-03-09.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailed information

if nargin<4 || isempty(nrmls)
    nrmls = om_normals(points,faces);
end

fid = fopen(filename,'w');
npoints = size(points,1);
fprintf(fid,'- %g\n',npoints);
fprintf(fid,'%g %g %g %g %g %g\n',[points , nrmls]');
nfaces = size(faces,1);
faces = faces-1;
faces = faces(:,[2 1 3]);  % flip triangles for openmeeg
fprintf(fid,'- %g %g %g\n', [nfaces nfaces nfaces]);
fprintf(fid,'%g %g %g\n',faces');

fclose(fid);

end %  function

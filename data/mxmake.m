function [A,v,w]=mxmake(file, nr, nc, rc)

if ( nargin < 4 ) rc = "real" ; end

A = randn(nr, nc) ;
v = randn(nc,1) ;
r = 1 ;

if ( strcmp(rc, "complex") )
  A += j*randn(nr, nc) ;
  v += j*randn(nc, 1) ;
  r = 2 ;
end

w = A*v ;

fid = fopen(file, "w") ;

mxwrite(A, fid) ;
vcwrite(v, fid) ;
vcwrite(w, fid) ;

fclose(fid) ;
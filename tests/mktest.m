nr = 13 ; nc = 22 ; cm = 1 ;

for cm=1:2
A = randn(nr,nc) ;
x = randn(nc,1) ;
if ( cm == 2 ) 
  A += j*randn(nr,nc) ;
  x += j*randn(nc,1) ;
end

b = A*x ;

fid = fopen(["vectors-" int2str(cm) ".dat"], "w") ;

fprintf(fid, "%d %d %d\n", nr, nc, cm) ;
dat = A.' ;
if ( cm == 2 ) 
  dat = [real(dat(:)) imag(dat(:))]' ;
end
fprintf(fid, "%1.16e\n", dat) ;

fprintf(fid, "%d %d\n", length(x), cm) ;
dat = x(:) ;
if ( cm == 2 ) 
  dat = [real(dat(:)) imag(dat(:))]' ;
end
fprintf(fid, "%1.16e\n", dat) ;

fprintf(fid, "%d %d\n", length(b), cm) ;
dat = b(:) ;
if ( cm == 2 ) 
  dat = [real(dat(:)) imag(dat(:))]' ;
end
fprintf(fid, "%1.16e\n", dat) ;

fclose(fid) ;

endfor

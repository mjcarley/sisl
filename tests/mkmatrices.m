nr = 13 ; nc = 15 ;

for cm = 1:2 ;

A = randn(nr, nc) ;
B = randn(nc, nr) ;

if ( cm == 2 ) 
  A += j*randn(nr,nc) ;
  B += j*randn(nc,nr) ;
end

D = A*B ;

fid = fopen(["matrices-" int2str(cm) ".dat"], "w") ;

fprintf(fid, "%d %d %d\n", nr, nc, cm) ;
dat = A.' ;
if ( cm == 2 ) 
  dat = [real(dat(:)) imag(dat(:))]' ;
end
fprintf(fid, "%1.16e\n", dat) ;

fprintf(fid, "%d %d %d\n", nc, nr, cm) ;
dat = B.' ;
if ( cm == 2 ) 
  dat = [real(dat(:)) imag(dat(:))]' ;
end
fprintf(fid, "%1.16e\n", dat) ;

fprintf(fid, "%d %d %d\n", nr, nr, cm) ;
dat = D.' ;
if ( cm == 2 ) 
  dat = [real(dat(:)) imag(dat(:))]' ;
end
fprintf(fid, "%1.16e\n", dat) ;

fclose(fid) ;

endfor

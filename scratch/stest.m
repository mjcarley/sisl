 nr = 64 ; nc = 128 ;

 i = floor(nr*rand(nr,1)) ; j = floor(nc*rand(nr,1)) ;

 x = rand(nr,1) ; y = rand(nr,1) ;
 z = rand(nc,1) + J*rand(nc,1) ;

 A = sparse(i+1, j+1, x+J*y, nr, nc, "unique") ;

 fid = fopen("test.dat", "w") ;

 fprintf(fid, "%d %d %1.16e %1.16e\n", [i(:) j(:) x(:) y(:)]') ;
% fprintf(fid, "%d %d %1.16e\n", [i(:) j(:) x(:)]') ;

 fprintf(fid, "0 0 %1.16e %1.16e\n", [real(z(:)), imag(z(:))]') ;

 fclose(fid) ;

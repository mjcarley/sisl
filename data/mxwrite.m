function mxwrite(A, f)

[nr, nc] = size(A) ;

if ( isreal(A) )
  fprintf(f, "%d %d %d\n", nr, nc, 1) ;
  for i=1:nr
    fprintf(f, "%1.16e ", A(i,:)) ;
  end
  fprintf(f, "\n") ;
else
  fprintf(f, "%d %d %d\n", nr, nc, 2) ;
  dat = zeros(1, 2*nc) ;
  for i=1:nr
    dat(:,1:2:end) = real(A(i,:)) ;
    dat(:,2:2:end) = imag(A(i,:)) ;    
    fprintf(f, "%1.16e ", dat) ;
  end
  fprintf(f, "\n") ;  
end

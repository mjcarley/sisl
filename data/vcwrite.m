function vcwrite(v, f)

[nr, nc] = size(v) ;
ne = max([nr nc]) ;

if ( isreal(v) )
  fprintf(f, "%d %d\n", ne, 1) ;
  fprintf(f, "%1.16e ", v) ;
  fprintf(f, "\n") ;
else
  fprintf(f, "%d %d\n", ne, 2) ;
  dat = zeros(1, 2*ne) ;
  dat(:,1:2:end) = real(v) ;
  dat(:,2:2:end) = imag(v) ;
  fprintf(f, "%1.16e ", dat) ;
  fprintf(f, "\n") ;  
end

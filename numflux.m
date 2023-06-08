function f=numflux(a,b)

LF=.25*(a.^2 + b.^2)- .5*max(abs(a),abs(b)).*(b-a);

f=.5*a.^2.* (a>0).*(b>0)  + LF .* (1 - (a>0).*(b>0) );

   
    
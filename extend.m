function y= extend(v, idx, xm, jumpidx,jumpsize, Lip,M,N,isdecreasing,dgn)
dx=xm(2)-xm(1);
umin=min(v);
umax=max(v);

y=zeros(1,N);


LM=max(max(max(Lip)),1);
Lip=LM*ones(size(Lip));
if idx==1
    y(1:jumpidx(1))=v(1:jumpidx(1));
    rlocslop=(v(jumpidx(1))-v(jumpidx(1)-1))/dx;
    
    if rlocslop<max(Lip(2:M+1))
        rspace=max(jumpsize(1)/(max(Lip(2:M+1))-rlocslop),0);
    else
        rspace=xm(end);
    end
    if isdecreasing(2) ==1
        rspace=.5*dgn;
    end
    rspace=min(rspace,1);
    %Radical measures
    %   rlocslop=LM;
    
    y(jumpidx(1)+1:end)=v(jumpidx(1))+rlocslop*(xm(jumpidx(1)+1:end)-xm(jumpidx(1)));
    y(jumpidx(1)+1:end)=y(jumpidx(1)+1:end)+(max(Lip(2:M+1))-rlocslop).*(xm(jumpidx(1)+1:end)-xm(jumpidx(1))-rspace).*(xm(jumpidx(1)+1:end)>xm(jumpidx(1))+rspace); 
    
    %rlocslop
   % max(Lip(2:M+1))
    %sum(jumpsize)
    y=min(y, umax+sum(jumpsize)+LM );
    
elseif idx==M+1
    y(jumpidx(M)+1:N)=v(jumpidx(M)+1:N);
    
    llocslop=(v(jumpidx(M)+2)-v(jumpidx(M)+1))/dx;
    if llocslop<max(Lip(1:M))
        lspace=max(jumpsize(M)/(max(Lip(1:M))-llocslop),0);
    else
        lspace=xm(end);
    end
     if isdecreasing(M) ==1
        lspace=.5*dgn;
    end
    lspace=min(lspace,1);
     %Radical measures
     %  llocslop=LM;
    
    y(1:jumpidx(M))=v(jumpidx(M)+1)+llocslop*(xm(1:jumpidx(M))-xm(jumpidx(M)+1));
    y(1:jumpidx(M))=y(1:jumpidx(M))+(max(Lip(2:M+1))-llocslop).*(xm(1:jumpidx(M))-xm(jumpidx(M)+1)+lspace).*(xm(1:jumpidx(M))<xm(jumpidx(M)+1)-lspace); 
   
    y=max(y, umin-sum(jumpsize)-LM );
    
    
    
    
else
    y(jumpidx(idx-1)+1:jumpidx(idx))=v(jumpidx(idx-1)+1:jumpidx(idx));
    
    ls=max(Lip(1:idx-1));
    rs=max(Lip(idx+1:end));
    
    llocslop=(v(jumpidx(idx-1)+2)-v(jumpidx(idx-1)+1))/dx;
    if llocslop<ls
        lspace=max(jumpsize(idx-1)/(ls-llocslop),0);
    else
        lspace=xm(end);
    end
    if isdecreasing(idx-1) ==1
        lspace=.5*dgn;
    end
    lspace=min(lspace,1);
    
    rlocslop=(v(jumpidx(idx))-v(jumpidx(idx)-1))/dx;
    if rlocslop<rs
        rspace=max(jumpsize(idx)/(rs-rlocslop),0);
    else
        rspace=xm(end);
    end
    
    if isdecreasing(idx+1) ==1
        rspace=.1;
    end
    rspace=min(rspace,1);
      %Radical measures
      % rlocslop=LM;
      % llocslop=LM;
    
    y(1:jumpidx(idx-1))=v(jumpidx(idx-1)+1)+llocslop*(xm(1:jumpidx(idx-1))-xm(jumpidx(idx-1)+1));
    
    y(1:jumpidx(idx-1))=y(1:jumpidx(idx-1))+(ls-llocslop)*(xm(1:jumpidx(idx-1))-xm(jumpidx(idx-1)+1)+lspace).*(xm(1:jumpidx(idx-1))<xm(jumpidx(idx-1)+1)-lspace); 
 
    y(jumpidx(idx)+1:end)=v(jumpidx(idx))+rlocslop*(xm(jumpidx(idx)+1:end)-xm(jumpidx(idx)));
    y(jumpidx(idx)+1:end)=y(jumpidx(idx)+1:end)+(rs-rlocslop).*(xm(jumpidx(idx)+1:end)-xm(jumpidx(idx))-rspace).*(xm(jumpidx(idx)+1:end)>xm(jumpidx(idx))+rspace); 
    
    y=min(y, umax+sum(jumpsize(idx:M))+LM );
    y=max(y, umin-sum(jumpsize(1:idx-1))-LM);
end
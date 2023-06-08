function uhat=recon(u,N,M,ul)
    uhat=zeros(2,N,M);
    uhat(1,1,:)=ul;
    for j=1:N-1
        uhat(2,j,:)=u(:,j);
        uhat(1,j+1,:)=u(:,j);
    end
    uhat(2,N,:)=u(:,N);
    
    
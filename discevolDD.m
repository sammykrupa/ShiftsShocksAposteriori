function hhat=discevolDD(hhat0,uhato,umid,ssp,uft,num,dt,lNbr,rNbr,idxset,x)
hhat=hhat0;

    for idx=idxset
            %find(x<hhat0(idx),1,'last');
            k=sum(ssp(lNbr(idx),1:num(lNbr(idx))-1)< hhat0(idx));
            uleft=uft(lNbr(idx),1+k);
            k=sum(ssp(rNbr(idx),1:num(rNbr(idx))-1)< hhat0(idx));
            uright=uft(rNbr(idx),1+k);
            
            hhat(idx)=hhat0(idx) + dt/2*(.5*(uleft+uright) );

            k=sum(ssp(lNbr(idx),1:num(lNbr(idx))-1)< hhat(idx));
            uleft=uft(lNbr(idx),1+k);
            k=sum(ssp(rNbr(idx),1:num(rNbr(idx))-1)< hhat(idx));
            uright=uft(rNbr(idx),1+k);
           
            hhat(idx)=hhat0(idx) + dt*(.5*(uleft+uright)  );
    end
    
    hhat=hhat(:,idxset);
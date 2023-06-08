function hhat=discevolLD(hhat0,uhato,umid,ssp,uft,num,dt,lNbr,rNbr,idxset,x)
hhat=hhat0;
dx=x(2)-x(1);
    for idx=idxset
            %find(x<hhat0(idx),1,'last');
            cell=find(x<hhat0(idx),1,'last');

            k=sum(  ssp(lNbr(idx),1:num(lNbr(idx))-1)< hhat0(idx)    );
            uleft=uft(lNbr(idx),1+k);
            uright=uhato(1,cell,rNbr(idx)) + (uhato(2,cell,rNbr(idx))-uhato(1,cell,rNbr(idx)))/dx*(hhat0(idx)-x(cell)) ;
            hhat(idx)=hhat0(idx) + dt/2*(.5*(uleft+uright) );

            cell=find(x<hhat(idx),1,'last');
            k=sum(ssp(lNbr(idx),1:num(lNbr(idx))-1)< hhat(idx));
            uleft=uft(lNbr(idx),1+k);
            uright=umid(1,cell,rNbr(idx)) + (umid(2,cell,rNbr(idx))-umid(1,cell,rNbr(idx)))/dx*(hhat(idx)-x(cell)) ;
            hhat(idx)=hhat0(idx) + dt*(.5*(uleft+uright)  );
    end
    
    hhat=hhat(:,idxset);
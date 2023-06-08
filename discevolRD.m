function hhat=discevolRD(hhat0,uhato,umid,ssp,uft,num,dt,lNbr,rNbr,idxset,x)
hhat=hhat0;
dx=x(2)-x(1);
    for idx=idxset
            %find(x<hhat0(idx),1,'last');
            cell=find(x<hhat0(idx),1,'last');
            uleft=uhato(1,cell,lNbr(idx)) + (uhato(2,cell,lNbr(idx))-uhato(1,cell,lNbr(idx)))/dx*(hhat0(idx)-x(cell)) ;
            k=sum(ssp(rNbr(idx),1:num(rNbr(idx))-1)< hhat0(idx));
            uright=uft(rNbr(idx),1+k);
            
            hhat(idx)=hhat0(idx) + dt/2*(.5*(uleft+uright) );

            cell=find(x<hhat(idx),1,'last');
            uleft=umid(1,cell,lNbr(idx)) + (umid(2,cell,lNbr(idx))-umid(1,cell,lNbr(idx)))/dx*(hhat(idx)-x(cell)) ;
            k=sum(ssp(rNbr(idx),1:num(rNbr(idx))-1)< hhat(idx));
            uright=uft(rNbr(idx),1+k);
           
            hhat(idx)=hhat0(idx) + dt*(.5*(uleft+uright)  );
    end
    
    hhat=hhat(:,idxset);
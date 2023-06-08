function hhat=discevol(hhat0,uhato,umid,dt,lNbr,rNbr,idxset,x)
hhat=hhat0;
dx=x(2)-x(1);
    for idx=idxset
            %find(x<hhat0(idx),1,'last');
            cell=find(x<hhat0(idx),1,'last');
            uleft=uhato(1,cell,lNbr(idx)) + (uhato(2,cell,lNbr(idx))-uhato(1,cell,lNbr(idx)))/dx*(hhat0(idx)-x(cell)) ;
            uright=uhato(1,cell,rNbr(idx)) + (uhato(2,cell,rNbr(idx))-uhato(1,cell,rNbr(idx)))/dx*(hhat0(idx)-x(cell)) ;
            hhat(idx)=hhat0(idx) + dt/2*(.5*(uleft+uright) );

            cell=find(x<hhat(idx),1,'last');
            uleft=umid(1,cell,lNbr(idx)) + (umid(2,cell,lNbr(idx))-umid(1,cell,lNbr(idx)))/dx*(hhat(idx)-x(cell)) ;
            uright=umid(1,cell,rNbr(idx)) + (umid(2,cell,rNbr(idx))-umid(1,cell,rNbr(idx)))/dx*(hhat(idx)-x(cell)) ;
            hhat(idx)=hhat0(idx) + dt*(.5*(uleft+uright)  );
    end
    
    hhat=hhat(:,idxset);
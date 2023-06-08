function hhat=discevolscal(hhat0,uhato,umid,dt,lNbr,rNbr,idxset,x)
dx=x(2)-x(1);
    
            find(x<hhat0,1,'last');
            wcell=find(x<hhat0,1,'last');
         
            uleft=uhato(1,wcell,lNbr) + (uhato(2,wcell,lNbr)-uhato(1,wcell,lNbr))/dx*(hhat0-x(wcell)) ;
            uright=uhato(1,wcell,rNbr) + (uhato(2,wcell,rNbr)-uhato(1,wcell,rNbr))/dx*(hhat0-x(wcell)) ;
            hhat=hhat0 + dt/2*(.5*(uleft+uright) );

            wcell=find(x<hhat,1,'last');
 
            uleft=umid(1,wcell,lNbr) + (umid(2,wcell,lNbr)-umid(1,wcell,lNbr))/dx*(hhat-x(wcell)) ;
            uright=umid(1,wcell,rNbr) + (umid(2,wcell,rNbr)-umid(1,wcell,rNbr))/dx*(hhat-x(wcell)) ;
            hhat=hhat0+ dt*(.5*(uleft+uright)  );
    
    
function R=compres(uold,unew,dx,dt)

        umid=.5*(uold+unew);
        Res(1,:)=(unew(1,:)-uold(1,:))/dt+uold(1,:).*(uold(2,:)-uold(1,:))/dx;
        Res(2,:)=(unew(2,:)-uold(2,:))/dt+uold(2,:).*(uold(2,:)-uold(1,:))/dx;
        Res(3,:)=4*((unew(1,:)-uold(1,:))/dt+umid(1,:).*(umid(2,:)-umid(1,:))/dx);
        Res(4,:)=4*((unew(2,:)-uold(2,:))/dt+umid(2,:).*(umid(2,:)-umid(1,:))/dx);
        Res(5,:)=(unew(1,:)-uold(1,:))/dt+unew(1,:).*(unew(2,:)-unew(1,:))/dx;
        Res(6,:)=(unew(2,:)-uold(2,:))/dt+unew(2,:).*(unew(2,:)-unew(1,:))/dx;
    
        
        
        R(:)=sum(Res(:,:).^2)/12;
      
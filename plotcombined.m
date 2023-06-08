function x=plotcombined(x,xm,hhat,uhat,M,uhatn,rev,pe)

% 
%clf

%subplot(1,2,1)
hold on 
for j=1:M+1
plot([x(1:end-1);x(2:end)],uhatn(:,:,j).*[rev(j,:);rev(j,:)],'r', 'LineWidth',3)
end

loc = find( x<=hhat(1));
plot(xm(loc),uhat(1,loc),'b','LineWidth',3)
for i=2:M 
    loc = find( (x(1:end-1)<=hhat(i) ) .* ( x(2:end) >= hhat(i-1) ) );
    plot(xm(loc),uhat(i,loc),'b','LineWidth',3)
end
   plot([x(1:end-1);x(2:end)],zeros(size([x(1:end-1);x(2:end)]  )  ),'k', 'LineWidth',3)
   loc = find(  x(2:end) >= hhat(M)-pe(M)  );
    plot(xm(loc),uhat(M+1,loc),'r','LineWidth',3)
    loc = find(  x(2:end) >= hhat(M)  );
    plot(xm(loc),uhat(M+1,loc),'b','LineWidth',3)
    x=1;
   
hold off



 
% subplot(1,2,2)
% hold on
% 
% for j=1:M+1
% plot(xm,R(j,:).*rev(j,:),'r')
% end
% 
% hold off
 
function [uhatn]=nwwn03(N,T)
 
%Input parameters of simulation
%N= number of cells
%T= final time


tic


%Number of cells
%N=3200;  is input argument

%Length of domain
L=1;


%sumsource=zeros(1,3);

ts=0; count=1;

%Compute mesh width
dx=L/N;

%Choose delta (potentially different from dx);
deltageneric=dx^(1/2);

%Create mesh and midpoints
x=linspace(0,L,N+1);
xm=.5*(x(1:end-1)+x(2:end));

%initialize list of collisions that have taken place in the computer
collist=zeros(1, 13);

%initialize list of collisions we know must have happened in \hat \psi
certlist=zeros(1, 2);


%nearly non-decreasing and rapidly decreasing piece (there might not be a discontinuity here because we have not yet discretized the rapidly decreasing pieces )
% starting point of h-curve is right boundary of cell jumpidx

jumpidx=[N/4 N/2 ceil(N*5/8)];

%M is the number of h_i curves (and \hat h_i curves) corresponding to
%`large' shocks and `boundary shocks'
M=length(jumpidx);

%Specify initial data
%Interesting initial data
% v=zeros(M+1,N);
% v(1,:)= xm;
% v(2,:)= 1/2- xm;
% v(3,:)= -1 + 2*xm;
% v(4,:)=-ones(1,length(xm));
% isdecreasing=[0, 1, 0, 0];


% second part is increasing initial data
%v1=@(z)  2.5 ;
%v2=@(z)  2.5 -5*(z-1/5);
%v3=@(z)  2.5 -5*(z-1/5)+6*(z-2/5);
%v4=@(z)  0;

%hand-curated list of which intervals are rapidly decreasing
%isdecreasing=[0, 1, 0, 0];

% all parts increasing initial data
v1=@(z) 3 ;
v2=@(z)  1+2*z;
v3=@(z)  z;
v4=@(z)  0;
%hand-curated list of which intervals are rapidly decreasing
isdecreasing=[0, 0, 0, 0];

v=zeros(M+1,N);
v(1,:)= v1(xm);
v(2,:)= v2(xm);
v(3,:)= v3(xm);
v(4,:)= v4(xm);




% build initial data as 'global' function
% nice for plotting
v0=zeros(1,N);
v0(1:jumpidx(1))=v(1,1:jumpidx(1));
v0(jumpidx(M)+1:end)=v(M+1,jumpidx(M)+1:end);
for i=2:M
  v0(jumpidx(i-1)+1:jumpidx(i))=v(i,jumpidx(i-1)+1:jumpidx(i));  
end
    

%plot(v0)
%pause


%Initilize total value of the L^2 error estimator (the squared space-time L^2 norm of
%\mathcal{R}_{\hat \psi})
Etot=0;

%the initial residual error from \hat\Lambda   
LhatE=0;

%-----------------------------------------------------------------------
%Process initial data
%Jumps adjacent to the first and last cell are currently not allowed


%At time t=0, we start with M h-curves (which are shocks, potentially of zero size). Here, we are not including the `front tracking' shocks.


%Determine initial size of jumps
jumpsize=zeros(1,M);

jumpsize(1)=v1(x(jumpidx(1)+1))-v2(x(jumpidx(1)+1));
jumpsize(2)=v2(x(jumpidx(2)+1))-v3(x(jumpidx(2)+1));
jumpsize(3)=v3(x(jumpidx(3)+1))-v4(x(jumpidx(2)+1));

%Check that there are no upward jumps
if min(min(jumpsize)) < 0
    disp('error upward jump')
  %  nc
end


%Create list with numbers of nearly non-decreasing pieces
incl=[];
for idx=1:M+1
    if isdecreasing(idx)==0
        incl=[incl,idx];
    end
end









%Compute Lipschitz constants of nearly non-decreasing pieces
Lip=zeros(1,M+1);

for i=1:M+1
    if i==1
    Lip(i) = max ((1-isdecreasing(i))*(v(i,2:jumpidx(i))-v(i,1:jumpidx(i)-1)))/dx;
    elseif i==M+1
    Lip(i) = max ((1-isdecreasing(i))*(v(i,jumpidx(i-1)+1:end)-v(i,jumpidx(i-1):end-1)))/dx;  
    else
    Lip(i) = max ((1-isdecreasing(i))*(v(i,jumpidx(i-1)+1:jumpidx(i))-v(i,jumpidx(i-1):jumpidx(i)-1)))/dx;
    end
end

%Compute maximal slope magnitude of nearly non-decreasing pieces:

MaxLip=max(max(Lip),1);



%Initialize vector of Lipschitz constants close to shocks  
LipM=  zeros(1,M);



%compute one additional solution for all the rapidly decreasings pieces to compute size of
%boundary shocks
Nfd=N;
xfd=linspace(-4*L,4*L,8*Nfd+1);
ufdo=MaxLip*(xfd(1:end-1)+xfd(2:end))/2;
dxfd=8*L/(8*Nfd);

%Compute extensions of nearly non-decreasing pieces
uM0=zeros(M+1,N);
 for idx=1:M+1    
     if isdecreasing(idx) == 0
        uM0(idx,:)= extend(v(idx,:), idx, xm, jumpidx,jumpsize, Lip,M,N,isdecreasing,deltageneric);     
     end
 end

 %Plot for verification purposes
 %plot(xm, uM0(1,:), 'b', xm, uM0(3,:),'r',xm, uM0(4,:), 'm')

 %nc

 
 
 
 

%Initialize positions of \hat h-curves
hhat=x(jumpidx+1);
hfrakl=x(jumpidx+1);  % Position of \hat h-curve based on neighbors with max index
hfrakr=x(jumpidx+1);  % Position of \hat h-curve based on neighbors with min index

ncc=zeros(1,M);  % number of collisions in computer for each curve.


%Initialze position error of each discontinuity
ipe=dx*ones(1,M);  % We might want to replace this by zero
pe=ipe; %position error at the beginning is the initial position error
pesource=zeros(1,M); %velocity difference caused by uncertainty of which shocks have collided in \hat u but not \hat \psi or vice-versa
abintns=ipe; %iteratively compute integral terms for error estimates instead of completing recomputing each integral at each time step
abintws=ipe; %iteratively compute integral terms for error estimates instead of completing recomputing each integral at each time step

pel=ipe;  %initialize position error for left frak curve
per=ipe;  %initialize position error for right frak curve





%Index of left and right neighbor (\hat v_i) of each discontiuity
lNbr= 1:M;   %left neighbor according to computer
minlNbr=lNbr; %minimal left neighbor
maxlNbr=lNbr; %maximal left neighbor
rNbr= lNbr+1;  %right neighbor according to computer
minrNbr=rNbr;  %minimal right neighbor
maxrNbr=rNbr;  %maximal right neighbor

%ul=2*uM0(:,1)-uM0(:,2); %linearly extrapolate to calculate boundary data for numerical solver
ur=uM0(:,end);

%Initialize "solution in previus time step" for numerical scheme
uMo=uM0; %the first `old value' is just the t=0 value (there is no history yet)
uMn=uMo; % This is just to ensure that a matrix of the correct size is allocated for uMn
t=0;



%Initialize front tracking for rapidly decreasing pieces
% Causes problem if first or last piece is rapidly decreasing (but this is
% not allowed, because we are assuming the initial data is in
% L^\infty(\mathbb R))
num=zeros(1,M+1);
delta=zeros(1,M+1);
for i=1:M+1
    if isdecreasing(i)==1
        dist=x(jumpidx(i))-x(jumpidx(i-1));
        num(i)=ceil(dist/deltageneric);
        snum(i)=num(i);
        delta(i)=dist/num(i);
       % setup positions of front tracking shocks
        ssp(i,1:num(i)-1)=x(jumpidx(i-1)+1)+delta(i)*(1:num(i)-1);
        %Set constant states
        %uft(i,1:num(i))= v(i,jumpidx(i-1)+1);
        
        uft(i,1)= v2(ssp(i,1)-delta(i)/2);
        
        for j=1:num(i)-1            
            uft(i,j+1)= v2(ssp(i,j)+delta(i)/2);
            %This line hard-codes wihich piece is decreasing ( ): )
        end
   
    end
end

%%Gather additional information for rapidly decreasing pieces
%collect information for the construction of the artificial extensions of the front tracking
%shocks using the \hat \Lambda and the \hat P
vim=zeros(1,M+1);vi=zeros(1,M+1);vikm=zeros(1,M+1);vik=zeros(1,M+1);
xim=zeros(1,M+1);xi=zeros(1,M+1);xikm=zeros(1,M+1);xik=zeros(1,M+1);


for i=1:M+1
    if isdecreasing(i)==1
    vim(i)=uft(i,1);
    xim(i)=hhat(i-1)-delta(i)/(4*MaxLip);
    vi(i)=uft(i,1);
    xi(i)=hhat(i-1)+delta(i)+.1;%delta(i)/(4*MaxLip);
    vikm(i)=uft(i,num(i));
    vik(i)=uft(i,num(i));
    xik(i)=hhat(i)+delta(i)/(4*MaxLip);
    xikm(i)=hhat(i)-delta(i)-.1;%delta(i)/(4*MaxLip);
    end
end

    




%Compute maximal time step based on maximal wave speed
dt=.5*dx/(max(max(max(abs(uM0))),max(abs(ufdo))));
dtfd=.5*dxfd/(max(max(max(abs(uM0))),max(abs(ufdo))));

ratio=dt/dtfd;

for k=1:M+1
    if isdecreasing(k)==1
          smin=  min(uft(k,1:num(k)-1) - uft(k,2:num(k)) ); %\min_alpha \frak{s}_\alpha
          smax=  max(uft(k,1:num(k)-1) - uft(k,2:num(k)) ); %\bar{\frak{s}}
    end
end


%Time evolution: time stepping starts
while t < T
    %time increment
   ulo=[3; 1; 0; 0 ];
   t=t+dt;
   uln=[3;1;0; 0];
   % reset pesource (just collecting new things)
   pesource=zeros(1,M);
   
   %Evolve numerical solution in nearly non-decreasing pieces
   
    fl=numflux([ulo,uMo(incl,1:end-1)], uMo(incl,:));
    fr=numflux(uMo(incl,:), [uMo(incl,2:end), ur(incl)]);
    uMn(incl,:)=uMo(incl,:)+dt/dx*(fl-fr);
    
    for k=1:ratio
    %for \hat\Lambda and \hat P
    flfd=numflux([-4*L*MaxLip,ufdo(1:end-1)], ufdo);
    frfd=numflux(ufdo, [ufdo(2:end),4*L*MaxLip]);
    
    ufdn=ufdo+dt/dx*(flfd-frfd);
   

    %Reconstruct \hat Lambda and \hat P
    ufdho(:,:,1)=recon(ufdo,4*N,1,ufdo(1));
    ufdhn(:,:,1)=recon(ufdn,4*N,1,ufdn(1));
    ufdhm=.5*(ufdo+ufdn);
    
    %Compute residuals for ufd (\hat Lambda and \hat P)
    %This is needed for computing error estimators that ensure that the
    %non-constant parts of the extensions of pieces between the front
    %tracking shocks are not revealed, cf. Section 5.1.1
    Rfd=compres(ufdho,ufdhn,dx,dt);
   
    ufdo=ufdn;
    end
    
    %plot(xm,uMn(1,:),'b', xm, uMn(3,:),'r',xm,uMn(4,:),'m')
    %pause
   
 

    
   
    
    %Reconstruct numerical solution
    uhato(:,:,incl)=recon(uMo(incl,:),N,length(incl),ulo);
    uhatn(:,:,incl)=recon(uMn(incl,:),N,length(incl),uln);
    umid=.5*(uhato+uhatn);
    

    
    
    
    

   
    %Compute residuals in nearly non-decreasing pieces
    %(indeed squares of residuals)
    R=zeros(length(incl),N);
    for idx=1:length(incl)
        R(idx,:)=compres(uhato(:,:,incl(idx)),uhatn(:,:,incl(idx)),dx,dt); %idx is index of nearly non-decreasing piece
    end
   %---------------------------------------------------------------------------------
   %Evolve front tracking in rapidly decreasing pieces
   % This does not work if the first or last piece is rapidly decreasing



   %Uft does not need to know about nearly non-decreasing piece
   %evolve shock positions for front tracking shocks
    for i=1:M+1
        if isdecreasing(i)==1
           ssp(i,1:num(i)-1)=ssp(i,1:num(i)-1) + dt*.5*(uft(i,1:num(i)-1)+uft(i,2:num(i)));
       end
    end

    %We write a simple scalar front tracking algorithm: collisions in rapidly decreasing piece
    for i=1:M+1
        if isdecreasing(i)==1
            j=1;
            while j < num(i)-1
                if ssp(i,j) > ssp(i,j+1)
                    s=1;
                    while j+s+1 < num(i) && ssp(i,j) > ssp(i,j+s+1)
                        s=s+1;
                        if j+s+1>num(i)-1
                            disp('all small shocks have collided: what do we do now?')
                            break
                        end
                    end
                    num(i)=num(i)-s;
                    % We take the average to avoid computing the exact time
                    % when the collision happened
                    ssp(i,j)= sum(ssp(i,j:j+s))/(1+s);
                    if j+1 < num(i)
                       ssp(i,j+1:num(i)-1)=ssp(i,j+s+1);
                    end
                      uft(i,j+1:num(i))=uft(i,j+1+s:num(i)+s);                                 
                end
                j=j+1;
            end
        end
    end
 
            
            
   
   %-----------------------------------------------------------------------------
    
    % On which cells might hatv(idx) be revealed, 1 means: might be revealed
    rev=zeros(M+1,N);
    rev(1,:)=(x(1:end-1) < hhat(1) + pe(1));
    for idx=2:M
        rev(idx,:)= (x(2:end) > hhat(idx-1) -pe(idx-1)).* (x(1:end-1)  < hhat(idx) + pe(idx));
        if size(certlist,1)>1
           for k=2:size(certlist,1)
               if idx > certlist(k,1) && idx <= certlist(k,2)
                  rev(idx,:)= zeros(1,N);
               end 
           end
        end
    end
    rev(M+1,:)= (x(2:end) > hhat(M) -pe(M));
    

    
    
    %Increment in L2-Norm^2 of residual
    
    Einkr=dx*dt*sum(sum(R.*rev(incl,:)));

    %for \hat\Lambda
    LhatEinkr=dt*max(Rfd(N:end-N))*L;
    
    %add both together 
    Etot=Etot+Einkr+LhatEinkr;
   
    
    
    LhatE=LhatE+LhatEinkr;
  
    
  %Evolve M \hat h-curves of discontinuity positions (using improved Euler) and
  %their \frak h curves (we are preemptively calculating \frak h curves)
  
  for i=1:M
      if isdecreasing(lNbr(i))==0 && isdecreasing (rNbr(i))==0
          hhat(i)=discevol(hhat,uhato,umid,dt,lNbr,rNbr,i,x);
      elseif isdecreasing(lNbr(i))==1 && isdecreasing (rNbr(i))==0
          hhat(i)=discevolLD(hhat,uhato,umid,ssp,uft,num,dt,lNbr,rNbr,i,x);
      elseif isdecreasing(lNbr(i))==0 && isdecreasing (rNbr(i))==1
          hhat(i)=discevolRD(hhat,uhato,umid,ssp,uft,num,dt,lNbr,rNbr,i,x);
      else
          hhat(i)=discevolDD(hhat,uhato,umid,ssp,uft,num,dt,lNbr,rNbr,i,x);
      end
      
      if isdecreasing(lNbr(i))==0 && isdecreasing (maxrNbr(i))==0
        hfrakl(i)=discevol(hfrakl,uhato,umid,dt,lNbr,maxrNbr,i,x);
      elseif isdecreasing(lNbr(i))==1 && isdecreasing (maxrNbr(i))==0
          hfrakl(i)=discevolLD(hhat,uhato,umid,ssp,uft,num,dt,lNbr,maxrNbr,i,x);
      elseif isdecreasing(lNbr(i))==0 && isdecreasing (maxrNbr(i))==1
          hfrakl(i)=discevolRD(hhat,uhato,umid,ssp,uft,num,dt,lNbr,maxrNbr,i,x);
      else 
          hfrakl(i)=discevolDD(hhat,uhato,umid,ssp,uft,num,dt,lNbr,maxrNbr,i,x);
      end
      
      if isdecreasing(minlNbr(i))==0 && isdecreasing (rNbr(i))==0
        hfrakr(i)=discevol(hfrakr,uhato,umid,dt,minlNbr,rNbr,i,x);
      elseif isdecreasing(minlNbr(i))==1 && isdecreasing (rNbr(i))==0
         hfrakr(i)=discevolLD(hhat,uhato,umid,ssp,uft,num,dt,minlNbr,rNbr,i,x);
      elseif isdecreasing(minlNbr(i))==0 && isdecreasing (rNbr(i))==1
         hfrakr(i)=discevolRD(hhat,uhato,umid,ssp,uft,num,dt,minlNbr,rNbr,i,x); 
      else
         hfrakr(i)=discevolDD(hhat,uhato,umid,ssp,uft,num,dt,minlNbr,rNbr,i,x);      
      end
  end
   
   

%-----------------------------------------------------------------------------------------------
    
   %Account for uncertainty due to h-curves having merged in \psi\hat but the corresponding \hat h curves not merging in computer
   %start with first shock furthest to the left, move to the right
   for i=1:M
       %i number of shock
       
       s= sum( (hfrakr(i) + per(i))>= (hfrakl(i+1:M)-pel(i+1:M)));

            %updating maximal right neighbor and minimal left neighbor
            %accordingly 
        for j=0:s-1
           maxrNbr(i+j)=max(maxrNbr(i+j),i+s+1);
           minlNbr(i+j+1)=min(i,minlNbr(i+j+1));
        end
       
       %is any one of the involved neighbors decreasing 
       oneNbrd=isdecreasing(rNbr(i))+isdecreasing(maxrNbr(i))+isdecreasing(minrNbr(i))+isdecreasing(lNbr(i))+isdecreasing(minlNbr(i))+ isdecreasing(maxlNbr(i));
       if oneNbrd == 0 
       % determine max left state
       
           cell=find(x<hhat(i),1,'last');
           uleftmax=uhato(1,cell,minlNbr(i)) + (uhato(2,cell,minlNbr(i))-uhato(1,cell,minlNbr(i)))/dx*(hhat(i)-x(cell)) ;
      
           cell=find(x<hhat(i),1,'last');
           uleftmin=uhato(1,cell,maxlNbr(i)) + (uhato(2,cell,maxlNbr(i))-uhato(1,cell,maxlNbr(i)))/dx*(hhat(i)-x(cell)) ;
      
       % determine computer left state
           cell=find(x<hhat(i),1,'last');
           uleft=uhato(1,cell,lNbr(i)) + (uhato(2,cell,lNbr(i))-uhato(1,cell,lNbr(i)))/dx*(hhat(i)-x(cell)) ;
       
        
        % determine max right state
           cell=find(x<hhat(i),1,'last');
           urightmax=uhato(1,cell,minrNbr(i)) + (uhato(2,cell,minrNbr(i))-uhato(1,cell,minrNbr(i)))/dx*(hhat(i)-x(cell)) ;
      
       %determine min right state 
           cell=find(x<hhat(i),1,'last');
           urightmin=uhato(1,cell,maxrNbr(i)) + (uhato(2,cell,maxrNbr(i))-uhato(1,cell,maxrNbr(i)))/dx*(hhat(i)-x(cell)) ;
       
       % determine computer right state
       
           cell=find(x<hhat(i),1,'last');
           uright=uhato(1,cell,rNbr(i)) + (uhato(2,cell,rNbr(i))-uhato(1,cell,rNbr(i)))/dx*(hhat(i)-x(cell)) ;
       
       veldiff1 = abs( .5* (uleftmax + urightmax) - .5* (uleft + uright));
       veldiff2 = abs( .5* (uleftmin + urightmin) - .5* (uleft + uright));
       
       
        
       else
           
          % determine computer left state
          %  if isdecreasing(lNbr(i)) == 0
          %     cell=find(x<hhat(i),1,'last');
          %     uleft=uhato(1,cell,lNbr(i)) + (uhato(2,cell,lNbr(i))-uhato(1,cell,lNbr(i)))/dx*(hhat(i)-x(cell)) ;
          %  else
           %  k=sum(  ssp(lNbr(i),1:num(lNbr(i))-1)< hhat(i)    );
          %   uleft=uft(lNbr(i),1+k);
       %     end  
           % determine computer right state
      %     if isdecreasing(rNbr(i)) == 0
       %       cell=find(x<hhat(i),1,'last');
      %        uright=uhato(1,cell,rNbr(i)) + (uhato(2,cell,rNbr(i))-uhato(1,cell,rNbr(i)))/dx*(hhat(i)-x(cell)) ;
       %    else
       %      k=sum(  ssp(rNbr(i),1:num(rNbr(i))-1)< hhat(i)    );
       %      uright=uft(rNbr(i),1+k);
        %   end
          %Speed of sound    
          sos=max(max(abs(uMn)));
          if maxrNbr(i) == rNbr(i) && minlNbr(i) == lNbr(i)
            veldiff1=0; 
            veldiff2=0;
          else
            veldiff1 = abs( sos - .5* (uleft + uright));
            veldiff2 = abs( -sos - .5* (uleft + uright));
          end
       end
      pesource(i)=pesource(i)+max(veldiff1,veldiff2);   
   end

  
    
   
        
    %Check whether `computer' \hat h  discontinuities have merged
  for i=1:M
      s = sum( hhat(i)>= hhat(i+1:M));
      if s > ncc(i)
            ncc(i)=s;
        %Length of time after tn before the merger has happend
        tstar = dt*( hhat0(i) - hhat0(i+s) )/(  hhat(i+s) - hhat0(i+s) - hhat(i) + hhat0(i)) ;
        hmeet=hhat0(i) + (hhat(i)-hhat0(i))*tstar/dt;
        dtloc=dt-tstar;
        
        
        if isdecreasing(rNbr(i))==1 && isdecreasing (lNbr(i+s))==1
           collist=[collist; i, s, lNbr(i) , rNbr(i) , lNbr(i+s), rNbr(i+s), pe(i), pe(i+s), tstar, 100, 100, 100, 1];
        else
           collist=[collist; i, s, hhat(i)-pe(i), hhat(i+s)+pe(i+s), pe(i), pe(i+s), maxlNbr(i), rNbr(i), lNbr(i+s), minrNbr(i+s),pe(i), pe(i+s),0 ];
        end
       for j=0:s
          %neighbors of new curve in the computer
          rNbr(i+j)=rNbr(i+s);
          lNbr(i+j)=lNbr(i);
       end  
              
  

       hloc=discevolscal(hmeet,uhato,umid,dtloc,lNbr(i),rNbr(i),i,x);
       for j=0:s  
         hhat(i+j)=hloc;
       end
       for j=0:s
          %position uncertainty of merged curve is minimum of  position
          %uncertainty of partners in merger
          %disp('merger')
          %pe(i+j)=min(pe(i:i+s))
          %t
       end
       
       
       
      end
      
      
   end
    
   
   %update maximal and minimal neighbors; when can we be sure collisions
   %have happened in \hat\psi ?
  for i=1:M
       
      if size(collist,1)>1
          
          for j=2:size(collist,1)
          if collist(j,13)==0    
              hfrak1= collist(j,3);
            if   isdecreasing(collist(j,7)) == 0 &&  isdecreasing(collist(j,8)) == 0
                cell=find(x<hfrak1,1,'last');
                uleft=uhato(1,cell,collist(j,7)) + (uhato(2,cell,collist(j,7))-uhato(1,cell,collist(j,7)))/dx*(hfrak1-x(cell)) ;
                uright=uhato(1,cell,collist(j,8)) + (uhato(2,cell,collist(j,8))-uhato(1,cell,collist(j,8)))/dx*(hfrak1-x(cell)) ;
                collist(j,3)=hfrak1 + dt*(.5*(uleft+uright) );
                LipL=max(max((uhatn(2,cell,collist(j,7))-uhatn(1,cell,collist(j,7)))/dx),max((uhatn(2,cell,collist(j,8))-uhatn(1,cell,collist(j,8)))/dx)) ;
                js=uleft-uright;
                collist(j,11)=collist(j,11)*exp(LipL*dt);
                collist(j,5) = collist(j,11) + sqrt(1./(4*LipL).*(exp(4*LipL*t)-1)).*(1./sqrt(js)*sqrt(Etot));
            elseif isdecreasing(collist(j,7)) == 1 &&  isdecreasing(collist(j,8)) == 0
                cell=find(hfrak1,1,'last');
                k=sum(  ssp(collist(j,7),1:num(collist(j,7))-1)< hfrak1    );
                uleft=uft(collist(j,7),1+k);
                uright=uhato(1,cell,collist(j,8)) + (uhato(2,cell,collist(j,8))-uhato(1,cell,collist(j,8)))/dx*(hfrak1-x(cell)) ;
                collist(j,3)=hfrak1 + dt*(.5*(uleft+uright) );
                collist(j,11)=collist(j,11)*exp(LipL*dt);
                collist(j,5) = collist(j,11) + sqrt(1./(4*LipL).*(exp(4*LipL*t)-1)).*(1./sqrt(js)*sqrt(Etot));
                collist(j,5) = collist(j,5) + Dinner(collist(j,7)) + (sqrt(2/eps(collist(j,7)))  + t* sqrt(eps(collist(j,7))/2) ) *sqrt(Gamma(collist(j,7)));
                
                 
             elseif isdecreasing(collist(j,7)) == 0 &&  isdecreasing(collist(j,8)) == 1
                 cell=find(x<hfrak1,1,'last');
                 uleft=uhato(1,cell,collist(j,7)) + (uhato(2,cell,collist(j,7))-uhato(1,cell,collist(j,7)))/dx*(hfrak1-x(cell)) ;
                 k=sum(ssp(collist(j,8),1:num(collist(j,8))-1)< hfrak1);
                 uright=uft(collist(j,8),1+k);
                 collist(j,3)=hfrak1 + dt*(.5*(uleft+uright) );
                 collist(j,11)=collist(j,11)*exp(LipL*dt);
                collist(j,5) = collist(j,11) + sqrt(1./(4*LipL).*(exp(4*LipL*t)-1)).*(1./sqrt(js)*sqrt(Etot));
                collist(j,5) = collist(j,5) + Dinner(collist(j,8)) + (sqrt(2/eps(collist(j,8)))  + t* sqrt(eps(collist(j,8))/2) ) *sqrt(Gamma(collist(j,8)));
            else
                k=sum(  ssp(collist(j,7),1:num(collist(j,7))-1)< hfrak1    );
                uleft=uft(collist(j,7),1+k);
                k=sum(ssp(collist(j,8),1:num(collist(j,8))-1)< hfrak1);
                uright=uft(collist(j,8),1+k);
                collist(j,3)=hfrak1 + dt*(.5*(uleft+uright) );
                collist(j,11)=collist(j,11)*exp(LipL*dt);
                collist(j,5) = collist(j,11) + sqrt(1./(4*LipL).*(exp(4*LipL*t)-1)).*(1./sqrt(js)*sqrt(Etot));
                collist(j,5) = collist(j,5) + Dinner(collist(j,8)) + (sqrt(2/eps(collist(j,8)))  + t* sqrt(eps(collist(j,8))/2) ) *sqrt(Gamma(collist(j,8)));
                collist(j,5) = collist(j,5) + Dinner(collist(j,7)) + (sqrt(2/eps(collist(j,7)))  + t* sqrt(eps(collist(j,7))/2) ) *sqrt(Gamma(collist(j,7)));
            end
           
          
            hfrak2= collist(j,4);
            if  isdecreasing(collist(j,9)) == 0 &&  isdecreasing(collist(j,10)) == 0
                cell=find(x<hfrak2,1,'last');
                uleft=uhato(1,cell,collist(j,9)) + (uhato(2,cell,collist(j,9))-uhato(1,cell,collist(j,9)))/dx*(hfrak2-x(cell)) ;
                uright=uhato(1,cell,collist(j,10)) + (uhato(2,cell,collist(j,10))-uhato(1,cell,collist(j,10)))/dx*(hfrak2-x(cell)) ;
                collist(j,4)=hfrak2 + dt*(.5*(uleft+uright) );
                js=uleft-uright;
                LipR=max(max((uhatn(2,cell,collist(j,9))-uhatn(1,cell,collist(j,9)))/dx),max((uhatn(2,cell,collist(j,10))-uhatn(1,cell,collist(j,10)))/dx)); 
                collist(j,12)=collist(j,12)*exp(LipR*dt);
                collist(j,6) = collist(j,12) + sqrt(1./(4*LipR).*(exp(4*LipR*t)-1)).*(1./sqrt(js)*sqrt(Etot));
            elseif isdecreasing(collist(j,9)) == 1 &&  isdecreasing(collist(j,10)) == 0
                cell=find(hfrak2,1,'last');
                k=sum(  ssp(collist(j,9),1:num(collist(j,10))-1)< hfrak2    );
                uleft=uft(collist(j,9),1+k);
                uright=uhato(1,cell,collist(j,10)) + (uhato(2,cell,collist(j,10))-uhato(1,cell,collist(j,10)))/dx*(hfrak2-x(cell)) ;
                collist(j,4)=hfrak2 + dt*(.5*(uleft+uright) );
                
                 
                collist(j,13)=collist(j,12)*exp(LipR*dt);
                collist(j,6) = collist(j,12) + sqrt(1./(4*LipR).*(exp(4*LipR*t)-1)).*(1./sqrt(js)*sqrt(Etot));
                collist(j,6) = collist(j,6) + Dinner(collist(j,9)) + (sqrt(2/eps(collist(j,9)))  + t* sqrt(eps(collist(j,9))/2) ) *sqrt(Gamma(collist(j,9))); 
                 
                 
             elseif isdecreasing(collist(j,9)) == 0 &&  isdecreasing(collist(j,10)) == 1
                 cell=find(x<hfrak2,1,'last');
                 uleft=uhato(1,cell,collist(j,9)) + (uhato(2,cell,collist(j,9))-uhato(1,cell,collist(j,9)))/dx*(hfrak2-x(cell)) ;
                 k=sum(ssp(collist(j,10),1:num(collist(j,10))-1)< hfrak2);
                 uright=uft(collist(j,10),1+k);
                 collist(j,4)=hfrak2 + dt*(.5*(uleft+uright) );
                
                 
                collist(j,13)=collist(j,12)*exp(LipR*dt);
                collist(j,6) = collist(j,12) + sqrt(1./(4*LipR).*(exp(4*LipR*t)-1)).*(1./sqrt(js)*sqrt(Etot));
                collist(j,6) = collist(j,6) + Dinner(collist(j,10)) + (sqrt(2/eps(collist(j,10)))  + t* sqrt(eps(collist(j,10))/2) ) *sqrt(Gamma(collist(j,10))); 
                 
                 
            else
                k=sum(  ssp(collist(j,9),1:num(collist(j,9))-1)< hfrak2    );
                uleft=uft(collist(j,9),1+k);
                k=sum(ssp(collist(j,10),1:num(collist(j,10))-1)< hfrak2);
                uright=uft(collist(j,10),1+k);
                collist(j,4)=hfrak2 + dt*(.5*(uleft+uright) );
                %Update for position error ?? 
                
                
                collist(j,12)=collist(j,12)*exp(LipR*dt);
                collist(j,6) = collist(j,12) + sqrt(1./(4*LipR).*(exp(4*LipR*t)-1)).*(1./sqrt(js)*sqrt(Etot));
                collist(j,6) = collist(j,6) + Dinner(collist(j,10)) + (sqrt(2/eps(collist(j,10)))  + t* sqrt(eps(collist(j,10))/2) ) *sqrt(Gamma(collist(j,10)));
                collist(j,6) = collist(j,6) + Dinner(collist(j,9)) + (sqrt(2/eps(collist(j,9)))  + t* sqrt(eps(collist(j,9))/2) ) *sqrt(Gamma(collist(j,9)));
            end
            
              
           

            
            if collist(j,3)- collist(j,5) > collist(j,4) + collist(j,6)
                minlNbr(collist(j,1):(collist(j,1)+collist(j,2)))=max(collist(j,1), minlNbr(collist(j,1):(collist(j,1)+collist(j,2))));
                maxlNbr(collist(j,1):(collist(j,1)+collist(j,2)))=collist(j,1);
                minrNbr(collist(j,1):(collist(j,1)+collist(j,2)))=collist(j,1)+collist(j,2)+1;
                maxrNbr(collist(j,1):(collist(j,1)+collist(j,2)))=min(collist(j,1)+collist(j,2)+1, maxrNbr(collist(j,1):(collist(j,1)+collist(j,2))));
                certlist=[certlist;collist(j,1), collist(j,2)];
                collist(j,:)=[];
    
            end
          else
               rho=max(uMn(collist(j,3),:)-uMn(collist(j,6)));
               leftdist=rho/2*(t-collist(j,9));
               rightdist=collist(j,7)+collist(j,8) + 2*max(ups(collist(j,4)),ups(collist(j,5)));
               if leftdist > rightdist
                  minlNbr(collist(j,1):(collist(j,1)+collist(j,2)))=max(collist(j,1), minlNbr(collist(j,1):(collist(j,1)+collist(j,2))));
                  maxlNbr(collist(j,1):(collist(j,1)+collist(j,2)))=collist(j,1);
                  minrNbr(collist(j,1):(collist(j,1)+collist(j,2)))=collist(j,1)+collist(j,2)+1;
                  maxrNbr(collist(j,1):(collist(j,1)+collist(j,2)))=min(collist(j,1)+collist(j,2)+1, maxrNbr(collist(j,1):(collist(j,1)+collist(j,2)))); 
                  certlist=[certlist;collist(j,1), collist(j,2)];
                  collist(j,:)=[]; 
               end
                   
          end
          
            
          end
          
          
      end
      
      
      
  end
  
   


 
   
   %Compute current minimal jump sizes and local Lipschitz constants around
   %shocks
   % Should account for boundary shocks as well
 
    for k=1:M    
         %  k number/label of shock
          loc=find( (xm+dx>= hhat(k) - pe(k)).*  (xm-dx<= hhat(k) + pe(k)) );
          LipM(k)=max(max((uhatn(2,loc,lNbr(k))-uhatn(1,loc,lNbr(k)))/dx),max((uhatn(2,loc,rNbr(k))-uhatn(1,loc,rNbr(k)))/dx)) ;
          if isdecreasing(maxlNbr(k))==0 && isdecreasing(minrNbr(k))==0
                jumpsize(k)=min(uhatn(2,loc,maxlNbr(k))-uhatn(2,loc,minrNbr(k)));
          elseif  isdecreasing(maxlNbr(k))==1 && isdecreasing(minrNbr(k))==0
                jumpsize(k)=1;
                for j = loc
                    NNN=maxlNbr(k);
                    xj=xm(j);
                    cell=find(xfd>xj-(xikm(NNN)-vikm(NNN)/MaxLip),1, 'first')-1;
                    Lxj= ufdn(cell);
                    cell=find(xfd>xj-(xik(NNN)-vik(NNN)/MaxLip),1, 'first')-1;
                    Pxj= ufdn(cell);
                    uft(NNN,num(NNN));
                    vln=max(Pxj,min(uft(NNN,num(NNN)),Lxj));   
                    vrn=uhatn(1,j,minrNbr(k));
                    jumpsize(k)= min(jumpsize(k),vln-vrn);
                    if jumpsize(k)<0
                        toc
                        disp('negative jumpsize')
                    %    pause
                    end
                end
          elseif isdecreasing(maxlNbr(k))==0 && isdecreasing(minrNbr(k))==1
                jumpsize(k)=1;
                for j = loc
                    NNN=minrNbr(k);
                    xj=xm(j);
                    cell=find(xfd>xj-(xim(NNN)-vim(NNN)/MaxLip),1, 'first')-1;
                    Lxj= ufdn(cell);
                    cell=find(xfd>xj-(xi(NNN)-vi(NNN)/MaxLip),1, 'first')-1;
                    Pxj= ufdn(cell);
                    
                    vrn=max(Pxj,min(uft(NNN,1),Lxj))  ; 
                    vln=uhatn(2,j,maxlNbr(k));
                    jumpsize(k)= min(jumpsize(k),vln-vrn);
                end
              
              
              %   jumpsize(k)=min(uhatn(2,loc,lNbr(k))-uft(rNbr(k),1));
          else
              % if both neighbors are rapidly decreasing use
              % \mathfrak{s}_min as mininmal jumpsize
                 jumpsize(k)=smin;
          end
    end

    
    
    %Specific position uncertainty for rapidly decreasing pieces
    %prepare for control on boundary shocks
    ups=zeros(1,M+1);
    Gamma=zeros(1,M+1);
    eps=zeros(1,M+1);
    Dinner=zeros(1,M+1);
    for k=1:M+1
        if isdecreasing(k)==1
          %compute the Upsilon
          ups(k)= sqrt(t/smin)*sqrt(2*delta(k)*smax^2+ Etot)*exp(t)+ t*smax*(smax + MaxLip*delta(k));
          ups(k)=ups(k)*(1 + 2*sum(isdecreasing)); % 
          
          Gamma(k)=(max(uft(k,1:num(k)))-min(uft(k,1:num(k)))+ snum(k)*MaxLip*ups(k))*sqrt(t)*sqrt(Etot)*exp(t);
          Gamma(k)=Gamma(k)+snum(k)*MaxLip^2*ups(k)^2;
          Gamma(k)=Gamma(k)+(max(uft(k,1:num(k)))-min(uft(k,1:num(k)))+ snum(k)*MaxLip*ups(k))*2*MaxLip*t*ups(k);
          
          eps(k)=smin/(2*delta(k));
          Dinner(k)=1/sqrt(eps(k))*sqrt(Gamma(k));
          
        end
    end
    
    
    
    
    
    
    
    
    
   
    
    
    
    
    

        
    %Update bound on position error
    
    abintws=(abintws+dt*pesource).*exp(LipM*dt);
    abintns=(abintns).*exp(LipM*dt);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %count=count+1;
    %sumsource(count,:)=sumsource(count-1)+dt*pesource;
    %petrack(count,:)=pe;
    
    
    %[t,pe,pesource]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % The following is needed, otherwise the case LipM=0 (locally) causes trouble
    for k=1:M
        if LipM(k)>1
          prefactor(k)=sqrt(1./(4*LipM(k)).*(exp(4*LipM(k)*t)-1));
        else
          prefactor(k)=  sqrt(3*LipM(k)*t);
        end
    end
            
    pe= abintws  + prefactor.*(sqrt(Etot)./sqrt(jumpsize));
    
    
    %Add additional position error term when \hat h-curve is boundary shock
    for k=1:M
        if isdecreasing(lNbr(k))==1
           pe(k)=  pe(k)+ Dinner(lNbr(k)) + (sqrt(2/eps(lNbr(k)))  + t* sqrt(eps(lNbr(k))/2) ) *sqrt(Gamma(lNbr(k)));
        end
         if isdecreasing(rNbr(k))==1
           pe(k)=  pe(k)+ Dinner(rNbr(k)) + (sqrt(2/eps(rNbr(k)))  + t* sqrt(eps(rNbr(k))/2) ) *sqrt(Gamma(rNbr(k)));
        end
    end
    
    
    %Dinner needs to be added in the two following lines
    pel= abintns  + prefactor.*(1./sqrt(jumpsize)*sqrt(Etot));
    
    
    %Add additional position error term when hat h frak-curve is boundary shock
    for k=1:M
        if isdecreasing(lNbr(k))==1
           pel(k)=  pel(k)+ Dinner(maxlNbr(k)) + (sqrt(2/eps(maxlNbr(k)))  + t* sqrt(eps(maxlNbr(k))/2) ) *sqrt(Gamma(maxlNbr(k)));
        end
         if isdecreasing(rNbr(k))==1
           pel(k)=  pel(k)+ Dinner(maxrNbr(k)) + (sqrt(2/eps(maxrNbr(k)))  + t* sqrt(eps(maxrNbr(k))/2) ) *sqrt(Gamma(maxrNbr(k)));
        end
    end
    
    
    
    per= abintns  + prefactor.*(1./sqrt(jumpsize)*sqrt(Etot));
     for k=1:M
        if isdecreasing(lNbr(k))==1
           per(k)=  per(k)+ Dinner(minlNbr(k)) + (sqrt(2/eps(minlNbr(k)))  + t* sqrt(eps(minlNbr(k))/2) ) *sqrt(Gamma(minlNbr(k)));
        end
         if isdecreasing(rNbr(k))==1
           per(k)=  per(k)+ Dinner(minrNbr(k)) + (sqrt(2/eps(minrNbr(k)))  + t* sqrt(eps(minrNbr(k))/2) ) *sqrt(Gamma(minrNbr(k)));
        end
    end
   
    %What was new is old again
    uMo=uMn;
    hhat0=hhat;
   
  
    
    if t>ts %plot not every time step, but every ts
      f=figure('visible','off');
      
      plotcombined(x,xm,hhat,uMo,M,uhatn,rev,pe);
       saveas(f,sprintf('FIG%03d.jpeg',count))
      ts=ts+T/6;
      count=count+1;
    end
end


%The value computed here is an estimate for L2 error UP TO SHIFTS
L2estup2s=sqrt(Etot)


%Now add the error due to uncertainty on what is revealed
L1est=L2estup2s; %use L^p nesting property on compact domains
i=1;
while i <=M
    loc=find((xm >= hhat(i)-pe(i) ).*(xm <= hhat(i)+pe(i) )); %use bluriness area on shocks
    Linftybound=0;
    for j=loc
       pmax=max(uMn(:,j).*rev(:,j));
       pmin=min(uMn(:,j).*rev(:,j)+1000*(rev(:,j)==0));
       Linftybound=max(pmax-pmin,Linftybound); %L^\infty difference on the region which might be revealed 
    end
    
    L1est=L1est+abs(Linftybound*pe(i));

    find(hhat(i:M) == hhat(i));
    
    %Next line: If shocks have merged in the computer, we only want to
    %count position uncertainty once
    i=i+find(hhat(i:M) == hhat(i),1,'last');
end


LamdaHatInfE=Etot^(1/3)*exp(t)*16*MaxLip;


%----------------------------------------------------------------------------------------
%Compare to "normal" numerical solution on extremely fine mesh, called reference
%solution

%Choose how much the reference solution is refined
Fac=2^5;


%Compute reference solution
Nref=Fac*N;
dxref=L/Nref;

xr=linspace(0,L,Nref+1);
dxr=xr(2)-xr(1);

xmr=.5*(xr(1:end-1)+xr(2:end));

dt=dt*N/Nref;




u0= v1(xmr).*(xmr<x(jumpidx(1)+1)) + v2(xmr).*(xmr >=x(jumpidx(1)+1) ).*(xmr <x(jumpidx(2)+1)) + v3(xmr).*(xmr>=x(jumpidx(2)+1)).*(xmr<x(jumpidx(3)+1)) + v4(xmr)*(xmr >=x(jumpidx(3)+1));
%u0=4*(xmr<-.8*L) +(2+4.5*(xmr+L).^2).*(xmr<-2*L/3).*(xmr>-.8*L)+3*(xmr+L).*(xmr<-L/3).*(xmr>-.8*L)+1*(xmr+L).*(xmr>-L/3).*(xmr<L/6)+1*(xmr>L/6);
uref=u0;
ul=u0(1);
ur=uref(end);
t=0;
while t<T

   t=t+dt;

   %Evolve numerical solution 
    fl=numflux([ul,uref(:,1:end-1)], uref);
    fr=numflux(uref, [uref(:,2:end), ur]);
    uref=uref+dt/dxr*(fl-fr);
  %  plot(uref)
  %  pause(.1)
end


%Project our solution to the reference (extremely fine) mesh
uproj=-100*ones(1,Nref);

loc = find( xmr<=hhat(1));
for j=loc 
    wcell=find(x<xmr(j),1,'last');
    uproj(j)=uhatn(1,wcell,1) + (uhatn(2,wcell,1)-uhatn(1,wcell,1))/dx*(xmr(j)-x(wcell));
end
for i=2:M 
    loc = find( (xmr<=hhat(i) ) .* ( xmr> hhat(i-1) ) );
    for j=loc 
        if isdecreasing(i)==0
        wcell=find(x<xmr(j),1,'last');
        uproj(j)=uhatn(1,wcell,i) + (uhatn(2,wcell,i)-uhatn(1,wcell,i))/dx*(xmr(j)-x(wcell));
        else
         if sum(find(ssp(i,1:num(i)-1)<xmr(j)))==0
            uproj(j)=uft(i,1)  ;
         else
            wcell=find(ssp(i,1:num(i)-1)<xmr(j),1,'last');     
            uproj(j)=uft(i,wcell)   ;
         end
        end
        
    end
end

    loc = find(  xmr > hhat(M)  );
 for j=loc 
        wcell=find(x<xmr(j),1,'last');
        uproj(j)= uhatn(1,wcell,M+1) + (uhatn(2,wcell,M+1)-uhatn(1,wcell,M+1))/dx*(xmr(j)-x(wcell));
 end

 L1est
 L1err=dxr*sum(abs(uref-uproj))
 
 toc
 
 save(['output' num2str(N) 'T' num2str(T) '.mat'])
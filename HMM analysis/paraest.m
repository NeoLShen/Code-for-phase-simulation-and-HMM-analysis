% This script was used for performing 2-state Hidden Malkov Model analysis.
%% Reforming trajectory data structure
ttt=tabulate(Trajsio(:,5)); 
Tralabel=ttt(:,1);
Tralength=ttt(:,2);
IN=Tralength>10 & Tralength<100; %set track length range for analysis
Tralabel(~IN)=[];
Tralength(~IN)=[];
R=cell(1,length(Tralength));
Dis=Trajsio(:,6);
%reforming selected displacement sequence into unit of micrometer
for i=1:length(Tralabel)
    r=Dis(Trajsio(:,5)==Tralabel(i))/1000;
    R{i}=r(2:end);
end
%% Set initial parameters for maximum likehood estimation
D1=0.1;D2=0.01;p12=0.1;p21=0.01; %initial parameters of 2-state model, D1 and D2 are the diffusion coefficients, p12 and p21 are the switching probabilities between 2 states.
%D1=0.1;Sig=0.01;p12=0.1;p21=0.01; %initial parameters of 2-state model, D1
%is the diffusion coefficient of mobile state.
%Sig is the standard deviation of detection error, when treat confined
%state as a pure detection error.
dt=0.03;%time interval between frames in unit of second
theta=[D1,D2,p12,p21]; %form the vector of parameters
LHt=likehood(R,theta,dt,Tralabel,Tralength); %calculate initial likelihood of selected trajectories
% LHt=likehood2(R,theta,dt,Tralabel,Tralength); %If consider the confined state as pure detection error,please use likehood2 function.
MCMC=30000; %total interations of refining parameters
Theta=zeros(MCMC+1,4); %create variable for store parameters' value in each interation
Theta(1,:)=theta;
%% Estimate parameters by maximum likehood estimation
tic
for i=1:MCMC/4
    for k=1:4
        theta=Theta(4*(i-1)+k,:);
        if k>2 %values of switching probabilities should between 0 to 1
            theta(k)=theta(k)+normrnd(0,0.003*min(theta(k),1-theta(k))); %random inreament for parameter
            while theta(k)>=1 || theta(k)<=0
                theta(k)=theta(k)+normrnd(0,0.003*min(theta(k),1-theta(k)));
            end
        else
            theta(k)=theta(k)+normrnd(0,0.003*theta(k)); %random inreament for parameter
        end
        %If consider the confined state as pure detection error,please use likehood2 function.
        %calculate current likelihood of seleted trajectories
        LHp=likehood(R,theta,dt);
        %determing whether accept this change on parameters
        if sum(LHp)>sum(LHt)
            LHt=LHp;
            Theta(4*(i-1)+k+1,:)=theta;
        else
            u=rand;
            if log(u)<(sum(LHp)-sum(LHt))
                LHt=LHp;
                Theta(4*(i-1)+k+1,:)=theta;
            else
                Theta(4*(i-1)+k+1,:)=Theta(4*(i-1)+k,:);
            end
        end
    end
    %reporting current iteration number for every 1000 interations
    if mod(i,1000)==0
        i
    end
end
toc
thetat=mean(Theta(MCMC-10000:MCMC,:)) %average value of parameters in last 10000 iterations
figure,plot(1:MCMC+1,Theta(:,1),1:MCMC+1,Theta(:,2)) %plot the changes of first two values along interations
figure,plot(1:MCMC+1,Theta(:,3),1:MCMC+1,Theta(:,4)) %plot the changes of last two values along interations
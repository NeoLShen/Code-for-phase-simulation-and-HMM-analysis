%% Determine boundary
CreateMask
%% Trajectory assign in condensed phase
Xm=Xm(MaskIN);Ym=Ym(MaskIN);Frame=Frame(MaskIN); %Extract condensed phase localizations first
d0=500; %using default maximum search range 500nm for first round trajectory assign
trajasin %assign trajectories in the condensed phase
% Calculate displacement and mean squared displacement for the first round trajectory assign.
maxt=max(trajs(:,5));
num=length(trajs(:,5));
MSD=zeros(num,1);
Xm=trajs(:,1);Ym=trajs(:,2);Tag=trajs(:,5);
temp=Tag(1);Dis=zeros(num,1);tag=1;
xt=Xm(1);yt=Ym(1);
for i=2:num
    if Tag(i)~=Tag(i-1)
        temp=Tag(i);
        tag=i;
        Dis(i)=0;MSD(i)=0;
        xt=Xm(i);yt=Ym(i);
    else
        Dis(i)=sqrt((Xm(i)-Xm(i-1))^2+(Ym(i)-Ym(i-1))^2);
        MSD(i)=(Xm(i)-xt)^2+(Ym(i)-yt)^2;
    end
end
Trajsio=[trajs,Dis,MSD];
save([name,'trajectoryinout-d500.mat'],'trajs','Dis','Trajsio','MSD')
%% Parameter estimation in the condensed phase with 2-state HMM
paraest
thetat=mean(Theta(20000:30000,:));
save('2state-est-in-d500.mat','Theta','Tralabel','Tralength','R','thetat');
figure,plot(Theta(:,1))
% Determing the optimized maximum search range
din=2.5*sqrt(4*thetat(1)*0.03)*1000;
p12=thetat(3);p21=thetat(4);
p1=p21/(p12+p21);
dout=2.5*sqrt(4*thetat(1)*p1*ratio*0.03)*1000;
%% trajectory assign simultaneously with optimized search range
load([name,'.mat'])
load([name,'-mask.mat'])
trajasinout
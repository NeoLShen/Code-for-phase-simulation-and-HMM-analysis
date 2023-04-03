% This script was used for generating trajectories from localization information recorded in each frame using minimum distance method.
% Trajectory merge and split will be eliminated and each localization will be only assigned to one trajectory.
% 2-state HMM analysis was performed followed by trajectory assgin.
maxf=max(Frame);
minf=min(Frame);
Dm=[];
num=length(Xm);
in0=(Frame==minf);
Tralabel=[1:1:sum(in0)]';
maxt=sum(in0);
% Assign trajectories.
trapre=Tralabel;
for i=minf:(maxf-1)
    inpre=(Frame==i);
    inpost=(Frame==(i+1));
    if (sum(inpost)&sum(inpre))==0
        if sum(inpost)>0
            lpost=sum(inpost);
            trapost=[maxt+1:1:maxt+lpost]';
            maxt=maxt+lpost;
            Tralabel=[Tralabel;trapost];
            trapre=trapost;
        end
    else
        Cpre=[Xm(inpre),Ym(inpre)];
        Cpost=[Xm(inpost),Ym(inpost)];
        Mpre=MaskIN(inpre);
        Mpost=MaskIN(inpost);
        [D,I]=pdist2(Cpre,Cpost,'euclidean','Smallest',1);
        INtemp=Mpre(I);
        dtemp=INtemp*din+~INtemp*dout;
        D1=(D<=dtemp');
        I1=I.*D1;
        It=I(D1);
        Nt=tabulate(It);
        I1t=I1*0;
        while sum(D1)&&sum(Nt(:,2)>1) %split
            St=Nt(:,2)>1;
            nt=Nt(:,1); %pre nodes list
            ntpre=nt(St); %pre-split nodes
            [D2,I2]=pdist2(Cpost,Cpre,'euclidean','Smallest',1);
            ntpost=I2(ntpre(1));
            I1t=(I1==ntpre(1));
            I1t(ntpost)=0;
            I1t=~I1t;%set other post split node as new node
            I1=I1.*I1t;%set other post split node as new node
            Iti=(I1>0);
            It=I1(Iti);
            Nt=tabulate(It);
        end
        lpost=length(I1);
        trapost=[];
        for j=1:lpost
            if I1(j)==0 %new trajectory
                maxt=maxt+1;
                trapost=[trapost;maxt];
            else %assign to old trajectory
                trapost=[trapost;trapre(I1(j))];
            end
        end
        Tralabel=[Tralabel;trapost];
        trapre=trapost;
        Dm=[Dm,D(D1)];
    end
end
traj=[Xm,Ym,Frame,MaskIN];
traj=sortrows(traj,3);
traj=[traj,Tralabel];
trajs=sortrows(traj,5);
Tralabel=trajs(:,5);
% Calculate displacement and mean squared displacement for the first round trajectory assign.
maxt=max(trajs(:,5));
num=length(trajs(:,1));
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
save([name,'trajectoryinout-final.mat'],'trajs','Dis','Trajsio','MSD')

%% Parameter estimation in the condensed phase with 2-state HMM
ttt=tabulate(Trajsio(:,5)); %get track length of each trajectory
time=ttt(:,2);
timein=(time>10)&(time<100); %set track length range for analysis
Tralabel=ttt(:,1);
Tralength=ttt(:,2);
tranum=Tralabel(timein);
li=length(tranum);
R=cell(1);
j=1;
Tralabelin=[];Tralengthin=[];
Dis=Trajsio(:,6);
% Reforming selected displacement sequence into unit of micrometer in condensed phase.
for i=1:li
    int=(Trajsio(:,5)==tranum(i));
    inin=Trajsio(:,4);
    inin=inin(int);
    if sum(~inin)==0 % condensed phase
        r=Dis(Trajsio(:,5)==tranum(i))/1000;
        R{j}=r(2:end);
        j=j+1;
        Tralabelin=[Tralabelin;tranum(i,1)];
        Tralengthin=[Tralengthin;sum(int)];
    end
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
%% Estimate motion parameters in condensed phase by maximum likehood estimation
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
thetat=mean(Theta(20000:30000,:));
p12=thetat(3);p21=thetat(4);
p1=p21/(p12+p21);
save([name,'-2state-est-in-final.mat'],'Theta','Tralabelin','Tralengthin','R','thetat');
%% Diffusion calculate in diluted phase
T=tabulate(trajs(:,5));
time=T(:,2);
timein=(time>5);
tranum=T(:,1);
tranum=tranum(timein);
li=length(tranum);
MSD10in=[];
MSD10out=[];
Endframe=[];
for i=1:li
    int=(Trajsio(:,5)==tranum(i));
    inin=Trajsio(:,4);
    inin=inin(int);
    if sum(inin)==0 % out phase
        msd=Trajsio(:,7);
        msd=msd(int);
        frame=Trajsio(:,3);
        frame=frame(int);
        Endframe=[Endframe;max(frame)];
        MSD10out=[MSD10out;msd(1:6)'];
    end
end
xtf=0.03:0.03:0.15;
MSD10out=MSD10out/10^6;
mean10out=mean(MSD10out);
mean10out=mean10out(2:6);
[a b]=size(MSD10out);
ste10out=std(MSD10out)/sqrt(a);
ste10out=ste10out(2:6);
msd6=MSD10out(:,6);
IN6=msd6>0.012;
MSD10outt=MSD10out(IN6,:);
mean10outt=mean(MSD10outt);
mean10outt=mean10outt(2:6);
ste10outt=std(MSD10outt)/sqrt(a);
ste10outt=ste10outt(2:6);
save([name,'MSD10out-final.mat'],'Endframe','MSD10out','mean10out','ste10out','MSD10outt','mean10outt','ste10outt')
% This script was used for generating simulated molecule motion followed
% 2-state motion switch model with mannually set phase mask in a simulated region with periodic boundaries.
% Generate [X,Y,P,M,t] data, contains x/y coordinates, labels indicate
% whether the molecule in condensed phase (P, 1 for in condensed phase and
% 0 for dilute phase)/motion state (M, 1 for mobile state and 0 for
% confined state)/time stemp (t).
% Parallel computation was required, if you don't want to use parallel
% computation, please replace 'parfor' with 'for'.
%% parameter setting
name=''; %name for storing simulated data, for example '1'
rng(seedn) %set random seed sequence
Din=0.2;Dout=100;
% Pm=1;%mobile ratio
fps=1000000;dt=1/fps;%simulation time step
Pcm30=0.01; %30ms switching probability from confined to mobile state
Pm=0.98; %mobile ratio
Pmc30=(1/Pm-1)*Pcm30; %theoretical swtiching probability from mobile to confined, calculated by mobile <--> immobile balance
Pmc=(Pmc30/(Pmc30+Pcm30))*(1-(1-Pmc30-Pcm30)^(dt/0.03)); %mobile to confined state probability
Pcm=(Pcm30/(Pmc30+Pcm30))*(1-(1-Pmc30-Pcm30)^(dt/0.03)); %confined to mobile state probability
EF=Dout/(Din*Pm);%enrichment fold
stepin=sqrt(2*Din*dt*100);stepout=sqrt(2*Dout*dt);%step size along x/y axis
n=50000;%number of molecules
T=100;%simulation time
%% Phase boundary setting
Xmax=40;Ymax=40; %simulated area: -Xmax to Xmax in x axis and -Ymax to Ymax in y axis
radius=[8,8,8,8,8]; %radius of condensed phase
XC=[-20,-20,0,20,20];YC=[20,-20,0,20,-20]; %center of the condensed phase
Sp=sum(pi*radius.^2); %total area of condensed phase region
Sd=4*Xmax*Ymax-Sp; %total area of dilute phase region
Nd=round(n*Sd/(EF*Sp+Sd)); %number of molecules in dilute phase reigon
Np=n-Nd; %number of molecules in condensed phase reigon
%% Initial position of N molecules
Ptag=ones(n,1);Ptag(randperm(n,Nd))=0;
Mtag=zeros(n,1);Xt=zeros(n,1);Yt=Xt;
for i=1:n
    if Ptag(i)>0 %condensed phase
        Mtag(i)=rand<Pm; %mobile state set
        xt=2*Xmax*rand-Xmax;
        yt=2*Ymax*rand-Ymax;
        while isinphase(xt,yt)==0 %check whether this molecule in the condensed phase, noted that you need to change the settings in this function once you changed the "Phase boundary setting"
            xt=2*Xmax*rand-Xmax;
            yt=2*Ymax*rand-Ymax;
        end
    else %diluted phase
        Mtag(i)=1;
        xt=2*Xmax*rand-Xmax;
        yt=2*Ymax*rand-Ymax;
        while isinphase(xt,yt)==1 %check whether this molecule in the condensed phase, noted that you need to change the settings in this function once you changed the "Phase boundary setting"
            xt=2*Xmax*rand-Xmax;
            yt=2*Ymax*rand-Ymax;
        end
    end
    Xt(i)=xt;Yt(i)=yt;
end
Data=[Xt,Yt,Ptag,Mtag];%[x,t,ptag,mtag]
%% Simulate T seconds
Datat=zeros(n,4,T+1);
Datat(:,:,1)=Data;
tic
for i=2:1+T %each second
    Xtt=Datat(:,1,i-1);Ytt=Datat(:,2,i-1);
    Ptt=Datat(:,3,i-1);Mtt=Datat(:,4,i-1);
%     Xtt=Xt;Ytt=Yt;Ptt=Pt;Mtt=Mt;
    parfor j=1:n
        xt=Xtt(j);yt=Ytt(j);pt=Ptt(j);mt=Mtt(j);
        t=0;
        while t<1
            if pt==1 %in condensed phase
                t=t+dt*100;
                if mt==1 %mobile
                    xt=xt+randn*stepin;yt=yt+randn*stepin;
                    % Periodic simulation boundary, if molecule exceed the
                    % simulation region, set it to periodic position.
                    if xt>Xmax
                        xt=xt-2*Xmax;
                    elseif xt<-Xmax
                        xt=xt+2*Xmax;
                    end
                    if yt>Ymax
                        yt=yt-2*Ymax;
                    elseif yt<-Ymax
                        yt=yt+2*Ymax;
                    end
                    pt=isinphase(xt,yt); %check whether this molecule in the condensed phase, noted that you need to change the settings in this function once you changed the "Phase boundary setting"
                    if pt==1 %still in phase
                        if rand<100*Pmc %switch
                            mt=0;
                        else
                            mt=1;
                        end
                    end
                else %confine
                    if rand<100*Pcm %switch
                        mt=1;
                    else
                        mt=0;
                    end
                end
            else %in diluted phase
                xt=xt+randn*stepout;yt=yt+randn*stepout;
                t=t+dt;
                % Periodic simulation boundary, if molecule exceed the
                % simulation region, set it to periodic position.
                if xt>Xmax
                    xt=xt-2*Xmax;
                elseif xt<-Xmax
                    xt=xt+2*Xmax;
                end
                if yt>Ymax
                    yt=yt-2*Ymax;
                elseif yt<-Ymax
                    yt=yt+2*Ymax;
                end
                pt=isinphase(xt,yt); %check whether this molecule in the condensed phase, noted that you need to change the settings in this function once you changed the "Phase boundary setting"
                if pt==1 %still in phase
                    if rand<Pmc %switch
                        mt=0;
                    else
                        mt=1;
                    end
                end
            end
        end
        Ptt(j)=pt;Mtt(j)=mt;
        Xtt(j)=xt;Ytt(j)=yt;
    end
    Datat(:,1,i)=Xtt;Datat(:,2,i)=Ytt;
    Datat(:,3,i)=Ptt;Datat(:,4,i)=Mtt;
    i
end
toc
save([name,'-',num2str(seedn),'.mat'],'Datat','Din','Dout','Pm','Pmc','Pcm','T','Xmax','Ymax','radius','XC','n','fps')
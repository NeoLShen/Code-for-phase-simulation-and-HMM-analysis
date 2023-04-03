% This script was used for generating simulated molecule motion followed
% 2-state motion switch model with experimental condensed phase mask in a simulated region with periodic boundaries.
% Generate [X,Y,P,M,t] data, contains x/y coordinates, labels indicate
% whether the molecule in condensed phase (P, 1 for in condensed phase and
% 0 for dilute phase)/motion state (M, 1 for mobile state and 0 for
% confined state)/time stemp (t).
%% Parameters setting
rng(seedn)
Areain=sum(sum(Mask)); %calculate total area of condensed phase region
Areaall=1500*3200; %total area of simulated region
Areaout=Areaall-Areain; %calculate the area of dilute phase region 
n=50000; % n is number of molecules
T=2000; % simulation time (Frames)
% Din=0.167; % diffusion coefficient in condensed phase
Din=0.17; % diffusion coefficient in condensed phase
Dout=0.47; % diffusion coefficient in diluted phase
Pm=0.044; % mobile fraction in condensed phase
Tml=0.5; % life time of mobile molecules in condensed phase
fps=10000; % recording fps
dt=1/fps; % simulation time step
Es=Tml/dt; % expectation of mobile steps
%Pmc30=0.03/Tml;
Pmc30=0.828; %30ms mobile to confine switch probability
Pcm30=0.038; %30ms confine to mobile switch probability
Pmc=(Pmc30/(Pmc30+Pcm30))*(1-(1-Pmc30-Pcm30)^(dt/0.03)); %mobile to confined state probability
Pcm=(Pcm30/(Pmc30+Pcm30))*(1-(1-Pmc30-Pcm30)^(dt/0.03)); %confined to mobile state probability
% Pmi=0;Pim=1; %no motion switch setting for all mobile conditions
Pm=Pcm/(Pcm+Pmc); %theoretical mobile ratio, calculated by mobile <--> immobile balance
P=P; % phase boundaries
Ef=Dout/(Din*Pm); % enrichment folds
ratio=Ef*Areain/Areaout;
stepin=1000*sqrt(2*Din*dt);stepout=1000*sqrt(2*Dout*dt); %average step size
%% Generate initial data set
Data=zeros(n,4); % [X,Y,P,M,0]
% Generate molecules with uniform distribution
for i=1:n
    xt=32000*rand;yt=15000*rand;
    px=ceil(xt/10);py=ceil(yt/10);
    ptag=Mask(py,px);
    if ptag==1
        if rand>Pm
            mtag=0;%immobile
        else
            mtag=1;%mobile
        end
    else
        mtag=1;
    end
    Data(i,:)=[xt,yt,ptag,mtag];
end
% Generate molecules by initial enrichment fold, please delete the '%' in
% this section and add '%' in previous section.
% for i=1:n
%     if rand<1/(ratio+1) % generate molecules by initial enrichment fold
%         ptag=0;%diluted phase
%     else
%         ptag=1;%condensed phase
%     end
%     xt=32000*rand;yt=15000*rand;
%     px=ceil(xt/10);py=ceil(yt/10);
%     while Mask(py,px)~=ptag
%         xt=32000*rand;yt=15000*rand;
%         px=ceil(xt/10);py=ceil(yt/10);
%     end
%     if ptag==1
%         if rand>Pm
%             mtag=0;%immobile
%         else
%             mtag=1;%mobile
%         end
%     else
%         mtag=1;
%     end
%     Data(i,:)=[xt,yt,ptag,mtag];
% end
%% Simulate for T seconds
Datat=zeros(n,4,T+1);
Datat(:,:,1)=Data;
tic
for j=1:n
    for i=2:1+T
        xt=Datat(j,1,i-1);yt=Datat(j,2,i-1);
        pt=Datat(j,3,i-1);mt=Datat(j,4,i-1);
        for k=1:0.03/dt
            if mt==0 %immobile
                if rand<Pcm %switch to mobile
                    mt=1;
                end
            else %mobile
                if pt==1 %in condensed phase
                    xt=xt+randn*stepin;yt=yt+randn*stepin;
                    % Periodic simulation boundary, if molecule exceed the
                    % simulation region, set it to periodic position.
                    if xt>32000
                        xt=xt-32000;
                    else
                        if xt<0
                            xt=xt+32000;
                        end
                    end
                    if yt>15000
                        yt=yt-15000;
                    else
                        if yt<0
                            yt=15000+yt;
                        end
                    end
                    px=ceil(xt/10);py=ceil(yt/10);
                    pt=Mask(py,px);
                    if pt==1
                        if rand<Pmc %switch to immobile
                            mt=0;
                        else 
                            mt=1;
                        end
                    end
                else %in diluted phase
                    xt=xt+randn*stepout;yt=yt+randn*stepout;
                    % Periodic simulation boundary, if molecule exceed the
                    % simulation region, set it to periodic position.
                    if xt>32000
                        xt=xt-32000;
                    else
                        if xt<0
                            xt=xt+32000;
                        end
                    end
                    if yt>15000
                        yt=yt-15000;
                    else
                        if yt<0
                            yt=15000+yt;
                        end
                    end
                    px=ceil(xt/10);py=ceil(yt/10);
                    pt=Mask(py,px);
                    if pt==1
                        if rand<Pmc %switch to immobile
                            mt=0;
                        else 
                            mt=1;
                        end
                    else
                        mt=1;
                    end
                end
            end
        end
        Datat(j,1,i)=xt;Datat(j,2,i)=yt;
        Datat(j,3,i)=pt;Datat(j,4,i)=mt;
    end
    if mod(j,100)==0
        j/100
    end
end
toc
save(['Datat-linear-2-',num2str(seedn),'.mat'],'Datat','Din','Dout','Pm','Pmc30','Pcm30','fps','P','Mask','Pcm','Pmc')
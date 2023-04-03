% This script was used for generating trajectories from localization information recorded in each frame using minimum distance method.
% Trajectory merge and split will be eliminated and each localization will be only assigned to one trajectory.
maxf=max(Frame);
minf=min(Frame);
Dm=[];
num=length(Xm);
in0=(Frame==minf);
Tralabel=[1:1:sum(in0)]';
maxt=sum(in0);
trapre=Tralabel;
% Assign trajectories.
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
        [D,I]=pdist2(Cpre,Cpost,'euclidean','Smallest',1);
        D1=(D<=d0);
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
%Generate new variable that contains trajectory labels and sort by trajectory labels.
traj=[Xm,Ym,Frame,MaskIN(MaskIN)];
traj=sortrows(traj,3);
traj=[traj,Tralabel];
trajs=sortrows(traj,5);
Tralabel=trajs(:,5);
%This script was used for creating mask for condensed phase based on local density of the localization got in the 2D single molecule tracking experiment.
%The imported data should at least contain three variables which are named "Xm", "Ym", "Frame".
%"Xm" and "Ym" are the vectors of x and y coordinates of the localizations with the unit of nanometer.
%"Frame" is a vector of every localization's frame number respectively.
%% Initial settings
name=''; %name of saved files, for example '1'
pixelsize=10; %pixel size (nm) to calculate the local density
extendsize=20; %extended pixels for morphological operation
DensityThreshold=1; %local density that higher than DensityThreshold times average local density will recognized as the condensed phase
mr=strel('disk',9); %morphological operation structures
mr2=strel('disk',5); 
mr3=strel('disk',2);
lx=ceil(max(Xm)/pixelsize)+20;ly=ceil(max(Ym)/pixelsize)+20;
%calculate local density of each pixel
I=zeros(ly,lx);
num=length(Xm);
for i=1:num
    x=Xm(i);y=Ym(i);
    xp=ceil(x/10);yp=ceil(y/10);
    I(yp,xp)=I(yp,xp)+1;
end
%% Morphological operation for creating condensed phase mask
I2=conv2(I,mr.Neighborhood,'same');
xpm=ceil(max(max(Xm))/10);xpn=ceil(min(min(Xm))/10);
ypm=ceil(max(max(Ym))/10);ypn=ceil(min(min(Ym))/10);
area=(xpm-xpn+1)*(ypm-ypn+1); %calculate total area
numt=sum(sum(I2))/area; %calculate average local density
I3=I2>DensityThreshold*numt; %create binary mask based on a density threshold
I4=imclose(I3,mr); %smooth the boundary
I4=imopen(I4,mr);
I4=imclose(I4,mr2);
I4=imopen(I4,mr2);
I4=imclose(I4,mr3);
I4=imopen(I4,mr3);
Mask=I4; %final mask for condensed phase
P=bwboundaries(I4,8); %get boundaries for each isolated region
% calculate total area of condensed phase region
L=bwlabel(I4,8);
Area=regionprops(L,'Area');
areain=0;
la=length(Area);
for i=1:la
    areain=areain+Area(i).Area;
end
%% Create label tags that whether the localization in the condensed phase (as 1) or not (as 0)
lp=length(P);
MaskIN=zeros(num,1);MaskIN=logical(MaskIN);
numin=0;
%generate demonstration figure for the phase boundary
figure,plot(Xm,Ym,'k.','markersize',4),axis equal
Ratio=[];
hold on
for i=1:lp
    bond=P{i};
    xt=bond(:,2)*10-5;yt=bond(:,1)*10-5;
    xt=[xt;xt(1)];yt=[yt;yt(1)];
    IN=inpolygon(Xm,Ym,xt,yt);
    MaskIN=MaskIN|IN;
    numin=numin+sum(IN);
    Area(i).num=sum(IN);
    Ratio(i)=(Area(i).num/Area(i).Area)/(numout/(area-areain));
    plot(xt,yt,'b-')
end
hold off
numout=num-numin;
EnrichmentFold=(numin/areain)/(numout/(area-areain)); %calculate the enrichment fold between condensed phase region and dilute phase region
%% Save the information of condensed phase mask and the demonstration figure
save([name,'-mask.mat'],'Mask','P','MaskIN','num','numin','numout','areain','area','EnrichmentFold')
saveas(1,[name,'-mask.fig'])
function LH=likehood(R,theta,dt,Tralabel,Tralength)
D1=theta(1);D2=theta(2);p12=theta(3);p21=theta(4);
p1=p21/(p12+p21);p11=1-p12;
p2=1-p1;p22=1-p21;
nt=length(Tralabel);
LH=zeros(nt,1);
for n=1:nt %calculate likehood for every trajectory
    r=R{n};
    lt=Tralength(n)-1;
    a(1,1)=log(p1)-log(D1*dt)-r(1)^2/(4*D1*dt);
    a(1,2)=log(p2)-log(D2*dt)-r(1)^2/(4*D2*dt);
    for j=2:lt
        a(j,1)=a(j-1,1)+log(p11+p21*exp(a(j-1,2)-a(j-1,1)))-log(D1*dt)-r(j)^2/(4*D1*dt);
        a(j,2)=a(j-1,1)+log(p12+p22*exp(a(j-1,2)-a(j-1,1)))-log(D2*dt)-r(j)^2/(4*D2*dt);
    end
    LH(n)=a(lt,1)+log(1+exp(a(lt,2)-a(lt,1)));
end
end
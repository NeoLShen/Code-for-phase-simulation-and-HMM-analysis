function id=isinphase(xt,yt)
id=0;
x=abs(xt);y=abs(yt);
if x^2+y^2<=64
    id=1;
elseif (x-20)^2+(y-20)^2<=64
    id=1;
end
end
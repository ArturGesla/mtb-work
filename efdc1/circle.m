function circle(x,y,r,str,fa)
phi=linspace(0,2*pi,100); 
plot(x+r*cos(phi),y+r*sin(phi),"Color",str);
a=char(str);
a1=a(2:3); a2=a(4:5); a3=a(6:7); 
cc=[hex2dec([a1])/256 hex2dec([a2])/256 hex2dec([a3])/256];
fill(x+r*cos(phi),y+r*sin(phi),cc,"facealpha",fa,"edgealpha",0);
end
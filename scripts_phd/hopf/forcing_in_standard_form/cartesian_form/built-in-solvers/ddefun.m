
function dydt = ddefun(t,y,Z,mu,gamma,om,k,beta) % equation being solved
  ylag1 = Z(:,1);
  ylag2 = Z(:,2);

  dydt = [ylag1(1); 
          ylag1(1)+ylag2(2); 
          y(2)];
  dydt=[y(1)*(mu+y(1)^2+y(2)^2-(y(1)^2+y(2)^2)^2)-y(2)*(om+gamma*(y(1)^2+y(2)^2))-k*cos(beta)*(y(1)-ylag1(1))+k*sin(beta)*(y(2)-ylag1(2)); 
      y(2)*(mu+y(1)^2+y(2)^2-(y(1)^2+y(2)^2)^2)+y(1)*(om+gamma*(y(1)^2+y(2)^2))-k*sin(beta)*(y(1)-ylag1(1))-k*cos(beta)*(y(2)-ylag1(2))];
end
%-------------------------------------------

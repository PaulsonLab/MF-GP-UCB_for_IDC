%%function [X] = Building_dynamics(xk,uk,wk)
% Simulates building response
function [X] = Building_dynamics(xk,uk,wk)

%Building dynamics function
C1 = 9.356*10^5; C2 = 2.970*10^6; C3 = 6.695*10^5;
K1 = 16.48; K2 = 108.5; K3 = 5; K4 = 30.5; K5 = 23.04;
eta_c = 2;
eta_h = 4;

%xk = [t1,t2,t3]; uk = [h,c]; wk=[Tod, DNI, IO]
%xk = sdpvar(3,1); uk = sdpvar(2,1); wk=sdpvar(3,1);
dx1 = 3600*(1/C1)*((K1+K2)*(xk(2)-xk(1)) + K3*(wk(1)-xk(1))+ K5*(xk(3)-xk(1))+ ...
            wk(2)+ uk(1)*eta_h- uk(2)*eta_c+ wk(3));
dx2 = 3600*(1/C2)*((K1+K2)*(xk(1)-xk(2)) + wk(2)  ); 
dx3 = 3600*(1/C3)*(K5*(xk(1)-xk(3))+ K4*(wk(1)-xk(3))) ;
X = [dx1;dx2;dx3] + xk;

end



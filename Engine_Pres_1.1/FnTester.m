% function [A_Orf, h_max] = OrficeArea(A_leak, Orf_diam, theta, h)
% h_max=(Orf_diam*(1-sqrt(1-cos(theta/2))))/sin(theta);
% 
% if (h<=0)
%     A_Orf= A_leak;
% elseif (h>0) && (h<h_max)
%     A_Orf=(pi*(Orf_diam-(h*sin(theta/2)*cos(theta/2)))*h*sin(theta/2))+A_leak;
% elseif (h>h_max)
%     A_Orf=(pi*(Orf_diam)^2/4)+A_leak;
% else
%     A_Orf=(pi*(Orf_diam)^2/4)+A_leak;
% end
% %A_Orf=(pi*(Orf_diam)^2/4)+A_leak;

% 
% 
% %
% 
% 
% 
% % function Cv = Cvcalc(u)
% % if (u<6)
% %     Cv = u/20;
% % elseif(u>6&&u<9)
% % Cv=((1.1169*(u-5.95)^2))/10+0.3;
% % else
% %   Cv=0;
% %     
% 
% function A_inj  = FnTester(m_dot_des)
% Cd=0.50;
% MW=44;
% P=3.62e6;
% Pe=2.413e6;
% Vhat_l=0.05;
% 
% A_inj=(m_dot_des)/(Cd*sqrt(abs(2/MW*(P - Pe)/Vhat_l))*MW);
% disp(A_inj);



function Z =FnTester(P,T)
Z=py.CoolProp.CoolProp.PropsSI('Z', 'P', P, 'T', T, gas);

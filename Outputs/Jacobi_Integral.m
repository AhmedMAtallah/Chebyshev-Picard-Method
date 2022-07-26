%================={ Copyright (c) Ahmad Bani-Younes 2013 }=================
%                Texas A&M University - AEROSPACE department
% 
% File name     : Jacobi_Integral.m
% Subject       : Computes the Hamiltonian.
% Description   : To validate orbit propagators accuracy.
%    
% Sub-files     : egm2008GPsphericalharmonic.m
%
% Compiler      : MATLAB 7.11.0 (R2010b)
% Date Modified : 04/10/2013
%==========================================================================
function Jacobi_Integral(time, soln, method)
global omega Deg;
r_B = soln(:,1:3);       % position in the rotating frame
v_B = soln(:,4:6);       % velocity in the rotating frame
term  = 0.5*omega*omega.*(r_B(:,1).^2 + r_B(:,2).^2);
KE    = 0.5*( v_B(:,1).^2 + v_B(:,2).^2 + v_B(:,3).^2 );
V      = -egm2008GPsphericalharmonic( r_B, Deg);        % Use for w_vec = [0 0 w]
% V = -gravitysphericalharmonic(r_B, Deg);   % Use for w_vec = [w1 w2 w3]
Enrgy  = V + KE - term ;
Enrgyo = Enrgy(1) ; 
dEnrgy = abs((Enrgy - Enrgyo)/Enrgyo);

b = time(end);
figure('color',[1 1 1]);
semilogy(time/b,dEnrgy,'r','LineWidth',2);grid;
set(gca,'fontsize',12)
set(findobj('type','axes'),'fontsize',14)
%title(['Hamiltonian: (', method, ')'])
xlabel('Integration Time (t/Tp)')
ylabel('| H - H_o | / | H_o |')









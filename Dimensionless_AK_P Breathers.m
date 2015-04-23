close all
clear all

t_cc=linspace(-150,150,500);
% x_cc=linspace(-50,50,2000);

x_plot = 160;

A = 0.45;
B = sqrt(8*A*(1-2*A));
a0_cc = 0.02;
T0 = 1; %s
g=9.81;
h = 1;

w = 2*pi/T0;

% k0, from w^2=g*k0*tanh(k0*h)
    syms k0_vpa   
    k0 = -vpasolve(k0_vpa*tanh(k0_vpa*h) == (w^2)/g, k0_vpa);
    k0 = double(k0) % convert sym to number
    % k0=w^2/g;   % k0=w^2/g is only valid for high h(deep water), difference shows at h<0.5

   
c_g=w/(2*k0)

% x_cc = c_g*t_cc; % so X=0
x_cc=linspace(0,80,200);
% x_cc=0;

alpha = -w/(8*k0^2);
beta = -w*k0^2/2;
w_mod=2*sqrt(1-2*A);

% Xprime = sqrt(beta/alpha)*a0_cc*(x_cc-c_g.*t_cc);
%  Xprime = 0;
Tprime = c_g^2*sqrt(beta/alpha)^2*a0_cc^2*alpha*t_cc;

% primes  
      j=1;
       for xi_cc=x_cc
            Xprime = c_g^2*sqrt(beta/alpha)*a0_cc*(xi_cc-c_g.*t_cc);
            
            q2(j,:) = (1+((2*(1-2*A)*cosh(2*B*Tprime)+1i*B*sinh(2*B*Tprime))...
            /(sqrt(2*A)*cos(w_mod*Xprime)-cosh(2*B*Tprime))))...
            *a0_cc.*exp(2*1i*Tprime);
        
            a2(j,:) = q2(j,:)./(sqrt(2)*k0^2); % rescaling no longer needed
            eta2(j,:)= -q2(j,:).*exp(1i*(k0*xi_cc-w*t_cc));
            
            j=j+1;
       end
  
%  dimensionless typical domain size
    %  T = -alpha*t_cc;
     X = linspace(-10,10,400);
     T = linspace(-10,10,800);
      j=1;
      for xi_cc=X
          if (A==0.5)
              C=4*(1+4*1i*T);
              D=1+4*xi_cc.^2+16*T.^2;
              q(j,:) = (1-(C./D)).*exp(2*1i*T);
              
              j=j+1;
          else
              C = (2*(1-2*A)*cosh(2*B*T)+1i*B*sinh(2*B*T));
              D = (sqrt(2*A)*cos(w_mod*xi_cc)-cosh(2*B*T));
              
              q(j,:) = (1+(C./D)).*exp(2*1i*T);
              a(j,:) = q(j,:)./(sqrt(2)*k0^2); % rescaling no longer needed
              
              j=j+1;
          end
      end

 % Dimensionless domain solver
     X_sol= linspace(-110,49,400);
     T_sol = linspace(-2,0.45,200);
      j=1;
       for xi_cc=X_sol   
           
           if (A==0.5)
               C=4*(1+4*1i*T_sol);
               D=1+4*xi_cc.^2+16*T_sol.^2;
               q_sol(j,:) = (1-(C./D)).*exp(2*1i*T_sol);
               
               j=j+1;
           else
               C = (2*(1-2*A)*cosh(2*B*T_sol)+1i*B*sinh(2*B*T_sol));
               D = (sqrt(2*A)*cos(w_mod*xi_cc)-cosh(2*B*T_sol));
               
               q_sol(j,:) = (1+(C./D)).*exp(2*1i*T_sol);               
               a_sol(j,:) = q_sol(j,:)./(sqrt(2)*k0^2); % rescaling no longer needed
               
               j=j+1;
           end
            
         
        end 
       
       
figure
plot(t_cc,abs(q2(x_plot,:)))
hold on
plot(t_cc,real(eta2(x_plot,:)),'--r')
title('prime coordinates, envelope and carrier wave in time @')

figure
plot(x_cc,abs(q2(:,10)))
hold on
plot(x_cc,real(eta2(:,19)),'--r')
title('prime coordinates, envelope and carrier wave in space @')

figure
imagesc(t_cc,x_cc,abs(q2))
xlabel('t')
ylabel('x')
colorbar
title('prime coordinates, abs, 2d plot')

figure
imagesc(t_cc,x_cc,real(q2))
xlabel('t')
ylabel('x')
colorbar
title('prime coordinates, real, 2d plot')

% Dimensionless, typical domain

        figure
        surf(T,X,abs(q))
        xlabel('v')
        ylabel('u')
        colorbar
        shading interp
        light
        lighting phong
         title('typical domain')
% Dimensionless, solver domain
        figure
        surf(T_sol,X_sol,abs(q_sol))
        xlabel('v')
        ylabel('u')
        colorbar
        shading interp
        light
        lighting phong
        title('domain solver')


        
        

function y=RungeKutta4_WB(Coeff_NL,a,dz,nt,tau,r,a_crit)   
F1 = zeros(1,nt);
F2 = zeros(1,nt);
F3 = zeros(1,nt);
F4 = zeros(1,nt);
aTemp = zeros(1,nt);
a2 = abs(a).^2;


%* Evaluate F1 = f(x,t).
WB = (-1/tau)*((a2./a_crit).^r-1).*heaviside(a2-a_crit);
F1 =  ((Coeff_NL*a2)+WB).*a;

%* Evaluate F2 = f( x+tau*F1/2, t+tau/2 ).
half_dz = 0.5*dz;
aTemp = a + half_dz*F1;
aTemp2 = abs(aTemp).^2;
WB_Temp = (-1/tau)*((aTemp2./a_crit).^r-1).*heaviside(a2-a_crit);

F2 = (Coeff_NL*aTemp2+WB_Temp).*aTemp;

%* Evaluate F3 = f( x+tau*F2/2, t+tau/2 ).
aTemp = a + half_dz*F2;
aTemp2 = abs(aTemp).^2;
WB_Temp = (-1/tau)*((aTemp2./a_crit).^r-1).*heaviside(a2-a_crit);

F3 = (Coeff_NL*aTemp2+WB_Temp).*aTemp;

%* Evaluate F4 = f( x+tau*F3, t+tau ).
aTemp = a + dz*F3;
aTemp2 = abs(aTemp).^2;
WB_Temp = (-1/tau)*((aTemp2./a_crit).^r-1).*heaviside(a2-a_crit);

F4 = (Coeff_NL*aTemp2+WB_Temp).*aTemp;

%* Return x(t+tau) computed from fourth-order R-K.
y = a + dz/6.*(F1 + F4 + 2.*(F2+F3));
return 

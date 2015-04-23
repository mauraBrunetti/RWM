function y=RungeKutta4(Coeff_NL,Envelope,dz,nt)   
F1 = zeros(1,nt);
F2 = zeros(1,nt);
F3 = zeros(1,nt);
F4 = zeros(1,nt);
EnvelopeTemp = zeros(1,nt);

%* Evaluate F1 = f(x,t).
F1 =  Coeff_NL*(abs(Envelope).^2).*Envelope;

%* Evaluate F2 = f( x+tau*F1/2, t+tau/2 ).
half_dz = 0.5*dz;
EnvelopeTemp = Envelope + half_dz*F1;
F2 = Coeff_NL*(abs(EnvelopeTemp).^2).*EnvelopeTemp;

%* Evaluate F3 = f( x+tau*F2/2, t+tau/2 ).
EnvelopeTemp = Envelope + half_dz*F2;
F3 = Coeff_NL*(abs(EnvelopeTemp).^2).*EnvelopeTemp;

%* Evaluate F4 = f( x+tau*F3, t+tau ).
EnvelopeTemp = Envelope + dz*F3;
F4 = Coeff_NL*(abs(EnvelopeTemp).^2).*EnvelopeTemp;

%* Return x(t+tau) computed from fourth-order R-K.
y = Envelope + dz/6.*(F1 + F4 + 2.*(F2+F3));
return 

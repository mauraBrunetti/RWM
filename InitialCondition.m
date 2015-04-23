 function InitialCondition(nc,mkFigsInitial);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NLS equation
% Akhmediew breather  
% March 17, 2015
% Maura Brunetti, Hubert Branger, Debbie Eeltink
%
% Generates an akhmediev breather initial condition for the NLSE solver
% First
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%
% nc = 1; % number to give filename
mkFigsInitial = 1; %(1= figures for the initial condition, 0= no figures)
close all 

%% Notes

%%% See Kibler et al. "The Peregrine soliton in nonlinear fibre optics", 
%%%     Nature Physics 6, 790-795 (2010)
%%% Note that: \xi -> 2T in our notation
%%%            in Peregrine there is a sign - 
%%%            \xi corresponds to 2T, \tau to X, beta = gamma, alpha =
%%%            beta_2/2

% maura: 
% a = envelope
% a_init  = surface elevation [note difference in the exponent with
% respect to my code: zetaT = real(y.*exp(-2*pi*1i*f0*t)); ]
%
% The following is used to define an elevation function that is 0 in t(1) 
% and in t(end) 
% Amplitude = real(a.*exp(1i.*(k0.*x-w.*t)));
% Put everywhere exp(-i omega t) (with - sign)  in RogueWaves_Hubert
% in order to have an increase of the sub-sideband
% (otherwise with opposite sign I obtain an increase of the super-sideband) !

%% Initialization
T0 = 1; 		% period of the wave carrier [s]
f0 = 1/T0;      % Peak of the spectrum [Hz]

g = 9.81;       % gravitation constant m/s^2
dt = 1/400;		% time step  (generation in ascii format each time step, 
                % 400 is frequency of wavemaker)
h = 0.85;       % water depth [m?] (Hubert: 0.85, Maura:100. 
ak = 0.07;	    % wave steepness, epsilon in Brunetti et al 2013

A = 0.25;	    % Akhmediev breather parameter. when aa decrease the period of the envelope decrease
                % 0.25 to have max instability (0.5 for a peregrine
                % breather) % aa = 0.4
x_f = -9.;      % distance from the focal point (x<0)
                % distance from the focal point (x<0) %x=-25;		
      
boundary = 0;   % zero is standard (?)
Nbramp= 0;  %20 (Hubert value);	 % for a delayed start

%Define type of inital condition
    % akhm = 0; no akhmedieve breather, but JONSWAB spectrum
    % akhm = 1; akhmediev breather
    % akhm = 2; this generated akhmediev input signal
    akhm = 2; 

% Define time domain 
    ntot = 10; % number of total periods/ wavelengths
    nbefore = 5; % number of periods before the soliton

% define space domain
   CalculateSpace = 1; % 1 = yes, 0 = no
   l_wavetank = 40; % m, length wavetank
   x=linspace(x_f,l_wavetank-(-x_f), 1200);
   %     x=linspace(0,140, 1200);
    

%% Calculate parameters

% time domain
    Ntot = 2 * 10 * ntot; % gives 10 periods, where do the factors come from?
    Nstart= 2 * 10 * nbefore; 
    Nend=Ntot-Nstart;  
    %
    prec = 1; %if I increase prec, spectrum in f is more accurate
    nf = 40000*prec; 

    rangeT = Nstart*T0+Nend*T0;
    rangeT = prec*rangeT;
    t = s_space(-Nstart*T0*prec,Nend*T0*prec,nf);  %evenly spaced sampled values from xstart to (xend-dx)

    dt = t(2)-t(1); % timestep
    df = 1/dt;      % frequency step

% carrier frequency omega
    w=2*pi/T0;  

% k0, from w^2=g*k0*tanh(k0*h)
    syms k0_vpa   
    k0 = -vpasolve(k0_vpa*tanh(k0_vpa*h) == (w^2)/g, k0_vpa);
    k0 = double(k0); % convert sym to number
    % k0=w^2/g;   % k0=w^2/g is only valid for high h(deep water), difference shows at h<0.5

% dispersion parameter
    q=k0*h;     

% group velocity: can be calculated in two ways, use c_g:
    c_g1=T0*g*(k0*h+sinh(k0*h)*cosh(k0*h))/4/pi/cosh(k0*h)^2;
    c_g1

    sigma=tanh(h.*k0);
    c_g= w/(2*k0*sigma)*(sigma + k0*h*(1.-sigma^2));
    c_g

%% Coefficents Nonlinear and Dispersion terms 

% alpha = dispersion coefficient in NLS eq (2) 
    % q = h*k0
    % alpha = -1/2*(d^2omega/dk^2) is only valid in deep water, st hk -> inf 
    % below is the general form:
    alpha = 0.5.*(- (2.*g.*h.*(tanh(q).^2 - 1) - 2.*g.*h.^2.*k0.*...
    tanh(q).*(tanh(q).^2 - 1))/(2.*(g.*k0.*tanh(q)).^(1/2)) -...
    (g.*tanh(q) - g.*h.*k0.*(tanh(q).^2 - 1)).^2/(4.*(g.*k0.*tanh(q)).^(3/2)));

% beta = nonlinear coefficient in NLS eq (2) ? plus minus can be constructed in two ways
    % beta1 gives the correct answer only for large h, such that kh-> inf
    % we do not use beta1, only beta
    %1.
    beta1 = -(w.*k0.^2)/(16.*(sinh(q)).^4).*(cosh(4.*q)+8-2.*(tanh(q)).^2)+...
    ((w/(2.*(sinh(2.*q)).^2)).*(2.*w.*(cosh(q)).^2+k0.*c_g).^2/(g.*h-c_g.^2));
    beta1 

    %2. See notes maura
    UU = 9-12*sigma^2 + 13*sigma^4 - 2*sigma^6;
    VV = 2 + c_g*(1-sigma^2)*k0/w;
    WW = (2*sigma^3)*VV/(sigma*c_g^2*k0^2/w^2-h*k0);
    beta = -w*g^2*(UU+VV*WW)/(16*sigma^2*w^4/k0^4);
    beta

    
%% Coordinate transformation to convert solution of the space NLSE to time NLSE

% additional normalization
    alpha = alpha/c_g^3; 
    beta = beta/c_g;
% go to the frame of ?
    t = t+x_f/c_g; % focus point/phase velocity
    a0 = ak/k0; % starting amplitude

% go to dimensionless NLSE (3)
    T=sqrt(beta/(2*alpha))*a0*(t-(1/c_g)*x_f);  
    X=sqrt(beta/(2*alpha))^2*a0^2*alpha*x_f;   

%% Calculate Initial Condition

% aa=0.5 	: Peregrine breather, the period in time and space is infinite
% aa<0.5	: Akhmediev-Peregrine, periodic in time
% solution of the amplitude envelope a is normalized by the initial
% envelope height a0

if A<0.5
    B=sqrt(8*A*(1-2*A)); % B in NLSE (3)
    w_mod=2*sqrt(1-2*A); % omega_mod in NLSE (3)
        
    a=-a0*(1+(2*(1-2*A)*cosh(2*B*X)+1i*B*sinh(2*B*X))./(sqrt(2*A)*cos(w_mod*T)-cosh(2*B*X))).*exp(2*1i*X);
else
    a= a0*(-1+(4*(1+4*1i*X))./(1+4*T.^2+16*X.^2)).*exp(2*1i*X); 
end;

% surface elevation initial condition, eqn (1): 
    eta = real(a.*exp(-1i.*(w.*t))); 
    % Hubert doest kx-wt, maura only wt, because FT is only on the time 
    % domain, so it should not matter.

% calculate anaylical solution for 2d domain:
if (CalculateSpace ==1)
    
   % use a coarser time vector for 2d picture
   t_coarse_init= linspace(min(t),max(t),400);
    
    j = 1;
    jj=1;
    a2d = zeros(size(x,2),size(t_coarse_init,2));
    if A<0.5
        B=sqrt(8*A*(1-2*A)); % B in NLSE (3)
        w_mod=2*sqrt(1-2*A); % omega_mod in NLSE (3)

        for (x_i = x)

        % go to dimensionless NLSE (3)
        t_coarse = t_coarse_init+(x_i./c_g);
        T_coarse = sqrt(beta/(2*alpha))*a0*(t_coarse-(1/c_g)*x_i);  
        X_coarse = sqrt(beta/(2*alpha))^2*a0^2*alpha*x_i;  
            
            
        a_coarse=-a0*(1+(2*(1-2*A)*cosh(2*B*X_coarse)+1i*B*sinh(2*B*X_coarse))./(sqrt(2*A)*cos(w_mod*T_coarse)-cosh(2*B*X_coarse))).*exp(2*1i*X_coarse);
        a2d (j,:) = a_coarse;
        j = j+1;
        if (ismember(j,[2;20;40;60;80;100;120;140;160;180;200]))         
            record(:,jj) = [x_i T min(T) max(T)];
            jj=jj+1;
    
        end
        end
    else
        for (x_i = x)
        t_coarse = t-(x_i./c_g);       
        T_coarse = sqrt(beta/(2*alpha))*a0*(t_coarse-(1/c_g)*x_i);  
        X_coarse = sqrt(beta/(2*alpha))^2*a0^2*alpha*x_i;  
        a_coarse= a0*(-1+(4*(1+4*1i*X_coarse))./(1+4*T_coarse.^2+16*X_coarse.^2)).*exp(2*1i*X_coarse);
        a2d (j,:) = a_coarse;
        j = j+1;
        end
    end;
    

end
    figure
    imagesc(t,x,abs(a2d))
    xlabel('t')
    ylabel('x')
    title('envelope absolute');
    colorbar;
    
    figure
    surf(t_coarse_init,x,abs(a2d))
    xlabel('t-x/c_g (s)')
    ylabel('x (s)')
    light
    shading interp
    lighting phong
    
          figure
    surf(t_coarse_init,x,abs(a2d)/a0)
    xlabel('t-x/c_g (s)')
    ylabel('x (s)')
    light
    shading interp
    lighting phong
    
    
    
% change orientation of vector and change imaginary part to negative?
    a = a'; 
%% make solution periodic on the time domain: trim edges
        a(1)
        a(end) 

        figure
        plot((repmat(real(a),[1 2])));
        
% trim to have minama on both sides
     min_a = min(real(a));
     tolerance = abs(max(real(a))-min_a)*1e-5; %depends on 'height' of a 
     %for T=0.1 a is in the range of 1e-5, for T=1 in the range of 1e-2 
     a_llim = find (abs(real(a)-min_a) < tolerance , 1, 'first'); %index fist minumum on the left side
     a_rlim =  find(abs(real(a)-min_a)< tolerance, 1, 'last'); %index fist minumum on the right side
     
     % clip a, and clip n
     a = a(a_llim:a_rlim);
     
     
     % clip a, eta and X
     eta = eta(a_llim:a_rlim);
     t = t(a_llim:a_rlim);
     T = T(a_llim:a_rlim);
     
% If this results in an odd vector, make it even by taking of 1 entry

      if(mod(size(a,1),2)==1)   % even number=0, odd number=1
        a = a(1:end-1);
        eta = eta(1:end-1);
        t = t(1:end-1);
        T = T(1:end-1);
      end
      
      
      figure
     plot((repmat(real(a),[2 1])));
     title('repeated matrix A');
%% Optional: boundary nonzero and/or Nramp nonzero

% if there is a nonzero boundary, shift things so the boundary matches.
if (boundary ~= 0) 
  if eta(1)>0
    t0 = find(eta<0,1,'first');
    t1 = find(eta<0,1,'last');
  else
    t0 = find(eta>0,1,'first');
    t1 = find(eta>0,1,'last');
  end
  %t = t; % can delete this?
  if t1==length(t)
    t(t1+1)=t(t1);
  end
  if (mod(length(eta(t0:t1)),2)==0) %%even
    a2=[0 eta(t0:t1) 0]';
    y2=[0 a(t0:t1) 0]';
    t2=t(t0-1:t1+1) ;
    t=t2-t2(1);
  else
    a2=[0 eta(t0:t1) 0 0]';
    y2=[0 a(t0:t1) 0 0]';
    t2=t(t0-1:t1+2) ;
    t=t2-t2(1);  
  end    
  eta=a2;
  a = y2;
end

% signal brut eta sans fonction de transfert
if (Nbramp ~= 0) 
    tramp=Nbramp*T0;
    ramp=ones(1,length(t));
    ramp(1:round(tramp*df))=sin((0:(round(tramp*df)-1))*pi/tramp/2/df);
    ramp(end-round(tramp*df)+1:end)=sin((((round(tramp*df)-1):-1:0))*pi/tramp/2/df);
    eta=(eta-mean(eta(1:round(tramp*df)))).*ramp';
end
    
%% Fourier transform: spectrum of surface elevation and envelope

nff = length(t);
f = s_space(0,df,nff);  %%fa = 400
nFreq = length(f);

%Spectrum Surface elevation 
eta_FT = FourierT(eta,dt); 

%Spectrum of the envelope/Amplitude
a_FT = FourierT(a, dt); %

%Fourier transform back: a_si should be the same as a
a_si = IFourierT(a_FT, dt);

%% Figures

if (mkFigsInitial == 1)
    
% Envelope/Amplitude modulation
    figure;
    hold on 
    plot(t,real(a),'b')
    plot(t,real(a_si),'r--');  %%to check fft
    plot(t,imag(a),'g--')
    xlabel('t (s)');
    ylabel('surface height (m)');
    legend('real(a)','real(a_si)','imag(a)')
    grid on
    title('Initial: Envelope (a) real time')

% Surface elevation and envelope
    figure
    hold on 
    plot(t,eta,'b','linewidth',1)
    plot(t,real(a),'m--','linewidth',2)
    plot(t,imag(a),'r','linewidth',2)
    plot(t,abs(a),'c','linewidth',2)
    xlabel('t')
    ylabel('surface height (m)');
    legend('eta','real(a)','imag(a)','abs(a)')
    title('Initial: surface elevation (eta) and envelope (a)')
    %axis([0 15 -6*10^(-3) 6*10^(-3)]) 
    
% Fourier Transform Surface Elevation 
    figure;
    plot(f(1:1000),abs(eta_FT(1:1000)/df),'r','linewidth',1);
    xlabel('f (Hz)');
    ylabel('fft eta');
    legend('eta_FT')
    title('Initial: Surface elevation Spectrum (eta_{FT})');
    grid on;axis([1.2 2.8 0 2*10^(-4)]);
    
% Fourier Transform of Envelope
    figure
    plot(f(1:50),abs(a_FT(1:50)/df),'linewidth',1)
    title('Initial: Spectrum envelope (a_{FT})')
    legend('a_{FT}')
    xlabel('f (Hz)')
    grid on
end

%% Significant wave height

% Both are not used for further calculation..

% First definition: Hs_sd0 = 4*sigma (sigma = standard deviation)
% this is the right definition
MeanH0 = mean(eta);
Hs_sd0 = 4*sqrt(mean((eta-MeanH0).^2)); % not used (?)
display(['Significant wave height Hs (from standard deviation): ' num2str(Hs_sd0) ' m']); 

%Second definition: integral of the surface-elevation spectrum
% sigma = \int S(f) df 
% Hs_spec0 = 4*sqrt(trapz(f,abs(eta_FT)/(df))), There is a factor of 2.92 
% of difference... sqrt(2*pi) = 2.5066, so this comes closer but is still not correct 
Hs_spec0 = 4*sqrt(trapz(f,abs(eta_FT)/(2*pi*df))); %%%???? normalization????
display(['Significant wave height Hs (from spectrum): ' num2str(Hs_spec0) ' m']); 

%% Save Data

% turn some vectors for processing in next script
    eta_FT = eta_FT.';
    a_FT = a_FT';
    a = a';

%save all relevant parameters in a matrix, clear the others
    save(sprintf('spectrumHubert%s.mat',num2str(nc)),'nc','a','t','f','dt','f0','nf','df','a_FT','Hs_spec0')
    save(sprintf('spectrumHubert%s.mat',num2str(nc)),'akhm','eta_FT','c_g','ak','g','x_f','-append')
    clear all
return
%exit


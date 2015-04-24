%function RogueWaves(nc) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NLS equation
% JONSWAP spectrum  (? not used, or is the option still in the code?
% October 24, 2013
% Maura Brunetti, Nadege Marchiando, Nicolas Berti
%
% test based on the script by Les Schmerr, Iowa Univ.
% in google: s_space.m matlab script schmerr
% see info in : http://www.public.iastate.edu/~e_m.350/FFT%205.pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
%% Set Parameters 

%1 = TRUE, 0 = FALSE

nonlinear = 1;      %?
RK = 4;             % level of ranga kutta
mkFigs = 1;         % make figures yes/no
mkFigsLive = 0;     % make the live figure of the wave propagation
mkFigsInitial = 1;  % make the figures for the initial condition

nc = 1; % filename for the initial condition matrix

FactorWO = 1;   %optional amplification factor weak wind (Wind Onorato) model
FactorWM = 1;   %optional amplification factor strong wind (Wind Maura) model

WaveBreaking = 0; % switch wavebreaking on/off

MoveWithGroup = 1; % Move with group velocity in t

windWOBlows = 1; % switch wind on/off for weak wind (Wind Onorato) model
windWMBlows = 1; % switch wind on/off for strong wind (Wind Maura) model

%% Load Initial Condtion File
tic % start timer
InitialCondition(nc,mkFigsInitial) % This will give the initial condition
% and its parameters, namely:

% ak: 
% g: gravitation constant

% fo: frequency carrier wave.
% t: time vector
% dt: timestep
% f: frequency vector
% df: frequency setep
% a: complex amplitude / envelope vector on t
% a_FT : Spectrum of the envelope on f
% eta_FT: Spectrum of the surface elevation on f
% Hs_spec0: Characteristic waveheight based on spectrum

load(sprintf('spectrumHubert%s.mat',num2str(nc)))
% only use this RogueWave solver for the Marseille Akhmediev breather initial
% condition: akhm = 2 
if (akhm ~= 2) 
    disp('Only valid for the Hubert case') 
    return
end 

%% Set Constants
nu = 1.05e-06; %Kinematic viscosity 
prec = 1e9; %//Check to 9 decimal points in mod function

Dist = 40.; %%20000 for the old simulations. propagation distance in meters
dz = 0.02;    %%1. distance step in calculations
outZ = 0.1;   %10  distance step for images

%wavebraking, see Trulsen and Dysthe 'Frequency downshift through self
%modulation'
tau = 1/8;
r = 2;
a_crit = 0.000035/ak; %k0 or k?

%% Initialization Hs, f*, viscosity
% Define characteristic maximum wave height
    HsMAX = 2.*Hs_spec0; %%1.2 for old sims

% Define f*: without f0 the pulse moves ? whats happening here?
    nFreq = length(f);
    f_mid = round(nFreq/2);
    fstar = f; 
    jj = 0;
    for ii=nFreq:-1:nFreq/2+1
        jj=jj+1;
        fstar(ii) = -f(jj);
    end   

% Viscosity
    visc = 4*nu*(2*pi)^5*f0^5/g^3; % where is this term in the equation?
    visc = 0.; %%%mau
    %only if I multiply by 1 million I start to observe the effect of
    %viscosity!!
    %visc = visc*1000000;

%% Define coefficients for Split Step Fourier
% initialize linear coefficients
    Coeff_L = zeros(1,nFreq);       % linear coefficient
    Coeff_LV = zeros(1,nFreq);      % linear coefficient plus viscoscity
    Coeff_LwindM = zeros(1,nFreq);  % linear coefficients in Strong wind model (maura)
  % Coeff_LwindO is a constant so no preallocation needed

    CL = zeros(1,nFreq);            % All linear coefficients NO wind
    CL_WO = zeros(1,nFreq);         % All linear coefficients Weak/Onorato wind
    CL_WM = zeros(1,nFreq);         % All linear coefficients Strong/ Maura wind
    
% linear coefficient no visc or wind
    Coeff_L(:) = 1i*4*pi^2*(fstar(:)-f0).^2/g;

% linear coefficient plus viscosity term
    %Coeff_LV(:) = Coeff_L(:) - visc; 
    Coeff_LV(1:nFreq/2) = Coeff_L(1:nFreq/2) - visc; 
    Coeff_LV(nFreq/2+1:nFreq) = Coeff_L(nFreq/2+1:nFreq) + visc; 

% linear coefficients wind term
    % GammaM (Miles growth rate) in eqn (4)
    % Wind as in Proment & Onorato 2012
    GammaMO = (ak)^2*f0; %%old sims: 0.0004*2*pi*f0 which gives 0.005 Hz = (0.05)^2*2
    Coeff_LwindO = 4*GammaMO*pi*f0/g;

    % GammaM, Wind as in paper Maura (strong wind)
    GammaMM = f0*ak/14.; %Hz %%old sims: 0.01;
    TermDeriv = zeros(1,nFreq);
    TermDeriv(1:nFreq/2) =  -8*pi*GammaMM*(fstar(1:nFreq/2)-f0)/g;
    TermDeriv(nFreq/2+1:nFreq) = 8*pi*GammaMM*(fstar(nFreq/2+1:nFreq)+f0)/g;
    Coeff_LwindM(1:nFreq) = - 2*1i*GammaMM^2/g + TermDeriv(1:nFreq);

%Nonlinear coefficient 
  %  Coeff_NL = -1i*(2*pi*f0)^6/g^3;
    Coeff_NL = -1i*(2*pi*f0)^6/g^3;

    
%% Initialize space steps and result vectors
% time/ space steps
    nZ = floor(Dist/outZ); % number of bigger z steps for total distance
    if (nZ==0) 
        nZ = 1;
    end    
    nt = length(t); % number of time steps

% Pre-allocate: vector Env2D is the envelope in t and z
    Env2D_ana=zeros(nZ,nt);         % no wind, test compare analytical
    Env2D=zeros(nZ,nt);         % no wind
    Env2D_WO=zeros(nZ,nt);      % weak wind (Wind Onorato)
    Env2D_WM=zeros(nZ,nt);      % strong wind (Wind Maura)
    
    if (MoveWithGroup==1)
        Env2D_cg=zeros(nZ,nt); 
    end
  
% Pre-allocate: Norm and Momentum Vector
    Momentum=zeros(nZ,1);        
    Norm=zeros(nZ,1);         
    
% Pre-allocate: vector Eta is the surface elevation 
    Eta=zeros(1,nZ);

% Pre-allocate: vectors for the the characteristic wave height Hs
    Hs_sd = zeros(1,nZ);
    Hs_sdWO = zeros(1,nZ);
    Hs_sdWM = zeros(1,nZ);

% a_FT is fourier transform of a(x,t), it is the spectrum of the envelope
    a_FT0 = a_FT; 
    a_FT_WO=a_FT;
    a_FT_WM=a_FT;
 
% Space steps
    Z=0.; % start position

    iz=1;  % first step bigger lengthscale
    izWO=1;
    izWM=1;

    Z0 = 0:outZ:Dist; % outer/bigger lengthscale steps (used for xx and xx plots)
    Zi = 0:dz:Dist; % inner/smaller caluclational step in z: dz
    tolerance = 0.001; % to determine if small step is close to big step

% live plot of surface elevation and envelope
    if (mkFigsLive == 1 && mkFigs==1)
       handle_LivePlotEta = figure;
       handle_LivePlotEnvelope = figure;
    end
    
% test for SSFT method
    a_new=a;
    
%% Split step Fourier Transform
for Z = Zi  
%% Split step calculation Envelope     
    %display(['Distance: ' num2str(Z) ' m']);
     CL = Coeff_LV; % linear coeffient without wind, but with with viscosity
     CL_WO = Coeff_L + FactorWO*Coeff_LwindO;  % linear with weak wind, viscosity inculeded in gamma factor Onorato.
     CL_WM = Coeff_LV + FactorWM*Coeff_LwindM;  % linear with strong wind
 
     if (nonlinear == 1)
      % if (STFT ==1 ) % calculate nonlinear part analytically
          
           %tic
           % linear part, encapuslation of NL
            a_new2=abs(a_new).^2;
            Env_NL = exp(-1i*Coeff_NL*dz.*a_new2).*a_new;     
            Env_NL_FT = FourierT(Env_NL,dt);            
            EnvTest = IFourierT(exp(CL*dz).*Env_NL_FT,dt);         
          % toc
        
      % else    %use Ranga Kutta       
           % Timestep linear part
             Envw = a_FT.*exp(CL*dz);
             Envw_WO = a_FT_WO.*exp(CL_WO*dz);
             Envw_WM = a_FT_WM.*exp(CL_WM*dz);

             Env = IFourierT(Envw,dt);  
             Env_WO = IFourierT(Envw_WO,dt); 
             Env_WM = IFourierT(Envw_WM,dt); 
           
           % Timestep nonlinear part 
           if (RK == 2)  
             K1 = Coeff_NL*(abs(Env).^2).*Env;
             K1_WO = Coeff_NL*(abs(Env_WO).^2).*Env_WO;
             K1_WM = Coeff_NL*(abs(Env_WM).^2).*Env_WM;
             %
             Env_temp = Env + K1*dz/2; 
             Env_temp_WO = Env_WO + K1_WO*dz/2;
             Env_temp_WM = Env_WM + K1_WM*dz/2;
             %estimate at the central point
             K1 = Coeff_NL*(abs(Env_temp).^2).*Env_temp; 
             K1_WO = Coeff_NL*(abs(Env_temp_WO).^2).*Env_temp_WO; 
             K1_WM = Coeff_NL*(abs(Env_temp_WM).^2).*Env_temp_WM; 
             %use the central-point estimate to perform the total time-step advancement
             Env = Env + K1*dz; 
             Env_WO = Env_WO + K1_WO*dz; 
             Env_WM = Env_WM + K1_WM*dz;
           elseif (RK == 4)
               if (WaveBreaking==1) % if wavebreaking is on or off
                   Env = RungeKutta4_WB(Coeff_NL,Env,dz,nt,tau,r,a_crit);
                   Env_WO = RungeKutta4_WB(Coeff_NL,Env_WO,dz,nt,tau,r,a_crit);
                   Env_WM = RungeKutta4_WB(Coeff_NL,Env_WM,dz,nt,tau,r,a_crit);
               else
                   Env = RungeKutta4(Coeff_NL,Env,dz,nt);
                   Env_WO = RungeKutta4(Coeff_NL,Env_WO,dz,nt);
                   Env_WM = RungeKutta4(Coeff_NL,Env_WM,dz,nt);             
               end
           end
           
           % Move with group velocity reference frame
           if (MoveWithGroup == 1)
             t_shift=Z/c_g;
             t_nshift = round(t_shift/dt);         
            Env_cg = circshift(Env,[0,t_nshift]);
           end   
          % toc
       %end   
%        close all
%        figure
%        title('Ranga Kutta Envelope v.s. analytical envelope');
%        plot(real(Env))
%        hold on
%        plot(real(EnvTest),'g')
%        legend('RK ','analytical');
%        
%        figure
%        title('Ranga Kutta Envelope v.s. analytical envelope');
%        plot(imag(Env))
%        hold on
%        plot(imag(EnvTest),'g')
%        legend('RK ','analytical');
       
     end
     
     %if (mod(floor(Z*prec), outZ*prec) == 0)
     %floor(Z/outZ*10.), floor(Z/outZ)*10
     
%% Live Plot 

     if (mkFigsLive == 1 && mkFigs==1)
     %Surface elevation    
          figure(handle_LivePlotEta)  % make the liveplot figure window active
          
          %plot surface elevation (see def in Onorato 2001)
          plot(t,real(Env.*exp(-2*1i*pi*f0*t))); %%changed sign
          plot(t,real(Env_WM.*exp(-2*1i*pi*f0*t)),'r');  %%changed sign$
          ylim([-9e-3,9e-3]) 
          xlabel('x [m] ')
          ylabel('surface elevation eta [m] ')
          grid on
          title([num2str(Z) ' m']);
          drawnow;      
          
      %Envelope    
          figure (handle_LivePlotEnvelope)
          
          subplot(1,2,1)
          plot(t,real(Env),'b');  
          title([num2str(Z) ' m']);
          ylim([-9e-3,9e-3]) 
          legend('Runga Kutta')
          drawnow;
          
               
          subplot(1,2,2)
          plot(t,real(EnvTest),'g');  
          ylim([-9e-3,9e-3]) 
          legend('analytical')
          title([num2str(Z) ' m Timeslice Envelope']);
          drawnow;
          pause(0.1)
    
     end  
%      if (Z == 28.5 || Z == 25 ||Z == 15 ||Z == 10.88 || Z == 6.04 || Z ==  1.3 || Z ==  0.5 || Z ==  0.1)
%              figure 
%              hold on
%                 plot(t,real(Env),'b');  
%                 plot(t,real(EnvTest),'g');      
%                ylim([-9e-3,9e-3]) 
%                xlabel('x [m] ')
%                ylabel('etha [m] ')
%                legend('RK ','analytical');
%                grid on
%                title([num2str(Z) ' m']);
%       end
      %set startingpoint for new step.
     a_new = EnvTest;
     a_FT = FourierT(Env,dt);
     a_FT_WO = FourierT(Env_WO,dt);
     a_FT_WM = FourierT(Env_WM,dt);

%% Surface elevation and significant wave height in z-direction
     % Calculate for steps zi (larger steps than z) in Z
   
     if (abs(Z - Z0(iz)) < tolerance)
         
     % Calculate Envelope in z direction for every t 
       Env2D_Test(iz,:) = EnvTest(:); %analytical SSFT
       Env2D(iz,:) = Env(:);
       Env2D_WO(iz,:) = Env_WO(:);
       Env2D_WM(iz,:) = Env_WM(:);
        
       if (MoveWithGroup == 1)
           Env2D_cg(iz,:) = Env_cg(:);
           
       end
       
       
     % Calculate surface elevation from Envelope  
       crossProd = Env2D(iz,:)'.*exp(-2*pi*1i*f0*t(:));  %%changed sign
       crossProdWO = Env2D_WO(iz,:)'.*exp(-2*pi*1i*f0*t(:)); %%changed sign
       crossProdWM = Env2D_WM(iz,:)'.*exp(-2*pi*1i*f0*t(:)); %%changed sign
      
     % Mean wave hight 
       MeanH0 = mean(mean(real(crossProd)));
       MeanH0_WO = mean(mean(real(crossProdWO)));
       MeanH0_WM = mean(mean(real(crossProdWM)));
       
     % Significant wave height Hs = 4x std(surface elevation)   
       Hs_sd(iz) = 4*sqrt(mean(mean((real(crossProd)-MeanH0).^2)));
       Hs_sdWO(iz) = 4*sqrt(mean(mean((real(crossProdWO)-MeanH0_WO).^2)));
       Hs_sdWM(iz) = 4*sqrt(mean(mean((real(crossProdWM)-MeanH0_WM).^2)));
       
     % Calculate Momentum and Norm, NOT YET FUNCTIONAL 
       a_FTshift = circshift(a_FT,[1 f_mid]);
       f_shift = f-f_mid; 
       Momentum(iz) = sum((f_shift.*abs(a_FTshift)),2); % momentum in  m^2 s^-1 
       Norm(iz) =  sum(abs(a_FT),2); % norm in m^2
    
     % Account for wavebreaking, OLD METHOD OF CUTTING AT 90% of max
       % Wave breaks is Amplitude > steepness * A0
%            if (windWOBlows == 1)
% 
%              if (Hs_sdWO(iz) > HsMAX)
%                   izWO = iz; 
%                   disp(['izWO : ' num2str(izWO)]) 
%                   FactorWO = 1; % continue to blow but amplitude is reduced              
%                   windWOBlows = 1;
%                   a_FT_WO = FourierT(Env_WO*0.9,dt);
%              end 
%            end
%            if (windWMBlows == 1)
% 
%              if (Hs_sdWM(iz) > HsMAX) %%% old sims: HsMAX*0.975) 
%                 izWM = iz; 
%                 disp(['izWM : ' num2str(izWM)]) 
%                 FactorWM = 1; % continue to blow but amplitude is reduced
%                 windWMBlows = 1;
%                 a_FT_WM = FourierT(Env_WM*0.9,dt);
%              end
%            end    
%            
      % Step     
        Eta(iz) = Z;
        iz = iz+1;

     end  
end

%% Significant wave height if akhm ~= 2 

%
%Find rogue waves
%[i,j]=find(Eta2D>7.3);
%%%
%

% find significant wave height in case of
    if (akhm ~= 2)
            %%%% SIGNIFICANT WAVE HEIGHT %%%%%%%%%
            %%%MeanH = mean(mean(abs(Eta2D)))
            %%%Hs = 2*MeanH/0.64
            nombre_RW = zeros(1,nZ+1);
            nombre_RW_WO = zeros(1,nZ+1);
            nombre_RW_WM = zeros(1,nZ+1);
            HTOT = [];
            HTOT_WO = [];
            HTOT_WM = [];
            save HsCheck Hs_sd Hs_sdWO Hs_sdWM HsMAX
        for k = izWO:min(izWO+1000,size(Env2D,1))
            HW1=HeightWave(Env2D(k,:).*exp(-2*1i*pi*f0*t),t); %%changed sign
            HTOT = [HTOT HW1];
            HW2=HeightWave(Env2D_WO(k,:).*exp(-2*1i*pi*f0*t),t); %%changed sign
            HTOT_WO = [HTOT_WO HW2];
            ii = 0; jj = 0; 
            [ii,jj] = hist(HW1,100);
            nombre_RW(k) = sum(ii(jj>=2.2*Hs_sd(k)));
            ii = 0; jj = 0; 
            [ii,jj] = hist(HW2,100);
            nombre_RW_WO(k) = sum(ii(jj>=2.2*Hs_sdWO(k)));
        end
        for k = izWM:min(izWM+1000,size(Env2D_WM,1))
        HW3=HeightWave(Env2D_WM(k,:).*exp(-2*1i*pi*f0*t),t); %%changed sign
        HTOT_WM = [HTOT_WM HW3];
        ii = 0; jj = 0; 
        [ii,jj] = hist(HW3,100);
        nombre_RW_WM(k) = sum(ii(jj>=2.2*Hs_sdWM(k)));    
        end   
        %%SAVE DATA
        save(sprintf('CommonOutput%s.mat',num2str(nc)),'t','f','Zeta','gamma','alpha','f0','Dist','nonlinear','dz','HsMAX','outZ')
        %save CommonOutput t f Zeta gamma alpha f0 Dist nonlinear dz HsMAX outZ
        %
        save(sprintf('Eta2dim%s.mat',num2str(nc)),'Eta2D','HTOT','Hs_sd','nombre_RW')
        %save Eta2dim Eta2D HTOT Hs_sd nombre_RW
        %
        save(sprintf('Eta2dimWO%s.mat',num2str(nc)),'Eta2D_WO','GammaMO','izWO','HTOT_WO','Hs_sdWO','nombre_RW_WO')
        %save Eta2dimWO Eta2D_WO GammaMO izWO HTOT_WO Hs_sdWO nombre_RW_WO
        %
        save(sprintf('Eta2dimWM%s.mat',num2str(nc)),'Eta2D_WM','GammaMM','izWM','HTOT_WM','Hs_sdWM','nombre_RW_WM')
        %save Eta2dimWM Eta2D_WM GammaMM izWM HTOT_WM Hs_sdWM nombre_RW_WM
        %
    end
    
   %
    if (akhm == 0) 
       psi = jonswap_spectrum(alpha,gamma,f0,f); 
       Hs_spec0 = 4*sqrt(trapz(f,psi)*0.5); % Significant wave height (integral of surface-elevation spectrum)
       display(['Hs from the initial spectrum : ' num2str(Hs_spec0) ' m']);
    end

%% PostProcessing: Figures and Information

% display significant wave height:
    display(['Significant wave height Hs: ' num2str(Hs_sd(1)) ' m']);
    display(['Significant wave height Hs WO: ' num2str(Hs_sdWO(1)) ' m']);
    display(['Significant wave height Hs WM: ' num2str(Hs_sdWM(1)) ' m']);

if (mkFigs == 1)
    
% % Surface elevation last Timeslice
%       figure
%       plot(t,real(Env.*exp(-2*1i*pi*f0*t)),'black');  %%changed sign$
%       hold on
%       plot(t,real(Env_WO.*exp(-2*1i*pi*f0*t)),'b--');  %%changed sign
%        plot(t,real(Env_WM.*exp(-2*1i*pi*f0*t)),'r--');  %%changed sign$
%       legend('no wind','weak wind (WO)', 'strong wind (WM)')
%       title('Surface elevation last timeslice'); 
%       
% % Envelopes at last time slice
%       figure
%       hold on
%       plot(t,real(a_new),'g');  %%changed sign
%       plot(t,real(Env),'b');  %%changed sign
%       legend('analytical','Runga Kutta')
%       title('Enelvipes at last time slice');
%       print(gcf, '-depsc2', 'Eta.eps')
  
% Characteristic wave height
      figure
      plot(Eta,Hs_sd,'b','linewidth',2)
      hold on
      plot(Eta,Hs_sdWO,'r--','linewidth',2)
      plot(Eta,Hs_sdWM,'k','linewidth',2)
      xlabel('x [m] ','fontsize',20)
      ylabel('Hs [m] ','fontsize',20)
      legend('no wind','wO', 'wM')
      set(gca,'fontsize',16)
%       print(gcf, '-depsc2', 'Hs.eps')
      title('Characteristic wave height different models')

% Envelope Spectrum
      f_mid = round(nFreq/2);
      a_FT0shift = circshift(a_FT0,[1 f_mid]);
      a_FTshift = circshift(a_FT,[1 f_mid]);

      figure
      plot(f(f_mid-200:f_mid+200)-f(f_mid), abs(a_FTshift(f_mid-200:f_mid+200)),'r','linewidth',2)
      hold on
      plot(f(f_mid-200:f_mid+200)-f(f_mid), abs(a_FT0shift(f_mid-200:f_mid+200)),'b+-','linewidth',1)
      % plot(f(1:300), abs(Y0(1:300)),'r') 
      % hold on
      % plot(f(1:300), abs(yf4(1:300)),'b+-') 
      title('Envelope Spectrum')
      xlabel('Frequency (Hz)')
      ylabel('|Y(f)|')
%       print(gcf, '-depsc2', 'Spectrum.eps')

% Momentum and Norm

    figure
    % momentum
        subplot(1,2,1)
        plot(Z0,Momentum*1e4) % go from m^2 to cm^2
        xlabel('distance z (m)')
        ylabel('Momentum (cm^2)')
    % norm
        subplot(1,2,2)
        plot(Z0,Norm*1e4) % go from m^2 to cm^2
        xlabel('distance z (m)')
        ylabel('Norm (cm^2 s^{-1})')

%% Surface Elevation Spectrum, for different models
      etaTFin = real(Env2D(end,:).*exp(-2*pi*1i*f0*t)); %%changed sign
      etaTFinO = real(Env2D_WO(end,:).*exp(-2*pi*1i*f0*t)); %%changed sign
      etaTFinM = real(Env2D_WM(end,:).*exp(-2*pi*1i*f0*t)); %%changed sign

      etawFin = FourierT(etaTFin,dt);
      etawFinO = FourierT(etaTFinO,dt);
      etawFinM = FourierT(etaTFinM,dt);

      figure   
      title('Single-Sided Amplitude Spectrum of a(t)')
      hold on
      
      % no wind
          subplot(3,1,1)    
          nm = 300; %%nm = 200; %%in old simulations
    %       plot(f(1:nm), abs(eta(1:nm)),'b')     
          plot(f(1:nm), abs(etawFin(1:nm)),'r--')        
          xlabel('Frequency (Hz)')
          ylabel('|Y(f)|')
          title('No Wind')
      
      % weak wind (Onorato)
          subplot(3,1,2)
    %       plot(f(1:nm), abs(zetaw(1:nm)),'b') 
          hold on
          plot(f(1:nm), abs(etawFinO(1:nm)),'r--') 
          xlabel('Frequency (Hz)')
          ylabel('|Y(f)|')
          title('Weak Wind')
      
      % strong wind (Maura)
          subplot(3,1,3)
    %       plot(f(1:nm), abs(zetaw(1:nm)),'b') 
          hold on
          plot(f(1:nm), abs(etawFinM(1:nm)),'r--') 
          xlabel('Frequency (Hz)')
          ylabel('|Y(f)|')
          title('Strong Wind')
%           print(gcf, '-depsc2', 'Spectrum.eps')

%% 2D plot envelope and line cuts Real Space
     
      % make coarser 2d envlope vectors for 3d plotting   
        Env2D_coarse = Env2D(:,1:20:end);
        Env2D_WO_coarse = Env2D_WO(:,1:20:end);
        Env2D_WM_coarse = Env2D_WM(:,1:20:end);
        Env2D_cg_coarse = Env2D_cg(:,1:20:end);
        t_coarse=linspace(min(t),max(t),size(Env2D_coarse,2));
      % NO WIND  
          % 2d envelope plots 
              figure
              subplot(1,2,1)
              imagesc(t,Z0,abs(Env2D_coarse))
              xlabel('time t (s)')
              ylabel('distance Z (m)')
              title('RK envelope [m]')
              colorbar

              subplot(1,2,2)
              surf(t_coarse,Z0,abs(Env2D_coarse))
              xlabel('time t (s)')
              ylabel('distance Z (m)')
              title('RK envelope [m]')
              light
              shading interp
              lighting phong
              colorbar

          % 2d evenlope NO WIND in cg frame
              figure
              subplot(1,2,1)
              imagesc(t,Z0,abs(Env2D_cg_coarse))
              xlabel('time t (s)')
              ylabel('distance Z (m)')
              title('RK envelope [m]')
              colorbar

               subplot(1,2,2)
               surf(t_coarse,Z0,abs(Env2D_cg_coarse))
              xlabel('time t (s)')
              ylabel('distance Z (m)')
              title('RK envelope [m]')
              light
              shading interp
              lighting phong
              colorbar
      
         % linecuts 
              % @  t=0 
              figure
              t_0 =  find(abs(t-(x_f/c_g))<dt, 1, 'first')
              subplot(1,2,1)
              plot(Z0,abs(Env2D(:,t_0)))
              xlabel('distance Z')
              ylabel('envelope (m)')
              title(['envelope surface height [m] @', num2str(t(t_0))])

              %@ Z=x (focus point)
              subplot(1,2,2)
              Z_x=  find(abs(Z0-abs(x_f))<outZ, 1, 'first');
              plot(t,abs(Env2D(Z_x,:)))
              xlabel('time t (s)')
              ylabel('envelope (m)')
                title('RK envelope height a [m] @ Z=x ')

         % 2d plots Analytical SSFT
              figure
              imagesc(t,Z0,abs(Env2D_Test))
              xlabel('time t (s)')
              ylabel('distance Z (m)')
              title('Analytical SSFT a  ')
              colorbar

         % linecuts Analytical SSFT
              % @  t=0 
              figure
              t_0 =  find(abs(t-(x_f/c_g))<dt, 1, 'first')
              subplot(1,2,1)
              plot(Z0,abs(Env2D_Test(:,t_0)))
              xlabel('distance Z')
              ylabel('envelope (m)')
              title(['envelope surface height [m] @', num2str(t(t_0))])

              %@ Z=x (focus point)
              subplot(1,2,2)
              Z_x=  find(abs(Z0-abs(x_f))<outZ, 1, 'first');
              plot(t,abs(Env2D_Test(Z_x,:)))
              xlabel('time t (s)')
              ylabel('envelope (m)')
              title('Analytical SSFT a [m] @ Z=x ')
      
      % STRONG WIND (Wind Maura
          % 2D plot envelope
              figure
              imagesc(t,Z0,abs(Env2D_WM))
              xlabel('time t (s)')
              ylabel('distance Z (m)')
              title('RK WM  ')
              colorbar
      
          % linecuts   
              figure
              t_0 =  find(abs(t-(x_f/c_g))<dt, 1, 'first')
              subplot(1,2,1)
              plot(Z0,abs(Env2D_WM(:,t_0)))
              xlabel('distance Z')
              ylabel('envelope (m)')
              title(['envelope surface height [m] @', num2str(t(t_0))])

              %@ Z=x (focus point)
              subplot(1,2,2)
              Z_x=  find(abs(Z0-abs(x_f))<outZ, 1, 'first');
              plot(t,abs(Env2D_WM(Z_x,:)))
              xlabel('time t (s)')
              ylabel('envelope (m)')
              title('Analytical SSFT a [m] @ Z=x ')

%% Non-Akhmediev Breather plots
      if (akhm ~=2)
      figure
      plot(Eta,nombre_RW,'b','linewidth',2)
      hold on
      plot(Eta,nombre_RW_WO,'r-.','linewidth',2)
      plot(Eta,nombre_RW_WM,'k--','linewidth',2)
      legend('no wind','wO', 'wM','Location','NorthWest')
      title('Number of RW as a function of distance')
      xlabel('x [m]')
      ylabel('Number of Rw')
      print(gcf, '-depsc2', 'RWdistance.eps')
      %
      figure
      hist(HTOT,100)
      h1 = findobj(gca,'Type','patch');
      hold on
      hist(HTOT_WO,100)
      h = findobj(gca,'type','patch');
      h2 = setdiff(h,h1);
      hist(HTOT_WM,100)
      h = findobj(gca,'type','patch');
      h3 = setdiff(h,h1);
      h4 = setdiff(h3,h2);
      set(h1,'FaceColor','g','EdgeColor','g');
      set(h2,'FaceColor','b','EdgeColor','b');
      set(h4,'FaceColor','r','EdgeColor','r');
      plot([HsMAX HsMAX], [0 2000],'k--')
      plot([2.2*HsMAX 2.2*HsMAX], [0 2000],'k:')
      h = findobj(gca,'Type','patch');
      set(h,'FaceAlpha',0.5);
      print(gcf, '-depsc2', 'histRW.eps')
      print(gcf, '-depsc2', 'histRW.fig')
      end
end

%% Calculate Final Surface elevation?

% ?? Final Zeta
%         zetaTFin = real(Eta2D(150,:).*exp(-2*pi*1i*f0*t));  %%changed sign
%         zetaTFinO = real(Eta2D_WO(150,:).*exp(-2*pi*1i*f0*t));  %%changed sign
%         zetaTFinM = real(Eta2D_WM(150,:).*exp(-2*pi*1i*f0*t));  %%changed sign
% 
%         T0 = 1/f0;
%         Nf=4096*2;
%         f1=3/T0;
%         df1=f1/Nf;
%         fi=df1:df1:f1;
%         c1=0; s=0; xx=0;
%         sWO=0; xxWO=0;
%         sWM=0; xxWM=0;
%         a=zeros(size(fi)); phi=a;
%         aWO=zeros(size(fi)); phiWO=a;
%         aWM=zeros(size(fi)); phiWM=a;
%         for ff=fi
%             omega=2*pi*ff;
%             c1=c1+1;
%             if mod(c1,100) == 0
%                 c1
%             end
%             %fft complexe
%             [a(c1),phi(c1)]=fftok1(zetaTFin' ,1/ff,1/dt,df1);
%             xx=xx+a(c1)*cos(omega*(t-t(1))-phi(c1));
%             s=s+a(c1)*cos(omega*(t-t(1))-phi(c1));
%             [aWO(c1),phiWO(c1)]=fftok1(zetaTFinO' ,1/ff,1/dt,df1);
%             xxWO=xxWO+aWO(c1)*cos(omega*(t-t(1))-phiWO(c1));
%             sWO=sWO+aWO(c1)*cos(omega*(t-t(1))-phiWO(c1));
%             [aWM(c1),phiWM(c1)]=fftok1(zetaTFinM' ,1/ff,1/dt,df1);
%             xxWM=xxWM+aWM(c1)*cos(omega*(t-t(1))-phiWM(c1));
%             sWM=sWM+aWM(c1)*cos(omega*(t-t(1))-phiWM(c1));
%         end
% 
%         figure
%         plot((fi-f0)/f0/ak,aWM*df*T0, 'g-+','linewidth',2)
%         hold on
%         plot((fi-f0)/f0/ak,aWO*df*T0, 'r','linewidth',2)
%         plot((fi-f0)/f0/ak,a*df*T0,'b','linewidth',2)
%         %axis([-10 10 0 0.037])
%         grid on
%         xlabel('(f-f0)/(f0 ak)','fontsize',20)
%         ylabel('fft elevation ','fontsize', 20)
%         set(gca,'fontsize',16)
%         print(gcf, '-depsc2', 'surfaceElevationSpectra_t150.eps')

%% end
%
%%%%%%
toc
%clear all
%close all
return
%exit

nt = 3;

filenametosave =strcat('SymbolicNSCavitDiffTTranspose_nt',num2str(nt))

%%
mode = 0;  % Need to think about this parameter. 

nvar = 4; 
cylindrical  = false;

%% Generate Symbolic Equations.
tic
ks =  0:(nt-1); % Spectral modes in "academic" ordering

% Declare Unknowns 
syms Unr Uwr Uor Uer Usr Unwr Uner User Uswr [1 nt]
syms Uni Uwi Uoi Uei Usi Unwi Unei Usei Uswi [1 nt]
syms Vnr Vwr Vor Ver Vsr Vnwr Vner Vser Vswr [1 nt]
syms Vni Vwi Voi Vei Vsi Vnwi Vnei Vsei Vswi [1 nt]
syms Pnr Pwr Por Per Psr Pnwr Pner Pser Pswr [1 nt]
syms Pni Pwi Poi Pei Psi Pnwi Pnei Psei Pswi [1 nt]
syms Wnr Wwr Wor Wer Wsr Wnwr Wner Wser Wswr [1 nt]
syms Wni Wwi Woi Wei Wsi Wnwi Wnei Wsei Wswi [1 nt]
syms Fq % 1/Period
%
Unknownsv = [...
             Unwr   Unwi    Unr    Uni    Uner    Unei    ...
             Uwr    Uwi     Uor    Uoi    Uer     Uei     ...
             Uswr   Uswi    Usr    Usi    User    Usei    ...
             Vnwr   Vnwi    Vnr    Vni    Vner    Vnei    ...
             Vwr    Vwi     Vor    Voi    Ver     Vei     ...
             Vswr   Vswi    Vsr    Vsi    Vser    Vsei    ...
             Pnwr   Pnwi    Pnr    Pni    Pner    Pnei    ...
             Pwr    Pwi     Por    Poi    Per     Pei     ...
             Pswr   Pswi    Psr    Psi    Pser    Psei    ...
             Wnwr   Wnwi    Wnr    Wni    Wner    Wnei    ...
             Wwr    Wwi     Wor    Woi    Wer     Wei     ...
             Wswr   Wswi    Wsr    Wsi    Wser    Wsei    ...
            ]; 
assume(Unknownsv,'real')
   
%
if nt > 1
    Unknownsf = [Fq];
    assume(Unknownsf,'real')
else
    Unknownsf = [];
    Fq = 0;
end
%
% Switch to complex to facilitate writing Equations.

Un = Unr + 1i*Uni; Uw =  Uwr + 1i*Uwi;        Uo =  Uor + 1i*Uoi;
Ue = Uer + 1i*Uei; Us =  Usr + 1i*Usi; 
Unw = Unwr + 1i*Unwi; Une = Uner + 1i*Unei; 
Use = User + 1i*Usei; Usw = Uswr + 1i*Uswi;

Vn = Vnr + 1i*Vni; Vw =  Vwr + 1i*Vwi;        Vo =  Vor + 1i*Voi;
Ve = Ver + 1i*Vei; Vs =  Vsr + 1i*Vsi; 
Vnw = Vnwr + 1i*Vnwi; Vne = Vner + 1i*Vnei; 
Vse = Vser + 1i*Vsei; Vsw = Vswr + 1i*Vswi;

Pn = Pnr + 1i*Pni; Pw =  Pwr + 1i*Pwi;        Po =  Por + 1i*Poi;
Pe = Per + 1i*Pei; Ps =  Psr + 1i*Psi; 
Pnw = Pnwr + 1i*Pnwi; Pne = Pner + 1i*Pnei; 
Pse = Pser + 1i*Psei; Psw = Pswr + 1i*Pswi;

Wn = Wnr + 1i*Wni; Ww =  Wwr + 1i*Wwi;        Wo =  Wor + 1i*Woi;
We = Wer + 1i*Wei; Ws =  Wsr + 1i*Wsi; 
Wnw = Wnwr + 1i*Wnwi; Wne = Wner + 1i*Wnei; 
Wse = Wser + 1i*Wsei; Wsw = Wswr + 1i*Wswi;


% Extra local variables but not unknowns 
syms xcw xco xce xuw xuo xue 
syms ycs yco ycn yvs yvo yvn 
extraxvars = [xcw xco xce xuw xuo xue];  % need to be coherent with vector xparams
extrayvars = [ycs yco ycn yvs yvo yvn];  % need to be coherent with vector yparams
assume([extraxvars extrayvars],'real');
% Variables for cylindrical geometry
syms rcw rco rce ruw ruo rue 
% 
Pi = sym(pi);
% Extra params
syms visc kappa buoy
extraparams = [visc kappa buoy];
assume(extraparams,'real');

%%   EQUATIONS AND BC

% Some extra definitions and tricks to switch to cylindrical coordinates 

dxco = xuo - xuw;
dyco = yvo - yvs;
dxuo = xce - xco;
dyvo = ycn - yco;
dyuo = dyco;
dxvo = dxco;

if cylindrical
    rcw = xcw; rco = xco; rce = xce;
    ruw = xuw; ruo = xuo; rue = xue;
else
    rcw = 1; rco = 1; rce = 1; 
    ruw = 1; ruo = 1; rue = 1;
end

% Interpolated quantities at interfaces
Uemidx = (Ue+Uo)/2;
Vnmidx = ( Ve*(xuo-xco) + Vo*(xce-xuo) )/(xce-xco);
Unmidy = ( Un*(yvo-yco) + Uo*(ycn-yvo) )/(ycn-yco);
Uwmidx = (Uo+Uw)/2;
Usmidy = ( Uo*(yvs-ycs) + Us*(yco-yvs)  )/(yco-ycs);
Vsmidx = ( Vse*(xuo-xco) + Vs*(xce-xuo) )/(xce-xco); 

Vnmidy  = (Vn+Vo)/2;
Uemidy  = ( Un*(yvo-yco)  + Uo*(ycn-yvo) )/(ycn-yco);
Vsmidy  = (Vs+Vo)/2;
Uwmidy  = ( Unw*(yvo-yco) + Uw*(ycn-yvo) )/(ycn-yco);          
Vwmidx  = ( Vo*(xuw-xcw)  + Vw*(xco-xuw) )/(xco-xcw);
Vemidx  = ( Ve*(xuo-xco)  + Vo*(xce-xuo) )/(xce-xco);

Wemidy  = ( Wn*(yvo-yco)  + Wo*(ycn-yvo) )/(ycn-yco);

Wwmidx  = ( Wo*(xuw-xcw) + Ww*(xco-xuw) )/(xco-xcw);
Wemidx  = ( We*(xuo-xco) + Wo*(xce-xuo) )/(xce-xco);
Wnmidy =  ( Wn*(yvo-yco) + Wo*(ycn-yvo) )/(ycn-yco);
Wsmidy =  ( Wo*(yvs-ycs) + Ws*(yco-yvs) )/(yco-ycs);


%% U Equation 

% Viscous coefficients
axup = 2*rce/(rce+rco)/(xue-xuo)/dxuo;
axum = 2*rco/(rce+rco)/(xuo-xuw)/dxuo;
ayup = 1/(ycn-yco)/dyuo;
ayum = 1/(yco-ycs)/dyuo;

% Skew symmetric scheme
EqU  = ... 
       Fq*2*1j*Pi*ks.*Uo ...
     + (2*rce/(rce+rco)*sspconv(Uemidx,Uemidx) - 2*rco/(rce+rco)*sspconv(Uwmidx,Uwmidx) )/dxuo ...
               + (sspconv(Unmidy,Vnmidx) - sspconv(Usmidy,Vsmidx))/dyuo  ...
                          + (Pe-Po)/dxuo  - buoy*Wemidx  ...
      -visc*( - (axup    + axum    + ayup    + ayum   )*Uo ...
                + (axup*Ue + axum*Uw + ayup*Un + ayum*Us)        ); 
     
if cylindrical
    EqU = EqU - sspconv(Wemidx,Wemidx)/ruo - visc*( -Uo/ruo^2 ) ;
    if mode >= 1
        EqU = EqU + 1i*mode*sspconv(Uo,Wemidx)/ruo ...
                  - visc*(-2i*mode*Wemidx/ruo^2 - mode^2*Uo/ruo^2);
    end
end
   
EqUn = (Uo + Us)/2 - 0; % No slip
%EqUn = (Uo - Us)/2 - 0;  % Free Slip
EqUs = (Un + Uo)/2 - 0; 
EqUe = Uo - 0;
EqUw = Uo - 0;

EqUfake = Uo - 0;

%% V Equation

% Viscous coefficients
axvp = ruo/rco/(xce-xco)/dxvo;
axvm = ruw/rco/(xco-xcw)/dxvo;
ayvp = 1/(yvn-yvo)/dyvo;
ayvm = 1/(yvo-yvs)/dyvo;

% Skew Sym

EqV  = ... 
       Fq*2*1j*Pi*ks.*Vo ...
   + ( sspconv(Uemidy,Vemidx)*ruo/rco - sspconv(Uwmidy,Vwmidx)*ruw/rco )/dxvo ...
                         + (sspconv(Vnmidy,Vnmidy) - sspconv(Vsmidy,Vsmidy))/dyvo  ... 
                         + (Pn-Po)/dyvo ... % - buoy*Wemidy ...
         -visc*( - (axvp    + axvm    + ayvp    + ayvm   )*Vo ...
                  + (axvp*Ve + axvm*Vw + ayvp*Vn + ayvm*Vs)        ); 
    if mode >= 1
        EqV = EqV + 1i*mode*sspconv(Vo,Wnmidy)/rco - visc*( - mode^2*Vo/rco^2);
    end
              
EqVn = Vo - 0.;
EqVs = Vo - 0.;
EqVe = (Vo+Vw)/2 - 0;
EqVw = (Ve+Vo)/2 - 0;

EqVfake = Vo - 0.;

%% P Equation  (div(V) = 0) 
EqP  = (Uo*ruo/rco-Uw*ruw/rco)/dxco  + (Vo-Vs)/dyco; 
    if mode >= 1
        EqP = EqP + 1i*mode*Wo/rco;
    end
 
EqPcte = Po - 0.; 
EqPn = Po - 0;
EqPs = Po - 0;
EqPe = Po - 0;
EqPw = Po - 0;

EqPsw = Po - 0;
EqPse = Po - 0;
EqPnw = Po - 0;
EqPne = Po - 0;


%% W Equation 

% Diffusive coefficients
axwp = ruo/rco/(xce-xco)/dxco;
axwm = ruw/rco/(xco-xcw)/dxco;
aywp = 1/(ycn-yco)/dyco;
aywm = 1/(yco-ycs)/dyco;

EqW  = ... 
       Fq*2*1j*Pi*ks.*Wo ...
  + (sspconv(Uo,Wemidx)*ruo-sspconv(Uw,Wwmidx)*ruw)/dxco/rco ...
         + (sspconv(Vo,Wnmidy) - sspconv(Vs,Wsmidy))/dyco ... 
         -kappa*( - (axwp    + axwm    + aywp    + aywm   )*Wo ...
                        + (axwp*We + axwm*Ww + aywp*Wn + aywm*Ws)   ); 

if cylindrical
    EqW = EqW + sspconv(Uwmidx,Wo)/rco - kappa*( -Wo/rco^2 ) ;
    if mode >= 1
        EqW = EqW + 1i*mode*sspconv(Wo,Wo)/rco +1i*mode*Po/rco ... 
            - kappa*( - mode^2*Wo/rco^2 + 2i*mode*Uwmidx/rco^2 );
    end
end
 
% EqWn = (Wo+Ws)/2 - (0.5 - xco).*fft(ones(1,nt))/nt;   % Conducting wall
% %EqWn = (Wo - Ws)/2 - 0*xco; % Zero flux
% EqWs = (Wn+Wo)/2 - (0.5 - xco).*fft(ones(1,nt))/nt; % Conducting wall
% EqWe = (Ww+Wo)/2 - (-0.5)*fft(ones(1,nt))/nt;  %
%EqWw = (Wo+We)/2 - (+0.5)*fft(ones(1,nt))/nt;
EqWn = (Wo+Ws)/2 - (-0.5).*fft(ones(1,nt))/nt;   % Conducting wall
EqWs = (Wn+Wo)/2 - (+0.5).*fft(ones(1,nt))/nt; % Conducting wall
EqWe = (Ww+Wo)/2 - (+0.5 - yco).*fft(ones(1,nt))/nt;  %
EqWw = (Wo+We)/2 - (+0.5 - yco).*fft(ones(1,nt))/nt;


EqWsw = Wo - 0;
EqWse = Wo - 0;
EqWnw = Wo - 0;
EqWne = Wo - 0;

%% Freq equation just impose some derivative in time to be 0 at t=0.
EqFq = sum(ks.*Uo) - sum(ks.*conj(Uo)) - 0 ; % This equation is purely imaginary.

%%         
% Make an array of Equations (must be coherent with the Nodes list see numerical code) 
Equationsvcplx = [...
                       EqUn,   ...
                EqUw,  EqU,    EqUe, ...
                       EqUs, ...
               EqUfake,...
                       EqVn,   ...
                EqVw,  EqV,    EqVe, ...
                       EqVs, ...
               EqVfake,...
               EqPnw,    EqPn,   EqPne,...
               EqPw,     EqP,    EqPe, ...
               EqPsw,    EqPs,   EqPse,...
               EqPcte ...
               EqWnw,    EqWn,   EqWne,...
               EqWw,     EqW,    EqWe, ...
               EqWsw,    EqWs,   EqWse...
               ];
% We now store in Equations real and imag of each complex equations 
Equations = [];
for neq = 1:length(Equationsvcplx)
    Equations = [Equations real(Equationsvcplx(neq))  imag(Equationsvcplx(neq)) ];
end
if nt>1 
    Equations = [Equations imag(EqFq)]; % Add the imaginary part of EqFq, the real part is zero anyway.
end

t=toc;disp(' ');disp(['Generating ',num2str(length(Equations)),' Symbolics Equations in ... ',num2str(t),' sec.']); 

%% Generate Numerical functions for F and Jacobian from symbolic functions
tic

Unknowns = [Unknownsv,Unknownsf];
nbEquations = length(Equations);
JEquationsv = jacobian(Equations,Unknownsv);
JEquationsf = jacobian(Equations,Unknownsf);
JEquations  = jacobian(Equations,Unknowns);
%%
%
%
% Clean-up all unused unknowns (Jacobian is 0) 
% Assign column values from Jacobian and stencil
% as well as the correct inputs for function.
jcoleqv = cell(nbEquations,1);
jcoleqf = cell(nbEquations,1);
xusedparam = cell(nbEquations,1);
yusedparam = cell(nbEquations,1);
usedparam  = cell(nbEquations,1);
finputs = cell(nbEquations,1);

for ind = 1:nbEquations
   jcoleqv{ind}    = ismember(Unknownsv,   symvar(Equations(ind))); % remove unused unknowns
   jcoleqf{ind}    = ismember(Unknownsf,   symvar(Equations(ind))); % remove unused unknowns
   xusedparam{ind} = ismember(extraxvars,  symvar(Equations(ind)));
   yusedparam{ind} = ismember(extrayvars,  symvar(Equations(ind)));
   usedparam{ind}  = ismember(extraparams, symvar(Equations(ind)));
   %
   finputs{ind} = [Unknowns(   ismember(Unknowns,    symvar(Equations(ind)))), ...
                   extraxvars( ismember(extraxvars,  symvar(Equations(ind)))), ...
                   extrayvars( ismember(extrayvars,  symvar(Equations(ind)))), ...
                   extraparams(ismember(extraparams, symvar(Equations(ind))))  ];
end

% Generate numerical functions in .m files from symbolic expressions
listF  =  cell(nbEquations,1);
listdF =  cell(nbEquations,1);

for ind = 1:nbEquations
   % ind
    listF{ind} = matlabFunction(Equations(ind),'vars',{finputs{ind}});
    tempo = JEquations(ind,:);tempo = tempo( ismember(Unknowns, symvar(Equations(ind))) ); % Remove zeros on each Jacobian line
    listdF{ind} = matlabFunction(tempo,'vars',{finputs{ind}});
end

save(filenametosave,'listF','listdF','usedparam', 'xusedparam', 'yusedparam', 'jcoleqv', 'jcoleqf' ,'nt','nvar', 'nbEquations', 'cylindrical')
%load('AllFdF.mat')

t=toc;disp(' ');disp(['Generating the numerical functions from the ',num2str(nbEquations), ...
                       ' symbolic equations in ... ',num2str(t),' sec.']);  

%return


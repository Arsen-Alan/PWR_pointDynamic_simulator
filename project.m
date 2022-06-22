%
% POINT DYNAMICS SIMULATOR
%

% FLAGS

promptFeedback  = 1; % 1 for ON 0 for OFF
delayedFeedback = 0; % 1 for ON 0 for OFF
source          = 0; % 1 for ON 0 for OFF
reactivity      = 1; % 0 OFF, 1 STEP, 2 TRIANG, 3 HARM

plotN           = 1; % 0 don't plot, 1 plot (neutron density)
plotRho         = 1; % 0 don't plot, 1 plot (reactivity)
plotT           = 1; % 0 don't plot, 1 plot (temperatures)
plotCon         = 1; % 0 don't plot, 1 plot (precursor concentrations)

plotRhoTot      = 1; % 0 don't plot, 1 plot (total reactivity)
plotRhoFeedF    = 1; % 0 don't plot, 1 plot (feedback reactivity)
plotRhoFeedC    = 1; % 0 don't plot, 1 plot (feedback reactivity)
plotRhoExt      = 1; % 0 don't plot, 1 plot (external reactivity)

plotTfuel       = 1; % 0 don't plot, 1 plot (fuel temp)
plotTCoolant    = 1; % 0 don't plot, 1 plot (coolant temp)

plotTypeN       = 0; % 0 normal, 1 semilogx, 2 semilogy, 3 loglog plot for n density
plotTypeRho     = 0; % 0 normal, 1 semilogx, 2 semilogy, 3 loglog plot for reactivity
plotTypeT       = 0; % 0 normal, 1 semilogx, 2 semilogy, 3 loglog plot for temps
plotTypeCon     = 0; % 0 normal, 1 semilogx, 2 semilogy, 3loglog plot for Prec concentrations

% TRANSIENT PARAMETERS 

global q
global react
global wSin
global tStart

q      = 200       ; % [1/m3s]  External source
react  = 0.006019  ; %          Reactivity (for triang is the peak)
wSin   = 0.25    ; % [1/s]    Harmonic period
tStart = 0         ; % [s]      Time of start of the transient

alphaF = 3.24e-5; % [1/K]    Fuel (Doppler) reactivity coefficient
alphaC = -21.3e-5; % [1/K]    Moderator/coolant reactivity coefficient

tspan  = [0 0.308]  ; % [s]      Time span in seconds


% REACTOR PARAMETERS

yield  =  [0.0001745, 0.001257, 0.0011405, 0.0023495, 0.000821, 0.0002765]; % yields of the precursors of delayed neutrons
dconst =  [0.01255, 0.0307, 0.1165, 0.3125, 1.19, 3.15]; % [1/s] decay constants of the precursors of delayed neutrons

LAMBDA =  2e-4    ; % [s]      Generation time
ne     =  200     ; % [1/m3]   Neutron density at equilibrium
cpF    =  657.5   ; % [j/kg/K] Fuel specific heat
cpC    =  4200    ; % [j/kg/K] Coolant specific heat
mF     =  40000   ; % [kg]     Mass of fuel
mC     =  17095.18; % [kg]     Mass of coolant
WCe    =  24285.71; % [kg]     Coolant flow at equilibrium
aF     =  12.5e6  ; % [jm3/s]  Power-to-neutron-density factor
h      =  6.734e6 ; % [j/K/s]  Fuel-coolant heat flux coefficient
TCine  =  563.15  ; % [K]      Coolant inlet temperature at equilibrium

% Calculated parameters: 

% ne might be included here
TFe    =  TCine + aF*ne*(1/(2*WCe*cpC) + 1/h);  % fuel temperature at equilibrium
TCe    =  TCine + aF*ne/(2*WCe*cpC);            % coolant temperature at equilibrium

%
% SOLVER
%

y0 = zeros(length(yield)+3,1);  % Initial conditions (all zero due to normalization)

% ODE solver
[t,y] = ode15s(@(t,y) PKModel(t,y,yield,dconst,LAMBDA,ne,cpF,cpC,mF,mC,WCe,aF,h,alphaF,alphaC,TCine,TCe,TFe,promptFeedback, delayedFeedback, source,reactivity),tspan, y0);

NPre = length(yield); % Get number of precursors
neq  = 3 + NPre;      % Total number of equations

n  = y(:,1)*ne  + ne;     % [1/cm3]  Neutron density
Tf = y(:,2)*TFe + TFe;    % [K]      Fuel temperature
Tc = y(:,3)*TCe + TCe;    % [K]      Coolant temperature

rho      = zeros(length(t),1); % initializing total reactivity array 
feedRhoF = zeros(length(t),1); % initializing Tf feedback reactivity array
feedRhoC = zeros(length(t),1); % initializing Tc feedback reactivity array
rhoExt   = zeros(length(t),1); % initializing external reactivity array 

switch reactivity % adding the right InReactivity shape to rho
    case 1
        for i = 1:length(t)
            rhoExt(i) = InReactivityStep(t(i));
        end
        
    case 2
        for i = 1:length(t)
            rhoExt(i) = InReactivityTriang(t(i));
        end
    case 3
        for i = 1:length(t)
            rhoExt(i) = InReactivityHarm(t(i));
        end
end

if promptFeedback == 1  % adding the promptFeedback to rho if needed
    feedRhoF = feedRhoF + alphaF*TFe*y(:,2);
end
   
if delayedFeedback == 1 % adding the delayedFeedback to rho if needed
    feedRhoC = feedRhoC + alphaC*TCe*y(:,3);
end  

rho = rhoExt + feedRhoF + feedRhoC; 

C = zeros(length(t),NPre); % initializing precursors matrix

for i = 1:NPre
        Ceq  = ne*yield(i)/(dconst(i)*LAMBDA);
        C(:,i) = y(:,i+neq-NPre)*Ceq + Ceq;
end

if plotN == 1
    
    figure(1)
    switch plotTypeN 
        case 0
            plot(t,n);
        
        case 1
            semilogx(t,n)
    
        case 2
        semilogy(t,n)
    
        case 3
            loglog(t,n)
    end
    
    title('Neutron population')
    ylabel('n [1/m^3]')
    xlabel('t [s]')
    grid on
    grid minor 
    hold on
    
end


if plotRho == 1
    
    figure(2)
    hold off
    switch plotTypeRho 
        case 0
            if plotRhoTot == 1                              
                plot(t,rho/sum(yield))
                hold on 
            end
               
            if plotRhoExt == 1
                plot(t,rhoExt/sum(yield))
            end
            
            if plotRhoFeedF == 1
                plot(t,feedRhoF/sum(yield))
                hold on
            end    
            
            if plotRhoFeedC == 1
                plot(t,feedRhoC/sum(yield))
                hold on
            end 
            
                    
        case 1
            if plotRhoTot == 1                              
                semilogx(t,rho/sum(yield))
                hold on 
            end
            
            if plotRhoExt == 1
                semilogx(t,rhoExt/sum(yield))
            end
               
            if plotRhoFeedF == 1
                semilogx(t,feedRhoF/sum(yield))
                hold on
            end 
            
            if plotRhoFeedC == 1
                semilogx(t,feedRhoC/sum(yield))
                hold on
            end
            
            
    
        case 2
            if plotRhoTot == 1                              
                semilogy(t,rho/sum(yield))
                hold on 
            end
            
            if plotRhoExt == 1
                semilogy(t,rhoExt/sum(yield))
            end
               
            if plotRhoFeedF == 1
                semilogy(t,feedRhoF/sum(yield))
                hold on
            end    
            
            if plotRhoFeedC == 1
                semilogy(t,feedRhoC/sum(yield))
                hold on
            end 
            
            
    
        case 3
            if plotRhoTot == 1                              
                loglog(t,rho/sum(yield))
                hold on 
            end
            
            if plotRhoExt == 1
                loglog(t,rhoExt/sum(yield))
            end 
               
            if plotRhoFeedF == 1
                loglog(t,feedRhoF/sum(yield))
                hold on
            end  
            
            if plotRhoFeedC == 1
                loglog(t,feedRhoC/sum(yield))
                hold on
            end
            
            
    end

    title('Reactivity')
    ylabel('Reactivity [\beta]')
    xlabel('t [s]')
    legend('Total reactivity','External reactivity','Doppler feedback reactivity','Coolant feedback reactiivy')
    grid on
    grid minor
end



if plotT == 1
    
    figure(3)
    %hold off    
    switch plotTypeT 
        case 0
            if plotTfuel == 1
                plot(t,Tf)
                hold on
            end
            if plotTCoolant == 1
                plot(t,Tc)
            end
        
        case 1
            if plotTfuel == 1
                semilogx(t,y(:,2))
                hold on
            end
            if plotTCoolant == 1
                semilogx(t,y(:,3))
            end
    
        case 2
            if plotTfuel == 1
                semilogy(t,Tf)
                hold on
            end
            if plotTCoolant == 1
                semilogy(t,Tc)
            end
    
        case 3
            if plotTfuel == 1
                loglog(t,Tf)
                hold on
            end
            if plotTCoolant == 1
                loglog(t,Tc)
            end 
    end

    title('Fuel & coolant temperature')
    ylabel('T [K]')
    xlabel('t [s]')
    legend('Fuel temperature','Coolant temperature')
    grid on
    grid minor
end

if plotCon == 1
    
    figure(4)
    hold off
    for i = 1:NPre
    
        switch plotTypeN 
            case 0
                plot(t,C(:,i));
        
            case 1
                semilogx(t,C(:,i));
    
            case 2
                semilogy(t,C(:,i));
    
            case 3
                loglog(t,C(:,i));
        end
        
        hold on
    end
    
    title('Precursors concentration')
    ylabel('Concentration [1/m3]')
    xlabel('t [s]')
    legend({'Group 1','Group 2','Group 3','Group 4','Group 5','Group 6'},'NumColumns',2,'Location','northwest')
    grid on 
    grid minor
end
  

function rho = InReactivityStep(t)

    global react
    global tStart
    
    if t < tStart
        rho = 0;
    else
        rho = react;
    end
end

function rho = InReactivityTriang(t)

    global react
    global tStart
    
    if t < (tStart)
        rho = 0;
    elseif t > (tStart) && t < (tStart + 0.5)
        rho = 2*react * (t-tStart);
    elseif t > (tStart + 0.5) && t < (tStart + 1)
        rho = 2*react - 2*react * (t-tStart);
    elseif t > (tStart+1)
        rho = 0;
    else
        rho = 0;
    end

end

function rho = InReactivityHarm(t)

    global react
    global wSin
    global tStart
    
    if t < tStart
        rho = 0;
    else
        rho = react*sin(wSin*t*2*pi());
    end
end

function S = sourceDef(t)
    
    global q
    global tStart
    
    if t < tStart
        S = 0;
    else
        S = q;
    end
end

function dy = PKModel(t,y,yield,dconst,LAMBDA,ne,cpF,cpC,mF,mC,WCe,aF,h,alphaF,alphaC,TCine,TCe,TFe, promptFeedback, delayedFeedback,source,reactivity)
    
    u   = 0;
    w   = 0;

    switch reactivity
    case 0
        rho = 0;
    case 1
        rho = InReactivityStep(t);
    case 2
        rho = InReactivityTriang(t);
    case 3
        rho = InReactivityHarm(t);
    end
    
    NPre = length(yield);                              % Get number of precursors
    neq  = 3 + NPre; % Total number of equations
    
    dy  = zeros(neq,1); % Initialize RHS
    
    if promptFeedback == 1  
        rho = rho + alphaF*TFe*y(2);
    end
   
    if delayedFeedback == 1
        rho = rho + alphaC*TCe*y(3);
    end
    
    if source == 1
        S = sourceDef(t);
    else
        S = 0;
    end
   
    dy(1) = (-sum(yield)*y(1) + yield*y(4:end) + rho + rho*y(1))/LAMBDA + S/ne;
    
    dy(2) = (aF*ne*y(1)/TFe - h*y(2) + h*TCe*y(3)/TFe)/(mF*cpF);
    
    dy(3) = (h*TFe*y(2)/TCe - (2*cpC*WCe+h)*y(3))/(mC*cpC) ...
            + 2*WCe*(TCine*u - (TCe-TCine)*w)/(mC*TCe) ...
            - 2*WCe*w*(y(3) - TCine*u/TCe)/mC;
    
    for i = 1:NPre
        dy(i+3) = dconst(i)*(y(1) - y(i+3));
    end

end


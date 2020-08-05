% calculates the molar volume of phases of interest, and adds excess volume to minerals that undergo phase transitions
function molar_volume = eos(P,T,datatable,phase)

% Initializing variables from thermodynamic_database and gas constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkbar=P;
k0=datatable{phase,'k0'};
k0_1=datatable{phase,'k0_1'};
k0_2=datatable{phase,'k0_2'};
v0=datatable{phase,'v0'};
s0=datatable{phase,'s0'}/1000;                  % convert entropy to kJ/K
alpha0=datatable{phase,'alpha0'}/10^5;
alpha1=datatable{phase,'alpha1'}/10^8;
lambda=datatable{phase,'lambda'};
sum_apfu=datatable{phase,'sum_apfu'};
crit_temp=datatable{phase,'crit_temp'};
max_s=datatable{phase,'max_s'}/1000;            % convert max entropy to kJ/K
max_v=datatable{phase,'max_v'};
delta_h=datatable{phase,'delta_h'};
delta_v=datatable{phase,'delta_v'};
w=datatable{phase,'w'};
wv=datatable{phase,'wv'};
n=datatable{phase,'n'};
sf=datatable{phase,'sf'};
einstein_temperature=datatable{phase,'einstein_T'};
thermal_term=datatable{phase,'thermal_term'};
eos_type=datatable{phase,'EoS_type'};
bulk_modulus_type=datatable{phase,'K_type'};
dK_dT=datatable{phase,'dK_dT'};
delta=datatable{phase,'delta_T'};
delta_prime=datatable{phase,'delta_T_Pr'};
R = 8.3144598;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EoS thermal model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thermal pressure Tait EoS
if thermal_term == 1                                                       
    
   einstein_temperature_check = isnan(einstein_temperature);
    
    if einstein_temperature_check == 0
        
        einstein_temperature;
        
    else
        
        einstein_temperature = 10636/(s0 * 1000/sum_apfu + 6.44);
        
    end
    
    u = einstein_temperature/T;
    
    u_ambient = einstein_temperature/298.15;
    
    einstein_function_ambient = ((u_ambient.^2) * exp(u_ambient))/((exp(u_ambient) - 1).^2);
        
    thermal_pressure = ((alpha0 * k0 * (einstein_temperature/einstein_function_ambient)) * (1/(exp(u) - 1) - 1/(exp(u_ambient) - 1)));
    
% Berman thermal model
elseif thermal_term == 2
    
    molar_volume_isobaric = v0 * (1 + alpha0 * (T - 298.15) + 0.5 * alpha1 * (T - 298.15).^2);       
    
% Kroll thermal model  
elseif thermal_term == 3                                                      
    
   einstein_temperature_check = isnan(einstein_temperature);
    
    if einstein_temperature_check == 0
        
        einstein_temperature;
        
    else
        
        einstein_temperature = 10636/(s0 * 1000/sum_apfu + 6.44);
        
    end
    
    u = einstein_temperature/T;
    
    u_ambient = einstein_temperature/298.15;
    
    einstein_function_ambient = ((u_ambient.^2) * exp(u_ambient))/((exp(u_ambient) - 1).^2);
    
    A = ((alpha0 * (einstein_temperature/einstein_function_ambient)) * (1/(exp(u) - 1) - 1/(exp(u_ambient) - 1)));
   
    
    if bulk_modulus_type == 2
              
        B = -1/(delta * (delta + 2)); % Hellfrich-Connolly approximation, delta substitution for K'
                
        molar_volume_isobaric = v0 * (-delta + (1 + delta) * (1 - ((delta * (delta + 2))/(delta + 1)) * A).^B); % Hellfrich-Connolly approximation, delta substitution for K' 
        
    else
        
        B = -1/(k0_1 * (k0_1 + 2));
        
        molar_volume_isobaric = v0 * (-k0_1 + (1 + k0_1) * (1 - ((k0_1 * (k0_1 + 2))/(k0_1 + 1)) * A).^B);
        
    end
         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bulk modulus variation with temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% linear approximation to bulk modulus variation with temperature 
if bulk_modulus_type == 1
    
    k0_T = k0 + dK_dT * (T - 298.15);                                          

% Hellfrich-Connolly approximation, delta substitution for K'     
elseif bulk_modulus_type == 2
              
    k0_T = k0 * (v0/molar_volume_isobaric).^delta;
    k0_1 = k0_1 * (molar_volume_isobaric/v0).^delta_prime;
    k0_2 = -1.0 * ((4 - k0_1) * (3 - k0_1) + 35/9) / k0_T;
               
else

% No bulk modulus variation with temperature
k0_T = k0;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% EoS type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tait EoS
if eos_type == 1                                                            
    
     a = (1 + k0_1)/(1 + k0_1 + k0_T * k0_2);                              % TEOS constant a
     b = k0_1/k0_T - k0_2/(1 + k0_1);                                      % TEOS constant b
     c = (1 + k0_1 + k0_T * k0_2)/(k0_1.^2 + k0_1 - k0_T * k0_2);          % TEOS constant c
     
     % thermal pressure Tait EoS 
     if thermal_term == 1
         
        molar_volume = v0 * (1 - a * (1 - (1 + b * (pkbar - thermal_pressure)).^-c));
     
     % Berman thermal model
     elseif thermal_term == 2
         
        molar_volume = molar_volume_isobaric * (1 - a * (1 - (1 + b * pkbar).^-c)); 
      
     % Kroll thermal model    
     elseif thermal_term == 3
         
        molar_volume = molar_volume_isobaric * (1 - a * (1 - (1 + b * pkbar).^-c)); 
     
     end

% Birch-Murnaghan EoS
elseif eos_type == 2                                                       
      
      options = optimoptions(@fsolve, 'display','off', 'tolfun', 1e-10);
      
      % thermal pressure Tait EoS        
      if thermal_term == 1                                             
                 
         Obj_func = @(molar_volume) 3 * k0_T * (((v0/molar_volume).^(2/3) - 1)/2)...
             * (1 + 2 * (((v0/molar_volume).^(2/3) - 1)/2)).^(5/2) * (1 + (3/2) * (k0_1 - 4)...
             * (((v0/molar_volume).^(2/3) - 1)/2) + (3/2)...
             * (k0_T * k0_2 + (k0_1 - 4) * (k0_1 - 3) + (35/9)) * (((v0/molar_volume).^(2/3) - 1)/2).^2) - (pkbar - thermal_pressure);

         molar_volume = fsolve(Obj_func, 0.001, options);     
      
      % Berman thermal model
      elseif thermal_term == 2                                                  
       
         Obj_func = @(molar_volume) 3 * k0_T * (((molar_volume_isobaric/molar_volume).^(2/3) - 1)/2)...
             * (1 + 2 * (((molar_volume_isobaric/molar_volume).^(2/3) - 1)/2)).^(5/2) * (1 + (3/2) * (k0_1 - 4)...
             * (((molar_volume_isobaric/molar_volume).^(2/3) - 1)/2) + (3/2)...
             * (k0_T * k0_2 + (k0_1 - 4) * (k0_1 - 3) + (35/9)) * (((molar_volume_isobaric/molar_volume).^(2/3) - 1)/2).^2) - pkbar;

         molar_volume = fsolve(Obj_func, 0.001, options);    
      
      % Kroll thermal model
      elseif thermal_term == 3
                                
         Obj_func = @(molar_volume) 3 * k0_T * (((molar_volume_isobaric/molar_volume).^(2/3) - 1)/2)...
             * (1 + 2 * (((molar_volume_isobaric/molar_volume).^(2/3) - 1)/2)).^(5/2) * (1 + (3/2) * (k0_1 - 4)...
             * (((molar_volume_isobaric/molar_volume).^(2/3) - 1)/2) + (3/2)...
             * (k0_T * k0_2 + (k0_1 - 4) * (k0_1 - 3) + (35/9)) * (((molar_volume_isobaric/molar_volume).^(2/3) - 1)/2).^2) - pkbar;

         molar_volume = fsolve(Obj_func, 0.001, options);  
        
      end
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% phase transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% holland and powell, 2011 landau theory excess volume calculation with angel 2017 curved boundary model beta quartz bare-phase parameters
if strcmp(phase, 'Quartz')
    
    aT = -1.884705 * 10.^-3;
    aTH = 0.706630 * 10.^-3;       
    beta = 0.5;                 
    beta_high = 0.25;           
    
    pgpa = pkbar/10;
        
    crit_temp_P = crit_temp + 270 * pgpa - 10 * pgpa.^2;  
              
            if T > crit_temp_P
                          
                Ev = aTH * abs((crit_temp_P - T)).^beta_high;
                                            
            else  
                             
                Ev = aT * abs((crit_temp_P - T)).^beta;
            
            end
        
                molar_volume = molar_volume * (1 + Ev);
                        
%% holland and powell, 2011 landau theory excess volume calculation
elseif lambda == 1
    
    crit_temp_P = crit_temp + (max_v/max_s) * pkbar;
    
    Q2_ambient = sqrt(1 - 298.15/crit_temp);
        
            if T > crit_temp_P
            
                Q2 = 0;
               
            else  
                
                Q2 = sqrt((crit_temp_P - T)/crit_temp);               
            
            end
         
            excess_v = (max_v) * (Q2_ambient - Q2);
            
            molar_volume = molar_volume + excess_v;
        
%% holland and powell, 1998 Bragg-Williams theory excess volume calculation          
elseif lambda == 2
    
    options = optimoptions(@fsolve, 'display','off', 'tolfun', 1e-10); 
    
    Obj_func = @(x) (sf * (n/(n + 1)) * R * T * (log(n - n * x) + log(1 - x) - log(1 + n * x) - log(n + x) ) + w * (2 * x - 1) + delta_h);
    
    Q = fsolve(Obj_func, 0.001, options);
    
    disorder_v = (1 - Q) * delta_v + wv * Q * (1 - Q);
    
    molar_volume = molar_volume + disorder_v;
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
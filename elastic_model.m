% elastic model
function mixed_objective = elastic_model(P_incl,P_entrapment,T_entrapment,inclusion_abv,host_abv,host_mole_fraction,inclusion_mole_fraction,hp_dataset)

%initialize ambient (or measurement) conditions
P_amb = 0.001;               % kbar (1 bar)
T_amb = 25 + 273.15;         % K  

options = optimoptions(@fsolve, 'display','off', 'tolfun', 1e-10); 

% solve for Pfoot in the function "elastic_model_P_foot
Obj_func = @(P_foot) elastic_model_P_foot(P_incl,P_foot,inclusion_abv,host_abv,host_mole_fraction,inclusion_mole_fraction,hp_dataset); 
P_footy = fsolve(Obj_func, 0.001, options); % 0.001 = start pressure value in kbars

% initialize inclusion variables
number_of_inclusion_phases = length(inclusion_mole_fraction);
molar_volume_inclusion_T_amb_P_incl = 0;
molar_volume_inclusion_T_entrapment_P_entrapment = 0;
molar_volume_inclusion_T_amb_P_footy = 0;

% calculate molar volumes
for i = 1:1:number_of_inclusion_phases
    molar_volume_inclusion_T_amb_P_incl = molar_volume_inclusion_T_amb_P_incl + eos(P_incl,T_amb,hp_dataset,inclusion_abv{i}) * inclusion_mole_fraction{i};
    molar_volume_inclusion_T_entrapment_P_entrapment = molar_volume_inclusion_T_entrapment_P_entrapment + eos(P_entrapment,T_entrapment,hp_dataset,inclusion_abv{i}) * inclusion_mole_fraction{i};
    molar_volume_inclusion_T_amb_P_footy = molar_volume_inclusion_T_amb_P_footy + eos(P_footy,T_amb,hp_dataset,inclusion_abv{i}) * inclusion_mole_fraction{i};
end

% initialize host variables
number_of_host_phases = length(host_mole_fraction);
molar_volume_host_ambient = 0;
molar_volume_host_T_entrapment_P_entrapment = 0;
shear_modulus = 0;

% calculate the host shear modulus
for i = 1:1:number_of_host_phases
      shear_modulus_check = isnan(hp_dataset{host_abv{i},'shear_modulus'});
      poisson_ratio_check = isnan(hp_dataset{host_abv{i},'poisson_ratio'});
      
      % check if user has input a shear modulus, if so, use input shear modulus
      if shear_modulus_check == 0
          molar_volume_host_ambient = molar_volume_host_ambient + eos(P_amb,T_amb,hp_dataset,host_abv{i}) * host_mole_fraction{i};
          molar_volume_host_T_entrapment_P_entrapment = molar_volume_host_T_entrapment_P_entrapment + eos(P_entrapment,T_entrapment,hp_dataset,host_abv{i}) * host_mole_fraction{i};
          shear_modulus_phase{i} = hp_dataset{host_abv{i},'shear_modulus'};
          shear_modulus = shear_modulus + shear_modulus_phase{i} * host_mole_fraction{i}; 
      
      % if no shear modulus has been input, check if user has input a poisson ratio, if so, use calculated shear modulus from poisson ratio and bulk modulus   
      else  
          if poisson_ratio_check == 0
                molar_volume_host_ambient = molar_volume_host_ambient + eos(P_amb,T_amb,hp_dataset,host_abv{i}) * host_mole_fraction{i};
                molar_volume_host_T_entrapment_P_entrapment = molar_volume_host_T_entrapment_P_entrapment + eos(P_entrapment,T_entrapment,hp_dataset,host_abv{i}) * host_mole_fraction{i};
                k0{i} = hp_dataset{host_abv{i},'k0'};
                poisson_ratio{i} = hp_dataset{host_abv{i},'poisson_ratio'};
                shear_modulus = shear_modulus + (3 * k0{i} * (1 - 2 * poisson_ratio{i}))/(2 * (1 + poisson_ratio{i})) * host_mole_fraction{i};
                
          elseif poisson_ratio_check == 1
                return
          end
          
      end
      
end

%output elastic model solution
mixed_objective = (molar_volume_inclusion_T_amb_P_incl./molar_volume_inclusion_T_entrapment_P_entrapment - ...
    molar_volume_host_ambient./molar_volume_host_T_entrapment_P_entrapment) * (molar_volume_inclusion_T_entrapment_P_entrapment./ ...
    molar_volume_inclusion_T_amb_P_footy) - (0.75/shear_modulus) * (P_incl - P_amb);






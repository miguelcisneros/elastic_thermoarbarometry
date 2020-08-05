clc
close all
clear all

%Read Tables
hp_dataset = readtable('thermodynamic_properties.xlsx','ReadRowNames',true);
data = readtable('calculations_input_entrapment_P.xlsx','ReadRowNames',true);

%Intialize Variables
P_entrapment_GPa = zeros(length(data{:,{'mean_stress_MPa'}}),1);
T_entrapment_celsius_2 = zeros(length(data{:,{'mean_stress_MPa'}}),1);

%ENTER YOUR PARAMETERS  BELOW HERE
%Set Inclusion and Mole Fraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inclusion_abv{1} = 'Quartz';
inclusion_mole_fraction{1} = 1;

%Second Inclusion Example
% inclusion_abv{1} = 'Apatite_Cl';
% inclusion_mole_fraction{1} = 0.7;
% inclusion_abv{2} = 'Apatite_F';
% inclusion_mole_fraction{2} = 0.2;
% inclusion_abv{3} = 'Apatite_OH';
% inclusion_mole_fraction{3} = 0.1;

%Set Host and Mole Fraction
host_abv{1} = 'Garnet_Almandine';
host_mole_fraction{1} = 0.7;
host_abv{2} = 'Garnet_Grossular';
host_mole_fraction{2} = 0.2;
host_abv{3} = 'Garnet_Pyrope';
host_mole_fraction{3} = 0.1;

%Second Host Example
% host_abv{1} = 'Epidote';
% host_mole_fraction{1} = 0.31;
% host_abv{2} = 'Clinozoisite';
% host_mole_fraction{2} = 0.69;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:length(data{:,{'mean_stress_MPa'}})

    raman_data = data{i,{'mean_stress_MPa','T_entrapment_celsius'}};
    
    mean_stress_MPa = raman_data(1,1);
    T_entrapment_celsius = raman_data(1,2);
    
    if mean_stress_MPa == 0
        
        P_entrapment_GPa(i,1) = 0;
        T_entrapment_celsius_2(i,1) = T_entrapment - 273.15;
    
    else
    
    P_incl = mean_stress_MPa/100;
    
    T_entrapment = T_entrapment_celsius + 273.15; %%celsius;
    
    options = optimoptions(@fsolve, 'tolfun', 1e-10); 
    
    Obj_func = @(P_entrapment) elastic_model(P_incl,P_entrapment,T_entrapment,inclusion_abv,host_abv,host_mole_fraction,inclusion_mole_fraction,hp_dataset);
    P_entrapment = fsolve(Obj_func, 0.001, options); % 0.001 = start pressure value in kbars
    fprintf('Entrapment P = %.2f GPa\n', P_entrapment/10);
    
    P_entrapment_GPa(i,1) = P_entrapment/10;
    T_entrapment_celsius_2(i,1) = T_entrapment - 273.15;
    
    end
    
end


%Creating Output table
P_entrapment_GPa = table(P_entrapment_GPa);
T_entrapment_celsius = table(T_entrapment_celsius_2);

new_table = [T_entrapment_celsius P_entrapment_GPa];
filename = 'calculations_output_entrapment_P.xlsx';
writetable(new_table,'calculations_output_entrapment_P.xlsx');

    
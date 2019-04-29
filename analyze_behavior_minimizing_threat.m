%% analyze_behavior_minimizing_threat
% This script performs the main behavioral analyses for the article 
%   "Minimizing threat via heuristic and optimal policies recruits 
%   hippocampus and medial prefrontal cortex" by Korn & Bach; 
%   for more information contact: christoph.w.korn@gmail.com
% subfunctions
% - BIC_fixed_plot_v1
% - PXP_random_plot_v1 (optional; requires SPM)
% - sigmoid
% - sigmoid_derive
% - fisher_z_calculate
% - fisher_z_calculate_back

clear
close all

%% load data
load('Behavioral_data_minimizing_threat.mat')
w_sample = 1; % w is for Which
if w_sample == 1
    USE = fMRI_sample;
elseif w_sample == 2
    USE = behavioral_sample;
end
n_sub = size( USE, 2 );


%% load variables related to (pseudo-)optimal policies
% consideration of predation AND starvation or of only one of the two
% J1: predation AND starvation: This is the true optimal policy
J1 = load('Optimal_policy_predation_starvation.mat');
% J2: ONLY starvation
J2 = load('Optimal_policy_starvation.mat');
% J3: ONLY predation
J3 = load('Optimal_policy_predation.mat');


%% delete trials in w participants "died" or did not press a key
% T is a structure for Trial-related variables
T.SUB_size_header = {'all ','alive','not_answered','final'};
for i_sub = 1:n_sub
       
    T.SUB_size( i_sub, 1 )  = size( USE{ i_sub }, 1 );
    T.w_alive            = USE{ i_sub }(:,5) > 0; 
    T.SUB_size( i_sub, 2 )  = sum( T.w_alive, 1);
    USE_1{ i_sub }          = USE{ i_sub }( T.w_alive, :);
    T.w_answered         = ~isnan( USE_1{ i_sub }(:,10) ); 
    T.SUB_size( i_sub, 3 )  = T.SUB_size( i_sub, 2 ) - sum( T.w_answered, 1);
    USE_2{ i_sub }          = USE_1{ i_sub }( T.w_answered, :); 
    % USE_2 contains the relevant data for further analyses    
    
    % here is a simple possibility to exclude the 1st or 2nd half of the runs
    T.w_gamles_not = []; % [1:5]; % [6,10];
    USE_2{ i_sub }( ismember( USE_2{ i_sub }(:,18), T.w_gamles_not ), : ) = [];    
    T.SUB_size( i_sub, 4 )  = size( USE_2{ i_sub }, 1 );   
   
end
T.SUB_size_mean(1,:) = mean( T.SUB_size );
T.SUB_size_mean(2,:) = std( T.SUB_size );


%% get overall choice percentages
% to check variability of choices
for i_sub = 1:n_sub     
    
    T.choice_all( i_sub, 1 ) = mean( USE_2{ i_sub }(:,10) );
    
end
T.choice_mean(1,1) = mean(T.choice_all);
T.choice_mean(1,2) = std(T.choice_all);
T.subs_below = T.choice_all < 0.85;
T.subs_above = T.choice_all > 0.15;


%% transform p_predator in USE_2
% the experimental script recorded 1_minus_p_predator
% i.e., the probability of not being attacked by the predator
% this is here transformed to p_predator for consistency
for i_sub = 1:n_sub    
         
    USE_2{ i_sub }(:,16) = 1 - USE_2{ i_sub }(:,16);
    
end
header_for_columns{1,16} = '16_p_predator';


%% add gain magnitudes to USE_2 
% these were not saved in experimental script
% but can be recovered from the variables used to calculate the optimal
% policy via the "3_index_forests"
for i_sub = 1:n_sub
    
    T.w_states = USE_2{ i_sub }(:, 3);
    T_gainss = J1.magn_gainss_take( T.w_states, :);
    % add gain magnitudes according to "6_weather_type"
    for i_f = 1:size( T_gainss, 1)
        USE_2{ i_sub }(i_f, 19 ) = T_gainss(i_f, USE_2{ i_sub }(i_f, 6));       
    end
    clear T_gainss
    
end
header_for_columns{1,19} = '19_gain_magnitude';

            
%% add (pseudo-)optimal policies to USE_2
% T.horiz_count contains the number of forests with the respective time
% horizon; this is used for the variable "pseudo-optimal: horizon-2.5"
T.horiz_count = [1,24; 2, 9; 3, 4; 4, 2; 5, 1; ];

% loop over all subjects
for i_sub = 1:n_sub
    
    % go through each trial individually
    for i_tri = 1:size( USE_2{ i_sub }, 1 ) 
        
        
        %% relevant variables to index value differences according to optimal policies
        T.w_gamble  = USE_2{ i_sub }(i_tri, 3);        
        T.w_energy  = USE_2{ i_sub }( i_tri, 5) + 1 + ( USE_2{ i_sub }( i_tri, 6)-1 ) * 6;
        % remaining time horizons
        % 6 - "4_order_days_in_forest" (i.e., days past)
        T.w_horiz_5 = 6 - USE_2{ i_sub }(i_tri, 4); % 5 days left
        T.w_horiz_4 = T.w_horiz_5 - 1;          % 4 days left
        T.w_horiz_3 = T.w_horiz_5 - 2;          % 3 days left    
        T.w_horiz_2 = T.w_horiz_5 - 3;          % 2 days left
        T.w_horiz_1 = T.w_horiz_5 - 4;          % 1 day  left
        % cap time horizons at 1
        if T.w_horiz_4 < 1; T.w_horiz_4 = 1; end
        if T.w_horiz_3 < 1; T.w_horiz_3 = 1; end
        if T.w_horiz_2 < 1; T.w_horiz_2 = 1; end
        if T.w_horiz_1 < 1; T.w_horiz_1 = 1; end
        
        
        %% actual addition of the (pseudo-)optimal policies to USE_2
        % rew_for REWard of FORaging
        % rew_wai REWard of WAIting
        % J1: predation AND starvation: This is the true optimal policy
        USE_2{ i_sub }(i_tri,20) = J1.rew_for( T.w_horiz_5 + 1, T.w_energy, T.w_gamble ) - J1.rew_wai( T.w_horiz_5 + 1, T.w_energy, T.w_gamble ); % normative horizon
        USE_2{ i_sub }(i_tri,21) = J1.rew_for( 1               + 1, T.w_energy, T.w_gamble ) - J1.rew_wai( 1               + 1, T.w_energy, T.w_gamble ); % horizon of 1 step
        % J2: ONLY starvation
        USE_2{ i_sub }(i_tri,22) = J2.rew_for( T.w_horiz_5 + 1, T.w_energy, T.w_gamble ) - J2.rew_wai( T.w_horiz_5 + 1, T.w_energy, T.w_gamble ); % normative horizon
        % J3: ONLY predation
        USE_2{ i_sub }(i_tri,23) = J3.rew_for( T.w_horiz_5 + 1, T.w_energy, T.w_gamble ) - J3.rew_wai( T.w_horiz_5 + 1, T.w_energy, T.w_gamble ); % normative horizon
        
        
        %% calculations for "pseudo-optimal: horizon-2.5"
        T.hori_pol{ i_sub }(i_tri,1) = J1.rew_for( T.w_horiz_4 + 1, T.w_energy, T.w_gamble ) - J1.rew_wai( T.w_horiz_4 + 1, T.w_energy, T.w_gamble );
        T.hori_pol{ i_sub }(i_tri,2) = J1.rew_for( T.w_horiz_3 + 1, T.w_energy, T.w_gamble ) - J1.rew_wai( T.w_horiz_3 + 1, T.w_energy, T.w_gamble );
        T.hori_pol{ i_sub }(i_tri,3) = J1.rew_for( T.w_horiz_2 + 1, T.w_energy, T.w_gamble ) - J1.rew_wai( T.w_horiz_2 + 1, T.w_energy, T.w_gamble );
        T.hori_pol{ i_sub }(i_tri,4) = J1.rew_for( T.w_horiz_1 + 1, T.w_energy, T.w_gamble ) - J1.rew_wai( T.w_horiz_1 + 1, T.w_energy, T.w_gamble );
        
        USE_2{ i_sub }(i_tri,24) = ...
          ( T.horiz_count(1,2) * T.hori_pol{ i_sub }(i_tri,4) + ... 
            T.horiz_count(2,2) * T.hori_pol{ i_sub }(i_tri,3) + ... 
            T.horiz_count(3,2) * T.hori_pol{ i_sub }(i_tri,2) + ... 
            T.horiz_count(4,2) * T.hori_pol{ i_sub }(i_tri,1) + ... 
            T.horiz_count(5,2) * USE_2{ i_sub }(i_tri,20) )/40;
            % division by 40 because there are 40 trials (aka days) per run
            
    end    
end
header_for_columns{1,20} = '20_optimal_policy';
header_for_columns{1,21} = '21_pseudo-optimal_horizon-1';
header_for_columns{1,22} = '22_pseudo-optimal_starvation-only';
header_for_columns{1,23} = '23_pseudo-optimal_predation-only';
header_for_columns{1,24} = '24_pseudo-optimal_horizon-2.5';


%% add other variables to USE_2
for i_sub = 1:n_sub

    
    %% binary energy
    USE_2{ i_sub }(:,25) = USE_2{ i_sub }(:,5) == 1;    
   
    
    %% expected energy & expected energy change
    % expected energy       = p_pred*0 + (1-p_pred)*( foraging_outcome )
    % foraging_outcome      = p_gain*foraging outcome_win + (1-p_gain)*foraging_outcome_loss 
    % foraging_outcome_win  = current_energy + m_gain
    % foraging_outcome_loss = current_energy - 2    
    
    % foraging_outcome_win
    USE_2{ i_sub }(:,26) =  USE_2{ i_sub }(:,5) + USE_2{ i_sub }(:,19); 
    USE_2{ i_sub }( USE_2{ i_sub }(:,26) > 5 ,26) = 5;
   
    % foraging_outcome_loss
    USE_2{ i_sub }(:,27) =  USE_2{ i_sub }(:,5) - 2; 
    USE_2{ i_sub }( USE_2{ i_sub }(:,27) < 0 ,27) = 0;
    
    % foraging_outcome
    USE_2{ i_sub }(:,28) = USE_2{ i_sub }(:,14) .* USE_2{ i_sub }(:,26) + ...
        (1-USE_2{ i_sub }(:,14)) .* USE_2{ i_sub }(:,27);
    % expected energy
    USE_2{ i_sub }(:,29) = (1-USE_2{ i_sub }(:,16)) .* USE_2{ i_sub }(:,28);
    % expected energy change = difference to current energy
    USE_2{ i_sub }(:,30) = USE_2{ i_sub }(:,29) - USE_2{ i_sub }(:,5); % v32: put as EV gained
    
    
    %% variables related to past trials or forests   
    % past energy change (has to be shifted by one trial)
    USE_2{ i_sub }(2:end,31) = USE_2{ i_sub }(1:end-1,11) - USE_2{ i_sub }(1:end-1,5);
    % win-stay-lose-shift (WSLS)
    USE_2{ i_sub }(:,32)     = USE_2{ i_sub }(:,31) >= 0;
    
    % death in past forest
    T.lfo_unique = unique( USE_2{ i_sub }(:,3) ); % lfo is for Last FOrest
    n_lfo = size( T.lfo_unique, 1 );
    for i_lfo = 1:n_lfo 
        T.lfo_find = find( USE_2{ i_sub }(:,3) == T.lfo_unique(i_lfo) );
        % death in previous forest ?
        T.lfo_prev = T.lfo_find(1,1)-1;
        if T.lfo_prev == 0
            T.lfo_dead = 0;
        else
            if USE_2{ i_sub }(T.lfo_prev,11) == 0
                T.lfo_dead = 1;
            else
                T.lfo_dead = 0;
            end
        end
        USE_2{ i_sub }(T.lfo_find,33) = T.lfo_dead;
    end 
    % outcome of current trial (same as past energy change but not shifted)
    USE_2{ i_sub }(:,34) = USE_2{ i_sub }(:,11) - USE_2{ i_sub }(:,5);
    
    %% log RT
    USE_2{ i_sub }(:,35) = log( USE_2{ i_sub }(:,9) );

end
header_for_columns{1,25} = '25_binary_energy';
header_for_columns{1,26} = '26_foraging_outcome_win';
header_for_columns{1,27} = '27_foraging_outcome_loss';
header_for_columns{1,28} = '28_foraging_outcome';
header_for_columns{1,29} = '29_expected_energy';
header_for_columns{1,30} = '30_expected_energy_30';
header_for_columns{1,31} = '31_past energy change';
header_for_columns{1,32} = '32_WSLS';
header_for_columns{1,33} = '33_death_in_past_forest';
header_for_columns{1,34} = '34_outcome_current_trial';
header_for_columns{1,35} = '35_log_RT';


%% FIRST model comparison 
% do simple regression loop for the relevant variables    
mod_indece_1 = [20,16,14,19,5,25,29,30,4,21,22,23,24,31,32,33];
mod_header_1 = header_for_columns( mod_indece_1 );
for i_sub = 1:n_sub    
    
    w_main      = 1:T.SUB_size( i_sub, 4 ); % select gambles
    y_regress   = USE_2{ i_sub }(w_main,10) + 1; % specify y 
    n_y         = size(y_regress,1);
    
    for i_mod = 1:length( mod_indece_1 )
        
        x_regress   = USE_2{ i_sub }( w_main, mod_indece_1(i_mod) );
        n_x         = size(x_regress,2)+1; % +1 because of intercept
        
        [SUB_B_1{i_mod}(i_sub, :), SUB_reg_dev_1(i_sub, i_mod), SUB_p_1{i_mod}(i_sub, :)] = mnrfit(x_regress, y_regress);
        BIC_1(i_sub, i_mod) = SUB_reg_dev_1(i_sub, i_mod) + n_x*log( n_y ); % Bayesian Information Criterion    
    
    end    
end


%% SECOND model comparison
mod_header_2 = mod_header_1( [1,3:end] );
mod_indece_2 = mod_indece_1( [1,3:end] );
for i_sub = 1:n_sub
    
    w_main      = 1:T.SUB_size( i_sub, 4 ); % select gambles
    y_regress   = USE_2{ i_sub }(w_main,10) + 1; % specify y 
    n_y         = size(y_regress,1);
    
    for i_mod = 1:length( mod_indece_2 )
        
        x_regress   = [ USE_2{ i_sub }( w_main, 16 ), USE_2{ i_sub }( w_main, mod_indece_2(i_mod) ) ];
        n_x         = size(x_regress,2)+1; % +1 because of intercept

        [SUB_B_2{i_mod}(i_sub, :), SUB_reg_dev_2(i_sub, i_mod), SUB_p_2{i_mod}(i_sub, :)] = mnrfit(x_regress, y_regress);      
        BIC_2(i_sub, i_mod) = SUB_reg_dev_2(i_sub, i_mod) + n_x*log( n_y ); % Bayesian Information Criterion
                
    end  
end


%% THIRD model comparison
mod_header_3 = mod_header_1( 3:end );
mod_indece_3 = mod_indece_1( 3:end );
for i_sub = 1:n_sub
    
    w_main      = 1:T.SUB_size( i_sub, 4 ); % select gambles
    y_regress   = USE_2{ i_sub }(w_main,10) + 1; % specify y 
    n_y         = size(y_regress,1);
    
    for i_mod = 1:length( mod_indece_3 )        

        x_regress   = [ USE_2{ i_sub }( w_main, 16 ), USE_2{ i_sub }( w_main, 20 ), USE_2{ i_sub }( w_main, mod_indece_3(i_mod) ) ];    
        n_x         = size(x_regress,2)+1; % +1 because of intercept

        [SUB_B_3{i_mod}(i_sub, :), SUB_reg_dev_3(i_sub, i_mod), SUB_p_3{i_mod}(i_sub, :)] = mnrfit(x_regress, y_regress);         
        BIC_3(i_sub, i_mod) = SUB_reg_dev_3(i_sub, i_mod) + n_x*log( n_y ); % Bayesian Information Criterion
        
    end
end


%% sumary statistics BIC 
BIC_sum_1 = sum( BIC_1 );
BIC_sum_2 = sum( BIC_2 );
BIC_sum_3 = sum( BIC_3 );
BIC_sum_1 = BIC_sum_1 - BIC_sum_1(1,1);
BIC_sum_2 = BIC_sum_2 - BIC_sum_2(1,1);
BIC_sum_3 = BIC_sum_3 - BIC_sum_3(1,1);
BIC_comp_1_2 = [BIC_1(:,2), BIC_2(:,1)];
BIC_sum_comp_1_2 = sum( BIC_comp_1_2 );
BIC_sum_comp_1_2 = BIC_sum_comp_1_2 - BIC_sum_comp_1_2(1,1);


%% model comparison plots
BIC_fixed_plot_v1( BIC_sum_1, 'FIRST models' )
PXP_1_H = PXP_random_plot_v1( -BIC_1/2, 'FIRST models' );

BIC_fixed_plot_v1( BIC_sum_2, 'SECOND models' )
PXP_2_H = PXP_random_plot_v1( -BIC_2/2, 'SECOND models' );

BIC_fixed_plot_v1( BIC_sum_3, 'THIRD models' )
PXP_3_H = PXP_random_plot_v1( -BIC_3/2, 'THIRD models' );

BIC_fixed_plot_v1( BIC_sum_comp_1_2, 'FIRST vs SECOND models' )
PXP_1_2_H = PXP_random_plot_v1( -BIC_comp_1_2/2, 'FIRST vs SECOND models' );


%% INTERACTION model comparison
for i_sub = 1:n_sub
    
    w_main      = 1:T.SUB_size( i_sub, 4 ); % select gambles
    y_regress   = USE_2{ i_sub }(w_main,10) + 1; % specify y
    n_y         = size(y_regress,1);
    
    % here slightly different logic for setting up x_regress: manually build larger x_regress and loop over it    
    x_regress_I{1} = [ USE_2{ i_sub }(w_main,16), USE_2{ i_sub }(w_main,14), USE_2{ i_sub }(w_main,16) .* USE_2{ i_sub }(w_main,14) ];    
    x_regress_I{2} = [ USE_2{ i_sub }(w_main,16), USE_2{ i_sub }(w_main,19), USE_2{ i_sub }(w_main,16) .* USE_2{ i_sub }(w_main,19) ];
    x_regress_I{3} = [ USE_2{ i_sub }(w_main,16), USE_2{ i_sub }(w_main, 5), USE_2{ i_sub }(w_main,16) .* USE_2{ i_sub }(w_main, 5) ];
    x_regress_I{4} = [ USE_2{ i_sub }(w_main,16), USE_2{ i_sub }(w_main,25), USE_2{ i_sub }(w_main,16) .* USE_2{ i_sub }(w_main,25) ];
    x_regress_I{5} = [ USE_2{ i_sub }(w_main,16), USE_2{ i_sub }(w_main, 4), USE_2{ i_sub }(w_main,16) .* USE_2{ i_sub }(w_main, 4) ];
    x_regress_I{6} = [ USE_2{ i_sub }(w_main,14), USE_2{ i_sub }(w_main,19), USE_2{ i_sub }(w_main,14) .* USE_2{ i_sub }(w_main,19) ];
    x_regress_I{7} = [ USE_2{ i_sub }(w_main,14), USE_2{ i_sub }(w_main, 5), USE_2{ i_sub }(w_main,14) .* USE_2{ i_sub }(w_main, 5) ];
    x_regress_I{8} = [ USE_2{ i_sub }(w_main,14), USE_2{ i_sub }(w_main,25), USE_2{ i_sub }(w_main,14) .* USE_2{ i_sub }(w_main,25) ];
    x_regress_I{9} = [ USE_2{ i_sub }(w_main,14), USE_2{ i_sub }(w_main, 4), USE_2{ i_sub }(w_main,14) .* USE_2{ i_sub }(w_main, 4) ];

    for i_mod = 1:length( x_regress_I );
        
        n_x         = size(x_regress,2)+1; % +1 because of intercept
        [SUB_B_I{i_mod}(i_sub, :), SUB_reg_dev_I(i_sub, i_mod), SUB_p_I{i_mod}(i_sub, :)] = mnrfit(x_regress_I{i_mod}, y_regress);        
        BIC_I(i_sub, i_mod) = SUB_reg_dev_I(i_sub, i_mod) + n_x*log( n_y ); % Bayesian Information Criterion
        
    end
end


%% BIC_I model comparisons & plots
BIC_I_with_win = [ BIC_2(:,1), BIC_I  ];
BIC_I_with_win_sum = sum( BIC_I_with_win, 1 );
BIC_I_with_win_sum = BIC_I_with_win_sum - BIC_I_with_win_sum(1,1);

BIC_fixed_plot_v1( BIC_I_with_win_sum, 'INTERACTION models' )
PXP_Inter_H = PXP_random_plot_v1( -BIC_I_with_win/2, 'INTERACTION models' );


%% mean parameters for all one- & two-variable models 
% these are needed for calculating uncertainties & discrepancies
% here they are calculated within the respective sample
% but they can also be taken from the behavioral sample for the fMRI sample
for i_mod = 1:size( SUB_B_1, 2 )
    
    SUB_B_1_mean(i_mod,:) = mean( SUB_B_1{i_mod} ); 
    
end
for i_mod = 1:size( SUB_B_2, 2 )
    
    SUB_B_2_mean(i_mod,:) = mean( SUB_B_2{i_mod} ); 
    
end


%% calculate uncertainties & discrepancies & add these to USE_2
for i_sub = 1:n_sub
       
   % take mean model fits
   T.fit_pol{ i_sub }(:, 1 )    = SUB_B_1_mean(2,1) + SUB_B_1_mean(2,2) * USE_2{ i_sub }(:,16); % One-parameter: 16_p_predator 
   T.fit_pol{ i_sub }(:, 2 )    = SUB_B_1_mean(1,1) + SUB_B_1_mean(1,2) * USE_2{ i_sub }(:,20); % One-parameter: 20_optimal_policy 
   T.fit_pol{ i_sub }(:, 3 )    = SUB_B_2_mean(1,1) + SUB_B_2_mean(1,2) * USE_2{ i_sub }(:,16) + SUB_B_2_mean(1,3) * USE_2{ i_sub }(:,20); % Two-parameter: 16_p_predator + 20_optimal_policy
   % convert these  
   T.fit_pol{ i_sub }(:, 4:6 )  = sigmoid( -T.fit_pol{ i_sub }(:, 1:3 ) ); 
   % uncertainties
   T.fit_pol{ i_sub }(:, 7:9 )  = sigmoid_derive( -T.fit_pol{ i_sub }(:, 1:3 ) );    
   % discrepancies
   T.fit_pol{ i_sub }(:, 10 )   = T.fit_pol{ i_sub }(:, 4 ) - T.fit_pol{ i_sub }(:, 5 ); % One-parameter: 16_p_predator & 20_optimal_policy
   T.fit_pol{ i_sub }(:, 11 )   = abs( T.fit_pol{ i_sub }(:, 10 ) );
   % transfer relevant variables to USE_2 
   USE_2{ i_sub }(:, 36:38 )    = T.fit_pol{ i_sub }(:, [7,8,11]);
   
end
header_for_columns{1,36} = '36_uncertainty_p_predator';
header_for_columns{1,37} = '37_uncertainty_optimal_policy';
header_for_columns{1,38} = '38_discrepancy_choice_probabilities';


%% do correlations
% ALSO: using fisher transformations to look at correlations
cor_mod = mod_indece_1;
cor_mri = [16,20,36:38,35,34];

for i_sub = 1:n_sub
    
    [C.cor_mod_r(:,:,i_sub), C.cor_mod_p(:,:,i_sub)] = corr( USE_2{ i_sub }(:,cor_mod) );   
    C.cor_mod_f(:,:,i_sub) = fisher_z_calculate( C.cor_mod_r(:,:,i_sub) );   
    
    [C.cor_mri_r(:,:,i_sub), C.cor_mri_p(:,:,i_sub)] = corr( USE_2{ i_sub }(:,cor_mri) );   
    C.cor_mri_f(:,:,i_sub) = fisher_z_calculate( C.cor_mri_r(:,:,i_sub) );   
   
end

C.cor_mod_f_mean       = mean( C.cor_mod_f, 3 );
C.cor_mod_f_mean_back  = fisher_z_calculate_back( C.cor_mod_f_mean );
C.cor_mod_f_mean_sqar  = C.cor_mod_f_mean_back .^2;

C.cor_mri_f_mean       = mean( C.cor_mri_f, 3 );
C.cor_mri_f_mean_back  = fisher_z_calculate_back( C.cor_mri_f_mean );
C.cor_mri_f_mean_sqar  = C.cor_mri_f_mean_back .^2;


%% get data for LME of RTs
LME.index_1 = [10,9,35,16,20,36,37,38]; 
LME.header_for_columns_RT = { 'id_sub', header_for_columns( LME.index_1 ) };
LME.total = [];

for i_sub = 1:n_sub
    
    LME.sub = [];
    LME.sub = [ zeros( size(USE_2{ i_sub }, 1), 1) + i_sub, ...        
        USE_2{ i_sub }(:, LME.index_1 ) ];
    LME.total = [LME.total; LME.sub];

end
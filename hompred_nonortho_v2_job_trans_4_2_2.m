%% new getonsets with ortho option
clear all
cwd = pwd;


%% Change according to your directory strucutre and scan parameters
fs                  = filesep;
dir_base            = 'D:\STUDIES_Zurich\HOMPRED_fMRI';
dir_analysis_base   = 'ANALYSIS';
dir_fMRI_base       = 'DATA_FMRI';
dir_epi             = 'functional_1';
% dir_multiple        = 'scanphys_3_multiple_pulse'; % 'scanphys_3_multiple_eeg'; % 

number_of_analysis_behaviour    = '4';
number_of_analysis_scanner      = '2';
number_of_analysis_model        = '2';

dir_analysis_behviour           = ['beh_v' number_of_analysis_behaviour];
dir_analysis_scanner            = ['stats_b_' number_of_analysis_behaviour '_s_' number_of_analysis_scanner '_m_' number_of_analysis_model];
cd([dir_base fs dir_analysis_base]);
mkdir(dir_analysis_scanner);


%% subjects
sub      = 240:241; % [233:236,238:241]; % [214:221,223:226,228,230:236,238:241]; 
estimate = 1;
ortho_yn = 0;


%% spm
spm fmri


for n = 1:length(sub)
    
    sub_string = num2str(sub(n));
       
    cd([dir_base fs dir_analysis_base fs dir_analysis_behviour]);
    filename_sub = ['beh_v' number_of_analysis_behaviour '_sub_' num2str(sub(n)) ];
    clear  'Z'
    load(filename_sub, 'Z');
    cd ..
    
    blockNo    = size( Z.which_ses, 2); % i.e. sessions
    
        
    cd(dir_analysis_scanner);
    dir_analysis_scanner_sub            = ['stats_b_' number_of_analysis_behaviour '_s_' number_of_analysis_scanner '_m_' number_of_analysis_model '_s_' sub_string];
    mkdir(dir_analysis_scanner_sub);
    
    outputDir = [dir_base fs dir_analysis_base fs dir_analysis_scanner fs dir_analysis_scanner_sub];
    
    %-----------------------------------------------------------------------
    % Job saved on 11-Mar-2016 11:20:25 by cfg_util (rev $Rev: 6134 $)
    % spm SPM - SPM12 (6225)
    % cfg_basicio BasicIO - Unknown
    %-----------------------------------------------------------------------
    matlabbatch{1}.spm.stats.fmri_spec.dir              = {outputDir}; % '<UNDEFINED>';
    matlabbatch{1}.spm.stats.fmri_spec.timing.units     = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT        = 2.1;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t    = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0   = 1;
    
    for k_i = 1:blockNo
        
        k = Z.which_ses( 1, k_i ); % this is new distinction between k and k_i
        
        %% scans
        dir_fMRI_sub        = ['sub_' sub_string ];
        dir_sess_specified  = ['session_' num2str(k)];
        epiDir              = [dir_base fs dir_fMRI_base fs dir_fMRI_sub fs dir_epi fs dir_sess_specified];
        f                   = spm_select('List', epiDir, '^swaumr.*\.img$');     % Select smoothed normalised images
        files               = cellstr([repmat([epiDir fs],size(f,1),1) f]);
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).scans                    = files; % '<UNDEFINED>';
        
        f = []; files = [];
        
        %% conditions
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).name             = 'forest_1';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).onset            = Z.forest{k}(:,1)'/1000; % '<UNDEFINED>';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).duration         = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).tmod             = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).pmod(1).name     = 'forest_2';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).pmod(1).param    = Z.forest{k}(:,2)'; % '<UNDEFINED>';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).pmod(1).poly     = 1;
        
%         matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).pmod(2).name     = 'forest_3';
%         matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).pmod(2).param    = Z.forest{k}(:,3)'; % '<UNDEFINED>';
%         matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).pmod(2).poly     = 1;

        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(1).orth             = ortho_yn;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).name             = 'choice_1';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).onset            = Z.choice{k}(:,1)'/1000; % '<UNDEFINED>';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).duration         = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).tmod             = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(1).name     = 'choice_2';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(1).param    = Z.choice{k}(:,2)'; % p_gain
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(1).poly     = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(2).name     = 'choice_3';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(2).param    = Z.choice{k}(:,3)'; % p_threat 
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(2).poly     = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(3).name     = 'choice_4';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(3).param    = Z.choice{k}(:,8)'; % CONFLICT p_threat p_gain
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(3).poly     = 1;
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(4).name     = 'choice_5';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(4).param    = Z.choice{k}(:,10)'; % optimal  
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(4).poly     = 1;
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(5).name     = 'choice_6';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(5).param    = Z.choice{k}(:,7); % CONFLICT threat & optimal
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(5).poly     = 1;
%         
%         matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(6).name     = 'choice_7';
%         matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(6).param    = abs(Z.choice{k}(:,12))'; % difficulty horizon_4
%         matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).pmod(6).poly     = 1;
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(2).orth             = ortho_yn;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(3).name             = 'outcome_1';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(3).onset            = Z.outcome{k}(:,1)'/1000; % '<UNDEFINED>';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(3).duration         = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(3).tmod             = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(3).pmod.name        = 'outcome_2';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(3).pmod.param       = Z.outcome{k}(:,2)'; % diff '<UNDEFINED>';
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(3).pmod.poly        = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).cond(3).orth             = ortho_yn;
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).multi                    = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).regress                  = struct('name', {}, 'val', {});
        
        
        %% motion, heart, breathing        
        motion_file     = spm_select('List', epiDir, '^.*\.txt$');
        multiregFile    = sprintf('regs%d_%d.mat',sub(n),k);        
        R               = textread([epiDir fs motion_file]);        
        multiregPath    = [outputDir fs multiregFile];
        save(multiregPath, 'R');        
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).multi_reg                = {multiregPath}; % {''}; %  
        R_read = []; R = [];
        
        
        %% hpf
        matlabbatch{1}.spm.stats.fmri_spec.sess(k_i).hpf                      = 128;
        
        
    end % block number
    
    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    
    %% run batch
    cfg_util('initjob', matlabbatch);
    cfg_util('run', matlabbatch);
    cfg_util('deljob', matlabbatch);
    
    
    if estimate
        clear matlabbatch
        outputDir_SPM = [outputDir fs 'SPM.mat'];
        %-----------------------------------------------------------------------
        % Job saved on 11-Mar-2016 12:14:08 by cfg_util (rev $Rev: 6134 $)
        % spm SPM - SPM12 (6225)
        % cfg_basicio BasicIO - Unknown
        %-----------------------------------------------------------------------
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {outputDir_SPM}; % '<UNDEFINED>';
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    end
    
    
    %% run batch
    cfg_util('initjob', matlabbatch);
    cfg_util('run', matlabbatch);
    cfg_util('deljob', matlabbatch);
    
    clear matlabbatch
    
end % subject
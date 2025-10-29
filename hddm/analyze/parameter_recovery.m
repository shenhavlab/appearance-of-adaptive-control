%% Run parameter recovery

clear; clc;

addpath(genpath('../util/ScientificColourMaps7'))
vik = load('vik.mat')
vik = vik.vik;

roma = load('roma.mat')
roma = roma.roma;

berlin = load('berlin.mat')
berlin = berlin.berlin;

batlow = load('batlow.mat')
batlow = batlow.batlow;

vec = @(x) x(:);


%% recover from which model


model_rec1 = 'VDOV_collapse_subj'; %'VDev_subj','VD_static_subj','VD_collapse_subj','VDOV_collapse_subj',...
model_rec2 = 'VDOV_collapse_subj'; %'VDev_subj','VD_static_subj','VD_collapse_subj','VDOV_collapse_subj',...



%% Get simualted parameters (traces)

p_name = {...
    'z_trans', 't', ...
    'v_maxVD', 'v_minVD',...
    'a_Intercept', 'a_absMaxVD', 'a_absMinVD', 'a_maxOV', 'a_minOV',...
    'theta_Intercept', 'theta_absMaxVD', 'theta_absMinVD', 'theta_maxOV', 'theta_minOV'}


n_params = length(p_name);

wild = sprintf('gen-%s__fit-%s*', model_rec1, model_rec2)

fn = dir(sprintf('../recovery-local/traces/%s/%s', wild, wild))
n_files = 40


% get parameters
params = cell(n_files, n_params);


parfor ww = 1:n_files % for each filname

    fit_t = readtable(fullfile(fn(ww).folder, fn(ww).name));
   
    for pp = 1:n_params
        params{ww,pp} = fit_t.(p_name{pp});
    end

end

disp('loaded')




%% Get original parameters (traces)

orig = readtable(sprintf('../fit-results/rep-models/results/%s_traces.csv', model_rec2));
sel_fin = @(x) x(isfinite(x));
disp('loaded')


%% recover parameters (traces)

plot_method = 'cat'


f=figure;
set(f,'defaultTextInterpreter','none')

for nn = 1:n_params


    nexttile; hold on;

    % plot orig
    orig_trace = orig.(p_name{nn});
    [fo,xo,bw]=ksdensity(orig_trace, 'Kernel','normal', 'Bandwidth','normal-approx');
    plt_orig=plot(xo,fo, '-k', 'LineWidth', 3);
    xline(median(orig_trace), '-k', 'LineWidth', 2)

    % plot fit
    switch plot_method
        case 'seperate'
            for ww = 1:size(params,1)
                [ff,xf]=ksdensity(params{ww,nn}, 'Bandwidth',bw);
                plt_fit=plot(xf,ff, ':g', 'LineWidth', 3);
            end

        case 'mean'
            trace_cat = [];
            for ww = 1:size(params,1)
                trace_cat = [trace_cat,params{ww,nn}];
            end

            [ff,xf]=ksdensity(median(trace_cat), 'Bandwidth',bw);
            plt_fit=plot(xf,ff, ':g', 'LineWidth', 3);
            xline(median(median(trace_cat)), ':g', 'LineWidth', 2)

        case 'cat'
            trace_cat = [];
            for ww = 1:size(params,1)
                trace_cat = [trace_cat,params{ww,nn}];
            end
            
            % [ff,xf]=ksdensity(vec(trace_cat), 'Bandwidth',bw);
            % plot(xf,ff, '--g', 'LineWidth', 3);

            centered_trace_cat = vec(trace_cat-mean(trace_cat)) + mean(trace_cat(:));
            [ff,xf]=ksdensity(centered_trace_cat, 'Bandwidth',bw);
            plt_fit=plot(xf,ff, ':g', 'LineWidth', 3);
            xline(median(centered_trace_cat), ':g', 'LineWidth', 2)

    end
    

    set(gca, 'TickDir', 'out', 'LineWidth', 1)
    title(p_name(nn))
    xlabel('parameter')
    ylabel('density')




    if nn==1

        legend([plt_orig, plt_fit],'generating', 'fit')

    end

end







%% Get simulated parameters (summary)
% 
% p_name = {...
%     'z_trans', 't', ...
%     'v_maxVD', 'v_minVD',...
%     'a_Intercept', 'a_absMaxVD', 'a_absMinVD', 'a_maxOV', 'a_minOV',...
%     'theta_Intercept', 'theta_absMaxVD', 'theta_absMinVD', 'theta_maxOV', 'theta_minOV'}
% 
% 
% np = length(p_name);
% 
% wild = sprintf('gen-%s__fit-%s*', model_rec1, model_rec2)
% 
% fn = dir(sprintf('../recovery-local/summary/%s/%s', wild, wild));
% 
% % get parameters
% params = nan(length(fn), np);
% 
% 
% for ww = 1:length(fn)
% 
%     fit_t = readtable(fullfile(fn(ww).folder, fn(ww).name));
% 
% 
%     name = split(fn(ww).name, {'-','_'});
%     sample_num = str2double(name{end-1});
% 
%     for pp = 1:np
% 
%         params(ww,pp) = fit_t{ismember(fit_t{:,1}, p_name(pp)),2};
%     end
% 
% end
% 






%% recover key parameters (summary)

% 
% f=figure;
% set(f,'defaultTextInterpreter','none')
% 
% for nn = 1:n_params
% 
% 
%     nexttile; hold on;
% 
%     [ff,xf,u]=ksdensity(params(:,nn));
% 
%     orig_trace = orig.(p_name{nn});
%     [fo,xo]=ksdensity(orig_trace, 'Bandwidth', u);
% 
% 
%     % orig
%     plot(xo,fo, '-g', 'LineWidth', 3);
%     plot(xf,ff, '--k', 'LineWidth', 3);
% 
%     % fit
%     set(gca, 'TickDir', 'out', 'LineWidth', 1)
%     title(p_name(nn))
% 
% 
% 
% 
%     if nn==1
% 
%         legend('generating', 'fit')
% 
%     end
% 
% end
% 



%% recover all parameters (don't)



% 
% 
% 
% 
% wild = sprintf('gen-%s__fit-%s*', model_rec1, model_rec2)
% fn = dir(sprintf('../recovery-local/summary/%s/%s', wild, wild));
% 
% % get parameters
% all_params = [];
% 
% for ww = 1:length(fn)
% 
% 
%     fit_t = readtable(fullfile(fn(ww).folder, fn(ww).name));
% 
%     name = split(fn(ww).name, {'-','_'});
%     sample_num = str2double(name{end-1});
% 
%     name = split(fn(ww).name, {'-','_'});
%     all_params(ww,:) = fit_t{2:end-1,2};
% 
% end
% 
% all_params(all_params==0) = nan;
% 
% 
% 
% f=figure;
% set(f,'defaultTextInterpreter','none')
% 
% for nn = 1:size(all_params,2)
% 
% 
%     nexttile; hold on;
% 
%     nn_name = strrep(fit_t{nn+1,1}{:}, '.', '_');
% 
% 
% 
%     [ff,xf,u]=ksdensity(all_params(:,nn), 'Bandwidth', 'normal-approx');
% 
%     orig_trace = orig.(nn_name);
%     [fo,xo]=ksdensity(orig_trace, 'Bandwidth', u);
%     title(nn_name)
% 
% 
%     % orig
%     plot(xo,fo, '-g', 'LineWidth', 2);
%     plot(xf,ff, '--k', 'LineWidth', 2);
%     set(gca, 'TickDir', 'out', 'LineWidth', 1)
% 
% 
% 
% end
% 
% 
% 












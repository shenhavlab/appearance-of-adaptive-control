%% make figures



%% ========= figure 1 =====================================================

%% posterior predictive check (VD ev) ==================================

clear; clc;

%% load

t_ppc = readtable('PATH_TO_VDev_subj_ppc.csv');
t_param = readtable('PATH_TO_VDev_subj_traces.csv');

disp('loaded')

%% plot

figure('units','inch','position',[0,0,15,5]); hold on;
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([1,2]); hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

yline(0, '-k', 'LineWidth',2)
xlim([.200, .800])
xlabel('RT')
ylabel('PDF')
title('PPC (Original)')

% easy
% sel = ismember(t_ppc.hard, 'easy');
sel = t_ppc.hard == 0 & (t_ppc.signVD ~= 0);

acc = mean(t_ppc.response(sel) == sign(t_ppc.signVD(sel)));
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response == sign(t_ppc.signVD) & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
feasy=plot(x, f*acc, '-g', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response ~= sign(t_ppc.signVD) & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '-g', 'LineWidth', 2);


acc = mean(t_ppc.response_sampled(sel)== sign(t_ppc.signVD(sel)));
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled == sign(t_ppc.signVD) & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, f*acc, '--g', 'LineWidth', 1.5);
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled ~= sign(t_ppc.signVD) & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '--g', 'LineWidth', 1.5);


% hard
% sel = ismember(t_ppc.hard, 'hard');
sel = (t_ppc.hard == 1) & (t_ppc.signVD ~= 0);

acc = mean(t_ppc.response(sel) == sign(t_ppc.signVD(sel)));
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response == sign(t_ppc.signVD) & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
fhard=plot(x, f*acc, '-m', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response ~= sign(t_ppc.signVD) & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '-m', 'LineWidth', 2);

acc = mean(t_ppc.response_sampled(sel)== sign(t_ppc.signVD(sel)));
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled == sign(t_ppc.signVD) & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, f*acc, '--m', 'LineWidth', 1.5);
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled ~= sign(t_ppc.signVD) & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '--m', 'LineWidth', 1.5);



legend([feasy, fhard], {'easy', 'hard'})






% drift
% nexttile; hold on;
% set(gca, 'TickDir', 'out', 'LineWidth', 1)
% 
% xlabel('drift')
% ylabel('PDF')
% title('drift by difficulty (Original)')
% 
% 
% A = t_param.v_signVD;
% A_lo = prctile(A, 2.5);
% A_hi = prctile(A, 97.5);
% 
% 
% [f,x]=ksdensity(A, 'NumPoints', 1000, 'Bandwidth', .001);
% plot(x, f, '-k', 'LineWidth', 1);
% plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-k', 'LineWidth', 3);
% % xlim([0 inf])

% threshold
nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

xlabel('threshold')
ylabel('PDF')
title('threshold by difficulty (Original)')



A = t_param.a_absVD;
A_lo = prctile(A, 2.5);
A_hi = prctile(A, 97.5);
a_absVD = mean(A>0)


[f,x]=ksdensity(A, 'NumPoints', 1000, 'Bandwidth', .001);
plot(x, f, '-k', 'LineWidth', 1);
plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-k', 'LineWidth', 3);
xline(0,'-k')







%% ========= figure 2 =====================================================


resp_t = readtable('PATH_TO_vcr_respCoded.csv');
orig_t = readtable('PATH_TO_vcr_orig.csv');
disp('loaded')

dEV = -2675.975059;
pEV = 59.634829;

dOrig = -2117.882290;
pOrig = 70.719346;

dic_VD_orig = -2047.162944;
dic_VDev_subj =  -2616.340230;


dic_VDev_subj - (dEV + pEV)
dic_VD_orig - (dOrig + pOrig)

orig_VDev = dic_VDev_subj - dic_VD_orig
orig_VDev_likNorm = (dEV/height(resp_t) + pEV) - (dOrig/height(orig_t) + pOrig)



%% model comparision ==================================

clear; clc;

dic_VD_orig = -2047.162944;
dic_VD_static = -2741.329206;
dic_VDcomb_static = -2736.005021;
dic_VD_collapse = -3877.349246;
dic_VD_collapseA =  -3879.871243;
dic_VD_collapseT =  -3877.410861;

dic_OV_static = -2775.329145;
dic_OV_collapse = -3939.723838;
dic_OV_collapseA = -3932.680471;
dic_OV_collapseT = -3901.452555;

dic_VDOV_static = -2776.673388;
dic_VDOV_collapse = -3934.211954;
dic_VDOV_collapseA = -3927.125863;
dic_VDOV_collapseT = -3907.447633;




dic_VDev_subj =  -2616.340230;
dic_VD_static_subj = -2735.160811;
dic_OV_static_subj = -2770.491619;
dic_VDOV_static_subj = -2842.494853;



dic_VD_collapse_subj = -3941.670806;
dic_VD_collapseA_subj =  -3946.283452;
dic_VD_collapseT_subj =  -3940.183433;
dic_OV_collapse_subj = -4013.616978;
dic_OV_collapseA_subj = -4017.501340;
dic_OV_collapseT_subj = -3986.664735;
dic_VDOV_collapse_subj = -4026.686943;
dic_VDOV_collapseA_subj = -4025.264462;
dic_VDOV_collapseT_subj = -4001.994055;








% static_dics = [...
%     dic_VD_static,...
%     dic_OV_static,...
%     dic_VDOV_static,...
%     ];


static_dics = [...
    dic_VDev_subj,...
    dic_VD_static_subj,...
    dic_OV_static_subj,...
    dic_VDOV_static_subj,...
    ];


dif_within_static = static_dics - static_dics(1)



% collapse_dics = [...
%     dic_VD_collapse,...
%     dic_VD_collapseA,...
%     dic_VD_collapseT,...
%     dic_OV_collapse,...
%     dic_OV_collapseA,...
%     dic_OV_collapseT,...
%     dic_VDOV_collapse,...
%     dic_VDOV_collapseA,...
%     dic_VDOV_collapseT,...
%     ];


collapse_dics = [...
    dic_VD_collapse_subj,...
    dic_VD_collapseA_subj,...
    dic_VD_collapseT_subj,...
    dic_OV_collapse_subj,...
    dic_OV_collapseA_subj,...
    dic_OV_collapseT_subj,...
    dic_VDOV_collapse_subj,...
    dic_VDOV_collapseA_subj,...
    dic_VDOV_collapseT_subj,...
    ];


collase_effect = fix([dic_VD_static_subj - dic_VD_collapse_subj, dic_OV_static_subj - dic_OV_collapse_subj, dic_VDOV_static_subj - dic_VDOV_collapse_subj])

vdov_effect = [dic_VD_static_subj - dic_VDOV_static_subj, dic_VD_collapse_subj - dic_VDOV_collapse_subj]



min_dic = min([static_dics, collapse_dics]);




figure('units','inch','position',[0,0,7,7]); hold on;
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% plot static
nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

ylim([0, length(static_dics)+1])

plot(static_dics-min_dic, length(static_dics):-1:1,  '+k', 'MarkerSize',12, 'LineWidth',1)

yticks(1:length(static_dics));
yticklabels(flip({'VD_{sum}', 'VD_{bias}', 'OV_{bias}', 'VDOV_{bias}'}))
xlabel('DIC - best')
title('fixed bound models')
xlim([min(static_dics-min_dic)-50, inf])

% plot collapse
nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

ylim([0, length(collapse_dics)+1])

plot(collapse_dics-min_dic, length(collapse_dics):-1:1,  '+k', 'MarkerSize',12, 'LineWidth',1)
plot(0, length(collapse_dics)-(find(collapse_dics==min_dic)-1),  '+g', 'MarkerSize',15, 'LineWidth',2)

xline(0)

yticks(1:length(collapse_dics));
yticklabels(flip({...
    'VD_{both}', 'VD_{init}', 'VD_{rate}',...
    'OV_{both}', 'OV_{init}', 'OV_{rate}',...
    'VDOV_{both}', 'VDOV_{init}', 'VDOV_{rate}',...
    }))
xlabel('DIC - best')
xlim([-5, inf])
title('collapsing bound models')



%% rhat (Collapse) ==================================
t = readtable('PATH_TO_VDOV_collapse_subj_rhat.csv');
t = t(2:end,:);
plot(sort(t{:,2}), 1:height(t), 'ok')
xline(1.1, '-k')
xline(max(t{:,2}), '-r')


%% posterior predictive check (Collapse) ==================================

clear; clc;

%% load

t_ppc_orig = readtable('PATH_TO_VDev_sAll_subj_ppc.csv');
t_ppc = readtable('PATH_TO_VDev_collapse_subj_ppc.csv');




disp('loaded')

%% plot



figure('units','inch','position',[0,0,15,5]); hold on;
% tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');
% nexttile([1,2]); hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

yline(0, '-k', 'LineWidth',2)
xlim([.200, .800])
xlabel('RT')
ylabel('PDF')
title('PPC (Best)')





% ====================== ORIG MODEL

corResp = sign(t_ppc_orig.signVD);
corResp(corResp==-1) = 0;

sel = ismember(t_ppc_orig.hard, 0) & (t_ppc_orig.signVD ~= 0);

% easy


acc = mean(t_ppc_orig.response_sampled(sel) == corResp(sel));
[f,x]=ksdensity(abs(t_ppc_orig.rt_sampled(t_ppc_orig.response_sampled == corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, f*acc, ':g', 'LineWidth', 1.5);
[f,x]=ksdensity(abs(t_ppc_orig.rt_sampled(t_ppc_orig.response_sampled ~= corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), ':g', 'LineWidth', 1.5);


% hard
sel = ismember(t_ppc_orig.hard, 1) & (t_ppc_orig.signVD ~= 0);

acc = mean(t_ppc_orig.response_sampled(sel) == corResp(sel));
[f,x]=ksdensity(abs(t_ppc_orig.rt_sampled(t_ppc_orig.response_sampled == corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, f*acc, ':m', 'LineWidth', 1.5);
[f,x]=ksdensity(abs(t_ppc_orig.rt_sampled(t_ppc_orig.response_sampled ~= corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), ':m', 'LineWidth', 1.5);












% ====================== BEST MODEL

corResp = sign(t_ppc.signVD);
% corResp(corResp==-1) = 0;


% easy
sel = ismember(t_ppc.hard, 0) & (t_ppc.signVD ~= 0);

acc = mean(t_ppc.response(sel) == corResp(sel));
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response == corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
feasy=plot(x, f*acc, '-g', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response ~= corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '-g', 'LineWidth', 2);


acc = mean(t_ppc.response_sampled(sel) == corResp(sel));
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled == corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, f*acc, '--g', 'LineWidth', 1.5);
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled ~= corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '--g', 'LineWidth', 1.5);


% hard
sel = ismember(t_ppc.hard, 1) & (t_ppc.signVD ~= 0);

acc = mean(t_ppc.response(sel) == corResp(sel));
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response == corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
fhard=plot(x, f*acc, '-m', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response ~= corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '-m', 'LineWidth', 2);


acc = mean(t_ppc.response_sampled(sel) == corResp(sel));
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled == corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, f*acc, '--m', 'LineWidth', 1.5);
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled ~= corResp & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '--m', 'LineWidth', 1.5);



legend([feasy, fhard], {'easy', 'hard'})




%% plot parameters


t_param = readtable('PATH_TO_VDOV_collapse_subj_traces.csv');
t_paramOrig = readtable('PATH_TO_VDev_subj_traces.csv');
disp('loaded')

%% plot key parameters


figure('units','inch','position',[0,0,10,5]); hold on;
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');
set(gca, 'TickDir', 'out', 'LineWidth', 1)
 


% initial bound
nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

xlabel('threshold')
ylabel('PDF')
title('inital bound')


[f,x]=ksdensity(A, 'NumPoints', 1000, 'Bandwidth', .001);
plot(x, f, '-b', 'LineWidth', 1);
plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-b', 'LineWidth', 3);


A = t_param.a_absMinVD;
A_lo = prctile(A, 2.5)
A_hi = prctile(A, 97.5)

a_absMinVD = mean(A>0)


[f,x]=ksdensity(A, 'NumPoints', 1000, 'Bandwidth', .001);
plot(x, f, '-k', 'LineWidth', 1);
plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-k', 'LineWidth', 3);
xline(0, '-k');



A = t_paramOrig.a_absVD;
A_lo = prctile(A, 2.5)
A_hi = prctile(A, 97.5)

a_absVD = mean(A>0)


[f,x]=ksdensity(A, 'NumPoints', 1000, 'Bandwidth', .001);
plot(x, f, ':k', 'LineWidth', 1);
plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), ':k', 'LineWidth', 2);
xline(0, '-k');




% initial bound
nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

xlabel('threshold')
ylabel('PDF')
title('overall value')


A = t_param.a_maxOV;
A_lo = prctile(A, 2.5);
A_hi = prctile(A, 97.5);

a_maxOV = mean(A>0)


[f,x]=ksdensity(A, 'NumPoints', 1000, 'Bandwidth', .001);
plot(x, f, '-b', 'LineWidth', 1);
plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-b', 'LineWidth', 3);


A = t_param.a_minOV;
A_lo = prctile(A, 2.5)
A_hi = prctile(A, 97.5)

a_minOV = mean(A>0)


[f,x]=ksdensity(A, 'NumPoints', 1000, 'Bandwidth', .001);
plot(x, f, '-k', 'LineWidth', 1);
plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-k', 'LineWidth', 3);
xline(0, '-k');









% rate
% nexttile; hold on;
% set(gca, 'TickDir', 'out', 'LineWidth', 1)
% 
% xlabel('threshold')
% ylabel('PDF')
% title('collapse rate')
% 
% 
% A = t_param.theta_absMaxVD
% A_lo = prctile(A, 2.5);
% A_hi = prctile(A, 97.5);
% 
% 
% [f,x]=ksdensity(A, 'NumPoints', 1000);
% plot(x, f, '-b', 'LineWidth', 1);
% plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-b', 'LineWidth', 3);
% 
% 
% A = t_param.theta_absMinVD
% A_lo = prctile(A, 2.5);
% A_hi = prctile(A, 97.5);
% 
% 
% [f,x]=ksdensity(A, 'NumPoints', 1000);
% plot(x, f, '-k', 'LineWidth', 1);
% plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-k', 'LineWidth', 3);





%% plot all parameters



tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
set(gca, 'TickDir', 'out', 'LineWidth', 1)



param_idx = find(~cellfun(@(x) contains(x, {'subj', 'Var', 'std'}, 'IgnoreCase', false), t_param.Properties.VariableNames));
for ii = 1:length(param_idx)

    param_name = t_param.Properties.VariableNames{param_idx(ii)};
  
    param_val = t_param{:, param_idx(ii)};
  param_lo = prctile(param_val, 2.5);
    param_hi = prctile(param_val, 97.5);

    % initial bound
    nexttile; hold on;
    set(gca, 'TickDir', 'out', 'LineWidth', 1)

    xlabel('')
    ylabel('PDF')
    title(param_name)



    [f,x]=ksdensity(param_val, 'NumPoints', 1000);
    plot(x, f, '-b', 'LineWidth', 1);
    plot(x(x>param_lo & x<param_hi), f(x>param_lo & x<param_hi), '-b', 'LineWidth', 3);
    xline(0, '-k');

end











%% ========= figure 3 =====================================================


%% threshold effect (LCA) ==================================

clear; clc;

%% load


% ORIG* MODEL (VDev) ----------------------


% NEW FIT LCA  *********
t_beh_orig = readtable('PATH_TO_VDev_subj_traces.csv');
disp('fitting LCA-fit orig')
t_lca_orig = readtable('PATH_TO_VDev_traces.csv');




% BEST MODEL (VDOV) ----------------------

t_beh_best = readtable('PATH_TO_VDOV_collapse_subj_traces.csv');

% NEW FIT LCA *********
disp('fitting LCA-fit best')
t_lca_best = readtable('PATH_TO_VDOV_collapse_traces.csv');


disp('loaded')

%% plot






figure('units','inch','position',[0,0,15,5]); hold on;
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');




nexttile; hold on;
xline(0, '-k', 'LineWidth',1)
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlabel('parameter')
title('parameter mimicry (VDOV collapse)')



% a ~ max OV
errorbar(mean(t_beh_best.a_maxOV), 1.55, ...
    abs(prctile(t_beh_best.a_maxOV, 2.5)-mean(t_beh_best.a_maxOV)), abs(prctile(t_beh_best.a_maxOV, 97.5)-mean(t_beh_best.a_maxOV)),...
    'horizontal','ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)

errorbar(mean(t_lca_best.a_maxOV),1.45,  ...
    abs(prctile(t_lca_best.a_maxOV, 2.5)-mean(t_lca_best.a_maxOV)), abs(prctile(t_lca_best.a_maxOV, 97.5)-mean(t_lca_best.a_maxOV)),...
   'horizontal', 'og', 'MarkerFaceColor', 'g', 'LineWidth', 1, 'MarkerSize', 10)



ylim([1.25,1.75]);
yticks([1.5]);
yticklabels({'a ~ OV'})
xlim([-.1, .1])


nexttile; hold on;
xline(0, '-k', 'LineWidth',1)
set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlabel('parameter')
title('parameter mimicry (VD fixed)')




errorbar(-mean(t_beh_orig.a_absVD),1.55,  ...
    abs(prctile(t_beh_orig.a_absVD, 2.5)-mean(t_beh_orig.a_absVD)), abs(prctile(t_beh_orig.a_absVD, 97.5)-mean(t_beh_orig.a_absVD)),...
     'horizontal', 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)


errorbar(-mean(t_lca_orig.a_absVD), 1.45, ...
    abs(prctile(t_lca_orig.a_absVD, 2.5)-mean(t_lca_orig.a_absVD)), abs(prctile(t_lca_orig.a_absVD, 97.5)-mean(t_lca_orig.a_absVD)),...
     'horizontal', 'og', 'MarkerFaceColor', 'g', 'LineWidth', 1, 'MarkerSize', 10)


ylim([1.25,1.75]);
yticks([1.5]);
yticklabels({'a ~ VD'})
xlim([-.015, .015])
















% --- ALL PARAMETERS


figure;



nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

xline(0, '-k', 'LineWidth',1)
xlabel('parameter')
title('parameter mimicry (VDOV collapse)')



% a ~ min OV
errorbar(mean(t_beh_best.a_minOV), 1.55,...
    abs(prctile(t_beh_best.a_minOV, 2.5)-mean(t_beh_best.a_minOV)), abs(prctile(t_beh_best.a_minOV, 97.5)-mean(t_beh_best.a_minOV)),...
    'horizontal', 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)

errorbar(mean(t_lca_best.a_minOV), 1.45,...
    abs(prctile(t_lca_best.a_minOV, 2.5)-mean(t_lca_best.a_minOV)), abs(prctile(t_lca_best.a_minOV, 97.5)-mean(t_lca_best.a_minOV)),...
    'horizontal', 'og', 'MarkerFaceColor', 'g', 'LineWidth', 1, 'MarkerSize', 10)

% a ~ max OV
errorbar(mean(t_beh_best.a_maxOV), 2.05,...
    abs(prctile(t_beh_best.a_maxOV, 2.5)-mean(t_beh_best.a_maxOV)), abs(prctile(t_beh_best.a_maxOV, 97.5)-mean(t_beh_best.a_maxOV)),...
    'horizontal', 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)

errorbar(mean(t_lca_best.a_maxOV), 1.95,...
    abs(prctile(t_lca_best.a_maxOV, 2.5)-mean(t_lca_best.a_maxOV)), abs(prctile(t_lca_best.a_maxOV, 97.5)-mean(t_lca_best.a_maxOV)),...
    'horizontal', 'og', 'MarkerFaceColor', 'g', 'LineWidth', 1, 'MarkerSize', 10)




% rate ~ min OV
errorbar(mean(t_beh_best.theta_minOV), 2.55,...
    abs(prctile(t_beh_best.theta_minOV, 2.5)-mean(t_beh_best.theta_minOV)), abs(prctile(t_beh_best.theta_minOV, 97.5)-mean(t_beh_best.theta_minOV)),...
    'horizontal', 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)

errorbar(mean(t_lca_best.theta_minOV), 2.45,...
    abs(prctile(t_lca_best.theta_minOV, 2.5)-mean(t_lca_best.theta_minOV)), abs(prctile(t_lca_best.theta_minOV, 97.5)-mean(t_lca_best.theta_minOV)),...
    'horizontal', 'og', 'MarkerFaceColor', 'g', 'LineWidth', 1, 'MarkerSize', 10)


% rate ~ max OV
errorbar(mean(t_beh_best.theta_maxOV), 3.05,...
    abs(prctile(t_beh_best.theta_maxOV, 2.5)-mean(t_beh_best.theta_maxOV)), abs(prctile(t_beh_best.theta_maxOV, 97.5)-mean(t_beh_best.theta_maxOV)),...
    'horizontal', 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)

errorbar(mean(t_lca_best.theta_maxOV), 2.95,...
    abs(prctile(t_lca_best.theta_maxOV, 2.5)-mean(t_lca_best.theta_maxOV)), abs(prctile(t_lca_best.theta_maxOV, 97.5)-mean(t_lca_best.theta_maxOV)),...
    'horizontal', 'og', 'MarkerFaceColor', 'g', 'LineWidth', 1, 'MarkerSize', 10)








% a ~ min VD
errorbar(mean(t_beh_best.a_absMinVD), 3.55,...
    abs(prctile(t_beh_best.a_absMinVD, 2.5)-mean(t_beh_best.a_absMinVD)), abs(prctile(t_beh_best.a_absMinVD, 97.5)-mean(t_beh_best.a_absMinVD)),...
    'horizontal', 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)

errorbar(mean(t_lca_best.a_absMinVD), 3.45,...
    abs(prctile(t_lca_best.a_absMinVD, 2.5)-mean(t_lca_best.a_absMinVD)), abs(prctile(t_lca_best.a_absMinVD, 97.5)-mean(t_lca_best.a_absMinVD)),...
    'horizontal', 'og', 'MarkerFaceColor', 'g', 'LineWidth', 1, 'MarkerSize', 10)

% a ~ max VD
errorbar(mean(t_beh_best.a_absMaxVD), 4.05,...
    abs(prctile(t_beh_best.a_absMaxVD, 2.5)-mean(t_beh_best.a_absMaxVD)), abs(prctile(t_beh_best.a_absMaxVD, 97.5)-mean(t_beh_best.a_absMaxVD)),...
    'horizontal', 'ok', 'MarkerFaceColor', 'k', 'LineWidth', 1, 'MarkerSize', 10)

errorbar(mean(t_lca_best.a_absMaxVD), 3.95,...
    abs(prctile(t_lca_best.a_absMaxVD, 2.5)-mean(t_lca_best.a_absMaxVD)), abs(prctile(t_lca_best.a_absMaxVD, 97.5)-mean(t_lca_best.a_absMaxVD)),...
    'horizontal', 'og', 'MarkerFaceColor', 'g', 'LineWidth', 1, 'MarkerSize', 10)









figure;
nexttile; hold on;
plot(t_beh_orig.v_signVD, t_beh_orig.a_absVD, 'ok')
corr(t_beh_orig.v_signVD, t_beh_orig.a_absVD)
title('behav')

nexttile; hold on;
plot(t_lca_orig.v_signVD, t_lca_orig.a_absVD, 'ok')
corr(t_lca_orig.v_signVD, t_lca_orig.a_absVD)
title('LCA')



%% sim LCA behavior
clear t_ppc t_ppc_beh
t_lca = readtable('PATH_TO_2025-08-22_14-38-02_newLCA-compare.csv')

t_lca=t_lca(t_lca.signVD~=0,:);

disp('loaded')


%%



% -- USING JUST LCA PPC DATA

figure; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

yline(0, '-k', 'LineWidth',2)
xlim([.200, .800])
xlabel('RT')
ylabel('PDF')
title('PPC (LCA) 2')




% LCA MODEL -----------
% t_lca.response(t_lca.response == 1) = 1;
% t_lca.response(t_lca.response == 2) = -1;

sel = t_lca.absVD > median(t_lca.absVD);
[f,x]=ksdensity(abs(t_lca.rt(t_lca.sim_acc)), 'Bandwidth', 1/60, 'NumPoints', 1000);
feasy=plot(x, f*mean(t_lca.sim_acc), '-g', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_lca.rt(~t_lca.sim_acc)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*mean(~t_lca.sim_acc), '-g', 'LineWidth', 2);


sel = t_lca.absVD < median(t_lca.absVD);

[f,x]=ksdensity(abs(t_lca.rt(t_lca.sim_acc)), 'Bandwidth', 1/60, 'NumPoints', 1000);
feasy=plot(x, f*mean(t_lca.sim_acc), '-m', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_lca.rt(~t_lca.sim_acc)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*mean(~t_lca.sim_acc), '-m', 'LineWidth', 2);




% Participant Behaviour -----------



mean(t_lca.orig_acc)
mean(t_lca.orig_resp == sign(t_lca.signVD))

% easy
sel = t_lca.absVD > median(t_lca.absVD);
[f,x]=ksdensity(abs(t_lca.orig_rt(t_lca.orig_acc)), 'Bandwidth', 1/60, 'NumPoints', 1000);
feasy=plot(x, f*mean(t_lca.orig_acc), '--g', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_lca.orig_rt(~t_lca.orig_acc)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*mean(~t_lca.orig_acc), '--g', 'LineWidth', 2);

% hard
sel = t_lca.absVD < median(t_lca.absVD);
[f,x]=ksdensity(abs(t_lca.orig_rt(t_lca.orig_acc)), 'Bandwidth', 1/60, 'NumPoints', 1000);
feasy=plot(x, f*mean(t_lca.orig_acc), '--m', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_lca.orig_rt(~t_lca.orig_acc)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*mean(~t_lca.orig_acc), '--m', 'LineWidth', 2);


















%% ========= SUPP figure 1 =====================================================

%% posterior predictive check (ORIG) ==================================

clear; clc;

%% load

t_ppc = readtable('PATH_TO_VD_orig_ppc.csv');
t_param = readtable('PATH_TO_VD_orig_traces.csv');


disp('loaded')

%% plot

figure('units','inch','position',[0,0,15,5]); hold on;
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile([1,2]); hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

yline(0, '-k', 'LineWidth',2)
xlim([.200, .800])
xlabel('RT')
ylabel('PDF')
title('PPC (Original)')

% easy
sel = ismember(t_ppc.hard, 'easy');
% sel = t_ppc.hard == 0;

acc = mean(t_ppc.response(sel));
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response == 1 & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
feasy=plot(x, f*acc, '-g', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response ~= 1 & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '-g', 'LineWidth', 2);


acc = mean(t_ppc.response_sampled(sel));
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled == 1 & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, f*acc, '--g', 'LineWidth', 1.5);
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled ~= 1 & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '--g', 'LineWidth', 1.5);


% hard
sel = ismember(t_ppc.hard, 'hard');
% sel = t_ppc.hard == 1;

acc = mean(t_ppc.response(sel));
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response == 1 & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
fhard=plot(x, f*acc, '-m', 'LineWidth', 2);
[f,x]=ksdensity(abs(t_ppc.rt(t_ppc.response ~= 1 & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '-m', 'LineWidth', 2);

acc = mean(t_ppc.response_sampled(sel));
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled == 1 & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, f*acc, '--m', 'LineWidth', 1.5);
[f,x]=ksdensity(abs(t_ppc.rt_sampled(t_ppc.response_sampled ~= 1 & sel)), 'Bandwidth', 1/60, 'NumPoints', 1000);
plot(x, -f*(1-acc), '--m', 'LineWidth', 1.5);



legend([feasy, fhard], {'easy', 'hard'})






% drift
% nexttile; hold on;
% set(gca, 'TickDir', 'out', 'LineWidth', 1)
% 
% xlabel('drift')
% ylabel('PDF')
% title('drift by difficulty (Original)')
% 
% 
% A = t_param.v_signVD;
% A_lo = prctile(A, 2.5);
% A_hi = prctile(A, 97.5);
% 
% 
% [f,x]=ksdensity(A, 'NumPoints', 1000, 'Bandwidth', .001);
% plot(x, f, '-k', 'LineWidth', 1);
% plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-k', 'LineWidth', 3);
% % xlim([0 inf])

% threshold
nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)

xlabel('threshold')
ylabel('PDF')
title('threshold by difficulty (Original)')



A = t_param.a_easy_ - t_param.a_hard_;
A_lo = prctile(A, 2.5);
A_hi = prctile(A, 97.5);

a_hard_ = mean(A>0)


[f,x]=ksdensity(A, 'NumPoints', 1000, 'Bandwidth', .001);
plot(x, f, '-k', 'LineWidth', 1);
plot(x(x>A_lo & x<A_hi), f(x>A_lo & x<A_hi), '-k', 'LineWidth', 3);
xline(0,'-k')





%% ========= SUPP figure 2 =====================================================

%% OV confound  ==================================


t_vd = readtable('PATH_TO_VDOV_collapse_subj_traces.csv');
t_dat = readtable('PATH_TO_speeded_behavioral_data.csv');

t_dat2 = readtable('PATH_TO_vcr_respCoded.csv');


disp('loaded')

%% plot

center = @(x) x-nanmean(x);

% t_dat = t_dat(t_dat.subjectID == 2, :);

maxVD =   abs(t_dat.MaxValue_left - t_dat.MaxValue_right);
minVD =   abs(t_dat.MinValue_left - t_dat.MinValue_right);
maxOV =   t_dat.MaxValue_left + t_dat.MaxValue_right;
minOV =   t_dat.MinValue_left + t_dat.MinValue_right;
OV    = (maxOV + minOV)/2;
VD = (maxVD + minVD)/2;
subj = t_dat.subjectID;







% vd_max = mean(t_vd.a_absMaxVD);
% vd_min = mean(t_vd.a_absMinVD);
vd_max = mean(t_vd.v_maxVD);
vd_min = mean(t_vd.v_minVD);
vd_sum = vd_max + vd_min;
vd_max = vd_max/vd_sum;
vd_min = vd_min/vd_sum;


% ov_max = mean(t_vd.a_maxOV);
% ov_min = mean(t_vd.a_minOV);
ov_max = mean(t_vd.v_maxVD);
ov_min = mean(t_vd.v_minVD);
ov_sum = ov_max + ov_min;
ov_max = ov_max/ov_sum;
ov_min = ov_min/ov_sum;



figure('units','inch','position',[0,0,9,3]); hold on;
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');


% orig
nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('objective - objective')
xlabel('objective VD')
ylabel('objective OV')

plot(VD, OV, 'ok'); lsline;
[r,p]=corr(VD, OV)
xlim([0, 60])
ylim([40, 150])


% confound
nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('subjective - subjective')
xlabel('subjective VD')
ylabel('subjective OV')

sub_VD = maxVD*vd_max + minVD*vd_min;
sub_OV = maxOV*ov_max + minOV*ov_min;  

plot(sub_VD, sub_OV, 'ok'); lsline;
[r,p]=corr(sub_VD, sub_OV)
xlim([0, 60])
ylim([40, 150])



% per-subj confound




subs = unique(t_dat.subjectID);
clear subj_r
for ii = 1:length(subs)

    vd_max = mean(t_vd.(sprintf('v_maxVD_subj_%d', subs(ii))));
    vd_min = mean(t_vd.(sprintf('v_minVD_subj_%d', subs(ii))));
    vd_sum = vd_max + vd_min;
    vd_max = vd_max/vd_sum;
    vd_min = vd_min/vd_sum;


    ov_max = mean(t_vd.(sprintf('v_maxVD_subj_%d', subs(ii))));
    ov_min = mean(t_vd.(sprintf('v_minVD_subj_%d', subs(ii))));
    ov_sum = ov_max + ov_min;
    ov_max = ov_max/ov_sum;
    ov_min = ov_min/ov_sum;



    sub_VD = maxVD(subj == subs(ii))*vd_max + minVD(subj == subs(ii))*vd_min;
    sub_OV = maxOV(subj == subs(ii))*ov_max + minOV(subj == subs(ii))*ov_min;

    subj_r(ii) = corr(sub_VD, sub_OV);


end

nexttile; hold on;
set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('value correlations')
xlabel('value correlation')
ylabel('density')
[f,x]=ksdensity(subj_r);
plot(x,f, '-k', 'LineWidth', 2)
% histogram(subj_r);
xline(0, '-k')








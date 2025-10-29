%% Fit LCA


clear;



addpath(genpath('./bads'))
addpath(genpath('./ibs'))

n_inits = 10;

max_bads_iter = 100
MaxFunEvals = 1000
use_targetSD = true

ibs_maxtime = 1.0

global_disc = false
n_discs = 10
dt = .001

center_rt = true


savename = 'ibs-20s'

%% prep data
PATH = 'PATH_TO_DATA';
D=readtable(sprintf('%sspeeded_behavioral_data.csv', PATH));
D = D(~D.responseID==0, :); % removing 3 trials where there was no response
D.responseID(D.responseID==-1) = 0;
tablesize = height(D);

pts = unique(D.subjectID);
n_pts = length(pts)

%% convert to accuracy code

% signVD = sign(D.ExpectedValue_left - D.ExpectedValue_right);

Vs = [ D.ExpectedValue_right, D.ExpectedValue_left];
[~,max_idx] = max(Vs,[],2);
[~,min_idx] = min(Vs,[],2);

D.acc = double((D.responseID+1) == max_idx);
D.rt = D.ResponseTime_seconds;
D.acc(isnan(D.rt))=0;% treat misses as errors

D.maxRT = zeros(size(D.rt)) + .750;
D.resp = D.responseID;



if center_rt
    % center = @(x) x-nanmean(x)
    grand_mean = nanmean(D.rt);

    for pp = 1:n_pts

        mRT = nanmean(D.rt(D.subjectID==pts(pp)));

        D.rt(D.subjectID==pts(pp)) = D.rt(D.subjectID==pts(pp)) - mRT;
        D.maxRT(D.subjectID==pts(pp)) = D.maxRT(D.subjectID==pts(pp)) - mRT;

    end

    D.rt = D.rt + grand_mean;
    D.maxRT = D.maxRT + grand_mean;

end




assert(mean(D.acc) > .5)

[D.CorrMaxVal, D.CorrMinVal, D.ErrMaxVal, D.ErrMinVal,D.CorrVal, D.ErrVal] = deal(D.ExpectedValue_right);
for tt = 1:tablesize
    if max_idx(tt) == 1

        D.CorrVal(tt) = D.ExpectedValue_right(tt);
        D.ErrVal(tt) = D.ExpectedValue_left(tt);

        D.CorrMaxVal(tt) = D.MaxValue_right(tt);
        D.CorrMinVal(tt) = D.MinValue_right(tt);
        D.ErrMaxVal(tt) = D.MaxValue_left(tt);
        D.ErrMinVal(tt) = D.MinValue_left(tt);

    else

        D.CorrVal(tt) = D.ExpectedValue_left(tt);
        D.ErrVal(tt) = D.ExpectedValue_right(tt);

        D.CorrMaxVal(tt) = D.MaxValue_left(tt);
        D.CorrMinVal(tt) = D.MinValue_left(tt);
        D.ErrMaxVal(tt) = D.MaxValue_right(tt);
        D.ErrMinVal(tt) = D.MinValue_right(tt);

    end

end

assert(all(D.CorrVal >= D.ErrVal))


%% Get VD and OV

D.OV = D.ExpectedValue_right+D.ExpectedValue_left;
D.discOV = D.OV > median(D.OV);

D.VD = abs(D.ExpectedValue_right-D.ExpectedValue_left);
D.discVD = D.VD > median(D.VD);


%% compute quantiles within each OV condition
cc =1;
quants = nan(4,n_discs+1);
mean_acc = nan(4,1);

[grp_subj, grp_OV, grp_VD, grp_maxRT, grp_RT, grp_rawRT, grp_acc,grp_resp,...
    grp_CorrMaxVal, grp_CorrMinVal, grp_ErrMaxVal,grp_ErrMinVal,...
    grp_LeftMaxVal, grp_LeftMinVal, grp_RightMaxVal, grp_RightMinVal] = deal([]);

for ov = 0:1
    for ac = 0:1

        Dg = D(D.discOV==ov & D.acc==ac,:);

        % create dmat
        grp_OV          = [grp_OV; Dg.discOV];
        grp_VD          = [grp_VD; Dg.discVD];
        grp_acc         = [grp_acc; Dg.acc];
        grp_resp        = [grp_resp; Dg.resp];

        grp_subj        = [grp_subj; Dg.subjectID];

        mean_acc(cc)    = nanmean(Dg.acc);

        quants(cc,:)    = quantile(Dg.rt, linspace(0,1,n_discs+1));
        % quants(cc,:)    = linspace(nanmin(D.rt), nanmax(D.rt),n_discs+1);
        grp_RT          = [grp_RT; discretize(Dg.rt, quants(cc,:))];
       
        grp_CorrMaxVal  = [grp_CorrMaxVal; Dg.CorrMaxVal];
        grp_CorrMinVal  = [grp_CorrMinVal; Dg.CorrMinVal];
        grp_ErrMaxVal   = [grp_ErrMaxVal; Dg.ErrMaxVal];
        grp_ErrMinVal   = [grp_ErrMinVal; Dg.ErrMinVal];


        grp_LeftMaxVal  = [grp_LeftMaxVal; Dg.MaxValue_left];
        grp_LeftMinVal  = [grp_LeftMinVal; Dg.MinValue_left];
        grp_RightMaxVal   = [grp_RightMaxVal; Dg.MaxValue_right];
        grp_RightMinVal   = [grp_RightMinVal; Dg.MinValue_right];


        grp_rawRT          = [grp_rawRT; Dg.rt];
        grp_maxRT       = [grp_maxRT; Dg.maxRT];


        cc=cc+1;
    end
end

% set missed responses to their own bin
grp_RT(isnan(grp_RT)) = 0;


if global_disc

    disp('RUNNING GLOBAL')
    quants = quantile(D.rt, linspace(0,1,n_discs+1));
    grp_RT = discretize(D.rt, quants);
    grp_RT(isnan(grp_RT)) = 0;

end



%% convert to matrix

designMat = [...
    grp_OV, ...1
    grp_acc,...2
    grp_CorrMaxVal/80, grp_CorrMinVal/80, grp_ErrMaxVal/80, grp_ErrMinVal/80, ... 3-6
    grp_VD,... 7
    grp_maxRT, ... 8
    grp_resp,... 9
    grp_LeftMaxVal, grp_LeftMinVal, grp_RightMaxVal, grp_RightMinVal,... 10-13
    grp_subj,... 14
    grp_RT,... 15
    grp_rawRT,... 16
    ];

% respMat = [...
%     grp_RT, grp_acc
%     ]; 

respMat = grp_RT; 


%% LCA parameters

badsopt = bads('defaults');
badsopt.UncertaintyHandling = true;
badsopt.SpecifyTargetNoise = use_targetSD;
badsopt.NoiseFinalSamples   = 30;
badsopt.MaxIter = max_bads_iter;
badsopt.Plot = 'off';
badsopt.MaxFunEvals=MaxFunEvals;
badsopt


param_names = {'t0', 'vin', 'vratio', 'leak', 'inhib', 'bound', 'collapse', 'sigma'};
x0 = [          0,     6,    .75,      1.75,   .03,     3,       1,          .8];
LB = [          0,     1,    .01,      0,      0,       2,       0,          .1];
UB = [          .300,  30,   .99,      5,      5,       20,      2,          10];
PLB = [         .200,  1,    .6,       .5,     .5,      3,       .5,         0.5];
PUB = [         .255,  10,   .9,       2,      2.0,     10,      1.5,        1.5];

n_param = length(param_names)

%% IBS options

ibsopt = ibslike('defaults');
ibsopt.Nreps = 10;
ibsopt.ReturnStd = true;
ibsopt.MaxTime = ibs_maxtime;
% ibsopt.MaxIter=5;
ibsopt

%% loglik options

llopt = struct;
llopt.quants = quants;
llopt.lates = mean(grp_RT==0);
llopt.accs = mean_acc;
llopt.global_disc = global_disc;
llopt.maxRT = nanmax(D.rt);
llopt.dt = dt;
llopt

%% create functions

inner_fun = @(x,dmat) lca_sim_ibs(x, dmat, llopt)
outer_fun = @(x) ibslike(inner_fun,x,respMat, designMat, ibsopt)


%% test inner


[resp,rt] = lca_sim_ibs(x0, designMat, llopt);

n_tests = 200
tic
for ii = 1:n_tests
    ll = inner_fun(x0, designMat);
end
inner_dur = toc;
fprintf('\ninner dur = %.4g s/iter\n', inner_dur/n_tests)


%% test outer


[NLOGL,NLOGLVAR,EXITFLAG,OUTPUT] = outer_fun(x0)



%% run fit

% [x,fval,exitflag,output,optimState,gpstruct] = bads(outer_fun,x0,LB,UB,PLB,PUB,[],badsopt)

rep_fval = nan(n_inits,1);
rep_x = nan(n_inits, n_param);
[rep_exitflag, rep_output, rep_optimState, rep_gpstruct] = deal(cell(n_inits,1));


time_on = string(datetime('now','Format','yyyy-MM-dd_HH-mm-ss'))
for rr = 1:n_inits
    tic
    
    if rr ==1
        [rep_x(rr,:),rep_fval(rr),rep_exitflag{rr},rep_output{rr},rep_optimState{rr},rep_gpstruct{rr}] = bads(outer_fun,x0,LB,UB,PLB,PUB,[],badsopt);
    else
        rep_x0 = unifrnd(PLB, PUB);
        [rep_x(rr,:),rep_fval(rr),rep_exitflag{rr},rep_output{rr},rep_optimState{rr},rep_gpstruct{rr}] = bads(outer_fun,rep_x0,LB,UB,PLB,PUB,[],badsopt);
    end

    % print
    tfit=toc;
    fprintf("\niter=%d // fval=%.2g // dur=%.2g min", rr, rep_fval(rr), tfit/60)

    % best
    [best_fval,best_fit]    = min(rep_fval)
    best_x          = rep_x(best_fit,:)
    best_exitflag   = rep_exitflag{best_fit};
    best_output     = rep_output{best_fit};
    best_optimState = rep_optimState{best_fit};
    best_gpstruct   = rep_gpstruct{best_fit};

    % save    
    save(sprintf('fit-results/%s_%s', time_on, savename))
   
end

best_x          
best_exitflag   
best_output     
best_optimState 
best_gpstruct   

%% slow eval at optima
center = @(x) x-nanmean(x)


ibsopt_final = ibsopt;
ibsopt_final.ReturnStd=false;
ibsopt_final.MaxTime=5.0;
ibsopt_final

% increase precision of dt
llopt.dt = 1e-3
inner_fun = @(x,dmat) lca_sim_ibs(x, dmat, llopt)
outer_fun = @(x) ibslike(inner_fun,x,respMat, designMat, ibsopt)

[NLOGL,NLOGLVAR] = deal(nan(length(rep_fval),1));


parfor gg = 1:length(rep_fval)

    if isnan(rep_fval(gg))
        continue;
    end

    [NLOGL(gg),NLOGLVAR(gg),EXITFLAG,OUTPUT] = ibslike(...
        inner_fun, ...
        rep_x(gg,:), ...
        respMat, ...
        designMat, ...
        ibsopt_final);

end


figure;hold on;

nexttile; hold on;
plot(rep_fval, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
plot(NLOGL, '-og', 'MarkerFaceColor', 'g', 'LineWidth', 2)

nexttile; hold on;
plot(center(rep_fval), '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 2)
yline(min(center(rep_fval)), '-k', 'LineWidth', 1)

plot(center(NLOGL), '-og', 'MarkerFaceColor', 'g', 'LineWidth', 2)
yline(min(center(NLOGL)), '-g', 'LineWidth', 1)


[best_orig_fval,best_orig_idx] = nanmin(rep_fval)
[best_final_fval,best_final_idx] = nanmin(NLOGL)

%% MAKE PLOTS ------------
vec = @(x) x(:);



% --------- LOAD -----------------------------------
load('PATH_TO_MODEL')


n_sims = 250

llopt.dt = 1e-4

max_norm = true

fprintf('\n ---- PLOTTING ----')
fprintf('\nmax time: %.2g // best model=%d // max fval=%.4g \n\n', ibs_maxtime, best_fit, best_fval)
disp('best parameters:')
disp(best_x)


plot_x = best_x;
% plot_x = rep_x(best_final_idx,:);




% fix grp_resp
D.acc(isnan(D.rt))=0;% treat misses as errors
D.resp = D.responseID;
[grp_resp,grp_rawRT,grp_RT,grp_subj] = deal([]);
for ov = 0:1
    for ac = 0:1

        Dg = D(D.discOV==ov & D.acc==ac,:);

        grp_subj    = [grp_subj; Dg.subjectID];
        grp_resp    = [grp_resp; Dg.resp];
        grp_rawRT   = [grp_rawRT; Dg.rt];

        quants_ii   = quantile(Dg.rt, linspace(0,1,n_discs+1));
        grp_RT      = [grp_RT; discretize(Dg.rt, quants_ii)];

    end
end
designMat(:,9) = grp_resp;
designMat(:,14) = grp_subj;
designMat(:,15) = grp_RT;
designMat(:,16) = grp_rawRT;








% --------- LOAD -----------------------------------
n_obs = length(D.rt)

[resp, rt, acc] = deal(nan(n_obs,n_sims));

designMat_stack = repmat(designMat, n_sims, 1);

parfor ii = 1:n_sims
    [resp_ii,rt_ii, acc_ii] = lca_sim_ibs(plot_x, designMat, llopt);
    
    resp(:,ii) = resp_ii';
    rt(:,ii) = rt_ii;
    acc(:,ii) = acc_ii;
end

resp = reshape(resp, [n_obs*n_sims,1]);
rt = reshape(rt, [n_obs*n_sims,1]);
acc = reshape(acc, [n_obs*n_sims,1]);

n_sim_obs = length(rt)




% --------- HISTOGRAM -----------------------------------
figure;
cc=1;

prob = @(x) x / nansum(x);
yl=[];
for ov = 0:1
    for ac = 0:1

        sel = D.discOV==ov & D.acc==ac;

        n(cc)=nexttile; hold on;
        % plot([0,quants(cc,2:end)], histcounts(grp_RT(sel,1), 'Normalization','probability'), '-ok', 'MarkerFaceColor', 'k', 'markersize', 8);
        % plot([-1,0,quants(cc,2:end)], histcounts(resp(sel,1), 'Normalization','probability'), '-or', 'MarkerFaceColor', 'r', 'markersize', 8);

        histogram(grp_RT(sel,1),'Normalization','probability');
        histogram(resp(sel,1),'Normalization','probability');
        legend({'obs','fit'})

        yl = [yl;ylim];

        title(sprintf('OV: %d - Accuracy: %d', ov, ac))
        cc=cc+1;
    end
end
for cc = 1:4
ylim(n(cc), [0, max(yl(:,2))]);
end


% --------- RT PPC -----------------------------------
figure;

D = D((D.ExpectedValue_left ~= D.ExpectedValue_right) & isfinite(D.rt),:);

% OV -----------------------------------
nexttile; hold on;
for ov = 0:1
    for ac = 0:1

        Dg = D(D.discOV==ov & D.acc==ac,:);
        scale = height(Dg)/height(D);

        if ac==1
            sgn=1;
        else
            sgn=-1;
        end

        if ov==1
            col='r';
        else
            col='b';
        end

        if max_norm

            % obs
            [f,x] = ksdensity(Dg.rt./Dg.maxRT, 'Bandwidth', 1/60, 'NumPoints', 1000);
            plot(x,scale*sgn*f, '-', 'LineWidth', 2, 'Color', col)

            % fit
            fit_sel = designMat_stack(:,1)==ov & designMat_stack(:,2)==ac;
            sim_rt = rt(fit_sel,:)./designMat_stack(fit_sel,8);
            sim_rt(sim_rt>1)=nan;
            [f,x] = ksdensity(sim_rt, 'Bandwidth', 1/60, 'NumPoints', 1000);
            plot(x,scale*sgn*f, '--', 'LineWidth', 2, 'Color', col);

        else


            % obs
            [f,x] = ksdensity(Dg.rt, 'Bandwidth', 1/60, 'NumPoints', 1000);
            plot(x,scale*sgn*f, '-', 'LineWidth', 2, 'Color', col)

            % fit
            [f,x] = ksdensity(rt(designMat_stack(:,1)==ov & designMat_stack(:,2)==ac,:), 'Bandwidth', 1/60, 'NumPoints', 1000);
            plot(x,scale*sgn*f, '--', 'LineWidth', 2, 'Color', col)

        end

    end
end

title('OV split')
yline(0, '-k', 'LineWidth', 2)
if max_norm
    xline(1, '-k', 'LineWidth', 1)
else
    xline(0.750, '-k', 'LineWidth', 1)
end
set(gca, 'TickDir', 'out', 'LineWidth', 1)




% VD -----------------------------------
nexttile; hold on;
for ov = 0:1
    for ac = 0:1

        Dg = D(D.discVD==ov & D.acc==ac,:);
        scale = height(Dg)/height(D);

        if ac==1
            sgn=1;
        else
            sgn=-1;
        end

        if ov==1
            col='r';
        else
            col='b';
        end

        if max_norm

            % obs
            [f,x] = ksdensity(Dg.rt./Dg.maxRT, 'Bandwidth', 1/60, 'NumPoints', 1000);
            plot(x,scale*sgn*f, '-', 'LineWidth', 2, 'Color', col)

            % fit
            fit_sel = designMat_stack(:,7)==ov & designMat_stack(:,2)==ac;
            sim_rt = rt(fit_sel,:)./designMat_stack(fit_sel,8);
            sim_rt(sim_rt>1)=nan;
            [f,x] = ksdensity(sim_rt, 'Bandwidth', 1/60, 'NumPoints', 1000);
            plot(x,scale*sgn*f, '--', 'LineWidth', 2, 'Color', col);


        else

            % obs
            [f,x] = ksdensity(Dg.rt, 'Bandwidth', 1/60, 'NumPoints', 1000);
            plot(x,scale*sgn*f, '-', 'LineWidth', 2, 'Color', col)

            % fit
            [f,x] = ksdensity(rt(designMat_stack(:,7)==ov & designMat_stack(:,2)==ac,:), 'Bandwidth', 1/60, 'NumPoints', 1000);
            plot(x,scale*sgn*f, '--', 'LineWidth', 2, 'Color', col)

        end


    end
end

title('VD split')
yline(0, '-k', 'LineWidth', 2)
set(gca, 'TickDir', 'out', 'LineWidth', 1)
if max_norm
    xline(1, '-k', 'LineWidth', 1)
else
    xline(0.750, '-k', 'LineWidth', 1)
end



%% create dataset for hddm


% --------- LOAD -----------------------------------
load('PATH_TO_MODEL')


center = @(x) x-nanmean(x);

version     = 'compare' % ddm, lca, compare
write_ddm   = true


% --- fix grp_resp
D.acc(isnan(D.rt))=0;% treat misses as errors
D.resp = D.responseID;
[grp_resp,grp_rawRT,grp_RT,grp_subj] = deal([]);
for ov = 0:1
    for ac = 0:1

        Dg = D(D.discOV==ov & D.acc==ac,:);

        grp_subj    = [grp_subj; Dg.subjectID];
        grp_resp    = [grp_resp; Dg.resp];
        grp_rawRT   = [grp_rawRT; Dg.rt];

        quants_ii   = quantile(Dg.rt, linspace(0,1,n_discs+1));
        grp_RT      = [grp_RT; discretize(Dg.rt, quants_ii)];

    end
end
designMat(:,9) = grp_resp;
designMat(:,14) = grp_subj;
designMat(:,15) = grp_RT;
designMat(:,16) = grp_rawRT;

% params ---------
llopt.dt = 1e-4
switch version
    case 'lca'
        n_sims = 3
    case 'ddm'
        n_sims = 1
    case 'compare'
        n_sims = 100
end

fprintf('\n ---- GENERATING HDDM ----')
fprintf('\nmax time: %.2g // best model=%d // max fval=%.4g \n\n', ibs_maxtime, best_fit, best_fval)
disp('best parameters:')
disp(best_x)


plot_x = best_x;
% plot_x = rep_x(best_final_idx,:);


% --------- GENERATE DATA -----------------------------------
n_obs = length(D.rt)


designMat_stack = repmat(designMat, n_sims, 1);
respMat_stack = repmat(respMat, n_sims, 1);

[resp, rt, acc] = deal(nan(n_obs,n_sims));

parfor ii = 1:n_sims
    [resp_ii,rt_ii, acc_ii] = lca_sim_ibs(plot_x, designMat, llopt);
    
    resp(:,ii) = resp_ii';
    rt(:,ii) = rt_ii;
    acc(:,ii) = acc_ii;
end

resp = reshape(resp, [n_obs*n_sims,1]);
rt = reshape(rt, [n_obs*n_sims,1]);
acc = reshape(acc, [n_obs*n_sims,1]);


n_sim_obs = length(rt)




% --------- FORMAT TABLE (RESPONSE CODED) -----------------------------------


tout = table;

tout.subj_idx   = ones(size(rt));
tout.sim_acc   = double(acc==1); % accuracy-coded
tout.rt         = rt;


tout.orig_acc   =  designMat_stack(:,2);
tout.orig_resp   = double(designMat_stack(:,9)==1);
tout.orig_resp(tout.orig_resp==0)   = -1;


tout.orig_subj   = designMat_stack(:,14);
tout.orig_discRT = designMat_stack(:,15);
tout.orig_rt     = designMat_stack(:,16);

% designMat = [...
%     1 grp_OV, 
%     2 grp_acc,...
%     3 grp_CorrMaxVal/80, 
%     4 grp_CorrMinVal/80, 
%     5 grp_ErrMaxVal/80, 
%     6 grp_ErrMinVal/80, ...
%     7 grp_VD, 
%     8 grp_maxRT...
%     ];

A_max  = designMat_stack(:,10);
A_min = designMat_stack(:,11);
A_ev    = A_max + A_min;

B_max  = designMat_stack(:,12);
B_min = designMat_stack(:,13);
B_ev    = B_max + B_min;


tout.signVD     = (A_ev - B_ev)./std(A_ev - B_ev);
tout.maxVD      = (A_max -B_max)./std(A_max - B_max);
tout.minVD      = (A_min - B_min)./std(A_min - B_min);

tout.absVD      = abs(A_ev - B_ev);
tout.absMaxVD   = abs(A_max - B_max);
tout.absMinVD   = abs(A_min - B_min);

tout.OV         = A_ev + B_ev;
tout.maxOV      = A_max + B_max;
tout.minOV      = A_min + B_min;


sign_acc = tout.sim_acc;
sign_acc(sign_acc==0) = -1;
tout.response = sign_acc .* sign(tout.signVD); % response coded




switch version

    case 'lca'
        tout = tout(isfinite(tout.rt), :);

    case 'ddm'
        tout = tout(isfinite(tout.orig_rt), :);

        % set behaviour to human
        tout.rt = tout.orig_rt;
        tout.response = tout.orig_resp;
        tout.sim_acc = tout.orig_acc;

    case 'compare'
        tout = tout(isfinite(tout.orig_rt), :);

end


pts = unique(tout.subj_idx);
for pp = 1:length(pts)
    
    tout.absVD(tout.subj_idx == pts(pp)) = zscore(tout.absVD(tout.subj_idx == pts(pp)));
    tout.absMaxVD(tout.subj_idx == pts(pp)) = zscore(tout.absMaxVD(tout.subj_idx == pts(pp)));
    tout.absMinVD(tout.subj_idx == pts(pp)) = zscore(tout.absMinVD(tout.subj_idx == pts(pp)));
    
    tout.OV(tout.subj_idx == pts(pp)) = zscore(tout.OV(tout.subj_idx == pts(pp)));
    tout.maxOV(tout.subj_idx == pts(pp)) = zscore(tout.maxOV(tout.subj_idx == pts(pp)));
    tout.minOV(tout.subj_idx == pts(pp)) = zscore(tout.minOV(tout.subj_idx == pts(pp)));

end

tout.sqrAbsVD       = tout.absVD.^2;
tout.sqrAbsMaxVD    = tout.absMaxVD.^2;
tout.sqrAbsMinVD    = tout.absMinVD.^2;


head(tout)
grpstats(tout,'subj_idx')


model_acc = corr(tout.response, tout.signVD)
human_acc = corr(tout.orig_resp, tout.signVD)




% ----- PLOT VD ----------
t_plt = tout(tout.signVD~=0 & isfinite(tout.orig_rt),:);
nexttile; hold on;
for ov = 0:1
    for ac = 0:1

        sel = ((t_plt.absVD>median(t_plt.absVD)) == ov) & t_plt.sim_acc== ac;
        scale = mean(sel);

        if ac==1
            sgn=1;
        else
            sgn=-1;
        end

        if ov==1
            col='r';
        else
            col='b';
        end


        % obs
        [f,x] = ksdensity(t_plt.orig_rt(sel,:), 'Bandwidth', 1/60, 'NumPoints', 1000);
        plot(x,scale*sgn*f, '-', 'LineWidth', 2, 'Color', col)

        % fit
        [f,x] = ksdensity(t_plt.rt(sel,:), 'Bandwidth', 1/60, 'NumPoints', 1000);
        plot(x,scale*sgn*f, '--', 'LineWidth', 2, 'Color', col)

        xline(nanmedian(t_plt.rt(sel,:)), 'Color', col)

    end
end

title('OV split')
yline(0, '-k', 'LineWidth', 2)
set(gca, 'TickDir', 'out', 'LineWidth', 1)


% ensure internally consistent
sel_resp = (tout.signVD~=0) & isfinite(tout.orig_rt);
assert(all(tout.sim_acc(sel_resp) == (tout.response(sel_resp) == sign(tout.signVD(sel_resp)))))
assert(all(tout.orig_acc(sel_resp) == (tout.orig_resp(sel_resp) == sign(tout.signVD(sel_resp)))))
assert(~all(tout.response(sel_resp) == sign(tout.signVD(sel_resp))))


if write_ddm

    switch version
        case 'lca'
            writetable(tout, sprintf('ddm-sims/%s_%s_fitLCA.csv', time_on, savename));
        case 'ddm'
            writetable(tout, sprintf('ddm-sims/%s_%s_fitDDM.csv', time_on, savename));
        case 'compare'
            writetable(tout, sprintf('ddm-sims/%s_%s_compare.csv', time_on, savename));
    end

    disp('saved!')

end






%% compare DDM and LCA fits (BIC)

% LOAD DATA

lca = readtable('LCA_FIT_PATH')

ddm_cell = cell(0);
ddm_cell{1} = readtable('DDM_FIT_PATH_ORIGINAL*')
ddm_cell{2} = readtable('DDM_FIT_PATH_VDOV-COLLAPSE')



ddm_name = {'Original*', 'VDOV_{collapse}'}

ddm_params = [5,14]





%% Run BIC comparison


n_disc = 10


% loop over iterations

assert(height(lca) == height(ddm_cell{1}))

n_iter = 100;
n_obs = height(lca)/n_iter;
sim_sel = repmat(1:n_iter, [n_obs,1]);
sim_sel = sim_sel(:);



[chi_lca_iter, G2_lca_iter, BIC_lca_iter,...
    chi_ddm_iter, G2_ddm_iter, BIC_ddm_iter,...
    chi_ddmVD_iter, G2_ddmVD_iter, BIC_ddmVD_iter] = deal(nan(n_iter,1));

[chi_ddm_iter,G2_ddm_iter,BIC_ddm_iter] = deal(cell(length(ddm_cell),1));
for dd = 1:length(ddm_cell)
    chi_ddm_iter{dd} = deal(nan(n_iter,1));
    G2_ddm_iter{dd} = deal(nan(n_iter,1));
    BIC_ddm_iter{dd} = deal(nan(n_iter,1));
end




for ii = 1:n_iter


    % get lca fit ------------------------------------------------


    lca_iter = lca(sim_sel==ii,:);

    split_VD = lca_iter.absVD > median(lca_iter.absVD);
    split_acc_lca = lca_iter.sim_acc;
    split_acc_hum = lca_iter.orig_acc;

    chi_lca = 0;
    G2_lca  = 0;
    BIC_lca = 0;

    for vv = 0:1
        for aa = 0:1

            lca_sel = split_VD==vv & split_acc_lca==aa;
            hum_sel = split_VD==vv & split_acc_hum==aa;

            hum_quants = quantile(lca_iter.orig_rt(hum_sel), linspace(0,1,n_disc+1));

            lca_discs = discretize(lca_iter.rt(lca_sel), hum_quants);
            hum_discs = discretize(lca_iter.orig_rt(hum_sel), hum_quants);

            for dd = 1:n_disc

                lca_P = nanmean(lca_discs==dd);
                hum_P = nanmean(hum_discs==dd);

                chi_lca = chi_lca + (n_obs*(hum_P - lca_P).^2)/lca_P;
                G2_lca = G2_lca + n_obs*hum_P*(log(hum_P) - log(lca_P));
                BIC_lca = BIC_lca + n_obs*hum_P*log(lca_P);

            end


        end
    end

    G2_lca = 2*G2_lca;
    BIC_lca = -2*BIC_lca + 8*log(n_obs);


    chi_lca_iter(ii) = chi_lca;
    G2_lca_iter(ii) = G2_lca;
    BIC_lca_iter(ii) = BIC_lca;







    % get ddm fit ------------------------------------------------
    for mm = 1:length(ddm_cell)

        ddm = ddm_cell{mm};

        ddm_iter = ddm(ddm.draw == (ii-1),:);

        split_VD = ddm_iter.absVD > median(ddm_iter.absVD);
        split_acc_ddm = ddm_iter.response_sampled == sign(ddm_iter.signVD);
        split_acc_hum = ddm_iter.response == sign(ddm_iter.signVD);

        chi_ddm = 0;
        G2_ddm = 0;
        BIC_ddm = 0;

        for vv = 0:1
            for aa = 0:1

                ddm_sel = split_VD==vv & split_acc_ddm==aa;
                hum_sel = split_VD==vv & split_acc_hum==aa;

                hum_quants = quantile(ddm_iter.rt(hum_sel), linspace(0,1,n_disc+1));

                ddm_discs = discretize(ddm_iter.rt_sampled(ddm_sel), hum_quants);
                hum_discs = discretize(ddm_iter.rt(hum_sel), hum_quants);

                for dd = 1:n_disc

                    ddm_P = nanmean(ddm_discs==dd);
                    hum_P = nanmean(hum_discs==dd);

                    chi_ddm = chi_ddm + (n_obs*(hum_P - ddm_P).^2)/ddm_P;
                    G2_ddm = G2_ddm + n_obs*hum_P*(log(hum_P) - log(ddm_P));
                    BIC_ddm = BIC_ddm + n_obs*hum_P*log(ddm_P);

                end





            end
        end


        G2_ddm = 2*G2_ddm;
        BIC_ddm = -2*BIC_ddm + ddm_params(mm)*log(n_obs);



        chi_ddm_iter{mm}(ii) = chi_ddm;
        G2_ddm_iter{mm}(ii) = G2_ddm;
        BIC_ddm_iter{mm}(ii) = BIC_ddm;


    end





end




% plot
figure; hold on;


b=80
[f,x,b] = ksdensity(BIC_lca_iter, 'NumPoints',1000, 'Bandwidth',b);
plt = [];
plt(1) = plot(x,f, '-g', 'LineWidth', 2);
% xline(mean(BIC_lca_iter), '-g', 'LineWidth', 1);

ddm_ls = {'-r', '-b'};
for dd = 1:length(ddm_cell)


    [f,x] = ksdensity(BIC_ddm_iter{dd}, 'NumPoints',1000, 'Bandwidth',b);
    plt(1+dd)=plot(x,f, ddm_ls{dd}, 'LineWidth', 2);

    mean_BIC = mean(BIC_ddm_iter{dd})
    % xline(mean_BIC, ddm_ls{dd}, 'LineWidth', 1);


end



set(gca, 'TickDir', 'out', 'LineWidth', 1)
title('LCA vs VDev')
xlabel('BIC (smaller is better)')
legend(plt,[{'lca'}, ddm_name])


LCA_VS_DDM = fix(mean(BIC_lca_iter) - mean(BIC_ddm_iter{2}))
LCA_VS_MAX_DDM = fix(min(BIC_lca_iter) - min(BIC_ddm_iter{2}))
mean(BIC_lca_iter < min(BIC_ddm_iter{2}))
mean(mean(BIC_lca_iter) > BIC_ddm_iter{2})

n_boot = 10000;
diff_boot = nan(n_boot,1);

for ii = 1:n_boot
    diff_boot(ii) = mean(BIC_lca_iter(randi(100,[100,1]))) - mean(BIC_ddm_iter{2}(randi(100,[100,1])));
end

boot_prct = prctile(diff_boot, [2.5, 50, 97.5])
max(boot_prct) - min(boot_prct)

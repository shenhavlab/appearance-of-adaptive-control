%% Run revovery

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

vec = @(x) x(:)


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





%% MODEL RECOVERY

% nn-based model recovery ---------------
fld = 'nn_summary'
model_name = '24-07-04_NN_'
models = {...
    'VDev_subj',...
    'VD_static_subj',...
    'VD_collapse_subj',...
    }


% standard model recovery ---------------
% fld = 'summary'
% model_name = ''
% 
% models = {...
%     'VDev_subj',...
%     'VD_static_subj',...
%     'VD_collapse_subj',...
%     }

% models = {...
%     'VDev_subj',...
%     'VD_static_subj',...
%     'VD_collapse_subj',...
%     'VDOV_collapseA_subj',...
%     }

n_models = length(models)

DIC = nan(n_models, n_models, 50);
sample_num = nan(n_models, n_models, 50);




for ff = 1:n_models % fits

    for gg = 1:n_models % gens



        wild = sprintf('gen-%s__fit-%s*', models{gg}, models{ff})

        fn = dir(sprintf('../recovery-local/%s/%s%s/%s',fld, model_name, wild, wild));


        for ww = 1:length(fn)

            fit = readtable(fullfile(fn(ww).folder, fn(ww).name));


            name = split(fn(ww).name, {'-','_'});
            sample_num(ff,gg,ww) = str2double(name{end-1});

            DIC(ff,gg,sample_num(ff,gg,ww)) = fit{end,2};

        end


    end


end






% plot

f=figure;
set(f,'defaultAxesTickLabelInterpreter','none')



% DIC vs diag
dics_mx = nanmean(DIC(:,:,:),3) - diag(nanmean(DIC(:,:,:),3))';

n1=nexttile;

imagesc(dics_mx, [-100, 100]); colorbar
colormap(n1,vik)
axis('square')


xticks(1:n_models)
xticklabels(models)
xlabel('generating')

yticks(1:n_models)
yticklabels(models)
ylabel('fit')
title('\DeltaDIC (relative to diagonal)')


% DIC observed vs null
dic_bias = vec(DIC(3,2,:)) - vec(DIC(2,2,:));
obs_diff = fix(dic_VD_collapse_subj - dic_VD_static_subj)
[min(dic_bias), nanmedian(dic_bias), max(dic_bias)]

n2=nexttile; hold on;
[f,x]=ksdensity(dic_bias);
plot(x,f, '-k', 'LineWidth', 2)
xline(dic_VD_collapse_subj - dic_VD_static_subj, '-g', 'LineWidth', 2)
xline(0, '-k', 'LineWidth', 2)

set(gca, 'TickDir', 'out', 'LineWidth', 1)
xlabel('\Delta DIC')
title('Observed \DeltaDIC much larger than bias')
yticks([])
xlim([-1250, 50])

% 
% % DIC vs min
% dics_mx = nanmean(DIC(:,:,:),3) - nanmin(nanmean(DIC(:,:,:),3));
% 
% n1=nexttile;
% 
% imagesc(dics_mx, [-100, 100]); colorbar
% colormap(n1,vik)
% axis('square')
% 
% 
% xticks(1:n_models)
% xticklabels(models)
% xlabel('generating')
% 
% yticks(1:n_models)
% yticklabels(models)
% ylabel('fit')
% title('\Delta DIC (relative to min)')
% 
% 
% 
% 
% 
% 
% 
% % number of simulations
% n2=nexttile;
% 
% imagesc(sum(isfinite(DIC),3)); colorbar;
% colormap(n2, 'hot')
% axis('square')
% 
% xticks(1:n_models)
% xticklabels(models)
% xlabel('generating')
% 
% yticks(1:n_models)
% yticklabels(models)
% ylabel('fit')
% title('N models')


%% plot per sample
f15=figure;
set(f15,'defaultAxesTickLabelInterpreter','none')

for ss = 1:size(DIC,3)
    
    dics_mx = nanmean(DIC(:,:,ss),3) - diag(nanmean(DIC(:,:,ss),3))';
    
    % DIC
    n1=nexttile;
    
    imagesc(dics_mx, [-100, 100]); colorbar
    colormap(n1,vik)
    
    
    xticks(1:n_models)
    xticklabels(models)
    xlabel('generating')
    
    yticks(1:n_models)
    yticklabels(models)
    ylabel('fit')
    title(['\Delta DIC ' num2str(ss)])
end





%% long-form model recovery plotting



% histogram


% histogram
f2=figure;
set(f2,'defaultAxesTickLabelInterpreter','none')


dic_l = reshape(DIC, [n_models*n_models, size(DIC,3)]);

for ii = 1:size(dic_l,1)
    nexttile;
    ksdensity(dic_l(ii,:)')

end




% plot long
f3=figure;
set(f3,'defaultAxesTickLabelInterpreter','none')


% fit
nf=nexttile;
dic_fit = imagesc(reshape(DIC, [n_models, n_models*50]));
colormap(nf,'hot')
title('fit')
yticks(1:n_models)
yticklabels(models)
ylabel('fit')
xlabel('gen x samples')


nf=nexttile;
dic_fit = imagesc(isfinite(reshape(DIC, [n_models, n_models*50])));
colormap(nf,'hot')
title('fit finite')
yticks(1:n_models)
yticklabels(models)
ylabel('fit')
xlabel('gen x samples')


% gen

ng=nexttile;
dic_fit = imagesc(reshape(permute(DIC,[2,1,3]), [n_models, n_models*50]));
colormap(ng,'hot')
title('generating')
yticks(1:n_models)
yticklabels(models)
ylabel('generating')
xlabel('fit x samples')

ng=nexttile;
dic_fit = imagesc(isfinite(reshape(permute(DIC,[2,1,3]), [n_models, n_models*50])));
colormap(ng,'hot')
title('generating finite')
yticks(1:n_models)
yticklabels(models)
ylabel('generating')
xlabel('fit x samples')


% sample
ng=nexttile;
dic_fit = imagesc(reshape(permute(DIC,[3,1,2]), [50, n_models*n_models]));
colormap(ng,'hot')
title('sample')
ylabel('sample')
xlabel('gen x fit')

ng=nexttile;
dic_fit = imagesc(isfinite(reshape(permute(DIC,[3,1,2]), [50, n_models*n_models])));
colormap(ng,'hot')
title('sample finite')
ylabel('sample')
xlabel('gen x fit')








%% check recovery issues

t = readtable('PATH_TO_VDev_subj_ppc.csv');
disp('loaded')

%% plot










rt_range = nan(50, 2);
perf_mean = nan(50,3);


figure;


for ss = 1:50

    sel = t.sample == ss;

    rt = t.rt_sampled(sel);
    resp = t.response_sampled(sel);
    acc = sign(resp) == sign(t.signVD(sel));
    macc = mean(acc);

    [f1,x1] = ksdensity(rt(acc==1));
    [f2,x2] = ksdensity(rt(acc==0));

    nexttile; hold on;
    plot(x1, f1*macc, '-g', 'LineWidth', 1.5);
    plot(x2, -f1*(1-macc), '-r', 'LineWidth', 1.5);
    yline(0, '-k')
    xlim([0,2])
    title(ss)


end




% plot subject means

sj = unique(t.subj_idx);
nsj = length(sj);


rt_subj = nan(50,nsj,6);
resp_subj = nan(50,nsj);
acc_subj = nan(50,nsj);

for pp = 1:nsj


    for ss = 1:50


        sel = (t.sample == ss) & (t.subj_idx == sj(pp));

        rt = t.rt_sampled(sel);
        resp = t.response_sampled(sel);
        acc = sign(resp) == sign(t.signVD(sel));

        rt_subj(ss,pp,:) = [...
            min(rt), max(rt), ...
            min(rt(acc==1)), max(rt(acc==1)),...
            min(rt(acc==0)), max(rt(acc==0)),...
            ];

        resp_subj(ss,pp) = mean(resp==1);
        acc_subj(ss,pp) = mean(acc);


    end

end


% response / rt
figure;

n1=nexttile;
imagesc(resp_subj, [0,1]); colorbar
colormap(n1, vik)
title('response bias')
xlabel('subjects')
ylabel('sims')

n2=nexttile;
imagesc(acc_subj, [0,1]);  colorbar
colormap(n2, vik)
title('accuracy')
xlabel('subjects')
ylabel('sims')


fprintf('\nmin resp = %f, max resp = %f\n', min(resp_subj(:)), max(resp_subj(:)) )
fprintf('min acc = %f, max acc = %f\n', min(acc_subj(:)), max(acc_subj(:)) )



figure;
for ss = 1:nsj

    nn=nexttile;
    imagesc(squeeze(rt_subj(:,ss,:)), [0,1])
    colormap(nn, batlow)

end


min(rt_subj(:))
max(rt_subj(:))


%% subj rt distributions


% plot subject means

sj = unique(t.subj_idx);
nsj = length(sj);



for pp = 1:nsj

    figure;


    for ss = 1:50


        sel = (t.sample == ss) & (t.subj_idx == sj(pp));

        rt = t.rt_sampled(sel);
        resp = t.response_sampled(sel);
        acc = sign(resp) == sign(t.signVD(sel));
        macc = mean(acc);

        [f1,x1] = ksdensity(rt(acc==1));
        [f2,x2] = ksdensity(rt(acc==0));

        nexttile; hold on;
        plot(x1, f1*macc, '-g', 'LineWidth', 1.5);
        plot(x2, -f1*(1-macc), '-r', 'LineWidth', 1.5);
        yline(0, '-k')
        xlim([0,2])
        title(sprintf('%d-%d', pp, ss))


    end

end









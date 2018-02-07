clear all
%org_data = load('E:\Box Sync\skinner\projects_analyses\Project Trust\data\model-derived\scan\actual\202200_cntr0_mltrun1_kappa2_censor0');
%org_data = load('E:\Box Sync\skinner\projects_analyses\Project Trust\data\model-derived\scan\f_trust_SVM1\202200_cntr0_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0');
%org_data = load('E:\Box Sync\skinner\projects_analyses\Project Trust\data\model-derived\scan\policy\202200_cntr2_mltrun1_kappa2_censor0');
org_data = load('E:\Box Sync\skinner\projects_analyses\Project Trust\data\model-derived\scan\regret\202200_cntr1_mltrun1_kappa2_censor0');


%Just use original u,y,f(g)_fname, options
y = org_data.y;
u = org_data.u;
f_fname = org_data.out.options.f_fname;
g_fname = org_data.out.options.g_fname;
dim = org_data.dim;
options = org_data.options;
[post_new,out_new] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);



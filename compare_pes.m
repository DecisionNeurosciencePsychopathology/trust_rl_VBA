function compare_pes
%Load original PEs
paper_pes=load('E:\trust_model_comparision\trust_rl_VBA\f_trust_Qlearn_original_PEs_paper\PEs\modelPEs_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');

%Load new kappa parameter model PEs
new_pes = load('E:\trust_model_comparision\trust_rl_VBA\f_trust_Qlearn_policy_censor0_mltrun1_kS_kT\PEs\modelPEs_f_trust_Qlearn_policy_censor0_mltrun1_kS_kT.mat');

%Create new plots and state correlation
for i = 1:length(new_pes.M)
   figure(i)
   clf;
   if paper_pes.M(i,1)==new_pes.M(i,1)
       %Plot
   plot(paper_pes.M(i,2:end),'b')
   hold on
   plot(new_pes.M(i,2:end),'r')
   
   %t-Test
   norm_original=norm_me(paper_pes.M(i,2:end));
   norm_new=norm_me(new_pes.M(i,2:end));
   h=ttest(norm_original,norm_new);
   
   rho=corr(paper_pes.M(1,2:end)',new_pes.M(1,2:end)');
   legend({'Paper pes','kappa2 pes'})
   title([num2str(paper_pes.M(i,1)) ' original PEs vs new PEs, ttest=' num2str(h)])
   if i==10
       waitforbuttonpress
   end
   else
       warning('Mismatch with index number %d',i)
       continue
   end
end

function norm_data=norm_me(bla)
norm_data = (bla - min(bla)) / ( max(bla) - min(bla) );
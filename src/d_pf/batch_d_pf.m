function batch_d_pf(x0,N_OL,N_CL,p,f,k_p,num_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,run_dir,Obsv_pf,Lu_pf,Ly_pf,Gu_pf,name_kvar)
system(append("echo p=,",num2str(p),', entered batch script'));
parfor k_e = 1:num_e %(k_e = 1:num_e, myCluster)
    seed_num = (k_p-1)*num_e+k_e+2520;
    loop_var(x0,N_OL,N_CL,p,f,k_p,k_e,plant,Ru,Re,ny,nu,nx,num_steps,Nbar,ref,Qk,Rk,dRk,num_c,Rdu,CL_sim_steps,run_dir,seed_num,Obsv_pf,Lu_pf,Ly_pf,Gu_pf,name_kvar);
end
end
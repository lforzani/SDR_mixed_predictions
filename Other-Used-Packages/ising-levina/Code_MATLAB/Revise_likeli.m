


%%%%%%%%%%%%%%%%%%%% Revise Simulation Scale %%%%%%%%%%%%%%%%%%%%%
n = 200;
q = 10;
p = 20;
Val = 2.^linspace(-1,6,8);
for j = 1:8
    j
name = sprintf('Value=%1.1f_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
load(name)
name = sprintf('ValueIdv=%1.1f_n=%d_q=%d_p=%d.mat',Val(j),n,q,p);
load(name);
[LambdaBest_jnt likeli_jnt]= Validtune(data_y, data_x, gamma_mat);
gamma_mat_max = select_gamma(gamma_mat_idv);
        [LambdaBest_idv likeli_idv]= Validtune(data_y, data_x, gamma_mat_max);
        for k = 1:K
            distIdv_1(:,k) = DistL2(gamma, gamma_mat_max(:,:,:,:,k), 1);
            distIdv_0(:,k) = DistL2(gamma, gamma_mat_max(:,:,:,:,k), 0);
        end
save(name,'gamma', 'gamma_mat_idv', 'lambda', 'sens_min', 'spec_min','sens_max', 'spec_max', ...
                    'likeli_jnt', 'likeli_idv','distIdv_1', 'distIdv_0', 'SensBest_jnt','SpecBest_jnt','SensBest_max','SpecBest_max');
end

%%%%%%%%%%%%%%%%%%%% Revise Simulation Sparsity %%%%%%%%%%%%%%%%%%%%%%%%%
Prop = [0.2 0.5 0.8];
nEdge = [10 20 30];
for i = 1:3     
        for j = 1:2
        edge = nEdge(i)    
        prop = Prop(j)
        name = sprintf('Sparse_prop=%.1f_nedge=%d.mat',prop, edge);
        load(name);
        name = sprintf('SparseIdv_prop=%.1f_nedge=%d.mat',prop, edge);
        load(name);
        [LambdaBest_jnt likeli_jnt]= Validtune(data_y, data_x, gamma_mat);
        gamma_mat_max = select_gamma(gamma_mat_idv);
        [LambdaBest_idv likeli_idv]= Validtune(data_y, data_x, gamma_mat_max);
        for k = 1:K
            distIdv_1(:,k) = DistL2(gamma, gamma_mat_max(:,:,:,:,k), 1);
            distIdv_0(:,k) = DistL2(gamma, gamma_mat_max(:,:,:,:,k), 0);
        end
        save(name,'gamma', 'gamma_mat_idv','lambda', 'sens_min', 'spec_min','sens_max', 'spec_max', ...
                    'likeli_jnt', 'likeli_idv','distIdv_1', 'distIdv_0',...
                    'SensBest_jnt','SpecBest_jnt','SensBest_max','SpecBest_max');
        end
end






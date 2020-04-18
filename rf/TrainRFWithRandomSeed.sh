#Train Models 

Rscript TrainRFWithRandomSeed.R --features_mat xin_feature_matrices/neil1_computational.features.csv.cleaned --output_dir xin_neil1 --n_iter 10

Rscript TrainRFWithRandomSeed.R --features_mat xin_feature_matrices/ajuba_bc_computational.features.csv.cleaned --output_dir xin_ajuba --n_iter 10

Rscript TrainRFWithRandomSeed.R --features_mat xin_feature_matrices/ttyh2_bc_computational.features.csv.cleaned --output_dir xin_ttyh2_bc_computational --n_iter 10

Rscript TrainRFWithRandomSeed.R --features_mat xin_feature_matrices/ttyh2_ecs_bc_computational.features.csv.cleaned --output_dir xin_ttyh2_ecs_bc_computational  --n_iter 10

Rscript TrainRFWithRandomSeed.R --features_mat xin_feature_matrices/ttyh2_computational.features.csv.cleaned --output_dir xin_ttyh2_computational  --n_iter 10

#Analyze Feature Correlations
Rscript feature_cor.R --feature_mat xin_feature_matrices/neil1_computational.features.csv.cleaned --output_prefix xin_neil1/neil1
Rscript feature_cor.R --feature_mat xin_feature_matrices/ajuba_bc_computational.features.csv.cleaned --output_prefix xin_ajuba/ajuba
Rscript feature_cor.R --feature_mat xin_feature_matrices/ttyh2_bc_computational.features.csv.cleaned --output_prefix xin_ttyh2_bc_computational/ttyh2_bc
Rscript feature_cor.R --feature_mat xin_feature_matrices/ttyh2_ecs_bc_computational.features.csv.cleaned --output_prefix xin_ttyh2_ecs_bc_computational/ttyh2_ecs_bc
Rscript feature_cor.R --feature_mat xin_feature_matrices/ttyh2_computational.features.csv.cleaned --output_prefix xin_ttyh2_computational/ttyh2_computational



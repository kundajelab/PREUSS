echo "NEIL1"
python generate_feature_matrix_for_rf_expanded.py --rna_lib_structure_summary_json /srv/scratch/annashch/adar_editing/data_from_xin/2018_07_30-Neil1_Structure_Summary/rna_lib-structure-summary-Neil-TG_1.json \
       --bpRNA_pickle /srv/scratch/annashch/adar_editing/bpRNA_wrapper/NEIL1.bpRNA.pkl \
       --outf NEIL1.features.txt \
       --approach computational
echo "DONE NEIL1"

echo "NEIL1 SHAPEMAPPER" 
python generate_feature_matrix_for_rf_expanded.py --rna_lib_structure_summary_json /srv/scratch/annashch/adar_editing/data_from_xin/2018_10_15-Neil1_Combined-Structure_Summary/rna_lib-structure_summary.json \
       --bpRNA_pickle /srv/scratch/annashch/adar_editing/bpRNA_wrapper/NEIL1.SHAPEMAPPER2.bpRNA.pkl \
       --outf NEIL1.SHAPEMAPPER2.features.txt \
       --approach computational
echo "DONE NEIL1 SHAPEMAPPER"

echo "TTYH2.BC" 
python generate_feature_matrix_for_rf_expanded.py --rna_lib_structure_summary_json /srv/scratch/annashch/adar_editing/data_from_xin/2018_11_13-TTYH2_BC-Structure_Summary-Reference_Only/rna_lib-structure_summary.json \
       --bpRNA_pickle /srv/scratch/annashch/adar_editing/bpRNA_wrapper/TTYH2.BC.bpRNA.pkl \
       --outf TTYH2.BC.features.txt \
       --approach computational
echo "DONE TTYH2.BC"

echo "TTYH2.ECS"
python generate_feature_matrix_for_rf_expanded.py --rna_lib_structure_summary_json /srv/scratch/annashch/adar_editing/data_from_xin/2018_11_13-TTYH2_ECS-Structure_Summary-Reference_Only/rna_lib-structure_summary.json \
       --bpRNA_pickle /srv/scratch/annashch/adar_editing/bpRNA_wrapper/TTYH2.ECS.bpRNA.pkl \
       --outf TTYH2.ECS.features.txt \
       --approach computational
echo "DONE TTYH2.ECS"

echo "AJUBA" 
python generate_feature_matrix_for_rf_expanded.py --rna_lib_structure_summary_json /srv/scratch/annashch/adar_editing/data_from_xin/2018_11_13-AJUBA_BC-Structure_Summary-Reference_Only/rna_lib-structure_summary.json \
       --bpRNA_pickle /srv/scratch/annashch/adar_editing/bpRNA_wrapper/AJUBA.bpRNA.pkl \
       --outf AJUBA.features.txt \
       --approach computational
echo "DONE AJUBA" 

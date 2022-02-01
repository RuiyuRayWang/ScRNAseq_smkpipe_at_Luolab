import pandas as pd

def wash_whitelist(path_to_whitelist,
                   path_to_ground_truth):

    # Read files
    whitelist = pd.read_csv(path_to_whitelist, sep='\t', names=["Barcode","Corrected","Count_BC","Count_Corrected"], 
                            header=None)
    BC_ground_truth_raw = pd.read_csv(path_to_ground_truth)
    
    # Parse ground truth list
    BC_ground_truth = BC_ground_truth_raw['Primer_sequence'].str.extract(r'TCAGACGTGTGCTCTTCCGATCT([ATCG]{8})', expand=False)
    
    # Keep entries whose Well BC and Ab BC are in their respective ground truth lists
    whitelist_washed = whitelist.loc[whitelist['Barcode'].isin(BC_ground_truth)]

    return(whitelist_washed)

whitelist_washed = wash_whitelist(path_to_whitelist = snakemake.input[0],
                                  path_to_ground_truth = snakemake.input[1])

whitelist_washed.to_csv(snakemake.output[0], index=False, sep="\t", header=False)

# # For debug
# whitelist = pd.read_csv('~/GITHUB_REPO/BLISAcounts/workflow/data/zhongshilin/210903_E00516_0748_BHG3T7CCX2/outs/AR005_S21_L006_whitelist.txt',
#                         sep='\t', names=["Barcode","Corrected","Count_BC","Count_Corrected"], header=None)
# well_ground_truth = pd.read_csv('~/GITHUB_REPO/BLISAcounts/resources/well_ground_truth.csv', names=['Well_BC'], header=None)
# ab_ground_truth = pd.read_csv('~/GITHUB_REPO/BLISAcounts/resources/ab_ground_truth.csv', names=['Ab_BC'], header=None)

# whitelist_washed.to_csv('~/GITHUB_REPO/BLISAcounts/workflow/data/zhongshilin/210903_E00516_0748_BHG3T7CCX2/outs/AR005_S21_L006_whitelist_washed.txt', index=False, sep="\t", header=False)
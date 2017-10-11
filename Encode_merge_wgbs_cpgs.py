import pandas as pd
import numpy as np
import os
import glob
import pickle
from os.path import join
from os.path import splitext
from os.path import basename

output_dir =  "/gpfs/gpfs1/home/schhetri/for_chris/batch_I"
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

cpg_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL*.cov"
cpg_file_list =  glob.glob(cpg_file_dir)


""" \nProcessing and creating a Cpg bed file with a strand cols from raw cpg file 
followed by pybedtool object formation...\n """
print " \nProcessing Cpg bed file\n "

for cpg_bed_file in cpg_file_list:
	os.environ["output_dir"] = output_dir
	os.environ["cpg_bed_file"] = cpg_bed_file
	cpg_bed_edited = splitext(basename(cpg_bed_file))[0] + "_edited" + splitext(cpg_bed_file)[1]
	print cpg_bed_edited
	os.environ["cpg_bed_edited"] = cpg_bed_edited

	if not os.path.exists(join(output_dir,cpg_bed_edited)):
	    CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3+1,$5,$6,"." }' $cpg_bed_file | sort -k1,1 -k2,2n | grep -v "^*" | grep -v "chrM" > $output_dir/$cpg_bed_edited'''
	    os.system(CMD) #subprocess.call(COMMAND, shell=True) # import subprocess

""" Find the intersect or common cpgs for merging the methylated reads and unmeth reads"""
os.system("/opt/BEDTools-2.19.1/bin/intersectBed -a SL83768_CpG_context_deduplicated.bismark_edited.cov -b SL83771_CpG_context_deduplicated.bismark_edited.cov -wa -wb > SL83768_SL83771_CpG_context_deduplicated.bismark.cov")

intersected_cpg_file = join(output_dir, "SL83768_SL83771_CpG_context_deduplicated.bismark_edited.cov" )

""" Read the intersect cpg file for adding up meth and unmeth reads """
combined_cpg_df =  pd.read_csv(intersected_cpg_file, sep="\t", header=None)
combined_cpg_df = combined_cpg_df.iloc[:, [0,1,2,3,4,9,10]]
combined_cpg_df.columns = ["chrom","start","end","meth1","unmeth1","meth2","unmeth2"]

combined_cpg_df["meth_sum"] = combined_cpg_df["meth1"]+combined_cpg_df["meth2"]
combined_cpg_df["unmeth_sum"] = combined_cpg_df["unmeth1"]+combined_cpg_df["unmeth2"]
combined_cpg_df["strand"] = "."
combined_cpg_df["coverage"] = combined_cpg_df["meth_sum"] + combined_cpg_df["unmeth_sum"]

select_cols = ["chrom", "start", "end", "meth_sum", "unmeth_sum", "strand"]
final_cpg_df = combined_cpg_df.loc[:,select_cols]
final_cpg_df.to_csv(join(output_dir, "SL83768_SL83771_merged_CpG_context_deduplicated.bismark.cov"), header=False, index=False, sep="\t")

final_cpg_df_300x = combined_cpg_df[~(combined_cpg_df["coverage"] > 300)]
select_cols = ["chrom", "start", "end", "meth_sum", "unmeth_sum", "strand"]
final_cpg_df_300x = final_cpg_df_300x.loc[:,select_cols]
final_cpg_df_300x.to_csv(join(output_dir, "SL83768_SL83771_300x_merged_CpG_context_deduplicated.bismark.cov"), header=False, index=False, sep="\t")

# """ Combine the edited CpG files """
# edited_cpg_file_list =  glob.glob("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL*edited.cov")

# concat_df_list = []
# for each_file in edited_cpg_file_list:
# 	cpg_df =  pd.read_csv(each_file, sep="\t", header=None)
# 	print cpg_df
# 	cpg_df.columns = ["chrom", "start", "end", "meth", "unmeth", "strand"]
# 	concat_df_list.append(cpg_df)

# combined_cpg_df = pd.concat(concat_df_list,ignore_index=True)


# """ Perform chromosome wise groupby operation """
# #chrom_list = combined_cpg_df["chrom"].unique().tolist()
# chrom_list = [ "chr"+ str(each) for each in range(1,23) + ["X"]]

# chrom_concat_df_list = []
# for chrom in chrom_list:
# 	chrom_df = combined_cpg_df[combined_cpg_df["chrom"] == chrom]
# 	chrom_gp_df = chrom_df.groupby(["chrom", "start", "end"]).apply(lambda x : (x["meth"].sum(), x["unmeth"].sum()))
# 	chrom_gp_reset = chrom_gp_df.reset_index()
# 	chrom_gp_reset[["meth", "unmeth"]] = chrom_gp_reset.loc[:,0].apply(pd.Series)
# 	chrom_gp_reset = chrom_gp_reset.drop([0], axis=1)
# 	chrom_gp_reset["strand"] =  "."
# 	chrom_concat_df_list.append(chrom_gp_reset)

# combined_cpg_final_df = pd.concat(chrom_concat_df_list, ignore_index=True)

# with open(join(output_dir,"combined_cpg_final_df.pkl"), "w") as outfile:
# 	pickle.dump(combined_cpg_final_df, outfile)

# with open(join(output_dir,"combined_cpg_final_df.pkl")) as infile:
#     combined_cpg_final_df_1 = pickle.load(infile)

# combined_cpg_final_df.to_csv(join(output_dir, "SL83768_SL83771_CpG_context_deduplicated.bismark_edited.cov"), header=False, index=False, sep="\t")










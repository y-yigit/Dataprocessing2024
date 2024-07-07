import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Directory containing VarScan data files
data_directory = 'out/varscan'

# Initialize lists to store aggregated counts, depth coverage, and filenames
all_snps_counts = []
de_novo_snps_counts = []
depth_coverages = []
sample_names = []  # To store file names

# Iterate through files in the directory
for filename in os.listdir(data_directory):
    if filename.endswith('.snp'):  # Adjust the file extension as necessary
        file_path = os.path.join(data_directory, filename)
        
        # Store the filename without the extension
        sample_name_without_extension = filename.replace(".varScan.snp", "")
        sample_names.append(sample_name_without_extension)
        
        # Read the data from the file
        df = pd.read_csv(file_path, sep='\t')
        
        # Extract relevant columns and calculate counts
        df['Reads1'] = df['Cons:Cov:Reads1:Reads2:Freq:P-value'].str.split(':').str[2].astype(int)
        df['Reads2'] = df['Cons:Cov:Reads1:Reads2:Freq:P-value'].str.split(':').str[3].astype(int)
        
        # Aggregate counts
        total_snps = df['Reads1'].sum() + df['Reads2'].sum()
        de_novo_snps = df[df['Reads1'] == 0]['Reads2'].sum() + df[df['Reads2'] == 0]['Reads1'].sum()
        
        # Append to lists
        all_snps_counts.append(total_snps)
        de_novo_snps_counts.append(de_novo_snps)
        
        # Extract depth coverage for legend
        depth_coverage = df['Cons:Cov:Reads1:Reads2:Freq:P-value'].str.split(':').str[1]
        depth_coverages.append(float(depth_coverage.iloc[0]))  # Convert to float assuming it's a numeric depth value

# Plotting
plt.figure(figsize=(8, 6))  # Adjust figure size if necessary

# Normalize depth_coverages to use as color values
normalized_depths = (np.array(depth_coverages) - np.min(depth_coverages)) / (np.max(depth_coverages) - np.min(depth_coverages))

# Scatter plot with color mapped by depth_coverage
scatter = plt.scatter(all_snps_counts, de_novo_snps_counts, marker='o', s=100, c=normalized_depths, cmap='bwr_r')

# Annotate each point with its corresponding filename (without extension)
for i, txt in enumerate(sample_names):
    plt.annotate(txt, (all_snps_counts[i], de_novo_snps_counts[i]), fontsize=8)

plt.colorbar(scatter, label='Depth Of Coverage')
plt.xlabel('All SNPs')
plt.ylabel('De Novo Variants')
plt.title('QC Matrix')
plt.grid(True)
plt.tight_layout()
plt.savefig("out/images/QC_plot.png")
plt.show()


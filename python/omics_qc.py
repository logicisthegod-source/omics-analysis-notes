"""Quick QC pipeline for an omics feature abundance table."""
import numpy as np
import pandas as pd

# Rows = samples, columns = features
df = pd.read_csv("data/metabolites.csv", index_col=0)

# Overview of missing values
missing = df.isna().sum().sort_values(ascending=False)
print("Features with most missing values:\n", missing.head(10))

# Drop features missing in more than 20% of samples
df = df.loc[:, df.isna().mean() <= 0.20]

# Impute remaining NAs with half of each feature minimum
df = df.fillna(df.min() / 2)

# Log2 transform then z-score scaling
log_df = np.log2(df)
z_scaled = (log_df - log_df.mean()) / log_df.std()

z_scaled.to_csv("results/metabolites_qc_scaled.csv")
print("QC table shape:", z_scaled.shape)

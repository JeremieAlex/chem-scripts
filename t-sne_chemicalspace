import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# Load your data with SMILES in a column named 'SMILES'
data = pd.read_csv('12hr_screens_Feb2024/Final_BIO2_Legend_Feb2024_Reads_12hrs.csv', encoding='latin1').dropna(subset=['SMILES'])  # Ignore empty rows

# Add a column for RDKit Mol objects
data['Mol'] = data['SMILES'].apply(Chem.MolFromSmiles)
data = data.dropna(subset=['Mol'])

# Calculate Morgan fingerprints
data['FP'] = data['Mol'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024))

# Convert fingerprints to a list
fingerprints = [list(fp.ToBitString()) for fp in data['FP']]

# Apply standard scaling
scaler = StandardScaler()
fingerprints_scaled = scaler.fit_transform(fingerprints)

# Apply t-SNE for dimensionality reduction
tsne = TSNE(n_components=2, random_state=42)
tsne_result = tsne.fit_transform(fingerprints_scaled)

# Define the cutoff value for the column you want to highlight

# Convert 'mean_norm' column to numeric, coercing non-numeric values to NaN
data['mean_norm'] = pd.to_numeric(data['mean_norm'], errors='coerce')

# Drop rows with NaN values in 'mean_norm' column

mean = data['mean_norm'].mean()
std_dev = data['mean_norm'].std()

# Create a boolean mask for points satisfying the condition
condition = data['mean_norm'] < (mean - (std_dev))

# Create a scatter plot with red points for hits and blue points otherwise
plt.scatter(tsne_result[condition, 0], tsne_result[condition, 1], color='red', alpha=0.4, label='hit', s=0.5)
plt.scatter(tsne_result[~condition, 0], tsne_result[~condition, 1], color='blue', alpha=0.4, label='non-hit', s=0.5)

# Add legend
plt.legend(title='', loc='upper right', fontsize='large')

plt.title('hit vs non-hit space')
plt.savefig('plot.png', dpi=800)

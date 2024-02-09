#Rank order

import pandas as pd
import matplotlib.pyplot as plt

# Load your CSV file
csv_file_path = '12hr_screens_Feb2024/Copy of Copy of BIO2_Legend_Feb2024_Reads_12hrs.csv'
df = pd.read_csv(csv_file_path)

# Convert 'normavg' column to numeric, dropping non-numeric values
df['normavg'] = pd.to_numeric(df['rep_mean'])

# Drop rows with non-numeric 'normavg' values
df = df.dropna(subset=['normavg'])

# Sort the DataFrame by 'normavg' column in descending order
df_sorted = df.sort_values(by='normavg', ascending=False)

# Plot the 'normavg' column
plt.plot(df_sorted['normavg'].values)
plt.title(f'{csv_file_path} - Rank Order Plot')
plt.xlabel('Index')
plt.ylabel('normavg')
plt.show()

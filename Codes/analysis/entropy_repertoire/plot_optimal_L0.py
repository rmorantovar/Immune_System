import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file into a DataFrame
df = pd.read_csv('../../out/repertoire_entropy/data_processed.csv')

# Plot the dependence of parameter p on L0 and l
plt.figure(figsize=(10, 6))
for l, group_df in df.groupby('l'):
    plt.plot(group_df['L0'], group_df['p'], label=f'l={l}')
plt.xlabel('L0')
plt.ylabel('p')
plt.title('Dependence of p on L0 for different values of l')
plt.legend()
plt.xscale('log')
plt.yscale('linear')
# plt.grid(True)
plt.show()

import pandas as pd

# Input and output file paths
input_xlsx = "input_file.xlsx"  # Change this to your actual file path
output_csv = "output_file.csv"  # The resulting CSV file

# Specify the sheet name and columns to extract (modify as needed)
sheet_name = "Photoactivation CGG"  # Adjust if different
columns_to_keep = ["Sequence ID", "CDR3:"]  # Modify column names if necessary

# Load the Excel file
df = pd.read_excel(input_xlsx, sheet_name=sheet_name, header=0, skiprows=1)

# Check if required columns exist
missing_columns = [col for col in columns_to_keep if col not in df.columns]
if missing_columns:
    print(f"❌ Missing columns in Excel file: {missing_columns}")
    exit(1)

# Select only the required columns
df_selected = df[columns_to_keep].dropna()  # Remove rows with missing values

# Remove duplicate sequences (based on the "CDR3:" column)
df_unique = df_selected.drop_duplicates(subset=["CDR3:"])

# Save as CSV
df_unique.to_csv(output_csv, index=False)

print(f"✅ Process completed! CSV file saved as '{output_csv}'")

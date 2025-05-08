import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

alphabet = list('ARNDCQEGHILKMFPSTWYV')
def to_float(s):
    try:
        return float(s)
    except:
        return 0

def read_triangular_matrices(filename):
    matrices = []
    matrices_names = []
    current_matrix = []
    l = 20
    elements_read = 0
    matrix_name = ""
    row = 0



    with open(filename, 'r') as f:
        matrix_begin = False
        triangular = True
        for line in f:
            if line.startswith("H"):
                matrix_name = line.strip().split()[1]
                continue
            if line.startswith("M"):
                matrix_begin = True
                continue
            if matrix_begin==True:
                elements = list(map(to_float, line.strip().split()))
                elements_read += len(elements)
                if row == 0:
                    if elements_read != 1:
                        triangular = False

                current_matrix.extend(elements)
                
                if triangular == True:
                    if elements_read == ((l * (l + 1)) // 2):
                        square_matrix = np.zeros((l, l))
                        idx = 0
                        for i in range(l):
                            for j in range(0, i+1):
                                square_matrix[i, j] = current_matrix[idx]
                                square_matrix[j, i] = current_matrix[idx]
                                idx += 1
                        df = pd.DataFrame(square_matrix)
                        df.to_csv("../Input_files/matrices/" + f"{matrix_name}.txt", sep="\t", index=False, header=False, float_format="%.6f")
                        print(f"Matrix '{matrix_name}' saved as {matrix_name}.txt")
                        matrices.append(df)
                        matrices_names.append(matrix_name)
                        current_matrix = []
                        #l = 1
                        elements_read = 0
                        matrix_begin = False
                        row = 0
                        triangular = True
                    else:
                        row += 1
                else:
                    if elements_read == ((l * l)):
                        square_matrix = np.zeros((l, l))
                        idx = 0
                        for i in range(l):
                            for j in range(l):
                                square_matrix[i, j] = current_matrix[idx]
                                #square_matrix[j, i] = current_matrix[idx]
                                idx += 1
                        df = pd.DataFrame(square_matrix)
                        df.to_csv("../Input_files/matrices/" + f"{matrix_name}.txt", sep="\t", index=False, header=False, float_format="%.6f")
                        print(f"Matrix '{matrix_name}' saved as {matrix_name}.txt")
                        matrices.append(df)
                        matrices_names.append(matrix_name)
                        current_matrix = []
                        #l = 1
                        elements_read = 0
                        matrix_begin = False
                        row = 0
                        triangular = True
                    else:
                        row += 1

    return matrices, matrices_names

def plot_heatmap(matrix, matrix_name):
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.set(font_scale=1)
    sns.heatmap(matrix, annot=True, annot_kws = {'fontsize':5}, fmt=".2f", cmap="coolwarm", linewidths=.5, square=True, cbar=True, ax = ax, xticklabels = alphabet, yticklabels = alphabet)
    ax.set_title(f"{matrix_name} Heatmap")
    fig.savefig('../../Figures/matrices/' + f"{matrix_name}.pdf")
    plt.close(fig)

filename = "../Input_files/aaindex3.txt"
matrices, matrices_names = read_triangular_matrices(filename)

for i, matrix in enumerate(matrices):
    plot_heatmap(matrix, matrices_names[i])








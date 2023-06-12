import pandas as pd
from sklearn.decomposition import PCA


def sc_pca(input_path, output_path, n_comp=5):
    input_data = pd.read_table(input_path, index_col=0)
    pca = PCA(n_components=n_comp)
    pca.fit(input_data)
    pcs = []
    for i in range(1, n_comp + 1):
        x = "PC{0}".format(i)
        pcs.append(x)
    pc_df = pd.DataFrame(pca.components_, columns=input_data.columns, index=pcs)
    pc_df.to_csv(output_path, sep="\t")
    return (pc_df)


import pandas as pd


def is_meta_valid(exp_path, meta_path):

    try:
        cells_exp = pd.read_table(exp_path, index_col=0).columns

        cells_meta = pd.read_table(meta_path, index_col=0).index

        if cells_meta != cells_exp:
            return False, "Cells in the meta data does not match with the cells in the expression data."

    except Exception as e:
        return False, e.message


def is_matrix_valid(exp_path):

    try:
        exp_mat = pd.read_table(exp_path, index_col=0)

        if (exp_mat < 0).values.any():
            return False, "Expression matrix contains negative values."
        if exp_mat.isnull().values.any():
            return False, "Expression matrix contains missing values."
        if not exp_mat.apply(lambda s: pd.to_numeric(s, errors='coerce').notnull().all()).all():
            return False, "Expression matrix is not numeric."
        if not exp_mat.index.is_unique:
            return False, "Genes are not unique."
        if not exp_mat.columns.is_unique:
            return False, "Cell names are not unique."
        else:
            return True, "Data can be processed."

    except Exception as e:
        return False, e.message



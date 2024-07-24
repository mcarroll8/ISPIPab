# Evan Edelstein
import pandas as pd
from sklearn.metrics import f1_score, matthews_corrcoef, precision_recall_curve, roc_curve, auc
from .compare_auc_delong_xu import delong_roc_test


def trainningsetstats():
    # df = pd.read_csv("/Users/evanedelstein/Desktop/MetaDPIv2/metadpi/input/input_data_all.csv")
    test_set = pd.read_csv(
        "/Users/evanedelstein/Desktop/MetaDPIv2/metadpi/input/test_all_prot.csv")
    proteins = test_set["protein"].tolist()
    proteins = [i.upper() for i in proteins]
    # df = pd.read_csv("/Users/evanedelstein/Desktop/Research_Evan/Raji_Summer2019_atom/Meta_DPI/Data/Test_data/meta-ppisp-results-comma-new.txt")
    df = pd.read_csv(
        "/Users/evanedelstein/Desktop/Research_Evan/Raji_Summer2019_atom/Meta_DPI/Data/vorffip_renumbered.txt")

    df["protein"] = [x.split('_')[1] for x in df['residue']]
    df.isnull().any()
    df = df.fillna(method='ffill')
    df = df[df['protein'].isin(proteins)]
    cutoff_dict = cutoff_file_parser(
        "/Users/evanedelstein/Desktop/MetaDPIv2/metadpi/input/cutoffs.csv")
    output = pd.DataFrame()
    for i in ['vorffip']:
        top = df.sort_values(by=[i], ascending=False).groupby((["protein"])).apply(
            lambda x: x.head(cutoff_dict[x.name])).index.get_level_values(1).tolist()

        df[f'{i}_bin'] = [
            1 if i in top else 0 for i in df.index.tolist()]

        output[f"{i}_fscore"] = df.groupby((["protein"])).apply(
            lambda x: fscore(x, 'annotated', i))
        output[f"{i}_mcc"] = df.groupby((["protein"])).apply(
            lambda x: mcc(x, 'annotated', i))

    # roc_and_pr_dic = roc_and_pr(df, "annotated", "meta-ppisp")
    roc_and_pr_dic = roc_and_pr(df, "annotated", "vorffip")
    f = fscore(df, "annotated", "vorffip")
    m = mcc(df, "annotated", "vorffip")
    print("fscore: ", f, "mcc: ", m)
    roc_df = pd.DataFrame()
    prdf = pd.DataFrame()
    roc_df["fpr"] = roc_and_pr_dic["fpr"]
    roc_df["tpr"] = roc_and_pr_dic["tpr"]

    prdf["precision"] = roc_and_pr_dic["precision"]
    prdf["recall"] = roc_and_pr_dic["recall"]
    roc_df.to_csv(
        "/Users/evanedelstein/Desktop/MetaDPIv2/metadpi/output/vorffip_roc.csv")
    prdf.to_csv(
        "/Users/evanedelstein/Desktop/MetaDPIv2/metadpi/output/vorffip_pr.csv")

    output.to_csv(
        "/Users/evanedelstein/Desktop/MetaDPIv2/metadpi/output/vorffip.csv")
    return


def fscore(x, annotated_col, pred) -> tuple:
    return f1_score(x[annotated_col], x[f'{pred}_bin'])


def mcc(x, annotated_col, pred) -> tuple:
    return matthews_corrcoef(x[annotated_col], x[f'{pred}_bin'])


def cutoff_file_parser(cutoff_frame) -> dict:
    cutoff_frame_df: pd.DataFrame = pd.read_csv(cutoff_frame)
    cutoff_frame_df.columns = cutoff_frame_df.columns.str.lower()
    cutoff_dict = dict(
        zip(cutoff_frame_df['protein'], cutoff_frame_df['cutoff res']))
    return cutoff_dict


def roc_and_pr(test_frame: pd.DataFrame, annotated_col, pred) -> dict:
    fpr, tpr, roc_thresholds = roc_curve(
        test_frame[annotated_col], test_frame[pred])
    roc_auc = round(auc(fpr, tpr), 3)
    precision, recall, pr_thresholds = precision_recall_curve(
        test_frame[annotated_col], test_frame[pred])
    pr_auc = round(auc(recall, precision), 3)
    roc_and_pr_dic: dict = {"fpr": fpr, "tpr": tpr, "roc_thresholds": roc_thresholds, "roc_auc": roc_auc,
                            "precision": precision, "recall": recall, "pr_thresholds": pr_thresholds, "pr_auc": pr_auc}
    print("Roc: ", roc_auc, "pr: ", pr_auc)
    return roc_and_pr_dic


trainningsetstats()

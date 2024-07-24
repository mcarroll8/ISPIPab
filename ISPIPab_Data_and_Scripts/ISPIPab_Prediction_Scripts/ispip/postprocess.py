# Evan Edelstein
from sklearn.metrics import auc, matthews_corrcoef, f1_score, precision_recall_curve, roc_curve
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from .compare_auc_delong_xu import delong_roc_test
from functools import reduce


def fscore_mcc(x, annotated_col, pred) -> tuple:
    return f1_score(x[annotated_col], x[f'{pred}_bin']), matthews_corrcoef(x[annotated_col], x[f'{pred}_bin'])


def compute_fscore(x, annotated_col, pred) -> tuple:
    return f1_score(x[annotated_col], x[f'{pred}_bin'])


def compute_mcc(x, annotated_col, pred) -> tuple:
    return matthews_corrcoef(x[annotated_col], x[f'{pred}_bin'])


def postprocess(test_frame, predicted_col,cutoff_frame, use_cutoff_from_file, annotated_col, autocutoff) -> tuple:
    proteins = test_frame.protein.unique()
    results = []
    roc_curve_data: list = []
    pr_curve_data: list = []
    fscore_mcc_dfs = []
    cutoff_dict = cutoff_file_parser(cutoff_frame) if (
        use_cutoff_from_file) else {protein: autocutoff for protein in proteins}
    # make testframe and cutoff dict this self in class
    params = [(pred, cutoff_dict, test_frame, annotated_col)
              for pred in predicted_col]

    with ProcessPoolExecutor(max_workers=4) as exe:
        return_vals = exe.map(analyses, params)
        for return_val in return_vals:
            test_frame[f'{return_val[0]}_bin'] = return_val[1]
            fscore_mcc_dfs.append(return_val[2])
            results.append(return_val[3])
            roc_curve_data.append(return_val[4])
            pr_curve_data.append(return_val[5])

    fscore_mcc_by_protein = reduce(lambda left, right: pd.merge(
        left, right,  on="protein", how='outer'), fscore_mcc_dfs)
    result_df = pd.DataFrame(
        results, columns=['predictor', 'f-score', 'mcc', 'roc_auc', 'pr_auc'])
    stats_df = pd.DataFrame(index=predicted_col, columns=predicted_col)
    test_frame = test_frame.sort_values(by=annotated_col, ascending=False)
    for index in predicted_col:
        for column in predicted_col:
            if index == column:
                stats_df.loc[index, column] = index
            else:
                pval, test, auc_diff = statistics(
                    test_frame, annotated_col, index, column)
                stats_df.loc[index, column] = pval  # below diagnol
                stats_df.loc[column, index] = auc_diff  # above diagnol

    return result_df, roc_curve_data, pr_curve_data, test_frame, fscore_mcc_by_protein, stats_df


def cutoff_file_parser(cutoff_frame) -> dict:
    cutoff_frame_df: pd.DataFrame = pd.read_csv(cutoff_frame)
    cutoff_frame_df.columns = cutoff_frame_df.columns.str.lower()
    cutoff_dict = dict(
        zip(cutoff_frame_df['protein'], cutoff_frame_df['cutoff res']))
    return cutoff_dict


def analyses(params) -> tuple:
    pred, cutoff_dict, test_frame, annotated_col = params
    top = test_frame.sort_values(by=[pred], ascending=False).groupby((["protein"])).apply(
        lambda x: x.head(cutoff_dict[x.name])).index.get_level_values(1).tolist()
    test_frame[f'{pred}_bin'] = [
        1 if i in top else 0 for i in test_frame.index.tolist()]

    fscores = test_frame.groupby(["protein"]).apply(
        lambda x: compute_fscore(x, annotated_col, pred))
    mccs = test_frame.groupby(["protein"]).apply(
        lambda x: compute_mcc(x, annotated_col, pred))

    fscores = pd.DataFrame({f"{pred}_fscore": fscores}).reset_index()
    mccs = pd.DataFrame({f"{pred}_mcc": mccs}).reset_index()
    fscore_mcc_by_protein = pd.merge(fscores, mccs, on="protein")

    fscore, mcc = fscore_mcc(test_frame, annotated_col, pred)
    roc_and_pr_dic = roc_and_pr(test_frame, annotated_col, pred)

    results_list = [pred, fscore, mcc,
                    roc_and_pr_dic["roc_auc"], roc_and_pr_dic["pr_auc"]]
    roclist = [pred, roc_and_pr_dic["fpr"], roc_and_pr_dic["tpr"],
               roc_and_pr_dic["roc_auc"], roc_and_pr_dic["roc_thresholds"]]

    prlist = [pred, roc_and_pr_dic["recall"], roc_and_pr_dic["precision"],
              roc_and_pr_dic["pr_auc"], roc_and_pr_dic["pr_thresholds"]]

    return pred, test_frame[f'{pred}_bin'], fscore_mcc_by_protein, results_list, roclist, prlist


def statistics(x, annotated_col, pred1, pred2) -> tuple:
    y_true = x[annotated_col]
    y1 = x[pred1]
    y2 = x[pred2]
    log10_pval, aucs = delong_roc_test(y_true, y1, y2)
    aucs = aucs.tolist()
    dauc = round(aucs[1] - aucs[0], 3)
    log10_pval = round(log10_pval.tolist()[0][0], 3)
    test = "signifigant" if log10_pval < -1.3 else "not significant"
    return log10_pval, test, dauc

# name better, type correct


def roc_and_pr(test_frame: pd.DataFrame, annotated_col, pred) -> dict:
    fpr, tpr, roc_thresholds = roc_curve(
        test_frame[annotated_col], test_frame[pred])
    roc_auc = round(auc(fpr, tpr), 3)
    precision, recall, pr_thresholds = precision_recall_curve(
        test_frame[annotated_col], test_frame[pred])
    pr_auc = round(auc(recall, precision), 3)
    roc_and_pr_dic: dict = {"fpr": fpr, "tpr": tpr, "roc_thresholds": roc_thresholds, "roc_auc": roc_auc,
                            "precision": precision, "recall": recall, "pr_thresholds": pr_thresholds, "pr_auc": pr_auc}
    return roc_and_pr_dic

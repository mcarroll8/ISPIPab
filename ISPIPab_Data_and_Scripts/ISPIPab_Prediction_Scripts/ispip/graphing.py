# Evan Edelstein
import subprocess
import matplotlib.pyplot as plt
import dtreeviz.trees
import pandas as pd
import pathlib
import os


def roc_viz(roc_curve_data, output_path_dir, model_name):
    roc_frames = []
    plt.figure()
    lw = 2
    for data in roc_curve_data:
        pred, fpr, tpr, roc_auc, thresholds = data
        new_roc_frame = pd.DataFrame()
        new_roc_frame[f"{pred}_fpr"] = fpr
        new_roc_frame[f"{pred}_tpr"] = tpr
        roc_frames.append(new_roc_frame)
        plt.plot(fpr, tpr, lw=lw, label=f"{pred} (area = {roc_auc})")

    roc_frame = pd.concat(roc_frames, axis=1)
    roc_frame.to_csv(f"{output_path_dir}/roc_{model_name}.csv")
    plt.plot([0, 1], [0, 1], color="gray", lw=lw, linestyle="--", alpha=0.5)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Receiver operating characteristic Curve")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(f"{output_path_dir}/ROC_{model_name}.png")
    # plt.show()
    return


def pr_viz(pr_curve_data, output_path_dir, model_name, df, annotated_col):
    pr_frames = [pd.DataFrame()]
    plt.figure()
    lw = 2
    for data in pr_curve_data:
        pred, recall, precision, pr_auc, thresholds = data
        pr_frame = pd.DataFrame()
        pr_frame[f"{pred}_fpr"] = recall
        pr_frame[f"{pred}_tpr"] = precision
        pr_frames.append(pr_frame)
        plt.plot(recall, precision, lw=lw, label=f"{pred} (area = {pr_auc})")

    pr_frame = pd.concat(pr_frames, axis=1)
    pr_frame.to_csv(f"{output_path_dir}/pr_{model_name}.csv")
    no_skill = len(df[df[annotated_col] == 1]) / len(df)
    plt.plot([0, 1], [no_skill, no_skill], linestyle='--', color='gray')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Precision-Recall curve")
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(f"{output_path_dir}/PR_{model_name}.png")
    # plt.show()
    return


def treeviz(tree, df, feature_cols, annotated_col, model_name, output_path_dir):
    try:
        viz = dtreeviz(tree,
                       df[feature_cols],
                       df[annotated_col],
                       target_name='Interface',
                       feature_names=feature_cols,
                       class_names=["non_interface", "interface"],
                       show_node_labels=True,
                       fancy=False
                       )
        path = f"{output_path_dir}/Rftree_{model_name}.svg"
        viz.save(path)
    except Exception as exception:
        print(exception)
        print("Please install Graphviz")
    return


def roc_to_csv(roc_curve_data, output_path_dir, model_name):
    roc_frame = pd.DataFrame()
    for data in roc_curve_data:
        pred, fpr, tpr, roc_auc, thresholds = data
        roc_frame[f"{pred}_fpr"] = fpr
        roc_frame[f"{pred}_tpr"] = tpr
    roc_frame.to_csv(f"{output_path_dir}/roc_{model_name}.csv")
    return


def pr_to_csv(pr_curve_data, output_path_dir, model_name):
    pr_frame = pd.DataFrame()
    for data in pr_curve_data:
        pred, recall, precision, pr_auc, thresholds = data
        pr_frame[f"{pred}_fpr"] = recall
        pr_frame[f"{pred}_tpr"] = precision
    pr_frame.to_csv(f"{output_path_dir}/pr_{model_name}.csv")
    return


def pymol_viz(bin_frame, proteins, predicted_col, annotated_col, pymolscriptpath, output_path_dir):
    new_folder = os.path.join(output_path_dir, "proteins")

    pathlib.Path(new_folder).mkdir(parents=True, exist_ok=True)
    for protein in proteins:
        protein_df = bin_frame[bin_frame['protein'] == protein]
        new_folder = os.path.join(output_path_dir, "proteins", protein)
        pathlib.Path(new_folder).mkdir(parents=True, exist_ok=True)

        for pred in predicted_col:

            pred_residues = protein_df[protein_df[f'{pred}_bin'] == 1].index.tolist()
            annotated_resiues = protein_df[protein_df[annotated_col] == 1].index.tolist()
            tp_residues = [i for i in pred_residues if i in annotated_resiues]
            fp_residues = [i for i in pred_residues if i not in annotated_resiues]
            tp_residues = [i.split("_")[0] for i in tp_residues]
            tp_residues = "+".join(tp_residues)
            fp_residues = [i.split("_")[0] for i in fp_residues]
            fp_residues = "+".join(fp_residues)
            annotated_resiues = [i.split("_")[0] for i in annotated_resiues]
            annotated_resiues = "+".join(annotated_resiues)

            subprocess.run(
                f"pymol -Q -c -q  {pymolscriptpath} -- {output_path_dir} {protein} {tp_residues} {fp_residues} {annotated_resiues} {pred}", shell=True)

    return

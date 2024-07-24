import os
from .argscontainer import ArgsContainer
from .cli_interface import userinterface
import pandas as pd
from .postprocess import postprocess
from .crossvalidation import hyperparamtertuning_and_crossvalidation
from .generate import generate
from .predict import predict
from .graphing import roc_viz, pr_viz, treeviz, pymol_viz
from .preprocess import data_preprocesss, data_split_auto,\
    data_split_from_file, cross_validation_set_generater
pd.options.mode.chained_assignment = None  # default='warn'


def main() -> None:
    # get Command line arguments and defaults.
    args_container: ArgsContainer = userinterface()

    # read input file containing residues, individual predictors and annotated columns.
    df: pd.DataFrame = pd.read_csv(args_container.input_frames_file)

    # preprocess data -> remove any null or missing data from the dataset and check that annoted is numbernulls.
    df, feature_cols, annotated_col, proteins = data_preprocesss(df)

    predicted_col = feature_cols + \
        ['logisticregresion', "linearregression",
            'randomforest'] + args_container.models_to_use

    # Mode 1: predict
    if args_container.mode == 'predict':
        df = predict(df, feature_cols, args_container.input_folder_path,
                     args_container.model_name, args_container.nn, args_container.xg)

        results_df, roc_curve_data, pr_curve_data, bin_frame, fscore_mcc_by_protein, stats_df = postprocess(
            df, predicted_col, args_container.cutoff_frame,args_container.use_cutoff_from_file, annotated_col, args_container.autocutoff)
        df_saver(results_df, "results", args_container.output_path_dir)
        df_saver(bin_frame, "bin_frame", args_container.output_path_dir)
        df_saver(fscore_mcc_by_protein, "fscore_mcc_by_protein",
                 args_container.output_path_dir)
        visualization(roc_curve_data, pr_curve_data, None, df, feature_cols,
                      annotated_col, predicted_col, df, bin_frame, args_container)

    # Mode 2: generate learned model

    elif args_container.mode == 'generate':
        models, tree = generate(df, feature_cols, annotated_col, args_container.output_path_dir,
                                args_container.model_name, args_container.rf_params, args_container.nn, args_container.xg)

        if (tree is not None) and args_container.save_tree:
            treeviz(tree, df, feature_cols, annotated_col,
                    args_container.model_name, args_container.output_path_dir)

    # Mode 3: test/train

    elif args_container.mode == 'test':

        if args_container.use_test_train_files:
            test_frame, train_frame = data_split_from_file(df, args_container)
        else:
            test_frame, train_frame = data_split_auto(df, proteins)

        print(f'length of test set: {len(test_frame)}',
              f"length of training set: {len(train_frame)}")

        # train
        models, tree = generate(train_frame, feature_cols, annotated_col, args_container.output_path_dir,
                                args_container.model_name, args_container.rf_params, args_container.nn, args_container.xg)  # train
    
        # test
        test_frame = predict(test_frame, feature_cols, args_container.input_folder_path,
                             args_container.model_name, args_container.nn, args_container.xg, models)

        results_df, roc_curve_data, pr_curve_data, bin_frame, fscore_mcc_by_protein, stats_df = postprocess(
            test_frame, predicted_col, args_container.cutoff_frame,args_container.use_cutoff_from_file, annotated_col, args_container.autocutoff)

        df_saver(results_df, "results", args_container.output_path_dir)
        df_saver(bin_frame, "bin_frame", args_container.output_path_dir)
        df_saver(fscore_mcc_by_protein, "fscore_mcc_by_protein",
                 args_container.output_path_dir)
        df_saver(stats_df, "pairtest", args_container.output_path_dir)

        visualization(roc_curve_data, pr_curve_data, tree, df, feature_cols,
                      annotated_col, predicted_col, test_frame, bin_frame, args_container)
        print(results_df)

    # Mode 4: Cross-validation:

    elif args_container.mode == 'cv':
        test_frame, cvs, train_proteins = cross_validation_set_generater(
            args_container.cvs_path, df)
        # print(test_frame, cvs, train_proteins)
        models = hyperparamtertuning_and_crossvalidation(
            df, train_proteins, feature_cols, annotated_col, args_container)

        model_param_writer(models, args_container.output_path_dir)
        test_frame = predict(test_frame, feature_cols, args_container.input_folder_path,
                             args_container.model_name, args_container.nn, args_container.xg, models)
        results_df, roc_curve_data, pr_curve_data, bin_frame, fscore_mcc_by_protein, stats_df = postprocess(
            test_frame, predicted_col, args_container.cutoff_frame,args_container.use_cutoff_from_file, annotated_col, args_container.autocutoff)

        df_saver(results_df, "results", args_container.output_path_dir)
        df_saver(bin_frame, "bin_frame", args_container.output_path_dir)
        df_saver(fscore_mcc_by_protein, "fscore_mcc_by_protein",
                 args_container.output_path_dir)
        df_saver(stats_df, "pairtest", args_container.output_path_dir)

        visualization(roc_curve_data, pr_curve_data, None, df, feature_cols,
                      annotated_col, predicted_col, test_frame, bin_frame, args_container)
        print(results_df)

    # Mode 5: visualization

    elif args_container.mode == 'viz':
        test_frame, cvs, train_proteins = cross_validation_set_generater(args_container.cvs_path, df)
        # bin_frame = pd.read_csv(os.path.join(
        #     args_container.output_path_dir, "bin_frame.csv"), index_col=0)
        bin_frame = pd.read_csv(args_container.results_df_input, index_col=0)
        protein_to_viz = test_frame["protein"].unique()

        pymol_viz(bin_frame, protein_to_viz, predicted_col, annotated_col,
                  args_container.pymolscriptpath, args_container.output_path_dir)
    # Mode 6: Reprocess: 
    elif args_container.mode == 'reprocess':
        results_df = pd.read_csv(args_container.results_df_input)
        results_df, roc_curve_data, pr_curve_data, bin_frame, fscore_mcc_by_protein, stats_df = postprocess(
            results_df, predicted_col, args_container.cutoff_frame,args_container.use_cutoff_from_file, annotated_col, args_container.autocutoff)
        
        df_saver(results_df, "results", args_container.output_path_dir)
        df_saver(bin_frame, "bin_frame", args_container.output_path_dir)
        df_saver(fscore_mcc_by_protein, "fscore_mcc_by_protein",
                 args_container.output_path_dir)
        visualization(roc_curve_data, pr_curve_data, None, df, feature_cols,
                      annotated_col, predicted_col, df, bin_frame, args_container)

    else:
        print("mode is set incorrectly")
    return


def df_saver(df, name, output_path_dir):
    out = os.path.join(output_path_dir, f'{name}.csv')
    df.to_csv(out)
    return


def model_param_writer(models, output_path_dir):
    rf_model, linear_model, logit_model, nn_model, xgb_model = models

    ensable_param: list = [f"{type(model).__name__}: {model.get_params()}\n"
                           for model in filter(None, [rf_model, nn_model, xgb_model])]

    regr_param: list = [f"{type(model).__name__}:    {model.coef_}\n"
                        for model in filter(None, [linear_model, logit_model])]

    out = os.path.join(output_path_dir, 'best_parameters.txt')
    with open(out, 'w+') as file:
        file.writelines(regr_param)
        file.writelines(ensable_param)
    return


def visualization(roc_curve_data, pr_curve_data, tree, df, feature_cols, annotated_col, predicted_col, test_frame, bin_frame, args_container):
    roc_viz(roc_curve_data, args_container.output_path_dir,
            args_container.model_name)
    pr_viz(pr_curve_data, args_container.output_path_dir,
           args_container.model_name, test_frame, annotated_col)
    if (tree is not None) and args_container.save_tree:
        treeviz(tree, df, feature_cols, annotated_col,
                args_container.model_name, args_container.output_path_dir)
    protein_to_viz = test_frame["protein"].unique()
    if args_container.usepymol:
        pymol_viz(bin_frame, protein_to_viz, predicted_col, annotated_col,
                  args_container.pymolscriptpath, args_container.output_path_dir)
    return

# if __name__ == "__main__":
#     main()


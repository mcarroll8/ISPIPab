# Evan Edelstein

import joblib
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingRegressor
# from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import GridSearchCV, cross_validate
from itertools import chain
import pandas as pd
import os

"""
currently scorer is roc_auc try adding average_precision as well
"""


def hyperparamtertuning_and_crossvalidation(df: pd.DataFrame, cvs, feature_cols, annotated_col, args_container) -> list:
    df.reset_index(level=0, inplace=True)
    cviterator = []
    for c, testk in enumerate(cvs):
        traink = cvs[:c] + cvs[c + 1:]  # train is k-1, test is k
        traink = list(chain.from_iterable(traink))
        trainindices = df[df["protein"].isin(traink)].index.values.astype(int)
        testindices = df[df["protein"].isin(testk)].index.values.astype(int)
        cviterator.append((trainindices, testindices))

    p_grid = {"n_estimators": [10,50,100, 200], "max_depth": [
        None, 5, 10, 15], "ccp_alpha": [0.0, 0.25, 0.5, 0.75], "bootstrap": [True, False]}

    rf_model = GridSearchCV(estimator=RandomForestClassifier(
    ), param_grid=p_grid, cv=cviterator, scoring="roc_auc").fit(df[feature_cols], df[annotated_col])
    hyperparam(rf_model, cvs, args_container.output_path_dir, "rf")
    rf_model = rf_model.best_estimator_

    logit_frame = pd.DataFrame(cross_validate(LogisticRegression(
    ), df[feature_cols], df[annotated_col], cv=cviterator, return_estimator=True, scoring="roc_auc"))
    print(logit_frame)
    logit_model = logit_frame.loc[
        logit_frame['test_score'].idxmax(), "estimator"]

    crossval_chart(logit_frame, args_container.output_path_dir, "logit")

    linmodel_frame = pd.DataFrame(cross_validate(LinearRegression(
    ), df[feature_cols], df[annotated_col], cv=cviterator, return_estimator=True, scoring="roc_auc"))
    print(linmodel_frame)
    linear_model = linmodel_frame.loc[linmodel_frame['test_score'].idxmax(
    ), "estimator"]
    crossval_chart(linmodel_frame, args_container.output_path_dir, "linear")

    # IN DEVELOPMENT 
    # if args_container.nn:
    #     param_grid = {'hidden_layer_sizes': [(50, 50, 50), (50, 100, 50), (100, 1)],
    #                   'alpha': [0.0001, 0.05],
    #                   'learning_rate': ['constant', 'adaptive'],
    #                   }
    #     nn_model = GridSearchCV(estimator=MLPRegressor(), param_grid=param_grid,
    #                             cv=cviterator, scoring="roc_auc").fit(df[feature_cols], df[annotated_col])
    #     hyperparam(nn_model, cvs, args_container.output_path_dir, "NN")
    #     nn_model = nn_model.best_estimator_
    # else:
    #     nn_model = None
    nn_model = None

    if args_container.xg:
        xparam_grid = {"loss": ["squared_error", "absolute_error", "poisson"]}
        xgb_model = GridSearchCV(estimator=HistGradientBoostingRegressor(), param_grid=xparam_grid, cv=cviterator, scoring="roc_auc").fit(df[feature_cols], df[annotated_col])
        hyperparam(xgb_model, cvs, args_container.output_path_dir, "xgb")
        xgb_model.best_estimator_
    else:
        xgb_model = None

    models = [rf_model, linear_model, logit_model, nn_model, xgb_model]

    _ = [joblib.dump(i, f"{args_container.output_path_dir}/{type(i).__name__}_{args_container.model_name}.joblib",
                     compress=3) for i in filter(None, models)]
    return models


def hyperparam(model, cvs, output_path_dir, name):
    cv_results = pd.DataFrame(model.cv_results_)
    chart1_list = []
    chart2_list = []
    for c, cv in enumerate(cvs):
        col = f"split{c}_test_score"
        mean = round(cv_results[col].mean(), 3)
        stdev = round(cv_results[col].std(), 3)
        print("\n\n", cv_results[col].idxmax(),
              "\n\n", cv_results[col], "\n\n")
        best_score = cv_results.loc[cv_results[col].idxmax(), [
            'params', col, ]].values
        chart1_list.append(
            [c, best_score[0], round(best_score[1], 3), mean, stdev])

    means = model.cv_results_["mean_test_score"]
    stds = model.cv_results_["std_test_score"]
    for mean, std, params in zip(means, stds, model.cv_results_["params"]):
        chart2_list.append([params, mean, std * 2])

    by_params = pd.DataFrame(chart2_list, columns=[
                             "params", "mean roc auc", "+/- std"])
    best_score = pd.DataFrame(chart1_list, columns=[
                              "CV", "best params", "best ROC_AUC", "mean", "stdev"])
    out = os.path.join(output_path_dir, f'{name}_cv_results.csv')
    cv_results.to_csv(out)
    out = os.path.join(output_path_dir, f'{name}_by_parametrs.csv')
    by_params.to_csv(out)
    out = os.path.join(output_path_dir, f'{name}_best_score.csv')
    best_score.to_csv(out)
    return


def crossval_chart(df, output_path_dir, name):
    df["coefs"] = [i.coef_ for i in df["estimator"].values.tolist()]
    out = os.path.join(output_path_dir, f'{name}_best_score.csv')
    df.to_csv(out)
    return

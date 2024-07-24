# Evan Edelstein
import pandas as pd
from sklearn import linear_model
import joblib
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingRegressor
# from sklearn.neural_network import MLPRegressor


def generate(df, feature_cols, annotated_col, output_path_dir, model_name, rf_params, nn, xg):
    rf_model, tree = randomforest_generate_model(
        df, feature_cols, annotated_col, rf_params, output_path_dir, model_name)
    linreg_model = linreg_generate_model(
        df, feature_cols, annotated_col, output_path_dir, model_name)
    logreg_model = logreg_generate_model(
        df, feature_cols, annotated_col, output_path_dir, model_name)
    NN_model = None
    # NN_model = NueralNet_generate_model(
    #     df, feature_cols, annotated_col, output_path_dir, model_name) if nn else None
    xgboost_model = xgboost_generate_model(
        df, feature_cols, annotated_col, output_path_dir, model_name) if xg else None
    models = [rf_model, linreg_model, logreg_model, NN_model, xgboost_model]
    return models, tree


def randomforest_generate_model(df: pd.DataFrame, feature_cols, annotated_col, rf_params, output_path_dir, model_name):
    trees, depth, ccp = rf_params
    X = df[feature_cols]
    y = df[annotated_col]
    model = RandomForestClassifier(
        n_estimators=trees, random_state=0, bootstrap=False, max_depth=depth, ccp_alpha=ccp).fit(X, y)
    tree = model.estimators_[0]

    joblib.dump(model, f"{output_path_dir}/RF_{model_name}.joblib", compress=3)
    return model, tree


def linreg_generate_model(df, feature_cols, annotated_col, output_path_dir, model_name):
    regr = linear_model.LinearRegression().fit(
        df[feature_cols], df[annotated_col])

    joblib.dump(
        regr, f"{output_path_dir}/LinRegr_{model_name}.joblib", compress=3)
    return regr


def logreg_generate_model(df, feature_cols, annotated_col, output_path_dir, model_name):
    regr = linear_model.LogisticRegression().fit(
        df[feature_cols], df[annotated_col])

    joblib.dump(
        regr, f"{output_path_dir}/Logit_{model_name}.joblib", compress=3)
    return regr

# IN DEVELOPMENT 
# def NueralNet_generate_model(df, feature_cols, annotated_col, output_path_dir, model_name): 
#     regr = MLPRegressor().fit(df[feature_cols], df[annotated_col])
#     joblib.dump(regr, f"{output_path_dir}/NN_{model_name}.joblib", compress=3)
#     return regr


def xgboost_generate_model(df, feature_cols, annotated_col, output_path_dir, model_name):
    regr = HistGradientBoostingRegressor().fit(
        df[feature_cols], df[annotated_col])

    joblib.dump(regr, f"{output_path_dir}/XGB_{model_name}.joblib", compress=3)
    return regr

# Evan Edelstein
import pandas as pd
import joblib


def predict(df, feature_cols, input_folder_path, model_name, nn, xg, models=None) -> pd.DataFrame:
    if models is not None:
        rf_model, linreg_model, logreg_model, nn_model, xgboost_model = models
    else:
        rf_model = joblib.load(f"{input_folder_path}/RF_{model_name}.joblib")
        logreg_model = joblib.load(
            f"{input_folder_path}/Logit_{model_name}.joblib")
        linreg_model = joblib.load(
            f"{input_folder_path}/LinRegr_{model_name}.joblib")
        # nn_model = joblib.load(
        #     f"{input_folder_path}/NN_{model_name}.joblib") if nn else None
        xgboost_model = joblib.load(
            f"{input_folder_path}/XGB_{model_name}.joblib") if xg else None
    out_df: pd.DataFrame = randomforest_predict_from_trained_model(
        df, feature_cols, rf_model)
    out_df: pd.DataFrame = logreg_predict_from_trained_model(
        out_df, feature_cols, logreg_model)
    out_df: pd.DataFrame = linreg_predict_from_trained_model(
        out_df, feature_cols, linreg_model)
    # out_df: pd.DataFrame = nueralnet_predict_from_trained_model(
    #     out_df, feature_cols, nn_model) if nn else df
    out_df: pd.DataFrame = xgboost_predict_from_trained_model(
        out_df, feature_cols, xgboost_model) if xg else df
    return out_df


def randomforest_predict_from_trained_model(df: pd.DataFrame, feature_cols, rf_model) -> pd.DataFrame:
    y_prob = rf_model.predict_proba(df[feature_cols])
    y_prob_interface = [p[1] for p in y_prob]
    df['randomforest'] = y_prob_interface
    return df


def logreg_predict_from_trained_model(df, feature_cols, logreg_model) -> pd.DataFrame:
    prediction = logreg_model.predict_proba(df[feature_cols])
    df["logisticregresion"] = [p[1] for p in prediction]
    return df


def linreg_predict_from_trained_model(df, feature_cols, linreg_model) -> pd.DataFrame:
    df["linearregression"] = linreg_model.predict(df[feature_cols])
    return df

# IN DEVELOPMENT 
# def nueralnet_predict_from_trained_model(df, feature_cols, nn_model) -> pd.DataFrame: 
#     # df["nueralnet"] = nn_model.predict(df[feature_cols])
#     return df


def xgboost_predict_from_trained_model(df, feature_cols, xgboost_model) -> pd.DataFrame:
    df["xgboost"] = xgboost_model.predict(df[feature_cols])
    return df

# Evan Edelstein
import argparse
import pathlib
import os
from .argscontainer import ArgsContainer
import sys

class DeprecateAction(argparse.Action):
    def __init__(self, *args, **kwargs):
        self.call_count = 0
        if 'help' in kwargs:
            kwargs['help'] = f'[DEPRECATED] {kwargs["help"]}'
        super().__init__(*args, **kwargs)


    def __call__(self, parser, namespace, values, option_string=None):
        if self.call_count == 0:
            sys.stderr.write(f"The option `{option_string}` is deprecated. It will be ignored.\n")
            sys.stderr.write(self.help + '\n')
            delattr(namespace, self.dest)
        self.call_count += 1

def userinterface() -> ArgsContainer:
    """
    See ReadMe for details
    """
    parser = argparse.ArgumentParser()
    parser.register('action', 'ignore', DeprecateAction)
    parser.add_argument('-i', '--inputfile',
                        default='input_data_all.csv', help='CSV Filename with columns: "residue","predus","ispred","dockpred","annotated". The column residue is of the form {residue number}_{PDB ID}.{chain}. The annotated column is 1 or interface residue and 0 for non-interface residue')
    parser.add_argument('-m', '--mode', dest="modeselection", choices=['predict', 'test', 'generate', 'cv', 'viz', "reprocess"], default='predict',
                        help="predict: Use pretrained model in input folder to predict on set.\ntest: generate a new rf model from a test set and train on a training set. Make sure --trainset and --testset is set.\ngenerate: generate a new rf model from a test set without predicting on any data.\ncv: perform cross validation. Use -cv to set test and train files.\nviz: perform data visualization with pymol use --results-df to specify output file from a previous prediction. reprocess: perform data analysis use --results-df to specify output file from a previous prediction")

    parser.add_argument('--trainset', default='train_set.csv', help='CSV Filename containing proteins for models to train on with columns: protein,size. The column protein is of the form {PDB ID}.{chain}')
    parser.add_argument('--testset', default='test_set.txt', help='Filename containing proteins for models to test on with columns: protein,size. The column protein is of the form {PDB ID}.{chain}')
    parser.add_argument("--rf-trees", dest='randomforest_parameter_trees', default=100, help='Scikit learn n_estimators parameter.')
    parser.add_argument("--rf-depth", dest='random_forest_parameter_depth',
                        default=None, help='Scikit learn max_depth parameter.')
    parser.add_argument("--rf-cpp", dest='random_forest_parameter_ccp', default=0.0, help='Scikit learn ccp_alpha parameter. ')
    parser.add_argument('-tv', '--tree_visualization',
                        action='store_true', help='Output Tree Visualization')
    parser.add_argument('-xg', '--xgboost',
                        action='store_true', help='Use XGboost')
    parser.add_argument('-nn', '--nuarelnet',
                        action='store_true', help='Use basic Neural network classifier Feature is IN DEVELOPMENT ')
    parser.add_argument("-p", '--pymol', dest='protein_visualization',
                        action='store_true', help='Output pymol images')
    parser.add_argument('--cutoffs', default='cutoffs.csv', help='CSV Filename containing length of interface or precalculated cutoff for each protein. File should have columns: Protein,surface res,cutoff res,annotated res ')
    parser.add_argument('--autocutoff', default='15', help='If no cutoff file is used this sets the default interface cutoff value.')
    parser.add_argument('--model-name', default='model', help='Name of models to import/export.')
    parser.add_argument('-of', '--outputfolder', default='output', help='- Directory to place output of ISPIP')
    parser.add_argument('-if', '--inputfolder', default='', help='Directory containing trained models. This folder should contain .joblib files to use as model inputs.')
    parser.add_argument('-cv', '--cvfoldername', default='cv', help='Directory containing test and train sets for cross-validation. Same csv format as train/test. Filenames should start with train and test.')
    # TODO add this argscontainer
    parser.add_argument('--plot', dest="plotselection", choices=[
                        'plot', 'csv', 'both'], default='both', help="output pr and roc curve as csv, png or both")
    parser.add_argument('--results-df', dest="results_df_input" ,help="csv output from predict mode to reprocess")

    # DEPRECIATED
    parser.add_argument('-mode', '--modeselection', action="ignore", help=" Pleas use -m or --mode")
    parser.add_argument('-trainset', action="ignore" , help="Please use --trainset")
    parser.add_argument('-testset', action="ignore" , help="Please use --testset")
    parser.add_argument('-randomforest_parameter_trees', action="ignore" , help="Please use  --rf-trees")
    parser.add_argument('-random_forest_parameter_ccp',  action="ignore" , help="Please use --rf-cpp ")
    parser.add_argument('-random_forest_parameter_depth',action="ignore" , help="Please use --rf-depth")
    parser.add_argument('-pymol', "--protein_visualization", action="ignore" , help="Please use -p or --pymol")
    parser.add_argument('-cutoffs', action="ignore" , help="Please use -p or --cutoffs")
    parser.add_argument('-autocutoff', action="ignore" , help="Please use  --autocutoff")
    parser.add_argument('-model_name', action="ignore" , help="Please use --model-name")
    parser.add_argument('-plot', "--plotselection", action="ignore" , help="Please use --plot")
    

    args: argparse.Namespace = parser.parse_args()
    args_container: ArgsContainer = parse(args)
    return args_container


def parse(args: argparse.Namespace) -> ArgsContainer:
    # folder_path: pathlib.Path = pathlib.Path(__file__).parent.parent
    # args_container = ArgsContainer(args, folder_path)
    args_container = ArgsContainer(args)

    os.makedirs(args_container.output_path_dir, exist_ok=True)

    if args_container.mode == 'test':
        if (not os.path.isfile(args_container.test_proteins_file)) or (not os.path.isfile(args_container.train_proteins_file)):
            print(
                "train and or test sets are not set, random 80/20 distribution will be used")
            args_container.use_test_train_files = False

    if not os.path.isfile(args_container.cutoff_frame):
        print(
            f"cutoffs not found, a global cutoff of {args.autocutoff} residues will be used (this value can be changed with the -autocutoff flag")
        args_container.use_cutoff_from_file = False

    elif not os.path.isfile(args_container.input_frames_file):
        print("please include an input csv file")
        sys.exit(0)
    if args_container.nn:
        args_container.models_to_use.append("nueralnet")
    if args_container.xg:
        args_container.models_to_use.append("xgboost")

    return args_container

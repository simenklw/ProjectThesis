import numpy as np
import pandas as pd

def island_split(pheno:str = "mass") -> pd.DataFrame:
    data = pd.read_feather(f"data/processed/{pheno}BV.feather")  

    table = pd.crosstab(data["hatchisland"],data["hatchisland"])
    diag = np.diag(table)

    # get the indexes of the 8 highest values in the diagonal
    idx = np.argpartition(diag, -8)[-8:]

    #islands with highest number of individuals
    islands = table.index[idx]

    # get the column hatchisland and ringnr for the 8 islands with the highest number of individuals
    x = data.loc[data["hatchisland"].isin(islands),["hatchisland","ringnr"]]
    del data

    #Now we will manually make folds for a 8-fold cross validation based on the hatchisland column
    #The idea is to make 8 folds where each fold contains all individuals from one of the 8 islands with the highest number of individuals.
    #This way we can make sure that the model is trained on all islands and tested on all islands. The size of the folds will be different.

    #We will make a new column in the data frame called "fold" which will contain the fold number for each individual.
    #The fold number will be between 1 and 8.


    #initialize the fold column with zeros
    x["fold"] = 0

    #initialize the fold number
    fold = 1

    #loop over the 8 islands with the highest number of individuals
    for island in islands:
        #get the indexes of the individuals from the current island
        idx = x.loc[x["hatchisland"] == island].index
        #set the fold number for the individuals from the current island
        x.loc[idx,"fold"] = fold
        #increment the fold number
        fold += 1
    
    #save the data frame with the fold column
    for i in range(1,9):
        fold = x[x["fold"] == i]
        fold = fold["ringnr"]
        fold.to_csv(f"data/interim/{pheno}BV_island_8fold_{i}.csv", index=False)

def random_split(pheno:str = "mass", num_folds:int = 5, seed: int = 42) -> pd.DataFrame:
    data = pd.read_feather(f"data/processed/{pheno}BV.feather")  
    x = data[["ringnr","hatchisland"]]
    del data
    np.random.seed(seed)
    #Split the data into num_folds folds with random sampling, each fold should have the same number of individuals
    #The folds should be saved in a new column called "fold" in the data frame

    #x = x.sample(frac=1, random_state=seed).reset_index(drop=True)
    # Number of folds
    n_folds = num_folds

    # Create an array with fold numbers
    folds = np.array([i % n_folds for i in range(len(x))])


    # Shuffle the fold assignment array
    np.random.shuffle(folds)

    # Assign the fold numbers to a new column in the dataframe
    x['fold'] = folds
    
    return x

    #for i in range(n_folds):
    #    fold = x[x["fold"] == i]
    #    fold = fold["ringnr"]
    #    fold.to_csv(f"data/interim/folds/random_{n_folds}fold_{i+1}_{pheno}.csv", index=False)



def prep_data_before_train(data: pd.DataFrame, phenotype: str) -> tuple:
    """
    prepare data for training, returns target vector and covariate-matrix and ringnrs for grouping
    :param data: all data from dataloader script
    :param phenotype: the phenotype to be predicted
    :return: X, Y, ringnrs, X contain covariates and Y contains the phenotype
    """
    ringnrs = data.ringnr
    mean_pheno = data.mean_pheno
    X = data.drop(
        columns=[
            "ID",
            "mass",
            "tarsus",
            "ringnr",
            "mean_pheno",
            "IID",
            "FID",
            "MAT",
            "PAT",
            "SEX",
            "PHENOTYPE",
        ],
        errors="ignore",
    )

    try:
        Y = data.loc[:, phenotype]
    except KeyError:
        try:
            Y = data.ID
        except AttributeError:
            Y = data.mass
    del data
    try:
        X.hatchyear = X.hatchyear.astype("int")
        X.island_current = X.island_current.astype("int")
    except AttributeError:
        pass
    Y = (Y - np.mean(Y))/np.std(Y)
    X = X.fillna(0)
    X = X.T
    X = X.astype("int")
    X = X.T 
    return X, Y, ringnrs, mean_pheno

def subset(X:pd.DataFrame, seed:int = 42, num_snps:int = 20000) -> pd.DataFrame:
    """
    Function to subset the data into a smaller dataset for testing purposes, 
    sampling only a subset of the features
    """
    np.random.seed(seed)
    sample_columns = np.random.choice(X.columns[:-5], num_snps, replace = False)
    sample_columns = np.append(sample_columns, ["ringnr","fold", "ID", "mean_pheno"])

    return X[sample_columns]

    
if __name__ == "__main__":
    #island_split()
    random_split("wing", 10, 42)
import optuna
import logging
from sklearn.model_selection import cross_val_score, KFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
import time


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Define objective function for Optuna
def objective(trial):
    
    logger.info(f"Starting trial {trial.number}")

    # Define hyperparameter search space
    n_estimators = trial.suggest_int('n_estimators', 50, 300)
    max_depth = trial.suggest_int('max_depth', 2, 20)
    min_samples_split = trial.suggest_int('min_samples_split', 2, 20)
    
    # Model to optimize
    model = RandomForestRegressor(
        n_estimators=n_estimators, 
        max_depth=max_depth, 
        min_samples_split=min_samples_split
    )
    
    # Cross-validation setup
    kfold = KFold(n_splits=5, shuffle=True, random_state=42)
    start_time = time.time()
    # Evaluate with cross-validation
    scores = cross_val_score(model, X, y, cv=kfold, scoring='neg_mean_squared_error')
    end_time = time.time()
    
    logger.info(f"Finished trial {trial.number} in {end_time - start_time:.2f} seconds")
    # Return the mean score (negative MSE, so Optuna minimizes)
    return scores.mean()


if __name__ == '__main__':
    # Load or generate your dataset
    X, y = make_regression(n_samples=1000, n_features=5)

    # Create Optuna study
    study = optuna.create_study(direction='maximize')

    n_jobs = 4  
    
    study.optimize(objective, n_trials=10, n_jobs=n_jobs)
    plot_slice = optuna.visualization.plot_slice(study)
    plot_slice.show()

    # Best hyperparameters
    print('Best hyperparameters:', study.best_params)
    print(X)
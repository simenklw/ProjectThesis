import optuna
import logging
from sklearn.model_selection import cross_val_score, KFold
from sklearn.datasets import make_regression
from multiprocessing import Pool
import time
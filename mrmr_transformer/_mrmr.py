import pandas as pd
from mrmr import mrmr_classif
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted


class MRMRTransformer(TransformerMixin, BaseEstimator):
    """MRMR feature selection for a classification task.

    Parameters
    ----------
    n_features : int
        The number of features of the data passed to :meth:`fit`.

    n_jobs : int, default=None
        The number of parallel jobs to run for neighbors search.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.
        Doesn't affect :meth:`fit` method.

    Attributes
    ----------
    selected_features_ : array of shape (n_features,)
        Selected features using MRMR algorithm.
    """

    def __init__(self, n_features=5, n_jobs=-1, verbose=None):
        self.n_features = n_features
        self.n_jobs = n_jobs
        self.verbose = verbose

    def fit(self, X, y=None):
        """A reference implementation of a fitting function for a transformer.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The training input samples.
        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        self : object
            Returns self.
        """
        # self._validate_params()
        _X = pd.DataFrame(X)
        _y = pd.DataFrame(y)

        if _X.empty | _y.empty:
            raise ValueError("Empty data is used to fit.")

        self.selected_features_ = mrmr_classif(
            X=_X,
            y=_y,
            K=self.n_features,
            n_jobs=self.n_jobs,
            show_progress=self.verbose,
        )

        # Return the transformer
        return self

    def transform(self, X):
        """A reference implementation of a transform function.

        Parameters
        ----------
        X : {array-like, sparse-matrix}, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        X_transformed : array, shape (n_samples, n_features)
            The array containing the element-wise square roots of the values
            in ``X``.
        """
        # Check is fit had been called
        check_is_fitted(self, "selected_features_")

        # Input validation
        # X = check_array(X, accept_sparse=True)

        # Check that the input is of the same shape as the one passed
        # during fit.
        # if X.shape[1] != self.n_features_:
        #     raise ValueError('Shape of input is different from what was seen'
        #                      'in `fit`')

        return X[self.selected_features_]

    # def get_features(self):
    #     return self.selected_features_

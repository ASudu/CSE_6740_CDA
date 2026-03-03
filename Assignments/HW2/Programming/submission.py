from typing import Tuple, Union, Callable
import numpy as np


def initialize_parameters(
    X: np.ndarray, k: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return initial values for training of the GMM
    Set component mean to a random
    pixel's value (without replacement),
    based on the mean calculate covariance matrices,
    and set each component mixing coefficient (PIs)
    to a uniform values
    (e.g. 4 components -> [0.25,0.25,0.25,0.25]).

    params:
    X = numpy.ndarray[numpy.ndarray[float]] - m x n
    k = int

    returns:
    (MU, SIGMA, PI)
    MU = numpy.ndarray[numpy.ndarray[float]] - k x n
    SIGMA = numpy.ndarray[numpy.ndarray[numpy.ndarray[float]]] - k x n x n
    PI = numpy.ndarray[float] - k x 1
    """
    # Mixing coefficients
    PI = np.full((k, 1), 1 / k)

    # Means
    # As per docstring, we set it to a random pixel's value, without replacement
    random_indices = np.random.choice(X.shape[0], size=k, replace=False)
    MU = X[random_indices]

    # Covariance matrices
    SIGMA = compute_sigma(X, MU)

    return MU, SIGMA, PI


def compute_sigma(X: np.ndarray, MU: np.ndarray) -> np.ndarray:
    """
    Calculate covariance matrix, based in given X and MU values

    params:
    X = numpy.ndarray[numpy.ndarray[float]] - m x n
    MU = numpy.ndarray[numpy.ndarray[float]] - k x n

    returns:
    SIGMA = numpy.ndarray[numpy.ndarray[numpy.ndarray[float]]] - k x n x n
    """
    m,n = X.shape
    k = MU.shape[0]

    # Initialize SIGMA with zeros
    SIGMA = np.zeros((k, n, n))
    for i in range(k):
        # Center the data around ith mean
        centered = X - MU[i]

        # Compute covariance matrix for the ith component
        SIGMA[i] = (1 / m) * np.dot(centered.T, centered)
    
    return SIGMA


def prob(x: np.ndarray, mu: np.ndarray, sigma: np.ndarray) -> Union[float, np.ndarray]:
    """Calculate the probability of x (a single
    data point or an array of data points,
    you have to take care both cases) under the
    component with the given mean and covariance.
    The function is intended to compute multivariate
    normal distribution, which is given by N(x;MU,SIGMA).

    params:
    x = numpy.ndarray[float] or numpy.ndarray[numpy.ndarray[float]]
    mu = numpy.ndarray[float]
    sigma = numpy.ndarray[numpy.ndarray[float]]

    returns:
    probability = float or numpy.ndarray[float]
    """
    
    if x.ndim == 1:
        x = x.reshape(1, -1)  # Convert to 2D array with shape (1, n)
    elif x.ndim == 2 and x.shape[1] != mu.shape[0]:
        raise ValueError("The number of features in x must match the length of mu.")

    num_samples, dim = x.shape

    # Normalization constant
    norm_const = 1.0 / (np.power(2 * np.pi, dim / 2) * np.sqrt(np.linalg.det(sigma)))

    def single_sample_prob(x_i):
        """
        Calculate the probability of a single data point x_i under the multivariate normal distribution

        Args:
            x_i (np.ndarray): A single data point (row vector: 1 x n)

        Returns:
            float: The probability of x_i under the multivariate normal distribution
        """
        # Center the data point
        centered = x_i - mu # shape: (1,n)

        # Exponent term
        exponent = -0.5 * np.dot(centered, np.linalg.solve(sigma, centered.T)) # shape: (1,1)

        return norm_const * np.exp(exponent)

    if num_samples == 1:
        return single_sample_prob(x[0])
    else:
        return np.array([single_sample_prob(x_i) for x_i in x])


def E_step(
    X: np.ndarray, MU: np.ndarray, SIGMA: np.ndarray, PI: np.ndarray, k: int
) -> np.ndarray:
    """
    E-step - Expectation
    Calculate responsibility for each
    of the data points, for the given
    MU, SIGMA and PI.

    params:
    X = numpy.ndarray[numpy.ndarray[float]] - m x n
    MU = numpy.ndarray[numpy.ndarray[float]] - k x n
    SIGMA = numpy.ndarray[numpy.ndarray[numpy.ndarray[float]]] - k x n x n
    PI = numpy.ndarray[float] - k x 1
    k = int

    returns:
    responsibility = numpy.ndarray[numpy.ndarray[float]] - k x m
    """
    num_samples, _ = X.shape
    resp = np.zeros((k, num_samples))  # shape: (k, m)

    for i in range(k):
        likelihoods = prob(X, MU[i], SIGMA[i])  # shape: (m,)
        resp[i] = PI[i] * likelihoods  # shape: (m,)

    # Normalize responsibilities
    resp = resp / np.sum(resp, axis=0)

    return resp

def M_step(
    X: np.ndarray, r: np.ndarray, k: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    M-step - Maximization
    Calculate new MU, SIGMA and PI matrices
    based on the given responsibilities.

    params:
    X = numpy.ndarray[numpy.ndarray[float]] - m x n
    r = numpy.ndarray[numpy.ndarray[float]] - k x m
    k = int

    returns:
    (new_MU, new_SIGMA, new_PI)
    new_MU = numpy.ndarray[numpy.ndarray[float]] - k x n
    new_SIGMA = numpy.ndarray[numpy.ndarray[numpy.ndarray[float]]] - k x n x n
    new_PI = numpy.ndarray[float] - k x 1
    """

    # New mixing coefficients
    new_PI = np.mean(r, axis=1)  # shape: (k, 1)'

    # New means
    new_MU = (r @ X) / np.sum(r, axis=1, keepdims=True)  # shape: (k, n)

    # New covariance matrices
    new_SIGMA = np.zeros((k, X.shape[1], X.shape[1]))  # shape: (k, n, n)
    for i in range(k):
        centered = X - new_MU[i]  # shape: (m, n)
        weighted_centered = r[i][:, np.newaxis] * centered  # shape: (m, n)
        new_SIGMA[i] = (weighted_centered.T @ centered) / np.sum(r[i])  # shape: (n, n)
    
    return new_MU, new_SIGMA, new_PI


def likelihood(
    X: np.ndarray, MU: np.ndarray, SIGMA: np.ndarray, PI: np.ndarray, k: int
) -> float:
    """Calculate a log likelihood of the
    trained model based on the following
    formula for posterior probability:

    log(Pr(X | mixing, mean, stdev)) = sum((n=1 to N), log(sum((k=1 to K),
                                      mixing_k * N(x_n | mean_k,stdev_k))))

    Make sure you are using natural log, instead of log base 2 or base 10.

    params:
    X = numpy.ndarray[numpy.ndarray[float]] - m x n
    MU = numpy.ndarray[numpy.ndarray[float]] - k x n
    SIGMA = numpy.ndarray[numpy.ndarray[numpy.ndarray[float]]] - k x n x n
    PI = numpy.ndarray[float] - k x 1
    k = int

    returns:
    log_likelihood = float
    """
    num_samples, _ = X.shape
    log_likelihood = 0.0

    for n in range(num_samples):
        likelihood_n = 0.0
        for i in range(k):
            likelihood_n += PI[i] * prob(X[n], MU[i], SIGMA[i])
        log_likelihood += np.log(likelihood_n)

    return log_likelihood


def default_convergence(prev_likelihood, new_likelihood, conv_ctr, conv_ctr_cap=10):
    """
    Default condition for increasing
    convergence counter:
    new likelihood deviates less than 10%
    from previous likelihood.

    params:
    prev_likelihood = float
    new_likelihood = float
    conv_ctr = int
    conv_ctr_cap = int

    returns:
    conv_ctr = int
    converged = boolean
    """
    increase_convergence_ctr = (
        abs(prev_likelihood) * 0.9 < abs(new_likelihood) < abs(prev_likelihood) * 1.1
    )

    if increase_convergence_ctr:
        conv_ctr += 1
    else:
        conv_ctr = 0

    return conv_ctr, conv_ctr > conv_ctr_cap


def train_model(
    X: np.ndarray,
    k: int,
    convergence_function: Callable = default_convergence,
    initial_values: Tuple[np.ndarray, np.ndarray, np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Train the mixture model using the
    expectation-maximization algorithm.
    E.g., iterate E and M steps from
    above until convergence.
    If the initial_values are None, initialize them.
    Else it's a tuple of the format (MU, SIGMA, PI).
    Convergence is reached when convergence_function
    returns terminate as True,
    see default convergence_function example
    above.

    params:
    X = numpy.ndarray[numpy.ndarray[float]] - m x n
    k = int
    convergence_function = func
    initial_values = None or (MU, SIGMA, PI)

    returns:
    (new_MU, new_SIGMA, new_PI, responsibility)
    new_MU = numpy.ndarray[numpy.ndarray[float]] - k x n
    new_SIGMA = numpy.ndarray[numpy.ndarray[numpy.ndarray[float]]] - k x n x n
    new_PI = numpy.ndarray[float] - k x 1
    responsibility = numpy.ndarray[numpy.ndarray[float]] - k x m
    """
    # Initialize parameters
    if initial_values is None:
        MU, SIGMA, PI = initialize_parameters(X, k)
    else:
        MU, SIGMA, PI = initial_values

    prev_likelihood = None
    conv_ctr = 0

    while True:
        # E-step
        responsibility = E_step(X, MU, SIGMA, PI, k)

        # M-step
        MU, SIGMA, PI = M_step(X, responsibility, k)

        # Check for convergence
        new_likelihood = likelihood(X, MU, SIGMA, PI, k)
        print(f"Log Likelihood: {new_likelihood}")

        if prev_likelihood is not None:
            conv_ctr, converged = convergence_function(prev_likelihood, new_likelihood, conv_ctr)
            if converged:
                print("Convergence reached.")
                break

        prev_likelihood = new_likelihood

    return MU, SIGMA, PI, responsibility

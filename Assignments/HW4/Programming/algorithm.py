import matplotlib.pyplot as plt
import numpy as np


def algo(q, Y):
    # init
    p = 0.0
    fig, ax = plt.subplots()

    # TODO implement your algorithm and return the (i) prob p and (ii) a matplotlib Figure object for the plot
    # 0 = bad, 1 = good for indexing hidden states
    # 0 = low, 1 = high for indexing observations
    pi = np.array([[0.8], [0.2]]) # initial probabilities
    a = np.array([[0.8, 0.2], [0.2, 0.8]]) # transition probabilities
    b = np.array([[q, 1-q], [1-q, q]]) # emission probabilities
    y = Y[:, 0] # observed price moves

    alphas = np.zeros((2, len(Y)))
    betas = np.zeros((2, len(Y)))

    # Initialize
    # (1+Y[0])//2 converts Y[0] from {-1, 1} to {0, 1}
    alphas[:, 0] = pi[:, 0] * b[:, (1+y[0])//2] # probability of starting in each state and observing Y[0]
    betas[:, -1] = 1 # probability of observing the empty sequence from any state is 1

    # Forward pass
    for t in range(1, len(Y)):
        for j in range(2): # for each state
            # Forward pass
            alphas[j, t] = b[j, (1+y[t])//2]*(alphas[:, t-1] @ a[:, j])

    # total probability of observing Y is the sum of probabilities of being in each state at the end and observing Y
    p_Y = np.sum(alphas[:, -1]) # total probability of observing Y

    # Backward pass
    for t in range(len(Y)-2, -1, -1):
        for i in range(2): # for each state
            # Backward pass
            betas[i, t] = a[i, :] @ (b[:, (1+y[t+1])//2] * betas[:, t+1])
    
    probs = [alphas[1,t]*betas[1,t]/p_Y for t in range(len(Y))] # probability of being in good state at time t given Y

    ax.plot(probs)
    ax.set_xlabel('Week (t)')
    ax.set_ylabel('P(Good State) in week t')
    ax.set_title(f'Probability of Being in Good State Over Time (q = {q})')
    plt.tight_layout()

    p = probs[-1] # probability of being in good state at the end of the sequence
    return p, fig
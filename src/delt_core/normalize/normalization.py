import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, pearsonr, spearmanr, halfnorm


def compute_probs(
        num_blocks: int,
        num_samples: int,
        affinities: np.ndarray,
) -> np.ndarray:
    epsilon = 1e-12
    probs = affinities / affinities.sum()
    distribution = np.random.choice(range(num_blocks), num_samples, p=probs)
    counts = np.bincount(distribution, minlength=num_blocks)
    return counts / num_samples + epsilon


def plot_distribution(
        probs: list,
        num_blocks: int,
        title: str,
        **kwargs,
) -> None:
    plt.bar(range(num_blocks), probs, **kwargs)
    plt.title(title)
    plt.xlabel('Building Block')
    plt.ylabel('Probability')


def compute_geometric_mean(
        data: list,
) -> float:
    return np.exp(np.log(data).mean())


def deseq(
        data: np.ndarray,
) -> np.ndarray:
    norm_counts = []
    for counts in data.T:
        geometric_mean = compute_geometric_mean(counts)
        count = counts[0]
        size_factor = np.median(count / geometric_mean)
        norm_counts += [count / size_factor]
    return np.array(norm_counts)


def main():

    num_blocks = int(1e5)
    num_samples = int(1e6)

    # Step 1
    uniform_affinities = np.ones(num_blocks)
    uniform_probs = compute_probs(num_blocks, num_samples, uniform_affinities)

    # Step 2
    neg_affinities = np.random.lognormal(size=num_blocks, sigma=0.1)
    neg_probs = compute_probs(num_blocks, num_samples, neg_affinities)

    # Step 3
    pos_affinities = halfnorm.rvs(size=num_blocks, scale=3) + neg_affinities
    pos_probs = compute_probs(num_blocks, num_samples, pos_affinities)

    # Normalization
    neg_counts = neg_probs * num_samples
    pos_counts = pos_probs * num_samples
    # counts = np.stack((pos_counts, neg_counts))
    # norm_counts = deseq(counts)
    norm_counts = pos_counts / uniform_probs / neg_probs
    norm_probs = norm_counts / norm_counts.sum()

    # Correlation
    rank_neg_affinities = rankdata(neg_affinities, method='ordinal')
    rank_pos_affinities = rankdata(pos_affinities, method='ordinal')
    rank_norm_probs = rankdata(norm_probs, method='ordinal')
    rank_pos_probs = rankdata(pos_probs, method='ordinal')

    pearson_corr_neg_pos, _ = pearsonr(neg_affinities, pos_probs)
    pearson_corr_neg_norm, _ = pearsonr(neg_affinities, norm_probs)
    pearson_corr_pos_pos, _ = pearsonr(pos_affinities, pos_probs)
    pearson_corr_pos_norm, _ = pearsonr(pos_affinities, norm_probs)

    spearman_corr_neg_pos, _ = spearmanr(rank_neg_affinities, rank_pos_probs)
    spearman_corr_neg_norm, _ = spearmanr(rank_neg_affinities, rank_norm_probs)
    spearman_corr_pos_pos, _ = spearmanr(rank_pos_affinities, rank_pos_probs)
    spearman_corr_pos_norm, _ = spearmanr(rank_pos_affinities, rank_norm_probs)
    
    print('neg aff - pos probs:\t', pearson_corr_neg_pos, spearman_corr_neg_pos)
    print('neg aff - norm probs:\t', pearson_corr_neg_norm, spearman_corr_neg_norm)
    print('pos aff - pos probs:\t', pearson_corr_pos_pos, spearman_corr_pos_pos)
    print('pos aff - norm probs:\t', pearson_corr_pos_norm, spearman_corr_pos_norm)
    
    exit()


    mean = np.mean(norm_probs)
    std = np.std(norm_probs)
    p = np.array(norm_probs)
    threshold = mean + 2 * std
    mask = p > threshold
    color = ['red' if significant else 'blue' for significant in mask]


    plt.figure(figsize=(18, 6))

    plt.subplot(1, 4, 1)
    plot_distribution(uniform_probs, num_blocks, 'Uniform Control Distribution')

    plt.subplot(1, 4, 2)
    plot_distribution(neg_probs, num_blocks, 'Negative Control Distribution')

    plt.subplot(1, 4, 3)
    plot_distribution(pos_probs, num_blocks, 'Positive Control Distribution')

    plt.subplot(1, 4, 4)
    plot_distribution(norm_probs, num_blocks, 'Normalized Distribution', color=color)
    plt.axhline(y=1/num_blocks, color='black', linestyle='--', linewidth=1)

    plt.tight_layout()
    # plt.savefig('normalization.png', dpi=200)
    plt.show()



if __name__ == '__main__':

    main()



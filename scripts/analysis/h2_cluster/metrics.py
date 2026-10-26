"""H2-only metrics: exactly two bonds per molecular pair."""
from .stats import ratio

def h2_metrics(h2, n, e, volume):
    if not (isinstance(h2, int) and 0 <= 2*h2 <= e <= 2*n and n > 0 and volume > 0):
        raise ValueError('Invalid H2, N, E, or volume')
    return dict(H2_per_volume=h2/volume,H2_per_bond=ratio(h2,e),
                H2_endpoint_fraction_bonded=ratio(2*h2,e),
                H2_endpoint_fraction_all=h2/n)



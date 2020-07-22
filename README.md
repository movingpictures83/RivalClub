# RivalClub
# Language: Python
# Input: CSV
# Output: CSV
# Tested with: PluMA 1.1, Python 3.6
# Dependency: PluMA AffinityPropagation plugin

Compute rival clubs within a network (Fernandez et al, 2015).
A rival club is a tightly connected component with a large amount
of negative connections with another tightly connected component.

The input network should be signed and weighted, in the format
of a CSV file where nodes correspond to rows and columns and entry
(i, j) is the weight of the edge from i to j.

Rival clubs are output to a CSV file, with each club delineated by
"","x".

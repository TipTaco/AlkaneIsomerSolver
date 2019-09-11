# Author Adrian Shedley,
# Date 11 Sep 2019
from timeit import default_timer as timer
import numpy as np


def multisets(n, k):
    """Multiset combinatory maths. Takes n possible ways and k repetitions"""
    return choose(n+k-1, k)


def choose(n, k):
    """N choose K"""
    c = 1
    for i in range(n-k+1, n+1):
        c *= i
    return int(c/factorial(k))


def factorial(n):
    """Factorial n is n!"""
    f = 1
    for i in range(2, n+1):
        f *= i
    return f


def next_partition(parts, sumN, maxN, numParts):
    """Takes the current partition parts and generates the next partition. Will give numParts parts and sum
    to sumN with limit maxN"""
    if parts[numParts - 1] == 0:
        for i in reversed(range(0, numParts-1+1)):
            parts[i] = min(sumN, maxN)
            sumN -= parts[i]
            if sumN <= 0: break
        return True

    for i in range(1, numParts):
        if parts[i] - parts[0] >= 2:
            for j in reversed(range(0, i-1+1)):
                if (parts[i]-parts[j]) >= 2:
                    parts[i] -= 1
                    parts[j] += 1
                    r = 0

                    for k in reversed(range(0, j-1+1)):
                        r += parts[k]
                        parts[k] = 0
                        if k < 0: break

                    k = j
                    while r > 0:
                        x = min(r, parts[k+1]-parts[k])
                        parts[k] += x
                        k -= 1
                        r -= x
                    return True
    return False


def count_occ(n, arr):
    """Count occurrences of n in arr"""
    count = 0
    for i in range(len(arr)):
        if n == arr[i]:
            count += 1
    return count


def r_solver(partition, r_table):
    """Returns possibilities for the sidechain R denoted by partition, given the current r_table"""
    completed = set()
    output = 1

    for part in partition:
        if part not in completed:
            n = r_table[part]
            k = count_occ(part, partition)
            if k == 1:
                output *= n
            else:
                output *= multisets(n, k)
            completed.add(part)

    return output


def calc_rtable(highest):
    """Populates the r_table and returns it up to highest+1 entries for the sidechain R combinations"""
    parts = 3
    arr = np.zeros(parts, dtype=np.uint16)
    r_table = list()
    r_table.append(1) # R0 = 1
    r_table.append(1) # R1 = 1

    for i in range(2, highest):
        sum = 0
        arr = np.zeros(parts, dtype=np.uint16)
        while next_partition(arr, i-1, highest, parts):
            sum += r_solver(arr, r_table)
        r_table.append(sum)

    return r_table


def num_isomers(n, r_table):
    """Calculates the number of isomers for an alkane of length n, given an appropriate precomputed r_table"""
    half_chain = int(n/2)
    sum = multisets(r_table[half_chain], 2)  # base case for R(n)-R(n) direct bonding

    if n % 2 == 1:
        # modified case for odd Carbon count. Can't use R(n)-R(n) symmetry so we must first convert to the form
        # R(n)-C-R(n) and then account for the set overlap with inclusion-exclusion principle
        sum = r_table[half_chain] * (r_table[half_chain+1]+1) - sum

    # Second half of this method generates the missing combinations that were not allowed with R(n)-R(n) formulation
    #          |            The number of carbons on either side of this split line, denoted |, is locked. The missing
    #  C-C-C-C-|-C-C-C-C    carbons are accounted for by regenerating all the partitions that sum to n-1 and have a
    #          |            maximum part size of half_chain-1. This will give a central carbon with FOUR bonds as shown
    #
    #          R1           This is generated as one of the partitions [1, 1, 2, 3] for C8 (octane)
    #       R3 C R2
    #          R1
    #
    parts = 4
    arr = np.zeros(parts, dtype=np.uint16)
    while next_partition(arr, n-1, half_chain-1, parts):
        sum += r_solver(arr, r_table)

    return sum


### BEGIN ###
ALKANE_SIZE = 20

start = timer()
r_table = calc_rtable(ALKANE_SIZE)
isomers = num_isomers(ALKANE_SIZE, r_table)
end = timer()

print('> Completed R Table:', r_table)
print('> Alkane of size', ALKANE_SIZE, 'has', isomers, 'isomers')
print(' > Time elapsed:', round((end - start)*1000, 2), 'ms')


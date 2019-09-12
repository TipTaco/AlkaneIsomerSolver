# Author Adrian Shedley,
# Date 11 Sep 2019
from timeit import default_timer as timer
import numpy as np
import math


def multisets(n, k):
    """Multiset combinatory maths. Takes n possible ways and k repetitions"""
    return choose(n+k-1, k)


def choose(n, k):
    """N choose K"""
    c = 1
    for i in range(n-k+1, n+1):
        c = i * c
    return c // math.factorial(k)


def trickle(parts, sumN, maxN, numParts, start):
    # sumN -= sum(parts[start+1:numParts])
    for i in range(start, -1, -1):
        if sumN > 0:
            parts[i] = min(sumN, maxN)
            sumN -= parts[i]
        else:
            parts[:i+1] = [0] * (i+1)
            break

    return True

def nxt_partition(parts, sumN, maxN, numParts):
    """Takes the current partition parts and generates the next partition. Will give numParts parts and sum
    to sumN with limit maxN"""
    if parts[numParts - 1] == 0:
        return trickle(parts, sumN, maxN, numParts, numParts-1)

    pSum = 0
    for i in range(1, numParts):
        pSum += parts[i - 1]
        if parts[i] - parts[0] < 2:
            continue
        for j in range(1, i+1):
            if parts[i] - parts[i-j] >= 2:
                parts[i] -= 1
                pSum += 1
                if i == 1:
                    parts[i-1] += 1
                else:
                    trickle(parts, pSum, parts[i], numParts, i-1)
                return True

    return False


def next_partition(parts, sumN, maxN, numParts):
    """Takes the current partition parts and generates the next partition. Will give numParts parts and sum
    to sumN with limit maxN"""
    if parts[numParts - 1] == 0:
        for i in range(numParts-1, -1, -1):
            parts[i] = min(sumN, maxN)
            sumN -= parts[i]
            if sumN <= 0: break
        return True

    for i in range(1, numParts):
        if parts[i] - parts[0] >= 2:
            for j in range(i-1, -1, -1):
                if (parts[i]-parts[j]) >= 2:
                    parts[i] -= 1
                    parts[j] += 1
                    r = 0

                    for k in range(0, j-1+1):
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


def partition_master(numN, parts):
    arr = [0] * parts
    partition_table = list()
    while next_partition(arr, numN, numN, parts):
        partition_table.append(arr.copy())
    return partition_table


def partition_table(master_table, numN, maxN):
    partition_table = master_table
    offset = master_table[0][2] - numN + 1
    new_table = list()

    for se in partition_table:
        se[2] -= offset
        #print(se, sum(se[:]), max(se[:]), se[2], se[1])
        if se[2] >= se[1]:
            new_table.append(se.copy())

    return new_table


def calc_rtable(highest):
    """Populates the r_table and returns it up to highest+1 entries for the sidechain R combinations"""
    parts = 3
    r_table = [1, 1]  # , 1, 2, 4, 8, 17, 39, 89, 211, 507, 1099, 2454, 5461]  # R(0) to R(13)

    for i in range(2, highest):
        sum = 0
        arr = [0] * parts
        while next_partition(arr, i-1, highest, parts):
            sum += r_solver(arr, r_table)
        r_table.append(sum)

    return r_table


def calc_rtable_mod(master_table, highest):
    """Populates the r_table and returns it up to highest+1 entries for the sidechain R combinations"""
    parts = 3
    r_table = [1, 1]  # , 1, 2, 4, 8, 17, 39, 89, 211, 507, 1099, 2454, 5461]  # R(0) to R(13)

    for i in range(2, highest):
        sum = 0
        p_tab = partition_table(master_table, i, highest)
        for arr in p_tab:
            sum += r_solver(arr, r_table)
        r_table.append(sum)

    return r_table


def r_solver(partition, r_table):
    """Returns possibilities for the sidechain R denoted by partition, given the current r_table"""
    completed = {0, 1, 2}  # The numbers from this partition that have been done
    output = 1  # Start with a 1 so we can multiply it
    adder = completed.add

    for part in partition:  # example '2' from [0,1,2,3]
        if part not in completed:  # if we havent done R(part)
            n = r_table[part]
            k = partition.count(part)
            if k == 1:
                output *= n
            else:
                output *= multisets(n, k)
            adder(part)

    return output


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
    arr = [0] * parts  # np.zeros(parts, dtype=np.uint32)
    while next_partition(arr, n-1, half_chain-1, parts):
        sum += r_solver(arr, r_table)

    return sum

def num_isomers_mod(master_table, n, r_table):
    """Calculates the number of isomers for an alkane of length n, given an appropriate precomputed r_table"""
    half_chain = int(n/2)
    sum = multisets(r_table[half_chain], 2)  # base case for R(n)-R(n) direct bonding

    if n % 2 == 1:
        sum = r_table[half_chain] * (r_table[half_chain+1]+1) - sum

    parts = 4
    partitions = partition_table(master_table, n-1, half_chain-1)
    for arr in partitions:
        sum += r_solver(arr, r_table)

    return sum

def task1():
    start = timer()
    r_table = calc_rtable(int(ALKANE_SIZE / 2 + 1))
    isomers = num_isomers(ALKANE_SIZE, r_table)

    end = timer()

    print('> Completed R Table:', r_table)
    print('> Alkane of size', ALKANE_SIZE, 'has', isomers, 'isomers')
    print(' > Time elapsed:', round((end - start) * 1000.0, 4), 'ms')
    print(" ")

    return isomers, r_table


def task1_mod():

    start = timer()
    m_tab = partition_master(ALKANE_SIZE//2, 3)
    #p_tab = partition_table(m_tab, 11, 11)
    # table_4way = partition_master(20+1, 4)

    r_table = calc_rtable_mod(m_tab, int(ALKANE_SIZE / 2 + 1))
    isomers = num_isomers(ALKANE_SIZE, r_table)  # num_isomers_mod(table_4way, ALKANE_SIZE, r_table)

    end = timer()

    print('> [MOD] R Table:', r_table)
    print('> Alkane of size', ALKANE_SIZE, 'has', isomers, 'isomers')
    print(' > Time elapsed:', round((end - start) * 1000.0, 4), 'ms')
    print(" ")

    return isomers, r_table


def task2():
    divis = 3

    start = timer()
    arr = [0] * divis  # np.zeros(parts, dtype=np.uint32)
    count1 = 0
    while next_partition(arr, 200, 100, divis):
        # print(arr)
        count1 += 1

    mid = timer()
    arr = [0] * divis  # np.zeros(parts, dtype=np.uint32)
    count2 = 0
    while nxt_partition(arr, 8, 11, divis):
        print(arr)
        count2 += 1

    end = timer()
    print('Alec', round((mid - start) * 1000.0, 4), 'ms vs Adrian', round((end - mid) * 1000.0, 4))
    print(count1, count2)


def task3():

    p1 = timer()
    m_tab = partition_master(11, 3)

    p2 = timer()
    p_tab = partition_table(m_tab, 11, 11)

    p3 = timer()

    divis = 3
    arr = [0] * divis  # np.zeros(parts, dtype=np.uint32)
    count2 = 0
    while nxt_partition(arr, 11, 11, divis):
        count2 += 1

    p4 = timer()

    print('Master table', m_tab)
    print(len(p_tab), p_tab)
    print('times', round((p2 - p1) * 1000.0, 4), 'ms and', round((p3 - p2) * 1000.0, 4), 'ms and ', round((p4 - p3) * 1000.0, 4), 'ms for std method')


### BEGIN ###
ALKANE_SIZE = 20
#task2()
#task3()

isomers, r_table = task1()
task1_mod()



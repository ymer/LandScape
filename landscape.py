from __future__ import division
import math
import gzip
from math import log, pow
import sys
import argparse
from collections import namedtuple
import pandas as pd
import plink
import random
import numpy as np

Landscape = namedtuple('Landscape', 'A s e Y position_Y sections section_indices maximal_segments')
MaximalSegment = namedtuple('Maximal', 's e sec_s sec_e Y independent')


def landscape_generator_plink():
    inds, snps, genos = plink.read_plink(settings.infile)
    random.seed(settings.random_seed)
    while 1:
        inds.pheno = np.random.permutation(inds.pheno)
        pvals = plink.assoc(inds, genos)
        yield get_landscape(pvals)


def landscape_generator_file():
    inp = open(settings.infile) if not settings.input_format == 'gzip' else gzip.open(settings.infile)
    next(inp)
    for line in inp:
        yield get_landscape([float(value) for value in line.split()])
    inp.close()


def write_results(landscape, pvalues, ntests, settings):

    columns = ['s', 'i', 'Y', 'type', 'p_value', 'adjusted_p_value']

    values = [[str(segment.s),
               str(segment.e),
               str(segment.Y),
               'independent' if segment.independent else 'dependent',
               str(round(min(pvalues[k], 1), 6)),
               str(min(pvalues[k] * ntests, 1))]
              for k, segment in enumerate(landscape.maximal_segments)]

    df = pd.DataFrame(values, columns=columns)
    path = settings.outfile + '_' + settings.score_evaluation + '.txt'
    df.to_csv(path, index=False, sep='\t')


def write_ak_text(snp_p, landscape, settings):
    columns = ['position', 'p_value', 'A_k']
    values = [[str(i),
               snp_p[i],
               landscape.A[i]]
              for i in range(len(snp_p))]

    df = pd.DataFrame(values, columns=columns)
    path = settings.outfile + '_A_k.txt'
    df.to_csv(path, index=False, sep='\t')


def get_pvalues_A0(landscape):
    if settings.lattice:
        lmbda = math.log((1 - settings.gamma) / settings.gamma)
        return [(1 - pow(math.e, -lmbda))*pow(math.e, -lmbda*segment.Y) for segment in landscape.maximal_segments]
    else:
        from scipy.optimize import minimize_scalar
        lmbda = minimize_scalar(lambda x: abs(log(1 - x) / x - log(settings.gamma)),
                                bounds=(0.001, 0.9999), method='bounded').x
        return [pow(math.e, -lmbda*segment.Y) for segment in landscape.maximal_segments]


def get_pvalues_A1(orig_landscape, landscape_generator):

    def add_hits(landscape):
        for orig_segment in orig_landscape.maximal_segments:
            hits_in_permutation = sum(1 for segment in landscape.maximal_segments if segment.Y >= orig_segment.Y)
            hits_per_segment = hits_in_permutation / len(landscape.maximal_segments) if landscape.maximal_segments else 0
            hits[orig_segment.s].append(hits_per_segment)

    hits = {segment.s: [0] for segment in orig_landscape.maximal_segments}

    M = [len(orig_landscape.maximal_segments)]
    add_hits(orig_landscape)

    for _ in range(settings.permutations):
        new_landscape = next(landscape_generator)
        M.append(len(new_landscape.maximal_segments))
        add_hits(new_landscape)

    pvalues = [sum(hits[segment.s]) / (settings.permutations + 1) for i, segment in enumerate(orig_landscape.maximal_segments)]
    EM = sum(M) / len(M)

    return pvalues, EM


def get_pvalues_A2(landscape, landscape_generator):

    tops = [k for k in range(len(landscape.position_Y)) if landscape.position_Y[k]]
    hits = [0 for _ in range(len(landscape.A))]
    perms = [0 for _ in range(len(landscape.A))]
    M = [len(landscape.maximal_segments)]

    for _ in range(settings.permutations):
        new_landscape = next(landscape_generator)
        if new_landscape.maximal_segments:
            M.append(len(new_landscape.maximal_segments))
            for k in tops:
                if new_landscape.position_Y[k]:
                    perms[k] += 1
                    if landscape.position_Y[k] <= new_landscape.position_Y[k]:
                        hits[k] += 1
    if any(perms[i] < 1 for segment in landscape.maximal_segments for i in range(segment.s, segment.e + 1)):
        print('More permutations are required for approach 2')
        assert not any(perms[i] < 1 for i in range(len(landscape.A)))

    pvalues = [max([hits[i] / perms[i] for i in range(segment.s, segment.e + 1)]) for segment in landscape.maximal_segments]
    EM = sum(M) / len(M)

    return pvalues, EM


def get_landscape(values):
    if settings.transformation == 'dichotomous':
        Z = [1 if value < settings.gamma else -1 for value in values]
    elif settings.transformation == 'log':
        Z = [log(settings.gamma / value) for value in values]
    else:
        Z = values

    A = [0]
    for n in range(len(Z)):
        A.append(max(0, A[n] + Z[n]))
    A.append(0)

    defined_values = range(1, len(A) - 1)
    s = [[n] for n in defined_values if A[n] > 0 and A[n - 1] == 0]
    t = [n for n in defined_values if A[n + 1] == 0 and A[n] > 0]

    section_indices = range(len(s))
    sections = [range(s[i][0], t[i] + 1) for i in section_indices]

    Y = [[max(A[n] for n in section)] for section in sections]
    e = [[min(n for n in sections[i] if A[n] == Y[i][0])] for i in section_indices]

    for i in section_indices:
        j = 0
        while 1:
            section_range = range(e[i][j] + 1, t[i] + 1)
            possible_starts = [n for n in section_range if A[n] > A[n - 1]]
            if not possible_starts:
                break
            start = min(possible_starts)
            possible_ends = [n for n in range(start + 1, t[i] + 2) if A[n] < A[start]]
            end = min(possible_ends)
            section = range(start, end + 1)
            if not section:
                break
            j += 1
            s[i].append(start)
            Y[i].append(max(A[n] for n in section) - A[s[i][j] -1])
            e[i].append(min(n for n in section if abs(A[n] - (Y[i][j] + A[s[i][j] -1])) < 0.00001))

    A = A[1:]
    position_Y = [0 for _ in range(len(A) + 1)]
    maximal_segments = []
    for i in section_indices:
        for j in range(len(s[i])):
            maximal_segments.append(MaximalSegment(s[i][j], e[i][j], sections[i][0], sections[i][-1], Y[i][j], j == 0))
            for k in range(s[i][j], e[i][j] + 1):
                position_Y[k] = Y[i][j]

    return Landscape(A, s, e, Y, position_Y, sections, section_indices, maximal_segments)


def run():

    if settings.input_format != 'plink':
        inp = open(settings.infile) if not settings.input_format == 'gzip' else gzip.open(settings.infile)
        line = inp.readline()
        orig_values = [float(value) for value in line.split()]
        if not settings.permutations:
            settings.permutations = sum(1 for _ in inp)
        inp.close()

        landscape = get_landscape(orig_values)
        landscape_generator = landscape_generator_file()
    else:
        inds, snps, genos = plink.read_plink(settings.infile)
        orig_values = plink.assoc(inds, genos)
        landscape = get_landscape(orig_values)
        landscape_generator = landscape_generator_plink()

    if settings.permutations > 0:
        if settings.score_evaluation == 'A1':
            pvalues, EM = get_pvalues_A1(landscape, landscape_generator)
        elif settings.score_evaluation == 'A2':
            pvalues, EM = get_pvalues_A2(landscape, landscape_generator)
        else:
            pvalues = get_pvalues_A0(landscape)
            EM = len(pvalues)
    else:
        pvalues = [1 for _ in landscape.maximal_segments]
        EM = 1

    write_results(landscape, pvalues, EM, settings)
    write_ak_text(orig_values, landscape, settings)


def get_settings(args):
    parser = argparse.ArgumentParser(prefix_chars='--', formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    parser.add_argument('--help', action='help', help='\nShow this help message and exit\n\n')
    parser.add_argument('--infile', help='\nThe input file path\n\n', required=True)
    parser.add_argument('--outfile', help='\nThe output file path\n\n', required=True)
    parser.add_argument('--input_format', choices=['text', 'gzip', 'plink'], default='text', help='Default: %(default)s\n\n')
    parser.add_argument('--transformation', choices=['dichotomous', 'log', 'none'], default='log', help='Transformation of the input data. dichotomous applies the 1/-1 transformation shown example 1, and log applies the log transformation shown in example 2.\nDefault: %(default)s\n\n')
    parser.add_argument('--score_evaluation', help='How to evaluate the scores. (See section 4.2-4.4).\n\n', choices=['A0', 'A1', 'A2'], required=True)
    parser.add_argument('--gamma', type=float, default=0.05, help='\nThe threshold of significance. (See Example 1 and 2).\nDefault: %(default)s\n\n')
    parser.add_argument('--no-lattice', dest='lattice', action='store_false', default=True, help= 'Set if Z is not a lattice variable. (See definition 8 and 9)\n\n')
    parser.add_argument('--permutations', type=int, help='The number of permutations that are performed if score_evaluation is A1 or A2.\nDefault: If input_format is text or gzip, the number of lines in the input file.\n\n')
    parser.add_argument('--random_seed', default=1, help='Random seed for permutations if input_format is plink.\nDefault: %(default)s\n\n')
    parser.add_argument('--no-draw', dest='draw', action='store_false', help='Do not draw the landscape')

    return parser.parse_args(args)


def draw():
    import matplotlib.pyplot as plt
    ak = pd.read_csv(settings.outfile + '_A_k.txt', sep='\s+')
    ak['logp'] = ak.p_value.apply(lambda x: -math.log10(x))
    plt.plot(ak.position, ak.logp)
    sig = -math.log10(0.05 / len(ak))
    plt.plot([min(ak.position), max(ak.position)], [sig, sig], color='0.75', linestyle='--')
    plt.ylabel('-log10(p-value)')
    plt.xlabel('position')

    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    fig.savefig(settings.outfile + '.png', dpi=100)


if __name__ == "__main__":

    settings = get_settings(sys.argv[1:])
    run()

    if settings.draw:
        draw()

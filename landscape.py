from __future__ import division
import random
import math
import gzip
from math import log
from math import pow
import sys
import argparse


class Landscape:
    def __init__(self, A, s, e, Y, position_Y, sections, section_indices, maximal_segments):
        self.A = A
        self.s = s
        self.e = e
        self.Y = Y
        self.position_Y = position_Y
        self.sections = sections
        self.section_indices = section_indices
        self.maximal_segments = maximal_segments


class MaximalSegment:
    def __init__(self, s, e, sec_s, sec_e, Y, independent):
        self.s = s
        self.e = e
        self.sec_s = sec_s
        self.sec_e = sec_e
        self.Y = Y
        self.independent = independent


def landscape_generator_plink(inds, snps, settings):
    def shuffle_affection(inds):
        aff = [ind.aff for ind in inds]
        random.shuffle(aff)
        for i in range(len(inds)):
            inds[i].aff = aff[i]
        return inds

    random.seed(settings["random_seed"])
    while 1:
        inds = shuffle_affection(inds)
        f.calc_p(inds, snps)
        yield get_landscape_from_plink(snps, settings)


def get_landscape_from_plink(snps, settings):
    assert settings["transformation"] in ["dichotomous", "log"]
    if settings["transformation"] == "dichotomous":
        return get_landscape([1 if snp.p < settings["gamma"] else -1 for snp in snps])
    elif settings["transformation"] == "log":
        return get_landscape([log(settings["gamma"] / snp.p) for snp in snps])


def landscape_generator_file(settings):
    if settings["input_type"] == "gzip":
        with gzip.open(settings["input"] + ".gz") as input:
            next(input)
            for line in input:
                line = line.split()
                yield get_landscape_from_file(settings, line)
    else:
        with open(settings["input"] + ".txt") as input:
            next(input)
            for line in input:
                line = line.split()
                yield get_landscape_from_file(settings, line)


def get_landscape_from_file(settings, line):
    if settings["transformation"] == "dichotomous":
        return get_landscape([1 if float(value) < settings["gamma"] else -1 for value in line])
    elif settings["transformation"] == "log":
        return get_landscape([log(settings["gamma"] / float(elem)) for elem in line])
    elif settings["transformation"] == "none":
        return get_landscape([float(elem) for elem in line])


def write_results(landscape, pvalues, ntests, settings):

    titles = ["s", "e", "Y", "type", "p_value", "adjusted_p_value"]

    values = [[str(segment.s), str(segment.e), str(segment.Y), "independent" if segment.independent else "dependent",
              str(round(min(pvalues[k], 1), 6)), str(min(pvalues[k] * ntests, 1))]
              for k, segment in enumerate(landscape.maximal_segments)]

    f.write(settings["output"] + "_" + settings["score_evaluation"] + ".txt", values, titles, settings["format_columns"])


def write_ak_plink(snps, snp_p, landscape, settings):
        titles = ["chr", "bp", "p_value", "A_k"]
        values = [[snps[i].chr, snps[i].pos, snp_p[i], landscape.A[i]] for i in range(len(snps))]
        f.write(settings["output"] + "_A_k.txt", values, titles, settings["format_columns"])


def write_ak_text(snp_p, landscape, settings):
        titles = ["position", "p_value", "A_k"]
        values = [[str(i), snp_p[i], landscape.A[i]] for i in range(len(snp_p))]
        f.write(settings["output"] + "_A_k.txt", values, titles, settings["format_columns"])


def get_pvalues_approach1(landscape, landscape_generator, settings):

    def get_hits_per_segment(landscape, y):
        hits_in_permutation = sum(1 for segment in landscape.maximal_segments if segment.Y >= y)
        if landscape.maximal_segments:
            return hits_in_permutation / len(landscape.maximal_segments)
        else:
            return 0

    M = [len(landscape.maximal_segments)]
    for segment in landscape.maximal_segments:
        segment.hits = [get_hits_per_segment(landscape, segment.Y)]

    for _ in range(settings["permutations"]):
        new_landscape = landscape_generator.next()
        M.append(len(new_landscape.maximal_segments))
        for segment in landscape.maximal_segments:
            segment.hits.append(get_hits_per_segment(new_landscape, segment.Y))

    pvalues = []
    for segment in landscape.maximal_segments:
        pvalues.append(sum(segment.hits) / (settings["permutations"] + 1))

    EM = sum(M) / len(M)

    return pvalues, EM


def get_pvalues_approach2(landscape, landscape_generator, settings):

    tops = [k for k in range(len(landscape.position_Y)) if landscape.position_Y[k]]
    hits = [0 for _ in range(len(landscape.A))]
    perms = [0 for _ in range(len(landscape.A))]
    M = [len(landscape.maximal_segments)]

    for _ in range(settings["permutations"]):
        new_landscape = landscape_generator.next()
        if new_landscape.maximal_segments:
            M.append(len(new_landscape.maximal_segments))
            for k in tops:
                if new_landscape.position_Y[k]:
                    perms[k] += 1
                    if landscape.position_Y[k] <= new_landscape.position_Y[k]:
                        hits[k] += 1
    if any(perms[i] < 1 for segment in landscape.maximal_segments for i in range(segment.s, segment.e + 1)):
        print "More permutations are required for approach 2"
        assert not any(perms[i] < 1 for i in range(len(landscape.A)))

    pvalues = []
    for segment in landscape.maximal_segments:
        pvalues_segment = [hits[i] / perms[i] for i in range(segment.s, segment.e + 1)]
        pvalues.append(max(pvalues_segment))

    EM = sum(M) / len(M)

    return pvalues, EM


def get_pvalues_independent(landscape, settings):
    pvalues = []
    if settings["lattice"]:
        lmbda = math.log((1 - settings["gamma"]) / settings["gamma"])
        for segment in landscape.maximal_segments:
            pvalue = (1 - pow(math.e, - lmbda)) * pow(math.e, -lmbda * (segment.Y - settings["minimum_score"]))
            pvalues.append(pvalue)
    else:
        from scipy.optimize import minimize_scalar

        def f(x):
            return abs(log(1 - x) / x - log(settings["gamma"]))

        lmbda = minimize_scalar(f, bounds=(0.001, 0.9999), method='bounded').x
        for segment in landscape.maximal_segments:
            pvalue = lmbda * pow(math.e, -lmbda * (segment.Y - settings["minimum_score"]))
            pvalues.append(pvalue)

    return pvalues


def get_landscape(Z):
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


def run(settings):
    if settings["input_type"] == "plink":
        inds, snps = f.read_plink(settings["input"])
        f.calc_p(inds, snps)
        orig_snp_p = [snp.p for snp in snps]
        landscape = get_landscape_from_plink(snps, settings)
        landscape_generator = landscape_generator_plink(inds, snps, settings)
    else:
        if settings["input_type"] == "gzip":
            with gzip.open(settings["input"] + ".gz") as input:
                line = input.readline()
                orig_snp_p = line.split()
        elif settings["input_type"] == "text":
            with open(settings["input"] + ".txt") as input:
                line = input.readline()
                orig_snp_p = line.split()
        landscape = get_landscape_from_file(settings, orig_snp_p)
        landscape_generator = landscape_generator_file(settings)

    if settings["permutations"] > 0:
        if settings["score_evaluation"] == "approach1":
            pvalues, EM = get_pvalues_approach1(landscape, landscape_generator, settings)

        if settings["score_evaluation"] == "approach2":
            pvalues, EM = get_pvalues_approach2(landscape, landscape_generator, settings)

        if settings["score_evaluation"] == "independent":
            pvalues = get_pvalues_independent(landscape, settings)
            EM = len(pvalues)
    else:
        pvalues = [1 for _ in landscape.maximal_segments]
        EM = 1

    write_results(landscape, pvalues, EM, settings)

    if settings["input_type"] == "plink":
        write_ak_plink(snps, orig_snp_p, landscape, settings)
    else:
        write_ak_text(orig_snp_p, landscape, settings)


def get_settings(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, prefix_chars='--')
    parser.add_argument('--infile', help="The input file path")
    parser.add_argument('--outfile', help="The output file path")
    parser.add_argument('--minimum_score', help="The minimum score if independent is chosen for multiple_testing_adjustment. (y_0 in formula 3.8 / 3.9)")
    parser.add_argument('--score_evaluation', help="How to evaluate the scores. Approaches explained in section 3.1 and section 3.2.")
    parser.add_argument('--permutations', help="The number of permutations that are performed if approach1 or approach2 is chosen for multiple_testing_adjustment. (Note, if input type 'text' is used, the program still expects the number of permutations to be set. It can be set to the full number of permutations in the file, or to a lower number, but not to a higher number.)")
    parser.add_argument('--gamma', help="The threshold of significance. (See Example 3.1 and 3.2).")
    parser.add_argument('--lattice', help="Set to true if Z is a lattice variable. (See example 3.1)")
    parser.add_argument('--random_seed', default=1)

    return parser.parse_args(args)


if __name__ == "__main__":

    settings = get_settings(sys.argv[1:])

    run(settings)

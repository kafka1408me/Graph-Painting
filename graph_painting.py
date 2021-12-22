import random
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import math

import graphs
import sys


count_nodes = 30
count_edges = 90 #count_nodes*(count_nodes-1)//2

k = count_nodes - 1
num_node = 0
list_adjacency = []
while k != 0:
    list_adjacency += [x for x in range(num_node*count_nodes + num_node + 1, count_nodes*(num_node + 1))]
#    print(num_node*(count_nodes + 1) + num_node + 1, count_nodes*(num_node + 1))
    num_node += 1
    k -= 1
    pass

random.shuffle(list_adjacency)
#print(list_adjacency)

list_edges = []
for x in list_adjacency[:count_edges]:
    list_edges += [(x // count_nodes, x % count_nodes)]
    pass

list_other_nodes = [x // count_nodes for x in list_adjacency[count_edges:]]

g = nx.from_edgelist(list_edges)
g.add_nodes_from(list_other_nodes)

# problem constants:
HARD_CONSTRAINT_PENALTY = 10  # the penalty factor for a hard-constraint violation

# Genetic Algorithm constants:
POPULATION_SIZE = 100
P_CROSSOVER = 0.9  # probability for crossover
P_MUTATION = 0.1   # probability for mutating an individual
MAX_GENERATIONS = 100

MAX_COLORS = 20

#g = nx.petersen_graph()

gcp = graphs.GraphColoringProblem(graph=g, hardConstraintPenalty=HARD_CONSTRAINT_PENALTY)






def gaAlgorithm(population_size, max_generations, crossover_percent, mutation_percent, max_colors, graph):
    gcp = graphs.GraphColoringProblem(graph=graph, hardConstraintPenalty=HARD_CONSTRAINT_PENALTY)
    count_chromosome = len(gcp)
    population = createPopulation(population_size, count_chromosome, max_colors)
    for i in range(max_generations):
        print('Generation', i + 1)
        # турнирный отбор
        selected_indxs = []
        for j in range(population_size*2):
            indx1 = random.randint(0, population_size-1)
            indx2 = random.randint(0, population_size-1)
            selected_indxs.append(indx1 if gcp.getCost(population[indx1]) < gcp.getCost(population[indx2]) else indx2)
            pass
        # скрещивание
        new_population = []
        for j in range(0, population_size, 2):
            parent1_indx = selected_indxs[j]
            parent2_indx = selected_indxs[j+1]
            child_1 = None
            child_2 = None
            if parent1_indx != parent2_indx and random.randint(0,100) >=crossover_percent:
                point_split = random.randint(1, count_nodes-1)
                child_1 = population[parent1_indx][:point_split] + population[parent2_indx][point_split:]
                child_2 = population[parent2_indx][:point_split] + population[parent1_indx][point_split:]
                pass
            else:
                child_1 = population[parent1_indx].copy()
                child_2 = population[parent2_indx].copy()
                pass
            new_population.append(child_1)
            new_population.append(child_2)
            pass
        population = new_population
        # mutation
        for chromosome in new_population:
            mutation(chromosome, mutation_percent, max_colors)
            pass
        pass

    min_cost = math.inf
    best = None
    for chromosome in population:
        cost = gcp.getCost(chromosome)
        if cost < min_cost:
            best = chromosome
            pass
        pass

    if gcp.getViolationsCount(best) == 0:
        return best
    return None

# best = gaAlgorithm(population_size=POPULATION_SIZE, max_generations=MAX_GENERATIONS, count_chromosome=count_nodes,
#             crossover_percent=P_CROSSOVER*100, mutation_percent=P_MUTATION*100, max_colors=MAX_COLORS)
#
#
# plot = gcp.plotGraph(best)
# plot.show()

def createIndividual(count_genes, max_colors):
    max_color_indx = max_colors - 1
    return [random.randint(0, max_color_indx) for _ in range(count_genes)]

def createPopulation(count_chromosomes, count_genes, max_colors):
    new_population = []
    for i in range(count_chromosomes):
        chromosome = []
        for j in range(count_genes):
            chromosome.append(random.randint(0, max_colors - 1))
            pass
        new_population.append(chromosome)
        pass
    return new_population


def immuneAlgorithm(population_size, max_generations, max_colors, p_stay_percent, updated_percent, graph):
    gcp = graphs.GraphColoringProblem(graph=graph, hardConstraintPenalty=HARD_CONSTRAINT_PENALTY)
    count_genes = len(gcp)
    population = createPopulation(population_size, count_genes, max_colors)
    count_stay_p = int(p_stay_percent*population_size)
    count_updated = int(updated_percent*population_size)
    count_deleted_p = population_size - count_stay_p

    average_clones_for_one = population_size / count_stay_p
    p = 1.0 / count_genes
    population_indxs = [i for i in range(population_size)]
    list_indxs = [i for i in range(count_genes)]

    utilities = [(i, gcp.getCost(population[i])) for i in range(population_size)]
    utilities = sorted(utilities, key=lambda el: el[1])[:count_stay_p]

    for generation in range(max_generations):
        print("Generation", generation)
        sum_utilities = 0
        for el in utilities:
            sum_utilities += el[1]

        new_population = []
        average_utility = sum_utilities / count_stay_p
        i = 0
        count_need_cloning = count_deleted_p
        utilities_2 = []
        while i < count_stay_p:
            count_clones = 0
            utility = utilities[i][1]
            if count_need_cloning > 0:
                utility_relative = (average_utility / utility) * average_clones_for_one
                count_clones = round(utility_relative)
                count_need_cloning -= count_clones

                if count_need_cloning < 0:
                    count_clones += count_need_cloning

            count_clones += 1
            utilities_2 += [(count_clones, utility)]

            new_population += [population[i] for _ in range(count_clones)]
            i += 1
            pass
        print('new_population size =', len(new_population))

        # Мутация
        i = 0
        for count_antibodies, utility in utilities_2:
            mut = math.exp(p * utility)
            for k in range(count_antibodies):
                count_changed_genes = math.ceil(mut)
                if count_changed_genes > count_genes:
                    count_changed_genes = count_genes
                random.shuffle(list_indxs)
                for t in range(count_changed_genes):
                    population[i][list_indxs[t]] = random.randint(0, max_colors - 1)
                    pass
                i += 1
                pass
            pass
        population = new_population

        utilities = [(i, gcp.getCost(population[i])) for i in range(population_size)]
        utilities = sorted(utilities, key=lambda el: el[1])[:count_stay_p]

        best = population[utilities[0][0]]
        if gcp.getViolationsCount(best) == 0:
            return best

        if count_updated != 0:
            updated_indxs = random.sample(population_indxs, count_updated)
            for indx in updated_indxs:
                population[indx] = createIndividual(count_genes, max_colors)


    return None


#best = immuneAlgorithm(POPULATION_SIZE, MAX_GENERATIONS, MAX_COLORS, 0.1, 0.1, g)
#
# if best is not None:
#     gcp = graphs.GraphColoringProblem(graph=g, hardConstraintPenalty=HARD_CONSTRAINT_PENALTY)
#     plot = gcp.plotGraph(best)
#     plot.show()






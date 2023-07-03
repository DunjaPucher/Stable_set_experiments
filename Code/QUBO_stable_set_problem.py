# --------------------------------------------------------------------------
#   Copyright (C) 2023 Dunja Pucher <dunja.pucher@aau.at>
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see https://www.gnu.org/licenses/.
# ==========================================================================

# This program solves the stable set problem by using D-Wave's QPU or hybrid 
# solver
# For QPU solver: use function compute_with_QPU_dwave
# For hybrid solver: use function compute_with_hybrid_dwave
# In order to perform computations, you need a token for D-Wave
# enter you token here:
    
token = 'secret'

# ==========================================================================

from dwave.system import LeapHybridSampler
import dimod
from dwave.system import DWaveSampler, EmbeddingComposite
import numpy as np

##################################################

def create_matrix_Q(Q, a):
    if a == 1:
        return Q
    if a > 1:
        n = Q.shape[1]
        for i in range(n):
            for j in range(n):
                if Q[i, j] == 1:
                    Q[i, j] = a
        return Q
            
##################################################

def compute_with_QPU_dwave(instance, a):
    
    # INPUT: instance, written as string (e.g. "Paley61"), value of the 
    # parameter beta denoted here as a
    # OUTPUT: output is written to files
    
    # prepare input for -DWave
    Q = np.loadtxt(instance+'_stable_set_matrix_Q.txt', dtype='i', delimiter=',')
    Q = create_matrix_Q(Q, a)
    bqm = dimod.BQM.from_qubo(Q)
    
    # computation with the QPU solver
    sampleset_QPU = EmbeddingComposite(DWaveSampler(token = token)).sample(bqm, num_reads = 1000)

    # saving files

    # all results
    R = sampleset_QPU.record
    f = open(instance+'_QPU_all_results_a_'+str(a)+'.txt','w')
    f.write('{}'.format(R)) 
    f.close()
    
    C = []
    E = []
    count_stable_sets = 0
    solution_count = len(R)
    for i in range(solution_count):
        c = sum(R[i][0])
        C.append(c)
        e = R[i][1]
        E.append(e)
        if c == - e:
            count_stable_sets = count_stable_sets + 1
    
    f = open(instance+'_QPU_count_solutions_a_'+str(a)+'.txt','w')
    f.write('{}'.format(solution_count)) 
    f.close()

    f = open(instance+'_QPU_count_stable_set_solutions_a_'+str(a)+'.txt','w')
    f.write('{}'.format(count_stable_sets)) 
    f.close()

    # best result
    B = sampleset_QPU.first
    f = open(instance+'_QPU_best_result_a_'+str(a)+'.txt','w')
    f.write('{}'.format(B)) 
    f.close()

    # best solution -- value and vector
    sol = B[1]
    f = open(instance+'_QPU_best_solution_energy_a_'+str(a)+'.txt','w')
    f.write('{}'.format(str(sol)))
    f.close()

    v = B[0]
    x = [0]*len(v)
    for i in range(len(v)):
        x[i] = v[i]

    f = open(instance+'_QPU_solution_vector_a_'+str(a)+'.txt','w')
    f.write('{}'.format(str(x)))
    f.close()
    
    c = sum(x)
    f = open(instance+'_QPU_solution_cardinality_a_'+str(a)+'.txt','w')
    f.write('{}'.format(str(c)))
    f.close()
    
    A = np.loadtxt(instance+'_stable_set_A_dense.txt', dtype='i', delimiter=',')
    count = 0
    for k in range(len(A)): # B[0]
        for l in range(len(A)): # B[0]
            if A[k, l] == 1:
                if x[k] + x[l] == 2:
                    count = count + 1
    f = open(instance+'_QPU_solution_type_a_'+str(a)+'.txt','w')
    if count == 0:
        f.write('{}'.format("This is a stable set"))
    else:
        f.write('{}'.format("This is not a stable set"))
    f.close()

    # time log 
    T = sampleset_QPU.info
    f = open(instance+'_QPU_dwave_times_a_'+str(a)+'.txt','w')
    f.write('{}'.format(T)) 
    f.close()
    
    return("Done")

##################################################

def compute_with_hybrid_dwave(instance, a):
    
    # INPUT: instance, written as string (e.g. "Paley61"), value of the 
    # parameter beta denoted here as a
    # OUTPUT: output is written to files
    
    # prepare input for D-Wave
    Q = np.loadtxt(instance+'_stable_set_matrix_Q.txt', dtype='i', delimiter=',')
    Q = create_matrix_Q(Q, a)
    bqm = dimod.BQM.from_qubo(Q)

    # computation with the hybrid solver
    sampleset_hybrid = LeapHybridSampler(token = token).sample(bqm)

    # saving files

    # result
    B = sampleset_hybrid.first
    f = open(instance+'_hybrid_result_a_'+str(a)+'.txt','w')
    f.write('{}'.format(B)) 
    f.close()

    # solution -- value and vector
    sol = B[1]
    f = open(instance+'_hybrid_solution_energy_a_'+str(a)+'.txt','w')
    f.write('{}'.format(str(sol)))
    f.close()

    v = B[0]
    x = [0]*len(v)
    for i in range(len(v)):
        x[i] = v[i]

    f = open(instance+'_hybrid_solution_vector_a_'+str(a)+'.txt','w')
    f.write('{}'.format(str(x)))
    f.close()
    
    c = sum(x)
    f = open(instance+'_hybrid_solution_cardinality_a_'+str(a)+'.txt','w')
    f.write('{}'.format(str(c)))
    f.close()
    
    A = np.loadtxt(instance+'_stable_set_A_dense.txt', dtype='i', delimiter=',')
    count = 0
    for k in range(len(A)): # B[0]
        for l in range(len(A)): # B[0]
            if A[k, l] == 1:
                if x[k] + x[l] == 2:
                    count = count + 1
    f = open(instance+'_hybrid_solution_type_a_'+str(a)+'.txt','w')
    if count == 0:
        f.write('{}'.format("This is a stable set"))
    else:
        f.write('{}'.format("This is not a stable set"))
    f.close()

    # time log 
    T = sampleset_hybrid.info
    f = open(instance+'_hybrid_dwave_times_a_'+str(a)+'.txt','w')
    f.write('{}'.format(T)) 
    f.close()
    
    return("Done")
   
##################################################


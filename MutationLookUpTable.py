#!/usr/bin/env python3

import math
import numpy as np

#*#*#*#*#*#*#*#*#*#*#*#*#*#*# NOTES #*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#                                                               #
# 1. 	DeltaGs from file are assumed to be +ve for             #
#       DeltaG_inactive and -ve for DeltaG_active.              #
#                                                               #
# 2.    IC50s from two sources A & B.                           #
#       WT, G250E, E255K, E255V, T315I are from Source A        #
#       T315M & Y253H-E255V are from Source B (rest of the      #
#       mutations aren't used in the nuclear approch            #
#       simulation                                              #
#                                                               #
#       A:  Three novel patient-derived BCR/ABL mutants show    #
#           different sensitivity to second and third           #
#           generation tyrosine kinase inhibitors               #
#           - Redaelli et al                                    #
#                                                               #
#       B:  Supplementary material to: BCR-ABL1 Compound        #
#           Mutations Combining Key Kinase Domain Positions     #
#           Confer Clinical Resistance to Ponatinib in Ph       #
#           Chromosome-Positive Leukemia                        #
#           - Zabriskie et al                                   #
#                                                               #
# 3.    Order and numbering of mutations:                       #
#       1.  Wild-type   2.  G250E       3. Y253H                #
#       4.  E255K       5.  E255V       6. T315I                #
#       7.  T315M       8.  G250E/T315I 9. Y253H/E255V          #
#       10. Y253H/T315I 11. E255V/T315I                         #
#							                                    #
# 4.    Drug numbering: 0. Ponatinib, 1. Imatinib, 2. Dasatinib #
#							                                    #
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#



#Vmax (pmol/min),kcatAct(/min),KM(microM),DeltaG(kcal/mol),Ponatinib-IC50(nanoM),Imatinib-IC50(nanoM),Dasatinib-IC50(nanoM), DelGErrMax
Mutant = np.array([[1.9,  66.0, 17.0,   0.0,    2.1,  527.0,  1.8, 2.6], \
                   [5.1, 175.2, 14.3,   -2.2,   12.5, 3613.0,  8.1, 3.3], \
                   [0.8,  26.9,  4.6, -24.0,   29.8,  17700,  5.9, 2.9], \
                   [1.8,  63.6, 15.6,   -0.5,   17.6,   3174, 10.3, 2.5],\
                   [0.2,   6.8, 22.1,   -0.4,   27.2,   8953,  6.3, 2.4], \
                   [0.4,  12.2,  7.2,   0.1,    6.3,   9221,137.3, 2.3], \
                   [0.2,   8.1,  1.9,   0.1,  577.5,  10240,  768, 4.6],\
                   [0.9,  30.0,  1.0,  19.1,  152.4,  10240,  768, 2.5], \
                   [0.4,  13.2,  0.7,   0.4,  203.5,  10240, 18.1, 2.5], \
                   [0.1,   3.8,  0.5, -29.3,  357.9,  10240,  768, 2.5], \
                   [0.1,   4.1,  0.6,  23.7,  659.5,  10240,  768, 2.5]],float)

DrugLabels = ["Ponatinib", "Imatinib", "Dasatinib"]
MutantLabels = ["Wild type", "G250E", "Y253H", "E255K", "E255V", "T315I", "T315M", "G250E-T315I", "Y253H-E255V", "Y253H-T315I", "E255V-T315I"]

#Halflife(hours), PeakTime(hours), AverageConc(nM), DayCountSteadyStateBegins
DrugDynamicTimings = np.array([ [24, 6,  200, 7],\
                                [18, 3, 2000, 5],\
                                [ 4, 1,  500, 1]],float)

#F_Bioavailibility, Dose([10**(-5)]mol), VolumeOfDistribution(L), EliminationRate([10**(-2)/hour), AbsorptionRate(/hour)
DrugDynamicsConstants = np.array([  [1, 8.449, 1223, 2.888, 1.302],\
                                    [1, 81.04,  435, 3.851, 0.940],\
                                    [1, 36.88, 2502, 17.32, 1.740]],float)


#kATrans, koffSubstrate, Ponatinib koffInhibitorInactive, Imatinib koffInhibitorInactive, Dasatinib koffInhibitorActive (all in /min)
Rates = np.array([[60,  33, 0.00488, 0.059, 0.00233],\
                  [55, 350, 0.00488, 0.059, 0.01010],\
                  [60,  75, 0.00488, 0.059, 1],\
                  [60, 127, 0.00488, 0.059, 0.01310],\
                  [60, 122, 0.00488, 0.059, 0.00861],\
                  [60,  29.7, 0.00488, 0.059, 0.136],\
                  [60,  16.2, 0.00488, 0.059, 0.354],\
                  [60,  75, 0.00488, 0.059, 1],\
                  [60,  26.4, 0.00488, 0.059, 0.00372],\
                  [60,  75, 0.00488, 0.059, 1],\
                  [60,  75, 0.00488, 0.059, 1]],float)
    


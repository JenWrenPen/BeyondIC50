#!/usr/bin/env python3

import numpy as np
import math
import random
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import stats
import sympy as sy

from MutationLookUpTable import * # imports values for Vmax, kCatAct, KM, IC50, DeltaG and Names

#*#*#*#*#*#*#*#*#*#*#*#*#*#*# NOTES #*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#                                                               #
# 1.    The timestep should be about 1/100th of a second        #
#                                                               #
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

################ Constant parameters ################

Temp = 310 # K
kB = 3.29982916 * 10**(-27) #in kcal/K
Beta = 1/(kB*Temp) # per kcal

SubConc = 10 * 10 ** (-6) # M
Etotal = 1*10**(-6) # M

WildTypeDeltaG = -1 # kcal/mol

MutationNumbersWanted = [0,1,3,4,5,6,8]
DrugsWanted = [0,1,2]

DrugName = ["Ponatinib","Imatinib","Dasatinib"]
MutationNames = ["wild-type", "G250E","E255K", "E255V", "T315I", "T315M", "Y253H-E255V"]

DayTotal = 10 # in days
TimeStep = 0.001 # in seconds
StepTotal = int(DayTotal * 24 * 3600 / TimeStep)
TimeTotal = DayTotal*24*3600

################ Some initial values ################

RateCatalysis = 0
RateTransActive = 0
RateOffSubstrate = 0
RateTransInactive = 0
RateOnSubstrate = 0
RateOnInhibitor = 0

PrevActUnb = 0
PrevActBou = 0
PrevActInh = 0
PrevInaUnb = 0
PrevInaInh = 0

PrevInhConc = 0

SwitchDay = 0
Gamma = 0
Epsilon = 0
Alpha = 0

ThisDrug = 0


################ Functions ################

#*# Convert per minute to per time step #*#

def convertRateToPerSec(RATE):
    return RATE/60

#*# Constants for concentration cycle equations #*#

def calculateGamma(FbIOAVAILABILITY,dOSE,aBSORPTIONrATE,eLIMINATIONrATE,vOLUMEoFdISTRIBUTION):
    return (FbIOAVAILABILITY * dOSE * aBSORPTIONrATE) / (vOLUMEoFdISTRIBUTION * (aBSORPTIONrATE - eLIMINATIONrATE))

def calculateEpsilon(eLIMINATIONrATE):
    return math.exp(eLIMINATIONrATE*24)

def calculateAlpha(aBSORPTIONrATE):
    return math.exp(aBSORPTIONrATE*24)

#*# Concentration cycle equations #*#

def findInhibitorConcentration(tIME,sWITCHdAY,gAMMA,ePSILON,aLPHA,aBSORPTIONrATE,eLIMINATIONrATE):

    DayCounter = math.floor(tIME/(3600*24))
    TimeInDay = (tIME / 3600) - (DayCounter * 24) #in hours
    TimeInHours = tIME / 3600

    if (DayCounter < sWITCHdAY):
        
        SumEpsilon = 0
        SumAlpha = 0

        for Day in range(DayCounter+1):

            SumEpsilon = SumEpsilon + (ePSILON ** Day)  # no units
            SumAlpha = SumAlpha + (aLPHA ** Day)        # no units

        return gAMMA * ((SumEpsilon * math.exp(-eLIMINATIONrATE * TimeInHours)) - (SumAlpha * math.exp(-aBSORPTIONrATE * TimeInHours)))

    else:
        return gAMMA * (((ePSILON * math.exp(-eLIMINATIONrATE * TimeInDay))/(ePSILON - 1)) - ((aLPHA * math.exp(-aBSORPTIONrATE * TimeInDay))/(aLPHA - 1)))

#*# Concentration equations #*#

#Decay#

def functionDecay(rate,time):
    return 1-rate*time

#Growth#

def functionGrowth(rate,time):
    return 1-functionDecay(rate,time)

#Enzyme States#

def findInhibitorBound():
    A = PrevInaInh * functionDecay(RateOffInhibitor,TimeStep)
    if ThisDrug==2:
        B = PrevActUnb * functionGrowth((RateTransActive+RateOnSubstrate*SubConc+RateOnInhibitor*PrevInhConc),TimeStep)*RateOnInhibitor*PrevInhConc/(RateTransActive+RateOnSubstrate*SubConc+RateOnInhibitor*PrevInhConc)
    else:
        B = PrevInaUnb * functionGrowth( (RateOnInhibitor * PrevInhConc + RateTransInactive),TimeStep ) * RateOnInhibitor * PrevInhConc / (RateOnInhibitor * PrevInhConc + RateTransInactive)
    return A + B

def findInactiveUnbound():
    if ThisDrug==2:
        A = PrevInaUnb*functionDecay(RateTransInactive,TimeStep)
        B = PrevActUnb * functionGrowth((RateTransActive+RateOnSubstrate*SubConc+RateOnInhibitor*PrevInhConc),TimeStep)*RateTransActive/(RateTransActive+RateOnSubstrate*SubConc+RateOnInhibitor*PrevInhConc)
        C = 0
    else:
        A = PrevInaUnb*functionDecay((RateOnInhibitor*PrevInhConc+RateTransInactive),TimeStep)
        B = PrevActUnb * functionGrowth( (RateTransActive + RateOnSubstrate*SubConc),TimeStep) * RateTransActive / (RateTransActive + RateOnSubstrate * SubConc)
        C = PrevInaInh * functionGrowth(RateOffInhibitor,TimeStep)
    return A + B + C

def findActiveUnbound():
    if ThisDrug==2:
        A = PrevActUnb * functionDecay((RateTransActive+RateOnSubstrate*SubConc+RateOnInhibitor*PrevInhConc),TimeStep)
        B = PrevActBou * functionGrowth((RateCatalysis+RateOffSubstrate),TimeStep)
        C = PrevInaUnb * functionGrowth(RateTransInactive,TimeStep)
        D = PrevInaInh * functionGrowth(RateOffInhibitor,TimeStep)
    else:
        A = PrevActUnb*functionDecay((RateTransActive+RateOnSubstrate*SubConc),TimeStep)
        B = PrevInaUnb*functionGrowth((RateOnInhibitor*PrevInhConc+RateTransInactive),TimeStep)*RateTransInactive/(RateOnInhibitor*PrevInhConc+RateTransInactive)
        C = PrevActBou*functionGrowth((RateCatalysis+RateOffSubstrate),TimeStep)
        D = 0
    return A + B + C + D

def findSubstrateBound():
    A = PrevActBou*functionDecay((RateCatalysis+RateOffSubstrate),TimeStep)
    if ThisDrug==2:
        B = PrevActUnb * functionGrowth((RateTransActive+RateOnSubstrate*SubConc+RateOnInhibitor*PrevInhConc),TimeStep)*RateOnSubstrate*SubConc/(RateTransActive+RateOnSubstrate*SubConc+RateOnInhibitor*PrevInhConc)
    else:
        B = PrevActUnb*functionGrowth((RateTransActive+RateOnSubstrate*SubConc),TimeStep)*RateOnSubstrate*SubConc/(RateTransActive+RateOnSubstrate*SubConc)

    return A + B

################ Main Calculations ################

plt.clf()

for ThisDrug in DrugsWanted:

    SwitchDay = DrugDynamicTimings[ThisDrug][3]
    FBioavailability = DrugDynamicsConstants[ThisDrug][0]               # no units
    Dose = DrugDynamicsConstants[ThisDrug][1] * 10 ** (-5)              # mol
    VolumeOfDistribution = DrugDynamicsConstants[ThisDrug][2]           # litre
    EliminationRate = DrugDynamicsConstants[ThisDrug][3] * 10 ** (-2)   # per hour
    AbsorptionRate = DrugDynamicsConstants[ThisDrug][4]                 # per hour

    Gamma = calculateGamma(FBioavailability,Dose,AbsorptionRate,EliminationRate,VolumeOfDistribution)   # M or mol/L
    Epsilon = calculateEpsilon(EliminationRate) # no units
    Alpha = calculateAlpha(AbsorptionRate)      # no units

    for ThisMutation in MutationNumbersWanted:

        kCat = Mutant[ThisMutation][1]                          # per minute
        KM = Mutant[ThisMutation][2] * 10 ** (-6)               # M
        DeltaG = (WildTypeDeltaG + Mutant[ThisMutation][3]) /(6.022 * 10 ** 23)  # kcal
        IC50 = Mutant[ThisMutation][ThisDrug + 4] * 10 ** (-9)  # M

        RateCatalysis = convertRateToPerSec(kCat)                                  # per second
        RateTransActive = convertRateToPerSec(Rates[ThisMutation][0])              # per second
        RateOffSubstrate = convertRateToPerSec(Rates[ThisMutation][1])             # per second

        RateOffInhibitor = convertRateToPerSec(Rates[ThisMutation][ThisDrug+2])    # per second
      
        ExpNegDelG = math.exp(-Beta * DeltaG)   # no unit
        ExpPosDelG = math.exp(Beta * DeltaG)    # no unit
        OnePlusSOverKM = 1 + (SubConc / KM)     # no unit

        if ThisDrug==2:
            RD = IC50 * ExpNegDelG / (1 + ExpNegDelG * OnePlusSOverKM) # M

        else:
            RD = IC50 / (1 + ExpNegDelG * OnePlusSOverKM) # M

        RateTransInactive = RateTransActive * ExpNegDelG
        RateOnSubstrate = (RateOffSubstrate + RateCatalysis) / KM
        RateOnInhibitor = RateOffInhibitor / RD

        WeightInaUnb = 1                            # no unit
        WeightActUnb = ExpNegDelG                   # no unit
        WeightActBou = ExpNegDelG * SubConc / KM    # no unit

        NonInhWeightTotal = WeightInaUnb + WeightActUnb + WeightActBou
        
        # Setting up initial conditions
        
        InitActUnb = Etotal * WeightActUnb / NonInhWeightTotal
        InitActBou = Etotal * WeightActBou / NonInhWeightTotal
        InitActInh = 0
        InitInaUnb = Etotal * WeightInaUnb / NonInhWeightTotal
        InitInaInh = 0

        PrevActUnb = InitActUnb
        PrevActBou = InitActBou
        PrevActInh = InitActInh
        PrevInaUnb = InitInaUnb
        PrevInaInh = InitInaInh

        PrevInhConc = 0

        ActUnb = 0
        ActBou = 0
        ActInh = 0 # gonna be lazy and code this as InaInh for Dasatinib as everything is taken care of in the functions
        InaUnb = 0
        InaInh = 0
        
        ProductInTimeStep = 0
        TotalProduct = 0

        Time=0
        RecordedTime = []
        InhConc = []
        ProductRate= []
        SubBound = []
        Inhibited = []
        Inactive = []
        Active = []
        Total = []

        CurrentTotal=0

        SecondsInMinute = 0
        
        ProductRateFile = 'ProductRateDrug{:02d}Mutant{:02d}.dat'.format(ThisDrug,ThisMutation)
        OUTPUT = open(ProductRateFile,"w+")

        while Time<TimeTotal :

            
            if Time==0:
                ActBou = PrevActBou
                ActUnb = PrevActUnb
                InaUnb = PrevInaUnb
                InaInh = PrevInaInh
            else:
                ActBou = findSubstrateBound()
                ActUnb = findActiveUnbound()
                InaUnb = findInactiveUnbound()
                InaInh = findInhibitorBound()

            PrevInhConc=findInhibitorConcentration(Time,SwitchDay,Gamma,Epsilon,Alpha,AbsorptionRate,EliminationRate)
            if SecondsInMinute>=60:
                RecordedTime.append(Time/(24*3600))
                SubBound.append(ActBou)
                ProductRate.append(ActBou * kCat)
                OUTPUT.write('%e  %e'%((Time/(24*3600)),(ActBou * kCat)))
                Inhibited.append(InaInh)
                Inactive.append(InaUnb)
                Active.append(ActUnb)
                InhConc.append(PrevInhConc*10**9)
                CurrentTotal=ActBou+ActUnb+InaUnb+InaInh
                Total.append(CurrentTotal)
                OUTPUT.write('  %e  %e  %e  %e  %e\n'%((PrevInhConc*10**9),ActBou,ActUnb,InaUnb,InaInh))
                SecondsInMinute=0
            else:
                SecondsInMinute = SecondsInMinute + TimeStep


            if (CurrentTotal>=(Etotal*1.01)):
                Time=TimeTotal
            
            PrevActUnb = ActUnb
            PrevActBou = ActBou
            PrevInaUnb = InaUnb
            PrevInaInh = InaInh

            Time=Time+TimeStep

        if ThisMutation==0:
            outputfilename = 'InhibitorConcentrationDrug{:02d}.png'.format(ThisDrug)
            plt.plot(RecordedTime,InhConc, label=DrugName[ThisDrug])
            plt.xlabel('Time (days)', fontsize = 12)
            plt.ylabel('Inhibitor Concentration (nM)', fontsize = 12)
            plt.legend('',frameon=False)
            plt.savefig(outputfilename)
            plt.clf()

        outputfilename = 'ProductRateDrug{:02d}Mutant{:02d}.png'.format(ThisDrug,ThisMutation)
        plt.plot(RecordedTime,ProductRate)
        plt.xlabel('Time (days)', fontsize = 12)
        plt.ylabel('Product Rate', fontsize = 12)
        plt.legend('',frameon=False)
        plt.savefig(outputfilename)

        plt.clf()

        outputfilename = 'EnzymeStatesDrug{:02d}Mutant{:02d}.png'.format(ThisDrug,ThisMutation)
        plt.plot(RecordedTime,SubBound, label='Substrate bound')
        plt.plot(RecordedTime,Active, label='Active')
        plt.plot(RecordedTime,Inactive, label='Inactive')
        plt.plot(RecordedTime,Inhibited, label='Inhibited')
        plt.plot(RecordedTime,Total, label='Total')
        plt.xlabel('Time (days)', fontsize = 12)
        plt.ylabel('Enzyme Concentration', fontsize = 12)
        plt.legend()
        plt.savefig(outputfilename)

        plt.clf()



# Use allel environment
import pandas as pd
import pdb

#--- run: python3 Assignment4_Pandas.py

## Making table 1 into dt1

def get_ifr(age):
    if 0 <= age <= 17:
        return 20 / 1_000_000
    elif 18 <= age <= 49:
        return 500 / 1_000_000
    elif 50 <= age <= 64:
        return 6000 / 1_000_000
    else:  # age >= 65
        return 90000 / 1_000_000

ages = list(range(0, 101))
ifr_values = [get_ifr(age) for age in ages]

dt1 = pd.DataFrame({'Age': ages, 'IFR': ifr_values})

## Reading in csv as dt2
dt2 = pd.read_csv("WorldDemographics.csv")


## Combine dt1 and dt2 together
dt3= pd.merge(dt1, dt2, on='Age')

## Calculating the estimated number of deaths for each age
dt3['#Deaths'] = dt3['#Alive']*dt3['IFR']

## Group dataframe by country and calculating total number alive and dead
dt4 = dt3.groupby('PopulationID').agg({'#Alive':'sum','#Deaths':'sum'})

## Refining final dataframe that has: Total Population, #Deaths, %Died
dt4['Total Population'] = dt4['#Alive'] + dt4['#Deaths']
dt4 = dt4.drop('#Alive', axis=1)
dt4['%Died'] = dt4['#Deaths']/dt4['Total Population']


dt4.to_csv('8802_asm4.csv', index=True)

print("Done")

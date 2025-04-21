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


''' CSV file looks like:
PopulationID,#Deaths,Total Population,%Died
Afghanistan,117963.90656000002,39046167.90656,0.0030211391510248913
Albania,42576.25422,2920216.25422,0.01457982920219457
Algeria,310409.63931999996,44160612.63932,0.007029106272941864
Angola,85247.39461999996,32951411.39462,0.0025870635281473363
Antigua and Barbuda,961.9778999999999,98839.9779,0.009732680241726358
Argentina,513107.6403199999,45705812.64032,0.01122631041171676
Armenia,36173.1498,2999307.1498,0.012060501973734868
Aruba,1589.5340199999996,108310.53402,0.014675710302623799
Australia,406443.91866,25902138.91866,0.015691519528033892
Austria,170543.43970000002,9174646.4397,0.01858855715268049
Azerbaijan,76576.10684,10215552.10684,0.007496032132098581
Bahamas,3299.3082600000002,396495.30826,0.008321178564454775
'''

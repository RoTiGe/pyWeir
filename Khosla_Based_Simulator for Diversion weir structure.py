
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 04:02:44 2019

@author: Dr.Robel Tilaye Geressu
@email: robtilaye@gmail.com
"""
from platypus import *
import numpy as np
import pandas as pd
import math
from pandas import ExcelWriter
from pandas import ExcelFile
from platypus import Problem, Real
from platypus import NSGAII, NSGAIII,OMOPSO, DTLZ2,GDE3, Hypervolume, experiment, calculate, display
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from multiprocessing import Pool, freeze_support



    
def ReadInputs(InputId):
    inputfileName = 'User input Simple_'+str(InputId)+'.xlsx'     
    df = pd.read_excel(inputfileName, index_col=None, header=0,sheet_name='Hydraulic')     

    dfT  = df;

    HydraulicInput = np.array(dfT)

    Head= HydraulicInput[0,0]
    SafeExitGradient = HydraulicInput[1,0]
    q = HydraulicInput[2,0] 
    LacySiltFactor = HydraulicInput[3,0] 
    print('-------Hydraulic Input-----')
    
    for i in range(len(HydraulicInput)):
        print(HydraulicInput[i][1]," = ",HydraulicInput[i][0])  
      

    inputfileName = 'User input Simple_'+str(InputId)+'.xlsx'     
    df = pd.read_excel(inputfileName, index_col=None, header=0,sheet_name='StructuralParameters')    
    dfT  = df;

    StructuralParameterInputs=np.array(dfT)

    ApronLevel = StructuralParameterInputs[0,0]
    
    print('ApronLevel',ApronLevel)
    
    print('-------Structural Parameter Inputs--------')    
    
    for i in range(len(StructuralParameterInputs)):
        print(StructuralParameterInputs[i][1]," = ",StructuralParameterInputs[i][0])  
        
    return [Head,LacySiltFactor,SafeExitGradient,q,ApronLevel] 

def ReadCost():
  
    dfc = pd.ExcelFile ('CostInputs.xlsx')     

    dfTc  = [];
    dfTc = dfc.parse('Sheet1')
    costData = np.array(dfTc)
   
    PileCost = costData[0,0]
    ApronCost = costData[1,0]

    
    print ('----Cost Data------')
    for i in range(2):
        print(costData[i][1]," = ",costData[i][0])  
    
        
    return [PileCost,ApronCost]#,ExcavationCostPile,ExcavationCostApron]


def calcScour(q,siltfactor):
    '''
    

    Parameters
    ----------
    q : TYPE
        DESCRIPTION. discharge intensity (High flood dischage/ River Width)
    siltfactor : TYPE
        DESCRIPTION. a factor representing the resistance of the river bed material to distabilizing force of the water movement

    Returns
    -------
    Scourdepth : TYPE
        DESCRIPTION. expected depth of scoure (borrowing the river bed)

    '''
 
    Scourdepth = 1.35*pow(pow(q,2)/siltfactor,(1/3))  

    
    return Scourdepth

def calcPileDepths(Scourdepth,ScaleUSPd,ScaleDSPd):
    '''
    

    Parameters
    ----------
    Scourdepth : TYPE
        DESCRIPTION. expected depth of scoure (borrowing the river bed)
    ScaleUSPd : TYPE
        DESCRIPTION. a factor by which the minimimum downstream pile depth is increased from the minimum value (calculated from emperical observations/recommendations)
    ScaleDSPd : TYPE
        DESCRIPTION.a factor by which the minimimum upstream pile depth is increased from the minimum value (calculated from emperical observations/recommendations)

    Returns
    -------
    list
        DESCRIPTION. the upstream and downstream pile depths to be implemented

    '''


    MinDSPd = 1.5*Scourdepth
    MinUSPd = 1.25*Scourdepth    

    DSPd = MinDSPd*ScaleDSPd
    USPd = MinUSPd*ScaleUSPd  
    

    return [round(USPd,2),round(DSPd,2)]

def calcMinTotalApronLen(DSPd,Head,SafeExitGradient):
    '''
    
    Parameters
    ----------
    DSPd : TYPE
        DESCRIPTION. the downstream pile depths to be implemented
    Head : TYPE
        DESCRIPTION. the difference between upstream and downstream water levels
    SafeExitGradient : TYPE
        DESCRIPTION. the exit gradient that can be sustained by the river bed material before the start of progressive dislodging of soil particles (onset of piping failure)

    Returns
    -------
    MinTotalApronLength : TYPE
        DESCRIPTION. The minimum total apron length required to be built to make the exit gradient less than the safe amount

    '''

    lambda_i= pow(Head/(SafeExitGradient*DSPd*math.pi),2)

    alphai = np.sqrt(pow((2*lambda_i-1),2)-1)

    MinTotalApronLength = alphai*DSPd

    
    return MinTotalApronLength

def calcApronLens(ScaleTotalApronLength,MinTotalApronLength):
    '''
    

    Parameters
    ----------
    ScaleTotalApronLength : TYPE
        DESCRIPTION.a factor by which the minimimum total apron length is increased from the minimum value

    MinTotalApronLength : TYPE
        DESCRIPTION. The minimum total apron length required to be built to make the exit gradient less than the safe amount


    Returns
    -------
    TotalApronLength : TYPE
        DESCRIPTION. the total apron length to be implemented

    '''

    TotalApronLength = ScaleTotalApronLength*MinTotalApronLength
    
        
    return TotalApronLength
   
def calcExitGradient(totalapronLength,dnpld,Head):
    '''
    
    Parameters
    ----------
    totalapronLength : float
        DESCRIPTION. the total apron length to be implemented
    dnpld : TYPE
        DESCRIPTION. the downstream pile depths to be implemented
    Head : float
        DESCRIPTION. the difference between upstream and downstream water levels

    Returns
    -------
    ExitGradiant : float
        DESCRIPTION. the exit gradient that tends to dislodges the river bed material due to the residual pressure on soil particels at downstream of the weir structure


    '''
    
    alpha = totalapronLength/dnpld
    lambda_ = (1+np.sqrt(1+pow(alpha,2)))/2
    ExitGradiant = Head/(dnpld*math.pi*np.sqrt(lambda_))
    

    return ExitGradiant



def CalcKholsa(TotalApronLength, USPd, DSPd):
    '''
    

    DESCRIPTION. 

    Parameters
    ----------
    TotalApronLength : TYPE
        DESCRIPTION. the total apron length to be implemented
    USPd : TYPE
        DESCRIPTION. the upstream and downstream pile depths to be implemented
    DSPd : TYPE
        DESCRIPTION. the downstream pile depths to be implemented

    Returns
    -------
    list
        DESCRIPTION. residual pressures  as percentage of head difference between upstream and downstream water levels

    '''
    # Upstream pile
    alpha_up = TotalApronLength/USPd
    lambda_up = 0.5*(1 + np.sqrt(1 + pow(alpha_up,2)))
    phiE = math.acos((lambda_up - 2)/lambda_up)/math.pi
    phiD = math.acos((lambda_up - 1)/lambda_up)/math.pi
    
    phiC1 = 1 - phiE
    phiD1 = 1 - phiD

    # Downstream pile 
    alpha_dp = TotalApronLength/DSPd
    lambda_dp = 0.5*(1 + np.sqrt(1 + pow(alpha_dp,2)))
    phiE2= math.acos(((lambda_dp-2)/lambda_dp))/math.pi
    phiD2 = math.acos(((lambda_dp-1)/lambda_dp))/math.pi


    return [phiC1,phiD1,phiE2,phiD2]

def CalcKholsaCorrections(Phis, t1,t2,TotalApronLength, USPd, DSPd,pileWidth):
    '''
    

    Parameters
    ----------
    Phis : TYPE
        DESCRIPTION.
    t1 : TYPE
        DESCRIPTION.
    t2 : TYPE
        DESCRIPTION.
    TotalApronLength : TYPE
        DESCRIPTION.
    USPd : TYPE
        DESCRIPTION.
    DSPd : TYPE
        DESCRIPTION.
    pileWidth : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    '''
 
    phiC = Phis[0]
    phiD1 = Phis[1]
    phiE = Phis[2]
    phiD2 = Phis[3]
    
    # corrections for thinkness
    ct1  = ((phiD1 - phiC)/USPd)*t1 
    ct2 =  ((phiE - phiD2)/DSPd)*t2 

        
    # corrections for pile interference
    
    b1 = TotalApronLength - pileWidth
    h1 = 0
    
    DC = h1 - t2 +  DSPd
    dc = USPd - t1
    DE = h1 - t1+ USPd
    de = DSPd - t2

    try:
        cp1 = 0.19*np.sqrt((DC/b1))*((DC+dc)/TotalApronLength)
        cp2 =  0.19*np.sqrt((DE/b1))*((DE+de)/TotalApronLength)
    except:
        print('error')
        print('t1 ,t2, USPd, DSPd')
        print(t1 ,t2, USPd,DSPd)

    # Overal correction
    phiC_c = phiC  - ct1   + cp1
    phiE_c = phiE - ct2 - cp2


    return [phiC_c,phiE_c]

def calcUpliftForce(pres1, pres2, totalapronLength,Head):
    '''
    

    Parameters
    ----------
    pres1 : TYPE
        DESCRIPTION.
    pres2 : TYPE
        DESCRIPTION.
    totalapronLength : TYPE
        DESCRIPTION.
    Head : TYPE
        DESCRIPTION.

    Returns
    -------
    Uplift : TYPE
        DESCRIPTION.

    '''


    Uplift = (pres1 + pres2)*totalapronLength*Head*WaterDensity*gravity

    return Uplift

def calcApronWeight(totalapronLength,thinkness1, thinkness2):
    '''
    

    Parameters
    ----------
    totalapronLength : TYPE
        DESCRIPTION.
    thinkness1 : TYPE
        DESCRIPTION.
    thinkness2 : TYPE
        DESCRIPTION.

    Returns
    -------
    ApronWeight : TYPE
        DESCRIPTION.

    '''
    ApronWeight1 = totalapronLength*(thinkness2)*ConcreteDensity*gravity
    ApronWeight2 = totalapronLength*abs(thinkness1 - thinkness2)*ConcreteDensity*gravity/2
    ApronWeight = ApronWeight1 + ApronWeight2
    
    # print(1)
    # print('ApronWeight,ApronWeight1,ApronWeight2,totalapronLength,thinkness1,thinkness2,ConcreteDensity,gravity')
    # print(ApronWeight,ApronWeight1, ApronWeight2,totalapronLength,thinkness1,thinkness2,ConcreteDensity,gravity)
    return ApronWeight

    
def calcWaterWeight(Head,lengthOfwater):
    '''
    

    Parameters
    ----------
    Head : TYPE
        DESCRIPTION.
    lengthOfwater : TYPE
        DESCRIPTION.

    Returns
    -------
    WaterWeight : TYPE
        DESCRIPTION.

    '''
    
    WaterWeight = Head*lengthOfwater*WaterDensity*gravity
    # print(7)
    # print('WaterWeight,lengthOfwater')
    # print(WaterWeight,lengthOfwater)
    return WaterWeight

def calcHorizontalForce(Head):
    '''
    

    Parameters
    ----------
    Head : TYPE
        DESCRIPTION.

    Returns
    -------
    HorizontalForce : TYPE
        DESCRIPTION.

    '''
    
    HorizontalForce = Head*Head*WaterDensity*gravity*1/2
    
    # print(2)
    # print('HorizontalForce,Head,WaterDensity')
    # print(HorizontalForce,Head,WaterDensity)    
    return HorizontalForce

######################
def calcUpliftForceMoment(pres1, pres2, totalapronLength,Head):
    '''
    

    Parameters
    ----------
    pres1 : TYPE
        DESCRIPTION.
    pres2 : TYPE
        DESCRIPTION.
    totalapronLength : TYPE
        DESCRIPTION.
    Head : TYPE
        DESCRIPTION.

    Returns
    -------
    Uplift : TYPE
        DESCRIPTION.

    '''
    
    Uplift1 = ((pres1 + pres2)*totalapronLength*totalapronLength*Head/2)*WaterDensity*gravity   
    Uplift2 = (abs(pres1 - pres2)*totalapronLength*totalapronLength*Head/(2*3))*WaterDensity*gravity
    Uplift = Uplift1 + Uplift2
    
    # print(3)
    # print('Uplift , Uplift1 , Uplift2, pres1 , pres2')
    # print(Uplift , Uplift1 , Uplift2, pres1 , pres2)   
    
    return Uplift

def calcApronWeightMoment(totalapronLength,thinkness1, thinkness2):
    
    ApronWeightMoment1 = (gravity*ConcreteDensity*totalapronLength*thinkness2)*(totalapronLength/2 )     
    ApronWeightMoment2 = (gravity*ConcreteDensity*totalapronLength*abs(thinkness2-thinkness1)/2)*(totalapronLength/3)  
    ApronWeightMoment =     ApronWeightMoment1 +   ApronWeightMoment2    
    
    # print(4)
    # print('ApronWeightMoment,ApronWeightMoment1,ApronWeightMoment2,thinkness2,thinkness1 ')
    # print(ApronWeightMoment,ApronWeightMoment1,ApronWeightMoment2,thinkness2,thinkness1 )  
    
    return ApronWeightMoment
 
def calcHorizontalForceMoment(Head):
    
    HorizontalForceMoment = (Head/3)*(Head*WaterDensity*gravity*1/2)
    
    return HorizontalForceMoment

    
def calcWaterWeightMoment(Head,lengthOfwater):
    
    WaterWeightMoment = (Head*lengthOfwater*WaterDensity*gravity)*(lengthOfwater/2)
    
    return WaterWeightMoment

def calcMomentratio(StablizingMoment, DestablizingMoment):  
    
    Momentratio = StablizingMoment/DestablizingMoment
    # print(8)
    # print('Momentratio,StablizingMoment,DestablizingMoment')
    # print(Momentratio,StablizingMoment,DestablizingMoment)    
    
    return  Momentratio

def calcSlidingRatio(vertcialForce, HorizontalForce):  
    
    VerticalForcetratio = vertcialForce/HorizontalForce
    
    return  VerticalForcetratio

def calcEccentricity(StablizingMoment, DestablizingMoment,netVerticalForce,TotalApronLength):  
    goodOrBad = 1

    netMoment = StablizingMoment-DestablizingMoment

    Eccentricity = (netMoment/netVerticalForce) - (TotalApronLength/2)
    Centroid = netMoment/netVerticalForce
    RhoMin = (netVerticalForce/2)*(1-6*Eccentricity/TotalApronLength)
    RhoMax = (netVerticalForce/2)*(1+6*Eccentricity/TotalApronLength)
    if TotalApronLength/3 < Centroid  and Centroid < 2*TotalApronLength/3:
        goodOrBad = 0

    elif RhoMin < 0 or RhoMax > 70:
        goodOrBad = 1
    Eccentricity = abs(netMoment/netVerticalForce - TotalApronLength/2)/TotalApronLength

    return  Eccentricity,goodOrBad

    
def calcVerticalforceratio(UpliftForce, ApronWeight,WaterWeight):  
    
    Verticalforceratio = (ApronWeight+WaterWeight)/UpliftForce
    
    return  Verticalforceratio     

def calcCost(pilecost,apronmaterialcost, totalapronLength, USPd, DSPd,thinkness1, thinkness2):

    cost = pilecost*(USPd+DSPd) + apronmaterialcost*totalapronLength*(thinkness1+thinkness2)
   
    return round(cost,2)

def Headworks(USPd, DSPd, ScaleTotalApronLength, ScaleDownstreamApronLength, Thickness_Calc_1, Thickness_Calc_2, inputs,CostInputs):
    # print('USPd, DSPd, ScaleTotalApronLength, ScaleDownstreamApronLength, Thickness_Calc_1, Thickness_Calc_2')
    # print(USPd, DSPd, ScaleTotalApronLength, ScaleDownstreamApronLength, Thickness_Calc_1, Thickness_Calc_2)
    '''

    Parameters
    ----------
    USPd : TYPE
        DESCRIPTION. upstream pile depth
    DSPd : TYPE
        DESCRIPTION. downstream pile depth
    ScaleTotalApronLength : TYPE
        DESCRIPTION.  factor multiplying the minimum total apron length  (which is calculated based on the safe exit gradient
                                                                          the Head difference and the depth of downstream pile)
    ScaleDownstreamApronLength : TYPE
        DESCRIPTION. factor for adding the remaining apron length to the minimum downstream apron length
    Thickness_Calc_1 : TYPE
        DESCRIPTION. user specified thickeness of apron at upstream of the apron
    Thickness_Calc_2 : TYPE
        DESCRIPTION.user specified thickeness of apron at downstream of the apron
    inputs : TYPE
        DESCRIPTION. Hyrologic and geologic factors specific to the site
    CostInputs : TYPE
        DESCRIPTION. cost of pile construction and apron construction at the specific site

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    '''
    Head  = inputs[0]
    LacySiltFactor  = inputs[1]
    SafeExitGradient = inputs[2]
    q = inputs[3]
    ApronLevel = inputs[4]
    Scourdepth = calcScour(q,LacySiltFactor)
    MinDownstreamApronLength = 10
    pileWidth = 1.0

    
    MinTotalApronLen = calcMinTotalApronLen(DSPd,Head,SafeExitGradient)
    totalapronLength = ScaleTotalApronLength*MinTotalApronLen

    DownstreamApronLength = MinDownstreamApronLength + (totalapronLength-MinDownstreamApronLength)*ScaleDownstreamApronLength
    UpstreamApronLength = totalapronLength -  DownstreamApronLength
    
    # print(5)
    # print('totalapronLength,UpstreamApronLength,DownstreamApronLength')
    # print(totalapronLength,UpstreamApronLength,DownstreamApronLength)
    exitGradient = calcExitGradient(totalapronLength, DSPd, Head) 
    Kholsa_Calc = CalcKholsa(totalapronLength, USPd, DSPd)

    
    Kholsa_Calc_ = CalcKholsaCorrections(Kholsa_Calc, Thickness_Calc_1, Thickness_Calc_2,totalapronLength, USPd, DSPd,pileWidth)
    
    pres1 = Kholsa_Calc_[0]
    pres2 = Kholsa_Calc_[1]
    

    
    WaterWeight = calcWaterWeight(Head,UpstreamApronLength)

    UpliftForce = calcUpliftForce(pres1, pres2, totalapronLength,Head)
    ApronWeight = calcApronWeight(totalapronLength,Thickness_Calc_1, Thickness_Calc_2)
    
    WaterWeight = calcWaterWeight(Head,UpstreamApronLength)
    HorizontalForce = calcHorizontalForce(Head)
    
    UpliftForceMoment = calcUpliftForceMoment(pres1, pres2, totalapronLength,Head)
    ApronWeightMoment = calcApronWeightMoment(totalapronLength,Thickness_Calc_1, Thickness_Calc_2)
    HorizontalForceMoment = calcHorizontalForceMoment(Head)
    WaterWeightMoment = calcWaterWeightMoment(Head,UpstreamApronLength)

    Momentratio = calcMomentratio((ApronWeightMoment + WaterWeightMoment), (HorizontalForceMoment + UpliftForceMoment))
    SlidingRatio = calcSlidingRatio(WaterWeight + ApronWeight - UpliftForce, HorizontalForce)
    #def calcEccentricity(StablizingMoment, DestablizingMoment,netVerticalForce,TotalApronLength)
    Eccentricity, goodOrBad = calcEccentricity( (ApronWeightMoment + WaterWeightMoment),(HorizontalForceMoment + UpliftForceMoment),WaterWeight + ApronWeight - UpliftForce,totalapronLength)
    # print(10)
    # print('Eccentricity')
    # print(Eccentricity)
    Verticalforceratio = calcVerticalforceratio(UpliftForce, ApronWeight, WaterWeight)
    

    pilecost= CostInputs[0]
    apronmaterialcost = CostInputs[1]

    Cost = calcCost(pilecost,apronmaterialcost, totalapronLength, USPd, DSPd,Thickness_Calc_1, Thickness_Calc_2)
    Verticalforceratio = calcVerticalforceratio(UpliftForce, ApronWeight,WaterWeight)

    res =  [Cost,exitGradient,Verticalforceratio,Momentratio,SlidingRatio,Eccentricity,totalapronLength, UpstreamApronLength, DownstreamApronLength, Thickness_Calc_1, Thickness_Calc_2,USPd,DSPd]

    res = [round(rr,2) if not r == 1 else round(rr,3) for r,rr in enumerate(res) ]

    return res


def Gri_Search_Optimisations():
    USPd = 3.4
    DSPd = 4.10


    ScaleTotalApronLength = 1
    ScaleDownstreamApronLength = 0
    
# =============================================================================
#     Minthinkness1 = 0.2
#     Maxthinkness1 = 3
#     Minthinkness2 = 0.2
#     Maxthinkness2 = 4   
# =============================================================================
     
    Minthinkness1 = 0.2
    Maxthinkness1 = 3
    Minthinkness2 = 0.2
    Maxthinkness2 = 4       
    



    MinUSPd = 3.42
    MaxUSPd = 10
    MinDSPd = 4.11
    MaxDSPd = 10   
    
    MinScaleTotalApronLength = 1
    MaxScaleTotalApronLength = 3
    
    MinScaleDownstreamApronLength = 0
    MaxScaleDownstreamApronLength = 1 
    
    thinkness1_divisions = 10
    thinkness2_divisions = 10    
    
    USPd_divisions = 10
    DSPd_divisions = 10
    
    ScaleTotalApronLength_divisions = 10
    ScaleDownstreamApronLength_divisions = 10    
  
    
    Results = []
    for ii in range(thinkness1_divisions + 1):
        for jj in range(thinkness2_divisions + 1):
            for kk in range(USPd_divisions + 1):  
                for pp in range(DSPd_divisions + 1):   
                    for oo in range(ScaleTotalApronLength_divisions + 1):  
                        for mm in range(ScaleDownstreamApronLength_divisions + 1):   
                            
# =============================================================================
#         for jj in range(thinkness2_divisions + 1):
#             for kk in range(USPd_divisions): # + 1):  
#                 for pp in range(DSPd_divisions): # + 1):   
#                     for oo in range(ScaleTotalApronLength_divisions): # + 1):  
#                         for mm in range(ScaleDownstreamApronLength_divisions): # + 1):    
# =============================================================================
                            
                            thinkness1 = Minthinkness1 + (Maxthinkness1 - Minthinkness1)*ii/(thinkness1_divisions)
                            thinkness2 = Minthinkness2 + (Maxthinkness2 - Minthinkness2)*jj/(thinkness2_divisions)
            
                            USPd = MinUSPd + (MaxUSPd - MinUSPd)*kk/(USPd_divisions)
                            DSPd = MinDSPd + (MaxDSPd - MinDSPd)*pp/(DSPd_divisions)
                            
                            ScaleTotalApronLength = MinScaleTotalApronLength + (MaxScaleTotalApronLength - MinScaleTotalApronLength)*oo/(ScaleTotalApronLength_divisions)
                            ScaleDownstreamApronLength = MinScaleDownstreamApronLength + (MaxScaleDownstreamApronLength - MinScaleDownstreamApronLength)*mm/(ScaleDownstreamApronLength_divisions)
                            
                            res = Headworks(USPd, DSPd, ScaleTotalApronLength, ScaleDownstreamApronLength, thinkness1, thinkness2, inputs,CostInputs)
                        
                            res[2] = -1*res[2]
                            res[3] = -1*res[3]    
                            res[4] = -1*res[4]
                            
                
                            res[6:] = [USPd,DSPd,ScaleTotalApronLength,ScaleDownstreamApronLength,thinkness1,thinkness2,999,999,999,ii,jj,kk,pp,oo,mm] 
                            Results.append(res)
                    
    resDf = pd.DataFrame(Results,columns = ['Cost','exitGradient','Verticalforceratio','Momentratio','SlidingRatio', 'Eccentricity', 'USPd','DSPd','ScaleTotalApronLength','ScaleDownstreamApronLength', 'thinkness1'	,'thinkness2','rs',	'k','j','ii','jj','kk','pp','oo','mm'])
    # resDf.to_csv('solns/Grid_searched_Traditional_thicknessOnlySearched100.csv') 
    resDf.to_csv('solns/Grid_searched_TraditionalAllVariablesMillionGrids.csv')    
    
    


class theHeadWorksSimulator(Problem):


    def __init__(self):
        super(theHeadWorksSimulator, self).__init__(6, 6)


        myList = []
        for (ResNo, row) in RangeOfDecisionVariables.iterrows():
            myList.append(Real(row['min'],row['max']))

        self.types[:] = myList
        

                  
    def evaluate(self, solution):
        x = solution.variables[:]

        # USPd, DSPd, ScaleTotalApronLength, ScaleDownstreamApronLength, thinkness1, thinkness2, inputs,CostInputs)
        a = Headworks(round(x[0],2),round(x[1],2),round(x[2],2),round(x[3],2), round(x[4],2),round(x[5],2),inputs,CostInputs)


        solution.objectives[:] = [round(a[0],3), round(a[1],3), -1*round(a[2],3), -1*round(a[3],3), -1*round(a[4],3), round(a[5],3)]

def optimiseMOsolns():
    Myspacing = [];
    MyHv = [];
    Solns = [];
    b = []
    NoOfFEEPerGen = 1000;


    for rs in range (31,32):
 
        for k in range(10,11):
            print('t,rs,k = ',t,rs,k)

            algorithms = [(NSGAII, {"epsilons":0.01,"population_size":25})]
        
            problems =  [theHeadWorksSimulator()]
            results = experiment(algorithms, problems,  nfe=(k+1)*NoOfFEEPerGen, seeds=(rs)) 
            hyp = Hypervolume(minimum=[0.0,0.0,0.0,0.0,0.0,0.0], maximum=[150.1,0.3,20.0,20.0,20.0,1])
            hyp_result = calculate(results, hyp)
            

            for j in range(len(hyp_result['NSGAII']['theHeadWorksSimulator']['Hypervolume'])):
  
                MyHv.append([t,j,(k+1)*NoOfFEEPerGen,(hyp_result['NSGAII']['theHeadWorksSimulator']['Hypervolume'][j])])
        
            # display the results
            display(hyp_result,ndigits=4)
            
        
            for i, algorithm in enumerate(six.iterkeys(results)):
                for j in range(len(results[algorithm]['theHeadWorksSimulator'])):
                    result = results[algorithm]["theHeadWorksSimulator"][j]      
                    # print(j,result)
        
                    for a in result:
                        c = [a.objectives[0],a.objectives[1],a.objectives[2],
                              a.objectives[3],a.objectives[4],a.objectives[5],
                              a.variables[0],
                             a.variables[1], a.variables[2],a.variables[3],
                             a.variables[4],a.variables[5],t*5,(k+1)*NoOfFEEPerGen,j]
                        b.append(c)
        
                   
        my_df = pd.DataFrame(b)
    my_df.columns = ['Cost','exitGradient','Verticalforceratio','Momentratio','SlidingRatio','Eccentricity',
                     'USPd','DSPd','ScaleTotalApronLength', 'ScaleDownstreamApronLength',  'thinkness1', 'thinkness2', 'rs','k','j']

    mystring = 'solns//HeadWorks_HV.csv';
    mystring3 = 'solns//HeadWorks_soln.csv';
    
    my_dfh = pd.DataFrame(MyHv,columns = ['RS','j','Gen','Hv'])        
    my_df.to_csv(mystring3, index=False)         
    my_dfh.to_csv(mystring, index=False)  
if __name__ == "__main__":
    
    ConcreteDensity = 2240
    WaterDensity = 1000
    gravity = 9.81

    InputId = 1
        
    inputs  = ReadInputs(InputId)
    CostInputs = ReadCost()
    RangeOfDecisionVariables = pd.read_csv('RangeOfDecisionVariables.csv')
    
    # Gri_Search_Optimisations()
    
    optimiseMOsolns()
    

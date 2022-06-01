#include <iostream>
using namespace std;
#include "../headers/grove.hpp"
#include "../headers/commodity.hpp"
#include "../headers/behavior.hpp"
#include "../headers/coord.hpp"
#include "../headers/parameterSet.hpp"
#include "../headers/bioABM.h"
#include <math.h>


//#define SURVIVAL_DEBUG

#ifdef SURVIVAL_DEBUG
#define sdebug(x) cout << x
#else
#define sdebug(x) 
#endif



/********************************************************************
* GLOBALS
*********************************************************************/
Grove agents[ParameterSet::gridLength][ParameterSet::gridWidth];
Behavior* behaviorPatterns[ParameterSet::numBehaviorPatterns];
const int numRows = 33;
const int rowLength = 75;
vector<double> baselineSurvival;
double rowPsyllidsWeight;
double colPsyllidsWeight;
double rowPsyllidsMean;
double colPsyllidsMean;
double rowPsyllidsScale;
double colPsyllidsScale;
vector<boost::tuple<double, double>> survivalTimeBuckets;
string harvestFilename;
string yieldFilename;
double commodityPrice;
double commodityCost;
int commodityMaxAge;
string baselineFilename;
string timeBucketsFilename;
ofstream outputFile;
string outputFilename;
int experimentID;



boost::random::mt19937 econ_rng(std::time(0));
boost::random::uniform_01<> econ_gen;


/***********************************************************
* Initialize Survival Model
* Initializes the survival model
***********************************************************/
void initializeSurvivalModel() {
    double input, input2;
    //Baseline survival
    ifstream bFile(baselineFilename);
    while (bFile >> input) {
        baselineSurvival.push_back(input);
    }

    //Time buckets
    ifstream tFile(timeBucketsFilename);
    while (tFile >> input >> input2) {
        survivalTimeBuckets.push_back(boost::tuple<double, double>(input, input2));
    }

}


/*******************************************************
* Plan Behavior
********************************************************/
void planBehavior(Grove* g, Behavior* b) {
    if (b == behaviorPatterns[0]) {
        g->clearPlannedActions();
    }
    else {
        int spray1Target = 16;
        int spray2Target = 30;
        int maxGen = ParameterSet::groupWindow;
        if (b == behaviorPatterns[1]) {
            maxGen = ParameterSet::individualWindow;
        }
        int gen_ub = floor(maxGen / 2);
        int gen_lb = -1 * gen_ub;
        boost::uniform_int<> shootGen(gen_lb, gen_ub);

        //Spring Flush
        int offset1 = shootGen(econ_rng);
        int offset2 = shootGen(econ_rng);
        int springSpray1 = bioABM::getSpringStart() + spray1Target + offset1;
        int springSpray2 = bioABM::getSpringStart() + spray2Target + offset2;
        //Summer flush
        offset1 = shootGen(econ_rng);
        offset2 = shootGen(econ_rng);
        int summerSpray1 = bioABM::getSummerStart() + spray1Target + offset1;
        int summerSpray2 = bioABM::getSummerStart() + spray2Target + offset2;
        //Fall Flush
        offset1 = shootGen(econ_rng);
        offset2 = shootGen(econ_rng);
        int fallSpray1 = bioABM::getFallStart() + spray1Target + offset1;
        int fallSpray2 = bioABM::getFallStart() + spray2Target + offset2;

        plan_func spray = &sprayGrove;
        g->clearPlannedActions();
        g->getPlanningQ()[springSpray1] = spray;
        g->getPlanningQ()[springSpray2] = spray;
        g->getPlanningQ()[summerSpray1] = spray;
        g->getPlanningQ()[summerSpray2] = spray;
        g->getPlanningQ()[fallSpray1] = spray;
        g->getPlanningQ()[fallSpray2] = spray;


    }
}
/**************************************************************
* Initialize CHMA 
* Take in our grove grid and uniform crop and populate the grid
* *************************************************************/
void InitialiseCHMA(Commodity crop) {
    //Populating manually due to lambda parameter limitation
    int rboundInterval = numRows / 3;
    int rbound0 = 0;
    int rbound1 = rbound0 + rboundInterval;
    int rbound2 = rbound1 + rboundInterval;
    int rbound3 = rbound2 + rboundInterval;
    int cboundInterval = rowLength / 3;
    int cbound0 = 0;
    int cbound1 = cbound0 + cboundInterval;
    int cbound2 = cbound1 + cboundInterval;
    int cbound3 = cbound2 + cboundInterval;
    
    agents[0][0] = Grove(crop, ParameterSet::g00_agency, behaviorPatterns[ParameterSet::g00_behavior], rbound0, rbound1, cbound0, cbound1, ParameterSet::g00_lambda, ParameterSet::g00_alpha, ParameterSet::sprayingPopEff, ParameterSet::planningLength);
    agents[0][1] = Grove(crop, ParameterSet::g01_agency, behaviorPatterns[ParameterSet::g01_behavior], rbound0, rbound1, cbound1, cbound2, ParameterSet::g01_lambda, ParameterSet::g01_alpha, ParameterSet::sprayingPopEff, ParameterSet::planningLength);
    agents[0][2] = Grove(crop, ParameterSet::g02_agency, behaviorPatterns[ParameterSet::g02_behavior], rbound0, rbound1, cbound2, cbound3, ParameterSet::g02_lambda, ParameterSet::g02_alpha, ParameterSet::sprayingPopEff, ParameterSet::planningLength);
    agents[1][0] = Grove(crop, ParameterSet::g10_agency, behaviorPatterns[ParameterSet::g10_behavior], rbound1, rbound2, cbound0, cbound1, ParameterSet::g10_lambda, ParameterSet::g10_alpha, ParameterSet::sprayingPopEff, ParameterSet::planningLength);
    agents[1][1] = Grove(crop, ParameterSet::g11_agency, behaviorPatterns[ParameterSet::g11_behavior], rbound1, rbound2, cbound1, cbound2, ParameterSet::g11_lambda, ParameterSet::g11_alpha, ParameterSet::sprayingPopEff, ParameterSet::planningLength);
    agents[1][2] = Grove(crop, ParameterSet::g12_agency, behaviorPatterns[ParameterSet::g12_behavior], rbound1, rbound2, cbound2, cbound3, ParameterSet::g12_lambda, ParameterSet::g12_alpha, ParameterSet::sprayingPopEff, ParameterSet::planningLength);
    agents[2][0] = Grove(crop, ParameterSet::g20_agency, behaviorPatterns[ParameterSet::g20_behavior], rbound2, rbound3, cbound0, cbound1, ParameterSet::g20_lambda, ParameterSet::g20_alpha, ParameterSet::sprayingPopEff, ParameterSet::planningLength);
    agents[2][1] = Grove(crop, ParameterSet::g21_agency, behaviorPatterns[ParameterSet::g21_behavior], rbound2, rbound3, cbound1, cbound2, ParameterSet::g21_lambda, ParameterSet::g21_alpha, ParameterSet::sprayingPopEff, ParameterSet::planningLength);
    agents[2][2] = Grove(crop, ParameterSet::g22_agency, behaviorPatterns[ParameterSet::g22_behavior], rbound2, rbound3, cbound2, cbound3, ParameterSet::g22_lambda, ParameterSet::g22_alpha, ParameterSet::sprayingPopEff, ParameterSet::planningLength);
    
    agents[0][0].setBehaviorCost(2, ParameterSet::g00_premium);
    agents[0][1].setBehaviorCost(2, ParameterSet::g01_premium);
    agents[0][2].setBehaviorCost(2, ParameterSet::g02_premium);
    agents[1][0].setBehaviorCost(2, ParameterSet::g10_premium);
    agents[1][1].setBehaviorCost(2, ParameterSet::g11_premium);
    agents[1][2].setBehaviorCost(2, ParameterSet::g12_premium);
    agents[2][0].setBehaviorCost(2, ParameterSet::g20_premium);
    agents[2][1].setBehaviorCost(2, ParameterSet::g21_premium);
    agents[2][2].setBehaviorCost(2, ParameterSet::g22_premium);

  

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            planBehavior(&agents[i][j], agents[i][j].getBehavior());
        }
    }
    

     
    return;
}


/*********************************************
* Get Noise
* Takes in a lower and an upper bound with a
* max of 2 digits of precision, and returns a
* noise value within those bounds
* ********************************************/
double getNoise(double lb, double ub) {
    int lb_i = lb * 100;
    int ub_i = lb * 100;
    double pull = econ_gen(econ_rng);
    int intPull = pull * 100;
    intPull = intPull % (ub_i - lb_i);
    intPull = intPull + lb_i;
    double res = intPull / 100;
    return res;
}

/************************************************************************
* Get Survival
* Gets the expected survival probability based on a pre-calibrated
* semi-parametric cox proportional hazards model
************************************************************************/
double getSurvival(double rowPsyllids, double columnPsyllids, double t) {
    // 1. Scale data
    //cout << "1 - S( " << rowPsyllids << ", " << columnPsyllids << " ) = ";
    rowPsyllids = (rowPsyllids - rowPsyllidsMean) / rowPsyllidsScale;
    columnPsyllids = (columnPsyllids - colPsyllidsMean) / colPsyllidsScale;
    sdebug("Scaled row/col: " << rowPsyllids << " " << columnPsyllids << endl);

    //2. Calculate phi
    double dotProduct = (rowPsyllids * rowPsyllidsWeight) + (columnPsyllids * colPsyllidsWeight);
    sdebug("dot product: " << dotProduct << endl);
    double phi = exp(dotProduct);
    sdebug("phi: " << phi << endl;);

    //3. Determine index
    int minIndex;
    double minValue;
    for (int i = 0; i < survivalTimeBuckets.size(); i++) {
        double element = abs(survivalTimeBuckets[i].get<0>() - t);
        if (i == 0) {
            minIndex = i;
            minValue = element;
        }
        else if (element < minValue) {
            minIndex = i;
            minValue = element;
        }
    }

    sdebug("Min index: " << minIndex << endl);
    //4. Calculate Survival
    double survivalProbability = pow(baselineSurvival[minIndex], phi);
    //cout << 1.0 - survivalProbability << endl;
    sdebug("Survival: " << survivalProbability << endl);
    sdebug("ERisk: " << 1.0 - survivalProbability << endl);
    //5. Complete
    return survivalProbability;

}

/*********************************************************************
* Get Expected Infection
* Gets the expected infection probability at (c_i, c_j). When extension
* is true, all data is used to calculate parameters. Non extension
* calls only use in-grove information + 1 cell over the border
*********************************************************************/
double getExpectedInfection(int c_i, int c_j, bool extension = false, double* rPsyllids = NULL, double* cPsyllids = NULL) {
    //Calculate parameters
    double rowPsyllids = 0;
    double columnPsyllids = 0;
    if (extension) {
        for (int i = 0; i < numRows; i++) {
            if (i == c_i) { continue; }
            columnPsyllids += ((double)bioABM::getPsyllidsAt(i, c_j) / (double)abs(i - c_i));
        }
        for (int j = 0; j < rowLength; j++) {
            if (j == c_j) { continue; }
            rowPsyllids += ((double)bioABM::getPsyllidsAt(c_i, j) / (double)abs(j - c_j));
        }
    } 
    else {
        int i_lb = (c_i / (numRows / ParameterSet::gridLength)) * (numRows / ParameterSet::gridLength) - 1;
        i_lb = max(i_lb, 0);
        int i_ub = i_lb + (numRows / ParameterSet::gridLength) + 1;
        i_ub = min(i_ub, numRows);
        int j_lb = (c_j / (rowLength / ParameterSet::gridWidth)) * (rowLength / ParameterSet::gridWidth) - 1;
        j_lb = max(j_lb, 0);
        int j_ub = j_lb + (rowLength / ParameterSet::gridWidth) + 1;
        j_ub = min(j_ub, rowLength);
        for (int i = i_lb; i < i_ub; i++) {
            if (i == c_i) { continue; }
            columnPsyllids += ((double)bioABM::getPsyllidsAt(i, c_j) / (double)abs(i - c_i));
        }
        for (int j = j_lb; j < j_ub; j++) {
            if (j == c_j) { continue; }
            rowPsyllids += ((double)bioABM::getPsyllidsAt(c_i, j) / (double)abs(j - c_j));
        }
    }

    if (rPsyllids != NULL) {
        *rPsyllids = rowPsyllids;
    }
    if (cPsyllids != NULL) {
        *cPsyllids = columnPsyllids;
    }

    //Calculate expectation
    double survivalProbability = getSurvival(rowPsyllids, columnPsyllids, bioABM::getModelDay());

    //Return 1 - expected survival probability
    return (1.0 - survivalProbability);
}

/**********************************************************
* Get Extension Expectation
* Given a set of cells defined by bounds (lower inclusive,
* upper exclusive), return the maximum expected probability
* of infection
************************************************************/
double getGroveExpectation(int* ibounds, int* jbounds, bool extension = false, Grove* g = NULL) {
    double maxE = 0;
    double* rowPsyllids = new double(0);
    double* colPsyllids = new double(0);
    double maxR = 0;
    double maxC = 0;
    for (int i = ibounds[0]; i < ibounds[1]; i++) {
        for (int j = jbounds[0]; j < jbounds[1]; j++) {
            //cout << "Checking (" << i << ", " << j << "): ";
            double expectation = getExpectedInfection(i, j, extension, rowPsyllids, colPsyllids);
            if (expectation >= maxE) {
                maxR = *rowPsyllids;
                maxC = *colPsyllids;
                maxE = expectation;
                if (g != NULL) {
                    g->maxE_i = i;
                    g->maxE_j = j;
                }
            }
        }
    }
    if (extension && g != NULL) {
        g->lastColPsyllids = maxC;
        g->lastRowPsyllids = maxR;
    }
    delete rowPsyllids;
    delete colPsyllids;
    return maxE;
}


/*************************************************************
* Phase 1
* Psyllid growth, movement, and infection. All handled by the
* biological model
*************************************************************/
void Phase1() {
    bioABM::advanceBiologicalModel();
}

/*************************************************************
* Phase 2
* Execution of planned actions
************************************************************/
void Phase2() {
    int period_t = bioABM::getModelDay();
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            int relativePeriod = period_t % 365;
            plan_func action = agents[i][j].getAction(relativePeriod);
            if (action != NULL) {
                action(&agents[i][j]);
            }
        }
    }
}

/*********************************************************
 * Weibull Survival
 * Returns the survival probability given 
 * vector of characteristsics
 * *******************************************************/
double weibullSurvival(int t, string strategy, double efficacy, string groveID, double alpha) {
    //cout << "Weibull(" << t << ", " << strategy << ", " << efficacy << ", " << groveID << ", " << alpha << "): ";
    double intercept = 6.345658;
    double indv_coef = 0.025527;
    double group_coef = 0.021299;
    double e75_coef = -0.000461;
    double e85_coef = 0.004604;
    double alpha1_coef = 0.396100;
    double g01_coef = -0.135702;
    double g02_coef = 0.001482;
    double g10_coef = -0.288761;
    double g11_coef = -0.510001;
    double g12_coef = -0.270179;
    double g20_coef = -0.082489;
    double g21_coef = -0.242598;
    double g22_coef = -0.067161;
    double threeterm_coef = 0.183670;
    double e75alpha1_coef = 0.279668;
    double e85alpha1_coef = 1.003784;
    double scale = 0.549;

    double lp_alpha0 = intercept;
    bool indv = (strategy == "Individual Action");
    bool group = (strategy == "Group Action");
    bool e75 = (efficacy == 0.75);
    bool e85 = (efficacy == 0.85);
    lp_alpha0 += (indv * indv_coef) + (group * group_coef) + (e75 * e75_coef) + (e85 * e85_coef);
    
    if (groveID == "g01") {
        lp_alpha0 += g01_coef;
    }
    else if (groveID == "g02") {
        lp_alpha0 += g02_coef;
    }
    else if (groveID == "g10") {
        lp_alpha0 += g10_coef;
    }
    else if (groveID == "g11") {
        lp_alpha0 += g11_coef;
    }
    else if (groveID == "g12") {
        lp_alpha0 += g12_coef;
    }
    else if (groveID == "g20") {
        lp_alpha0 += g20_coef;
    }
    else if (groveID == "g21") {
        lp_alpha0 += g21_coef;
    }
    else if (groveID == "g22") {
        lp_alpha0 += g22_coef;
    }
    
    double lp_alpha1 = lp_alpha0;
    bool threeTerm = (e85 && (indv || group));
    lp_alpha1 += alpha1_coef + (e75 * e75alpha1_coef) + (e85 * e85alpha1_coef) + (threeTerm * threeterm_coef);

    double survival_alpha0 = exp(-1 * pow(t / exp(lp_alpha0), 1/scale));
    double survival_alpha1 = exp(-1 * pow(t / exp(lp_alpha1), 1/scale));
    double scaled = survival_alpha0 + alpha*(survival_alpha1 - survival_alpha0);
    //cout << scaled << endl;
    return scaled;
}
/**********************************************************
* Get Expected Risk
* Calculates the growers expected risk of infection based on
* their expectation, the extension agents expectation, and
* their trust in the extension agent
***********************************************************/
double getExpectedRisk(Grove*g, int growerI, int growerJ) {
    string growerID = string("g") + to_string(growerI) + to_string(growerJ);
    string behaviorName;
    if (g->getBehavior() == behaviorPatterns[0]) {
        behaviorName = "No Action";
    }
    else if (g->getBehavior() == behaviorPatterns[1]) {
        behaviorName = "Individual Action";
    }
    else {
        behaviorName = "Group Action";
    }
    double extensionExpectation = 1 - weibullSurvival(bioABM::getModelDay(), behaviorName, 
                                                  g->getSprayEfficacy(), growerID,
                                                  g->getAlpha());
    double adjustedExpectation = g->getLambda() * (extensionExpectation - g->lastAdjustedRisk) + g->lastAdjustedRisk;
    g->lastGrowerRisk = adjustedExpectation;
    g->lastExtensionRisk = extensionExpectation;
    g->lastAdjustedRisk = adjustedExpectation;
    return adjustedExpectation;
}


double getMeanHLB(Grove g) {
    int* ibounds = g.getIBounds();
    int* jbounds = g.getJBounds();
    double totalCells = 0.0;
    double totalHLB = 0.0;
    for (int i = ibounds[0]; i < ibounds[1]; i++) {
        for (int j = jbounds[0]; j < jbounds[1]; j++) {
            totalCells += 1.0;
            totalHLB += bioABM::getSeverityAt(i, j);
        }
    }
    return totalHLB / totalCells;
    
}
/********************************************************
* Phase 3
* Behavior Determination
********************************************************/
void Phase3() {
    //Only determine during a planning period
    int period_t = bioABM::getModelDay();
    if ( (period_t % ParameterSet::planningLength) != 0 ) { return; }

    

    //Determine behavior for agents with agency
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            if (!agents[i][j].hasAgency()) { continue; }
            //Gather info for assessing risk
            //Assess risk
            double risk = getExpectedRisk(&agents[i][j], i, j);
            double meanHLB = getMeanHLB(agents[i][j]);
            bool findsHLB = econ_gen(econ_rng) <= meanHLB;
            if (findsHLB || agents[i][j].foundHLB) {
                risk = 1.0;
                if (!agents[i][j].foundHLB) {
                    agents[i][j].foundHLB = true;
                    agents[i][j].foundHLB_day = bioABM::getModelDay();
                }
            }
            agents[i][j].lastAdjustedRisk = risk;
            double maxExpectedValue = 0;
            int maxEVIndex = -1;
            int totalProjectingDays = ParameterSet::projectionLength + bioABM::getModelDuration();
            for (int k = 0; k < ParameterSet::numBehaviorPatterns; k++) {
                //additional costs are per planning period, given 6 sprays
                double additionalCosts = agents[i][j].behaviorCosts[k] * 6 / 4;
                double EV = behaviorPatterns[k]->getExpectedValue(agents[i][j], risk, totalProjectingDays, 
                                                                period_t, ParameterSet::planningLength,
                                                                agents[i][j].getSprayEfficacy(), agents[i][j].getAlpha(),
                                                                additionalCosts);
                if ( (EV > maxExpectedValue) || k == 0 ) {
                    maxExpectedValue = EV;
                    maxEVIndex = k;
                }
                switch (k) {
                    case 0:
                        agents[i][j].lastNAEV = EV;
                        break;
                    case 1:
                        agents[i][j].lastISEV = EV;
                        break;
                    case 2:
                        agents[i][j].lastGSEV = EV;
                        break;
                    default:
                        break;
                }
            }
            if (agents[i][j].getBehavior() != behaviorPatterns[maxEVIndex]) {
                planBehavior(&agents[i][j], behaviorPatterns[maxEVIndex]);
                agents[i][j].setBehavior(behaviorPatterns[maxEVIndex]);
       
            }
        }
    }

}

/********************************************************
* Phase 4
* Plan Actions
*********************************************************/
void Phase4() {
    //only plan manually during the new year
    int period_t = bioABM::getModelDay();
    if ( (period_t % 365) != 0 ) { return; }
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            planBehavior(&agents[i][j], agents[i][j].getBehavior());
        }
    }
}

/***********************************************
* Phase 5: Accounting
* Economic accounting for crops and mitigation
***********************************************/
void Phase5() {
    int t = bioABM::getModelDay();
    int rel_t = t % 365;
    
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            int* ibounds = agents[i][j].getIBounds();
            int* jbounds = agents[i][j].getJBounds();
            int numCrops = (ibounds[1] - ibounds[0]) * (jbounds[1] - jbounds[0]);

            if (rel_t == 0) {
                //VC
                //agents[i][j].costs += numCrops * agents[i][j].getCrop()->getVariableCost();
                //FC
                //agents[i][j].costs += agents[i][j].getFixedCosts();
                agents[i][j].costs += agents[i][j].getCrop()->costs;
            }

            //Harvest
            if (agents[i][j].getCrop()->isHarvestPeriod(rel_t)) {
                for (int k = ibounds[0]; k < ibounds[1]; k++) {
                    for (int l = jbounds[0]; l < jbounds[1]; l++) {
                        //Age of crop at time in projection
                        //int age = bioABM::getAgeAt(k, l);
                        //Projected severity based on days since initial infection
                        double severity = bioABM::getSeverityAt(k, l);
                        //Yield of crop at projected age
                        double returns = agents[i][j].getCrop()->getReturns();
                        //Infected yield: Units yielded times projected decay
                        double adjustedReturns = returns * getInfectedYield(severity);
                        agents[i][j].returns += adjustedReturns;

                    }
                }
            }

            //Mitigation costs
            if (rel_t % ParameterSet::planningLength == 0) {
                agents[i][j].costs += agents[i][j].getBehavior()->getVariableCosts();
            }
        }
    }
    
    
    return;
}

/*
void initializeSQL() {
  driver = get_driver_instance();
  con = driver->connect("tcp://127.0.0.1:3306","sam","CitrusABM21");
  con->setSchema("citrus");
}*/
/*************************************************************
* Write CSV Line
***************************************************************/
/*
void writeSQLLine() {
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            sql::Statement *stmt;
            std::stringstream cmd;
            cmd << "INSERT INTO econ VALUES(";
            string behaviorID;
            if (agents[i][j].getBehavior() == behaviorPatterns[0]) {
                behaviorID = "1";
            }
            else if (agents[i][j].getBehavior() == behaviorPatterns[1]) {
                behaviorID = "2";
            }
            else {
                behaviorID = "3";
            }
            cmd << bioABM::getModelDay() << ","; 
            cmd << i << ",";
            cmd << j << ",";
            cmd << behaviorID << ",";
            cmd << agents[i][j].costs << ",";
            cmd << agents[i][j].returns << ",";
            cmd << (agents[i][j].returns - agents[i][j].costs) << ",";
            cmd << agents[i][j].getLambda() << ",";
            cmd << agents[i][j].lastNAEV << ",";
            cmd << agents[i][j].lastISEV << ",";
            cmd << agents[i][j].lastGSEV << ",";
            cmd << agents[i][j].lastAdjustedRisk << ",";
            cmd << agents[i][j].getAlpha() << ",";
            cmd << experimentID << ");";
            stmt = con->createStatement();
            stmt->execute(cmd.str());
            delete stmt;
        }
    }
}*/

/*************************************************************
* Write CSV Line
***************************************************************/
void writeCSVLine() {
    //compute total alpha
    double alphaCount = 0;
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            if (agents[i][j].getBehavior() == behaviorPatterns[2]) {
                alphaCount += 1;
            }
        }
    }
    for (int i = 0; i < ParameterSet::gridLength; i++) {
        for (int j = 0; j < ParameterSet::gridWidth; j++) {
            string behaviorID;
            if (agents[i][j].getBehavior() == behaviorPatterns[0]) {
                behaviorID = "1";
            }
            else if (agents[i][j].getBehavior() == behaviorPatterns[1]) {
                behaviorID = "2";
            }
            else {
                behaviorID = "3";
            }

            double totalSeverity = 0;
            for (int bio_i = 11*i; bio_i < 11*(i+1); bio_i++) {
                for (int bio_j = 25*j; bio_j < 25*(j+1); bio_j++) {
                    totalSeverity += bioABM::getSeverityAt(bio_i,bio_j);
                }
            }
            double spraying = (behaviorID == "3");
            double meanSeverity = totalSeverity / (275);

            outputFile << bioABM::getModelDay() << ","; 
            outputFile << i << ",";
            outputFile << j << ",";
            outputFile << "g" << i << j << ",";
            outputFile << behaviorID << ",";
            outputFile << agents[i][j].costs << ",";
            outputFile << agents[i][j].returns << ",";
            outputFile << (agents[i][j].returns - agents[i][j].costs) << ",";
            outputFile << agents[i][j].getLambda() << ",";
            outputFile << agents[i][j].lastNAEV << ",";
            outputFile << agents[i][j].lastISEV << ",";
            outputFile << agents[i][j].lastGSEV << ",";
            outputFile << agents[i][j].lastAdjustedRisk << ",";
            outputFile << agents[i][j].getAlpha() << ",";
            outputFile << (alphaCount - spraying) / 8.0 << ",";
            outputFile << meanSeverity << ",";
            outputFile << (totalSeverity > 0) << ",";
            outputFile << experimentID << ",";
            outputFile << agents[i][j].foundHLB << ",";
            outputFile << agents[i][j].foundHLB_day << ",";
            outputFile << agents[i][j].behaviorCosts[2] << endl;
        }
    }
}
/*************************************************************
* Run Model
* Runs through the phases of the model until the simulation is
* complete
**************************************************************/
void runModel() {
    outputFile.open(outputFilename);
    outputFile
       << fixed << "t,i,j,grove_id,behavior,costs,returns,profit,lambda,lastNAEV,lastISEV,lastGSEV,lastAdjustedRisk,alphaP,alphaA,hlb_severity,infected,experiment_id,foundHLB,foundHLB_day,premium" << endl;
    while (bioABM::getModelDay() <= bioABM::getModelDuration()) {
        // Stage 1: Psyllid Growth and Movement
        Phase1();
        //cout << "Period " << bioABM::getModelDay() << endl;
        if (!ParameterSet::biologicalRun) {
            if (bioABM::getModelDay() % ParameterSet::planningLength == 0) {
              //cout << "Planning period!\n";
            }
            // Stage 2: Execution of Planned Actions
            Phase2();

            // Stage 3: Behavior determination
            Phase3();

            //Stage 4: Planning
            Phase4();

            //Stage 5: Accounting
            Phase5();
            writeCSVLine();
        }
    }
    outputFile.close();
}

/******************************************************
* Parse Parameter File
* Parses a parameter file 
* ****************************************************/
void parseParameterFile(string fileName) {
    ifstream is(fileName);
    cereal::JSONInputArchive archive(is);
    try {
        archive(ParameterSet::planningLength, ParameterSet::sprayingPopEff, ParameterSet::freshYield, ParameterSet::juiceYield,
            ParameterSet::freshPrice, ParameterSet::juicePrice, ParameterSet::costs, ParameterSet::biologicalRun,
            ParameterSet::numIndividualSprays, ParameterSet::numGroupSprays, ParameterSet::projectionLength, ParameterSet::sprayCost,
            ParameterSet::groupWindow, ParameterSet::individualWindow,
            ParameterSet::g00_lambda, ParameterSet::g01_lambda, ParameterSet::g02_lambda,
            ParameterSet::g10_lambda, ParameterSet::g11_lambda, ParameterSet::g12_lambda,
            ParameterSet::g20_lambda, ParameterSet::g21_lambda, ParameterSet::g22_lambda,
            ParameterSet::g00_alpha, ParameterSet::g01_alpha, ParameterSet::g02_alpha,
            ParameterSet::g10_alpha, ParameterSet::g11_alpha, ParameterSet::g12_alpha,
            ParameterSet::g20_alpha, ParameterSet::g21_alpha, ParameterSet::g22_alpha,
            ParameterSet::g00_behavior, ParameterSet::g01_behavior, ParameterSet::g02_behavior,
            ParameterSet::g10_behavior, ParameterSet::g11_behavior, ParameterSet::g12_behavior,
            ParameterSet::g20_behavior, ParameterSet::g21_behavior, ParameterSet::g22_behavior,
            ParameterSet::g00_agency, ParameterSet::g01_agency, ParameterSet::g02_agency,
            ParameterSet::g10_agency, ParameterSet::g11_agency, ParameterSet::g12_agency,
            ParameterSet::g20_agency, ParameterSet::g21_agency, ParameterSet::g22_agency,
            ParameterSet::g00_premium, ParameterSet::g01_premium, ParameterSet::g02_premium,
            ParameterSet::g10_premium, ParameterSet::g11_premium, ParameterSet::g12_premium,
            ParameterSet::g20_premium, ParameterSet::g21_premium, ParameterSet::g22_premium,
            harvestFilename, baselineFilename, timeBucketsFilename, rowPsyllidsWeight, colPsyllidsWeight,
            rowPsyllidsMean, colPsyllidsMean, rowPsyllidsScale, colPsyllidsScale, outputFilename, experimentID);
    }
    catch (exception e) {
        cout << "ERROR WITH ECON JSON" << endl;
        exit(-1);
    }
}   

Commodity getCommodity() {
    //get yield profile
   

    //get harvests
    vector<int> harvests;
    ifstream hFile(harvestFilename);
    int harvest;
    while (hFile >> harvest) {
        harvests.push_back(harvest);
    }

    Commodity retVal(ParameterSet::freshYield, ParameterSet::juiceYield, ParameterSet::freshPrice, ParameterSet::juicePrice, ParameterSet::costs, harvests);
    return retVal;
}


void testSurvivalFunction() {
    string input;
    string input2;
    string input3;
    while (true) {
        cin >> input;
        cin >> input2;
        cin >> input3;
        double rowPsyllids = atof(input.c_str());
        double colPsyllids = atof(input2.c_str());
        int time = atoi(input3.c_str());
        double result = getSurvival(rowPsyllids, colPsyllids, time);
        cout << result << endl;
    }
}


double getProfitCurves(double premium) {
    GroupAction groupSpray = GroupAction(21, ParameterSet::sprayCost,  premium, ParameterSet::planningLength);
    behaviorPatterns[2] = &groupSpray;
    //ofstream pOut;
    //pOut.open("C:/dev/ECNR.csv");
    int projectionLength = 3650;
    int planningLength = 91;
    //pOut << "efficacy,risk,alpha,t,strategy,ev,er,c,foundHLB,foundHLB_day" << endl;
    //No found HLB
    int totalGroupChoices = 0;
    int totalChoices = 0;
    for (double efficacy = 0.85; efficacy <= 0.85; efficacy += 0.1) {
        for (double risk = 0.1; risk <= 1; risk += 0.1) {
            for (double alpha=0; alpha <=1; alpha +=0.25) {
                for (int t = 91; t <= 1825; t+=91) {
                    bool choseGroup = false;
                    double ev = 0;
                    double strategy = 0;
                    double* na = behaviorPatterns[0]->getExpectedValueTester(agents[0][0], risk, projectionLength, t, planningLength, efficacy, alpha);
                    ev = na[0];
                    strategy = 0;
                    double* ia = behaviorPatterns[1]->getExpectedValueTester(agents[0][0], risk, projectionLength, t, planningLength, efficacy, alpha);
                    if (ia[0] >= ev) {
                        ev = ia[0];
                        strategy = 1;
                    }
                    double* ga = behaviorPatterns[2]->getExpectedValueTester(agents[0][0], risk, projectionLength, t, planningLength, efficacy, alpha);
                    if (ga[0] >= ev) {
                        ev = ga[0];
                        strategy = 2;
                        totalGroupChoices++;
                    }
                    /*
                    //na
                    pOut    << efficacy << ","
                            << risk << ","
                            << alpha << ","
                            << t << ","
                            << 0 << ","
                            << na[0] << ","
                            << na[2] << ","
                            << na[1] << ","
                            << false << ","
                            << -1 << endl;
                    //ia
                    pOut    << efficacy << ","
                            << risk << ","
                            << alpha << ","
                            << t << ","
                            << 1 << ","
                            << ia[0] << ","
                            << ia[2] << ","
                            << ia[1] << ","
                            << false << ","
                            << -1 << endl;
                    //ga
                    pOut    << efficacy << ","
                            << risk << ","
                            << alpha << ","
                            << t << ","
                            << 2 << ","
                            << ga[0] << ","
                            << ga[2] << ","
                            << ga[1] << ","
                            << false << ","
                            << -1 << endl;*/
                    totalChoices++;
                }
            }
        }
    }

    //Found HLB
    agents[0][0].foundHLB = true;
    for (double efficacy = 0.85; efficacy <= 0.85; efficacy += 0.1) {
        for (double risk = 0.1; risk <= 1; risk += 0.1) {
            for (double alpha=0; alpha <=1; alpha +=0.25) {
                for (int foundDay = 91; foundDay <= 1825; foundDay += 91) {
                    agents[0][0].foundHLB_day = foundDay;
                    for (int t = foundDay + 91; t <= 1825; t+=91) {
                        double ev = 0;
                        double strategy = 0;
                        double* na = behaviorPatterns[0]->getExpectedValueTester(agents[0][0], risk, projectionLength, t, planningLength, efficacy, alpha);
                        ev = na[0];
                        strategy = 0;
                        double* ia = behaviorPatterns[1]->getExpectedValueTester(agents[0][0], risk, projectionLength, t, planningLength, efficacy, alpha);
                        if (ia[0] >= ev) {
                            ev = ia[0];
                            strategy = 1;
                        }
                        double* ga = behaviorPatterns[2]->getExpectedValueTester(agents[0][0], risk, projectionLength, t, planningLength, efficacy, alpha);
                        if (ga[0] >= ev) {
                            ev = ga[0];
                            strategy = 2;
                            totalGroupChoices++;
                        }
                        
                        //na
                        /*
                        pOut    << efficacy << ","
                                << risk << ","
                                << alpha << ","
                                << t << ","
                                << 0 << ","
                                << na[0] << ","
                                << na[2] << ","
                                << na[1] << ","
                                << true << ","
                                << foundDay << endl;
                        //ia
                        pOut    << efficacy << ","
                                << risk << ","
                                << alpha << ","
                                << t << ","
                                << 1 << ","
                                << ia[0] << ","
                                << ia[2] << ","
                                << ia[1] << ","
                                << true << ","
                                << foundDay << endl;
                        //ga
                        pOut    << efficacy << ","
                                << risk << ","
                                << alpha << ","
                                << t << ","
                                << 2 << ","
                                << ga[0] << ","
                                << ga[2] << ","
                                << ga[1] << ","
                                << true << ","
                                << foundDay << endl;*/
                        totalChoices++;
                    }
                }
            }
        }
    }
    
    //pOut.close();
    return (double)totalGroupChoices / (double)totalChoices;
}

void findProfitBounds() {
    ofstream out;
    out.open("C:/dev/profitBounds.csv");
    out << "premium,groupPct" << endl; 
    double groupChoice = -1;
    double premium = -1;
    while (groupChoice != 0) {
        premium += 1;
        groupChoice = getProfitCurves(premium);
        cout << premium << ": " << groupChoice << endl;
        out << premium << "," << groupChoice << endl;
    }
    out.close();
}

/*****************************************************
* Main
*****************************************************/
int main(int argc, char ** argv) {
    string econConfigFile;
    string bioConfigFile;

    if (argc == 3) {
        econConfigFile = argv[1];
        bioConfigFile = argv[2];
        
    } 
    else {
        cout << "Using default filenames\n";
        econConfigFile = "C:/dev/EconABM/configs/econConfig.json";
        bioConfigFile = "C:/dev/EconABM/configs/bioConfig.json";
    }
    parseParameterFile(econConfigFile);
    bioABM::parseParameterFile(bioConfigFile);


    // Set up the behavior patterns, parameters arent used
    IndividualAction individualSpray = IndividualAction(60, ParameterSet::sprayCost, ParameterSet::planningLength);
    
    GroupAction groupSpray = GroupAction(21, ParameterSet::sprayCost,  0, ParameterSet::planningLength);

    NoAction noAction = NoAction();

    //List of possible behaviors
    behaviorPatterns[0] = &noAction;
    behaviorPatterns[1] = &individualSpray;
    behaviorPatterns[2] = &groupSpray;
    bioABM::setExperimentID(experimentID);
    InitialiseCHMA(getCommodity());
    runModel();
    //findProfitBounds();

    return 0;
}

#include "../headers/behavior.hpp"




/***************************************************************
* Get Infected yield
* Relates HLB severity to yield. From Bassanezi et. al (2011)
******************************************************************/
double getInfectedYield(double hlbSeverity) {
    return exp(-1.85 * hlbSeverity);
}

/***************************************************************
* Spray Grove
***************************************************************/
void sprayGrove(Grove* g) {
    vector<boost::tuple<int, int>> coords;
    for (int i = g->getIBounds()[0]; i < g->getIBounds()[1]; i++) {
        for (int j = g->getJBounds()[0]; j < g->getJBounds()[1]; j++) {
            coords.push_back(boost::tuple<int, int>(i, j));
        }
    }
    bioABM::sprayTrees(g->getSprayEfficacy(), coords);
}

/*********************************************************
* Behavior : Get Expected Value
* Returns the expected value of a behavior pattern
**********************************************************/
double Behavior::getExpectedValue(Grove g, double risk, int projectionLength, int startingPeriod, 
                                  int planningLength, double sprayEfficacy, double alpha, double additionalCosts) {
    double EV = 0;
    double ui_outcome = 0;
    double i_outcome = 0;
    int* ibounds = g.getIBounds();
    int* jbounds = g.getJBounds();
   
    for (int t = startingPeriod; t < projectionLength; t++) {

        //Add variable costs from mitigation every planning period
        if (t % planningLength == 0) {
            ui_outcome -= this->getVariableCosts();
            i_outcome -= this->getVariableCosts();
            ui_outcome -= additionalCosts;
            i_outcome -= additionalCosts;
        }

        //Add crop income every harvest
        if (g.getCrop()->isHarvestPeriod(t % 365)) {
            int numCrops = (g.getIBounds()[1] - g.getIBounds()[0]) * (g.getJBounds()[1] - g.getJBounds()[0]);
            double severity;
            if (!g.foundHLB) {
                severity = this->hlbSpread(t - startingPeriod, sprayEfficacy, alpha);
            }
            else {
                severity = this->hlbSpread(t - g.foundHLB_day, sprayEfficacy, alpha);
            }
         
            double returns = g.getCrop()->getReturns();
            double infectedReturns = returns * getInfectedYield(severity);
            ui_outcome += returns * numCrops;
            i_outcome += infectedReturns * numCrops;
        }

        //Add variable costs from crops and fixed costs on land at the end of the year
        if (t % 365 == 0) {
            ui_outcome -= g.getCrop()->costs;
            i_outcome -= g.getCrop()->costs;
        }
    }
    EV = (risk * i_outcome) + ((1 - risk) * ui_outcome);
    return EV;
}

double* Behavior::getExpectedValueTester(Grove g, double risk, int projectionLength, int startingPeriod, 
                                  int planningLength, double sprayEfficacy, double alpha) {
    double EV = 0;
    double ui_outcome = 0;
    double i_outcome = 0;
    double costs = 0;
    int* ibounds = g.getIBounds();
    int* jbounds = g.getJBounds();
   
    for (int t = startingPeriod; t < projectionLength; t++) {

        //Add variable costs from mitigation every planning period
        if (t % planningLength == 0) {
            ui_outcome -= this->getVariableCosts();
            i_outcome -= this->getVariableCosts();
            costs += this->getVariableCosts();
        }

        //Add crop income every harvest
        if (g.getCrop()->isHarvestPeriod(t % 365)) {
            int numCrops = (g.getIBounds()[1] - g.getIBounds()[0]) * (g.getJBounds()[1] - g.getJBounds()[0]);
            double severity;
            if (!g.foundHLB) {
                severity = this->hlbSpread(t - startingPeriod, sprayEfficacy, alpha);
            }
            else {
                severity = this->hlbSpread(t - g.foundHLB_day, sprayEfficacy, alpha);
            }
         
            double returns = g.getCrop()->getReturns();
            double infectedReturns = returns * getInfectedYield(severity);
            ui_outcome += returns * numCrops;
            i_outcome += infectedReturns * numCrops;
        }

        //Add variable costs from crops and fixed costs on land at the end of the year
        if (t % 365 == 0) {
            ui_outcome -= g.getCrop()->costs;
            i_outcome -= g.getCrop()->costs;
            costs += g.getCrop()->costs;
        }
    }
    EV = (risk * i_outcome) + ((1 - risk) * ui_outcome);
    i_outcome += costs;
    ui_outcome += costs;
    double ER = (risk * i_outcome) + ((1 - risk) * ui_outcome);

    double* retval = new double[3];
    retval[0] = EV;
    retval[1] = costs;
    retval[2] = ER;
    return retval;
}

/********************************************************
 * Beta Spread
 * Measures HLB Spread using beta regression results
 * *****************************************************/
double betaSpread(int relT, double efficacy, string strategy, double alpha) {
    // Coefficients
    double intercept = -2.44978330;
    double e75_coef = -0.20743718;
    double e85_coef = -0.07299259;
    double alpha1_coef = 0.01494436;
    double indv_coef = -0.15463509;
    double group_coef = -0.10212705;
    double maxT_coef = 0.00228886;
    double e75alpha1_coef = -0.06266322;
    double e85alpha1_coef = -0.79121443;
    double e75indv_coef = -0.16353347;
    double e85indv_coef = -0.64362937;
    double e75group_coef = -0.13874484;
    double e85group_coef = -0.85018937;
    double e75maxT_coef = 0.00018147;
    double e85maxT_coef = 0.00020873;
    double indvmaxT_coef = -0.00009667;
    double groupmaxT_coef = -0.00013745;
    double e75indvmaxT_coef = -0.00003253;
    double e85indvmaxT_coef = -0.00017724;
    double e75groupmaxT_coef = -0.00005853;
    double e85groupmaxT_coef = -0.00014676;
    
    //covariates
    int e75 = (efficacy == 0.75);
    int e85 = (efficacy == 0.85);
    int indv = (strategy == "Individual Action");
    int group = (strategy == "Group Action");

    //sum em up
    double sum_alpha1 = 0;
    double sum_alpha0 = 0;
    //sum the alpha0 case
    sum_alpha0 += intercept + 
                  (e75_coef * e75) + 
                  (e85_coef * e85) + 
                  (indv_coef * indv) + 
                  (group_coef * group) +
                  (maxT_coef * relT) + 
                  (e75indv_coef * (e75 * indv)) +
                  (e85indv_coef * (e85 * indv)) + 
                  (e75group_coef * (e75 * group)) +
                  (e85group_coef * (e85 * group)) +
                  (e75maxT_coef * (e75 * relT)) +
                  (e85maxT_coef * (e85 * relT)) +
                  (indvmaxT_coef * (indv * relT)) +
                  (groupmaxT_coef * (group * relT)) +
                  (e75indvmaxT_coef * (e75 * indv * relT)) +
                  (e85indvmaxT_coef * (e85 * indv * relT)) +
                  (e75groupmaxT_coef * (e75 * group * relT)) +
                  (e85groupmaxT_coef * (e85 * group * relT));
    // add the additional alpha terms
    sum_alpha1 = sum_alpha0 + 
                 alpha1_coef + 
                 (e75alpha1_coef * e75) +
                 (e85alpha1_coef * e85);
    //Take the weighted sum based on alpha
    double sum = sum_alpha0 + alpha*(sum_alpha1 - sum_alpha0);
    //Apply the logistic transformation 
    double transformed = 1 / (1 + exp(-1*sum));
    return transformed;
}

/********************************************************
* No Action : HLB Spread
* Returns the mean hlb spread t days after infection
********************************************************/
double NoAction::hlbSpread(int t, double efficacy, double alpha) {
    double res = betaSpread(t, efficacy, "No Action", alpha);
    return res;
}

/*******************************************************
* No Action : Plan Actions
* Fills a planning queue with empty actions
********************************************************/
void NoAction::PlanActions(vector<plan_func> * v) {
        plan_func doNothing = NULL;
        for (int i = 0; i < v->size(); i++) {
            (*v)[i] = doNothing;
        }
}

/********************************************************
 * Individual Action: HLB Spread
 * Returns the mean hlb spread t days after infection
*********************************************************/
double IndividualAction::hlbSpread(int t, double efficacy, double alpha) {
    double res = betaSpread(t, efficacy, "Individual Action", alpha);
    return res;
}

/*******************************************************
* Individual Action : Plan Actions
* CURRENTLY DEPRECATED PLANS TO REIMPLEMENT
********************************************************/
void IndividualAction::PlanActions(vector<plan_func> * v) {
        plan_func doNothing = NULL;
        for (int i = 0; i < v->size(); i++) {
            (*v)[i] = doNothing;
        }
}

/********************************************************
 * Individual Action: HLB Spread
 * Returns the mean hlb spread t days after infection
*********************************************************/
double GroupAction::hlbSpread(int t, double efficacy, double alpha) {
    double res = betaSpread(t, efficacy, "Group Action", alpha);
    return res;
}

/*******************************************************
* Group Action : Plan Actions
* CURRENTLY DEPRECATED PLANS TO REIMPLEMENT
********************************************************/
void GroupAction::PlanActions(vector<plan_func> * v) {
        plan_func doNothing = NULL;
        for (int i = 0; i < v->size(); i++) {
            (*v)[i] = doNothing;
        }
}














/*******************************************************
* Spray_K : HLB Spread
* Returns the mean HLB spread at time t
********************************************************/
double Spray_K::hlbSpread(int t) {
    string strategy; 
    double res;
    return res;
}

/********************************************************
* Spray_K : Plan Actions
* Fills the planning queue with evenly spaced spray events
********************************************************/
void Spray_K::PlanActions(vector<plan_func> * v) {
        plan_func spray = &sprayGrove;
        for (int t = 0; t < v->size(); t++) {
            if (t % (v->size() / numSprays) == 0) {
                (*v)[t] = spray;
            }
            else {
                (*v)[t] = NULL;
            }
        }
}


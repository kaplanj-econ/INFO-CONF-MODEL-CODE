#ifndef GROVE_HPP
#define GROVE_HPP
#include "commodity.hpp"
#include "planningFunc.hpp"
#include <vector>
using namespace std;

class Behavior;

class Grove {
private:
    Commodity crop; // The crop that is being grown
    bool agency; // 1: This grove makes decisions, 0: The grove has a fixed behavior pattern
    plan_func plannedActions[365]; //A vector of length planningLength that holds the planned actions for that period
    double fixedCosts; //Fixed cost per year associated with the grove
    Behavior* behaviorPattern;
    int ibounds[2]; // Lower inclusive, Upper Exclusive
    int jbounds[2]; // Lower inclusive, Upper Exclusive
    double lambda; // Trust in extension agent
    double alpha; // Expectation of neighbors coordination
    double sprayEfficacy;

public:
    double costs = 0;
    double returns = 0;
    double lastExtensionRisk = 0;
    double lastGrowerRisk = 0;
    double lastAdjustedRisk = 0;
    double lastRowPsyllids = 0;
    double lastColPsyllids = 0;
    double maxE_i = 0;
    double maxE_j = 0;
    double lastNAEV = 0;
    double lastISEV = 0;
    double lastGSEV = 0;
    bool foundHLB = false;
    int foundHLB_day = -1;
    double behaviorCosts[3];

    Grove();
    Grove(Commodity crop, bool agency, Behavior* behaviorPattern, int i_lb, int i_ub, int j_lb, int j_ub, double lambda, double alpha, double sprayE, int planningLength);
    //Grove();
    // Getters 
    Behavior* getBehavior() { return this->behaviorPattern; }
    bool hasAgency() { return this->agency; }
    plan_func* getPlanningQ() { return this->plannedActions; }
    Commodity* getCrop() { return &crop; }
    double getFixedCosts() { return this->fixedCosts; }
    int* getIBounds();
    int* getJBounds();
    //Setters
    void setAgency(bool);
    void setBehavior(Behavior* b) { this->behaviorPattern = b; }
    //Get a planned action
    plan_func getAction(int relativePeriod);
    double getLambda() { return this->lambda; }
    double getAlpha() { return this->alpha; }
    double getSprayEfficacy() { return this->sprayEfficacy; }
    void clearPlannedActions() {
        for (int i = 0; i < 365; i++) {
            plannedActions[i] = NULL;
        }
    }
    void setBehaviorCost(int idx, double value) {
        behaviorCosts[idx] = value;
    }
  
};

#endif
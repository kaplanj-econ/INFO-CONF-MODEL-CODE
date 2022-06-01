#ifndef BEHAVIOR_HPP
#define BEHAVIOR_HPP
#include "grove.hpp"
#include "math.h"
#include "planningFunc.hpp"
#include "bioABM.h"
#include<vector>
using namespace std;

void sprayGrove(Grove * g);
double getInfectedYield(double severity);
class Behavior {
public:
    //Fill up the planned actions vector for the next season
    virtual void PlanActions(vector<plan_func> * v) = 0;
    //Simulate outcomes of a grove under this behavior pattern
    double getExpectedValue(Grove g, double risk, int simulationLength, 
                            int startingPeriod, int planningLength,
                            double sprayEfficacy, double alpha, double additionalCosts);

    double* getExpectedValueTester(Grove g, double risk, int simulationLength, 
                            int startingPeriod, int planningLength,
                            double sprayEfficacy, double alpha);    
    //Get the expected mean infection at time t
    virtual double hlbSpread(int t, double efficacy, double alpha) = 0;

    //Returns the variable costs of behavior per planning period
    virtual double getVariableCosts() = 0;
};

class NoAction : public Behavior {
public:
    //Fills planning Q with NULLs
    void PlanActions(vector<plan_func> * v);

    //Returns expected value of this behavior pattern if continued until the end of the simulation
    double SimulateOutcome(Grove * g, double risk, int simulationLength,int startingPeriod);

    //Returns the expected mean infection with no mitigation
    double hlbSpread(int t, double efficacy, double alpha);

    //Returns the variable costs per planning period
    double getVariableCosts() { return 0; }
};

class IndividualAction: public Behavior {
    private:
        int windowSize;
        double vc;
    public:
        IndividualAction(int window, double sprayCost, int planningLength) { 
            windowSize = window; 
            this->vc = sprayCost *6 / 4; 
        }
        // //Fills planning Q with sprays in a window
        void PlanActions(vector<plan_func> * v);

        //Returns expected value of this behavior pattern if continued until the end of the simulation
        double SimulateOutcome(Grove * g, double risk, int simulationLength,int startingPeriod, double efficacy);

        //Returns the expected mean infection
        double hlbSpread(int t, double efficacy, double alpha);

        //Returns the variable costs per year
        double getVariableCosts() { return this->vc; }
};

class GroupAction: public Behavior {
 private:
    int windowSize;
    double vc;
 public:
        GroupAction(int window, double sprayCost, double groupCost, int planningLength) { 
            windowSize = window; 
            this->vc = (sprayCost + groupCost) * 6 / 4; }
        // //Fills planning Q with sprays in a window
        void PlanActions(vector<plan_func> * v);

        //Returns expected value of this behavior pattern if continued until the end of the simulation
        double SimulateOutcome(Grove * g, double risk, int simulationLength,int startingPeriod, double efficacy);

        //Returns the expected mean infection
        double hlbSpread(int t, double efficacy, double alpha);

        //Returns the variable costs per year
        double getVariableCosts() { return this->vc; }
};

class Spray_K : public NoAction {
public:
    static double sprayEfficacy;
    int numSprays;
    double sprayCost;

    int getNumSprays() { return numSprays; }

    void setNumSprays(int sprays) { numSprays = sprays; }

    void PlanActions(vector<plan_func>* v);

    double SimulateOutcome(Grove* g, double risk, int simulationLength, int startingPeriod, double yieldDecay);
    double hlbSpread(int t);

    double getVariableCosts() { return sprayCost * numSprays; }
};

class Individual_Spray : public Spray_K {

};

class Group_Spray : public Spray_K {

};
#endif

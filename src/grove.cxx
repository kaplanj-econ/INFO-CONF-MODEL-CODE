#include "../headers/grove.hpp"

/*
Grove::Grove() {
    //Please dont use this
    setAge(0);
}*/



Grove::Grove() {
}

Grove::Grove(Commodity crop, bool agency, Behavior * behaviorPattern, int i_lb, int i_ub, int j_lb, int j_ub, double lambda, double alpha, double sprayE, int planningLength) {
    this->crop = crop;
    setAgency(agency);
    this->behaviorPattern = behaviorPattern;
    ibounds[0] = i_lb;
    ibounds[1] = i_ub;
    jbounds[0] = j_lb;
    jbounds[1] = j_ub;
    this->lambda = lambda;
    this->alpha = alpha;
    this->sprayEfficacy = sprayE;
    this->clearPlannedActions();
    for (int i = 0; i < 3; i++) {
        behaviorCosts[i] = 0;
    }
}

int* Grove::getIBounds() {
    return this->ibounds;
}

int* Grove::getJBounds() {
    return this->jbounds;
}

//set agency
void Grove::setAgency(bool agency) {
    this->agency = agency;
}

//Get planned action at relative period
plan_func Grove::getAction(int relativePeriod) {
    plan_func p = NULL;
    if (this->getPlanningQ() == NULL) {
        return p;
    }
    if ( (relativePeriod >= 0) && (relativePeriod < 365) ) {
        p = this->plannedActions[relativePeriod];    
    }
    return p;
}

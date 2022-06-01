#include "../headers/bioABM.h"
using namespace std;
//using namespace concurrency;
namespace bioABM {

//#define _DEBUG

#ifdef _DEBUG
#define debug_log(s) cout << s << endl
#else 
#define debug_log(s)
#endif

//#include "../headers/tbb/concurrent_vector.h"
//tbb::concurrent_vector<int> testVector;
/***************************************************
 * TO-DO
 * 110 trees per acre gotten from orange extension
 * docs
 * *************************************************/
/****************************************************
 * PAPER PARAMETERS
 * **************************************************/
int maxFlushAge = 30; // Maximum Flush Age
int flushEmerging = 20; // Flush Shoots Emerging
int eggAdultTransition_25C = 17;
int eggAdultTransition_28C = 14;
int eggAdultTransition = 17; //Egg to Adult Transition
int durationYoungFlush_25C = 13;
int durationYoungFlush_28C = 16;
int durationYoungFlush = 13; //Duration of young flush
float proportionMigrating = 0.4; //Proportion of Migrating Adults
float withinRowP = 0.95; //Within-row probability
float betweenRowP = 0.5; //Between-row probability
int eggDuration_25C = 4;
int eggDuration_28C = 3;
int eggDuration = 4; //Egg Duration
int nymphDuration_25C = 13;
int nymphDuration_28C = 11;
int nymphDuration = 13; //Nymph Duration
int shootCapacity = 40; //Maximum shoots a flush can support
int shootEggCapacity = 40; // Daily flush shoot egg capacity
int eggsPerFemaleAdult = 10; //Eggs per Female Adult
float transmissionFlushNymph = 0.083; //Transmission from flushes to nymphs
float transmissionAdultFlush = 0.3; //Transmission from adults to flushes
int latentPeriod = 15; //Latent Period, T=15 (see Supplement 1)
float eggSurvivalP_25C = 0.8614;
float eggSurvivalP_28C = 0.8343;
float eggSurvivalP = 0.8614; // Egg/nymph survival probability
float adultSurvivalP_25C = 0.9847;
float adultSurvivalP_28C = 0.9659;
float adultSurvivalP = 0.9847; //Adult survival probability OG: 0.9847
int nymphMinAgeToInfect = eggAdultTransition; // Minimum age to transmit HLB to flush, paper says adult??
int nymphMinAgeToBeInfected_25C = 9;
int nymphMinAgeToBeInfected_28C = 8;
int nymphMinAgeToBeInfected = 9; //Minimum age to receive HLB from flush
int modelDuration = 366 * 2; //Duration of model in days
int springFlushStart = 80;
int springFlushEnd = 140;
int summerFlushStart = 180;
int summerFlushEnd = 210;
int fallFlushStart = 250;
int fallFlushEnd = 265;
int invasionDay = 80;
int carryingCapacity = 40000;
int invasionModality = 1;
int invasionGrove = 0;
double borderCrossingP = 0.01;

// Lattice dimensions
const int rowLength = 75;
const int numRows = 33;
const int hBorders = 2;
const int vBorders = 2;


// SQL
int experimentID;

void setExperimentID(int id) {
    experimentID = id;
}

/**********************************************************
 * OTHER
 * ********************************************************/
bool isFlushingPeriod = false;
bool modelStarted = false;
int modelDay = -1;

ofstream csvFile;
string csvName = "vs_output.csv";
boost::random::lagged_fibonacci607 rng(std::random_device{}());
boost::random::uniform_01<> gen;
struct FlushPatch;
typedef boost::tuple<int, int> coord;

int intRand(const int& min, const int& max) {
    static thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
}

double doubleRand(double min, double max) {
    static thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution(min, max);
    //return gen(rng);
    return distribution(generator);
}


struct Psyllid {
public:
    int age = 0;
    bool infected = false;
    bool female = false;
    coord pos = coord(-1,-1);
    bool alive = true;
    bool moved = false;
};

struct FlushShoot {
    int age = 0;
    bool infected = false;
    bool asymptomatic = false;
    int numEggs = 0;
    int daysAsymptomatic = 0;
    //atomic<int> numPsyllids_male;
    //atomic<int> numPsyllids_female;
    //atomic<int> numInfectedPsyllids_male;
    //atomic<int> numInfectedPsyllids_female;
    //array<int,17> numNymphs;
    //array<int,17> numInfectedNymphs;
    bool alive = true;
    bool bark = false;
    int numPsyllids = 0;
    int numInfectedPsyllids = 0;
    vector<Psyllid*> psyllids;

    
    FlushShoot() {
        /*numPsyllids_male.store(0);
        numPsyllids_female.store(0);
        numInfectedPsyllids_female.store(0);
        numInfectedPsyllids_male.store(0);*/
    }

    FlushShoot operator=(const FlushShoot&) {
        
        FlushShoot fs;
        /*
        fs.numPsyllids_male.store(numPsyllids_male.load());
        fs.numPsyllids_female.store(numPsyllids_female.load());
        fs.numInfectedPsyllids_male.store(numInfectedPsyllids_male.load());
        fs.numInfectedPsyllids_female.store(numInfectedPsyllids_female.load());
        fs.bark = bark;
        fs.infected = infected;
        fs.asymptomatic = asymptomatic;
        fs.daysAsymptomatic = daysAsymptomatic;
        fs.age = age;*/
        return fs;
    }

    FlushShoot(const FlushShoot& f) {
        /*
        numPsyllids_male.store(f.numPsyllids_male.load());
        numPsyllids_female.store(f.numPsyllids_female.load());
        numInfectedPsyllids_male.store(f.numInfectedPsyllids_male.load());
        numInfectedPsyllids_female.store(f.numInfectedPsyllids_female.load());
        bark = f.bark;
        infected = f.infected;
        asymptomatic = f.asymptomatic;
        daysAsymptomatic = f.daysAsymptomatic;
        age = f.age;*/
        return;
    }
    
    void placePsyllid(Psyllid * p) {
        if (p->infected) {
            numInfectedPsyllids++;
        }
        else {
            numPsyllids++;
        }


    }

    /*
    void placePsyllid(bool female, bool adult, bool infected, int age = -1) {
        if (female && adult && infected) {
            numInfectedPsyllids_female++;
        }
        else if (female && adult && !infected) {
            numPsyllids_female++;
        }
        else if (!female && adult && infected) {
            numInfectedPsyllids_male++;
        }
        else if (!female && adult && !infected) {
            numPsyllids_male++;
        }
        else if (!adult && infected) {
            numInfectedNymphs[age]++;
        }
        else if (!adult && !infected) {
            numNymphs[age]++;
        }
    }*/

    

    //Deprecated
    void kill() {
        auto it = psyllids.begin();
        while (it != psyllids.end()) {
            delete (*it);
            it = psyllids.erase(it);
        }
    }

    void increaseAge() {
        age += 1;
        if (!bark) {
            if (asymptomatic) {
                daysAsymptomatic += 1;
                if (daysAsymptomatic == latentPeriod) {
                    infected = true;
                    asymptomatic = false;
                }
            }
        }
    }
};

struct FlushPatch {
    int age = 0;

    int oldInfectedShoots = 0;
    int oldUninfectedShoots = 0;
    int numPsyllids_male = 0;
    int numPsyllids_female = 0;
    int numInfectedPsyllids_male = 0;
    int numInfectedPsyllids_female = 0;
    array<int,17> numNymphs;
    array<int,17> numInfectedNymphs;
    array<int, 30> numShoots;
    array<int, 30> numInfectedShoots;

    FlushPatch() {      
        for (int i = 0; i < 17; i++) {
            numNymphs[i] = 0;
            numInfectedNymphs[i] = 0;
        }
        for (int i = 0; i < 30; i++) {
            numShoots[i] = 0;
            numInfectedShoots[i] = 0;
        }
    }

    bool validate() {
        for (int i = 0; i < 17; i++) {
            if (numNymphs[i] < 0 || numInfectedNymphs[i] < 0) {
                return false;
            }
        }
        for (int i = 0; i < 30; i++) {
            if (numShoots[i] < 0 || numInfectedShoots[i] < 0) {
                return false;
            }
        }
        return true;
    }
    int getAge() {
        return age;
    }

    int getNumPsyllids() {
        int numPsyllids = 0;
        numPsyllids += numPsyllids_male;
        numPsyllids += numPsyllids_female;
        numPsyllids += numInfectedPsyllids_male;
        numPsyllids += numInfectedPsyllids_female;
        numPsyllids += accumulate(numNymphs.begin(), numNymphs.end(), 0);
        numPsyllids += accumulate(numInfectedNymphs.begin(), numInfectedNymphs.end(), 0);
        return numPsyllids;
    }


    double getHLBSeverity() {
        int uninfected = accumulate(numShoots.begin(), numShoots.end(), 0);
        int infected = accumulate(numInfectedShoots.begin(), numInfectedShoots.end(), 0);
        double hlbNum = (double)infected + (double)oldInfectedShoots;
        double hlbDenom = (double)uninfected + (double)oldUninfectedShoots + (double)hlbNum;
       
        if (hlbDenom == 0) {
            return 0;
        }
        else {
            assert((hlbNum / hlbDenom) >= 0 && (hlbNum / hlbDenom) <= 1);
            double severity = hlbNum / hlbDenom;
            return severity;
        }

    }

    int getTotalPsyllids() {
        return getNumPsyllids();
    }

    void placePsyllid(bool female, bool adult, bool infected, int age = -1) {
        if (adult && female && infected) {
            numInfectedPsyllids_female++;
        }
        else if (adult && female && !infected) {
            numPsyllids_female++;
        }
        else if (adult && !female && infected) {
            numInfectedPsyllids_male++;
        }
        else if (adult && !female && !infected) {
            numPsyllids_male++;
        }
        else if (!adult && !infected) {
            numNymphs[age]++;
        }
        else if (!adult && infected) {
            numInfectedNymphs[age]++;
        }
    }
};



typedef array<FlushPatch, rowLength> FlushRow;

// Connect signal to mark dead

array<FlushRow, numRows> lattice; //11x25 lattice that represents flush patches
vector<Psyllid*> homelessPsyllids;

/****************************************************
 * CUSTOM DATA STRUCTURES
 * **************************************************/


enum PositionType {MIDDLE, EDGE, CORNER};


/********************************************************
 * MODEL PARAMETERS
 * ******************************************************/
float initialInfectedPortion = 0.18;
int initialNumPsyllids = 300;
typedef boost::tuple<vector<coord>, vector<double>> TransitionMap;
//Transition possibilities for "Middle" cells
vector<coord> middleCellDifferentials = boost::assign::list_of(coord(0,1))(coord(0,-1))(coord(1,0))(coord(-1,0));
vector<double> middleCellProbabilities = boost::assign::list_of(0.45)(0.45)(0.05)(0.05);
TransitionMap middleCellTransitions(
    middleCellDifferentials,
    middleCellProbabilities
);
//Transition possibilities for "Edge" cells
vector<coord> edgeCellDifferentials = boost::assign::list_of(coord(0,1))(coord(0,-1))(coord(1,0))(coord(-1,0));
vector<double> edgeCellProbabilities = boost::assign::list_of(0.45)(0.45)(0.05)(0.05);
TransitionMap edgeCellTransitions(
    edgeCellDifferentials,
    edgeCellProbabilities
);
//Transition possibilities for "Corner" cells
vector<coord> cornerCellDifferentials = boost::assign::list_of(coord(0,1))(coord(1,0));
vector<double> cornerCellProbabilities = boost::assign::list_of(0.9)(0.1);
TransitionMap cornerCellTransitions(
    cornerCellDifferentials,
    cornerCellProbabilities
);
//Setup Map from PositionType to it's associated transition map, MAKE SURE IT MATCHES ENUM DECLARATION
array<TransitionMap, 3> positionProbabilityMap{
    middleCellTransitions,
    edgeCellTransitions,
    cornerCellTransitions
};


/************************************************************
* Get Flushing Bounds
************************************************************/
int getFallStart() {
    return fallFlushStart;
}
int getSpringStart() {
    return springFlushStart;
}
int getSummerStart() {
    return summerFlushStart;
}
/*************************************************************
* Get Age at
* Returns flush patch age at coordinates
*************************************************************/
int getAgeAt(int i, int j, int differential) {
    return (lattice[i][j].getAge() + differential) / 365;
}


/************************************************************
* Get Severity At
* Returns flush patch severity at coordinates
*************************************************************/
double getSeverityAt(int i, int j) {
    return lattice[i][j].getHLBSeverity();
}

/**************************************************************
* Get Psyllids At
* Returns number of psyllids at coordinates
**************************************************************/
double getPsyllidsAt(int i, int j) {
    if (i > numRows || j > rowLength || i < 0 || j < 0) {
        return 0;
    }
    else {
        return lattice[i][j].getTotalPsyllids();
    }
}

/*************************************************************
 * Update Seasonal Parameters
 * Updates parameters that vary based on operative temperature
 * ***********************************************************/
void updateSeasonalParameters(int relativeT) {
    if ( (relativeT >= springFlushStart && relativeT <= springFlushEnd) ||
         (relativeT >= fallFlushStart && relativeT <= fallFlushEnd) ) {
             eggAdultTransition = eggAdultTransition_25C;
             durationYoungFlush = durationYoungFlush_25C;
             eggDuration = eggDuration_25C;
             nymphDuration = nymphDuration_25C;
             eggSurvivalP = eggSurvivalP_25C;
             adultSurvivalP = adultSurvivalP_25C;
         }
    else if (relativeT >= summerFlushStart && relativeT <= summerFlushEnd) {
        eggAdultTransition = eggAdultTransition_28C;
        durationYoungFlush = durationYoungFlush_28C;
        eggDuration = eggDuration_28C;
        nymphDuration = nymphDuration_28C;
        eggSurvivalP = eggSurvivalP_28C;
        adultSurvivalP = adultSurvivalP_28C;
    }
}

/**************************************************************
* Get Model Day
* For econ wrapper
***************************************************************/
int getModelDay() {
    return modelDay;
}

/***************************************************************
* Get Model Duration
* For econ wrapper
***************************************************************/
int getModelDuration() {
    return modelDuration;
}

/**************************************************************
 * Is Flushing Period
 * ************************************************************/
void setFlushingPeriod(int period) {
    int t = period % 365;
    if ( (t >= springFlushStart && t <= springFlushEnd) ||
         (t >= summerFlushStart && t <= summerFlushEnd) || 
         (t >= fallFlushStart   && t <= fallFlushEnd)) {
        isFlushingPeriod = true;
        updateSeasonalParameters(t);
    }
    else {
        isFlushingPeriod = false;
    }
}

/*****************************************************************
* rogueTreeAt
******************************************************************/

void rogueTreeAt(coord pos) {
    /*
    //Create a new tree
    FlushPatch fp;
    FlushShoot bark;
    bark.bark = true;
    fp.shoots.push_back(bark);

    //Remove bugs on desired patch
    FlushPatch* dying = &lattice[pos.get<0>()][pos.get<1>()];
    auto it = dying->shoots.begin();
    while (it != dying->shoots.end()) {
        auto pIt = (*it).psyllids.begin();
        while (pIt != (*it).psyllids.end()) {
            delete (*pIt);
            pIt = (*it).psyllids.erase(pIt);
        }
    }

    //Replace patch
    lattice[pos.get<0>()][pos.get<1>()] = fp;*/
    return;
}

/****************************************************************
* Spray trees
******************************************************************/

void sprayTrees(double efficacy, vector<coord> locations) {
    int psyllidsRemoved = 0;
    for (int i = 0; i < locations.size(); i++) {
        int r = locations[i].get<0>();
        int c = locations[i].get<1>();
        int beforePsyllids = lattice[r][c].getTotalPsyllids();
        /*lattice[r][c].numPsyllids_male = ceil((1.0 - efficacy) * (double)lattice[r][c].numPsyllids_male);
        lattice[r][c].numPsyllids_female = ceil((1.0 - efficacy) * (double)lattice[r][c].numPsyllids_female);
        lattice[r][c].numInfectedPsyllids_male = ceil((1.0 - efficacy) * (double)lattice[r][c].numInfectedPsyllids_male);
        lattice[r][c].numInfectedPsyllids_female = ceil((1.0 - efficacy) * (double)lattice[r][c].numInfectedPsyllids_female);*/
        int startingMales = lattice[r][c].numPsyllids_male;
        for (int j = 0; j < startingMales; j++) {
            if (doubleRand(0, 1) <= efficacy) {
                lattice[r][c].numPsyllids_male--;
            }
        }
        int startingFemales = lattice[r][c].numPsyllids_female;
        for (int j = 0; j < startingFemales; j++) {
            if (doubleRand(0, 1) <= efficacy) {
                lattice[r][c].numPsyllids_female--;
            }
        }
        int startingMales_i = lattice[r][c].numInfectedPsyllids_male;
        for (int j = 0; j < startingMales_i; j++) {
            if (doubleRand(0, 1) <= efficacy) {
                lattice[r][c].numInfectedPsyllids_male--;
            }
        }
        int startingFemales_i = lattice[r][c].numInfectedPsyllids_female;
        for (int j = 0; j < lattice[r][c].numInfectedPsyllids_female; j++) {
            if (doubleRand(0, 1) <= efficacy) {
                lattice[r][c].numInfectedPsyllids_female--;
            }
        }
        for (int j = 0; j < 17; j++) {
            int startingNymphs = lattice[r][c].numNymphs[j];
            for (int k = 0; k < startingNymphs; k++) {
                if (doubleRand(0, 1) <= efficacy) {
                    lattice[r][c].numNymphs[j]--;
                }
            }
            int startingNymphs_i = lattice[r][c].numInfectedNymphs[j];
            for (int k = 0; k < startingNymphs_i; k++) {
                if (doubleRand(0, 1) <= efficacy) {
                    lattice[r][c].numInfectedNymphs[j]--;
                }
            }
            //lattice[r][c].numNymphs[j] = ceil((1.0 - efficacy) * (double)lattice[r][c].numNymphs[j]);
            //lattice[r][c].numInfectedNymphs[j] = ceil((1.0 - efficacy) * (double)lattice[r][c].numInfectedNymphs[j]);
        }
        int afterPsyllids = lattice[r][c].getTotalPsyllids();
        psyllidsRemoved += beforePsyllids - afterPsyllids;
    }
    //cout << "Spray removed " << psyllidsRemoved << endl;
}

/**************************************************************
 * crossesBorder
 * Determines if movement between the two cells is a border 
 * crossing
 * ************************************************************/
bool crossesBorder(coord a, coord b) {
    int aQuadrant_i = a.get<0>() / (numRows / (hBorders+1));
    int aQuadrant_j = a.get<1>() / (rowLength / (vBorders+1));
    int bQuadrant_i = b.get<0>() / (numRows / (hBorders+1));
    int bQuadrant_j = b.get<1>() / (rowLength / (vBorders+1));
    if (aQuadrant_i != bQuadrant_i || aQuadrant_j != bQuadrant_j) {
        return true;
    }
    else {
        return false;
    }
}

/***************************************************************
 * Is Valid Coordinate
 * ************************************************************/
bool isValidCoordinate(coord pos) {
    int row = pos.get<0>();
    int col = pos.get<1>();
    if (row >= 0 && row < numRows && col >= 0 && col < rowLength) {
        return true;
    }
    else {
        return false;
    }
}

/**************************************************************
 * Determine position type
 * Determines if a coordinate is a middle, edge, or corner
 * ************************************************************/
PositionType determinePositionType(coord pos) {
    int row = pos.get<0>();
    int col = pos.get<1>();
    array<coord, 4> neighbors = { coord(row,col+1), coord(row, col-1), coord(row+1, col), coord(row-1, col)};
    int numNeighbors = 0;
    for (int i = 0; i < neighbors.size(); i ++) {
        if (isValidCoordinate(neighbors[i])) {
            numNeighbors += 1;
        }
    }
    if (numNeighbors < 2 || numNeighbors > 4) {
        debug_log("Error in determine position type");
        assert(numNeighbors >= 2 && numNeighbors <= 4);
    }
    switch(numNeighbors) {
        case 2:
            return PositionType::CORNER;
        case 3:
            return PositionType::EDGE;
        case 4:
            return PositionType::MIDDLE;
        default:
            assert(false); //ERROR
            return PositionType::MIDDLE;
    }
}



/**************************************************************
 * Discrete Probability Match
 * Takes in a discrete vector of probabilities, and chooses the
 * corresponding index based on a uniform pull 
 * and those probabilities
 * ***********************************************************/
int discreteProbabilityMatch(vector<double> probabilities) {
    double pull;
    pull = doubleRand(0, 1);
    double cumSum = 0;
    int resultIdx = -1;
    for (int i = 0; i < probabilities.size(); i++) {
        if (pull >= cumSum && pull <= (cumSum + probabilities[i])) {
            resultIdx = i;
            break;
        }
        else {
            cumSum += probabilities[i];
        }
    }
    if (resultIdx < 0) {
        if (pull == 1) {
            resultIdx = probabilities.size() - 1;
        }
        else {
            cout << "Well there's ya problem\n";
            assert(resultIdx >= 0);
        }

    }
    return resultIdx;
}




/*************************************************
 * Initialize Lattice
 * ***********************************************/
void initializeLattice() {
    for (int i = 0; i < numRows; i++) {
        FlushRow row;
        for (int j = 0; j < rowLength; j++) {
            FlushPatch patch;
            row[j] = patch;
        }
        lattice[i] = row;
    }
}




/***************************************************
* getGroveBounds
* Returns the row and column upper/lower bounds
* based on a two digit identifier, indices are
* 0: row bounds
* 1: column bounds
* Bounds are lower inclusive, upper exclusive
****************************************************/
vector<coord> getGroveBounds(int identifier) {
    int rowID = identifier / 10;
    int colID = identifier % 10;
    
    //First row bounds
    int gamma_r = numRows / (hBorders + 1);
    int rowLB = gamma_r * rowID;
    int rowUB = gamma_r * (rowID + 1);
    coord rowBounds(rowLB, rowUB);

    //Same for column
    int gamma_c = rowLength / (vBorders + 1);
    int colLB = gamma_c * colID;
    int colUB = gamma_c * (colID + 1);
    coord colBounds(colLB, colUB);

    //Package and return
    vector<coord> bounds;
    bounds.push_back(rowBounds);
    bounds.push_back(colBounds);
    return bounds;
}
/***************************************************
* uniformPsyllidDistribution
* Distributes a number of uninfected psyllids evenly 
* amongst cells not in the specified vector
* *************************************************/
void uniformPsyllidDistribution(double percent, int numPsyllids, vector<coord> occupied, int groveID) {
    //Determine number of trees to infect
    int rowsPerGrove = numRows / (hBorders + 1);
    int rowLPerGrove = rowLength / (vBorders + 1);
    int totalTrees = rowsPerGrove * rowLPerGrove;
    int availableTrees = totalTrees - occupied.size();
    int treesToInfect = floor((double)availableTrees * percent);

    //Set up coordinates available
    vector<coord> groveTrees = getGroveBounds(groveID);
    int psyllidsPerTree = floor((double)numPsyllids / (double)treesToInfect);

    //For using std::find
    
    vector<coord> cells;
    for (int i = 0; i < treesToInfect; i++) {
        //Choose a tree
        coord c;
        c = coord(intRand(groveTrees[0].get<0>(), groveTrees[0].get<1>() - 1), intRand(groveTrees[1].get<0>(), groveTrees[1].get<1>() - 1));
        //Reroll until it's not one of the off-limit trees or one we already picked
        while (find(occupied.begin(), occupied.end(), c) != occupied.end() ||
               find(cells.begin(), cells.end(), c) != cells.end()) {
            c = coord(intRand(groveTrees[0].get<0>(), groveTrees[0].get<1>() - 1), intRand(groveTrees[1].get<0>(), groveTrees[1].get<1>() - 1));
        }
        cells.push_back(c);
  
    }
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < psyllidsPerTree; j++) {
            bool infected = doubleRand(0, 1) < initialInfectedPortion;
            bool female = doubleRand(0, 1) < 0.5;
            lattice[cells[i].get<0>()][cells[i].get<1>()].placePsyllid(female, true, infected);
        }
    }
}


/***************************************************
* invasionModality1
* 200 psyllids, evenly distributed 
* amongst 4 trees in the SW corner
****************************************************/
vector<coord> invasionModality1(int groveID) {
    int psyllidsPerTree = 50;
    vector<coord> bounds = getGroveBounds(groveID);
    vector<coord> cells;

    //Collect the 4 trees in the SW corner
    cells.push_back(coord(bounds[0].get<1>() - 1, bounds[1].get<0>()));
    cells.push_back(coord(bounds[0].get<1>() - 2, bounds[1].get<0>()));
    cells.push_back(coord(bounds[0].get<1>() - 1, bounds[1].get<0>() + 1));
    cells.push_back(coord(bounds[0].get<1>() - 2, bounds[1].get<0>() + 1));

    //Place psyllids
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < psyllidsPerTree; j++) {
            bool infected = doubleRand(0,1) < initialInfectedPortion;
            bool female = doubleRand(0,1) < 0.5;
            lattice[cells[i].get<0>()][cells[i].get<1>()].placePsyllid(female, true, infected);
        }
    }
    //cout << "Psyllids placed: " << numPsyllids << endl;
    return cells;

}

/***************************************************************
* invasionModality2
*   Modality 1 AND 35% of randomly selected trees from
*   remaining patches are occupied by 200 uninfected 
*   ACP distributed evenly
****************************************************************/
void invasionModality2(int groveID) {
    vector<coord> cells = invasionModality1(groveID);
    uniformPsyllidDistribution(0.35, 200, cells, groveID);
}

/****************************************************************
* invasionModality3
* 200 psyllids, placed on 25% of trees on the southern edge
* and 100% trees on the eastern edge
****************************************************************/
vector<coord> invasionModality3(int groveID) {
    vector<coord> cells;
    //Get bounds
    vector<coord> bounds = getGroveBounds(groveID);

    //Determine number of psyllids per tree
    int southernTrees = bounds[1].get<1>() - bounds[1].get<0>();
    southernTrees = floor((double)southernTrees * 0.25);
    int easternTrees = bounds[0].get<1>() - bounds[0].get<0>();
    int totalTrees = southernTrees + easternTrees;
    int psyllidsPerTree = floor(200 / (double)totalTrees);

    //Add eastern trees
    for (int i = bounds[0].get<0>(); i < bounds[0].get<1>(); i++) {
        coord c(i, bounds[1].get<1>() - 1);
        cells.push_back(c);
    }

    //Add southern trees
    for (int i = 0; i < southernTrees; i++) {
        coord c(bounds[0].get<1>() - 1, intRand(bounds[1].get<0>(), bounds[1].get<1>() - 2));
        while (find(cells.begin(), cells.end(), c) != cells.end()) {
            c = coord(bounds[0].get<1>() - 1, intRand(bounds[1].get<0>(), bounds[1].get<1>() - 2));
        }
        cells.push_back(c);
    }

    //Distribute psyllids
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < psyllidsPerTree; j++) {
            bool infected = doubleRand(0, 1) < initialInfectedPortion;
            bool female = doubleRand(0, 1) < 0.5;
            lattice[cells[i].get<0>()][cells[i].get<1>()].placePsyllid(female, true, infected);
        }
    }
    return cells;
}

/***************************************************************
* invasionModality4
* Modality 3 AND 35% of randomly selected trees from 
* remaining patches are occupied by 200 uninfected
* ACP distributed evenly
* **************************************************************/
void invasionModality4(int groveID) {
    vector<coord> cells = invasionModality3(groveID);
    uniformPsyllidDistribution(0.35, 200, cells, groveID);
}

/****************************************************************
* invasionModality5
* 10 trees distributed around the center of the grove
* are occupied by 200 psyllids
* ***************************************************************/
vector<coord> invasionModality5(int groveID) {
    //Collect the ten middle cells
    vector<coord> cells;
    vector<coord> bounds = getGroveBounds(groveID);
    int center_i = ceil(  ( (double)bounds[0].get<1>() - 1 - (double)bounds[0].get<0>() ) / 2) + bounds[0].get<0>();
    int center_j = ceil(  ( (double)bounds[1].get<1>() - 1 - (double)bounds[1].get<0>() ) / 2) + bounds[1].get<0>();
    for (int i = center_i - 1; i <= center_i + 1; i++) {
        for (int j = center_j - 1; j <= center_j + 1; j++) {
            coord c(i, j);
            cells.push_back(c);
        }
    }
    //Distribute psyllids
    int psyllidsPerTree = 20;
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < psyllidsPerTree; j++) {
            bool infected = doubleRand(0, 1) < initialInfectedPortion;
            bool female = doubleRand(0, 1) < 0.5;
            lattice[cells[i].get<0>()][cells[i].get<1>()].placePsyllid(female, true, infected);
        }
    }

    return cells;

}

/***************************************************************
* invasionModality6
* Modality 5 AND 35% of randomly selected trees from
* remaining patches are occupied by 200 uninfected 
* ACP distributed evenly
****************************************************************/
void invasionModality6(int groveID) {
    vector<coord> cells = invasionModality5(groveID);
    uniformPsyllidDistribution(0.35, 200, cells, groveID);
}


/***************************************************************
 * Place Initial Psyllids
 * INVASION SPATIAL MODALITIES:
 *      1: 200 psyllids, evenly distributed 
 *          amongst 4 trees in the SW corner
 *      2: Modality 1 AND 35% of randomly selected trees from
 *         remaining patches are occupied by 200 uninfected 
           ACP distributed evenly
 *      3: 200 psyllids, placed on 25% of trees on the southern edge
 *         and 100% of trees on the eastern edge
 *      4: Modality 3 AND 35% of randomly selected trees from 
 *         remaining patches are occupied by 200 uninfected
 *         ACP distributed evenly
 *      5: 10 trees distributed around the center of the grove
 *         are occupied by 200 psyllids
 *      6: Modality 5 AND 35% of randomly selected trees from
 *         remaining patches are occupied by 200 uninfected 
 *         ACP distributed evenly
 * ************************************************************/
void placeInitialPsyllids(int invasionModality, int groveID) {
    switch (invasionModality) {
        case 1:
            invasionModality1(groveID);
            break;
        case 2:
            invasionModality2(groveID);
            break;
        case 3:
            invasionModality3(groveID);
            break;
        case 4:
            invasionModality4(groveID);
            break;
        case 5:
            invasionModality5(groveID);
            break;
        case 6:
            invasionModality6(groveID);
            break;
        default:
            invasionModality1(groveID);
    }
    return;
}




/************************************************
 * Initialize Model
 * Desc
 * **********************************************/
void initializeModel() {
    //cout << "Initializing lattice...";
    initializeLattice();
    //cout << "Plaincg initial psyllids...\n";
    //placeInitialPsyllids();
    return;
}

/************************************************
 * Birth New Flush
 * Emergence of new flush: activity 1
 * - For each flush patch, birth 20 shoots
 * **********************************************/
void birthNewFlush() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            double birthInfectChance = lattice[i][j].getHLBSeverity();
            //assert(birthInfectChance < 1);
            if (birthInfectChance > 0) {
                for (int k = 0; k < flushEmerging; k++) {
                    if (doubleRand(0, 1) <= birthInfectChance) {
                        lattice[i][j].numInfectedShoots[0]++;
                    }
                    else {
                        lattice[i][j].numShoots[0]++;
                    }
                }
            }
            else {
                lattice[i][j].numShoots[0] += flushEmerging;
            }
        }
    }
}

void migratePsyllidAt(coord pos, bool infected, bool female) {
    
    PositionType posType = determinePositionType(pos);
    TransitionMap transitionMap = positionProbabilityMap[posType];
    int stateIdx = discreteProbabilityMatch(transitionMap.get<1>());
    
    
    //Associated position differential
    coord positionDiff = transitionMap.get<0>()[stateIdx];
    coord destination = coord(positionDiff.get<0>() + pos.get<0>(), positionDiff.get<1>() + pos.get<1>());

    bool borderCrossing = false;
    bool successfulCrossing = false;
    if (hBorders > 0 || vBorders > 0) {
        if (crossesBorder(pos, destination)) {
            double bCrossPull = doubleRand(0, 1);
            bCrossPull = doubleRand(0, 1);
            borderCrossing = true;
            successfulCrossing = (bCrossPull <= borderCrossingP);
        }
    }
    if (!isValidCoordinate(destination) || (borderCrossing && !successfulCrossing)) {
        destination.get<0>() = -1 * positionDiff.get<0>() + pos.get<0>();
        destination.get<1>() = -1 * positionDiff.get<1>() + pos.get<1>();
    }
    //Just in case it's not
    assert(isValidCoordinate(destination));
    
    //Removal accounting
    if (infected && female) {
        lattice[pos.get<0>()][pos.get<1>()].numInfectedPsyllids_female--;
    }
    else if (infected && !female) {
        lattice[pos.get<0>()][pos.get<1>()].numInfectedPsyllids_male--;
    }
    else if (!infected && female) {
        lattice[pos.get<0>()][pos.get<1>()].numPsyllids_female--;
    }
    else if (!infected && !female) {
        lattice[pos.get<0>()][pos.get<1>()].numPsyllids_male--;
    }

    //Add psyllid at new location
    lattice[destination.get<0>()][destination.get<1>()].placePsyllid(female, true, infected);

}

void migration_parallel() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            int malesMigrating = floor(lattice[i][j].numPsyllids_male * proportionMigrating);
            int femalesMigrating = floor(lattice[i][j].numPsyllids_female * proportionMigrating);
            int malesMigrating_i = floor(lattice[i][j].numInfectedPsyllids_male * proportionMigrating);
            int femalesMigrating_i = floor(lattice[i][j].numInfectedPsyllids_female * proportionMigrating);

            for (int l = 0; l < malesMigrating; l++) {
                migratePsyllidAt(coord(i, j), false, false);
            }
            //#pragma omp parallel for
            for (int l = 0; l < femalesMigrating; l++) {
                migratePsyllidAt(coord(i, j), false, true);
            }
            //#pragma omp parallel for
            for (int l = 0; l < malesMigrating_i; l++) {
                migratePsyllidAt(coord(i, j), true, false);
            }
            //#pragma omp parallel for
            for (int l = 0; l < femalesMigrating_i; l++) {
                migratePsyllidAt(coord(i, j), true, true);
            }
        }
    }
}

void psyllidAging_parallel() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            //nymphs
            double modifier = 1;
            if (lattice[i][j].getTotalPsyllids() > carryingCapacity) {
                modifier = carryingCapacity / (double)lattice[i][j].getTotalPsyllids();
            }
            double nymphSurvivalChance = eggSurvivalP * modifier;
            double adultSurvivalChance = adultSurvivalP * modifier;
            //Age nymphs
            for (int k = 16; k >= 0; k--) {
                if (k == 16) {
                    int startingNymphs_ui = lattice[i][j].numNymphs[k];
                    for (int l = 0; l < startingNymphs_ui; l++) {
                        if (doubleRand(0, 1) <= nymphSurvivalChance) {
                            if (doubleRand(0, 1) <= 0.5) {
                                lattice[i][j].numPsyllids_female++;
                            }
                            else {
                                lattice[i][j].numPsyllids_male++;
                            }
                        }
                    }
                    int startingNymphs_i = lattice[i][j].numInfectedNymphs[k];
                    for (int l = 0; l < startingNymphs_i; l++) {
                        if (doubleRand(0, 1) <= nymphSurvivalChance) {
                            if (doubleRand(0, 1) <= 0.5) {
                                lattice[i][j].numInfectedPsyllids_female++;
                            }
                            else {
                                lattice[i][j].numInfectedPsyllids_male++;
                            }
                        }
                    }
                }
                
                if (k == 0) {
                    lattice[i][j].numNymphs[k] = 0;
                    lattice[i][j].numInfectedNymphs[k] = 0;
                }
                else {
                    lattice[i][j].numNymphs[k] = 0;
                    lattice[i][j].numInfectedNymphs[k] = 0;
                    int startingNymphs_ui = lattice[i][j].numNymphs[k-1];
                    for (int l = 0; l < startingNymphs_ui; l++) {
                        if (doubleRand(0, 1) < nymphSurvivalChance) {
                            lattice[i][j].numNymphs[k]++;
                        }
                    }
                    int startingNymphs_i = lattice[i][j].numInfectedNymphs[k - 1];
                    for (int l = 0; l < startingNymphs_i; l++) {
                        if (doubleRand(0, 1) < nymphSurvivalChance) {
                            lattice[i][j].numInfectedNymphs[k]++;
                        }
                    }
                }
            }
            
           // adults
            int startingMales = lattice[i][j].numPsyllids_male;
            for (int k = 0; k < startingMales; k++) {
                if (doubleRand(0, 1) > adultSurvivalChance) {
                    lattice[i][j].numPsyllids_male--;
                }
            }
            int startingFemales = lattice[i][j].numPsyllids_female;
            for (int k = 0; k < startingFemales; k++) {
                if (doubleRand(0, 1) > adultSurvivalChance) {
                    lattice[i][j].numPsyllids_female--;
                }
            }
            int startingMales_i = lattice[i][j].numInfectedPsyllids_male;
            for (int k = 0; k < startingMales_i; k++) {
                if (doubleRand(0, 1) > adultSurvivalChance) {
                    lattice[i][j].numInfectedPsyllids_male--;
                }
            }
            int startingFemales_i = lattice[i][j].numInfectedPsyllids_female;
            for (int k = 0; k < startingFemales_i; k++) {
                if (doubleRand(0, 1) > adultSurvivalChance) {
                    lattice[i][j].numInfectedPsyllids_female--;
                }
            }
          
        }
    }
}

/*****************************************************
 * Egg management
 * Egg laying and egg survival: Activity 4 originally,
 * Activity 3 in this implementation. Since the model
 * treats eggs and nymphs interchangebly, egg/nymph
 * mortality is handled in the aging function
 * ****************************************************/

void eggManagement_parallel() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            int numMothers = lattice[i][j].numPsyllids_female + lattice[i][j].numInfectedPsyllids_female;
            int numViableShoots = accumulate(lattice[i][j].numShoots.begin(), lattice[i][j].numShoots.end(), 0);
            numViableShoots += accumulate(lattice[i][j].numInfectedShoots.begin(), lattice[i][j].numInfectedShoots.end(), 0);
            int totalNymphs = min(numMothers * eggsPerFemaleAdult, numViableShoots * shootEggCapacity);
            lattice[i][j].numNymphs[0] += totalNymphs;
        }
    }
}

/**************************************************
 * Disease transmission
 * Transmission between ACP and flush. A combination
 * of activities 5 and 6. First from f->p, then 
 * from p->f
 * ************************************************/
void diseaseTransmission_parallel() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {        
            //f->p
            int numUninfectedShoots = accumulate(lattice[i][j].numShoots.begin(), lattice[i][j].numShoots.end(), 0);
            int numInfectedShoots = accumulate(lattice[i][j].numInfectedShoots.begin(), lattice[i][j].numInfectedShoots.end(), 0);
            if (numInfectedShoots != 0) {
                double infectedproportion = (double)numInfectedShoots / ((double)numInfectedShoots + (double)numUninfectedShoots);
                for (int k = nymphMinAgeToBeInfected; k < 17; k++) {
                    int ui_nymphs = lattice[i][j].numNymphs[k];
                    ui_nymphs = floor(ui_nymphs * infectedproportion);
                    int nymphsInfected = floor(ui_nymphs * transmissionFlushNymph);
                    lattice[i][j].numNymphs[k] -= nymphsInfected;
                    lattice[i][j].numInfectedNymphs[k] += nymphsInfected;                 
                }
            }
         
            //p->f
            if (numUninfectedShoots != 0) {
                int numInfectedPsyllids = lattice[i][j].numInfectedPsyllids_female + lattice[i][j].numInfectedPsyllids_male;
                vector<double> probabilities;
                probabilities.reserve(60);
                vector<int> indices;
                indices.reserve(60);
                int totalShoots = 0;
                for (int l = 0; l < 30; l++) {
                    if (lattice[i][j].numShoots[l] > 0) {
                        probabilities.push_back(lattice[i][j].numShoots[l]);
                        indices.push_back(l);
                        totalShoots += lattice[i][j].numShoots[l];
                    }
                    if (lattice[i][j].numInfectedShoots[l] > 0) {
                        probabilities.push_back(lattice[i][j].numInfectedShoots[l]);
                        indices.push_back(l + 30);
                        totalShoots += lattice[i][j].numInfectedShoots[l];
                    }
                }
                for (int l = 0; l < probabilities.size(); l++) {
                    probabilities[l] = probabilities[l] / (double)totalShoots;
                }

                for (int l = 0; l < numInfectedPsyllids; l++) {
                    int shootIdx = discreteProbabilityMatch(probabilities);
                    if (shootIdx < 30) {
                        if (lattice[i][j].numShoots[shootIdx] > 0) {
                            if (doubleRand(0, 1) <= transmissionAdultFlush) {
                                lattice[i][j].numShoots[shootIdx]--;
                                lattice[i][j].numInfectedShoots[shootIdx]++;
                            }
                        }
                    }
                }
            }
        }
    }
}



void ageFlush_parallel() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            for (int k = 29; k >= 0; k--) {
                //Oldest shoots graduate to old shoots
                if (k == 29) {
                    lattice[i][j].oldUninfectedShoots += lattice[i][j].numShoots[k];
                    lattice[i][j].oldInfectedShoots += lattice[i][j].numInfectedShoots[k];
                }

                // All ages shift forward 1
                if (k == 0) {
                    lattice[i][j].numShoots[k] = 0;
                    lattice[i][j].numInfectedShoots[k] = 0;
                }
                else {
                    lattice[i][j].numShoots[k] = lattice[i][j].numShoots[k - 1];
                    lattice[i][j].numInfectedShoots[k] = lattice[i][j].numInfectedShoots[k - 1];
                }
            }
        }
    }
}


void initializeCSV() {
    try {
        csvFile.open(csvName);
    }
    catch (exception e) {
        cout << e.what() << endl;
    }
    csvFile << "t,i,j,numPsyllids,numInfectedPsyllids,hlbSeverity\n";
}

void write_csv_batch(int t) {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            int numPsyllids = lattice[i][j].numPsyllids_male + 
                lattice[i][j].numPsyllids_female + 
                accumulate(lattice[i][j].numNymphs.begin(), lattice[i][j].numNymphs.end(), 0);
            int numInfectedPsyllids = lattice[i][j].numInfectedPsyllids_male 
                + lattice[i][j].numInfectedPsyllids_female +
                accumulate(lattice[i][j].numInfectedNymphs.begin(), lattice[i][j].numInfectedNymphs.end(), 0);
            double severity = lattice[i][j].getHLBSeverity();
            if (severity > 1) {
                cout << "WARNING: HLB ERROR AT (" << i << ", " << j << ")\n";
            }
            csvFile << t << "," << i << "," << j << "," << numPsyllids << "," 
                << numInfectedPsyllids << "," << std::fixed 
                << setprecision(5) << severity << "," << experimentID << "\n";
        }
    }
}

/*
void writeSQLBatch(int t) {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            int numPsyllids = lattice[i][j].numPsyllids_male + 
                lattice[i][j].numPsyllids_female + 
                accumulate(lattice[i][j].numNymphs.begin(), lattice[i][j].numNymphs.end(), 0);
            int numInfectedPsyllids = lattice[i][j].numInfectedPsyllids_male 
                + lattice[i][j].numInfectedPsyllids_female +
                accumulate(lattice[i][j].numInfectedNymphs.begin(), lattice[i][j].numInfectedNymphs.end(), 0);
            double severity = lattice[i][j].getHLBSeverity();
            if (severity > 1) {
                cout << "WARNING: HLB ERROR AT (" << i << ", " << j << ")\n";
            }
            sql::Statement* stmt;
            stmt = con->createStatement();
            std::stringstream cmd;
            cmd << "INSERT INTO bio VALUES("
                << t << "," << i << "," << j << "," << numPsyllids << "," 
                << numInfectedPsyllids << "," << std::fixed 
                << setprecision(5) << severity << "," << experimentID << ");";
            stmt->execute(cmd.str());
            delete stmt;
        }
    }
} 
*/
void parseParameterFile(string fileName) {
    ifstream is(fileName);
    cereal::JSONInputArchive archive(is);
    try {
        archive(maxFlushAge, flushEmerging, eggAdultTransition,
            durationYoungFlush, proportionMigrating, withinRowP,
            betweenRowP, eggDuration, nymphDuration, shootCapacity,
            shootEggCapacity, eggsPerFemaleAdult, transmissionFlushNymph,
            transmissionAdultFlush, latentPeriod, eggSurvivalP, adultSurvivalP,
            nymphMinAgeToInfect, nymphMinAgeToBeInfected,
            modelDuration, csvName, initialInfectedPortion, initialNumPsyllids,
            invasionDay, carryingCapacity, borderCrossingP, springFlushStart,
            springFlushEnd, summerFlushStart, summerFlushEnd, fallFlushStart,
            fallFlushEnd, invasionModality, invasionGrove);
    }
    catch (exception e) {
        cout << "ERROR WITH BIO JSON" << endl;
        exit(-1);
    }
}


int countPsyllids() {
    int num = 0;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            num += lattice[i][j].getTotalPsyllids();
        }
    }
    return num;
}

bool validateLattice() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            if (!lattice[i][j].validate()) {
                cout << "Failed validation at (" << i << ", " << j << ")\n";
                return false;
            }
        }
    }
    return true;
}


void advanceBiologicalModel() {
    if (modelStarted && modelDay == -1) {
        throw("ERROR: advanceBiologicalModel and runModel cannot be used together.");
    }
    if (!modelStarted && modelDay == -1) {
        modelStarted = true;
        initializeModel();
        initializeCSV();
        modelDay = 0;
    }
    if (modelDay == invasionDay) {
        placeInitialPsyllids(invasionModality, invasionGrove);
    }
    if (modelDay >= invasionDay) {
        //cout << countFlush() << " flush\n";
        if (modelDay >= invasionDay) {
            setFlushingPeriod(modelDay);
            if (isFlushingPeriod) {
                //cout << "Birthing Flush\n";
                birthNewFlush();
            }
            migration_parallel();
            eggManagement_parallel();
            diseaseTransmission_parallel();
            ageFlush_parallel();
            psyllidAging_parallel();
        }
    }
    write_csv_batch(modelDay);
    modelDay++;
}


void runModelTest() {
    for (int i = 0; i < modelDuration; i++) {
        advanceBiologicalModel();
    }
}


void finishRun() {
    //csvFile.close();
}


int main(int argc, char *argv[]) {
    if (argc == 2) {
        parseParameterFile(string(argv[1]));
    }
    runModelTest();
    return 0;
}

} //end namespace
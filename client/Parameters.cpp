/*
 *===============================================
 *  Created by Andreas Pedersen on 10/30/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#include "Parameters.h"
#include "INIFile.h"

//using namespace constants;

Parameters::Parameters(){

    // Default values
    randomSeed_ = -1;
    reactantStateTag_ = 0;
    potentialTag_ = 1;
    potentialNoTranslation_ = 0;
    getPrefactorsTag_ = 0;
    typePerturbation_SP_ = 1;
    refine_SP_=false;
    lowestEigenmodeDetermination_SP_ = 1;
    minimize_only_ = 0;
    minimize_box_ = 0;

    job_Type_ = PROCESS_SEARCH;
    
    
    // Tweakable parameters, default values if not read in from parameters_passed.dat

    // Value used in the Relaxation   
    converged_Relax_ = 0.005;

    cgCurvatureStep_ = 0.001;
    cgMaxMoveFullRelax_ = 0.2;
    qmTimeStep_ = 0.1;
    maxDifferencePos_ = 0.1;
    neighborCutoff_ = 0.33;

    // Values used in the Saddle Point determination   
    converged_SP_ = 0.025;
    maxJumpAttempts_SP_ = 0;
    nrOfTriesToDetermineSaddlePoint_SP_ = 1;
    maxStepSize_SP_ = 0.2;
//    maxStepSizeConcave_SP_ = 0.2;
//    maxStepSizeConvex_SP_ = 0.1;
    maxEnergy_SP_ = 20.0;
    normPerturbation_SP_ = 0.1;
    withinRadiusPerturbated_SP_ = 4.0;
    maxSinglePerturbation_SP_ = 0.1;
    maximumIterations_ = 512;

    maxIterations_ = 256;
    maxIterationsConcave_ = 256;

    perpendicularForceRatio_=0.0;

    // Values used in the Hessian determination   
	maxSize_Hessian_ = 0;
	minDisplacement_Hessian_ = 0.25;
    withinRadiusDisplaced_Hessian_ = 5.0;

    prefactorMax_ = 10e20;
    prefactorMin_ = 10e8;

    // Values used in the Dimer method
    rotations_Dimer_ = 1;
    //rotationsNewSearch_Dimer_ = 10; // GH: not using this any more
    separation_Dimer_ = 0.0001;
    rotationAngle_Dimer_ = 0.005;
    
    // Defaults for Lanczos
    iterationLimit_Lanczos_=50;
    convergenceLimit_Lanczos_=1e-4;

    // Initializing the cummulative output
    forceCalls_ = 0;
    forceCallsSaddlePointConcave_ = 0;
    forceCallsSaddlePointConvex_ = 0;
    forceCallsPrefactors_ = 0;
    displacement_saddle_distance_ = 0;

    // Constants used in the client
    maxDifferencePos_ = 0.1; // The distance criterion for comparing geometries

    // Constants used in prefactor determination
    prefactorMax_ = 10e20; // max prefactor allowed
    prefactorMin_ = 10e8; // min prefactor allowed

    // Constants used by the optimizers
    cgCurvatureStep_ = 0.001; // finite difference step size used in conjugate gradients
    cgMaxMoveFullRelax_ = 0.2; // maximum displacement vector for a step during minimization
    qmTimeStep_ = 0.1; // time step size used in Quickmin.


    return;
}


Parameters::~Parameters(){
    return;
}

void Parameters::load(string filename)
{
    FILE *parametersFile;

    parametersFile = fopen(filename.c_str(), constants::READ.c_str());
    load(parametersFile);
    fclose(parametersFile);
}


void Parameters::load(FILE *file){
    
    //If we fail to parse the file as an INI, we will need to rewind the
    //file. So, we store the current stream position.
    fpos_t pos;
    fgetpos(file, &pos); 

    CIniFile ini;
    ini.CaseInsensitive();
    if(ini.ReadFile(file))
    {
        //If we succesfully read the file, then parse it as an INI
        randomSeed_ = ini.GetValueL("Default", "RANDOM_SEED", randomSeed_);
        reactantStateTag_ = ini.GetValueL("Default", "REACTANT_STATE_TAG", reactantStateTag_);
        potentialTag_ = ini.GetValueL("Default", "POTENTIAL_TAG", potentialTag_);
        potentialNoTranslation_ = ini.GetValueL("Default", "POTENTIAL_NO_TRANSLATION", potentialNoTranslation_);
        minimize_only_ = ini.GetValueL("Default", "MINIMIZE_ONLY", minimize_only_);
        minimize_box_ = ini.GetValueL("Default", "MINIMIZE_BOX", minimize_box_);
        getPrefactorsTag_ = ini.GetValueL("Default", "GET_PREFACTORS_TAG", getPrefactorsTag_);
        converged_Relax_ = ini.GetValueF("Default", "CONVERGED_RELAX", converged_Relax_);
        maximumIterations_ = ini.GetValueL("Default", "MAXIMUM_ITERATIONS", maximumIterations_);
        perpendicularForceRatio_ = ini.GetValueL("Default", "PERPENDICULAR_FORCE_RATIO", perpendicularForceRatio_);

        string job_type_string;
        job_type_string = ini.GetValue("Default", "JOB_TYPE", "ProcessSearch");

        if (job_type_string == "ProcessSearch") {
            job_Type_ = PROCESS_SEARCH;
        }else if (job_type_string == "SaddleSearch") {
            job_Type_ = SADDLE_SEARCH;
        }else if (job_type_string == "Minimization") {
            job_Type_ = MINIMIZATION;
        }else{
            job_Type_ = UNKNOWN_JOBTYPE;
        }
        
        
        typePerturbation_SP_ = ini.GetValueL("Saddle_Point", "TYPE_PERTURBATION", typePerturbation_SP_);
        lowestEigenmodeDetermination_SP_ = ini.GetValueL("Saddle_Point", "LOWEST_EIGENMODE_DETERMINATION", lowestEigenmodeDetermination_SP_);
        refine_SP_ = ini.GetValueB("Saddle_Point", "REFINE", refine_SP_);
        converged_SP_ = ini.GetValueF("Saddle_Point", "CONVERGED", converged_SP_);
        maxJumpAttempts_SP_ = ini.GetValueL("Saddle_Point", "MAX_JUMP_ATTEMPTS", maxJumpAttempts_SP_);
        nrOfTriesToDetermineSaddlePoint_SP_ = ini.GetValueL("Saddle_Point", "NR_OF_TRIES_TO_DETERMINE_SADDLE_POINT", nrOfTriesToDetermineSaddlePoint_SP_);
        maxStepSize_SP_ = ini.GetValueF("Saddle_Point", "MAX_STEP_SIZE", maxStepSize_SP_);
        maxEnergy_SP_ = ini.GetValueF("Saddle_Point", "MAX_ENERGY", maxEnergy_SP_);
        normPerturbation_SP_ = ini.GetValueF("Saddle_Point", "NORM_PERTURBATION", normPerturbation_SP_);
        maxSinglePerturbation_SP_ = ini.GetValueF("Saddle_Point", "MAX_SINGLE_PERTURBATION", maxSinglePerturbation_SP_);
        withinRadiusPerturbated_SP_ = ini.GetValueF("Saddle_Point", "WITHIN_RADIUS_PERTURBATED", withinRadiusPerturbated_SP_);


        maxSize_Hessian_ = ini.GetValueL("Hessian", "MAX_SIZE", maxSize_Hessian_);
        withinRadiusDisplaced_Hessian_ = ini.GetValueF("Hessian", "WITHIN_RADIUS_DISPLACED", withinRadiusDisplaced_Hessian_);
        minDisplacement_Hessian_ = ini.GetValueF("Hessian", "MIN_DISPLACEMENT", minDisplacement_Hessian_);

        
        rotations_Dimer_ = ini.GetValueL("Dimer", "ROTATIONS", rotations_Dimer_);
        separation_Dimer_ = ini.GetValueF("Dimer", "SEPARATION", separation_Dimer_);

        convergenceLimit_Lanczos_ = ini.GetValueF("Lanczos", "CONVERGENCE", convergenceLimit_Lanczos_);
        iterationLimit_Lanczos_ = ini.GetValueF("Lanczos", "ITERATIONS", iterationLimit_Lanczos_);        
    }
    else
    {
        //Otherwise, parse it as an old-style parameters file
        printf("Reading old-style parameters file\n"); 
        //Rewind the file to before CIniFile::ReadFile was called
        fsetpos(file, &pos);

        char **parms;
        double *values;
        long nLines, i;
        
        nLines = linesInFile(file);
        
        values = new double[nLines];
        parms = new char*[nLines];
        for(int i=0; i<nLines; i++)
            parms[i] = new char[STRING_SIZE];
        i = 0;

        // Important that the variable value is float
        // The used function sscanf only support this type 
        float value;
        char parm[STRING_SIZE];
        char lineAll[STRING_SIZE];
        for (i=0; i<nLines; i++){
            fgets(lineAll, STRING_SIZE, file);
            std::sscanf(lineAll, "%s %e", parm, &value);            
            for(int j = 0; parm[j]; j++)
                parm[j] = toupper(parm[j]);
                
            values[i] = value;
            strcpy(parms[i], parm);
        }
        
        for(i=0; i<nLines; i++){
            // Note strcmp() return 0 if strings are equal
            if(!strcmp(parms[i], "RANDOM_SEED"))
                randomSeed_ = long(values[i]);
            else if(!strcmp(parms[i], "REACTANT_STATE_TAG"))
                reactantStateTag_ = long(values[i]);
            else if(!strcmp(parms[i], "POTENTIAL_TAG"))
                potentialTag_ = long(values[i]);
            else if(!strcmp(parms[i], "POTENTIAL_NO_TRANSLATION"))
                potentialNoTranslation_ = long(values[i]);
            else if(!strcmp(parms[i], "MINIMIZE_ONLY"))
                minimize_only_ = long(values[i]);
            else if(!strcmp(parms[i], "MINIMIZE_BOX"))
                minimize_box_ = long(values[i]);
            else if(!strcmp(parms[i], "GET_PREFACTORS_TAG"))
                getPrefactorsTag_ = long(values[i]);
            else if(!strcmp(parms[i], "TYPE_PERTURBATION_SP"))
                typePerturbation_SP_ = (long) values[i];
            else if(!strcmp(parms[i], "LOWEST_EIGENMODE_DETERMINATION_SP"))
                lowestEigenmodeDetermination_SP_ = long(values[i]);
            else if(!strcmp(parms[i], "REFINE_SP")) {
                refine_SP_ = (bool) values[i];
            }

            // Tweakable parameters

            // Relaxation related
            else if(!strcmp(parms[i], "CONVERGED_RELAX"))
                converged_Relax_ = values[i];

            // Saddle Point related
            else if(!strcmp(parms[i], "CONVERGED_SP"))
                converged_SP_ = values[i];
            else if(!strcmp(parms[i], "MAX_JUMP_ATTEMPTS_SP"))
                maxJumpAttempts_SP_ = (long) values[i];
            else if(!strcmp(parms[i], "NR_OF_TRIES_TO_DETERMINE_SADDLE_POINT_SP"))
                nrOfTriesToDetermineSaddlePoint_SP_ = long(values[i]);
            else if(!strcmp(parms[i], "MAX_STEP_SIZE_SP"))
                maxStepSize_SP_ = values[i];
/*            else if(!strcmp(parms[i], "MAX_STEP_SIZE_CONCAVE_SP"))
                maxStepSizeConcave_SP_ = values[i];
            else if(!strcmp(parms[i], "MAX_STEP_SIZE_CONVEX_SP"))
                maxStepSizeConvex_SP_ = values[i]; */
            else if(!strcmp(parms[i], "MAX_ENERGY_SP"))
                maxEnergy_SP_ = values[i];
            else if(!strcmp(parms[i], "NORM_PERTURBATION_SP"))
                normPerturbation_SP_ = values[i];
            else if(!strcmp(parms[i], "WITHIN_RADIUS_PERTURBATED_SP"))
                withinRadiusPerturbated_SP_ = values[i];
            else if(!strcmp(parms[i], "MAX_SINGLE_PERTURBATION_SP"))
                maxSinglePerturbation_SP_ = values[i];
            else if(!strcmp(parms[i], "MAXIMUM_ITERATIONS"))
                maximumIterations_ = (long)values[i];
            else if(!strcmp(parms[i], "PERPENDICULAR_FORCE_RATIO"))
                perpendicularForceRatio_ = (double)values[i];

            // Hessian related
            else if(!strcmp(parms[i], "MAX_SIZE_HESSIAN"))
                maxSize_Hessian_ = values[i];
            else if(!strcmp(parms[i], "MIN_DISPLACEMENT_HESSIAN"))
                minDisplacement_Hessian_ = values[i];
            else if(!strcmp(parms[i], "WITHIN_RADIUS_DISPLACED_HESSIAN"))
                withinRadiusDisplaced_Hessian_ = values[i];

            // Dimer related
            else if(!strcmp(parms[i], "ROTATIONS_DIMER"))
                rotations_Dimer_ = (long) values[i];
/*            else if(!strcmp(parms[i], "ROTATIONS_NEW_SEARCH_DIMER"))
                rotationsNewSearch_Dimer_ = (long) values[i];*/
            else if(!strcmp(parms[i], "SEPARATION_DIMER"))
                separation_Dimer_ = (double) values[i];

            else if(!strcmp(parms[i], "LANCZOS_CONVERGENCE"))
                convergenceLimit_Lanczos_ = (double) values[i];
            else if(!strcmp(parms[i], "LANCZOS_ITERATION"))
                iterationLimit_Lanczos_ = (double) values[i];
            // Lines with user comment are started with #
            else if(parms[i][0]=='#'){}
            else
                std::cout<<"Unknown property: "<<parms[i]<<"\n";
        }
        delete [] values;
        for(i=0; i<nLines; i++)
            delete [] parms[i];
        delete [] parms;
    }
    return;
}

// Passing the input parameters
long Parameters::getRandomSeed(){
    return randomSeed_;
}
void Parameters::setRandomSeed(long randomSeed){
    randomSeed_=randomSeed;
}
long Parameters::getReactantStateTag(){
    return reactantStateTag_;
}
long Parameters::getPotentialTag(){
    return potentialTag_;
}
long Parameters::getPotentialNoTranslation(){
    return potentialNoTranslation_;
}
bool Parameters::getMinimizeOnly(){
    return minimize_only_;
}
bool Parameters::getMinimizeBox(){
    return minimize_box_;
}
long Parameters::getPrefactorsTag(){
    return getPrefactorsTag_;
}
// Tweakable parameters

// Relaxation related
double Parameters::getConverged_Relax(){
    return converged_Relax_;
}

// Optimizer related
double  Parameters::getCgCurvatureStep(){
    return cgCurvatureStep_;
}
double  Parameters::getCgMaxMoveFullRelax(){
    return cgMaxMoveFullRelax_;
}
double  Parameters::getQmTimeStep(){
    return qmTimeStep_;
}

// Matter related
double  Parameters::getMaxDifferencePos(){ 
    return maxDifferencePos_;
}

// EpiCenter related
double  Parameters::getNeighborCutoff(){ 
    return neighborCutoff_;
}

// Saddle Point determination related
double Parameters::getConverged_SP(){
    return converged_SP_;
}
long Parameters::getMaxJumpAttempts_SP(){
    return maxJumpAttempts_SP_;
}
long Parameters::getTypePerturbation_SP(){
    return typePerturbation_SP_;
}
bool Parameters::getRefineSP(){
    return refine_SP_;
}
long Parameters::getLowestEigenmodeDetermination_SP(){
    return lowestEigenmodeDetermination_SP_;
}
long Parameters::getNrOfTriesToDetermineSaddlePoints_SP(){
    return nrOfTriesToDetermineSaddlePoint_SP_;
}
double Parameters::getMaxStepSize_SP(){
    return maxStepSize_SP_;
}
/*double Parameters::getMaxStepSizeConcave_SP(){
    return maxStepSizeConcave_SP_;
}
double Parameters::getMaxStepSizeConvex_SP(){
    return maxStepSizeConvex_SP_;
}*/
double Parameters::getMaxEnergy_SP(){
    return maxEnergy_SP_;
}
double Parameters::getNormPerturbation_SP(){
    return normPerturbation_SP_;
}
double Parameters::getWithinRadiusPerturbated_SP(){
    return withinRadiusPerturbated_SP_;
}
double Parameters::getMaxSinglePerturbation_SP(){
    return maxSinglePerturbation_SP_;
}
/// Limit on the number of iterations that may be performed by the saddle point searches and minimization
long Parameters::getMaximumIterations(){
    return maximumIterations_;
}
double Parameters::getPerpendicularForceRatio(){
    return perpendicularForceRatio_;
}
long Parameters::getMaxIterations(){
    return maxIterations_;
}
long Parameters::getMaxIterationsConcave(){
    return maxIterationsConcave_;
}

// Hessian related
long Parameters::getMaxSize_Hessian(){
    return maxSize_Hessian_;
}
double Parameters::getMinDisplacement_Hessian(){
    return minDisplacement_Hessian_;
}
double Parameters::getWithinRadiusDisplaced_Hessian(){
    return withinRadiusDisplaced_Hessian_;
}
double Parameters::getPrefactorMax(){
    return prefactorMax_;
}
double Parameters::getPrefactorMin(){
    return prefactorMin_;
}

// Dimer related
long Parameters::getRotations_Dimer(){
    return rotations_Dimer_;
}
/*long Parameters::getRotationsNewSearch_Dimer(){
    return rotationsNewSearch_Dimer_;
}*/
double Parameters::getSeparation_Dimer(){
    return separation_Dimer_;
}
double Parameters::getRotationAngle_Dimer(){
    return rotationAngle_Dimer_;
}

// Lanczos related
double Parameters::getConvergenceLimit_Lanczos(){
      return convergenceLimit_Lanczos_;
}
double Parameters::getIterationLimit_Lanczos(){
      return iterationLimit_Lanczos_;
}

// Setting results values
void Parameters::setTerminationReason(long terminationReason){
    terminationReason_ = terminationReason;
    return;
}
void Parameters::setPotentialEnergySP(double potentialEnergySP){
    potentialEnergySP_ = potentialEnergySP;
    return;
}
void Parameters::setPotentialEnergyMin1(double potentialEnergyMin1){
    potentialEnergyMin1_ = potentialEnergyMin1;
    return;
}
void Parameters::setPotentialEnergyMin2(double potentialEnergyMin2){
    potentialEnergyMin2_ = potentialEnergyMin2;
    return;
}
void Parameters::setPrefactorReac_Prod(double prefactorReac_Prod){
    prefactorReac_Prod_ = prefactorReac_Prod;
    return;
}
void Parameters::setPrefactorProd_Reac(double prefactorProd_Reac){
    prefactorProd_Reac_ = prefactorProd_Reac;
    return;
}
void Parameters::setBarrierReac_Prod(double barrierReac_Prod){
    barrierReac_Prod_ = barrierReac_Prod;
    return;
}
void Parameters::setBarrierProd_Reac(double barrierProd_Reac){
    barrierProd_Reac_ = barrierProd_Reac;
    return;
}
void Parameters::setDisplacementSaddleDistance(double displacement_saddle_distance){
    displacement_saddle_distance_ = displacement_saddle_distance;
    return;
}

long Parameters::linesInFile(FILE *file){
    int nLines = 0;
    char ch, prev = '\n'; 
    // read all chars in the file
    while((ch = fgetc(file)) != EOF){
        if(ch == '\n'){
            nLines = nLines+1;
        }
        prev = ch; // keep a copy to later test whether
    }
    if(prev!='\n') // ...the last line did not end in a newline.
        nLines = nLines+1;

    rewind(file);
    return (nLines);
}


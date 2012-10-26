//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "PointJob.h"
#include "Matter.h"

PointJob::PointJob(Parameters *params)
{
    parameters = params;
}

PointJob::~PointJob(){ }

std::vector<std::string> PointJob::run(void)
{
    std::vector<std::string> returnFiles;
    string posInFilename("pos.con");
    string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);

    Matter *pos = new Matter(parameters);
    pos->con2matter(posInFilename);

    printf("Energy:         %f\n", pos->getPotentialEnergy());
    printf("Max atom force: %g\n", pos->maxForce());

    FILE *fileResults = fopen(resultsFilename.c_str(), "wb");
    fprintf(fileResults, "%f Energy\n", pos->getPotentialEnergy());
    fprintf(fileResults, "%f Max_Force\n", pos->maxForce());
    fclose(fileResults);

    return returnFiles;
}

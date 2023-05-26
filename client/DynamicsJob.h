#ifndef DYNAMICSJOB_H
#define DYNAMICSJOB_H

#include "Matter.h"
#include "Parameters.h"
#include "Job.h"

class DynamicsJob : public Job 
{

    public:

        DynamicsJob(std::unique_ptr<Parameters> parameters)
            : Job(std::move(parameters)) {}
        ~DynamicsJob(void) = default;
        std::vector<std::string> run(void);
        std::vector<std::string> returnFiles;
};

#endif

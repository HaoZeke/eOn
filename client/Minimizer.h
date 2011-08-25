//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#ifndef MINIMIZER_H
#define MINIMIZER_H

class Minimizer
{
    
    public:
        virtual ~Minimizer(){};
        void virtual setOutput(int level) = 0;
        void virtual oneStep() = 0;
        long virtual fullRelax() = 0;
        bool virtual isItConverged(double convergeCriterion) = 0; 
        long totalForceCalls;

        enum {
            STATUS_GOOD, //0
            STATUS_MAX_ITERATIONS, //1
            STATUS_POTENTIAL_FAILED, //2
            STATUS_CHECKPOINT, //3
        };
};

#endif
